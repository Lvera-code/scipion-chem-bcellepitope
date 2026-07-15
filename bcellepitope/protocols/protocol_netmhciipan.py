"""Protocolo Scipion-chem: Fase 5 (ultima) del pipeline B-Cell-Epitope-Prediction.

Evalua promiscuidad T-helper (MHC-II) via NetMHCIIpan-4.3 local sobre los ROI
con ``status == 'Segura'`` (salida del protocolo de filtro BLASTp de este
plugin), enrutando cada peptido a modo peptido exacto o modo proteina segun
su longitud, y traza cada candidato de vuelta a su ROI padre. Ver
``bcellepitope.utils.netmhciipan`` para el parseo del .xls y el traceback
(portado de ``B-Cell-Epitope-Prediction/src/engines/netmhciipan_engine.py``).
"""

import os

import pandas as pd
from pwchem.objects import Sequence, SequenceROI, SetOfSequenceROIs
from pwem.protocols import EMProtocol
from pyworkflow.object import Float, Integer, String
from pyworkflow.protocol import params

from .. import Plugin as bcellPlugin
from ..utils.exceptions import NetMHCIIpanExecutionError
from ..utils.netmhciipan import (
    IEDB_REFERENCE_PANEL, MAX_PEPTIDE_MODE_LENGTH, MIN_PEPTIDE_LENGTH,
    build_traceback_report, parse_xls,
)


class ProtBCellEpitopeNetMHCIIpan(EMProtocol):
    """Promiscuidad T-helper (MHC-II) via NetMHCIIpan-4.3 local sobre los ROI
    'Segura' (filtro BLASTp), con traceback a la region de origen."""

    _label = 'netmhciipan th promiscuity'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSequence', params.PointerParam, pointerClass='Sequence',
                       label='Input protein sequence: ',
                       help='La secuencia de la que provienen los ROI a evaluar.')
        form.addParam('inputROIs', params.PointerParam, pointerClass='SetOfSequenceROIs',
                       label='BLAST-annotated epitope ROIs: ',
                       help="Salida del protocolo de filtro BLASTp de este plugin. Solo se "
                            "evaluan los ROI con status == 'Segura'.")
        form.addParam('allelePanel', params.StringParam, default=IEDB_REFERENCE_PANEL,
                       label='HLA-DR/DQ/DP allele panel: ',
                       help='Alelos separados por coma SIN espacios (formato NetMHCIIpan). Por '
                            'defecto, el panel de referencia IEDB de 27 alelos mas representativos '
                            'de cobertura poblacional. Para anexar un alelo especifico, edita este '
                            'campo directamente (ej. agregar ",DRB1_1602" al final).')

        gGroup = form.addGroup('Promiscuity')
        gGroup.addParam('rankWeak', params.FloatParam, label='Weak binder %Rank threshold: ', default=5.0,
                         help='Umbral de %Rank_EL de NetMHCIIpan por debajo (inclusive) del cual un '
                              'alelo cuenta como aglutinador (SB o WB) para la promiscuidad.')
        gGroup.addParam('minPromiscuousAlleles', params.IntParam, label='Min. promiscuous alleles: ', default=3,
                         help="Numero minimo de alelos del panel (en orientacion normal, ver ayuda del "
                              "modulo de utils) que deben clasificar como aglutinadores para reportar "
                              "un peptido como 'Candidato Valido'.")
        gGroup.addParam('timeoutSeconds', params.IntParam, label='Timeout (s): ', default=600,
                         expertLevel=params.LEVEL_ADVANCED)

    def _insertAllSteps(self):
        self._insertFunctionStep(self.netmhciipanStep)

    # ---------------------------------- Steps -----------------------------------

    def netmhciipanStep(self):
        # Ver ADR en protocol_blast.py: iterar un SetOfXXX de Scipion reutiliza
        # el mismo objeto Python por fila, hay que clonar al materializar.
        rois = [roi.clone() for roi in self.inputROIs.get()]
        segura = [roi for roi in rois if hasattr(roi, '_status') and roi._status.get() == 'Segura']

        parentRecords = []
        peptides = []
        for roi in segura:
            seq = roi.getROISequence()
            if len(seq) < MIN_PEPTIDE_LENGTH:
                continue
            peptides.append(seq)
            parentRecords.append({
                'sequence': seq,
                'start': roi.getROIIdx(),
                'origen': roi._origen.get() if hasattr(roi, '_origen') else '',
                'bepipred_score': roi._bepipredScore.get() if hasattr(roi, '_bepipredScore') else None,
                'epidope_score': roi._epidopeScore.get() if hasattr(roi, '_epidopeScore') else None,
            })

        if not peptides:
            return

        shortPeptides = [p for p in peptides if len(p) <= MAX_PEPTIDE_MODE_LENGTH]
        longPeptides = [p for p in peptides if len(p) > MAX_PEPTIDE_MODE_LENGTH]

        allelePanel = self.allelePanel.get()
        nAlleles = len([a for a in allelePanel.split(',') if a])
        binary = bcellPlugin.getNetMHCIIpanBinaryPath()

        reportFrames = []
        if shortPeptides:
            reportFrames.append(self._runMode(
                binary, allelePanel, nAlleles,
                modeArgs=['-p', '-f', self._writePeptideFile(shortPeptides)],
                xlsPath=self._getExtraPath('peptide_mode_output.xls'),
                modeDesc='modo peptido exacto',
            ))
        if longPeptides:
            reportFrames.append(self._runMode(
                binary, allelePanel, nAlleles,
                modeArgs=['-f', self._writeFragmentsFasta(longPeptides)],
                xlsPath=self._getExtraPath('protein_mode_output.xls'),
                modeDesc='modo proteina (ventana deslizante)',
            ))

        reportDf = pd.concat(reportFrames, ignore_index=True) if reportFrames else pd.DataFrame()
        tracebackDf = build_traceback_report(reportDf, parentRecords)

        inputSeq = self.inputSequence.get()
        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for row in tracebackDf.itertuples(index=False):
            roiId = f'ROI_{row.start}-{row.end}'
            roiSeq = Sequence(sequence=row.sequence_f5, name=roiId, id=roiId,
                               description='NetMHCIIpan T-helper candidate')
            newRoi = SequenceROI(sequence=inputSeq, seqROI=roiSeq, roiIdx=row.start, roiIdx2=row.end)
            newRoi._origen = String(row.origen)
            newRoi._core9aa = String(row.core_9aa)
            newRoi._nAlelosPromiscuos = Integer(int(row.n_alelos_promiscuos))
            newRoi._nAlelosEvaluados = Integer(int(row.n_alelos_evaluados))
            newRoi._minRankEl = Float(row.min_rank_el)
            if pd.notna(row.bepipred_score):
                newRoi._bepipredScore = Float(row.bepipred_score)
            if pd.notna(row.epidope_score):
                newRoi._epidopeScore = Float(row.epidope_score)
            outROIs.append(newRoi)

        if len(outROIs) > 0:
            self._defineOutputs(outputROIs=outROIs)
            self._defineSourceRelation(self.inputROIs, outROIs)

    def _writePeptideFile(self, peptides):
        pepPath = self._getExtraPath('peptides.pep')
        with open(pepPath, 'w') as fh:
            fh.write('\n'.join(peptides) + '\n')
        return pepPath

    def _writeFragmentsFasta(self, peptides):
        fastaPath = self._getExtraPath('fragments.fasta')
        with open(fastaPath, 'w') as fh:
            for i, seq in enumerate(peptides):
                fh.write(f'>candidato_{i}\n{seq}\n')
        return fastaPath

    def _runMode(self, binary, allelePanel, nAlleles, modeArgs, xlsPath, modeDesc):
        args = ' '.join(modeArgs) + f' -a {allelePanel} -xls -xlsfile {xlsPath}'
        self.runJob(binary, args)

        if not os.path.isfile(xlsPath):
            raise NetMHCIIpanExecutionError(
                f"NetMHCIIpan-4.3 ({modeDesc}) termino sin generar el archivo de salida esperado "
                f"en '{xlsPath}'. Causas conocidas: un peptido de entrada excede el limite del modo "
                "usado, o la linea NMHOME dentro del script wrapper apunta a una ruta desactualizada."
            )
        return parse_xls(xlsPath, nAlleles, self.rankWeak.get(), self.minPromiscuousAlleles.get())

    # ---------------------------------- Validation -------------------------------

    def _validate(self):
        return bcellPlugin.validateNetMHCIIpanInstallation()

    def _summary(self):
        summary = []
        if self.isFinished():
            outROIs = getattr(self, 'outputROIs', None)
            n = len(outROIs) if outROIs is not None else 0
            nAlleles = len([a for a in self.allelePanel.get().split(',') if a])
            summary.append(f'{n} candidato(s) T-helper promiscuo(s) (>= '
                            f'{self.minPromiscuousAlleles.get()} de {nAlleles} alelos) tras traceback '
                            'y deduplicacion.')
        return summary
