"""Protocolo Scipion-chem: Fase 4 (filtro de tolerancia inmunologica, BLASTp
local) del pipeline B-Cell-Epitope-Prediction.

Anota CADA ROI de entrada (no filtra: mismo criterio que blast_engine.py) con
'blast_task'/'blast_evalue'/'max_pident'/'status' ('Segura' o
'Autoinmunidad'), agrupando los peptidos por tramo de longitud (mismo
'-task'/E-value) para minimizar invocaciones de 'blastp' -ver
``bcellepitope.utils.blast`` para la seleccion dinamica de tramo, portada tal
cual de ``blast_engine.py::_select_task``/``_select_evalue``.
"""

import os
from pathlib import Path

import pandas as pd
from pwchem.objects import Sequence, SequenceROI, SetOfSequenceROIs
from pwem.protocols import EMProtocol
from pyworkflow.object import Float, String
from pyworkflow.protocol import params

from .. import Plugin as bcellPlugin
from ..utils.blast import (
    DEFAULT_EVALUE_LONG, DEFAULT_EVALUE_MEDIUM, DEFAULT_EVALUE_SHORT,
    DEFAULT_MEDIUM_PEPTIDE_MAX_LEN, DEFAULT_SHORT_PEPTIDE_MAX_LEN,
    select_evalue, select_task,
)
from ..utils.exceptions import BlastExecutionError

_OUTFMT6_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]


class ProtBCellEpitopeBlastFilter(EMProtocol):
    """Filtro de tolerancia inmunologica: descarta por BLASTp local los ROI con
    alta homologia al proteoma humano (los anota, no los elimina del set)."""

    _label = 'blastp immune tolerance filter'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSequence', params.PointerParam, pointerClass='Sequence',
                       label='Input protein sequence: ',
                       help='La misma secuencia de la que provienen los ROI a evaluar.')
        form.addParam('inputROIs', params.PointerParam, pointerClass='SetOfSequenceROIs',
                       label='Epitope ROIs: ',
                       help='SetOfSequenceROIs a evaluar contra el proteoma humano (tipicamente '
                            'la salida del protocolo de consenso, o de un motor individual).')

        eGroup = form.addGroup('Immune tolerance')
        eGroup.addParam('identityThreshold', params.FloatParam, label='Identity threshold (%): ', default=75.0,
                         help='Porcentaje de identidad (exclusivo) por encima del cual un ROI se '
                              'descarta por riesgo de autoinmunidad.')

        tGroup = form.addGroup('Length-dependent BLAST tiers', expertLevel=params.LEVEL_ADVANCED)
        tGroup.addParam('shortPeptideMaxLen', params.IntParam, label='Short peptide max length: ',
                         default=DEFAULT_SHORT_PEPTIDE_MAX_LEN, expertLevel=params.LEVEL_ADVANCED,
                         help="Longitud (aa) hasta la cual se usa '-task blastp-short'.")
        tGroup.addParam('mediumPeptideMaxLen', params.IntParam, label='Medium peptide max length: ',
                         default=DEFAULT_MEDIUM_PEPTIDE_MAX_LEN, expertLevel=params.LEVEL_ADVANCED,
                         help="Longitud (aa) hasta la cual se usa el E-value 'medium' (por encima, 'long').")
        tGroup.addParam('evalueShort', params.FloatParam, label='E-value (short): ',
                         default=DEFAULT_EVALUE_SHORT, expertLevel=params.LEVEL_ADVANCED,
                         help='E-value laxo para peptidos cortos: la estadistica de BLAST penaliza '
                              'su longitud y descartaria hits identicos reales.')
        tGroup.addParam('evalueMedium', params.FloatParam, label='E-value (medium): ',
                         default=DEFAULT_EVALUE_MEDIUM, expertLevel=params.LEVEL_ADVANCED)
        tGroup.addParam('evalueLong', params.FloatParam, label='E-value (long): ',
                         default=DEFAULT_EVALUE_LONG, expertLevel=params.LEVEL_ADVANCED,
                         help='E-value estricto para peptidos largos, evita ruido de homologias irrelevantes.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.blastStep)

    # ---------------------------------- Steps -----------------------------------

    def blastStep(self):
        inputSeq = self.inputSequence.get()
        # Scipion reutiliza el mismo objeto Python por cada fila al iterar un
        # SetOfXXX (cursor de sqlite subyacente): hay que clonar cada item al
        # materializar la lista, o las N referencias terminan apuntando todas
        # al ultimo estado del cursor.
        rois = [roi.clone() for roi in self.inputROIs.get()]

        tiers = {}
        for idx, roi in enumerate(rois):
            seq = roi.getROISequence()
            task = select_task(len(seq), self.shortPeptideMaxLen.get())
            evalue = select_evalue(
                len(seq), self.shortPeptideMaxLen.get(), self.mediumPeptideMaxLen.get(),
                self.evalueShort.get(), self.evalueMedium.get(), self.evalueLong.get(),
            )
            tiers.setdefault((task, evalue), []).append((idx, seq))

        maxPident = {}
        for (task, evalue), entries in tiers.items():
            hits = self._runBlastpBatch(entries, task, evalue)
            if not hits.empty:
                best = hits.groupby('qseqid')['pident'].max()
                for qseqid, pident in best.items():
                    maxPident[int(qseqid.split('_')[1])] = float(pident)

        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for idx, roi in enumerate(rois):
            task = select_task(roi.getROILength(), self.shortPeptideMaxLen.get())
            evalue = select_evalue(
                roi.getROILength(), self.shortPeptideMaxLen.get(), self.mediumPeptideMaxLen.get(),
                self.evalueShort.get(), self.evalueMedium.get(), self.evalueLong.get(),
            )
            pident = maxPident.get(idx, 0.0)
            status = 'Autoinmunidad' if pident > self.identityThreshold.get() else 'Segura'

            roiSeq = Sequence(sequence=roi.getROISequence(), name=roi.getROIName(), id=roi.getROIId(),
                               description=roi._ROISequence.getDescription())
            newRoi = SequenceROI(sequence=inputSeq, seqROI=roiSeq,
                                  roiIdx=roi.getROIIdx(), roiIdx2=roi.getROIIdx2())
            for attrName in ('_meanScore', '_origen', '_bepipredScore', '_epidopeScore',
                              '_bepipredRegion', '_epidopeRegion'):
                if hasattr(roi, attrName):
                    setattr(newRoi, attrName, getattr(roi, attrName).clone())
            newRoi._blastTask = String(task)
            newRoi._blastEvalue = Float(evalue)
            newRoi._maxPident = Float(pident)
            newRoi._status = String(status)
            outROIs.append(newRoi)

        if len(outROIs) > 0:
            self._defineOutputs(outputROIs=outROIs)
            self._defineSourceRelation(self.inputROIs, outROIs)

    def _runBlastpBatch(self, entries, task, evalue):
        blastp_bin = bcellPlugin.getBlastpBin()
        db = bcellPlugin.getBlastHumanDb()

        tier_dir = self._getExtraPath(f'{task}_{evalue}'.replace('.', '_'))
        os.makedirs(tier_dir, exist_ok=True)
        query_path = os.path.join(tier_dir, 'candidates.fasta')
        out_path = os.path.join(tier_dir, 'blast_results.tsv')

        with open(query_path, 'w') as fh:
            for idx, seq in entries:
                fh.write(f'>peptide_{idx}\n{seq}\n')

        args = f'-task {task} -query {query_path} -db {db} -outfmt 6 -evalue {evalue} -out {out_path}'
        self.runJob(blastp_bin, args)

        if not os.path.isfile(out_path) or os.path.getsize(out_path) == 0:
            return pd.DataFrame(columns=_OUTFMT6_COLUMNS)
        return pd.read_csv(out_path, sep='\t', names=_OUTFMT6_COLUMNS)

    # ---------------------------------- Validation -------------------------------

    def _validate(self):
        return bcellPlugin.validateBlastInstallation()

    def _summary(self):
        summary = []
        if self.isFinished():
            outROIs = getattr(self, 'outputROIs', None)
            if outROIs is not None:
                nSafe = sum(1 for roi in outROIs if roi._status.get() == 'Segura')
                nRejected = sum(1 for roi in outROIs if roi._status.get() == 'Autoinmunidad')
                summary.append(f'{nSafe} segura(s) / {nRejected} rechazada(s) por homologia con '
                                f'el proteoma humano (umbral {self.identityThreshold.get()}%).')
        return summary
