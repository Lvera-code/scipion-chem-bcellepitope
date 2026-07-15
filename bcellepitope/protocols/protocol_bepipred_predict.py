"""Protocolo Scipion-chem: Fases 2+3 (BepiPred-3.0) del pipeline
B-Cell-Epitope-Prediction.

Ejecuta BepiPred-3.0 localmente contra el CLI oficial via ``runJob`` (mismo
patron de subprocess que ``B-Cell-Epitope-Prediction/src/engines/bepipred_engine.py``,
sin activacion de entorno conda: se invoca directamente el interprete
configurado en ``BEPIPRED_PYTHON_BIN``), y aplica sobre el ``raw_output.csv``
resultante la misma ventana deslizante tolerante a gaps ya validada en el
pipeline original (``bcellepitope.utils.epitope_mapping``), sin depender de
los propios ficheros de prediccion de BepiPred (que se ignoran, igual que en
el pipeline original). El resultado se expone como ``SetOfSequenceROIs``,
compatible con el resto de protocolos de Scipion-chem.
"""

import os
from pathlib import Path

import pandas as pd
from pwchem.objects import Sequence, SequenceROI, SetOfSequenceROIs
from pwem.protocols import EMProtocol
from pyworkflow.protocol import params

from .. import Plugin as bcellPlugin
from ..utils.epitope_mapping import extract_epitope_regions
from ..utils.exceptions import BepiPredExecutionError

# Columnas confirmadas empiricamente en el raw_output.csv real de BepiPred-3.0
# (ver B-Cell-Epitope-Prediction/src/engines/bepipred_engine.py).
ACCESSION_COLUMN = "Accession"
SCORE_COLUMN = "BepiPred-3.0 score"
RESIDUE_COLUMN_CANDIDATES = ("Residue", "residue", "Residues", "AA", "aa")
RAW_OUTPUT_FILENAME = "raw_output.csv"


class ProtBCellEpitopeBepiPredPredict(EMProtocol):
    """Prediccion de antigenicidad B-cell con BepiPred-3.0 (ejecucion local) y
    extraccion de regiones de epitopo por ventana deslizante tolerante a gaps."""

    _label = 'bepipred antigenicity prediction'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSequence', params.PointerParam, pointerClass='Sequence',
                       label='Input protein sequence: ',
                       help='Protein sequence to run BepiPred-3.0 antigenicity prediction on.')

        form.addParam('predMode', params.EnumParam, label='Prediction mode: ',
                       choices=['vt_pred', 'mjv_pred'], default=0,
                       expertLevel=params.LEVEL_ADVANCED,
                       help="Flag obligatorio del CLI de BepiPred-3.0 (variable threshold vs. "
                            "majority vote). No afecta la extraccion de epitopos de este "
                            "protocolo, que corre sobre 'raw_output.csv' y no sobre los "
                            "propios ficheros de prediccion de BepiPred (que se ignoran).")

        eGroup = form.addGroup('Epitope extraction')
        eGroup.addParam('threshold', params.FloatParam, label='Score threshold: ', default=0.1512,
                         help='Umbral de score de BepiPred-3.0 para considerar un residuo/ventana '
                              'como parte de un epitopo candidato.')
        eGroup.addParam('minLength', params.IntParam, label='Minimum epitope length: ', default=9,
                         help='Longitud minima (aa) de una region fusionada para reportarse como epitopo.')
        eGroup.addParam('windowSize', params.IntParam, label='Sliding window size: ', default=9,
                         expertLevel=params.LEVEL_ADVANCED,
                         help='Footprint minimo de reconocimiento de celula B (9 aa por defecto).')
        eGroup.addParam('maxGapResidues', params.IntParam, label='Max gap residues per window: ', default=2,
                         expertLevel=params.LEVEL_ADVANCED,
                         help='Residuos individuales por debajo del umbral tolerados dentro de una '
                              'misma ventana, para no descartar un epitopo real por un unico residuo debil.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.bepipredStep)
        self._insertFunctionStep(self.createOutputStep)

    # ---------------------------------- Steps -----------------------------------

    def writeInputFasta(self):
        faFile = self._getExtraPath('inputSequence.fa')
        self.inputSequence.get().exportToFile(faFile)
        return os.path.abspath(faFile)

    def bepipredStep(self):
        bepipred_home = bcellPlugin.getBepipredHome()
        python_bin = bcellPlugin.getBepipredPythonBin()
        cli_script = bcellPlugin.getBepipredCliScript()

        fasta_file = self.writeInputFasta()
        out_dir = os.path.abspath(self._getExtraPath('bepipred_raw'))
        os.makedirs(out_dir, exist_ok=True)

        args = (
            f'{cli_script} -i {fasta_file} -o {out_dir} '
            f'-pred {self.getEnumText("predMode")} -t {self.threshold.get()}'
        )
        self.runJob(python_bin, args, cwd=bepipred_home)

    def createOutputStep(self):
        out_dir = self._getExtraPath('bepipred_raw')
        csv_path = self._locateRawOutput(out_dir)
        raw_df = pd.read_csv(csv_path)

        missing = {ACCESSION_COLUMN, SCORE_COLUMN} - set(raw_df.columns)
        if missing:
            raise BepiPredExecutionError(
                f"El CSV de salida '{csv_path}' no contiene las columnas confirmadas "
                f"{sorted(missing)}. Columnas encontradas: {list(raw_df.columns)}."
            )

        epitopes_df = extract_epitope_regions(
            raw_df,
            accession_col=ACCESSION_COLUMN,
            score_col=SCORE_COLUMN,
            residue_col_candidates=RESIDUE_COLUMN_CANDIDATES,
            threshold=self.threshold.get(),
            min_length=self.minLength.get(),
            window_size=self.windowSize.get(),
            max_gap_residues=self.maxGapResidues.get(),
        )

        inputSeq = self.inputSequence.get()
        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for row in epitopes_df.itertuples(index=False):
            roiId = f'ROI_{row.start}-{row.end}'
            roiSeq = Sequence(sequence=row.sequence, name=roiId, id=roiId,
                               description='BepiPred-3.0 epitope')
            seqROI = SequenceROI(sequence=inputSeq, seqROI=roiSeq, roiIdx=row.start, roiIdx2=row.end)
            outROIs.append(seqROI)

        if len(outROIs) > 0:
            self._defineOutputs(outputROIs=outROIs)
            self._defineSourceRelation(self.inputSequence, outROIs)

    # ---------------------------------- Utils -----------------------------------

    @staticmethod
    def _locateRawOutput(out_dir: str) -> Path:
        result_dir = Path(out_dir)
        primary = list(result_dir.rglob(RAW_OUTPUT_FILENAME))
        if primary:
            return primary[0]

        candidates = sorted(result_dir.rglob('*.csv'))
        if not candidates:
            found = sorted(p.name for p in result_dir.rglob('*') if p.is_file())
            raise BepiPredExecutionError(
                f"No se encontro ningun CSV de salida en '{result_dir}'. "
                f"Archivos generados por BepiPred-3.0: {found or '<ninguno>'}."
            )
        return candidates[0]

    # ---------------------------------- Validation -------------------------------

    def _validate(self):
        return bcellPlugin.validateInstallation()

    def _summary(self):
        summary = []
        if self.isFinished():
            outROIs = getattr(self, 'outputROIs', None)
            n = len(outROIs) if outROIs is not None else 0
            summary.append(f'{n} epitope region(s) found above threshold {self.threshold.get()}.')
        return summary
