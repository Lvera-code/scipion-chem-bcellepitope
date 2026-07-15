"""Protocolo Scipion-chem: Fases 2+3 (EpiDope) del pipeline B-Cell-Epitope-Prediction.

Mismo patron que protocol_bepipred_predict.py, pero para el segundo motor de
antigenicidad (EpiDope, codigo abierto MIT). Diferencias clave frente a
BepiPred, portadas de epidope_engine.py:

* El binario 'epidope' (shim autocontenido en '<EPIDOPE_HOME>/bin/epidope')
  se invoca DIRECTAMENTE, sin 'conda run' ni activacion de entorno: el ADR en
  epidope_engine.py documenta que 'conda run --no-capture-output' producia
  fallos espurios intermitentes con stdout/stderr capturados por pipe.
* La salida cruda no es un unico CSV: EpiDope escribe un CSV por-accession en
  '<outdir>/epidope/<accession>.csv' (separado por tabs, columnas
  'position'/'aminoacid'/'score'), que hay que concatenar.
* El umbral por defecto de EpiDope (0.818) NO es comparable en escala al de
  BepiPred (0.1512): son motores independientes con sus propios parametros.
"""

import os
from pathlib import Path

import pandas as pd
from pwchem.objects import Sequence, SequenceROI, SetOfSequenceROIs
from pwem.protocols import EMProtocol
from pyworkflow.object import Float
from pyworkflow.protocol import params

from .. import Plugin as bcellPlugin
from ..utils.epitope_mapping import extract_epitope_regions
from ..utils.exceptions import EpidopeExecutionError

ACCESSION_COLUMN = "Accession"
RESIDUE_COLUMN = "Residue"
SCORE_COLUMN = "EpiDope score"

# Columnas confirmadas en <outdir>/epidope/<accession>.csv (ver
# epidope_engine.py, epidope/epidope2.py::output_results).
_RAW_COLUMNS = ("position", "aminoacid", "score")


class ProtBCellEpitopeEpiDopePredict(EMProtocol):
    """Prediccion de antigenicidad B-cell con EpiDope (ejecucion local) y
    extraccion de regiones de epitopo por ventana deslizante tolerante a gaps."""

    _label = 'epidope antigenicity prediction'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSequence', params.PointerParam, pointerClass='Sequence',
                       label='Input protein sequence: ',
                       help='Protein sequence to run EpiDope antigenicity prediction on.')

        eGroup = form.addGroup('Epitope extraction')
        eGroup.addParam('threshold', params.FloatParam, label='Score threshold: ', default=0.818,
                         help='Umbral de score de EpiDope para considerar un residuo/ventana como '
                              'epitopo. NO comparable en escala al umbral de BepiPred: cada motor '
                              'conserva el suyo.')
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
        self._insertFunctionStep(self.epidopeStep)
        self._insertFunctionStep(self.createOutputStep)

    # ---------------------------------- Steps -----------------------------------

    def writeInputFasta(self):
        faFile = self._getExtraPath('inputSequence.fa')
        self.inputSequence.get().exportToFile(faFile)
        return os.path.abspath(faFile)

    def epidopeStep(self):
        epidope_bin = bcellPlugin.getEpidopeBin()

        fasta_file = self.writeInputFasta()
        out_dir = os.path.abspath(self._getExtraPath('epidope_raw'))
        os.makedirs(out_dir, exist_ok=True)

        args = f'-i {fasta_file} -o {out_dir} -t {self.threshold.get()}'
        self.runJob(epidope_bin, args)

    def createOutputStep(self):
        out_dir = self._getExtraPath('epidope_raw')
        raw_df = self._loadRawScores(out_dir)

        epitopes_df = extract_epitope_regions(
            raw_df,
            accession_col=ACCESSION_COLUMN,
            score_col=SCORE_COLUMN,
            residue_col_candidates=(RESIDUE_COLUMN,),
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
                               description='EpiDope epitope')
            seqROI = SequenceROI(sequence=inputSeq, seqROI=roiSeq, roiIdx=row.start, roiIdx2=row.end)
            seqROI._meanScore = Float(row.mean_score)
            outROIs.append(seqROI)

        if len(outROIs) > 0:
            self._defineOutputs(outputROIs=outROIs)
            self._defineSourceRelation(self.inputSequence, outROIs)

    # ---------------------------------- Utils -----------------------------------

    @staticmethod
    def _loadRawScores(out_dir: str) -> pd.DataFrame:
        """Concatena los CSV por-accession que escribe EpiDope en '<out_dir>/epidope/'."""
        per_gene_dir = Path(out_dir) / 'epidope'
        csv_files = sorted(per_gene_dir.glob('*.csv'))
        if not csv_files:
            found = sorted(p.name for p in Path(out_dir).rglob('*') if p.is_file())
            raise EpidopeExecutionError(
                f"No se encontro ningun CSV de scores por-accession en '{per_gene_dir}'. "
                f"Archivos generados por EpiDope: {found or '<ninguno>'}."
            )

        frames = []
        for csv_path in csv_files:
            df = pd.read_csv(csv_path, sep='\t')
            missing = set(_RAW_COLUMNS) - set(df.columns)
            if missing:
                raise EpidopeExecutionError(
                    f"El CSV de salida '{csv_path}' no contiene las columnas confirmadas "
                    f"{sorted(missing)}. Columnas encontradas: {list(df.columns)}."
                )
            df = df.sort_values('position').reset_index(drop=True)
            df.insert(0, ACCESSION_COLUMN, csv_path.stem)
            frames.append(df)

        combined = pd.concat(frames, ignore_index=True)
        combined = combined.rename(columns={'aminoacid': RESIDUE_COLUMN, 'score': SCORE_COLUMN})
        return combined[[ACCESSION_COLUMN, RESIDUE_COLUMN, SCORE_COLUMN]]

    # ---------------------------------- Validation -------------------------------

    def _validate(self):
        return bcellPlugin.validateEpidopeInstallation()

    def _summary(self):
        summary = []
        if self.isFinished():
            outROIs = getattr(self, 'outputROIs', None)
            n = len(outROIs) if outROIs is not None else 0
            summary.append(f'{n} epitope region(s) found above threshold {self.threshold.get()}.')
        return summary
