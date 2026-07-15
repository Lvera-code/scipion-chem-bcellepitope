"""Protocolo Scipion-chem: Fase 3 (union logica anotada) del pipeline
B-Cell-Epitope-Prediction.

Toma los SetOfSequenceROIs producidos por ProtBCellEpitopeBepiPredPredict y/o
ProtBCellEpitopeEpiDopePredict (sobre la MISMA secuencia de entrada) y los une
(no interseca), fusionando regiones solapadas y etiquetando cada region final
como 'Consenso' (ambos motores), 'BepiPred' o 'EpiDope' (un solo motor). Ver
``bcellepitope.utils.consensus`` para el algoritmo (portado de
``B-Cell-Epitope-Prediction/src/engines/consensus.py``).
"""

from pwchem.objects import Sequence, SequenceROI, SetOfSequenceROIs
from pwem.protocols import EMProtocol
from pyworkflow.object import Float, String
from pyworkflow.protocol import params

from ..utils.consensus import MIN_FINAL_PEPTIDE_LENGTH, merge_annotated_intervals


class ProtBCellEpitopeConsensus(EMProtocol):
    """Union logica anotada de regiones de epitopo entre BepiPred-3.0 y EpiDope."""

    _label = 'bepipred + epidope consensus'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSequence', params.PointerParam, pointerClass='Sequence',
                       label='Input protein sequence: ',
                       help='La misma secuencia usada como entrada de los protocolos BepiPred-3.0 '
                            'y/o EpiDope cuyos ROI se van a unir.')
        form.addParam('inputROIsBepiPred', params.PointerParam, pointerClass='SetOfSequenceROIs',
                       allowsNull=True, label='BepiPred-3.0 epitope ROIs: ',
                       help='SetOfSequenceROIs producido por el protocolo BepiPred-3.0 de este '
                            'plugin. Opcional si solo se dispone de resultados de EpiDope.')
        form.addParam('inputROIsEpiDope', params.PointerParam, pointerClass='SetOfSequenceROIs',
                       allowsNull=True, label='EpiDope epitope ROIs: ',
                       help='SetOfSequenceROIs producido por el protocolo EpiDope de este plugin. '
                            'Opcional si solo se dispone de resultados de BepiPred-3.0.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.consensusStep)

    # ---------------------------------- Steps -----------------------------------

    def consensusStep(self):
        bpSet = self.inputROIsBepiPred.get()
        edSet = self.inputROIsEpiDope.get()

        bpIntervals = [(roi.getROIIdx(), roi.getROIIdx2(), roi._meanScore.get()) for roi in bpSet] if bpSet is not None else []
        edIntervals = [(roi.getROIIdx(), roi.getROIIdx2(), roi._meanScore.get()) for roi in edSet] if edSet is not None else []

        inputSeq = self.inputSequence.get()
        fullSequence = inputSeq.getSequence()

        records = merge_annotated_intervals(
            bpIntervals, edIntervals, fullSequence, min_length=MIN_FINAL_PEPTIDE_LENGTH
        )

        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for rec in records:
            roiId = f"ROI_{rec['start']}-{rec['end']}"
            roiSeq = Sequence(sequence=rec['sequence'], name=roiId, id=roiId,
                               description=f"{rec['origen']} epitope")
            seqROI = SequenceROI(sequence=inputSeq, seqROI=roiSeq, roiIdx=rec['start'], roiIdx2=rec['end'])
            seqROI._origen = String(rec['origen'])
            seqROI._bepipredRegion = String(rec['bepipred_region'])
            seqROI._epidopeRegion = String(rec['epidope_region'])
            if rec['bepipred_score'] is not None:
                seqROI._bepipredScore = Float(rec['bepipred_score'])
            if rec['epidope_score'] is not None:
                seqROI._epidopeScore = Float(rec['epidope_score'])
            outROIs.append(seqROI)

        if len(outROIs) > 0:
            self._defineOutputs(outputROIs=outROIs)
            if bpSet is not None:
                self._defineSourceRelation(self.inputROIsBepiPred, outROIs)
            if edSet is not None:
                self._defineSourceRelation(self.inputROIsEpiDope, outROIs)

    # ---------------------------------- Validation -------------------------------

    def _validate(self):
        errors = []
        if self.inputROIsBepiPred.get() is None and self.inputROIsEpiDope.get() is None:
            errors.append('Debes proporcionar al menos un SetOfSequenceROIs (BepiPred-3.0 y/o EpiDope).')
        return errors

    def _summary(self):
        summary = []
        if self.isFinished():
            outROIs = getattr(self, 'outputROIs', None)
            if outROIs is None:
                summary.append(f'No se encontraron regiones de epitopo de al menos '
                                f'{MIN_FINAL_PEPTIDE_LENGTH} aa tras la union.')
            else:
                nConsenso = sum(1 for roi in outROIs if roi._origen.get() == 'Consenso')
                nBepipred = sum(1 for roi in outROIs if roi._origen.get() == 'BepiPred')
                nEpidope = sum(1 for roi in outROIs if roi._origen.get() == 'EpiDope')
                summary.append(f'{len(outROIs)} region(es) tras la union: {nConsenso} Consenso, '
                                f'{nBepipred} solo BepiPred, {nEpidope} solo EpiDope.')
        return summary
