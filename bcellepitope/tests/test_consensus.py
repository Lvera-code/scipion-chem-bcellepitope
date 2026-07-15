from pyworkflow.tests import setupTestProject, BaseTest

from pwem.protocols import ProtImportSequence

from ..protocols import (
    ProtBCellEpitopeBepiPredPredict,
    ProtBCellEpitopeConsensus,
    ProtBCellEpitopeEpiDopePredict,
)


class TestBCellEpitopeConsensus(BaseTest):
    NAME = 'LISOZIMA_P00698'
    DESCRIPTION = 'Lisozima C Gallus gallus (UniProt P00698), control positivo de B-Cell-Epitope-Prediction'
    AMINOACIDSSEQ = (
        'MRSLLILVLCFLPLAALGKVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGST'
        'DYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDV'
        'QAWIRGCRL'
    )

    # Referencia real: fasta_outputs/lisozima_bepipred_union_epitopes.csv del
    # repo B-Cell-Epitope-Prediction. EpiDope no encuentra epitopos para esta
    # secuencia con su umbral por defecto, asi que la union coincide con los
    # ROI de BepiPred sin fusiones (todos 'BepiPred', epidope_score ausente).
    EXPECTED = [
        (59, 72, 'QATNRNTDGSTDYG', 'BepiPred'),
        (115, 127, 'KIVSDGNGMNAWV', 'BepiPred'),
        (129, 140, 'WRNRCKGTDVQA', 'BepiPred'),
    ]

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        setupTestProject(cls)

        cls._runImportSeq()
        cls._waitOutput(cls.protImportSeq, 'outputSequence', sleepTime=5)

    @classmethod
    def _runImportSeq(cls):
        kwargs = {
            'inputSequenceName': cls.NAME,
            'inputSequenceDescription': cls.DESCRIPTION,
            'inputRawSequence': cls.AMINOACIDSSEQ,
        }
        cls.protImportSeq = cls.newProtocol(ProtImportSequence, **kwargs)
        cls.proj.launchProtocol(cls.protImportSeq, wait=False)

    def _runBepiPred(self):
        prot = self.newProtocol(ProtBCellEpitopeBepiPredPredict)
        prot.inputSequence.set(self.protImportSeq)
        prot.inputSequence.setExtended('outputSequence')
        self.launchProtocol(prot, wait=True)
        return prot

    def _runEpiDope(self):
        prot = self.newProtocol(ProtBCellEpitopeEpiDopePredict)
        prot.inputSequence.set(self.protImportSeq)
        prot.inputSequence.setExtended('outputSequence')
        self.launchProtocol(prot, wait=True)
        return prot

    def test(self):
        protBepiPred = self._runBepiPred()
        protEpiDope = self._runEpiDope()

        protConsensus = self.newProtocol(ProtBCellEpitopeConsensus)
        protConsensus.inputSequence.set(self.protImportSeq)
        protConsensus.inputSequence.setExtended('outputSequence')
        protConsensus.inputROIsBepiPred.set(protBepiPred)
        protConsensus.inputROIsBepiPred.setExtended('outputROIs')
        # EpiDope no genero outputROIs (0 epitopos): se deja sin fijar, igual
        # que se dejaria en la GUI si ese protocolo no produjo salida.
        self.launchProtocol(protConsensus, wait=True)

        outROIs = getattr(protConsensus, 'outputROIs', None)
        self.assertIsNotNone(outROIs)
        got = sorted(
            (roi.getROIIdx(), roi.getROIIdx2(), roi.getROISequence(), roi._origen.get())
            for roi in outROIs
        )
        self.assertEqual(got, self.EXPECTED)
