from pyworkflow.tests import setupTestProject, BaseTest

from pwem.protocols import ProtImportSequence

from ..protocols import ProtBCellEpitopeBepiPredPredict


class TestBCellEpitopeBepiPredPrediction(BaseTest):
    NAME = 'LISOZIMA_P00698'
    DESCRIPTION = 'Lisozima C Gallus gallus (UniProt P00698), control positivo de B-Cell-Epitope-Prediction'
    AMINOACIDSSEQ = (
        'MRSLLILVLCFLPLAALGKVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGST'
        'DYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDV'
        'QAWIRGCRL'
    )

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

    def _runBepiPredPrediction(self):
        protBepiPred = self.newProtocol(ProtBCellEpitopeBepiPredPredict)
        protBepiPred.inputSequence.set(self.protImportSeq)
        protBepiPred.inputSequence.setExtended('outputSequence')
        self.proj.launchProtocol(protBepiPred, wait=False)
        return protBepiPred

    def test(self):
        protBepiPred = self._runBepiPredPrediction()
        self._waitOutput(protBepiPred, 'outputROIs', sleepTime=10)
        outROIs = getattr(protBepiPred, 'outputROIs', None)
        self.assertIsNotNone(outROIs)
        self.assertGreater(len(outROIs), 0)
