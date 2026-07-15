from pyworkflow.tests import setupTestProject, BaseTest

from pwem.protocols import ProtImportSequence

from ..protocols import ProtBCellEpitopeEpiDopePredict


class TestBCellEpitopeEpiDopePrediction(BaseTest):
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

    def _runEpiDopePrediction(self):
        protEpiDope = self.newProtocol(ProtBCellEpitopeEpiDopePredict)
        protEpiDope.inputSequence.set(self.protImportSeq)
        protEpiDope.inputSequence.setExtended('outputSequence')
        return protEpiDope

    def test(self):
        protEpiDope = self._runEpiDopePrediction()
        # Con la lisozima y el umbral por defecto de EpiDope (0.818), el
        # pipeline original (fasta_outputs/lisozima_bepipred_epidope_epitopes.csv)
        # no encuentra ninguna region de epitopo: el protocolo no define
        # 'outputROIs' en absoluto si el DataFrame de epitopos queda vacio (ver
        # createOutputStep), asi que esperamos sincronicamente a que el
        # protocolo termine (en vez de a un output concreto) y solo
        # comprobamos que no haya fallado.
        self.launchProtocol(protEpiDope, wait=True)
