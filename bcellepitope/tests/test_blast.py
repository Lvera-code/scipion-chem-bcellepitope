from pyworkflow.tests import setupTestProject, BaseTest

from pwem.protocols import ProtImportSequence

from ..protocols import (
    ProtBCellEpitopeBepiPredPredict,
    ProtBCellEpitopeBlastFilter,
    ProtBCellEpitopeConsensus,
)


class TestBCellEpitopeBlastFilter(BaseTest):
    NAME = 'LISOZIMA_P00698'
    DESCRIPTION = 'Lisozima C Gallus gallus (UniProt P00698), control positivo de B-Cell-Epitope-Prediction'
    AMINOACIDSSEQ = (
        'MRSLLILVLCFLPLAALGKVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGST'
        'DYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDV'
        'QAWIRGCRL'
    )

    # Referencia real: fasta_outputs/lisozima_bepipred_blast_report.csv del
    # repo B-Cell-Epitope-Prediction. Los 3 peptidos de BepiPred son <= 30 aa
    # (tramo 'corto'), por lo que usan blastp-short/evalue=50.0, y los 3 dan
    # 100% de identidad contra el proteoma humano -> 'Autoinmunidad'.
    EXPECTED = [
        (59, 72, 'blastp-short', 50.0, 100.0, 'Autoinmunidad'),
        (115, 127, 'blastp-short', 50.0, 100.0, 'Autoinmunidad'),
        (129, 140, 'blastp-short', 50.0, 100.0, 'Autoinmunidad'),
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

    def _runConsensus(self, protBepiPred):
        prot = self.newProtocol(ProtBCellEpitopeConsensus)
        prot.inputSequence.set(self.protImportSeq)
        prot.inputSequence.setExtended('outputSequence')
        prot.inputROIsBepiPred.set(protBepiPred)
        prot.inputROIsBepiPred.setExtended('outputROIs')
        self.launchProtocol(prot, wait=True)
        return prot

    def test(self):
        protBepiPred = self._runBepiPred()
        protConsensus = self._runConsensus(protBepiPred)

        protBlast = self.newProtocol(ProtBCellEpitopeBlastFilter)
        protBlast.inputSequence.set(self.protImportSeq)
        protBlast.inputSequence.setExtended('outputSequence')
        protBlast.inputROIs.set(protConsensus)
        protBlast.inputROIs.setExtended('outputROIs')
        self.launchProtocol(protBlast, wait=True)

        outROIs = getattr(protBlast, 'outputROIs', None)
        self.assertIsNotNone(outROIs)
        got = sorted(
            (roi.getROIIdx(), roi.getROIIdx2(), roi._blastTask.get(),
             roi._blastEvalue.get(), roi._maxPident.get(), roi._status.get())
            for roi in outROIs
        )
        self.assertEqual(got, self.EXPECTED)
