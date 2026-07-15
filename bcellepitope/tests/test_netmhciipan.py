from pyworkflow.tests import setupTestProject, BaseTest

from pwem.protocols import ProtImportSequence

from ..protocols import (
    ProtBCellEpitopeBepiPredPredict,
    ProtBCellEpitopeBlastFilter,
    ProtBCellEpitopeConsensus,
    ProtBCellEpitopeEpiDopePredict,
    ProtBCellEpitopeNetMHCIIpan,
)


class TestBCellEpitopeNetMHCIIpan(BaseTest):
    NAME = 'GP120_P03377'
    DESCRIPTION = 'ENV_HV1BR (GP120), UniProt P03377'
    AMINOACIDSSEQ = (
        'MRVKEKYQHLWRWGWKWGTMLLGILMICSATEKLWVTVYYGVPVWKEATTTLFCASDAKAYDTEVHNVW'
        'ATHACVPTDPNPQEVVLVNVTENFNMWKNDMVEQMHEDIISLWDQSLKPCVKLTPLCVSLKCTDLGNAT'
        'NTNSSNTNSSSGEMMMEKGEIKNCSFNISTSIRGKVQKEYAFFYKLDIIPIDNDTTSYTLTSCNTSVIT'
        'QACPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGPCTNVSTVQCTHGIRPVVSTQLLLNGSLAEEEVV'
        'IRSANFTDNAKTIIVQLNQSVEINCTRPNNNTRKSIRIQRGPGRAFVTIGKIGNMRQAHCNISRAKWNA'
        'TLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHSFNCGGEFFYCNSTQLFNSTWFNSTWSTEGSNNTE'
        'GSDTITLPCRIKQFINMWQEVGKAMYAPPISGQIRCSSNITGLLLTRDGGNNNNGSEIFRPGGGDMRDN'
        'WRSELYKYKVVKIEPLGVAPTKAKRRVVQREKRAVGIGALFLGFLGAAGSTMGARSMTLTVQARQLLSG'
        'IVQQQNNLLRAIEAQQHLLQLTVWGIKQLQARILAVERYLKDQQLLGIWGCSGKLICTTAVPWNASWSN'
        'KSLEQIWNNMTWMEWDREINNYTSLIHSLIEESQNQQEKNEQELLELDKWASLWNWFNITNWLWYIKI'
        'FIMIVGGLVGLRIVFAVLSIVNRVRQGYSPLSFQTHLPTPRGPDRPEGIEEEGGERDRDRSIRLVNGSL'
        'ALIWDDLRSLCLFSYHRLRDLLLIVTRIVELLGRRGWEALKYWWNLLQYWSQELKNSAVSLLNATAIAV'
        'AEGTDRVIEVVQGACRAIRHIPRRIRQGLERILL'
    )

    # Referencia real (regenerada corriendo pipeline.py --input fasta_inputs/GP120.fasta):
    # fasta_outputs/candidatos_finales.csv del repo B-Cell-Epitope-Prediction.
    # 17 candidatos tras traceback + deduplicacion de ventanas de modo proteina.
    EXPECTED = sorted([
        (88, 123, 'WKNDMVEQM', 3, 27),
        (309, 323, 'IRIQRGPGR', 7, 27),
        (316, 330, 'RAFVTIGKI', 4, 27),
        (317, 331, 'FVTIGKIGN', 4, 27),
        (318, 332, 'FVTIGKIGN', 6, 27),
        (319, 333, 'FVTIGKIGN', 5, 27),
        (338, 352, 'KWNATLKQI', 3, 27),
        (339, 353, 'KWNATLKQI', 5, 27),
        (340, 354, 'WNATLKQIA', 4, 27),
        (349, 363, 'KLREQFGNN', 3, 27),
        (350, 364, 'LREQFGNNK', 3, 27),
        (639, 653, 'NNYTSLIHS', 3, 27),
        (640, 654, 'YTSLIHSLI', 4, 27),
        (642, 656, 'LIHSLIEES', 5, 27),
        (643, 657, 'LIHSLIEES', 8, 27),
        (644, 658, 'IHSLIEESQ', 7, 27),
        (656, 670, 'KNEQELLEL', 3, 27),
    ])

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

    def _runProt(self, protClass, **extraPointers):
        prot = self.newProtocol(protClass)
        prot.inputSequence.set(self.protImportSeq)
        prot.inputSequence.setExtended('outputSequence')
        for paramName, (sourceProt, extended) in extraPointers.items():
            getattr(prot, paramName).set(sourceProt)
            getattr(prot, paramName).setExtended(extended)
        self.launchProtocol(prot, wait=True)
        return prot

    def test(self):
        protBepiPred = self._runProt(ProtBCellEpitopeBepiPredPredict)
        protEpiDope = self._runProt(ProtBCellEpitopeEpiDopePredict)

        protConsensus = self.newProtocol(ProtBCellEpitopeConsensus)
        protConsensus.inputSequence.set(self.protImportSeq)
        protConsensus.inputSequence.setExtended('outputSequence')
        protConsensus.inputROIsBepiPred.set(protBepiPred)
        protConsensus.inputROIsBepiPred.setExtended('outputROIs')
        protConsensus.inputROIsEpiDope.set(protEpiDope)
        protConsensus.inputROIsEpiDope.setExtended('outputROIs')
        self.launchProtocol(protConsensus, wait=True)

        protBlast = self._runProt(ProtBCellEpitopeBlastFilter, inputROIs=(protConsensus, 'outputROIs'))
        protNetMHCIIpan = self._runProt(ProtBCellEpitopeNetMHCIIpan, inputROIs=(protBlast, 'outputROIs'))

        outROIs = getattr(protNetMHCIIpan, 'outputROIs', None)
        self.assertIsNotNone(outROIs)
        got = sorted(
            (roi.getROIIdx(), roi.getROIIdx2(), roi._core9aa.get(),
             roi._nAlelosPromiscuos.get(), roi._nAlelosEvaluados.get())
            for roi in outROIs
        )
        self.assertEqual(got, self.EXPECTED)
