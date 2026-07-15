"""Plugin Scipion-chem para el pipeline B-Cell-Epitope-Prediction.

BepiPred-3.0 no se instala automaticamente (licencia academica DTU Health
Tech, no redistribuible): el usuario lo descarga manualmente y configura
BEPIPRED_HOME / BEPIPRED_PYTHON_BIN en ``scipion.conf``, exactamente el mismo
patron ya validado en B-Cell-Epitope-Prediction/src/config/settings.py (un
interprete Python dedicado en vez de activacion de un entorno conda inline:
BepiPred-3.0 fija versiones antiguas de dependencias -torch==1.12.0,
numpy==1.20.2- incompatibles con el entorno principal de Scipion).
"""

import os

from pwchem import Plugin as pwchemPlugin

from .constants import (
    BEPIPRED_DIC, NOINSTALL_WARNING,
    EPIDOPE_DIC, EPIDOPE_NOINSTALL_WARNING,
)


class Plugin(pwchemPlugin):

    @classmethod
    def _defineVariables(cls):
        cls._defineVar(BEPIPRED_DIC['home'], '')
        cls._defineVar(BEPIPRED_DIC['python_bin'], '')
        cls._defineVar(EPIDOPE_DIC['home'], '')

    @classmethod
    def defineBinaries(cls, env):
        """No-op: ni BepiPred-3.0 ni EpiDope se instalan automaticamente (ver
        validateBepipredInstallation/validateEpidopeInstallation)."""
        pass

    @classmethod
    def validateInstallation(cls):
        """Agregado de todos los requisitos del plugin (llamado por el gestor de
        plugins de Scipion). Los protocolos individuales validan solo el motor
        que necesitan via validateBepipredInstallation/validateEpidopeInstallation."""
        return cls.validateBepipredInstallation() + cls.validateEpidopeInstallation()

    @classmethod
    def validateBepipredInstallation(cls):
        """Comprueba que BEPIPRED_HOME/BEPIPRED_PYTHON_BIN esten configurados y
        que el CLI de BepiPred-3.0 exista, devolviendo una lista de errores
        accionables (lista vacia = instalacion correcta)."""
        errors = []

        bepipred_home = cls.getVar(BEPIPRED_DIC['home'])
        if not bepipred_home or not os.path.isdir(bepipred_home):
            errors.append(
                f"BEPIPRED_HOME no configurado o no existe: '{bepipred_home}'."
            )
        else:
            cli_script = os.path.join(bepipred_home, BEPIPRED_DIC['cli_script'])
            if not os.path.isfile(cli_script):
                errors.append(f"No se encontro '{cli_script}' dentro de BEPIPRED_HOME.")

        python_bin = cls.getVar(BEPIPRED_DIC['python_bin'])
        if not python_bin:
            errors.append("BEPIPRED_PYTHON_BIN no configurado.")

        if errors:
            errors.append(NOINSTALL_WARNING)
        return errors

    @classmethod
    def validateEpidopeInstallation(cls):
        """Comprueba que EPIDOPE_HOME este configurado y que el binario
        'epidope' exista dentro de ese prefijo de entorno conda."""
        errors = []

        epidope_bin = cls.getEpidopeBin()
        if not epidope_bin or not os.path.isfile(epidope_bin):
            errors.append(f"No se encontro el binario de EpiDope en '{epidope_bin}'.")

        if errors:
            errors.append(EPIDOPE_NOINSTALL_WARNING)
        return errors

    # ---------------------------------- Utils -----------------------------------

    @classmethod
    def getBepipredHome(cls):
        return cls.getVar(BEPIPRED_DIC['home'])

    @classmethod
    def getBepipredPythonBin(cls):
        return cls.getVar(BEPIPRED_DIC['python_bin'])

    @classmethod
    def getBepipredCliScript(cls):
        return os.path.join(cls.getBepipredHome(), BEPIPRED_DIC['cli_script'])

    @classmethod
    def getEpidopeHome(cls):
        return cls.getVar(EPIDOPE_DIC['home'])

    @classmethod
    def getEpidopeBin(cls):
        home = cls.getEpidopeHome()
        if not home:
            return None
        return os.path.join(home, 'bin', EPIDOPE_DIC['binary'])
