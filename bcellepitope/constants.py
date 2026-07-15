# Common constants

DEFAULT_VERSION = '3.0'

# BepiPred-3.0: se ejecuta 100% en local via subprocess, sin instalacion
# automatica por parte de Scipion (licencia academica DTU Health Tech, no
# redistribuible). El usuario descarga el paquete manualmente y crea un
# entorno Python dedicado (ver README.rst), igual que en el pipeline
# B-Cell-Epitope-Prediction original (src/config/settings.py).
BEPIPRED_DIC = {
    'name': 'BepiPred',
    'version': '3.0',
    'home': 'BEPIPRED_HOME',
    'python_bin': 'BEPIPRED_PYTHON_BIN',
    'cli_script': 'bepipred3_CLI.py',
}

READ_URL = 'https://github.com/Lvera-code/scipion-chem-bcellepitope'

NOINSTALL_WARNING = (
    'BepiPred-3.0 installation could not be found. Please check the '
    'scipion-chem-bcellepitope README file to see how to download BepiPred-3.0 '
    f'and configure BEPIPRED_HOME / BEPIPRED_PYTHON_BIN. See {READ_URL}'
)
