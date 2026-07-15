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

# EpiDope: codigo abierto (licencia MIT, github.com/rnajena/EpiDope), instalable
# via conda sin solicitud academica, pero con un stack de dependencias muy
# antiguo (Python 3.6, TensorFlow 1.13, ELMo/AllenNLP) que exige un entorno
# conda dedicado creado EXACTAMENTE con el 'epidope.yml' oficial del repo. El
# binario 'epidope' que instala conda en '<prefix>/bin/epidope' es un shim
# autocontenido (shebang apuntando al python del propio entorno): se invoca
# directamente, SIN pasar por 'conda run' (ver ADR en
# B-Cell-Epitope-Prediction/src/engines/epidope_engine.py::_build_command:
# 'conda run --no-capture-output' demostro ser poco fiable con stdout/stderr
# capturados por pipe, fallos espurios intermitentes).
EPIDOPE_DIC = {
    'name': 'EpiDope',
    'home': 'EPIDOPE_HOME',
    'binary': 'epidope',
}

EPIDOPE_DOWNLOAD_URL = 'https://github.com/rnajena/EpiDope'

EPIDOPE_NOINSTALL_WARNING = (
    'EpiDope installation could not be found. A diferencia de BepiPred-3.0, '
    f'EpiDope es codigo abierto y no requiere solicitud academica: {EPIDOPE_DOWNLOAD_URL} . '
    'Instalalo en un entorno conda dedicado con el epidope.yml oficial del repo '
    '(NO instales los paquetes a mano) y configura EPIDOPE_HOME apuntando al '
    'prefijo de ese entorno. Ver README.rst.'
)
