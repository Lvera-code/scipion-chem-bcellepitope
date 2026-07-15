"""Jerarquia de excepciones especificas del plugin, portada de
B-Cell-Epitope-Prediction/src/utils/exceptions.py (mismo autor, mismo criterio:
nunca dejar escapar un FileNotFoundError/CalledProcessError crudo hacia la GUI
de Scipion sin un mensaje accionable).
"""


class BCellEpitopeExecutionError(Exception):
    """Clase base para todos los errores controlados de este plugin."""


class BepiPredExecutionError(BCellEpitopeExecutionError):
    """Fallo al ejecutar BepiPred-3.0 localmente (instalacion ausente, subprocess
    fallido, timeout o salida en un formato inesperado)."""


class EpidopeExecutionError(BCellEpitopeExecutionError):
    """Fallo al ejecutar EpiDope localmente (instalacion ausente, subprocess
    fallido, timeout o salida en un formato inesperado)."""


class BlastExecutionError(BCellEpitopeExecutionError):
    """Fallo al ejecutar el filtro de tolerancia inmunologica (BLASTp local):
    binario ausente, base de datos no indexada, o subprocess fallido/timeout."""


class NetMHCIIpanExecutionError(BCellEpitopeExecutionError):
    """Fallo al ejecutar NetMHCIIpan-4.3 localmente: instalacion ausente,
    subprocess fallido/timeout, o el .xls de salida no se genero/no tiene el
    formato esperado."""
