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
