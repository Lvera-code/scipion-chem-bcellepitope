"""Seleccion dinamica de tarea/E-value de BLASTp por longitud de peptido.

Portado de B-Cell-Epitope-Prediction/src/engines/blast_engine.py
(``_select_task``/``_select_evalue``): BLAST recomienda un modo de busqueda y
una sensibilidad estadistica distintos segun la longitud de la secuencia
consultada, y la estadistica por defecto penaliza a los peptidos cortos (un
match identico de 9-25 aa contra un proteoma completo puede descartarse como
"no significativo" con el E-value estandar de blastp), arruinando el filtro
de autoinmunidad justo donde mas importa.
"""

DEFAULT_SHORT_PEPTIDE_MAX_LEN = 30
DEFAULT_MEDIUM_PEPTIDE_MAX_LEN = 100
DEFAULT_EVALUE_SHORT = 50.0
DEFAULT_EVALUE_MEDIUM = 0.1
DEFAULT_EVALUE_LONG = 0.05


def select_task(sequence_length: int, short_max_len: int = DEFAULT_SHORT_PEPTIDE_MAX_LEN) -> str:
    """'blastp-short' si sequence_length <= short_max_len, si no 'blastp'."""
    return "blastp-short" if sequence_length <= short_max_len else "blastp"


def select_evalue(
    sequence_length: int,
    short_max_len: int = DEFAULT_SHORT_PEPTIDE_MAX_LEN,
    medium_max_len: int = DEFAULT_MEDIUM_PEPTIDE_MAX_LEN,
    evalue_short: float = DEFAULT_EVALUE_SHORT,
    evalue_medium: float = DEFAULT_EVALUE_MEDIUM,
    evalue_long: float = DEFAULT_EVALUE_LONG,
) -> float:
    """E-value segun el tramo de longitud del peptido (corto/medio/largo)."""
    if sequence_length <= short_max_len:
        return evalue_short
    if sequence_length <= medium_max_len:
        return evalue_medium
    return evalue_long
