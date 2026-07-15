"""Parseo de la salida de NetMHCIIpan-4.3 y traceback a la region padre.

Portado de B-Cell-Epitope-Prediction/src/engines/netmhciipan_engine.py
(``_parse_xls``, ``build_traceback_report``,
``_deduplicate_protein_mode_windows``), simplificado: el CLI original agrupa
por ``accession`` porque procesa FASTA multi-registro; en Scipion cada
protocolo corre sobre un unico ``Sequence`` de entrada (mismo criterio ya
aplicado en ``consensus.py``), asi que no hace falta esa columna.
"""

from typing import List, Optional

import numpy as np
import pandas as pd

# Panel de referencia de 27 alelos HLA-DR/DQ/DP mas representativos usado por
# el IEDB para estimar cobertura poblacional amplia en el diseno de epitopos
# T-helper (MHC-II). NUNCA se le agregan espacios entre comas: NetMHCIIpan lo
# pasa tal cual a su parser de '-a' y un espacio rompe el parseo del alelo
# siguiente.
IEDB_REFERENCE_PANEL = (
    "DRB1_0101,DRB1_0301,DRB1_0401,DRB1_0405,DRB1_0701,DRB1_0802,DRB1_0901,"
    "DRB1_1101,DRB1_1201,DRB1_1302,DRB1_1501,DRB3_0101,DRB3_0202,DRB4_0101,DRB5_0101,"
    "HLA-DQA10501-DQB10201,HLA-DQA10501-DQB10301,HLA-DQA10301-DQB10302,"
    "HLA-DQA10401-DQB10402,HLA-DQA10101-DQB10501,HLA-DQA10102-DQB10602,"
    "HLA-DPA10201-DPB10101,HLA-DPA10103-DPB10201,HLA-DPA10103-DPB10401,"
    "HLA-DPA10301-DPB10402,HLA-DPA10201-DPB10501,HLA-DPA10201-DPB11401"
)

# Footprint minimo del core de union a MHC-II: NetMHCIIpan descarta (o
# calcula sobre un core mas corto que el peptido, degradando la prediccion)
# peptidos mas cortos que esto. Deliberadamente NO configurable (piso
# biologico fijo, mismo criterio que MIN_FINAL_PEPTIDE_LENGTH en consensus.py).
MIN_PEPTIDE_LENGTH = 9

# Longitud maxima segura para el modo peptido exacto ('-p'): el binario
# NetMHCIIpan-4.3 (Linux_x86_64) revienta con "buffer overflow" (SIGABRT)
# para entradas mas largas, ver ADR completo en netmhciipan_engine.py. NO
# configurable: es un limite de seguridad del binario, no un parametro de
# usuario.
MAX_PEPTIDE_MODE_LENGTH = 40

_OUTPUT_COLUMNS = [
    "sequence", "core_9aa", "n_alelos_evaluados", "n_alelos_promiscuos", "min_rank_el", "veredicto",
]

_TRACEBACK_COLUMNS = [
    "sequence_f5", "core_9aa", "start", "end", "origen",
    "n_alelos_promiscuos", "n_alelos_evaluados", "min_rank_el",
    "bepipred_score", "epidope_score",
]


class NetMHCIIpanParseError(Exception):
    """El .xls de salida de NetMHCIIpan no tiene el formato esperado."""


def parse_xls(xls_path: str, n_alleles: int, rank_weak: float, min_promiscuous_alleles: int) -> pd.DataFrame:
    """Parsea el .xls de NetMHCIIpan y evalua la promiscuidad de cada peptido.

    Alelos invertidos (columna ``Inverted``) se excluyen POR COMPLETO antes de
    cualquier calculo: ni cuentan para la promiscuidad ni pueden ser el alelo
    "ganador" que determina ``core_9aa``/``min_rank_el`` (ver ADR completo en
    ``netmhciipan_engine.py::_parse_xls`` del pipeline original). Esto
    garantiza que ``core_9aa`` siempre es una subcadena literal del peptido de
    entrada.

    Returns:
        DataFrame con columnas ``sequence``, ``core_9aa``,
        ``n_alelos_evaluados``, ``n_alelos_promiscuos``, ``min_rank_el`` y
        ``veredicto`` ('Candidato Valido' / 'Rechazado').
    """
    try:
        raw = pd.read_csv(xls_path, sep="\t", skiprows=2)
    except Exception as exc:
        raise NetMHCIIpanParseError(f"No se pudo parsear la salida de NetMHCIIpan en '{xls_path}': {exc}") from exc

    rank_cols = [c for c in raw.columns if c == "Rank_EL" or c.startswith("Rank_EL.")]
    core_cols = [c for c in raw.columns if c == "Core" or c.startswith("Core.")]
    inverted_cols = [c for c in raw.columns if c == "Inverted" or c.startswith("Inverted.")]
    if (
        len(rank_cols) != n_alleles
        or len(core_cols) != n_alleles
        or len(inverted_cols) != n_alleles
        or "Peptide" not in raw.columns
    ):
        raise NetMHCIIpanParseError(
            f"El formato de salida .xls de NetMHCIIpan no coincide con lo esperado: "
            f"se encontraron {len(rank_cols)} columna(s) 'Rank_EL', {len(core_cols)} "
            f"columna(s) 'Core' y {len(inverted_cols)} columna(s) 'Inverted' para "
            f"{n_alleles} alelo(s) evaluado(s). Columnas encontradas: {list(raw.columns)}."
        )

    rank_matrix = raw[rank_cols].to_numpy()
    core_matrix = raw[core_cols].to_numpy()
    is_inverted = raw[inverted_cols].to_numpy().astype(bool)
    row_idx = np.arange(len(raw))

    rank_matrix_normal = np.where(is_inverted, np.inf, rank_matrix)
    best_allele_idx = rank_matrix_normal.argmin(axis=1)
    best_core = core_matrix[row_idx, best_allele_idx]

    is_binder_normal = (rank_matrix_normal <= rank_weak)
    n_alelos_promiscuos = is_binder_normal.sum(axis=1)

    result = pd.DataFrame(
        {
            "sequence": raw["Peptide"],
            "core_9aa": best_core,
            "n_alelos_evaluados": n_alleles,
            "n_alelos_promiscuos": n_alelos_promiscuos,
            "min_rank_el": rank_matrix_normal.min(axis=1),
        }
    )
    result["veredicto"] = result["n_alelos_promiscuos"].apply(
        lambda n: "Candidato Valido" if n >= min_promiscuous_alleles else "Rechazado"
    )
    return result[_OUTPUT_COLUMNS]


def build_traceback_report(report_df: pd.DataFrame, parent_records: List[dict]) -> pd.DataFrame:
    """Cruza los 'Candidato Valido' con su ROI padre (Fase 3/4) y deduplica ventanas.

    Args:
        report_df: Salida de :func:`parse_xls` (incluye rechazados; aqui se
            filtra a ``veredicto == 'Candidato Valido'``).
        parent_records: Lista de dicts, uno por ROI evaluado, con
            ``start``, ``sequence``, ``origen``, ``bepipred_score``,
            ``epidope_score`` (estos dos ultimos ``None`` si el motor
            correspondiente no contribuyo a esa region).

    Returns:
        DataFrame con columnas ``_TRACEBACK_COLUMNS``, ventanas redundantes
        del modo proteina ya colapsadas (mismo core_9aa + misma
        n_alelos_promiscuos -> se conserva solo la de menor min_rank_el).
    """
    if report_df.empty or not parent_records:
        return pd.DataFrame(columns=_TRACEBACK_COLUMNS)

    valid_df = report_df[report_df["veredicto"] == "Candidato Valido"]
    if valid_df.empty:
        return pd.DataFrame(columns=_TRACEBACK_COLUMNS)

    records = []
    for candidate in valid_df.itertuples(index=False):
        for parent in parent_records:
            if candidate.sequence not in parent["sequence"]:
                continue
            offset = parent["sequence"].find(candidate.sequence)
            start_real = parent["start"] + offset
            end_real = start_real + len(candidate.sequence) - 1
            records.append(
                {
                    "sequence_f5": candidate.sequence,
                    "core_9aa": candidate.core_9aa,
                    "start": start_real,
                    "end": end_real,
                    "origen": parent["origen"],
                    "n_alelos_promiscuos": candidate.n_alelos_promiscuos,
                    "n_alelos_evaluados": candidate.n_alelos_evaluados,
                    "min_rank_el": candidate.min_rank_el,
                    "bepipred_score": parent["bepipred_score"],
                    "epidope_score": parent["epidope_score"],
                }
            )

    traceback_df = pd.DataFrame.from_records(records, columns=_TRACEBACK_COLUMNS)
    return _deduplicate_protein_mode_windows(traceback_df)


def _deduplicate_protein_mode_windows(traceback_df: pd.DataFrame) -> pd.DataFrame:
    """Colapsa ventanas redundantes del modo proteina (mismo core_9aa + misma
    n_alelos_promiscuos -> se conserva solo la de menor min_rank_el)."""
    if traceback_df.empty:
        return traceback_df

    best_idx = traceback_df.groupby(["core_9aa", "n_alelos_promiscuos"], sort=False)["min_rank_el"].idxmin()
    return traceback_df.loc[best_idx].sort_index().reset_index(drop=True)
