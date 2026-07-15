"""Ventana deslizante tolerante a gaps para mapear regiones de epitopo B-cell.

Portado tal cual (mismo autor, misma logica validada) desde
B-Cell-Epitope-Prediction/src/engines/epitope_mapping.py, para no depender en
tiempo de ejecucion de ese repo (los plugins de Scipion se instalan de forma
independiente, no como submodulo del pipeline original).
"""

import warnings
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd


def resolve_residue_column(df: pd.DataFrame, candidates: Sequence[str]) -> Optional[str]:
    for candidate in candidates:
        if candidate in df.columns:
            return candidate
    warnings.warn(
        f"No se encontro ninguna columna de residuo entre {candidates}. "
        f"Columnas disponibles: {list(df.columns)}. Las regiones de epitopo se "
        "reportaran sin secuencia de aminoacidos."
    )
    return None


def find_valid_windows(
    scores: List[float], threshold: float, window_size: int, max_gap_residues: int
) -> List[Tuple[int, int]]:
    """Desliza una ventana de ``window_size`` (paso=1) y devuelve los rangos validos.

    Una ventana ``[i, i + window_size - 1]`` (0-indexada, inclusive) es valida
    si, a la vez: (a) a lo sumo ``max_gap_residues`` de sus residuos tienen un
    score individual por debajo de ``threshold``, y (b) el score medio de la
    ventana completa es ``>= threshold``.
    """
    n = len(scores)
    valid_windows = []
    for i in range(0, n - window_size + 1):
        window = scores[i : i + window_size]
        below_count = sum(1 for score in window if score < threshold)
        if below_count <= max_gap_residues and (sum(window) / window_size) >= threshold:
            valid_windows.append((i, i + window_size - 1))
    return valid_windows


def merge_overlapping_windows(windows: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Fusiona ventanas validas solapadas o adyacentes en regiones continuas."""
    if not windows:
        return []

    merged = [windows[0]]
    for start, end in windows[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + 1:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def extract_epitope_regions(
    raw_scores_df: pd.DataFrame,
    accession_col: str,
    score_col: str,
    residue_col_candidates: Sequence[str],
    threshold: float,
    min_length: int,
    window_size: int,
    max_gap_residues: int,
) -> pd.DataFrame:
    """Mapea regiones de epitopo con una ventana deslizante tolerante a gaps.

    Returns:
        DataFrame con columnas ``accession``, ``start``, ``end``, ``length``,
        ``mean_score``, ``max_score`` y ``sequence``.
    """
    missing = {accession_col, score_col} - set(raw_scores_df.columns)
    if missing:
        raise ValueError(
            f"El DataFrame de entrada no contiene las columnas requeridas {sorted(missing)}. "
            f"Columnas encontradas: {list(raw_scores_df.columns)}."
        )

    residue_col = resolve_residue_column(raw_scores_df, residue_col_candidates)
    records = []

    for accession, group in raw_scores_df.groupby(accession_col, sort=False):
        group = group.reset_index(drop=True)
        scores = group[score_col].tolist()

        valid_windows = find_valid_windows(scores, threshold, window_size, max_gap_residues)
        merged_regions = merge_overlapping_windows(valid_windows)

        for start, end in merged_regions:
            length = end - start + 1
            if length < min_length:
                continue

            block = group.iloc[start : end + 1]
            sequence = "".join(block[residue_col].astype(str)) if residue_col else ""
            records.append(
                {
                    "accession": accession,
                    "start": start + 1,
                    "end": end + 1,
                    "length": length,
                    "mean_score": float(block[score_col].mean()),
                    "max_score": float(block[score_col].max()),
                    "sequence": sequence,
                }
            )

    return pd.DataFrame.from_records(
        records,
        columns=["accession", "start", "end", "length", "mean_score", "max_score", "sequence"],
    )
