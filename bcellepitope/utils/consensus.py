"""Union logica anotada de regiones de epitopo entre BepiPred-3.0 y EpiDope.

Portado de B-Cell-Epitope-Prediction/src/engines/consensus.py::_merge_accession_intervals,
simplificado: el pipeline CLI original opera sobre FASTA multi-accession (de
ahi la normalizacion ``accession_id`` para cruzar cabeceras entre motores);
en Scipion cada protocolo de motor (BepiPred/EpiDope) corre sobre un unico
``Sequence`` de entrada, asi que ambos SetOfSequenceROIs de un mismo consenso
ya pertenecen a la misma secuencia por construccion y no hace falta agrupar
por accession.

Misma directiva de "Union Logica Anotada" que el original:

* Preservacion de datos: TODA region de BepiPred y/o EpiDope avanza al
  resultado final, no solo las que coinciden entre ambos motores.
* Fusion de solapamientos: regiones que comparten al menos un residuo se
  fusionan (start minimo, end maximo), incluida fusion transitiva
  (A-B-C encadenadas por solapes sucesivos se fusionan en una sola region
  aunque A y C por si solas no se solapen).
* Etiquetado: cada region final queda marcada como 'Consenso' (ambos
  motores), 'BepiPred' o 'EpiDope' (un solo motor).
"""

from typing import List, Optional, Tuple

MIN_FINAL_PEPTIDE_LENGTH = 9


def merge_annotated_intervals(
    bepipred_intervals: List[Tuple[int, int, float]],
    epidope_intervals: List[Tuple[int, int, float]],
    full_sequence: str,
    min_length: int = MIN_FINAL_PEPTIDE_LENGTH,
) -> List[dict]:
    """Fusiona por solapamiento los intervalos de BepiPred y EpiDope de una secuencia.

    Args:
        bepipred_intervals: Lista de ``(start, end, mean_score)`` (1-indexados,
            inclusive) de las regiones detectadas por BepiPred-3.0.
        epidope_intervals: Igual, para EpiDope.
        full_sequence: Secuencia completa de aminoacidos, usada para
            reconstruir la subsecuencia de regiones fusionadas cuyo span
            puede exceder el de cualquiera de los dos motores por separado.
        min_length: Filtro de longitud inquebrantable (9 aa por defecto,
            mismo footprint minimo de union a MHC-II que exige NetMHCIIpan;
            deliberadamente NO expuesto como parametro de formulario, ver
            docstring del modulo original).

    Returns:
        Lista de dicts: ``start``, ``end``, ``length``, ``sequence``,
        ``origen`` ('Consenso'/'BepiPred'/'EpiDope'), ``bepipred_score``/
        ``epidope_score`` (``None`` si ese motor no contribuyo a la region) y
        ``bepipred_region``/``epidope_region`` (coordenadas de origen,
        separadas por ';' si se fusiono mas de una region del mismo motor).
    """
    intervals = []
    for start, end, score in bepipred_intervals:
        intervals.append((start, end, "BepiPred", score))
    for start, end, score in epidope_intervals:
        intervals.append((start, end, "EpiDope", score))

    if not intervals:
        return []

    intervals.sort(key=lambda iv: (iv[0], iv[1]))

    def _new_bucket_map(start: int, end: int, source: str, score: float) -> dict:
        return {source: {"scores": [score], "regions": [f"{start}-{end}"]}}

    first_start, first_end, first_source, first_score = intervals[0]
    merged_groups = [[first_start, first_end, _new_bucket_map(first_start, first_end, first_source, first_score)]]

    for start, end, source, score in intervals[1:]:
        group = merged_groups[-1]
        if start <= group[1]:  # solapamiento: comparte al menos un residuo con el grupo abierto
            group[1] = max(group[1], end)
            bucket = group[2].setdefault(source, {"scores": [], "regions": []})
            bucket["scores"].append(score)
            bucket["regions"].append(f"{start}-{end}")
        else:
            merged_groups.append([start, end, _new_bucket_map(start, end, source, score)])

    records = []
    for group_start, group_end, members in merged_groups:
        length = group_end - group_start + 1
        if length < min_length:
            continue

        sources = set(members.keys())
        if sources == {"BepiPred", "EpiDope"}:
            origen = "Consenso"
        elif sources == {"BepiPred"}:
            origen = "BepiPred"
        else:
            origen = "EpiDope"

        bp_info = members.get("BepiPred")
        ed_info = members.get("EpiDope")
        records.append(
            {
                "start": group_start,
                "end": group_end,
                "length": length,
                "sequence": full_sequence[group_start - 1 : group_end] if full_sequence else "",
                "origen": origen,
                "bepipred_score": (sum(bp_info["scores"]) / len(bp_info["scores"])) if bp_info else None,
                "epidope_score": (sum(ed_info["scores"]) / len(ed_info["scores"])) if ed_info else None,
                "bepipred_region": ";".join(bp_info["regions"]) if bp_info else "",
                "epidope_region": ";".join(ed_info["regions"]) if ed_info else "",
            }
        )
    return records
