from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from genoprotein.core.sequence import translate_dna, is_start_codon, is_stop_codon


@dataclass(frozen=True)
class OrfResult:
    start: int
    end: int
    frame: int
    strand: str
    length: int
    nucleotide_sequence: str
    protein_sequence: str

    def __repr__(self) -> str:
        return (
            f"OrfResult(start={self.start}, end={self.end}, "
            f"frame={self.frame}, strand={self.strand}, "
            f"length={self.length}aa, protein={self.protein_sequence[:20]}...)"
        )


def find_orfs(
    sequence: str,
    min_length: int = 30,
    include_nested: bool = False,
) -> list[OrfResult]:
    seq = sequence.upper()
    results: list[OrfResult] = []

    for strand, seq_actual in [("+", seq), ("-", _revcomp(seq))]:
        for frame in range(3):
            i = frame
            while i < len(seq_actual) - 2:
                codon = seq_actual[i : i + 3]
                if len(codon) < 3:
                    break
                if not is_start_codon(codon):
                    i += 3 if include_nested else 1
                    continue

                orf_start = i
                for j in range(i + 3, len(seq_actual) - 2, 3):
                    codon2 = seq_actual[j : j + 3]
                    if len(codon2) < 3:
                        break
                    if is_stop_codon(codon2):
                        orf_end = j + 3
                        orf_len = (orf_end - orf_start) // 3
                        if orf_len >= min_length / 3:
                            nuc = seq_actual[orf_start:orf_end]
                            prot = translate_dna(nuc)
                            results.append(
                                OrfResult(
                                    start=orf_start
                                    if strand == "+"
                                    else len(seq) - orf_end,
                                    end=orf_end
                                    if strand == "+"
                                    else len(seq) - orf_start,
                                    frame=frame,
                                    strand=strand,
                                    length=orf_len,
                                    nucleotide_sequence=nuc,
                                    protein_sequence=prot,
                                )
                            )
                        i = j + 3
                        break
                else:
                    i = len(seq_actual)
                    continue
                if include_nested:
                    i += 3
            else:
                continue

    results.sort(key=lambda r: r.length, reverse=True)
    return results


def _revcomp(seq: str) -> str:
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq))
