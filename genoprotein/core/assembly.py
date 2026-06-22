from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from genoprotein.core.orf import find_orfs, OrfResult
from genoprotein.core.sequence import translate_dna


@dataclass
class AssemblyResult:
    query_sequence: str
    assembled_cds: str
    protein_sequence: str
    orfs_found: list[OrfResult] = field(default_factory=list)
    confidence: float = 0.0
    warnings: list[str] = field(default_factory=list)


class ProteinReconstructor:
    def __init__(
        self,
        sequence: str,
        sequence_type: str = "auto",
        min_orf_length: int = 30,
    ):
        self.raw_sequence = sequence.upper().replace(" ", "").replace("\n", "")
        self.sequence_type = sequence_type
        self.min_orf_length = min_orf_length
        self._result: Optional[AssemblyResult] = None
        self._reconstruct()

    @property
    def result(self) -> AssemblyResult:
        if self._result is None:
            raise RuntimeError("Reconstruction not run")
        return self._result

    @property
    def protein_sequence(self) -> str:
        return self.result.protein_sequence

    @property
    def assembled_cds(self) -> str:
        return self.result.assembled_cds

    def _reconstruct(self) -> None:
        seq = self.raw_sequence
        if not seq:
            raise ValueError("Empty sequence provided for reconstruction")
        warnings: list[str] = []
        seq_type = self.sequence_type

        if seq_type == "auto":
            seq_type = self._detect_type(seq)

        if seq_type == "protein":
            self._result = AssemblyResult(
                query_sequence=seq,
                assembled_cds="",
                protein_sequence=seq,
                confidence=0.8,
            )
            return

        orfs = find_orfs(seq, min_length=self.min_orf_length)

        if not orfs:
            orfs = find_orfs(seq, min_length=1)
            if orfs:
                longest = orfs[0]
                if longest.length < 10:
                    warnings.append("No significant ORF found. Sequence may be non-coding.")
                else:
                    warnings.append(f"Only short ORF found ({longest.length}aa)")
            else:
                warnings.append("No ORF found in any frame")
                self._result = AssemblyResult(
                    query_sequence=seq,
                    assembled_cds=seq,
                    protein_sequence=translate_dna(seq),
                    warnings=warnings,
                    confidence=0.1,
                )
                return

        longest = orfs[0]
        cds = longest.nucleotide_sequence
        protein = longest.protein_sequence

        confidence = min(1.0, 0.3 + longest.length / 500)
        if protein.startswith("M"):
            confidence = min(1.0, confidence + 0.1)
        if protein.endswith("*"):
            confidence = min(1.0, confidence + 0.1)

        self._result = AssemblyResult(
            query_sequence=seq,
            assembled_cds=cds,
            protein_sequence=protein.rstrip("*"),
            orfs_found=orfs,
            confidence=round(confidence, 2),
            warnings=warnings,
        )

    @staticmethod
    def _detect_type(seq: str) -> str:
        valid_nt = set("ATGCNatgcn")
        valid_aa = set("ARNDCQEGHILKMFPSTWYVBZX*-arndcqeghilkmfpstwyvbz")
        nt_ratio = sum(1 for c in seq if c in valid_nt) / max(len(seq), 1)
        if nt_ratio > 0.85:
            return "nucleotide"
        aa_ratio = sum(1 for c in seq if c in valid_aa) / max(len(seq), 1)
        if aa_ratio > 0.85:
            return "protein"
        return "nucleotide"

    def splice_in(self, insert_sequence: str, position: int) -> str:
        cds = self.assembled_cds
        if position < 0 or position > len(cds):
            raise ValueError(f"Position {position} out of range [0, {len(cds)}]")
        return cds[:position] + insert_sequence.upper() + cds[position:]

    def splice_out(self, start: int, end: int) -> str:
        cds = self.assembled_cds
        if start < 0 or end > len(cds) or start >= end:
            raise ValueError(f"Invalid range ({start}, {end}) for CDS of length {len(cds)}")
        removed = cds[start:end]
        warnings = self._result.warnings if self._result else []
        warnings.append(f"Spliced out [{start}:{end}] ({len(removed)}bp)")
        return cds[:start] + cds[end:]
