from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from genoprotein.core.sequence import translate_dna


class Orientation(Enum):
    N_TO_C = "N_to_C"
    C_TO_N = "C_to_N"


class LinkerType(Enum):
    FLEXIBLE = "GGTGGCGGTGGCTCT"
    RIGID = "GCGGCGGCGGAGGCGGA"
    HELICAL = "GAGGCTGCCGCCGAGGCAGCGGCTGAG"
    CUSTOM = ""


FLEXIBLE_LINKERS: dict[str, str] = {
    "GS": "GGTTCT",
    "G4S": "GGTGGCGGTGGCTCT",
    "G4S_3x": "GGTGGCGGTGGCTCTGGTGGCGGTGGCTCTGGTGGCGGTGGCTCT",
    "P2A": "GCCACCAACTTCAGCCTGCTGAAGCAGGCCGGCGACGTGGAGGAGAACCCCGGCCCT",
    "T2A": "GAGGGCAGAGGAAGTCTGCTAACATGCGGTGACGTCGAGGAGAATCCTGGCCCA",
}


COMMON_TAGS: dict[str, str] = {
    "His6": "CATCACCATCACCATCAC",
    "FLAG": "GACTACAAGGACGACGACGACAAA",
    "HA": "TACCCCTACGACGTGCCCGACTACGCC",
    "Myc": "GAGCAGAAACTGATCTCCGAGGAAGACCTG",
    "StrepII": "TGGAGCCACCCGCAGTTCGAAAAA",
    "GFP": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA",
    "mCherry": "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCGATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGTAA",
}


@dataclass
class FusionDesign:
    gene_a_name: str
    gene_b_name: str
    nucleotide_sequence: str
    protein_sequence: str
    linker_used: str
    orientation: Orientation
    linker_type: LinkerType

    def __repr__(self) -> str:
        return (
            f"FusionDesign({self.gene_a_name} -> {self.gene_b_name}, "
            f"len={len(self.nucleotide_sequence)}bp, "
            f"protein={len(self.protein_sequence)}aa"
            f"{' [linker]' if self.linker_used else ''})"
        )


@dataclass
class DesignSpec:
    sequence: str
    description: str = ""
    gc_content: float = 0.0
    warnings: list[str] = field(default_factory=list)


def design_fusion(
    gene_a: str,
    gene_b: str,
    linker: str = "GS",
    orientation: Orientation = Orientation.N_TO_C,
    gene_a_name: str = "GeneA",
    gene_b_name: str = "GeneB",
    remove_stop_a: bool = True,
    remove_start_b: bool = True,
) -> FusionDesign:
    a = gene_a.upper()
    b = gene_b.upper()

    if linker in FLEXIBLE_LINKERS:
        linker_seq = FLEXIBLE_LINKERS[linker]
        linker_type = LinkerType.FLEXIBLE
    else:
        linker_seq = linker.upper()
        linker_type = LinkerType.CUSTOM

    if remove_stop_a:
        for codon in ["TAA", "TAG", "TGA"]:
            if a.endswith(codon):
                a = a[: -len(codon)]
                break

    if remove_start_b and b.startswith("ATG"):
        b = b[3:]

    if orientation == Orientation.N_TO_C:
        fused_nt = a + linker_seq + b
    else:
        fused_nt = b + linker_seq + a

    fused_protein = translate_dna(fused_nt, to_stop=True)

    return FusionDesign(
        gene_a_name=gene_a_name,
        gene_b_name=gene_b_name,
        nucleotide_sequence=fused_nt,
        protein_sequence=fused_protein,
        linker_used=linker_seq,
        orientation=orientation,
        linker_type=linker_type,
    )


def add_linker(
    sequence: str,
    linker: str = "G4S",
    position: str = "end",
) -> str:
    if linker in FLEXIBLE_LINKERS:
        linker_seq = FLEXIBLE_LINKERS[linker]
    else:
        linker_seq = linker.upper()
    seq = sequence.upper()
    if position == "end":
        return seq + linker_seq
    elif position == "start":
        return linker_seq + seq
    else:
        raise ValueError(f"Position must be 'start' or 'end', got {position}")


def add_tag(
    sequence: str,
    tag: str = "His6",
    position: str = "C_terminal",
    linker: str = "",
) -> str:
    seq = sequence.upper()
    if tag in COMMON_TAGS:
        tag_seq = COMMON_TAGS[tag]
    else:
        tag_seq = tag.upper()

    linker_seq = ""
    if linker:
        if linker in FLEXIBLE_LINKERS:
            linker_seq = FLEXIBLE_LINKERS[linker]
        else:
            linker_seq = linker.upper()

    if position.lower() in ("c_terminal", "c", "cterm"):
        return seq + linker_seq + tag_seq if linker_seq else seq + tag_seq
    elif position.lower() in ("n_terminal", "n", "nterm"):
        return (tag_seq + linker_seq + seq) if linker_seq else tag_seq + seq
    else:
        raise ValueError(f"Position must be C_terminal or N_terminal, got {position}")
