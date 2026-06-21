from __future__ import annotations

import re
from typing import Optional

STANDARD_CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate_dna(
    sequence: str,
    codon_table: Optional[dict[str, str]] = None,
    frame: int = 0,
    to_stop: bool = False,
) -> str:
    if codon_table is None:
        codon_table = STANDARD_CODON_TABLE
    seq = sequence[frame:].upper()
    seq = re.sub(r"[^ATGCN]", "N", seq)
    amino_acids: list[str] = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        if "N" in codon:
            amino_acids.append("X")
            continue
        aa = codon_table.get(codon, "X")
        if to_stop and aa == "*":
            break
        amino_acids.append(aa)
    return "".join(amino_acids)


DNA_COMPLEMENT = str.maketrans("ATGCNatgcn", "TACGNtacgn")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(DNA_COMPLEMENT)[::-1]


AA_CODES = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D",
    "Cys": "C", "Glu": "E", "Gln": "Q", "Gly": "G",
    "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
    "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "Asx": "B", "Glx": "Z", "Xaa": "X",
}
AA_CODES.update({v: v for v in "ARNDCQEGHILKMFPSTWYVBZX"})


def gc_content(sequence: str) -> float:
    seq = sequence.upper()
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq) * 100


def is_start_codon(codon: str) -> bool:
    return codon.upper() == "ATG"


def is_stop_codon(codon: str) -> bool:
    return codon.upper() in ("TAA", "TAG", "TGA")
