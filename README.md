# GenoProtein

Reconstruct, analyze, and manipulate proteins from partial genomic data with gene splicing operations.

## Installation

```bash
pip install genoprotein
```

Or from source:

```bash
git clone https://github.com/YOUR_USERNAME/genoprotein.git
cd genoprotein
pip install -e ".[dev]"
```

## Quick Start

### Reconstruct a protein from a partial genomic sequence

```python
from genoprotein import ProteinReconstructor

# Reconstruct from a partial CDS sequence
reconstructor = ProteinReconstructor("ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA")
print(reconstructor.protein_sequence)
```

### Splice a gene into another

```python
from genoprotein.splice import design_fusion

# Create a fusion of GFP and a target gene
fusion = design_fusion(
    gene_a="ATGGTGAGCAAGGGCGAGGAG...",
    gene_b="ATGGCCGAGCGCGACCTGATC...",
    linker="GGATCC",  # BamHI linker
    orientation="N_to_C",
)
print(fusion.nucleotide_sequence)
```

### Find open reading frames

```python
from genoprotein.core.orf import find_orfs

orfs = find_orfs("ATGCGTAGCGTGACGTAGCGATGCGTGA...", min_length=50)
for orf in orfs:
    print(f"ORF at {orf.start}-{orf.end}: {orf.protein_sequence}")
```

## Features

- **Protein Reconstruction** — Reconstruct full protein sequences from partial genomic data
- **ORF Detection** — Find open reading frames in any nucleotide sequence
- **Gene Splicing** — Splice genes in, out, fuse, or replace them
- **Designer Fusions** — Design fusion proteins with custom linkers
- **Format Support** — FASTA, GenBank, raw sequences
- **Database Lookup** — Fetch reference sequences from NCBI/UniProt

## API Overview

| Module | Description |
|--------|-------------|
| `genoprotein.core` | Core sequence analysis and protein reconstruction |
| `genoprotein.splice` | Gene splicing operations (in/out/fusion/replace) |
| `genoprotein.io` | Sequence format parsing |
| `genoprotein.utils` | Database and helper utilities |

## License

All Rights Reserved
