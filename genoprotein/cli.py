from __future__ import annotations

import sys
from pathlib import Path

import click

from genoprotein import __version__
from genoprotein.core.assembly import ProteinReconstructor
from genoprotein.core.orf import find_orfs
from genoprotein.core.sequence import translate_dna, reverse_complement, gc_content
from genoprotein.io.formats import read_fasta, write_fasta, SequenceRecord
from genoprotein.splice.operations import splice_in, splice_out, splice_replace, splice_fusion
from genoprotein.splice.design import design_fusion, add_tag


@click.group()
@click.version_option(version=__version__, prog_name="genoprotein")
def cli():
    """Reconstruct, analyze, and manipulate proteins from genomic data."""


@cli.command()
@click.argument("sequence_or_file")
@click.option("--min-orf", default=30, help="Minimum ORF length in aa", type=int)
@click.option("--file", "is_file", is_flag=True, help="Treat input as FASTA file path")
@click.option("--json", "json_output", is_flag=True, help="JSON output (terminal)")
def reconstruct(sequence_or_file, min_orf, is_file, json_output):
    """Reconstruct protein from a genomic sequence."""
    if is_file:
        records = read_fasta(sequence_or_file)
        for rec in records:
            reconstructor = ProteinReconstructor(rec.sequence, min_orf_length=min_orf)
            _print_reconstruction(reconstructor, rec.id, json_output)
    else:
        reconstructor = ProteinReconstructor(sequence_or_file, min_orf_length=min_orf)
        _print_reconstruction(reconstructor, "input", json_output)


def _print_reconstruction(reconstructor, label, json_output):
    r = reconstructor.result
    if json_output:
        import json
        data = {
            "label": label,
            "protein_sequence": r.protein_sequence,
            "cds_length": len(r.assembled_cds),
            "protein_length": len(r.protein_sequence),
            "confidence": r.confidence,
            "orfs_found": [
                {"start": o.start, "end": o.end, "length": o.length}
                for o in r.orfs_found[:5]
            ],
            "warnings": r.warnings,
        }
        click.echo(json.dumps(data, indent=2))
        return
    click.echo(f"--- {label} ---")
    click.echo(f"Protein ({len(r.protein_sequence)}aa): {r.protein_sequence}")
    click.echo(f"CDS length: {len(r.assembled_cds)}bp")
    click.echo(f"Confidence: {r.confidence}")
    if r.warnings:
        for w in r.warnings:
            click.echo(f"Warning: {w}")
    if r.orfs_found:
        click.echo(f"Top ORFs: {r.orfs_found[0]}")


@cli.command()
@click.argument("sequence")
@click.option("--min-len", default=30, help="Minimum ORF length", type=int)
def orfs(sequence, min_len):
    """Find open reading frames in a sequence."""
    results = find_orfs(sequence, min_length=min_len)
    if not results:
        click.echo("No ORFs found.")
        return
    for i, orf in enumerate(results[:10]):
        click.echo(
            f"ORF {i+1}: [{orf.start}-{orf.end}] "
            f"frame={orf.frame} strand={orf.strand} "
            f"{orf.length}aa"
        )
        click.echo(f"  Protein: {orf.protein_sequence[:60]}")


@cli.command()
@click.argument("sequence")
@click.option("--frame", default=0, help="Reading frame (0,1,2)", type=int)
@click.option("--to-stop", is_flag=True, help="Stop at first stop codon")
def translate(sequence, frame, to_stop):
    """Translate a DNA sequence to protein."""
    result = translate_dna(sequence, frame=frame, to_stop=to_stop)
    click.echo(result)


@cli.command()
@click.argument("sequence")
def revcomp(sequence):
    """Reverse complement of a DNA sequence."""
    click.echo(reverse_complement(sequence))


@cli.command()
@click.argument("sequence")
def gc(sequence):
    """Calculate GC content of a sequence."""
    pct = gc_content(sequence)
    click.echo(f"GC content: {pct:.1f}%")


@cli.command()
@click.argument("sequence")
@click.argument("insert")
@click.argument("position", type=int)
@click.option("--frame/--no-frame", default=True, help="Maintain reading frame")
def splice_in_cmd(sequence, insert, position, frame):
    """Splice a sequence into another at a position."""
    result = splice_in(sequence, insert, position, in_frame=frame)
    click.echo(result)


@cli.command()
@click.argument("sequence")
@click.argument("start", type=int)
@click.argument("end", type=int)
@click.option("--frame/--no-frame", default=True, help="Maintain reading frame")
def splice_out_cmd(sequence, start, end, frame):
    """Remove a region from a sequence."""
    result = splice_out(sequence, start, end, in_frame=frame)
    click.echo(result)


@cli.command()
@click.argument("gene_a")
@click.argument("gene_b")
@click.option("--linker", default="GS", help="Linker type (GS, G4S, P2A, T2A)")
@click.option("--name-a", default="GeneA", help="Name of gene A")
@click.option("--name-b", default="GeneB", help="Name of gene B")
def fuse(gene_a, gene_b, linker, name_a, name_b):
    """Design a gene fusion between two sequences."""
    design = design_fusion(gene_a, gene_b, linker=linker, gene_a_name=name_a, gene_b_name=name_b)
    click.echo(f"Fusion: {design.gene_a_name} -> {design.gene_b_name}")
    click.echo(f"Nucleotide ({len(design.nucleotide_sequence)}bp):")
    click.echo(design.nucleotide_sequence)
    click.echo(f"Protein ({len(design.protein_sequence)}aa):")
    click.echo(design.protein_sequence)


@cli.command()
@click.argument("input_file")
@click.argument("output_file")
@click.option("--format", "fmt", default="fasta", help="Output format (fasta)")
@click.option("--action", type=click.Choice(["translate", "revcomp", "filter"]), default=None)
@click.option("--min-len", default=0, type=int, help="Minimum sequence length filter")
@click.option("--gene", help="Gene name to annotate")
def process(input_file, output_file, fmt, action, min_len, gene):
    """Process a FASTA file with various operations."""
    records = read_fasta(input_file)
    processed: list[SequenceRecord] = []

    for rec in records:
        if min_len and len(rec.sequence) < min_len:
            continue

        if action == "translate":
            new_seq = translate_dna(rec.sequence)
            new_id = f"{rec.id}_translated"
            processed.append(SequenceRecord(id=new_id, description="Translation", sequence=new_seq))
        elif action == "revcomp":
            new_seq = reverse_complement(rec.sequence)
            new_id = f"{rec.id}_revcomp"
            processed.append(SequenceRecord(id=new_id, description="Reverse complement", sequence=new_seq))
        else:
            processed.append(rec)

    if gene:
        from genoprotein.utils.databases import GENE_LOOKUP
        if gene.upper() in GENE_LOOKUP:
            info = GENE_LOOKUP[gene.upper()]
            for rec in processed:
                rec.description += f" | Gene: {gene} ({info['description']})"

    write_fasta(processed, output_file)
    click.echo(f"Wrote {len(processed)} records to {output_file}")


def main():
    cli()


if __name__ == "__main__":
    main()
