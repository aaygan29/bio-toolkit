from __future__ import annotations

import json

import click

from genoprotein import __version__
from genoprotein.core.assembly import ProteinReconstructor
from genoprotein.core.orf import find_orfs
from genoprotein.core.sequence import translate_dna, reverse_complement, gc_content
from genoprotein.io.formats import read_fasta, write_fasta, SequenceRecord
from genoprotein.splice.operations import splice_in, splice_out
from genoprotein.splice.design import design_fusion
from genoprotein.repository.store import ProteinRepository
from genoprotein.decoder.matcher import match_partial_sequence
from genoprotein.security.screener import SequenceScreener


@click.group()
@click.version_option(version=__version__, prog_name="genoprotein")
def cli():
    """Reconstruct, analyze, and manipulate proteins from genomic data."""


@cli.command()
@click.argument("sequence_or_file")
@click.option("--min-orf", default=30, help="Minimum ORF length in aa", type=int)
@click.option("--file", "is_file", is_flag=True, help="Treat input as FASTA file path")
@click.option("--json-output", is_flag=True, help="JSON output")
@click.option("--decode", is_flag=True, help="Decode against known protein database")
def reconstruct(sequence_or_file, min_orf, is_file, json_output, decode):
    """Reconstruct protein from a genomic sequence."""
    if is_file:
        records = read_fasta(sequence_or_file)
        for rec in records:
            _do_reconstruct(rec.sequence, rec.id, min_orf, json_output, decode)
    else:
        _do_reconstruct(sequence_or_file, "input", min_orf, json_output, decode)


def _do_reconstruct(sequence, label, min_orf, json_output, decode):
    recon = ProteinReconstructor(sequence, min_orf_length=min_orf)
    r = recon.result

    if json_output:
        data = {
            "label": label,
            "protein_sequence": r.protein_sequence,
            "cds_length": len(r.assembled_cds),
            "protein_length": len(r.protein_sequence),
            "confidence": r.confidence,
            "factors": r.enhanced_confidence.factors if r.enhanced_confidence else {},
            "warnings": r.warnings,
        }
        if decode:
            dec = recon.decode()
            data["decoder"] = {
                "best_gene": dec.best_gene, "candidates": dec.candidates[:5],
            }
        click.echo(json.dumps(data, indent=2))
        return

    click.echo(f"--- {label} ---")
    click.echo(f"Protein ({len(r.protein_sequence)}aa): {r.protein_sequence}")
    click.echo(f"CDS length: {len(r.assembled_cds)}bp")
    if r.enhanced_confidence:
        f = r.enhanced_confidence.factors
        click.echo(f"Confidence: {r.confidence:.4f}  [{', '.join(f'{k}={v:.2f}' for k,v in f.items())}]")
    if r.warnings:
        for w in r.warnings:
            click.echo(f"Warning: {w}")
    if r.orfs_found:
        click.echo(f"Top ORF: {r.orfs_found[0]}")
    if decode:
        dec = recon.decode()
        click.echo(f"Best gene match: {dec.best_gene}")
        if dec.ambiguity:
            click.echo(f"Ambiguity: best={dec.ambiguity.best_candidate} resolved={dec.ambiguity.is_resolved}")
            for d in dec.ambiguity.covered_diagnostic_positions[:5]:
                click.echo(f"  Pos {d['position']}: {d['expected']}->{d['observed']} ({d['description']})")


@cli.command()
@click.argument("sequence")
@click.option("--min-len", default=30, type=int)
def orfs(sequence, min_len):
    """Find open reading frames."""
    for i, orf in enumerate(find_orfs(sequence, min_length=min_len)[:10]):
        click.echo(f"ORF {i+1}: [{orf.start}-{orf.end}] frame={orf.frame} strand={orf.strand} {orf.length}aa  {orf.protein_sequence[:60]}")
    else:
        click.echo("No ORFs found.")


@cli.command()
@click.argument("sequence")
@click.option("--frame", default=0, type=int)
@click.option("--to-stop", is_flag=True)
def translate(sequence, frame, to_stop):
    """Translate DNA to protein."""
    click.echo(translate_dna(sequence, frame=frame, to_stop=to_stop))


@cli.command()
@click.argument("sequence")
def revcomp(sequence):
    """Reverse complement."""
    click.echo(reverse_complement(sequence))


@cli.command()
@click.argument("sequence")
def gc(sequence):
    """GC content."""
    click.echo(f"GC content: {gc_content(sequence):.1f}%")


@cli.command()
@click.argument("sequence")
@click.argument("insert")
@click.argument("position", type=int)
@click.option("--frame/--no-frame", default=True)
def splice_in_cmd(sequence, insert, position, frame):
    click.echo(splice_in(sequence, insert, position, in_frame=frame))


@cli.command()
@click.argument("sequence")
@click.argument("start", type=int)
@click.argument("end", type=int)
@click.option("--frame/--no-frame", default=True)
def splice_out_cmd(sequence, start, end, frame):
    click.echo(splice_out(sequence, start, end, in_frame=frame))


@cli.command()
@click.argument("gene_a")
@click.argument("gene_b")
@click.option("--linker", default="GS")
@click.option("--name-a", default="GeneA")
@click.option("--name-b", default="GeneB")
def fuse(gene_a, gene_b, linker, name_a, name_b):
    """Design a gene fusion."""
    d = design_fusion(gene_a, gene_b, linker=linker, gene_a_name=name_a, gene_b_name=name_b)
    click.echo(f"Fusion: {d.gene_a_name} -> {d.gene_b_name}")
    click.echo(f"Nucleotide ({len(d.nucleotide_sequence)}bp): {d.nucleotide_sequence}")
    click.echo(f"Protein ({len(d.protein_sequence)}aa): {d.protein_sequence}")


@cli.command()
@click.argument("input_file")
@click.argument("output_file")
@click.option("--action", type=click.Choice(["translate", "revcomp"]))
@click.option("--min-len", default=0, type=int)
@click.option("--gene")
def process(input_file, output_file, action, min_len, gene):
    """Process a FASTA file."""
    records = read_fasta(input_file)
    processed = []

    for rec in records:
        if min_len and len(rec.sequence) < min_len:
            continue
        if action == "translate":
            processed.append(SequenceRecord(id=f"{rec.id}_translated", description="Translation", sequence=translate_dna(rec.sequence)))
        elif action == "revcomp":
            processed.append(SequenceRecord(id=f"{rec.id}_revcomp", description="RevComp", sequence=reverse_complement(rec.sequence)))
        else:
            processed.append(rec)

    if gene:
        with ProteinRepository() as db:
            for entry in db.get_by_gene(gene):
                for rec in processed:
                    rec.description += f" | {gene}: {entry.full_name}"

    write_fasta(processed, output_file)
    click.echo(f"Wrote {len(processed)} records to {output_file}")


@cli.group()
def repo():
    """Manage protein repository."""


@repo.command("list")
def repo_list():
    with ProteinRepository() as db:
        for e in db.list_all():
            click.echo(f"{e.accession}  {e.gene_name:8s}  {e.full_name[:50]}")
        click.echo(f"Total: {len(db.list_all())} entries")


@repo.command("info")
@click.argument("gene_or_accession")
def repo_info(gene_or_accession):
    with ProteinRepository() as db:
        entry = db.get_by_accession(gene_or_accession) or (db.get_by_gene(gene_or_accession) or [None])[0]
    if not entry:
        click.echo(f"No entry for {gene_or_accession}")
        return
    click.echo(f"{entry.accession}  {entry.gene_name}  {entry.full_name}  ({len(entry)}aa)")
    if entry.variants:
        click.echo("Variants:")
        for v in entry.variants:
            click.echo(f"  {v['id']:15s} pos={v['position']:5d} {v['ref']}->{v['alt']}  {v.get('description','')}")
    if entry.domains:
        click.echo("Domains:")
        for d in entry.domains:
            click.echo(f"  {d['name']:25s} [{d['start']}-{d['end']}]")


@cli.command()
@click.argument("query")
@click.option("--min-id", default=0.3, type=float)
@click.option("--json-output", is_flag=True)
def decode(query, min_id, json_output):
    """Decode partial sequence against known proteins."""
    result = match_partial_sequence(query, min_identity=min_id)
    if json_output:
        click.echo(json.dumps({
            "query_length": result.query_length, "query_type": result.query_type,
            "matches": [{"gene": m.entry.gene_name, "identity": m.identity, "coverage": m.coverage_fraction}
                        for m in result.matches[:10]],
            "best_gene": result.best_gene, "ambiguity": result.ambiguity,
        }, indent=2))
        return

    if not result.matches:
        click.echo("No matches.")
        return
    click.echo(f"Top: {result.top_match.entry.gene_name} ({result.top_match.identity:.2%}, cov={result.top_match.coverage_fraction:.2%})")
    if result.ambiguity:
        click.echo("Ambiguity: " + "; ".join(result.ambiguity))


@cli.group()
def security():
    """Biosecurity screening."""


@security.command("screen")
@click.argument("sequence")
@click.option("--customer", help="Customer ID")
@click.option("--institution")
@click.option("--use-case")
@click.option("--json-output", is_flag=True)
def screen_cmd(sequence, customer, institution, use_case, json_output):
    """Screen sequence against pathogen database."""
    result = SequenceScreener(customer_id=customer, institution=institution, use_case=use_case).screen(sequence)
    if json_output:
        click.echo(json.dumps({
            "level": result.level.value, "score": result.score,
            "allowed": result.is_allowed, "warnings": result.warnings,
            "matches": result.matches,
        }, indent=2))
        return
    click.echo(f"[{result.level.value.upper()}] score={result.score:.2f} allowed={result.is_allowed}")
    for w in result.warnings:
        click.echo(f"  {w}")
