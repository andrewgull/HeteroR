#!/usr/bin/env python3
"""
Script to create BED files for a given assembly and extract flanking regions.
Requires: biopython, bcbio-gff, pandas, numpy
"""

import os
import sys
import logging
from pathlib import Path
from typing import List, Tuple, Optional

import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from BCBio import GFF

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s", handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

# Standard BED columns
BED_COLUMNS = ["chrom", "range_start", "range_end", "name", "score", "strand"]


def make_bed(collection: pd.DataFrame, score: int = 0) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Creates BED-formatted DataFrames for normal, 5'-end crossing, and 3'-end crossing spans.
    """

    def create_bed_df(df, start_col, end_col):
        bed = pd.DataFrame()
        if not df.empty:
            bed["chrom"] = df["chrom"]
            bed["range_start"] = df[start_col]
            bed["range_end"] = df[end_col]
            bed["name"] = df["gene_id"]
            bed["score"] = str(score)
            bed["strand"] = df["strand"]
        return bed

    # normal coordinates
    bed_normal = create_bed_df(collection, "span_start", "span_end")

    # spans crossing 5-end
    bed_5 = create_bed_df(collection, "span_over_5_start", "span_over_5_end")
    bed_5 = bed_5[bed_5["range_start"].notnull()]

    # spans crossing 3-end
    bed_3 = create_bed_df(collection, "span_over_3_start", "span_over_3_end")
    bed_3 = bed_3[bed_3["range_end"].notnull()]

    return bed_normal, bed_5, bed_3


def make_bed_file_for_rg(
    gff_record, rgi_df: pd.DataFrame, dna_len: int, span_len: int, is_polypolish: bool, chr_name: str
) -> Tuple[List[pd.DataFrame], str]:
    """
    Makes BED files for genomic ranges with resistance genes.
    """
    genes = [feature for feature in gff_record.features if feature.type == "gene"]
    item_id = gff_record.id

    # find median gene length
    gene_lengths = [gene.location.end.position - gene.location.start.position for gene in genes]
    median_gene = pd.Series(gene_lengths).median() if genes else 0

    # adjust span_len if record is too short
    record_length = len(gff_record.seq)
    msg_len = ""
    if record_length < span_len * 2:
        span_len = round((record_length - median_gene) / 2)
        msg_len = f"span length was adjusted to {span_len} in record {item_id}\n"

    # Match RGI results to GFF gene IDs
    resistance_genes_coords = []
    rgi_ids = set(orf.split(" ")[0] for orf in rgi_df["ORF_ID"])
    for gene in genes:
        if any(rid in gene.id for rid in rgi_ids):
            resistance_genes_coords.append(gene)

    msg_count = f"In record {item_id} {len(resistance_genes_coords)} of {len(rgi_df)} resistance genes found\n"

    if not resistance_genes_coords:
        return [pd.DataFrame(columns=BED_COLUMNS)], msg_len + msg_count

    rg_collection = []
    for gene in resistance_genes_coords:
        chrom_name = chr_name
        if not is_polypolish:
            # Unicycler/SPAdes naming convention
            chrom_name = item_id.split("_")[-1]

        start = int(gene.location.start)
        end = int(gene.location.end)
        span_start = start - span_len - 1
        span_end = end + span_len

        rg_collection.append([chrom_name, gene.id, start, end, span_start, span_end, int(gene.location.strand)])

    rg_ranges = pd.DataFrame(
        rg_collection, columns=["chrom", "gene_id", "gene_start", "gene_end", "span_start", "span_end", "strand"]
    )

    # Handle circular overlaps
    rg_ranges["span_over_5_start"] = np.where(rg_ranges["span_start"] < 0, rg_ranges["span_start"] + dna_len, np.nan)
    rg_ranges["span_over_5_end"] = np.where(rg_ranges["span_over_5_start"].isna(), np.nan, dna_len)
    rg_ranges["span_over_3_end"] = np.where(rg_ranges["span_end"] > dna_len, rg_ranges["span_end"] - dna_len, np.nan)
    rg_ranges["span_over_3_start"] = np.where(rg_ranges["span_over_3_end"].isna(), np.nan, 0)

    bed_dfs = list(make_bed(rg_ranges, score=0))

    # Final coordinate clipping and type conversion for BED compliance
    for df in bed_dfs:
        if not df.empty:
            df["range_start"] = np.clip(df["range_start"].astype(float), 0, dna_len).astype("int64")
            df["range_end"] = np.clip(df["range_end"].astype(float), 0, dna_len).astype("int64")

    return bed_dfs, msg_len + msg_count


def load_inputs(
    assembly_path: Path, gff_dir: Path, rgi_path: Path, min_plasmid_size: int, strain: str
) -> Tuple[List[SeqRecord], List[SeqRecord], pd.DataFrame, bool]:
    """
    Loads and filters assembly, GFF, and RGI data.
    """
    # Check assembler type from first record
    try:
        first_record = next(SeqIO.parse(assembly_path, "fasta"))
        is_polypolish = "contig" in first_record.id and "polypolish" in first_record.id
    except StopIteration:
        raise ValueError(f"Assembly file {assembly_path} is empty")

    # Load and filter RGI
    rgi = pd.read_csv(rgi_path, sep="\t")
    rgi_strict = rgi[rgi["Cut_Off"] != "Loose"]

    # Load and filter GFF
    gff_path = gff_dir / f"{strain}_genomic.gff"
    if not gff_path.exists():
        raise FileNotFoundError(f"GFF file not found: {gff_path}")

    with open(gff_path) as f:
        gff_records = [rec for rec in GFF.parse(f) if len(rec.seq) >= min_plasmid_size]

    # Load assembly records
    all_assembly_records = list(SeqIO.parse(assembly_path, "fasta"))

    return all_assembly_records, gff_records, rgi_strict, is_polypolish


def match_assembly_gff(
    assembly_records: List[SeqRecord], gff_records: List[SeqRecord]
) -> List[Tuple[SeqRecord, SeqRecord]]:
    """
    Pairs assembly contigs with GFF records by length.
    """
    # Filter assembly to match GFF records by length
    assembly_by_len = {len(rec): rec for rec in assembly_records}
    matched_pairs = []

    for gff_rec in gff_records:
        rec_len = len(gff_rec.seq)
        if rec_len in assembly_by_len:
            matched_pairs.append((assembly_by_len[rec_len], gff_rec))
        else:
            logger.warning(f"No assembly record found matching GFF record length {rec_len}")

    # Sort by length for consistency (descending)
    return sorted(matched_pairs, key=lambda x: len(x[0].seq), reverse=True)


def save_results(bed_lol: List[List[pd.DataFrame]], output_paths: List[Path], messages: List[str], log_path: Path):
    """
    Concatenates BED DataFrames and saves them along with the log messages.
    """
    for i, output_path in enumerate(output_paths):
        if bed_lol[i]:
            final_df = pd.concat(bed_lol[i], ignore_index=True)
            final_df.to_csv(output_path, sep="\t", index=False, header=False)
        else:
            # Ensure output file exists even if empty
            pd.DataFrame(columns=BED_COLUMNS).to_csv(output_path, sep="\t", index=False, header=False)

    with open(log_path, "w") as log:
        log.writelines(messages)


def main():
    # Snakemake parameters
    in_assembly = Path(snakemake.input[0])
    # Extract strain from results/assemblies/STRAIN/assembly.fasta
    strain = in_assembly.parent.name

    in_gff_dir = Path(snakemake.input[1])
    in_rgi = Path(snakemake.input[2])
    range_len = int(snakemake.params[0])
    min_plasmid_size = int(snakemake.params[1])

    output_files = [Path(p) for p in snakemake.output]
    log_file = Path(snakemake.log[0])

    # 1. Load data
    assembly_records, gff_records, rgi_strict, is_polypolish = load_inputs(
        in_assembly, in_gff_dir, in_rgi, min_plasmid_size, strain
    )

    # 2. Match records
    matched_pairs = match_assembly_gff(assembly_records, gff_records)

    # 3. Process records
    bed_lol = [[], [], []]
    all_messages = []

    for assembly_rec, gff_rec in matched_pairs:
        record_len = len(assembly_rec.seq)
        record_id = assembly_rec.id

        record_bed_dfs, message = make_bed_file_for_rg(
            gff_record=gff_rec,
            rgi_df=rgi_strict,
            dna_len=record_len,
            span_len=range_len,
            is_polypolish=is_polypolish,
            chr_name=record_id,
        )

        for i, df in enumerate(record_bed_dfs):
            bed_lol[i].append(df)
        all_messages.append(message)

    # 4. Save outputs
    save_results(bed_lol, output_files, all_messages, log_file)


if __name__ == "__main__":
    main()
