#!/usr/bin/env python3
"""
Parses GRF output (grf-main) and creates GFF files for visualization of repeats.
This script processes perfect and imperfect repeat identifications from GRF
and converts them into GFF format for genomic browsers.
"""

import logging
import os
import sys
from pathlib import Path
from typing import List, Tuple

import pandas as pd
from BCBio import GFF
from Bio import SeqIO
from Bio.SeqFeature import ExactPosition, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s", handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)


def parse_grf_output_w_ranges(header: str) -> List:
    """
    Parses a GRF header that contains range coordinates.

    Format example: '>1:0-190543:951:190482:14m'

    Args:
        header (str): The GRF header string to parse.

    Returns:
        List: A list containing [record_id, rep1_start, rep1_end, rep2_start, rep2_end, repeat_len].
    """
    first_split = header.split(":")
    record_id = first_split[0][1:]
    
    # Extract repeat length from the last part (e.g., '14m' -> 14)
    repeat_len = int(first_split[-1].split("m")[0])
    
    # Extract range start/end from the second part (e.g., '0-190543')
    range_coords = first_split[1].split("-")
    range_start = int(range_coords[0])
    
    # Repeat offsets within the range
    repeat_1_start_in_range = int(first_split[2])
    repeat_2_end_in_range = int(first_split[3])

    # Calculate absolute coordinates on the chromosome/contig
    repeat_1_start_in_chrom = range_start + repeat_1_start_in_range
    repeat_1_end_in_chrom = repeat_1_start_in_chrom + repeat_len

    repeat_2_end_in_chrom = range_start + repeat_2_end_in_range
    repeat_2_start_in_chrom = repeat_2_end_in_chrom - repeat_len

    return [
        record_id,
        repeat_1_start_in_chrom,
        repeat_1_end_in_chrom,
        repeat_2_start_in_chrom,
        repeat_2_end_in_chrom,
        repeat_len,
    ]


def parse_grf_output_no_ranges(header: str) -> List:
    """
    Parses a GRF header that does not contain range coordinates.

    Used for newer GRF versions where headers are relative to the record itself.
    Format example: '>IPFHMEHC_00036_gene:95122:101284:29m'

    Args:
        header (str): The GRF header string to parse.

    Returns:
        List: A list containing [record_id, rep1_start, rep1_end, rep2_start, rep2_end, repeat_len].
    """
    first_split = header.split(":")
    record_id = first_split[0][1:]
    repeat_len = int(first_split[-1].split("m")[0])
    
    repeat_1_start = int(first_split[1])
    repeat_1_end = repeat_1_start + repeat_len - 1  # Inclusive coordinates
    
    repeat_2_end = int(first_split[2])
    repeat_2_start = repeat_2_end - repeat_len + 1  # Inclusive coordinates
    
    return [record_id, repeat_1_start, repeat_1_end, repeat_2_start, repeat_2_end, repeat_len]


def create_gff_record_with_features(
    features_df: pd.DataFrame, record: SeqRecord, strand: int = 1, feature_type: str = "direct_repeat"
) -> SeqRecord:
    """
    Appends repeat features to a SeqRecord from a coordinates DataFrame.

    Args:
        features_df (pd.DataFrame): DataFrame with columns [record_id, start_1, end_1, start_2, end_2, length].
        record (SeqRecord): The biopython SeqRecord to which features will be added.
        strand (int, optional): The strand for the GFF features. Defaults to 1.
        feature_type (str, optional): The GFF feature type label. Defaults to "direct_repeat".

    Returns:
        SeqRecord: The record with appended features.
    """
    # The record_id in features_df is expected to have a "_gene" suffix relative to the assembly record.id
    features_in_record = features_df[features_df.record_id == f"{record.id}_gene"]
    
    for i, (_, row) in enumerate(features_in_record.iterrows(), 1):
        repeat_id = str(i)
        
        # Add the first repeat of the pair
        record.features.append(
            SeqFeature(
                FeatureLocation(ExactPosition(row["start_1"]), ExactPosition(row["end_1"]), strand=strand),
                type=feature_type,
                id=repeat_id,
            )
        )
        # Add the second repeat of the pair
        record.features.append(
            SeqFeature(
                FeatureLocation(ExactPosition(row["start_2"]), ExactPosition(row["end_2"]), strand=strand),
                type=feature_type,
                id=repeat_id,
            )
        )

    return record


def main_process(input_assembly: Path, input_grf: Path, min_len: int) -> List[SeqRecord]:
    """
    Coordinates the loading of data, parsing of repeats, and generation of GFF records.

    Args:
        input_assembly (Path): Path to the assembly FASTA file.
        input_grf (Path): Path to the GRF output file.
        min_len (int): Minimum repeat length to include.

    Returns:
        List[SeqRecord]: A list of SeqRecords populated with repeat features.
    """
    if not input_grf.exists():
        logger.warning(f"GRF input file not found: {input_grf}")
        return []

    with open(input_grf) as f:
        repeat_headers = [line.rstrip() for line in f if line.strip()]

    if not repeat_headers:
        logger.info(f"No repeats found in {input_grf}")
        return []

    # Parse headers and create DataFrame
    gff_rows = [parse_grf_output_no_ranges(line) for line in repeat_headers]
    cols = ["record_id", "start_1", "end_1", "start_2", "end_2", "length"]
    gff_df = pd.DataFrame(data=gff_rows, columns=cols).drop_duplicates()

    # Filter by minimum length
    gff_df = gff_df[gff_df["length"] >= min_len]

    # Load assembly records
    assembly_records = list(SeqIO.parse(input_assembly, "fasta"))

    # Generate GFF records by matching features to assembly records
    gff_records = [create_gff_record_with_features(gff_df, rec) for rec in assembly_records]
    
    # Only return records that actually have features
    return [rec for rec in gff_records if rec.features]


def run_snakemake():
    """Wrapper function to handle Snakemake execution context."""
    with open(snakemake.log[0], "w") as f:
        sys.stderr = sys.stdout = f
        
        # Setup inputs from Snakemake context
        input_dir = Path(snakemake.input[0])
        assembly_dir = Path(snakemake.input[1])
        strain = assembly_dir.name
        
        in_perfect = input_dir / "perfect.spacer.id"
        in_imperfect = input_dir / "imperfect.id"
        in_assembly = assembly_dir / f"{strain}_genomic.faa"
        
        out_perfect = Path(snakemake.output[0])
        out_imperfect = Path(snakemake.output[1])
        
        min_repeat_length = int(snakemake.params[0])

        logger.info(f"Processing strain: {strain}")

        # Process Perfect repeats
        logger.info("Processing perfect repeats...")
        perfect_recs = main_process(in_assembly, in_perfect, min_repeat_length)
        with open(out_perfect, "w") as out:
            GFF.write(recs=perfect_recs, out_handle=out, include_fasta=False)

        # Process Imperfect repeats
        logger.info("Processing imperfect repeats...")
        imperfect_recs = main_process(in_assembly, in_imperfect, min_repeat_length)
        with open(out_imperfect, "w") as out:
            GFF.write(recs=imperfect_recs, out_handle=out, include_fasta=False)


if __name__ == "__main__":
    if "snakemake" in globals():
        run_snakemake()
    else:
        logger.error("This script is intended to be run via Snakemake.")
        sys.exit(1)
