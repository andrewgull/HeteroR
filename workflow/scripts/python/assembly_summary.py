import subprocess
import os
import sys
import logging
from pathlib import Path
from typing import Optional, List

import pandas as pd
import numpy as np
from Bio import SeqIO

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

# Standard columns for the assembly summary
COLUMNS = [
    "Component", "Segments", "Links", "Length", "N50", 
    "Longest_component", "Status", "Strain", "Repeat", 
    "Mult", "Coverage", "Type"
]

def parse_flye_info(info_path: Path, strain: str) -> pd.DataFrame:
    """Parses Flye/Medaka assembly_info.txt."""
    logger.info(f"Parsing Flye info from {info_path}")
    df = pd.read_csv(info_path, sep="\t")
    
    # Rename and drop columns to match standard format
    # Standard Flye headers often include '#' and '.'
    # Headers can be: #seq_name, length, cov., circ., repeat, mult.
    df = df.rename(columns={
        "#seq_name": "Component",
        "seq_name": "Component",
        "length": "Length",
        "cov.": "Coverage",
        "cov": "Coverage",
        "circ.": "Status",
        "circ": "Status",
        "repeat": "Repeat",
        "mult.": "Mult",
        "mult": "Mult"
    })
    df.drop(["alt_group", "graph_path"], axis=1, inplace=True, errors="ignore")
    
    # Map status: Flye uses Y/N, some versions might use C/N
    if "Status" in df.columns:
        df["Status"] = df["Status"].str.replace("Y", "complete").str.replace("C", "complete").str.replace("N", "incomplete")
    
    # Add missing standard columns
    df["Strain"] = strain
    df["Type"] = ["Chromosome"] + ["Plasmid"] * (len(df) - 1)
    for col in ["Segments", "Links", "N50", "Longest_component"]:
        df[col] = np.nan
        
    return df[COLUMNS]

def parse_unicycler_log(assembly_dir: Path, strain: str) -> pd.DataFrame:
    """Parses Unicycler log file to extract summary table."""
    log_path = assembly_dir / "unicycler.log"
    logger.info(f"Parsing Unicycler log from {log_path}")
    
    # Use sed to extract the component table
    # This matches the original logic but encapsulated
    cmd = f"sed -n '/^Component/,/^Polishing/{{p;/^Polishing/q}}' {log_path} | head -n -3"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    lines = result.stdout.splitlines()
    if not lines:
        logger.warning(f"No summary found in Unicycler log: {log_path}")
        return pd.DataFrame(columns=COLUMNS)

    # Process lines: skip header and 'total' if present
    data = [line.split() for line in lines]
    if len(data) > 2:
        # Header is data[0], second line usually 'total' or first component
        # Original logic skipped first two elements
        summary_data = data[2:]
    else:
        # Just header and one component
        summary_data = data[1:]

    col_names = ["Component", "Segments", "Links", "Length", "N50", "Longest_component", "Status"]
    df = pd.DataFrame(summary_data, columns=col_names)
    
    # Add missing standard columns
    df["Strain"] = strain
    df["Type"] = ["Chromosome"] + ["Plasmid"] * (len(df) - 1)
    for col in ["Repeat", "Mult", "Coverage"]:
        df[col] = np.nan
        
    return df[COLUMNS]

def parse_plasmid_assembly(plasmid_dir: Path, strain: str) -> Optional[pd.DataFrame]:
    """Parses plasmid assembly scaffolds.fasta if it exists."""
    plasmid_fasta = plasmid_dir / "scaffolds.fasta"
    if not (plasmid_fasta.exists() and plasmid_fasta.stat().st_size > 0):
        logger.info(f"No plasmid assembly found at {plasmid_fasta}")
        return None

    logger.info(f"Parsing plasmid assembly from {plasmid_fasta}")
    lengths = [len(seq) for seq in SeqIO.parse(plasmid_fasta, "fasta")]
    
    df = pd.DataFrame({
        "Length": lengths,
        "Longest_component": lengths,
        "Segments": 1,
        "Status": "scaffold",
        "Strain": strain,
        "Type": "Plasmid"
    })
    
    # Add other columns as NaNs/Nones
    for col in ["Component", "Links", "N50", "Repeat", "Mult", "Coverage"]:
        df[col] = np.nan
        
    return df[COLUMNS]

def main():
    # Setup paths from Snakemake
    input_assembly = Path(snakemake.input[0])
    input_plasmid = Path(snakemake.input[1])
    strain_pos = int(snakemake.params[0])
    output_file = Path(snakemake.output[0])
    log_file = Path(snakemake.log[0])

    # Redirect stderr/stdout to log file
    with open(log_file, "w") as f:
        sys.stderr = sys.stdout = f
        
        # Determine strain name
        # Path parts: results/assemblies/STRAIN -> index strain_pos
        # Note: parts is 0-indexed. If strain_pos was 2 for results/assemblies/DA..., parts is ('results', 'assemblies', 'DA...')
        strain = input_assembly.parts[strain_pos]
        logger.info(f"Processing strain: {strain}")

        # 1. Parse main assembly
        flye_info = input_assembly / "assembly_info.txt"
        if flye_info.exists():
            main_df = parse_flye_info(flye_info, strain)
        else:
            main_df = parse_unicycler_log(input_assembly, strain)

        # 2. Parse plasmid assembly
        plasmid_df = parse_plasmid_assembly(input_plasmid, strain)

        # 3. Join and Finalize
        if plasmid_df is not None:
            final_df = pd.concat([main_df, plasmid_df], ignore_index=True)
            # Re-index components to be sequential across all
            final_df["Component"] = range(1, len(final_df) + 1)
        else:
            final_df = main_df

        # Write output
        logger.info(f"Writing summary to {output_file}")
        final_df.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":
    main()
