import pytest
import pandas as pd
import numpy as np
import os
import sys
from unittest.mock import MagicMock, patch

# Mock dependencies before imports to avoid ImportError during collection
mock_bio = MagicMock()
mock_seqio = MagicMock()
mock_bcbio = MagicMock()
mock_gff = MagicMock()
mock_seqrecord = MagicMock()

sys.modules["Bio"] = mock_bio
sys.modules["Bio.SeqIO"] = mock_seqio
sys.modules["Bio.SeqRecord"] = mock_seqrecord
sys.modules["BCBio"] = mock_bcbio
sys.modules["BCBio.GFF"] = mock_gff

# Allow Bio to be treated as a package for 'from Bio.SeqRecord import SeqRecord'
mock_bio.__path__ = []

# Add repository root to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

from workflow.scripts.python import flanking_regions

def test_make_bed():
    """Test the make_bed function with various span types."""
    data = {
        "chrom": ["chr1", "chr1", "chr1"],
        "gene_id": ["g1", "g2", "g3"],
        "strand": [1, -1, 1],
        "span_start": [100, 200, 300],
        "span_end": [200, 300, 400],
        "span_over_5_start": [900, np.nan, np.nan],
        "span_over_5_end": [1000, np.nan, np.nan],
        "span_over_3_start": [np.nan, np.nan, 0],
        "span_over_3_end": [np.nan, np.nan, 100]
    }
    collection = pd.DataFrame(data)
    
    bed_normal, bed_5, bed_3 = flanking_regions.make_bed(collection, score=100)
    
    # Check normal BED
    assert len(bed_normal) == 3
    assert bed_normal.iloc[0]["range_start"] == 100
    assert bed_normal.iloc[0]["score"] == "100"
    
    # Check 5'-end crossing BED
    assert len(bed_5) == 1
    assert bed_5.iloc[0]["name"] == "g1"
    assert bed_5.iloc[0]["range_start"] == 900
    
    # Check 3'-end crossing BED
    assert len(bed_3) == 1
    assert bed_3.iloc[0]["name"] == "g3"
    assert bed_3.iloc[0]["range_start"] == 0

def test_make_bed_file_for_rg_basic():
    """Test make_bed_file_for_rg with a simple single contig case."""
    # Mock GFF record
    mock_gene = MagicMock()
    mock_gene.type = "gene"
    mock_gene.id = "gene_001"
    
    # Configure location mocks with .position and int() support
    mock_gene.location.start.position = 500
    mock_gene.location.end.position = 600
    mock_gene.location.start.__int__.return_value = 500
    mock_gene.location.end.__int__.return_value = 600
    mock_gene.location.strand = 1
    
    mock_record = MagicMock()
    mock_record.id = "contig_1"
    mock_record.features = [mock_gene]
    mock_record.seq = "A" * 2000
    
    # Mock RGI results
    rgi_df = pd.DataFrame({
        "ORF_ID": ["gene_001 [some metadata]"],
        "Cut_Off": ["Strict"]
    })
    
    bed_dfs, message = flanking_regions.make_bed_file_for_rg(
        gff_record=mock_record,
        rgi_df=rgi_df,
        dna_len=2000,
        span_len=100,
        is_polypolish=True,
        chr_name="contig_1"
    )
    
    # Expectation: 1 gene found
    assert "1 of 1 resistance genes found" in message
    assert len(bed_dfs) == 3
    
    bed_normal = bed_dfs[0]
    assert len(bed_normal) == 1
    assert bed_normal.iloc[0]["chrom"] == "contig_1"
    assert bed_normal.iloc[0]["range_start"] == 500 - 100 - 1
    assert bed_normal.iloc[0]["range_end"] == 600 + 100

def test_make_bed_file_for_rg_unicycler_naming():
    """Test that chrom name is correctly parsed for Unicycler."""
    mock_gene = MagicMock()
    mock_gene.type = "gene"
    mock_gene.id = "g1"
    mock_gene.location.start.position = 10
    mock_gene.location.end.position = 20
    mock_gene.location.start.__int__.return_value = 10
    mock_gene.location.end.__int__.return_value = 20
    
    mock_record = MagicMock()
    mock_record.id = "abc_def_12" # Expected format for split("_")[-1]
    mock_record.features = [mock_gene]
    mock_record.seq = "A" * 1000
    
    rgi_df = pd.DataFrame({"ORF_ID": ["g1"], "Cut_Off": ["Strict"]})
    
    bed_dfs, _ = flanking_regions.make_bed_file_for_rg(
        mock_record, rgi_df, 1000, 10, False, "ignore_me"
    )
    
    assert bed_dfs[0].iloc[0]["chrom"] == "12"

def test_make_bed_file_for_rg_circular_overlaps():
    """Test handling of spans that cross the 0 boundary."""
    mock_gene = MagicMock()
    mock_gene.type = "gene"
    mock_gene.id = "edge_gene"
    mock_gene.location.start.position = 50
    mock_gene.location.end.position = 100
    mock_gene.location.start.__int__.return_value = 50
    mock_gene.location.end.__int__.return_value = 100
    
    mock_record = MagicMock()
    mock_record.id = "circ"
    mock_record.features = [mock_gene]
    mock_record.seq = "A" * 1000
    
    rgi_df = pd.DataFrame({"ORF_ID": ["edge_gene"], "Cut_Off": ["Strict"]})
    
    # span_len = 100 => span_start = 50 - 100 - 1 = -51
    bed_dfs, _ = flanking_regions.make_bed_file_for_rg(
        mock_record, rgi_df, 1000, 100, True, "circ"
    )
    
    # bed_normal should be clipped to 0
    assert bed_dfs[0].iloc[0]["range_start"] == 0
    
    # bed_5 should contain the wrapped part: -51 + 1000 = 949 to 1000
    assert len(bed_dfs[1]) == 1
    assert bed_dfs[1].iloc[0]["range_start"] == 949
    assert bed_dfs[1].iloc[0]["range_end"] == 1000
