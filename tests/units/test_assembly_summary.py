import pytest
import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

# Mock Bio dependency before any imports
mock_bio = MagicMock()
mock_seqio = MagicMock()
mock_bio.SeqIO = mock_seqio
sys.modules["Bio"] = mock_bio
sys.modules["Bio.SeqIO"] = mock_seqio

# Add repository root to sys.path to allow importing from workflow.scripts
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

from workflow.scripts.python import assembly_summary

@pytest.fixture
def temp_assembly_dir(tmp_path):
    """Fixture for a temporary assembly directory."""
    d = tmp_path / "assembly"
    d.mkdir()
    return d

@pytest.fixture
def temp_plasmid_dir(tmp_path):
    """Fixture for a temporary plasmid directory."""
    d = tmp_path / "plasmid"
    d.mkdir()
    return d

def test_parse_flye_info(temp_assembly_dir):
    """Test parsing of Flye assembly_info.txt."""
    info_file = temp_assembly_dir / "assembly_info.txt"
    # Create mock Flye info content
    content = "seq_name\tlength\tcov\tcirc\trepeat\tmult\talt_group\tgraph_path\n" \
              "contig_1\t1000\t50.5\tC\tN\t1\t*\t*\n" \
              "contig_2\t500\t20.0\tN\tY\t2\t*\t*\n"
    info_file.write_text(content)
    
    df = assembly_summary.parse_flye_info(info_file, "STRAIN1")
    
    assert len(df) == 2
    assert list(df["Component"]) == ["contig_1", "contig_2"]
    assert list(df["Status"]) == ["complete", "incomplete"]
    assert list(df["Type"]) == ["Chromosome", "Plasmid"]
    assert df.iloc[0]["Coverage"] == 50.5
    assert df.iloc[0]["Strain"] == "STRAIN1"
    assert all(col in df.columns for col in assembly_summary.COLUMNS)

@patch("subprocess.run")
def test_parse_unicycler_log(mock_run, temp_assembly_dir):
    """Test parsing of Unicycler log."""
    log_file = temp_assembly_dir / "unicycler.log"
    log_file.write_text("dummy") # Just needs to exist
    
    # Mock subprocess output for the sed command
    # Header + Total + 2 Components
    mock_stdout = "Component Segments Links Length N50 Longest_component Status\n" \
                  "Total 3 2 1500 1000 1000 -\n" \
                  "1 1 0 1000 1000 1000 complete\n" \
                  "2 2 2 500 500 500 incomplete\n"
    mock_run.return_value = MagicMock(stdout=mock_stdout)
    
    df = assembly_summary.parse_unicycler_log(temp_assembly_dir, "STRAIN1")
    
    assert len(df) == 2
    assert list(df["Component"]) == ["1", "2"]
    assert list(df["Type"]) == ["Chromosome", "Plasmid"]
    assert df.iloc[0]["Length"] == "1000"
    assert df.iloc[1]["Status"] == "incomplete"
    assert all(col in df.columns for col in assembly_summary.COLUMNS)

@patch("workflow.scripts.python.assembly_summary.SeqIO")
def test_parse_plasmid_assembly(mock_seqio_patch, temp_plasmid_dir):
    """Test parsing of plasmid assembly scaffolds.fasta."""
    fasta_file = temp_plasmid_dir / "scaffolds.fasta"
    fasta_file.write_text("dummy content") # Must exist and have size > 0
    
    # Configure mock for SeqIO.parse
    mock_seqio_patch.parse.return_value = [
        MagicMock(__len__=MagicMock(return_value=4)),
        MagicMock(__len__=MagicMock(return_value=8))
    ]
    
    df = assembly_summary.parse_plasmid_assembly(temp_plasmid_dir, "STRAIN1")
    
    assert len(df) == 2
    assert list(df["Length"]) == [4, 8]
    assert all(df["Type"] == "Plasmid")
    assert all(df["Status"] == "scaffold")
    assert all(col in df.columns for col in assembly_summary.COLUMNS)

def test_parse_plasmid_assembly_not_exists(temp_plasmid_dir):
    """Test plasmid parsing when file is missing."""
    df = assembly_summary.parse_plasmid_assembly(temp_plasmid_dir, "STRAIN1")
    assert df is None

def test_parse_plasmid_assembly_empty(temp_plasmid_dir):
    """Test plasmid parsing when file is empty."""
    fasta_file = temp_plasmid_dir / "scaffolds.fasta"
    fasta_file.write_text("")
    df = assembly_summary.parse_plasmid_assembly(temp_plasmid_dir, "STRAIN1")
    assert df is None
