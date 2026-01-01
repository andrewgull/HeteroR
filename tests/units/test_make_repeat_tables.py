import pytest
import pandas as pd
import os
import sys
from unittest.mock import MagicMock, patch

# Add repository root and script directory to sys.path
repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, repo_root)
script_dir = os.path.join(repo_root, "workflow/scripts/python")
sys.path.insert(0, script_dir)

# Mock missing dependencies before any imports that use them
mock_bcbio = MagicMock()
sys.modules["BCBio"] = mock_bcbio
sys.modules["BCBio.GFF"] = MagicMock()

mock_bio = MagicMock()
sys.modules["Bio"] = mock_bio
sys.modules["Bio.SeqIO"] = MagicMock()
sys.modules["Bio.SeqFeature"] = MagicMock()
sys.modules["Bio.SeqRecord"] = MagicMock()

from workflow.scripts.python.make_repeat_tables import make_repeats_df

def test_make_repeats_df(tmp_path):
    # Setup mock data for GRF output
    # Format according to parse_grf_output_no_ranges: '>record_id:rep1_start:rep2_end:lengthm'
    mock_grf_content = ">rec1:100:200:10m\n>rec1:150:250:15m\n"
    
    # We need to make sure the strain_index correctly extracts the strain from the path
    # Path: /dummy/strain_name/suffix/perfect.spacer.id
    # split("/") -> ['', 'dummy', 'strain_name', 'suffix', 'perfect.spacer.id']
    # index 2 would be 'strain_name'
    
    spacer_file = tmp_path / "strain1" / "suffix" / "perfect.spacer.id"
    spacer_file.parent.mkdir(parents=True)
    spacer_file.write_text(mock_grf_content)
    
    # Call the function
    # Using absolute path string to ensure split("/") works as expected on this system
    df = make_repeats_df(str(spacer_file), strain_index=-3)
    
    assert len(df) == 2
    assert df.loc[0, "record_id"] == "rec1"
    assert df.loc[0, "start_1"] == 100
    assert df.loc[0, "end_1"] == 109 # 100 + 10 - 1
    assert df.loc[0, "start_2"] == 191 # 200 - 10 + 1
    assert df.loc[0, "end_2"] == 200
    assert df.loc[0, "length"] == 10
    assert df.loc[0, "strain"] == "strain1"
    
    assert df.loc[1, "length"] == 15
    assert df.loc[1, "strain"] == "strain1"

def test_make_repeats_df_different_strain_index(tmp_path):
    # Test that strain_index parameter works
    mock_grf_content = ">rec1:100:200:10m\n"
    spacer_file = tmp_path / "some_strain" / "file.id"
    spacer_file.parent.mkdir(parents=True)
    spacer_file.write_text(mock_grf_content)
    
    # If path is .../some_strain/file.id, then .split("/")[-2] is some_strain
    df = make_repeats_df(str(spacer_file), strain_index=-2)
    assert df.loc[0, "strain"] == "some_strain"

def test_make_repeats_df_duplicates_handling(tmp_path):
    # make_repeats_df itself doesn't drop duplicates, the main block does
    # But we can verify it parses duplicates if they exist in the file
    mock_grf_content = ">rec1:100:200:10m\n>rec1:100:200:10m\n"
    spacer_file = tmp_path / "strain" / "file.id"
    spacer_file.parent.mkdir(parents=True)
    spacer_file.write_text(mock_grf_content)
    
    df = make_repeats_df(str(spacer_file), strain_index=-2)
    assert len(df) == 2
    
    # Now verify that drop_duplicates handles it correctly (as it does in main)
    df.drop_duplicates(inplace=True)
    assert len(df) == 1
