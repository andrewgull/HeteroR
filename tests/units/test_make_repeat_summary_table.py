import pytest
import pandas as pd
import os
import sys
import numpy as np

# Add repository root to sys.path to allow importing from workflow.scripts.python
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

from workflow.scripts.python.make_repeat_summary_table import process_repeat_summary

def test_process_repeat_summary_basic(tmp_path):
    # Setup mock data for repeat summary
    # Using a small example where end_1=10, start_2=20
    # AR_length should be 20 - 10 + 1 = 11
    repeat_data = {
        "Unnamed: 0": [0],
        "record_id": ["rec1"],
        "end_1": [10],
        "start_2": [20]
    }
    repeat_df = pd.DataFrame(repeat_data)
    repeat_file = tmp_path / "repeat.csv"
    repeat_df.to_csv(repeat_file, index=False)

    # Setup mock data for BED file
    # Format: chrom, start, end, record_id
    # Using start=0, end=100
    # gene_center = (100 - 0 + 1) / 2 = 50.5
    bed_file = tmp_path / "bed.bed"
    with open(bed_file, "w") as f:
        f.write("chr1\t0\t100\trec1\n")

    # Run the function
    result_df = process_repeat_summary([str(repeat_file)], [str(bed_file)])

    # Assertions
    assert result_df.loc[0, "AR_length"] == 11
    assert result_df.loc[0, "gene_center"] == 50.5
    
    # spans_center logic: (end_1 <= gene_center) & (start_2 >= gene_center)
    # (10 <= 50.5) & (20 >= 50.5) -> True & False -> False
    assert result_df.loc[0, "spans_center"] == "no"
    
    # Check that Unnamed: 0 was dropped
    assert "Unnamed: 0" not in result_df.columns

def test_spans_center_yes(tmp_path):
    # Case where it should span center
    # gene_center = (100 - 0 + 1) / 2 = 50.5
    # end_1 = 40, start_2 = 60
    # (40 <= 50.5) & (60 >= 50.5) -> True & True -> True
    repeat_data = {
        "record_id": ["rec2"],
        "end_1": [40],
        "start_2": [60]
    }
    repeat_df = pd.DataFrame(repeat_data)
    repeat_file = tmp_path / "repeat2.csv"
    repeat_df.to_csv(repeat_file, index=False)

    bed_file = tmp_path / "bed2.bed"
    with open(bed_file, "w") as f:
        f.write("chr1\t0\t100\trec2\n")

    result_df = process_repeat_summary([str(repeat_file)], [str(bed_file)])
    
    assert result_df.loc[0, "spans_center"] == "yes"

def test_multiple_inputs(tmp_path):
    # Test concatenation of multiple repeat files
    repeat_df1 = pd.DataFrame({"record_id": ["rec1"], "end_1": [10], "start_2": [20]})
    repeat_df2 = pd.DataFrame({"record_id": ["rec2"], "end_1": [30], "start_2": [40]})
    r1 = tmp_path / "r1.csv"; repeat_df1.to_csv(r1, index=False)
    r2 = tmp_path / "r2.csv"; repeat_df2.to_csv(r2, index=False)

    bed_file = tmp_path / "bed.bed"
    with open(bed_file, "w") as f:
        f.write("chr1\t0\t100\trec1\n")
        f.write("chr1\t0\t100\trec2\n")

    result_df = process_repeat_summary([str(r1), str(r2)], [str(bed_file)])
    
    assert len(result_df) == 2
    assert set(result_df["record_id"]) == {"rec1", "rec2"}

def test_missing_record_id(tmp_path):
    # Test outer join behavior when a record_id is missing in one of the files
    repeat_df = pd.DataFrame({"record_id": ["rec1"], "end_1": [10], "start_2": [20]})
    repeat_file = tmp_path / "repeat.csv"
    repeat_df.to_csv(repeat_file, index=False)

    bed_file = tmp_path / "bed.bed"
    with open(bed_file, "w") as f:
        f.write("chr1\t0\t100\trec2\n") # mismatch

    result_df = process_repeat_summary([str(repeat_file)], [str(bed_file)])
    
    # Due to how="outer", we should have both rec1 and rec2
    assert len(result_df) == 2
    assert set(result_df["record_id"]) == {"rec1", "rec2"}
    
    # rec1 will have NaN for gene_center, rec2 will have NaN for end_1, start_2 etc.
    # We should check how np.where handles NaNs
    # (NaN <= NaN) is False in numpy normally? 
    # Actually np.where with NaNs in conditions:
    # (NaN <= values) or (values <= NaN) will be False.
    # So spans_center should be 'no' for both.
    
    rec1_row = result_df[result_df["record_id"] == "rec1"].iloc[0]
    rec2_row = result_df[result_df["record_id"] == "rec2"].iloc[0]
    
    assert rec1_row["spans_center"] == "no"
    assert rec2_row["spans_center"] == "no"
