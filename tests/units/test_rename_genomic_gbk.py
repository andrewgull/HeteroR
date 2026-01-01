import pytest
from unittest.mock import MagicMock
import sys
import os

# Add repository root to sys.path
repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, repo_root)

# Mock Bio.SeqIO and Bio.SeqRecord before importing the script
mock_bio = MagicMock()
mock_seqio = MagicMock()
mock_seqrecord = MagicMock()

sys.modules["Bio"] = mock_bio
sys.modules["Bio.SeqIO"] = mock_seqio
sys.modules["Bio.SeqRecord"] = mock_seqrecord

from workflow.scripts.python.rename_genomic_gbk import rename_gbk_records

def test_rename_gbk_records():
    # Setup mock records
    rec1 = MagicMock()
    rec1.id = "LOOOHHPE_1"
    
    rec2 = MagicMock()
    rec2.id = "ABCDEFG_2"
    
    rec3 = MagicMock()
    rec3.id = "PLAIN_ID"
    
    records = [rec1, rec2, rec3]
    
    # Call the function
    renamed_records = rename_gbk_records(records)
    
    # Assertions
    assert renamed_records[0].id == "1"
    assert renamed_records[1].id == "2"
    assert renamed_records[2].id == "ID" # split("_")[-1] on "PLAIN_ID" is "ID"

def test_rename_gbk_records_multiple_underscores():
    # Case with multiple underscores
    rec = MagicMock()
    rec.id = "CONTIG_ID_123"
    
    renamed = rename_gbk_records([rec])
    assert renamed[0].id == "123"

def test_rename_gbk_records_empty():
    assert rename_gbk_records([]) == []
