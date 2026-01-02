import pytest
from unittest.mock import MagicMock, patch
import sys
import os

# Mock Biopython before any imports that might use it
mock_bio = MagicMock()
mock_seqio = MagicMock()
mock_seqrecord_mod = MagicMock()
mock_seq = MagicMock()

sys.modules["Bio"] = mock_bio
sys.modules["Bio.SeqIO"] = mock_seqio
sys.modules["Bio.SeqRecord"] = mock_seqrecord_mod
sys.modules["Bio.Seq"] = mock_seq

# Mock classes
mock_seq_class = MagicMock()
mock_seq.Seq = mock_seq_class
mock_seqrecord_class = MagicMock()
mock_seqrecord_mod.SeqRecord = mock_seqrecord_class

# Add repository root to sys.path to allow importing from workflow.scripts.python
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

def test_join_ends_function():
    # Setup test data
    from workflow.scripts.python.join_ends import join_ends
    
    # Configure mock SeqRecord to return identifiable objects
    def seq_record_side_effect(seq=None, id=None, name=None, description=None):
        m = MagicMock()
        m.seq = seq
        m.id = id
        return m
    mock_seqrecord_class.side_effect = seq_record_side_effect
    
    side_rec = MagicMock()
    side_rec.seq = "AAAA"
    side_rec.id = "rec1"
    
    normal_rec = MagicMock()
    normal_rec.seq = "TTTT"
    normal_rec.id = "rec1"
    
    side_dict = {"rec1": side_rec}
    normal_dict = {"rec1": normal_rec}
    
    # Test left=True (side + normal)
    joined_left = join_ends(side_dict, normal_dict, left=True)
    assert len(joined_left) == 1
    # Check that addition happened: side.seq + normal.seq
    # Since we can't easily mock __add__ on a string unless it's a mock
    # Let's ensure seq is a mock that supports addition
    side_rec.seq = MagicMock()
    normal_rec.seq = MagicMock()
    side_rec.seq.__add__.return_value = "AAAATTTT"
    
    joined_left = join_ends(side_dict, normal_dict, left=True)
    assert joined_left[0].seq == "AAAATTTT"
    side_rec.seq.__add__.assert_called_with(normal_rec.seq)

def test_main_logic_case1():
    """CASE 1: normal is non-empty, left and right are empty"""
    from workflow.scripts.python import join_ends
    
    # Mock snakemake inside the module
    join_ends.snakemake = MagicMock()
    join_ends.snakemake.input = {"within": "n.fasta", "five_end": "l.fasta", "three_end": "r.fasta"}
    join_ends.snakemake.output = {"within": "out_n.fasta", "five_end": "out_l.fasta", "three_end": "out_r.fasta"}
    join_ends.snakemake.log = ["test.log"]
    
    rec_n = MagicMock(id="n1")
    
    with patch("workflow.scripts.python.join_ends.SeqIO") as mock_seqio, \
         patch("workflow.scripts.python.join_ends.open", create=True) as mock_open:
        
        mock_seqio.to_dict.side_effect = [
            {"n1": rec_n}, # normal_records
            {},            # left_records
            {}             # right_records
        ]
        
        join_ends.run_snakemake()
                
        assert mock_seqio.write.call_count == 3
        args, _ = mock_seqio.write.call_args_list[0]
        assert list(args[0])[0].id == "n1"
        assert args[1] == "out_n.fasta"

def test_main_logic_case2():
    """CASE 2: normal and left are non-empty, right is empty"""
    from workflow.scripts.python import join_ends
    
    join_ends.snakemake = MagicMock()
    join_ends.snakemake.input = {"within": "n.fasta", "five_end": "l.fasta", "three_end": "r.fasta"}
    join_ends.snakemake.output = {"within": "out_n.fasta", "five_end": "out_l.fasta", "three_end": "out_r.fasta"}
    join_ends.snakemake.log = ["test.log"]
    
    rec_n = MagicMock(id="key1")
    rec_n.seq = MagicMock()
    rec_l = MagicMock(id="key1")
    rec_l.seq = MagicMock()
    rec_l.seq.__add__.return_value = "AAAATTTT"
    
    with patch("workflow.scripts.python.join_ends.SeqIO") as mock_seqio, \
         patch("workflow.scripts.python.join_ends.open", create=True) as mock_open:
        
        mock_seqio.to_dict.side_effect = [
            {"key1": rec_n, "key2": MagicMock(id="key2")}, # normal_records
            {"key1": rec_l}, # left_records
            {}  # right_records
        ]
        
        join_ends.run_snakemake()
                
        # Check output_5_end (second call)
        args, _ = mock_seqio.write.call_args_list[1]
        joined_recs = args[0]
        assert len(joined_recs) == 1
        assert joined_recs[0].seq == "AAAATTTT"
        
        # Check output_normal (first call) - key1 should be popped
        args_normal, _ = mock_seqio.write.call_args_list[0]
        normal_recs_out = list(args_normal[0])
        assert len(normal_recs_out) == 1
        assert normal_recs_out[0].id == "key2"
