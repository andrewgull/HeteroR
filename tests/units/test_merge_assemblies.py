import os
import sys
from unittest.mock import MagicMock, patch

# Add repository root to sys.path
repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, repo_root)

# Mock Bio components
mock_bio = MagicMock()
sys.modules["Bio"] = mock_bio
mock_seqio = MagicMock()
sys.modules["Bio.SeqIO"] = mock_seqio
mock_seqrecord = MagicMock()
sys.modules["Bio.SeqRecord"] = mock_seqrecord

# Import the functions to test
from workflow.scripts.python.merge_assemblies import joiner

def test_joiner_both_files_exist(tmp_path):
    file1 = "file1.fasta"
    file2 = "file2.fasta"
    
    rec1 = MagicMock()
    rec1.id = "rec1"
    rec2 = MagicMock()
    rec2.id = "rec2"
    
    # Mocking os.path.isfile and os.path.getsize
    with patch("os.path.isfile", return_value=True), \
         patch("os.path.getsize", return_value=100), \
         patch("workflow.scripts.python.merge_assemblies.SeqIO.parse") as mock_parse:
        
        def side_effect(filename, format):
            if filename == file1:
                return [rec1]
            if filename == file2:
                return [rec2]
            return []
        
        mock_parse.side_effect = side_effect
        
        joined = joiner(file1, file2)
        
        assert len(joined) == 2
        assert joined[0].id == "rec1"
        assert joined[1].id == "rec2"

def test_joiner_only_first_exists(tmp_path):
    file1 = "file1.fasta"
    file2 = "non_existent.fasta"
    
    rec1 = MagicMock()
    rec1.id = "rec1"
    
    with patch("os.path.isfile") as mock_isfile, \
         patch("os.path.getsize", return_value=100), \
         patch("workflow.scripts.python.merge_assemblies.SeqIO.parse") as mock_parse:
        
        def isfile_side_effect(filename):
            return filename == file1
        
        mock_isfile.side_effect = isfile_side_effect
        mock_parse.return_value = [rec1]
        
        joined = joiner(file1, file2)
        
        assert len(joined) == 1
        assert joined[0].id == "rec1"

def test_joiner_second_is_empty(tmp_path):
    file1 = "file1.fasta"
    file2 = "empty.fasta"
    
    rec1 = MagicMock()
    rec1.id = "rec1"
    
    with patch("os.path.isfile", return_value=True), \
         patch("os.path.getsize") as mock_getsize, \
         patch("workflow.scripts.python.merge_assemblies.SeqIO.parse") as mock_parse:
        
        def getsize_side_effect(filename):
            if filename == file1: return 100
            return 0
        
        mock_getsize.side_effect = getsize_side_effect
        mock_parse.return_value = [rec1]
        
        joined = joiner(file1, file2)
        
        assert len(joined) == 1
        assert joined[0].id == "rec1"

def test_joiner_second_is_none(tmp_path):
    file1 = "file1.fasta"
    
    rec1 = MagicMock()
    rec1.id = "rec1"
    
    with patch("os.path.isfile", return_value=True), \
         patch("os.path.getsize", return_value=100), \
         patch("workflow.scripts.python.merge_assemblies.SeqIO.parse") as mock_parse:
        
        mock_parse.return_value = [rec1]
        
        joined = joiner(file1, None)
        
        assert len(joined) == 1
        assert joined[0].id == "rec1"
