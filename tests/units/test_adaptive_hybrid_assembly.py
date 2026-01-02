import pytest
from unittest.mock import MagicMock, patch
import os
import sys

# Add repository root to sys.path
repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, repo_root)

# Import the functions to test
from workflow.scripts.python.adaptive_hybrid_assembly import calculate_coverage, main_logic, run_unicycler, run_fmp

def test_calculate_coverage():
    # Setup mock subprocess output for seqkit stats
    # Format: file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len
    # sum_len is index 4
    mock_run = MagicMock()
    mock_run.stdout = "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\nreads.fasta\tFASTA\tDNA\t100\t1000000\t1000\t10000\t50000\n"
    
    with patch("subprocess.run", return_value=mock_run) as mock_subprocess:
        # 1,000,000 / 50,000 = 20
        coverage, out = calculate_coverage("reads.fasta", 50000)
        assert coverage == 20
        mock_subprocess.assert_called_once_with("seqkit stats reads.fasta -T", shell=True, capture_output=True, text=True, check=True)

def test_main_logic_shallow_coverage():
    # Mocking functions inside main_logic to test the decision flow
    with patch("workflow.scripts.python.adaptive_hybrid_assembly.calculate_coverage") as mock_cov, \
         patch("workflow.scripts.python.adaptive_hybrid_assembly.run_unicycler") as mock_uni, \
         patch("workflow.scripts.python.adaptive_hybrid_assembly.run_fmp") as mock_fmp:
        
        mock_cov.return_value = (10, MagicMock()) # coverage 10
        mock_uni.return_value = (MagicMock(), "ok")
        
        main_logic("sr1", "sr2", "lr", "assembly", "draft", "polish", 8, "caller", "size", 5000000, 20)
        
        mock_uni.assert_called_once()
        mock_fmp.assert_not_called()

def test_main_logic_deep_coverage():
    with patch("workflow.scripts.python.adaptive_hybrid_assembly.calculate_coverage") as mock_cov, \
         patch("workflow.scripts.python.adaptive_hybrid_assembly.run_unicycler") as mock_uni, \
         patch("workflow.scripts.python.adaptive_hybrid_assembly.run_fmp") as mock_fmp:
        
        mock_cov.return_value = (50, MagicMock()) # coverage 50
        mock_fmp.return_value = [MagicMock()]
        
        main_logic("sr1", "sr2", "lr", "assembly", "draft", "polish", 8, "caller", "size", 5000000, 20)
        
        mock_uni.assert_not_called()
        mock_fmp.assert_called_once()

def test_main_logic_unicycler_fallback():
    # If Unicycler fails, it should run FMP
    with patch("workflow.scripts.python.adaptive_hybrid_assembly.calculate_coverage") as mock_cov, \
         patch("workflow.scripts.python.adaptive_hybrid_assembly.run_unicycler") as mock_uni, \
         patch("workflow.scripts.python.adaptive_hybrid_assembly.run_fmp") as mock_fmp:
        
        mock_cov.return_value = (10, MagicMock()) # coverage 10 (triggers Unicycler)
        mock_uni.return_value = (MagicMock(), "error")
        mock_fmp.return_value = [MagicMock()]
        
        main_logic("sr1", "sr2", "lr", "assembly", "draft", "polish", 8, "caller", "size", 5000000, 20)
        
        mock_uni.assert_called_once()
        mock_fmp.assert_called_once()

def test_run_unicycler(tmp_path):
    # Test directory creation when Unicycler succeeds
    assembly_dir = tmp_path / "assembly"
    draft_dir = tmp_path / "draft"
    polish_dir = tmp_path / "polish"
    
    mock_run = MagicMock()
    mock_run.stderr = ""
    
    with patch("subprocess.run", return_value=mock_run):
        out, status = run_unicycler("sr1", "sr2", "lr", 8, str(assembly_dir), str(draft_dir), str(polish_dir))
        
        assert status == "ok"
        assert draft_dir.exists()
        assert polish_dir.exists()

def test_run_fmp_side_effects(tmp_path):
    # Test directory creation and file operations in FMP
    assembly_dir = tmp_path / "assembly"; assembly_dir.mkdir()
    draft_dir = tmp_path / "draft"
    polish_dir = tmp_path / "polish"
    
    # Create a mock assembly.fasta to satisfy os.rename
    (assembly_dir / "assembly.fasta").write_text("dummy")
    
    mock_run = MagicMock()
    mock_run.stdout = "done"
    mock_run.stderr = ""
    
    with patch("subprocess.run", return_value=mock_run) as mock_run_call, \
         patch("os.symlink") as mock_symlink:
        
        # Define a side effect to create draft_dir when medaka is called
        def side_effect(cmd, **kwargs):
            if "medaka_consensus" in cmd:
                os.mkdir(draft_dir)
            return mock_run
        
        mock_run_call.side_effect = side_effect
        
        out = run_fmp("sr1", "sr2", "lr", 8, str(assembly_dir), str(draft_dir), str(polish_dir), "caller", "size", 50)
        
        assert len(out) == 7 # flye, medaka, bwa index, bwa mem 1, bwa mem 2, plp filter, plp polish
        assert draft_dir.exists()
        assert polish_dir.exists()
        assert (assembly_dir / "assembly_flye_raw.fasta").exists()
        mock_symlink.assert_called_once()
