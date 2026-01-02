import pandas as pd
import io
import os
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

# Mock dependencies before imports
mock_bio = MagicMock()
mock_seqio = MagicMock()
mock_seqrecord_mod = MagicMock()
mock_seqfeature_class = MagicMock()
mock_featurelocation_class = MagicMock()
mock_exactposition_class = MagicMock()

# Configure Bio.SeqFeature module
mock_seqfeature_mod = MagicMock()
mock_seqfeature_mod.SeqFeature = mock_seqfeature_class
mock_seqfeature_mod.FeatureLocation = mock_featurelocation_class
mock_seqfeature_mod.ExactPosition = mock_exactposition_class

mock_bcbio = MagicMock()
mock_gff = MagicMock()

sys.modules["Bio"] = mock_bio
sys.modules["Bio.SeqIO"] = mock_seqio
sys.modules["Bio.SeqRecord"] = mock_seqrecord_mod
sys.modules["Bio.SeqFeature"] = mock_seqfeature_mod
sys.modules["BCBio"] = mock_bcbio
sys.modules["BCBio.GFF"] = mock_gff

# Allow Bio to be treated as a package for nested imports
mock_bio.__path__ = []

# Add repository root to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

from workflow.scripts.python import GRF_parser

def test_parse_grf_output_w_ranges():
    """Test parsing headers with ranges (legacy version)."""
    header = ">1:0-1000:100:900:10m"
    result = GRF_parser.parse_grf_output_w_ranges(header)
    
    # expected: [record_id, rep1_start, rep1_end, rep2_start, rep2_end, repeat_len]
    # record_id: 1
    # range_start: 0
    # repeat_len: 10
    # rep1_start_in_range: 100 => rep1_start_in_chrom: 0 + 100 = 100
    # rep1_end_in_chrom: 100 + 10 = 110
    # rep2_end_in_range: 900 => rep2_end_in_chrom: 0 + 900 = 900
    # rep2_start_in_chrom: 900 - 10 = 890
    
    assert result == ["1", 100, 110, 890, 900, 10]

def test_parse_grf_output_no_ranges():
    """Test parsing headers without ranges (newer version)."""
    header = ">contig_1_gene:100:900:20m"
    result = GRF_parser.parse_grf_output_no_ranges(header)
    
    # expected: [record_id, rep1_start, rep1_end, rep2_start, rep2_end, repeat_len]
    # record_id: contig_1_gene
    # repeat_len: 20
    # rep1_start: 100
    # rep1_end: 100 + 20 - 1 = 119
    # rep2_end: 900
    # rep2_start: 900 - 20 + 1 = 881
    
    assert result == ["contig_1_gene", 100, 119, 881, 900, 20]

def test_create_gff_record_with_features():
    """Test adding features to a SeqRecord."""
    features_df = pd.DataFrame({
        "record_id": ["contig_1_gene", "contig_1_gene"],
        "start_1": [100, 300],
        "end_1": [120, 320],
        "start_2": [800, 1000],
        "end_2": [820, 1020],
        "length": [20, 20]
    })
    
    mock_record = MagicMock()
    mock_record.id = "contig_1"
    mock_record.features = []
    
    # Configure SeqFeature class to return a mock with an ID
    def side_effect(*args, **kwargs):
        m = MagicMock()
        m.id = kwargs.get("id")
        return m
    mock_seqfeature_class.side_effect = side_effect
    
    result_record = GRF_parser.create_gff_record_with_features(features_df, mock_record)
    
    # Each row adds 2 features (pair of repeats)
    assert len(result_record.features) == 4
    assert result_record.features[0].id == "1"
    assert result_record.features[1].id == "1"
    assert result_record.features[2].id == "2"
    assert result_record.features[3].id == "2"

@patch("workflow.scripts.python.GRF_parser.SeqIO.parse")
def test_main_process(mock_seqio_parse):
    """Test the main procedural flow."""
    # Mock inputs
    input_assembly = Path("dummy.fasta")
    input_grf = Path("dummy.id")
    min_len = 10
    
    # Mock SeqIO records
    mock_rec = MagicMock()
    mock_rec.id = "contig_1"
    mock_seqio_parse.return_value = [mock_rec]
    
    # Mock GRF file content
    grf_content = ">contig_1_gene:100:900:20m\n"
    
    with patch("pathlib.Path.exists", return_value=True):
        with patch("builtins.open", return_value=io.StringIO(grf_content)):
            with patch("workflow.scripts.python.GRF_parser.create_gff_record_with_features") as mock_gff_func:
                mock_rec.features = [MagicMock()] # Ensure it has features to be returned
                mock_gff_func.return_value = mock_rec
                
                result = GRF_parser.main_process(input_assembly, input_grf, min_len)
                
                assert len(result) == 1
                assert result[0] == mock_rec
                mock_gff_func.assert_called_once()
