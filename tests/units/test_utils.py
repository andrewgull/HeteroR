import pytest
import pandas as pd
from unittest.mock import MagicMock, patch
import os
import sys
from io import StringIO

# Add the repository root to sys.path to allow importing from workflow.scripts
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

from workflow.scripts import utils

@pytest.fixture
def mock_snakemake_env():
    """
    Fixture to mock Snakemake globals (expand, strains, parents) in the utils module.
    """
    # Store original values if they exist (unlikely in this context but good practice)
    original_expand = getattr(utils, 'expand', None)
    original_strains = getattr(utils, 'strains', None)
    original_parents = getattr(utils, 'parents', None)
    
    # Mock expand function
    utils.expand = MagicMock()
    # Mock strains and parents data
    utils.strains = {"strains": ["S1", "S2"]}
    utils.parents = {"parents": ["P1"]}
    
    yield utils
    
    # Cleanup: restore or delete mocked attributes
    if original_expand:
        utils.expand = original_expand
    else:
        del utils.expand
        
    if original_strains:
        utils.strains = original_strains
    else:
        del utils.strains
        
    if original_parents:
        utils.parents = original_parents
    else:
        del utils.parents

# --- Tests for final_files ---

@pytest.mark.parametrize("workflow_type, expected_wildcard, expected_values", [
    ("assembly", "strain", ["S1", "S2"]),
    ("annotation", "strain", ["S1", "S2"]),
    ("mutants", "parent", ["P1"]),
    ("phylogeny", "strain", ["S1", "S2"]),
])
def test_final_files_success(mock_snakemake_env, workflow_type, expected_wildcard, expected_values):
    """
    Test successful expansion for valid workflow types in final_files.
    """
    utils.final_files(workflow_type)
    
    utils.expand.assert_called_with(
        f"results/final/{{{expected_wildcard}}}_{workflow_type}_all.done", 
        **{expected_wildcard: expected_values}
    )

def test_final_files_unknown_type(mock_snakemake_env):
    """
    Test that an unknown workflow type raises a ValueError in final_files.
    """
    with pytest.raises(ValueError, match="Unknown workflow type: unknown"):
        utils.final_files("unknown")

# --- Tests for get_sample_path ---

def test_get_sample_path_success():
    """
    Test that get_sample_path returns the correct path for a valid sample.
    """
    df = pd.DataFrame({"sample": ["A", "B"], "path": ["path/A", "path/B"]})
    wildcards = MagicMock()
    wildcards.sample = "A"
    
    result = utils.get_sample_path(wildcards, df)
    assert result == "path/A"

def test_get_sample_path_failure():
    """
    Test that get_sample_path raises ValueError for a missing sample.
    """
    df = pd.DataFrame({"sample": ["A", "B"], "path": ["path/A", "path/B"]})
    wildcards = MagicMock()
    wildcards.sample = "C"
    
    with pytest.raises(ValueError, match="Sample 'C' not found in the provided DataFrame"):
        utils.get_sample_path(wildcards, df)

# --- Tests for tsv2dict ---

def test_tsv2dict(tmp_path):
    """
    Test that tsv2dict correctly converts a TSV file to a dictionary.
    """
    d = tmp_path / "subdir"
    d.mkdir()
    p = d / "config.tsv"
    p.write_text("param\tvalue\nparam1\tval1\nparam2\tval2\n")
    
    result = utils.tsv2dict(str(p))
    assert result == {"param1": "val1", "param2": "val2"}
