import pytest
import pandas as pd
import os
import sys
from unittest.mock import MagicMock

# Add repository root to sys.path
repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, repo_root)

# Mock Bio components to avoid dependency issues in test environment
mock_bio = MagicMock()
sys.modules["Bio"] = mock_bio
sys.modules["Bio.SeqIO"] = MagicMock()
sys.modules["Bio.SeqFeature"] = MagicMock()
sys.modules["Bio.SeqRecord"] = MagicMock()

from workflow.scripts.python.rgi2gff import parse_rgi_output, process_gbk_records

def test_parse_rgi_output(tmp_path):
    # Setup mock RGI TSV data
    data = {
        "ORF_ID": ["gene1 [partial]", "gene2 [partial]", "gene3 [partial]"],
        "Col1": [0,0,0], "Col2": [0,0,0], "Col3": [0,0,0], "Col4": [0,0,0],
        "Cut_Off": ["Perfect", "Strict", "Loose"],
        "Col6": [0,0,0], "Col7": [0,0,0], "Col8": [0,0,0], "Col9": [0,0,0], "Col10": [0,0,0], "Col11": [0,0,0],
        "SNP": ["No_SNP", "S83L", "No_SNP"],
        "Col13": [0,0,0],
        "Drug Class": ["class1", "class2", "class3"],
        "Col15": [0,0,0],
        "AMR Gene Family": ["family1", "family2", "family3"]
    }
    # Ensure indices match: 
    # 0: ORF_ID, 5: Cut_Off, 12: SNP, 14: Drug Class, 16: AMR Gene Family
    # The dict above has 17 keys, so indices should work if ordered correctly.
    # However, row[1][idx] in iterrows() uses positional index if it's a Series.
    
    rgi_df = pd.DataFrame(data)
    rgi_file = tmp_path / "rgi.tsv"
    rgi_df.to_csv(rgi_file, sep="\t", index=False)
    
    ids, metadata = parse_rgi_output(str(rgi_file), "Loose")
    
    assert "gene1" in ids
    assert "gene2" in ids
    assert "gene3" not in ids
    
    assert metadata["gene1"] == ["No_SNP", "class1", "family1"]
    assert metadata["gene2"] == ["S83L", "class2", "family2"]

def test_process_gbk_records():
    # Setup mock SeqRecord and features
    # Since we mocked Bio at the module level, we can just use MagicMocks directly.
    
    def create_mock_feature(locus_tag):
        f = MagicMock()
        f.type = "gene"
        f.qualifiers = {"locus_tag": [locus_tag]}
        return f

    rec1 = MagicMock()
    rec1.id = "PREFIX_1"
    feat1 = create_mock_feature("gene1")
    feat2 = create_mock_feature("gene2")
    rec1.features = [feat1, feat2]
    
    rec2 = MagicMock()
    rec2.id = "PREFIX_2"
    feat3 = create_mock_feature("gene3")
    rec2.features = [feat3]
    
    records = [rec1, rec2]
    rgi_ids = ["gene1", "gene3"]
    rgi_dict = {
        "gene1": ["n/a", "class1", "family1"],
        "gene3": ["n/a", "class3", "family3"]
    }
    
    processed = process_gbk_records(records, rgi_ids, rgi_dict)
    
    # Assertions
    assert len(processed) == 2 # both rec1 and rec2 have at least one matching gene
    
    # rec1 should now only have gene1
    assert processed[0].id == "1"
    assert len(processed[0].features) == 1
    assert processed[0].features[0].qualifiers["locus_tag"][0] == "gene1"
    assert processed[0].features[0].qualifiers["SNP"] == "n/a"
    
    # rec2 should now only have gene3
    assert processed[1].id == "2"
    assert len(processed[1].features) == 1
    assert processed[1].features[0].qualifiers["locus_tag"][0] == "gene3"
    assert processed[1].features[0].qualifiers["AMR_gene_family"] == "family3"
