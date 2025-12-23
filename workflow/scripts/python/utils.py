import pandas as pd
from typing import Any
from snakemake.io import expand


def get_sample_path(wildcards: Any, data_frame: pd.DataFrame) -> str:
    """
    Get the path to the fastq file using its sample name.
    arguments:
        wildcards: wildcards object.
        data_frame: A DataFrame object containing sample names and paths.
    returns:
        str: The path to the fastq file.
    """
    row = data_frame[data_frame["sample"] == wildcards.sample]
    if not row.empty:
        return row["path"].values[0]
    else:
        raise ValueError(f"Sample '{wildcards.sample}' not found in the provided DataFrame.")


def tsv2dict(file_path: str) -> dict:
    """
    Read a TSV file and convert it to a dictionary.
    arguments:
        file_path: Path to the TSV file.
    returns:
        dict: A dictionary with parameters as keys and their values.
    """
    df = pd.read_csv(
        file_path,
        sep="\t",
        dtype={"param": str, "value": str}
    )
    data_dict = df.set_index("param")["value"].to_dict()
    return data_dict


def get_strains(config: dict) -> list:
    """Get the list of strains from the config."""
    return pd.read_csv(config.get("strains", ""), dtype={"strains": str})["strains"].tolist()


def get_parents(config: dict) -> list:
    """Get the list of parents from the config."""
    return pd.read_csv(config.get("parents", ""), dtype={"parents": str})["parents"].tolist()


def final_files(config: dict, workflow_type: str) -> Any:
    """
    Get the final files for a given workflow type.
    arguments:
        workflow_type: The type of workflow.
    returns:
        Any: The final files for the given workflow type.
    """
    if workflow_type == "assembly":
        return expand(f"results/final/{{strain}}_{workflow_type}_all.done", strain=get_strains(config))
    elif workflow_type == "annotation":
        return expand(f"results/final/{{strain}}_{workflow_type}_all.done", strain=get_strains(config))
    elif workflow_type == "mutants":
        return expand(f"results/final/{{parent}}_{workflow_type}_all.done", parent=get_parents(config))
    elif workflow_type == "phylogeny":
        return expand(f"results/final/{{strain}}_{workflow_type}_all.done", strain=get_strains(config))
    else:
        raise ValueError(f"Unknown workflow type: {workflow_type}")


if __name__ == "__main__":
    print("hello, this is utils.py")