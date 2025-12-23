import os

# --- Config Validation ---

def validate_config(config):
    """
    Validates the Snakemake configuration to ensure all required parameters
    and files are present.
    """
    
    # Check for mutual exclusivity of workflows
    run_flags = ["run_assembly", "run_annotation", "run_mutants", "run_phylogeny"]
    active_workflows = [flag for flag in run_flags if config.get(flag, False)]
    
    if len(active_workflows) > 1:
        raise ValueError(f"Config error: Multiple workflows enabled ({', '.join(active_workflows)}). Only one workflow can be active at a time.")

    # Check for required paths if specific parts of the pipeline are enabled
    
    # Assembly, Annotation, and Phylogeny workflow requirements
    if config.get("run_assembly", False) or config.get("run_annotation", False) or config.get("run_phylogeny", False):
        if "strains" not in config:
            raise ValueError("Config error: 'strains' path is required for enabled workflows (assembly/annotation/phylogeny).")
        if not os.path.exists(config["strains"]):
            raise FileNotFoundError(f"Config error: 'strains' file not found at '{config['strains']}'.")

    # Assembly and Annotation requirements
    if config.get("run_assembly", False) or config.get("run_annotation", False):
        if "reads_path" not in config:
             raise ValueError("Config error: 'reads_path' is required for enabled workflows (assembly/annotation).")

    # Mutant workflow requirements
    if config.get("run_mutants", False):
        if "parents" not in config:
            raise ValueError("Config error: 'parents' path is required for mutant workflow.")
        if not os.path.exists(config["parents"]):
            raise FileNotFoundError(f"Config error: 'parents' file not found at '{config['parents']}'.")
        if "mutants_path" not in config:
            raise ValueError("Config error: 'mutants_path' is required for mutant workflow.")
        if "parents_path" not in config:
            raise ValueError("Config error: 'parents_path' is required for mutant workflow.")

# Execute validation
validate_config(config)
