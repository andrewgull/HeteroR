#!/bin/bash
cd /home/andrei/Data/HeteroR/results

directories=( assemblies assemblies_joined annotations direct_repeats mapping plasmids qualcheck_assembly qualcheck_reads resistance_genes )

print_help() {
   echo "this script removes directories in 'results' for a particular strain"
   echo "required argument: list of strains, one per line"
}