#!/bin/bash
cd /home/andrei/Data/HeteroR/results
directories=( assemblies assemblies_joined annotations direct_repeats mapping plasmids qualcheck_assembly qualcheck_reads resistance_genes )

print_help() {
   echo "required argument: list of strains, one per line"
}

if [ "$#" -eq 0 ]
then
    print_help
    exit 1
fi


while IFS= read -r line; do
  for dir in "${directories[@]}"; do
    sudo rsync -av "$dir/$line" "/home/andrei/Data/Argos/imb_sal_wrk/500_Sepsis_Eco/$dir"
  done
done < "$1"
