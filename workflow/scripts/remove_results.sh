#!/bin/bash
cd /home/andrei/Data/HeteroR/results

directories=( annotations assemblies assemblies_joined data_filtered direct_repeats mapping plasmids qualcheck_assembly qualcheck_reads resistance_genes )

print_help() {
   echo "this script removes directories in 'results' for a particular strain"
   echo "required argument: list of strains, one per line"
}


if [ "$#" -eq 0 ]
then
    print_help
    exit 1
fi

# 'line' is 'strain'
while IFS= read -r line; do
  for dir in "${directories[@]}"; do
    #echo "$dir/$line"
    rm -r "$dir/$line"
  done
  # rm final
  #echo "final/$line""_all.done"
  rm "final/$line""_all.done"
done < "$1"