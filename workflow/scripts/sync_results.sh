#!/bin/bash
cd /home/andrei/Data/HeteroR/results
directories=( assemblies_joined annotations direct_repeats mapping plasmids qualcheck_assembly qualcheck_reads )

while IFS= read -r line; do
  for dir in "${directories[@]}"; do
    sudo rsync -av "$dir/$line" "/home/andrei/Data/Argos/imb_sal_wrk/500_Sepsis_Eco/$dir"
  done
done < "$1"
