#!/bin/bash

print_help() {
   printf "required arguments:\n sorted list of strains available in resources\n sorted list of strain available on Argos"
}

if [ "$#" -eq 0 ]
then
    print_help
    exit 1
fi

# sort lists first
list_old_sorted=$1  # list of strains currently in ./resources
list_new_sorted=$2  # list of all strains available on Argos

# update sample list
comm -13 $1 $2 > strains_to_add.txt  # find difference between two lists
cat list_old_sorted strains_to_add > strains_all.txt
