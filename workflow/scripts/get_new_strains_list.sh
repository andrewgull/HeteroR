#!/bin/bash

print_help() {
   printf "This script takes strains available in ./resources,
compares them to a provided list of strains
and returns strains in the list but not in resources
Required arguments:
    list of strain available on Argos\n"
}

if [ "$#" -eq 0 ]
then
    print_help
    exit 1
fi

## VARIABLES ##
RESOURCES_PATH="/home/andrei/Data/HeteroR/resources/data_raw"
STRAINS_ALL_ARGOS=$1  # list of all strains available on Argos

# GET SORTED LIST OF STRAINS FROM RESOURCES
ls $RESOURCES_PATH | sort > strains_in_resources.tmp

# SORT PROVIDED STRAINS LIST
cat $STRAINS_ALL_ARGOS | sort > strains_on_argos.tmp

# FIND DIFFERENCE BETWEEN TWO LISTS
comm -13 strains_in_resources.tmp strains_on_argos.tmp
rm strains_in_resources.tmp
rm strains_on_argos.tmp
