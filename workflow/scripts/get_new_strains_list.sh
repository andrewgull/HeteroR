#!/bin/bash

# sort lists first
list_old_sorted=$1
list_new_sorted=$2

# update sample list
comm -13 $1 $2 > strains_to_add.txt
cat list_old_sorted strains_to_add > strains_all.txt
