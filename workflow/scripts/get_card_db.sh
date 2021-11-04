#!/usr/bin/env bash

address="https://card.mcmaster.ca/latest/data"
dest="."

conda activate rgi-env
wget $address &&
tar -xvf $dest"/data"  $dest"/card.json" &&
rgi load --card_json $dest"/card.json" --local &&
# clean up
rm  $dest"/data"
rm  $dest"/card.json"
echo "done"
