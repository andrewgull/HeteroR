#!/bin/bash
MACHINE=$1
GENOME_SIZE="5131220"
for strain in data_raw/DA*; do
  # get total length
  L=`seqkit stats $strain/$MACHINE/*.gz -T | cut -f 5 | sed '1d' | awk '{s+=$1}END{print s}'`
  # get mean coverage
  let "C = $L / $GENOME_SIZE"
  echo "$strain - $C"
done
