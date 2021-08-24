#!/bin/bash

mkdir qualcheck
mkdir qualcheck/Illumina; mkdir qualcheck/Nanopore

# Exception in thread "Thread-1" java.lang.OutOfMemoryError: GC overhead limit exceeded
# makes fastqc unresponsive

fastqc -t 10 -q -o qualcheck/Illumina data_raw/DA*/Illumina/*gz
multiqc -f -o qualcheck/Illumina qualcheck/Illumina

fastqc -t 10 -q -o qualcheck/Nanopore data_raw/DA*/Nanopore/*gz
multiqc -f -o qualcheck/Nanopore qualcheck/Nanopore
