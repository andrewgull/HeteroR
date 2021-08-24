#!/usr/bin/env bash

# a script to map reads back to assembly
DIR=$1
CPU=$2

bwa index "$DIR"/assembly/hybrid_uni/assembly.fasta;
bwa mem -t 12 "$DIR"/assembly/hybrid_uni/assembly.fasta "$DIR"/Illumina/"$DIR"_1.fq.gz "$DIR"/Illumina/"$DIR"_2.fq.gz > "$DIR"/assembly/hybrid_uni/assembly.sam;
samtools view -b "$DIR"/assembly/hybrid_uni/assembly.sam > "$DIR"/assembly/hybrid_uni/assembly.bam;
samtools sort -o "$DIR"/assembly/hybrid_uni/assembly.sorted.bam -O BAM -@ $CPU "$DIR"/assembly/hybrid_uni/assembly.bam;
samtools index "$DIR"/assembly/hybrid_uni/assembly.sorted.bam;
rm "$DIR"/assembly/hybrid_uni/assembly.sam;
rm "$DIR"/assembly/hybrid_uni/assembly.bam
