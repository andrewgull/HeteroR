DRAFT="/home/andrei/Data/HeteroR/results/assemblies/DA63960/pooled/flye/medaka/consensus.fasta";
R1="/home/andrei/Data/HeteroR/resources/data_raw/DA63960/Illumina/renamed/DA63960_1.fq.gz";
R2="/home/andrei/Data/HeteroR/resources/data_raw/DA63960/Illumina/renamed/DA63960_2.fq.gz";

# activate bwa-env
bwa index $DRAFT
bwa mem -t 10 -a $DRAFT $R1 > alignments_1.sam
bwa mem -t 10 -a $DRAFT $R2 > alignments_2.sam

# activate polypolish-env
polypolish_insert_filter.py --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
polypolish $DRAFT filtered_1.sam filtered_2.sam > assembly_polished.fasta

