N="DA63960_all_pooled_filt.fastq.gz"
ASS_DIR="flye"
POL_DIR="medaka"
R1="/home/andrei/Data/HeteroR/resources/data_raw/DA63960/Illumina/renamed/DA63960_1.fq.gz";
R2="/home/andrei/Data/HeteroR/resources/data_raw/DA63960/Illumina/renamed/DA63960_2.fq.gz";

# assemble with flye
flye --nano-raw $N --threads 8 --out-dir $ASS_DIR -g 5m --asm-coverage 50

# polish with medaka
medaka_consensus -i $N -d $ASS_DIR'/assembly.fasta' -o $POL_DIR -t 8 -m r941_min_fast_g507

DRAFT=$POL_DIR'/consensus.fasta';


# map Illumina on medaka consensus
bwa index $DRAFT
bwa mem -t 10 -a $DRAFT $R1 > alignments_1.sam
bwa mem -t 10 -a $DRAFT $R2 > alignments_2.sam

# polish medaka consensus
polypolish_insert_filter.py --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
polypolish $DRAFT filtered_1.sam filtered_2.sam > assembly_polished.fasta

