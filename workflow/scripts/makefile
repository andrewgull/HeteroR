# make file for bash scripts
SHELL=/bin/bash

#### retrieve particular lines from the full strain list ####
### line X to Y
sed -n '11,20p' strain_list_full.txt > strain_list11-20.txt


#### file preparation/processing ####
join_nano:
	for S in DA*; do zcat $S/Nanopore/"$S"_*.fastq.gz > $S/Nanopore/"$S"_total.fastq && pigz -p 10 "$S"_total.fastq; done

#### quality control ####
fastqc:
	fastqc -t 10 -q -o qualcheck/Illumina DA*/Illumina/*gz
	fastqc -t 10 -q -o qualcheck/Nanopore DA*/Nanopore/*gz

multiqc:
	multiqc -f -o qualcheck/Illumina qualcheck/Illumina
	multiqc -f -o qualcheck/Nanopore qualcheck/Nanopore

#### AUTOMATED QC ####
# default filters:
# quality filter: enabled
# length filter: enabled
# low complexity filter: disabled
# adapter trimming: enabled
# overlap base correction: disabled (use --correction to enable it)
# polyG tail trimming: enabled for NextSeq/NoveSeq data
# polyX trimming: disabled
fastp:
  fastp --in1 file_1.fq.gz --in2 file_2.fq.gz --out1 file1_1_fastp_out.fq.gz  \
	--out2 file1_2_fastp_out.fq.gz --thread 6 --html report.html --qualified_quality_phred 20 \
	--cut_window_size 4 --cut_right 20 --length_required 50

#### length QC ####
soft_qc_nano:
		filtlong --min_length 1000 --keep_percent 95 DAXXX_total.fastq.gz > DAXXX_total_filt.fastq;
		pigz -p 10 DAXXX_total_filt.fastq

#### coverage calculations ####
coverage_nano:
	./scripts/get_coverage.sh Nanopore

coverage_illu:
	./scripts/get_coverage.sh Illumina

map_back:
	for D in DA*; do ./scripts/map_back.sh $D 12; done

backmapped_coverage:
	bedtools genomecov -ibam assembly.sorted.bam -bga > assembly.coverage.bga
	grep -w "0$" assembly.coverage.bga
