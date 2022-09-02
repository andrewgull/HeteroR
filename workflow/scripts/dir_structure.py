import glob
import re

# n = 435
Fastq = glob.glob("/home/andrei/Data/Argos/imb_sal_raw/Sequenced_reference_strains/Sequencing/Strains/DA*/Nanopore/Fastq",
                  recursive=True)

# n = 508
Nanopore = glob.glob("/home/andrei/Data/Argos/imb_sal_raw/Sequenced_reference_strains/Sequencing/Strains/DA*/Nanopore/",
                     recursive=True)


def get_strain_name(string):
    m = re.search('DA[0-9].*/Nanopore', string)
    return m.group(0)


nanopore_strains = list(map(get_strain_name, Nanopore))
fastq_strains = list(map(get_strain_name, Fastq))

# n = 73
not_fastq_strains = set(nanopore_strains) - set(fastq_strains)

# n = 60
fastq_pass = glob.glob("/home/andrei/Data/Argos/imb_sal_raw/Sequenced_reference_strains/Sequencing/Strains/DA*/Nanopore/fastq_pass",
                  recursive=True)

fastq_pass_strains = list(map(get_strain_name, fastq_pass))

set(not_fastq_strains) - set(fastq_pass_strains)

# in the end I have three types of names:
# Fastq, fastq_pass/barcode..., fast5q_pass/barcode...,
