#!/usr/bin/env Rscript
library(fastqcr)
library(stringr)
library(ggplot2)

keyword <- function(file_out){
  # this function returns keyword depending on the provided utput file path
  if (grepl("Illumina", file_out) & !grepl("trimmed", file_out)){
    key_word <- "Illumina"
    }
  else if (grepl("Illumina_trimmed", file_out)){
    key_word <- "Illumina_trimmed"
  }
  else if (grepl("Nanopore", file_out) & !grepl("filtered", file_out)){
    key_word <- "Nanopore"
  }
  else if (grepl("Nanopore_filtered", file_out)){
    key_word <- "Nanopore_filtered"
  }
  return(key_word)
}

quality_reports <- function(file_in, file_out, fastqc_exe, cpus){
  # TODO: check beahaviour with the new 'renamed' subdir in Illumina data
  # run fastqc, write summary table and length distribution for Nanopore reads
  # ARGUMENT PARSING SECTION
  key_word <- keyword(file_out = file_out)
  if (key_word == "Illumina"){
    # file_in: "data_raw/{strain}/Illumina/{strain}_1.fq.gz"
    input_dir <- paste0(str_replace(file_in, "Illumina/DA.*_1.fq.gz", ""), key_word)  # "data_raw/DA62920/" + Illumina
    output_dir <- paste0(str_replace(file_out, "Illumina/DA.*_1.fastqc.zip", ""), key_word)  # "qualcheck/DA62920/" + Illumina
  } else if (key_word == "Nanopore"){
    input_dir <- paste0(str_replace(file_in, "Nanopore/DA.*_all.fastq.gz", ""), key_word)
    output_dir <- paste0(str_replace(file_out, "Nanopore/DA.*_all.fastqc.zip", ""), key_word)
  } else if (key_word == "Illumina_trimmed"){
    # input: "data_filtered/{strain}/Illumina/{strain}_1.fq.gz",
    # output: "qualcheck/{strain}/Illumina_trimmed/{strain}_1.fastqc.zip"
    input_dir <- paste0(str_replace(file_in, "Illumina/DA.*_1.fq.gz", ""), "Illumina")
    output_dir <- paste0(str_replace(file_out, "Illumina_trimmed/DA.*_1.fastqc.zip", ""), key_word)  # this creates new dir
  } else if (key_word == "Nanopore_filtered"){
    # input: "data_filtered/{strain}/Nanopore/{strain}_all.fastq.gz"
    # output: "qualcheck/{strain}/Nanopore/{strain}_all.fastqc.zip"
    input_dir <- paste0(str_replace(file_in, "Nanopore/DA.*_all.fastq.gz"), "Nanopore")
    output_dir <- paste0(str_replace(file_out, "Nanopore_filtered/DA.*_all.fastqc.zip"), key_word)  # this creates new dir
  }

  # FASTQC RUN
  fastqc(fq.dir=input_dir, qc.dir=output_dir, fastqc.path = fastqc_exe, threads = cpus)

  # QC SUMMARY
  qc <- qc_aggregate(output_dir)
  qc_summary <- summary(qc)
  samples_warn_failed <- qc_problems(qc)

  # WRITING TO FILES
  strain <- str_extract(file_in, "DA[0-9]*")
  summary_name <- paste0("/", strain, "_summary.tsv")
  warning_name <- paste0("/", strain, "_warnings.tsv")
  write.table(qc_summary, paste0(output_dir, summary_name), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(samples_warn_failed, paste0(output_dir, warning_name), sep="\t", row.names = FALSE, quote = FALSE)

  # PLOTS FOR NANOPORE
  if (grepl("Nanopore", key_word)){
    # collect fastqc zip files
    nanopore_zip <- dir(output_dir, pattern = "*_all_.*.zip")
      # speed is not an issue here, so use a for loop
      for (filename in nanopore_zip){
        qc_results <- qc_read(paste0(output_dir, "/", filename))
        #qc_results <- qc_read(file_out)
        length_plot <- ggplot(data = qc_results$sequence_length_distribution, aes(Length, Count)) +
          geom_bar(stat="identity") +
          theme(axis.text.x = element_text(angle=45, size=8, vjust=0.5)) +
          ggtitle(filename)
        # save the plot, cut off '_fastqc.zip' from the filename
        out_filename <- str_replace(filename, "_fastqc.zip", "")
        png(paste0(output_dir, "/", out_filename, "_length_distribution.png"))
        print(length_plot)
        dev.off()
      }
  }
}

quality_reports(file_in = snakemake@input[[1]],
                    file_out = snakemake@output[[1]],
                    fastqc_exe = "/home/andrei/miniconda3/bin/fastqc",
                    cpus = 4)