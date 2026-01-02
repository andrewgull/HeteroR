##########################################################################
# Script for find over-covered regions based on the provided depth file
# find normalized depth per window of specified size
# make bed file with coords of windows and contig names
# make and save plots
# save some summary along with bed file
# depth.txt is made by samtools depth/coverage
#########################################################################

#### OPEN LOG ####
sink(snakemake@log[[1]])

#### LIBRARIES ####
library(data.table)
library(ggplot2)
library(purrr)
library(dplyr)
library(ggpubr)

### FUNCTIONS ####
# calculate raw depth in windows of selected size
# and normalize it
normalize_depth <- function(depth_dt, contig_name, window_size = 1000) {
  depth_contig <- depth_dt[V1 == contig_name, ][, window := 1 + (V2 - 1) %/% window_size][, coverage := mean(V3), by = window][, z := (coverage - mean(coverage)) / sd(coverage)]
  return(unique(depth_contig, by = c("V1", "window")))
}

# make a BED file from normalized depth
make_bed <- function(depth_dt, window_size) {
  bed <- depth_dt[, start := V2 - 1][, end := start + window_size][, V3 := NULL][, coverage := NULL][, V2 := NULL]
  # set names
  setnames(bed, "V1", "contig")
  # order columns
  setcolorder(bed, c("contig", "start", "end", "window", "z"))
  return(bed)
}

# make a plot of coverage per window
plot_coverage <- function(df, title, z_line) {
  ggplot(df, aes(window, z)) +
    geom_point(size = 0.2) +
    geom_hline(yintercept = z_line, color = "red") +
    ggtitle(title)
}

# get normalized depth per contig, filter them & convert to BED
depth2bed <- function(depth_dt, window_len = 1000, z_threshold = 2) {
  # available contig names
  contigs <- unique(depth_dt, by = "V1")[, V1]
  # get normalized depth for each contig
  depth_norm_list <- map(
    contigs,
    ~ normalize_depth(depth_dt, ., window_size = window_len)
  )
  # make list of coverage plots
  depth_plots <- map2(
    depth_norm_list,
    contigs,
    ~ plot_coverage(.x, .y, z_threshold)
  )
  # filter normalized depth
  depth_norm_filtered <- bind_rows(depth_norm_list)[z >= z_threshold, ]
  # turn the filtered depth into a bed file (keep in mind the 0-based indexing!)
  bed <- make_bed(depth_norm_filtered, window_size = window_len)
  return(list(bed, depth_plots))
}

# apply the functions above
main <- function(input_depth_file, win_len, z_score,
                 output_bed_file, output_plot_file) {
  # read depth file
  depth <- fread(input_depth_file, sep = "\t")
  beds_and_plots <- depth2bed(depth,
    window_len = win_len,
    z_threshold = z_score
  )
  # save BED files
  fwrite(beds_and_plots[[1]],
    output_bed_file,
    sep = "\t",
    col.names = FALSE
  )
  # arrange and save plots
  n_plots <- length(beds_and_plots[[2]])
  plots <- ggarrange(
    plotlist = beds_and_plots[[2]],
    ncol = 1,
    nrow = n_plots
  )
  ggsave(
    filename = output_plot_file,
    plot = plots, width = 14,
    height = 2 * n_plots
  )
}

#### RUN #####
main(
  input_depth_file = snakemake@input[["depth"]],
  win_len = snakemake@params[["w"]],
  z_score = snakemake@params[["z"]],
  output_bed_file = snakemake@output[["bed"]],
  output_plot_file = snakemake@output[["plot"]]
)

print("Script finished, no errors.")

#### CLOSE LOG ####
sink()
