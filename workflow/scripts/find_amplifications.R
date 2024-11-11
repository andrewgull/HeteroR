##########################################################################
#
# A script for find over-covered regions based on the provided depth file
# find normalized depth per window of specified size
# make bed file with coords of windows and contig names
# make and save plots
# save some summary along with bed file
# depth.txt is made by samtools depth/coverage
#
#########################################################################

library(optparse)

# CLI parsing
option_list <- list(
    make_option(c("-i", "--input"),
                type = "character",
                default = NULL,
                help = "TSV file with per position sequence depth",
                metavar = "path"),
    make_option(c("-b", "--output_bed"),
                type = "character",
                default = NULL,
                help = "BED output file with coordinates of over-covered regions",
                metavar = "path"),
    make_option(c("-l", "--output_plot"),
                type = "character",
                default = NULL,
                help = "PNG output file with coverage plots",
                metavar = "path"),
    make_option(c("-z", "--threshold"),
                type = "integer",
                default = 2,
                help = "threshold to filter candidate regions (in SD units)",
                metavar = "int"
    ),
    make_option(c("-w", "--window"),
                type = "integer",
                default = 1000,
                help = "Window size to calculate coverage",
                metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)){
 print_help(opt_parser)
 stop("Input file must be provided", call. = FALSE)
}
if (is.null(opt$output_bed)){
 print_help(opt_parser)
 stop("Output file (BED) must be provided", call. = FALSE)
}
if (is.null(opt$output_plot)){
 print_help(opt_parser)
 stop("Output file (PNG) must be provided", call. = FALSE)
}

input_depth_file <- opt$input
output_bed_file <- opt$output_bed
output_plot_file <- opt$output_plot
z_threshold <- opt$threshold
win_len <- opt$window

#### LIBRARIES ####
library(data.table)
library(ggplot2)
library(purrr)
library(dplyr)
library(ggpubr)



# function to calculate raw depth in windows of selected size and normalize it
normalize_depth <- function(depth_dt, contig_name, window_size = 1000) {
  depth_contig <- depth_dt[V1 == contig_name,
  ][, window := 1 + (V2 - 1) %/% window_size
  ][, coverage := mean(V3), by = window
  ][, z := (coverage - mean(coverage)) / sd(coverage)]
  return(unique(depth_contig, by = c("V1", "window")))
}

# function to make BED file from normalized depth
make_bed <- function(depth_dt, window_size) {
  bed <- depth_dt[, start := V2 - 1
  ][, end := start + window_size
  ][, V3 := NULL
  ][, coverage := NULL
  ][, V2 := NULL]
  # set names
  setnames(bed, "V1", "contig")
  # order columns
  setcolorder(bed, c("contig", "start", "end", "window", "z"))
  return(bed)
}

# function to make a plot of coverage per window
plot_coverage <- function(df, title, z_line) {
  ggplot(df, aes(window, z)) +
  geom_point(size = 0.2) +
  geom_hline(yintercept = z_line, color = "red") +
  ggtitle(title)
}

# main function
main <- function(depth_dt, window_length = 1000, z_threshold = 2) {
  # available contig names
  contigs <- unique(depth_dt, by = "V1")[, V1]
  # get normalized depth for each contig
  depth_norm_list <- map(contigs, ~ normalize_depth(depth, ., window_size = window_length))
  # make list plots
  depth_plots <- map2(depth_norm_list, contigs, ~ plot_coverage(.x, .y, z_threshold))
  # filter normalized depth
  depth_norm_filtered <- bind_rows(depth_norm_list)[z >= z_threshold, ]
  # turn the filtered depth into a bed file (keep in mind the 0-based indexing!)
  bed <- make_bed(depth_norm_filtered, window_size = window_length)

  return(list(bed, depth_plots))
}

# Apply the functions above
depth <- fread(input_depth_file, sep = "\t")
bed_and_plots <- main(depth, window_length = win_len, z_threshold = z_threshold)
# write files on the disc
fwrite(bed_and_plots[[1]], output_bed_file, sep = "\t", col.names = FALSE)
# arrange plots
n_plots <- length(bed_and_plots[[2]])
plots <- ggarrange(plotlist = bed_and_plots[[2]], ncol = 1, nrow = n_plots)
ggsave(filename = output_plot_file, plot = plots, width = 14, height=2*n_plots)
print("Script finished, no errors.")
