# script for find over-covered regions based on the provided depth file
library(optparse)

# CLI parsing
option_list <- list(
    make_option(c("-i", "--input"),
                type = "character",
                default = NULL,
                help = "TSV file with per position sequence depth",
                metavar = "path"),
    make_option(c("-o", "--output"),
                type = "character",
                default = NULL,
                help = "BED output file with coordinates of over-covered regions",
                metavar = "path")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)){
 print_help(opt_parser)
 stop("Input file must be provided", call. = FALSE)
}
if (is.null(opt$output)){
 print_help(opt_parser)
 stop("Output file must be provided", call. = FALSE)
}

input_depth_file <- opt$input
output_bed_file <- opt$output

library(data.table)
library(ggplot2)

# find normalized depth per window of specified size
# make bed file with coords of windows and contig names
# make and save plots
# save some summary along with bed file
# depth.txt is made by samtools depth/coverage

# function to calculate raw depth in windows of selected size and normalize it
normalize_depth <- function(depth_dt, contig_name, window_size=1000){
 depth_contig <- depth_dt[V1 == contig_name,
 ][, window := 1 + (V2 - 1) %/% window_size
 ][, coverage := mean(V3), by = window
 ][, z := (coverage - mean(coverage)) / sd(coverage)]
 return(unique(depth_contig, by = c("V1", "window")))
}

# function to make BED file from normalized depth
make_bed <- function(depth_dt, window_size){
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

main <- function(depth_dt, window_length=1000, z_threshold=2){
 # available contig names
 contigs <- unique(depth_dt, by = "V1")[, V1]
 # get normalized depth for each contig
 depth_norm <- purrr::map_dfr(contigs, ~ normalize_depth(depth, ., window_size = window_length))
 # filter normalized depth
 depth_norm_filtered <- depth_norm[z >= z_threshold,]
 # turn the filtered depth into a bed file (keep in mind the 0-based indexing!)
 bed <- make_bed(depth_norm_filtered, window_size = window_length)

 return(bed)
}

# Apply the functions above
depth <- fread(input_depth_file, sep = "\t")
candidate_regions <- main(depth)
fwrite(candidate_regions, output_bed_file, sep = "\t", col.names = FALSE)

