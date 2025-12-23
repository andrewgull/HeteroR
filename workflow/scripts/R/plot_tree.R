###########################################################
# script to plot and save phylogenetic tree
# based on core genome alignment
# input: newick tree, outgroup taxa name, HR/non-HR labels
# the exact look of the tree is hard-coded
###########################################################

#### OPEN LOG ####
sink(snakemake@log[[1]])

#### PACKAGES ####
library(ggtree)
library(treeio)
library(ggplot2)

#### FUNCTIONS ####
build_tree <- function(tree_file, outgroup_name, res_labels) {
  # read the tree file and re-root it
  coli_tree <- read.newick(tree_file, node.label = "label")
  coli_tree_rooted <- root(coli_tree, outgroup = outgroup_name)

  # read in the HR classification labels
  hr_labels <- read.csv(res_labels, row.names = 1)
  # tree plot object
  tree_plot <- ggtree(coli_tree_rooted,
                      branch.length = "",
                      layout = "rectangular", ladderize = TRUE, size = .6) +
    geom_strip("DA63750", "DA63034",
               barsize = 4, color = "darkred",
               label = "B1",
               offset.text = 0.005, offset = -0.12, fontsize = 10) +
    geom_strip("DA63090", "DA62946",
               barsize = 4, color = "blue",
               label = "A", offset.text = 0.005,
               offset = -0.12, fontsize = 10) +
    geom_strip("DA63836", "DA63628",
               barsize = 4, color = "magenta2",
               label = "D1",
               offset.text = 0.005, offset = -0.12, fontsize = 10) +
    geom_strip("DA63974", "DA63328",
               barsize = 4, color = "darkorange",
               label = "D2/F", offset.text = 0.005, offset = -0.12,
               fontsize = 10, extend = TRUE) +
    geom_strip("DA63268", "DA63010",
               barsize = 4, color = "darkgreen",
               label = "B2", offset.text = 0.005,
               offset = -0.12, fontsize = 10) +
    geom_strip("O157_H7_Sakai", "O55_H7_CB9615",
               barsize = 4, color = "purple",
               label = "E", offset.text = 0.005,
               offset = -0.12, fontsize = 7)

  # add annotation
  tree_annotated <- gheatmap(tree_plot, hr_labels,
                             offset = .05, width = 0.05,
                             font.size = 3, hjust = 0, colnames = FALSE) +
    scale_fill_manual(breaks = c("nonHR", "HR", "R"),
                      values = c("steelblue", "firebrick", "gold"),
                      name = "phenotype")
  return(tree_annotated)
}

#### RUN ####
# build tree filename becuse the inptu is directory
tree_file_name <- paste(c(snakemake@input[[1]],
                          snakemake@config[["filename"]]),
                        sep = "/")
tree_plot_obj <- build_tree(tree_file = tree_file_name,
                            res_labels = snakemake@input[[2]],
                            outgroup_name = snakemake@config[["outgroup"]])
# save it to file
ggsave(file = snakemake@output[[1]],
       plot = tree_plot_obj,
       width = snakemake@config[["width"]],
       height = snakemake@config[["height"]],
       units = snakemake@config[["units"]])
print("The tree was saved. No errors.")

#### CLOSE LOG ####
sink()
