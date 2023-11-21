# testing tools and approaches for variant calling

library(GenomicRanges)

# beta-lac TEM-1 copy number increase in DA62886 mutant
bam_path <-"~/Data/HeteroR/results/snippy/DA62886/snps.bam"
bed_path <- "~/Data/HeteroR/results/snippy/DA62886/snps.bed"

gr <- GRanges("BNIFGIJN_1", IRanges(681022, 681023))

bfls <- dir(".", "bam$")

t1 <- autoplot(TxDb.Mmusculus.UCSC.mm10.knownGene, which = gr)

t2 <- lapply(bfls, autoplot, which = gr)

names(t2) <- bfls

t2 <- tracks(t2)

t1 <- tracks(t1)

c(t1,t2)