library(BSgenome.Hsapiens.UCSC.hg19.masked)
genome <- BSgenome.Hsapiens.UCSC.hg19
out_file <- file.path(snakemake@output[["hg19"]])
export(genome, out_file)
