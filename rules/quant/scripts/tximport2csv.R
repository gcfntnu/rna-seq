#/usr/bin/env Rscript


library(argparse)
library(tximport)
library(readr)


parser <- ArgumentParser(description="tximport obj csv export")

parser$add_argument("input", help="tximport obj RDS file")

parser$add_argument("--txinfo", required=TRUE,
                    help="Transcript info (required). Tab delimited file, needs columns with `gene_id`, `transcript_id`")

parser$add_argument("-t", "--type", type="character", default="gene",
                    help="csv count output option (gene, gene_tpm, gene_tpm_scaled, gene_tpm_length_scaled, tx, tx_tpm, tx_tpm_scaled, gene_length, variances)")

parser$add_argument("-o", "--output", required=TRUE, help="Output tsv file")

parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="Print extra output")

args <- parser$parse_args(args=commandArgs(TRUE))

tx.info <- readr::read <- tsv(args$txinfo)
tx2gene <- tx.info[,c("transcript_id", "gene_id")]

txi.tx <- readRDS(args$input)

if (args$type == "gene"){
    out <- summarizeToGene(txi.tx, tx2gene)$counts
} else if (args.type == "gene_tpm"){
    out <- summarizeToGene(txi.tx, tx2gene)$abundance
} else if (args.type == "gene_tpm_scaled"){
    out <- summarizeToGene(txi.tx, tx2gene, countsFromAbundance="scaledTPM")$counts
} else if (args.type == "gene_tpm_length_scaled"){
    out <- summarizeToGene(txi.tx, tx2gene, countsFromAbundance="lengthScaledTPM")$counts
} else if (args.type == "gene_length"){
    out <- summarizeToGene(txi.tx, tx2gene)$length
} else if (args.type == "tx"){
    out <- txi.tx$counts
} else if (args.type == "tx_tpm"){
    out <- txi.tx$abundance
} else if (args.type == "tx_tpm_scaled"){
    length4CFA <- tximport:::medianLengthOverIsoform(txi.tx$length, tx2gene)
    out <- tximport:::makeCountsFromAbundance(countsMat = txi.tx$counts, 
                                              abundanceMat = txi.tx$abundance,
                                              lengthMat=length4CFA,
                                              countsFromAbundance = "lengthScaledTPM") 
} else if (args.type == "tx_length"){
    out <- txi.tx$length
} else if (args.type == "variances"){
    if (! args.type in c("salmon", "kallisto")){
        stop("only salmon, kallisto estimate variances!")
        }
    out <- summarizeToGene(txi.tx, tx2gene, varReduce=TRUE)$variance
}



write_tsv(out, args$output, delim="\t")






