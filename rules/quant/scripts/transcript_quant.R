#/usr/bin/env Rscript

if (!require(argparse)) {
    install.packages("argparse", repos="http://cran.rstudio.com")
    library("argparse")
}



parser <- ArgumentParser(description="Transcript quantification from Kallisto/Salmon/Sailfish/RSEM")

parser$add_argument("input", nargs="+", help="Kallisto/Salmon/Sailfish/RSEM files")

parser$add_argument("--txinfo", required=TRUE,
                    help="Transcript info (required). Tab delimited file, needs columns with `gene_id`, `transcript_id`")

parser$add_argument("--ginfo", nargs="?",
                    help="Gene info (optional). Tab delimited file, needs column with `gene_id`")

parser$add_argument("-t", "--type", type="character", default="kallisto",
                    help="Data origins (kallisto, salmon, sailfish, rsem)")

parser$add_argument("-o", "--output", default="data/processed/salmon",
                    help="Output directory")

parser$add_argument("--output-scaled-tpm", action="store_true",
                    default=FALSE,
                    help="Output count (.quant) estimates derived from TPM scaled by library size")

parser$add_argument("--output-length-scaled-tpm", action="store_true",
                    default=FALSE,
                    help="Output count (.quant) estimates derived from TPM scaled by library size and gene lengths")

parser$add_argument("--output-genelength", action="store_true",
                    default=FALSE,
                    help="Additional output of (effective) genelength table")

parser$add_argument("--output-transcripts", action="store_true",
                    default=FALSE,
                    help="Additional output Counts/TPM per transcript.")

parser$add_argument("-v", "--verbose", action="store_true",
                    default=FALSE, help="Print extra output")

args <- parser$parse_args(args=commandArgs(TRUE))


if (!require(tximport)) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("tximport")
    library("tximport")
}
if (!require(data.table)) {
    install.packages("data.table", repos="http://cran.rstudio.com")
    library("data.table")
}


## add sample names to input files
nn <- as.character(sapply(args$input, function(x) {a <- strsplit(x, "/")[[1]]; a[[length(a)-1]]}))
names(args$input) <- nn

if (args$verbose == TRUE){
    print(args)
    options(echo=TRUE)
}

main <- function(args){
    tx.info <- fread(args$txinfo)
    keep <- colSums(is.na(tx.info)) != dim(tx.info)[1]
    tx.info <- tx.info[,..keep]
    tx2gene <- tx.info[,c("transcript_id", "gene_id")]
    
    reader <- function(x, ...) fread(x)
    importer <- function(x, ...) reader(x)
    countsFromAbundance <- "no"
    if (args$output_scaled_tpm == TRUE){
        countsFromAbundance <- "scaledTPM"
        if (args$output_length_scaled_tpm == TRUE){
            stop("Set *one* of output-scaled-tpm, output-length-scaled-tpm")
        }
    }
    if (args$output_length_scaled_tpm == TRUE){
        countsFromAbundance <- "lengthScaledTPM"
    }

    
    txi.tx <- tximport(args$input, type=args$type, tx2gene=tx2gene,
                       importer=importer, txOut=TRUE)
    txi <- summarizeToGene(txi.tx, tx2gene, countsFromAbundance=countsFromAbundance)
    
    gene.quant.fn <- file.path(args$output, "genes.quant")
    write.table(txi$counts, file=gene.quant.fn, sep="\t", quote=FALSE)
    tpm.fn <- file.path(args$output, "genes.tpm")
    write.table(txi$abundance, file=tpm.fn, sep="\t", quote=FALSE)

    if (args$output_transcripts == TRUE){
        tx.quant.fn <- file.path(args$output, "transcripts.quant")
        write.table(txi.tx$counts, file=tx.quant.fn, sep="\t", quote=FALSE)
        cat("Wrote file: ", tx.quant.fn)
        tpm.fn <- file.path(args$output, "transcripts.tpm")
        write.table(txi.tx$abundance, file=tpm.fn, sep="\t", quote=FALSE)
        tx.info <- as.data.frame(tx.info)
        rownames(tx.info) <- tx.info[,"transcript_id"]
        tx.info <- tx.info[rownames(txi.tx$counts),]
        tx.info.fn <- file.path(args$output, "transcript_info.tsv")
        write.table(tx.info, file=tx.info.fn, sep="\t", quote=FALSE, row.names=FALSE)
    }
    
    if (args$output_genelength == TRUE){
        ## transcript
        tpm.fn <- file.path(args$output, "transcripts.length")
        write.table(txi.tx$length, file=tpm.fn, sep="\t", quote=FALSE)
        
        ## gene
        tpm.fn <- file.path(args$output, "genes.length")
        write.table(txi$length, file=tpm.fn, sep="\t", quote=FALSE)
    }

    if (!is.null(args$ginfo)){
        gene.info <- read.delim(args$ginfo, sep=",", row.names=1)
        keep <- colSums(is.na(gene.info)) != dim(gene.info)[1]
        gene.info <- gene.info[,keep]
        gene.info <- gene.info[rownames(txi$counts),]
        
        gene.info.fn <- file.path(args$output, "gene_info.tsv")
        write.table(gene.info, file=gene.info.fn, sep="\t", quote=FALSE, row.names=FALSE)
    }

}

if (args$verbose == TRUE){
    debug(main)
}


main(args)
