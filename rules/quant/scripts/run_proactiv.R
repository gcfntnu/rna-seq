#!/usr/bin/env Rscript

## List of STAR junction files as input
##files <- list.files(system.file('extdata/vignette/junctions',package = 'proActiv'), full.names = TRUE)
## Vector describing experimental condition
#condition <- rep(c('A549','HepG2'), each=3)
## Promoter annotation for human genome GENCODE v34

#promoterAnnotation <- promoterAnnotation.gencode.v34.subset

library(argparse)
library(stringr)
library(proActiv)

parser <- ArgumentParser(description="Promoter quantification from STAR junction file")
parser$add_argument("input", nargs="+", help="Star junction files")
parser$add_argument("--sample-info", default=NULL, help="Sample info (required). Tab delimited file, needs columns with `sample_id`")
parser$add_argument("--condition", default="Sample_Group", help="Condition column name in sample_info file")
parser$add_argument("--ref-level", help="Condition reference level.")
parser$add_argument("--gtf", required=TRUE, help="Feature gtf file")
parser$add_argument("--species", default="Homo_sapiens", help="Feature gtf file")
parser$add_argument("-o", "--output", default="proactiv.rds", help="Output rds file")
parser$add_argument("-v", "--verbose", action="store_true",
                    default=FALSE, help="Print extra output")
args <- parser$parse_args(args=commandArgs(TRUE))
args$species <- stringr::str_to_title(args$species)
sample.names <- as.character(sapply(args$input, function(x) {a <- strsplit(x, "/")[[1]]; a[[length(a)-1]]}))
names(args$input) <- sample.names

if (args$verbose == TRUE){
    print(args)
    options(echo=TRUE)
}

anno <- preparePromoterAnnotation(file=args$gtf, species = args$species)

condition <- NULL
if (!is.null(args$sample_info)){
    tab <- read.delim(args$sample_info, sep="\t", row.names=1)
    tab <- tab[sample.names,]
    condition <- factor(tab[,args$condition])
    if (!is.null(args$ref_level)){
        condition <- relevel(condition, ref=args$ref_level)
    }
}


result <- proActiv(files = args$input, promoterAnnotation=input$gtf, condition=condition)
result <- result[complete.cases(assays(result)$promoterCounts),]
result$$countData <- data.frame(assays(result)$promoterCounts, rowData(result))

saveRDS(result, args$output, compress=FALSE)
