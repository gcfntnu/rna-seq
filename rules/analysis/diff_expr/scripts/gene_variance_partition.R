################################################################################
### R script to compare several complex desings using the variancePartition package
###
### Requires installed packages: optparse, variancePartition, doParallel
################################################################################

rm(list=ls())
library(optparse)
library(variancePartition)
library(doParallel)

# options list with associated default value.
option_list <- list(
    
    make_option(c("-P", "--projectName"),
                default=basename(getwd()),
                dest="projectName",
                help="name of the project used for the report [default: name of the current directory]."),
    
    make_option(c("-A", "--author"),
                default=Sys.info()[7],
                dest="author",
                help="name of the report author [default: %default]."),
    
    make_option(c("-t", "--targetFile"),
                default="target.txt",
                dest="targetFile",
                help="path to the design/target file [default: %default]."),

    make_option(c("-r", "--countsFile"),
                default="data/processed/rnaseq/quant/salmon/gene.quant",
                dest="countsFile",
                help="path to the count matrix file [default: %default]."),
    
    make_option(c("-m", "--metaFile"),
                default="data/processed/rnaseq/quant/salmon/gene_info.tsv",
                dest="metaFile",
                help="path to the features info file [default: %default]."),
    
    make_option(c("-R", "--templateFile"),
                default="src/rna-seq/rules/analysis/diff_expr/scripts/GCF_DESeq2.rmd",
                dest="templateFile",
                help="path to the directory R markdown template [default: %default]."),
    
    make_option(c("-F", "--featuresToRemove"),
                default="alignment_not_unique,ambiguous,no_feature,not_aligned,too_low_aQual",
                dest="FTR",
                help="names of the features to be removed, more than once can be specified [default: %default]"),
    
    make_option(c("-v", "--varInt"),
                default="group",
                dest="varInt", 
                help="factor of interest [default: %default]"),
    
    make_option(c("-c", "--condRef"),
                default="WT",
                dest="condRef",
                help="reference biological condition [default: %default]"),
    
    make_option(c("-b", "--batch"),
                default=NULL,
                dest="batch",
                help="blocking factor [default: %default] or \"batch\" for example"),

    make_option(c("-f", "--fitType"),
                default="parametric",
                dest="fitType", 
                help="mean-variance relationship: [default: %default],local or mean"),

    make_option(c("-o", "--cooksCutoff"),
                default=TRUE,
                dest="cooksCutoff", 
                help="perform the outliers detection (default is TRUE)"),
    
    make_option(c("-i", "--independentFiltering"),
                default=TRUE,
                dest="independentFiltering",
                help="perform independent filtering (default is TRUE)"),

    make_option(c("-a", "--alpha"),
                default=0.05,
            dest="alpha", 
            help="threshold of statistical significance [default: %default]"),
    
    make_option(c("-p", "--pAdjustMethod"),
                default="BH",
                dest="pAdjustMethod", 
                help="p-value adjustment method: \"BH\" or \"BY\" [default: %default]"),
    
    make_option(c("-T", "--typeTrans"),
                default="VST",
                dest="typeTrans", 
                help="transformation for PCA/clustering: \"VST\" ou \"rlog\" [default: %default]"),
    
    make_option(c("-l", "--locfunc"),
                default="median",
                dest="locfunc", 
                help="median or shorth to estimate the size factors [default: %default]"),
    
    make_option(c("-C", "--colors"),
                default="dodgerblue,firebrick1,MediumVioletRed,SpringGreen,chartreuse,cyan,darkorchid,darkorange",
                dest="cols",
                help="colors of each biological condition on the plots\n\t\t\"col1,col2,col3,col4\"\n\t\t[default: %default]"),
    
    make_option(c("-O", "--output"),
                default="data/processed/rnaseq/sartools",
                dest="output",
                help="output directory [default: %default]"),
    
    make_option(c("--forceCairoGraph"),
                action="store_true",
                default=FALSE,
                dest="forceCairoGraph",
                help="activate cairo type")
    
)

# now parse the command line to check which option is given and get associated values
parser <- OptionParser(usage="usage: %prog [options]",
					   option_list=option_list, 
					   description="Compare two or more biological conditions in a RNA-Seq framework with DESeq2.",
					   epilogue="For comments, bug reports etc... please contact Hugo Varet <hugo.varet@pasteur.fr>")
opt <- parse_args(parser, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options
