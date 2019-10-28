#################################################################
# Function: To make distance matrices for a data matrix
# Call: Rscript PM_makedist.R -i matrix_file -o output
# R packages used: reshape,vegan,plyr,optparse
# Last update: 2018-03-22, Shi Huang
#################################################################

## clean R environment
rm(list = ls())
setwd('./')
## install necessary libraries
p <- c("reshape","ggplot2","vegan","plyr","vegan","optparse")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="http://mirrors.opencas.cn/cran/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

## parsing arguments
args <- commandArgs(trailingOnly=TRUE)
# make option list and parse command line
option_list <- list(
    make_option(c("-i", "--matrix_file"), type="character", help="Input data matrix table [Required]"),
    make_option(c("-o", "--outdir"), type="character", default='dist_matrix', help="Output directory [default %default]"),
    make_option(c("-p", "--prefix"), type="character",default='Out', help="Output file prefix [Optional, default %default]")
    )
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# Error checking
if(is.null(opts$matrix_file)) stop('Please supply a data matrix table.')

dir.create(outpath<-paste(opts$outdir,"/",sep=""),showWarnings=FALSE, recursive=TRUE)

filename<-opts$matrix_file
#dm_name<-opts$dist_name
prefix<-opts$prefix
#--------------------------------
JSD<-function(object, sweep=FALSE, eps=10^-4, overlap=TRUE,...)
{
    if(!is.numeric(object))
    stop("object must be a numeric matrix\n")
    
    z <- matrix(NA, nrow=ncol(object), ncol=ncol(object))
    colnames(z) <- rownames(z) <- colnames(object)
    
    w <- object < eps
    if (any(w)) object[w] <- eps
    ## If you takes as input a matrix of density values
    ## with one row per observation and one column per
    ## distribution, add following statement below.
    if(sweep){ object <- sweep(object, 2, colSums(object) , "/") }
    
    for(k in seq_len(ncol(object)-1)){
        for(l in 2:ncol(object)){
            ok <- (object[, k] > eps) & (object[, l] > eps)
            if (!overlap | any(ok)) {
                m=0.5*(object[,k]+object[,l])
                z[k,l] <- sqrt(0.5*sum(object[,k] *(log(object[,k]) - log(m)))+0.5*sum(object[,l] *(log(object[,l]) - log(m))))
                z[l,k] <- sqrt(0.5*sum(object[,l] *(log(object[,l]) - log(m)))+0.5*sum(object[,k] *(log(object[,k]) - log(m))))
            }
        }
    }
    diag(z)<-0
    z
}

#--------------------------------
data<-read.table(filename,header=T,row.names=1)
tdata<-t(data)
#--------------------------------------------------
dm_name<-c("JSD","bray","jaccard","euclidean")
#--------------------------------------------------
for(k in 1:length(dm_name)){
    if(k==1){ dm<-as.dist(JSD(tdata))
    }else{
        dm<-vegdist(data, method = dm_name[k])
    }
    sink(paste(outpath,prefix,".",dm_name[k],".dm.txt",sep="")); cat("\t"); write.table(as.matrix(dm),sep="\t",quote=FALSE); sink()
}




