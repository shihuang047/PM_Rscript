#################################################################
# Function: Multivariate statistical analysis based on distance matrix for a certain sample category
# Call: Rscript DistBoxplot.R -m map_file -d dist_file -i HostID -o output
# R packages used: reshape,ggplot2,pheatmap,combinat,plyr,vegan,optparse
# Last update: 2018-03-22, Shi Huang
#################################################################

## install necessary libraries
p <- c("reshape2","ggplot2","pheatmap","combinat","plyr","vegan","optparse")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="http://cran.us.r-project.org/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

## clean R environment
rm(list = ls())
setwd('./')

## parsing arguments
args <- commandArgs(trailingOnly=TRUE)
args <- commandArgs(trailingOnly=TRUE)
sourcedir <- Sys.getenv("PM_Rscript")
source(sprintf('%s/util.R',sourcedir))

# make option list and parse command line
option_list <- list(
    make_option(c("-d", "--dist_file"), type="character", help="Input distance matrix table [required]."),
    make_option(c("-n", "--dist_name"), type="character", default='Default', help="The distance measure used such as Jensen-Shannon, Bray-Curtis, Euclidean et al.  [default %default]"),
    make_option(c("-m", "--map_file"), type="character", help="Input metadata mapping file [required]."),
    make_option(c("-p", "--prefix"), type="character",default='Out', help="Output file prefix [Optional, default %default]"),
    make_option(c("-c", "--category"), type="character", help="Metadata categories you want to investigate in distance matrix, e.g. \"Timepoint,Status,HostID\" [required]"),
    make_option(c("-i", "--individual_id"), type="character", default=NULL, help="The subject categories to specify all within-subject distance required, e.g. \"HostID\"  [default %default]"),
    make_option(c("-o", "--outdir"), type="character", default='DistBoxplot', help="Output directory [default %default]")
    )
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# Error checking
if(is.null(opts$map_file)) stop('Please supply a mapping file.')
if(is.null(opts$dist_file)) stop('Please supply an distance matrix table.')
if(is.null(opts$category)) stop('Please supply a category in the mapping file.')

# create output directory if needed
#if(opts$outdir != ".") dir.create(opts$outdir,showWarnings=FALSE, recursive=TRUE)
dir.create(outpath<-paste(opts$outdir,"/",sep=""),showWarnings=FALSE, recursive=TRUE)

filename<-opts$dist_file
metadata.filename<-opts$map_file 
individual_id<-opts$individual_id
dm_name<-opts$dist_name          
group<-opts$category
prefix_name<-opts$prefix

con <- file(paste(opts$outdir,"/",prefix_name,".DistBoxplot.",group,".Ind_",individual_id,".log",sep=""))
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")
#--------------------------------
dm<-read.table(filename,header=T,row.names=1,sep="\t")
dm<-dm[order(rownames(dm)),order(colnames(dm))]
allmetadata<-read.table(metadata.filename,header=T,sep="\t",row.names=1)
if(length(allmetadata)==1){metadata<-data.frame(allmetadata[order(rownames(allmetadata)),])
    all_group<-colnames(metadata)<-colnames(allmetadata)
}else{
    metadata<-allmetadata[order(rownames(allmetadata)),which(sapply(allmetadata,var)!=0)]
    all_group<-colnames(metadata)
    all_group_f<-colnames(metadata)[sapply(metadata,class)=="factor"]
    all_group_n<-colnames(metadata)[sapply(metadata,class)!="factor"]
}
#--------------------------------Data Check
if(any((colnames(dm)==rownames(dm))==FALSE))
  {cat("The column names do not exactly match the row names! Please revise!")}
if(any((rownames(metadata)==rownames(dm))==FALSE)) 
  {cat("The row names in Map file do not exactly match the row names in the distance matrix! Please revise!\n")}
if(!is.element(group, colnames(metadata))) stop('Please make sure that category existed in the mapping file.')

if(is.factor(metadata[,group])){
if(is.null(individual_id)){
d<-DistBoxplot(dm,dm_name=dm_name,group=metadata[,group],group_name=group,IndividualID=NULL,outpath=outpath)}else{
d<-DistBoxplot(dm,dm_name=dm_name,group=metadata[,group],group_name=group,IndividualID=metadata[,individual_id],outpath=outpath)
}

ano.P<-anosim(dm,metadata[,group])$signif
ano.R<-anosim(dm,metadata[,group])$statistic
    print(paste(group," P value (ANOSIM)=", ano.P,sep=""))
    print(paste(group," R value (ANOSIM)=", ano.R,sep=""))
    cat("\n")

ado<-adonis(dm~metadata[,group])
ado.P<-ado$aov.tab$P[1]
ado.F<-ado$aov.tab$F.Model[1]
    ado$aov.tab
    print(paste(group," P value (ADONIS/PERMANOVA)=", ado.P,sep=""))
    print(paste(group," F value (ADONIS/PERMANOVA)=", ado.F,sep=""))
    cat("\n")
}else{
    cat("The category you specified should be a categorical variable in metadata! If NOT, please revise!\n")}

#--------------------------------
# Correlation with continuous variable(s) in metadata
#--------------------------------    
if(length(all_group_n)>=1){

for(n in all_group_n) {
    dir.create(outpath1<-paste(outpath,dm_name,".Btw_",group,"_VS._",n,"/",sep=""))
    dm_n<-data.matrix(dist(metadata[,n])); colnames(dm_n)<-rownames(dm_n)<-rownames(metadata)
    if(is.null(individual_id)){
        d_n<-DistBoxplot(dm_n,dm_name="Euclidean",group=metadata[,group],group_name=group,IndividualID=NULL,outpath=outpath1)}else{
        d_n<-DistBoxplot(dm_n,dm_name="Euclidean",group=metadata[,group],group_name=group,IndividualID=metadata[,individual_id],outpath=outpath1)
        }
    if(var(d_n$Dist)!=0){
    d_comp<-data.frame(d,Dist_n=d_n$Dist)
    corr<-with(d_comp,cor.test(Dist,Dist_n,method="spearman"))
    p<-ggplot(data=d_comp, aes(x=Dist_n,y=Dist))+geom_point(alpha=0.2)+ 
       annotate("text", x=max(d_comp$Dist_n)*0.9, y=max(d_comp$Dist)*0.9, label= paste("Rho=",round(corr$estimate,2),"\n","P=",round(corr$p.value,4),"\n",sep="")) +
       ylab(paste(dm_name," Distance",sep=""))+ 
       xlab(bquote(paste(~Delta, .(n),sep=" ")))+
       theme_bw()
       #if(corr$estimate>0.5 && corr$p.value<0.01){p<-p+geom_smooth(method = "loess", se=TRUE, span=1)}
    
    ggsave(filename=paste(outpath1,"/",prefix_name,".",dm_name,".Btw_",group,"_VS._",n,".Scatterplot.pdf",sep=""), plot=p, limitsize=TRUE, width=4, height=4)
    
    sink(paste(outpath1, dm_name, ".Btw_",group, "_VS._", n ,".dm_values.xls",sep=""))
    write.table(d_comp,quote=FALSE,sep="\t",row.names=FALSE)
    sink()
    
    }else{
    cat("The variance of Delta",n, "= 0\n",sep=" ")
    }
}
}

sink()
sink(type="message")

