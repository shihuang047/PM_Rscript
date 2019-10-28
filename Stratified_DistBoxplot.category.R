#!/usr/bin/env Rscript
#################################################################
# Function: Multivariate statistical analysis based on distance matrix
# Call: Rscript PM_Bdiversity.R -m map_file -d dist_file -o output
# R packages used: reshape,ggplot2,pheatmap,pROC,combinat,plyr,vegan,optparse
# Last update: 2018-05-18, Shi Huang
#################################################################
options(warn=1)
#Rprof()
## install necessary libraries
p <- c("reshape2","ggplot2","pheatmap","pROC","combinat","plyr","vegan","optparse")
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
sourcedir <- Sys.getenv("PM_Rscript")
source(sprintf('%s/util.R',sourcedir))
# make option list and parse command line
option_list <- list(
make_option(c("-d", "--dist_file"), type="character", help="Input distance matrix table [Required]"),
make_option(c("-m", "--meta_data"), type="character", help="Input meta data file [Required]"),
make_option(c("-o", "--out_dir"), type="character", default='Out', help="Output directory [default %default]"),
make_option(c("-p", "--prefix"), type="character",default='Out', help="Output file prefix [Optional, default %default]"),
make_option(c("-s", "--strata"), type="character", help="The strata you specified should be found in the metadata. [Required]"),
make_option(c("-c", "--category_list"), type="character", help="The categories you specified should be found in the metadata. [Required]"),
make_option(c("-n", "--dist_name"), type="character", default='Default', help="The distance metrics name such as Meta-Storms, Jensen-Shannon, Euclidean et al. [Optional, default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# Error checking
if(is.null(opts$meta_data)) stop('Please input a meta data file')
if(is.null(opts$dist_file)) stop('Please input a distance matrix table')

# create output directory if needed
#if(opts$out_dir != ".")
#dir.create(outpath1<-paste(opts$out_dir,"/",sep=""),showWarnings=FALSE, recursive=TRUE)

prefix<-opts$prefix
filename<-opts$dist_file
metadata.filename<-opts$meta_data
dm_name<-opts$dist_name
strata<-opts$strata
category_list<-opts$category_list; category_list<-unlist(strsplit(category_list, ","))


con <- file(paste("./",prefix,".Stratified_DistBoxplot4Categories.log",sep=""))
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")

dir.create(outpath1<-paste(opts$out_dir,"_DistStratifiedBy_",strata,"/",sep=""),showWarnings=FALSE, recursive=TRUE)
#--------------------------------
dm<-read.table(filename, header=T, row.names=1); dm<-dm[order(rownames(dm)),order(colnames(dm))]

#--------------------------------metadata checking
allmetadata<-read.table(metadata.filename,header=T,sep="\t",row.names=1); 
metadata<-data.frame(allmetadata[order(rownames(allmetadata)), ]) 
colnames(metadata)<-colnames(allmetadata)
all_group<-colnames(metadata)
all_group_f<-colnames(metadata)[sapply(metadata,class)=="factor"]; try(metadata_f<-metadata[, all_group_f])
all_group_n<-colnames(metadata)[sapply(metadata,class)!="factor"]; try(metadata_n<-metadata[, all_group_n])

cat("All the sample metadata: ",all_group, "\n\n",sep=" ")
#--------------------------------Data Check
if(any((colnames(dm)==rownames(dm))==FALSE)) 
  {cat("The column names do not exactly match the row names! Please revise!")}
if(any((rownames(metadata)==rownames(dm))==FALSE)) 
  {cat("The row names in Map file do not exactly match the row names in the distance matrix! Please revise!\n")}
#--------------------------------To form a split_list
s<-metadata[,strata]
if(length(strata)==1){
   names(s)<-rownames(metadata)
   }else{
   split<-factor(apply(s, 1, paste,collapse="."))
}
split_list<-split(s, s)
#--------------------------------To delete sample category, of which the number of levels <2 in the splited metadata
sub_metadata_list<-list()
for(i in 1:length(split_list)){
sub_metadata_list[[i]]<-metadata[names(split_list[[i]]), category_list]
}
names(sub_metadata_list)<-names(split_list)
nlevels_sub_metadata<-sapply(sub_metadata_list,function(x) apply(x,2,function(x) nlevels(as.factor(x))))
del_category_list<-rownames(nlevels_sub_metadata)[which(nlevels_sub_metadata==1)%/%ncol(nlevels_sub_metadata)]
if(length(del_category_list)>0) {cat("The sample categories can NOT be analyzed in one of the splited data:", del_category_list, "\n")
category_list<-category_list[which(!(category_list %in% del_category_list))]}
cat("The sample categories finally analyzed in the splited data:", category_list, "\n")

#--------------------------------

all_dm_v_list<-list(); length(all_dm_v_list)<-length(category_list)
for(g in 1:length(category_list)){
   group<-category_list[g]
#--------------------------------For each group, we extract the dm values.
   dm_v_list<-list(); length(dm_v_list)<-length(split_list)
   stat_summ<-matrix(NA, nrow=6, ncol=length(split_list)); colnames(stat_summ)<-names(split_list); rownames(stat_summ)<-c("df", "Adonis_Pseudo-F","Adonis_R2","Adonis_P","Anosim_R","Anosim_P")
     #--------------------------------
     for(i in 1:length(split_list)){
     sub_dm<-dm[names(split_list[[i]]),names(split_list[[i]])]
     sub_dm_name<-names(split_list)[i]
     sub_metadata<-metadata[names(split_list[[i]]),]
     #--------------------------------
     dir.create(outpath2<-paste(outpath1,dm_name, ".", group, ".", sub_dm_name, "/",sep=""))
     #--------------------------------
     d<-DistBoxplot(sub_dm, sub_metadata[, group], dm_name=dm_name, group_name=group, outpath=outpath2)
     #--------------------------------AllBetween and AllWithin BoxPlot
     dm_value<-d$dm_value
     dm_value_AB<-subset(dm_value, DistType=="AllBetween")
     write.table(dm_value_AB, paste(outpath2, dm_name, '.', group, '.AllBetween.xls',sep=''), quote=F, sep="\t", row.names=F)
     p<-qplot(x=GroupPair, y=Dist,  data=dm_value_AB, geom='boxplot',main='',xlab="Group pair", ylab=paste(dm_name,'_Distance',sep='')) + coord_flip() + theme_bw()
     suppressMessages(ggsave(filename=paste(outpath2, dm_name, '.', group, '.AllBetween.boxplot.ggplot.pdf',sep=''),plot=p, height=ifelse(nlevels(metadata[,group])>2,nlevels(metadata[,group]),2), width=5))
     #----------------
     dm_value_AB_med<-ddply(dm_value_AB, .(GroupPair), summarize, MedianDist=median(Dist))
     GroupPairs<-do.call(rbind, strsplit(as.character(dm_value_AB_med$GroupPair), "_VS._"));colnames(GroupPairs)<-c("Group1","Group2")
     dm_value_AB_med<-data.frame(GroupPairs, dm_value_AB_med)
     p <- ggplot(dm_value_AB_med, aes(Group1, Group2)) + geom_tile(aes(fill = MedianDist), colour = "white") + 
     ggtitle(paste("Median", dm_name, "Distance", sep=" ")) + theme(plot.title = element_text(hjust = 0.5))+ 
     scale_fill_gradient(low = "white", high = "steelblue")
     suppressMessages(ggsave(filename=paste(outpath2, dm_name, '.', group, '.AllBetween.Median_heatmap.ggplot.pdf',sep=''),plot=p))
     #--------------------------------
     dm_value_AW<-subset(dm_value, DistType=="AllWithin")
     p<-qplot(x=GroupPair, y=Dist,  data=dm_value_AW, geom='boxplot',main='',xlab="Group pair", ylab=paste(dm_name,'_Distance',sep='')) + coord_flip() + theme_bw()
     suppressMessages(ggsave(filename=paste(outpath2, dm_name, '.', group, '.AllWithin.boxplot.ggplot.pdf',sep=''),plot=p, height=ifelse(nlevels(metadata[,group])>2,nlevels(metadata[,group]),2), width=5))
     #--------------------------------
     cat("\n\nstrata: ",sub_dm_name,"\n--------------------------------\n")
     #print(paste(sub_dm_name," All_between ",group," VS All_within ",group," P value (T-test)=", d$p_t,sep=""))
     #print(paste(sub_dm_name," All_between ",group," VS All_within ",group," P value (Wilcox-test)=", d$p_w,sep=""))
     stat_summ[1,i]<-length(split_list[[i]])-1
     #--------------------------------
     ano<-anosim(sub_dm,sub_metadata[,group])
     stat_summ[6,i]<-ano.P<-ano$signif
     stat_summ[5,i]<-ano.R<-round(ano$statistic,3)
     cat(paste(group, ".ANOSIM: \n", sep=""))
     cat("--------------------------------")
     print(ano)
     #--------------------------------
     ado<-adonis(sub_dm~sub_metadata[,group])
     stat_summ[4,i]<-ado.P<-ado$aov.tab$P[1]
     stat_summ[3,i]<-ado.R<-round(ado$aov.tab$R2[1],3)
     stat_summ[2,i]<-ado.F<-round(ado$aov.tab$F.Model[1],3)
     cat(paste(group, ".ADNOIS: \n", sep=""))
     cat("--------------------------------\n")
     print(ado$aov.tab)
     cat("--------------------------------\n\n")
     strata<-sub_dm_name
     dm_v_list[[i]]<-sub_dm_v<-data.frame(Group=rep(group, nrow(d$dm_value)), strata=rep(strata, nrow(d$dm_value)),d$dm_value)
     sink(paste(outpath2,prefix, ".", sub_dm_name,".txt",sep=""));cat("\t"); write.table(sub_dm,quote=FALSE,sep="\t");sink()
     sink(paste(outpath2,prefix,".", sub_dm_name,"_Map.txt",sep=""));cat("SampleID\t"); write.table(sub_metadata,quote=FALSE,sep="\t");sink()
     }
   all_dm_v_list[[g]]<-dm_v<-do.call(rbind, dm_v_list)
   #--------------------------------
   sink(paste(outpath1, prefix, ".", group, ".stat_summ.xls",sep=""));cat("\t");write.table(stat_summ,quote=FALSE,sep="\t");sink()
   #--------------------------------

}
    all_dm_v<-do.call(rbind, all_dm_v_list) 
         
        
p<-ggplot(all_dm_v, aes(x=Group, y=Dist)) + geom_boxplot(aes(fill=DistType), position="dodge") + ylab(paste(dm_name," Distance",sep=""))+ 
   coord_flip()+
   theme_bw()+
   facet_wrap(~strata)
all_groups<-paste(category_list, collapse="-")
ggsave(filename=paste(outpath1, prefix, ".", dm_name, ".", all_groups, ".ggplot.pdf",sep=""), plot=p, limitsize=TRUE, width=6, height=3)


sink()
sink(type="message")
