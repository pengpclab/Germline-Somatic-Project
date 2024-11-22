
# Title: ActiveDriverWGS Analysis Script
# Author: Justin Ping-Hung Lai, Pei-Chen Peng
# Date: 11-21-2024
# Description: Analyze genomic data for enriched genomic regions using ActiveDriverWGS.


#
# Load required libraries
#
library(ActiveDriverWGS)
library(tidyverse)
library(GenomicRanges)


#
# Function to clean raw genomic data
#
clean <- function(rawdata){
  validChr = c(paste("chr",1:22,sep=""),"chrX","chrY")
  if (length(grep("chr",rawdata$chr))==0){
    rawdata$chr = paste("chr",rawdata$chr,sep="")
  }
  rawdata <- rawdata %>%
    filter(chr%in%validChr) %>%
    mutate_if(is.factor, as.character)
  return(rawdata)
}


#
# Variables
#
args = commandArgs(trailingOnly=TRUE)
cohort = args[1]
annoFile = args[2]
annoName = args[3]
genomeBuild="hg19"


#
# Load cohort-specific mutation data
#
if (cohort=="POCROC"){ ## patient with multiple samples will only be counted once
  data.mutations = read.table("/path/to/POCROC_concat.CMG_pos_head.vcf",header=T)
  data.mutations$patient = gsub("\\-.*","",data.mutations$patient)
}else if (cohort=="PCAWG.ovary"){
  data.mutations = read.table("/path/to/PCAWG_OvCa_merged.vcf",header=T)
}

#
# Clean mutation data 
#
data.mutations = data.mutations %>% distinct()
data.mutations = clean(data.mutations)
data.mutations = data.mutations[!(data.mutations$ref == ""), ]
data.mutations$pos1 = as.numeric(data.mutations$pos1)
data.mutations$pos2 = as.numeric(data.mutations$pos2)


#
# Load and clean annotation data
#
data.anno = read.table(annoFile)
data.anno$id = paste(data.anno$V1,"_",data.anno$V2,"_",data.anno$V3,":",annoName,sep="")
colnames(data.anno) = c("chr","start","end","id")
data.anno = clean(data.anno)
data.anno$start = as.numeric(data.anno$start)
data.anno$end = as.numeric(data.anno$end)


#
# Run ActiveDriverWGS analysis
#
finalresults = c()
results = c()
validChr = c(paste("chr",1:22,sep=""),"chrX","chrY")	
for (selectedchr in validChr){
	print(paste("Analyzing:", selectedchr))
	data.anno.sub = data.anno %>% filter(chr==selectedchr)
	tryCatch({
	  results = c()
	  results = ActiveDriverWGS(mutations = data.mutations, elements= data.anno.sub, ref_genome= genomeBuild)
	}, error=function(e){cat(selectedchr)})
	if (is.null(results)){next}	
	finalresults = rbind(finalresults, results %>% filter(element_enriched==TRUE& pp_element<=0.05))
}


#
# Save results
#
OutFile=paste("/path/to/result/folder/",cohort,"_",annoName,"_activeDriver.txt",sep="")
write.table(finalresults,OutFile,sep="\t",row.names=F, col.names=T,quote=F)
