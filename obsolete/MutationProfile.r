q# NOTES:
# * adjusted for analysing reads mapped to the spike gene
# * lofreq options are outdated

# clear workspace
rm(list = ls())

# set working directory
#setwd("dir")

# library install and import
list.of.packages <- c("optparse", "pheatmap", "RColorBrewer", "dplyr", "fs", "tools", "stringr", "data.table", "vcfR")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)


# INPUT:
# 1. path to folder
# 2. pattern to cath vcf files in folder

# options 
args <- commandArgs(TRUE) 
cat("input files:\t", args[1], "/", args[2], "\n", sep="")

outfile <- args[3]
POS_lower <- 21563                  ## adjust the variant frequency
POS_upper <- 25384                  
percent <- F                        ## should variants be shown as percent?
rounded <- F                        ## should rounded values be displayed in the map
number <- T                         ## do the file names contain a number and you want to sort
clustering <- F                     ## should the samples be clustered?



# data preparation
## import data
files<-list.files(path = args[2], pattern = args[3], recursive = F, full.names = TRUE)
print(files)
## create data.frame for first dataset
if (type == 'ivar'){
  variants<-read.table(files[1], header = F, col.names=c('REGION','POS','REF','ALT', 'REF_DP','REF_RV','REF_QUAL','ALT_DP','ALT_RV','ALT_QUAL','ALT_FREQ','TOTAL_DP','PVAL','PASS','GFF_FEATURE','REF_CODON','REF_AA','ALT_CODON','ALT_AA'), sep = "\t")
  AF <- variants$INFO %>% strsplit(';') %>% sapply("[", 2)
  AF <- AF %>% strsplit('=') %>% sapply("[", 2)
  variants$AF <- AF
}
# medaka
if (type == 'medaka'){
  variants<-read.table(files[1], header = F, col.names=c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE', '2:SAMPLE'), sep = "\t", colClasses=c("character", "integer", rep("character",3), "numeric", rep("character", 5)))
  AF <- rep(1, nrow(variants))
  variants$AF <- AF
  variants$POS = variants$POS + POS_lower
}

variants <- variants[variants$POS > POS_lower, ]
variants <- variants[variants$POS < POS_upper, ]
print(variants)
final<- paste0(variants$REF, variants$POS, variants$ALT)
final<- data.frame(cbind(final, variants$AF))
Sample_number<-path_file(files[1])
Sample_number <- sapply(strsplit(Sample_number, ".", fixed=T), '[', 1)
colnames(final) <- c("Mutation", Sample_number)
print(final)

## loop over all other datasets and join by Mutation
for(i in 2:length(files)) {
  # lofreq variant calls
  if (type == 'lofreq'){
    variants<-read.table(files[i], header = F, col.names=c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'), sep = "\t")
    AF <- variants$INFO %>% strsplit(';') %>% sapply("[", 2)
    AF <- AF %>% strsplit('=') %>% sapply("[", 2)
    variants$AF <- AF
  }
  # medaka variant calls
  if (type == 'medaka'){
    variants<-read.table(files[i], header = F, col.names=c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','SAMPLE', '2:SAMPLE'), sep = "\t", colClasses=c("character", "integer", rep("character",3), "numeric", rep("character",5)))
    AF <- rep(1, nrow(variants))
    variants$AF <- AF
    variants$POS = variants$POS + POS_lower
  }
  variants <- variants[variants$POS > POS_lower, ]
  variants <- variants[variants$POS < POS_upper, ]
  varpos<- paste0(variants$REF, variants$POS, variants$ALT)
  varpos<- cbind(varpos, variants$AF)
  Sample_number<-path_file(files[i])
  Sample_number <- sapply(strsplit(Sample_number, ".", fixed=T), '[', 1)
  colnames(varpos) <- c("Mutation", Sample_number)
  final <- full_join(final, varpos, by="Mutation", copy=T)
}

## sort
final <- final[with(final, order(as.numeric(gsub('\\D', '', Mutation)))), ]
## set row names
row.names(final)<-final$Mutation
final<-final[,-1]
## order samples after number in name
if(number){final<-final[,str_order(colnames(final), numeric = T)]}
final <- t(final)
final[is.na(final)] <- 0
class(final) <- "numeric"
# convert values to percent
if(percent){final<-final*100}
# show rounded values in pheatmap
if(rounded){final_rounded <- round(final, 0)} else {final_rounded <- F}

print(final)

# visualize heatmap
## colormangment heatmap
my_colors <- colorRampPalette(c("grey93","brown","black"))

## The pdf scales with the number of variants and samples
png(paste0(outfile, ".png"))
if (type == 'lofreq'){ title <- "LoFreq variants" }
if (type == 'medaka'){ title <- "Medaka variants plotted as present (1) or absent (0)" }
pheatmap(final, 
         main = title,
         color = my_colors(100),
         cluster_rows = clustering,
         cluster_cols = clustering,
         display_numbers = final_rounded,
         legend = T,
         fontsize = 10
         )
dev.off()
