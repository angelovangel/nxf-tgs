#!/usr/bin/env Rscript

# make valid bed files from feature_table.txt
# because the bed files produced by wf-clone-validation are not correct
# arg[1] is feature_table.txt, output is written to file
# example feature_table.txt

#Sample_name,Feature,Database,Identity,Match Length,Description,Start Location,End Location,Length,Strand,Plasmid length
#pGEM2,f1 ori,snapgene,100.0%,100.0%,f1 bacteriophage origin of replication; arrow indicates direction of (+) strand synthesis ,1619,2048,429,1,3197
#pGEM2,AmpR promoter,snapgene,100.0%,100.0%,bla,2406,2511,105,1,3197
#pGEM2,penA,swissprot,68.8%,10.2%,PENA_BURM1 - Experimental evidence at transcript level: Swiss-Prot protein existence level 2. Upon expression in E.coli enables the latter to utilize penicillin as a carbon source. From Burkholderia multivorans (strain ATCC 17616 / 249).,1207,1303,96,1,3197
#pGEM3,pVS1 RepA,snapgene,100.0%,100.0%,"replication protein from the Pseudomonas plasmid pVS1 (Heeb et al., 2000)",3109,4183,1074,1,18989

#library(readr)
#library(dplyr)
#library(stringr)

arg <- commandArgs(trailingOnly = T)

df <- utils::read.delim(arg[1], sep = ',', check.names = F)
# identity % as numbers
df$Identity <- as.numeric(gsub(pattern = "%", "", df$Identity))
# change strand to +/-
df$Strand <- replace(df$Strand, df$Strand == 1, "+")
df$Strand <- replace(df$Strand, df$Strand == -1, "-")
# select only required columns
df <- df[c('Sample_name', 'Start Location', 'End Location', 'Feature', 'Identity', 'Strand', 'Plasmid length')]
  #dplyr::select(Sample_name, `Start Location`, `End Location`, Feature, Identity, Strand, `Plasmid length`)
if (nrow(df) < 1) {
  #quit(status = 0)
  stop('No features in feature_table.txt')
}

# actual fix of features spanning origin
df1 <- df[df$`Start Location` < df$`End Location`, ]
# these have to be regenerated in 2 rows - 1..end and start..len
df2 <- df[df$`Start Location` > df$`End Location`, ]
df2a <- df2
df2b <- df2
if(nrow(df2) > 0) {
  df2a$`Start Location` <- 1
  df2b$`End Location` <- df2$`Plasmid length`  
}


finaldf <- rbind(df1, df2a, df2b)
finaldf <- subset(finaldf, select = -c(`Plasmid length`))

# check all is ok
if(nrow(finaldf[finaldf$`Start Location` > finaldf$`End Location`, ]) != 0) {
  stop('Check the feature_table.txt file!')
}

# save bed files
s <- split.data.frame(finaldf, finaldf$Sample_name)
snames <- names(s)
mapply(function(x,y) 
  utils::write.table(x, file = paste0(y, '.annotations2.bed'), sep = '\t', col.names = F, row.names = F, quote = F), s, snames
  )



  
  
