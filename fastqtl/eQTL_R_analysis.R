setwd('/Users/yocelyn/eQTLs/Actual_eQTLs')
library("VennDiagram")
#------------------------------------------------------------------------#-----------------------------------------------------------------------------------------------------------------------------------------------
#Europeans eQTLs chromosome 1  ##threshold for 0.00001 through fastqtl 
#read table
EUR_eqtl_1 <- read.table("EUR_eQTL_chr1.results.txt")
#add colnames
colnames(x = EUR_eqtl_1) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
EUR_eqtl_ordered_1 <- subset(x = EUR_eqtl_1, subset = pvalue <= 0.0000001)
#order p-value
EUR_eqtl_ordered_1 <- EUR_eqtl_ordered_1[order(EUR_eqtl_ordered_1$pvalue),]
#save SNP column 
EUR_eqtls_SNP_1 <- EUR_eqtl_ordered_1[,c(2,4), drop =FALSE]
#write csv for ordered EUR chr
#write.csv(EUR_eqtl_ordered_1, "EUR_eqtl_ordered_1.csv")
#extra commands 
chr1_MP_EUR = plot(-log10(x = EUR_eqtl_1$pvalue) , main='EUR_chr1')                 
#------------------------------------------------------------------------#-----------------------------------------------------------------------------------------------------------------------------------------------
#AA eQTLs chromosome 1 ##threshold for 0.00001 through fastqtl 
#read table
AA_eqtl_1 <- read.table("AA_eQTL_chr1.results.txt")

#rename columns 
colnames(x = AA_eqtl_1) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
AA_eqtl_ordered_1 <- subset(x = AA_eqtl_1, subset = pvalue <= 0.0000001)
#ordered by pvalue 
AA_eqtl_ordered_1 <- AA_eqtl_ordered_1[order(AA_eqtl_ordered_1$pvalue),]
##filter african american for population specific ones 
#subset eqtls not in europeans
AA_eqtl_ordered_specific_1 <- AA_eqtl_ordered_1[!(AA_eqtl_ordered_1$snpid %in% EUR_eqtl_ordered_1$snpid),]
#get names of eqtls ordered pvalue 
AA_eqtl_ordered_specific_names_1 <- AA_eqtl_ordered_specific_1[order(AA_eqtl_ordered_specific_1$pvalue),] 
#subset eqtls not in african americans and only in europeans
EUR_eqtl_ordered_specific_1 <- EUR_eqtl_ordered_1[!(EUR_eqtl_ordered_1$snpid %in% AA_eqtl_ordered_1$snpid),]

#Data for eQTLs number less than 0.00005
EUR_1 <- nrow(EUR_eqtl_ordered_1)
AA_1 <- nrow(AA_eqtl_ordered_1)
EUR_specific_1 <- nrow(EUR_eqtl_ordered_specific_1)
AA_specific_1 <- nrow(AA_eqtl_ordered_specific_names_1) 
shared_1 <- EUR_1 - EUR_specific_1

  


AA_eqtl_ordered_only_0.00001 <- subset(x = AA_eqtl_ordered_only, subset = pvalue <= 0.00001)
write.csv(AA_eqtl_ordered_specific_names_1, "AA_eqtl_ordered_specific_names_1.csv")     

#manhattan plot for chr 1
chr1_MP_AA = plot(-log10(x = AA_eqtl$pvalue), main='AA_chr1')

#EUR eQTLs chromosome 2 

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#read table
EUR_eqtl_2 <- read.table("EUR_eQTL_chr2.results.txt")
#add colnames
colnames(x = EUR_eqtl_2) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
EUR_eqtl_ordered_2 <- subset(x = EUR_eqtl_2, subset = pvalue <= 0.0000001)
#order p-value
EUR_eqtl_ordered_2 <- EUR_eqtl_ordered_2[order(EUR_eqtl_ordered_2$pvalue),]
#save SNP column 
EUR_eqtls_SNP_2 <- EUR_eqtl_ordered_2[,2, drop =FALSE]
#write csv for ordered EUR chr
#write.csv(EUR_eqtl_ordered_2, "EUR_eqtl_ordered_2.csv")
#extra commands 
chr2_MP_EUR = plot(-log10(x = EUR_eqtl_2$pvalue) , main='chr2')                 
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#AA eQTLs chromosome 2 ##threshold for 0.00001 through fastqtl 
#read table
AA_eqtl_2 <- read.table("AA_eQTL_chr2.results.txt")

#rename columns 
colnames(x = AA_eqtl_2) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
AA_eqtl_ordered_2 <- subset(x = AA_eqtl_2, subset = pvalue <= 0.0000001)
#ordered by pvalue 
AA_eqtl_ordered_2 <- AA_eqtl_ordered_2[order(AA_eqtl_ordered_2$pvalue),]
##filter african american for population specific ones 
#subset eqtls not in europeans
AA_eqtl_ordered_specific_2 <- AA_eqtl_ordered_2[!(AA_eqtl_ordered_2$snpid %in% EUR_eqtl_ordered_2$snpid),]
#get names of eqtls ordered pvalue 
AA_eqtl_ordered_specific_names_2 <- AA_eqtl_ordered_specific_2[order(AA_eqtl_ordered_specific_2$pvalue),] 
#subset eqtls not in african americans and only in europeans
EUR_eqtl_ordered_specific_2 <- EUR_eqtl_ordered_2[!(EUR_eqtl_ordered_2$snpid %in% AA_eqtl_ordered_2$snpid),]

#Data for eQTLs number less than 0.00005
EUR_2 <- nrow(EUR_eqtl_ordered_2)
AA_2 <- nrow(AA_eqtl_ordered_2)
EUR_specific_2 <- nrow(EUR_eqtl_ordered_specific_2)
AA_specific_2 <- nrow(AA_eqtl_ordered_specific_names_2) 
shared_2 <- EUR_2 - EUR_specific_2




#AA_eqtl_ordered_only_0.00001 <- subset(x = AA_eqtl_ordered_only, subset = pvalue <= 0.00001)
write.csv(AA_eqtl_ordered_specific_names_2, "AA_eqtl_ordered_specific_names_2.csv")
#-----------------------------------------------------------------------------------------------------------------------------------------------#------------------------------------------------------------------------

#EUR eQTLs chromosome 3 


#EUR eQTLs chromosome 3 

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#read table
EUR_eqtl_3 <- read.table("EUR_eQTL_chr3.results.txt")
#add colnames
colnames(x = EUR_eqtl_3) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
EUR_eqtl_ordered_3 <- subset(x = EUR_eqtl_3, subset = pvalue <= 0.0000001)
#order p-value
EUR_eqtl_ordered_3 <- EUR_eqtl_ordered_3[order(EUR_eqtl_ordered_3$pvalue),]
#save SNP column 
EUR_eqtls_SNP_3 <- EUR_eqtl_ordered_3[,3, drop =FALSE]
#write csv for ordered EUR chr
#write.csv(EUR_eqtl_ordered_3, "EUR_eqtl_ordered_3.csv")
#extra commands 
chr3_MP_EUR = plot(-log10(x = EUR_eqtl_3$pvalue) , main='chr3')                 
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#AA eQTLs chromosome 3 ##threshold for 0.00001 through fastqtl 
#read table
AA_eqtl_3 <- read.table("AA_eQTL_chr3.results.txt")

#rename columns 
colnames(x = AA_eqtl_3) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
AA_eqtl_ordered_3 <- subset(x = AA_eqtl_3, subset = pvalue <= 0.0000001)
#ordered by pvalue 
AA_eqtl_ordered_3 <- AA_eqtl_ordered_3[order(AA_eqtl_ordered_3$pvalue),]
##filter african american for population specific ones 
#subset eqtls not in europeans
AA_eqtl_ordered_specific_3 <- AA_eqtl_ordered_3[!(AA_eqtl_ordered_3$snpid %in% EUR_eqtl_ordered_3$snpid),]
#get names of eqtls ordered pvalue 
AA_eqtl_ordered_specific_names_3 <- AA_eqtl_ordered_specific_3[order(AA_eqtl_ordered_specific_3$pvalue),] 
#subset eqtls not in african americans and only in europeans
EUR_eqtl_ordered_specific_3 <- EUR_eqtl_ordered_3[!(EUR_eqtl_ordered_3$snpid %in% AA_eqtl_ordered_3$snpid),]

#Data for eQTLs number less than 0.00005
EUR_3 <- nrow(EUR_eqtl_ordered_3)
AA_3 <- nrow(AA_eqtl_ordered_3)
EUR_specific_3 <- nrow(EUR_eqtl_ordered_specific_3)
AA_specific_3 <- nrow(AA_eqtl_ordered_specific_names_3) 
shared_3 <- EUR_3 - EUR_specific_3




#AA_eqtl_ordered_only_0.00001 <- subset(x = AA_eqtl_ordered_only, subset = pvalue <= 0.00001)
write.csv(AA_eqtl_ordered_specific_names_3, "AA_eqtl_ordered_specific_names_3.csv")

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#EUR eQTLs chromosome 4 

#read table
EUR_eqtl_4 <- read.table("EUR_eQTL_chr4.results.txt")
#add colnames
colnames(x = EUR_eqtl_4) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
EUR_eqtl_ordered_4 <- subset(x = EUR_eqtl_4, subset = pvalue <= 0.0000001)
#order p-value
EUR_eqtl_ordered_4 <- EUR_eqtl_ordered_4[order(EUR_eqtl_ordered_4$pvalue),]
#save SNP column 
EUR_eqtls_SNP_4 <- EUR_eqtl_ordered_4[,4, drop =FALSE]
#write csv for ordered EUR chr
#write.csv(EUR_eqtl_ordered_4, "EUR_eqtl_ordered_4.csv")
#extra commands 
chr4_MP_EUR = plot(-log10(x = EUR_eqtl_4$pvalue) , main='chr4')                 
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#AA eQTLs chromosome 4 ##threshold for 0.00001 through fastqtl 
#read table
AA_eqtl_4 <- read.table("AA_eQTL_chr4.results.txt")

#rename columns 
colnames(x = AA_eqtl_4) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
AA_eqtl_ordered_4 <- subset(x = AA_eqtl_4, subset = pvalue <= 0.0000001)
#ordered by pvalue 
AA_eqtl_ordered_4 <- AA_eqtl_ordered_4[order(AA_eqtl_ordered_4$pvalue),]
##filter african american for population specific ones 
#subset eqtls not in europeans
AA_eqtl_ordered_specific_4 <- AA_eqtl_ordered_4[!(AA_eqtl_ordered_4$snpid %in% EUR_eqtl_ordered_4$snpid),]
#get names of eqtls ordered pvalue 
AA_eqtl_ordered_specific_names_4 <- AA_eqtl_ordered_specific_4[order(AA_eqtl_ordered_specific_4$pvalue),] 
#subset eqtls not in african americans and only in europeans
EUR_eqtl_ordered_specific_4 <- EUR_eqtl_ordered_4[!(EUR_eqtl_ordered_4$snpid %in% AA_eqtl_ordered_4$snpid),]

#Data for eQTLs number less than 0.00005
EUR_4 <- nrow(EUR_eqtl_ordered_4)
AA_4 <- nrow(AA_eqtl_ordered_4)
EUR_specific_4 <- nrow(EUR_eqtl_ordered_specific_4)
AA_specific_4 <- nrow(AA_eqtl_ordered_specific_names_4) 
shared_4 <- EUR_4 - EUR_specific_4




#AA_eqtl_ordered_only_0.00001 <- subset(x = AA_eqtl_ordered_only, subset = pvalue <= 0.00001)
write.csv(AA_eqtl_ordered_specific_names_4, "AA_eqtl_ordered_specific_names_4.csv")

#------------------------------------------------------------------------
#EUR eQTLs chromosome 5 

#------------------------------------------------------------------------

#read table
EUR_eqtl_5 <- read.table("EUR_eQTL_chr5.results.txt")
#add colnames
colnames(x = EUR_eqtl_5) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
EUR_eqtl_ordered_5 <- subset(x = EUR_eqtl_5, subset = pvalue <= 0.0000001)
#order p-value
EUR_eqtl_ordered_5 <- EUR_eqtl_ordered_5[order(EUR_eqtl_ordered_5$pvalue),]
#save SNP column 
EUR_eqtls_SNP_5 <- EUR_eqtl_ordered_5[,5, drop =FALSE]
#write csv for ordered EUR chr
#write.csv(EUR_eqtl_ordered_5, "EUR_eqtl_ordered_5.csv")
#extra commands 
chr5_MP_EUR = plot(-log10(x = EUR_eqtl_5$pvalue) , main='chr5')                 
#-----------------------------------------------------------------------------------------------------------------------------------------------
#AA eQTLs chromosome 5 ##threshold for 0.00001 through fastqtl 
#read table
AA_eqtl_5 <- read.table("AA_eQTL_chr5.results.txt")

#rename columns 
colnames(x = AA_eqtl_5) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
AA_eqtl_ordered_5 <- subset(x = AA_eqtl_5, subset = pvalue <= 0.0000001)
#ordered by pvalue 
AA_eqtl_ordered_5 <- AA_eqtl_ordered_5[order(AA_eqtl_ordered_5$pvalue),]
##filter african american for population specific ones 
#subset eqtls not in europeans
AA_eqtl_ordered_specific_5 <- AA_eqtl_ordered_5[!(AA_eqtl_ordered_5$snpid %in% EUR_eqtl_ordered_5$snpid),]
#get names of eqtls ordered pvalue 
AA_eqtl_ordered_specific_names_5 <- AA_eqtl_ordered_specific_5[order(AA_eqtl_ordered_specific_5$pvalue),] 
#subset eqtls not in african americans and only in europeans
EUR_eqtl_ordered_specific_5 <- EUR_eqtl_ordered_5[!(EUR_eqtl_ordered_5$snpid %in% AA_eqtl_ordered_5$snpid),]

#Data for eQTLs number less than 0.00005
EUR_5 <- nrow(EUR_eqtl_ordered_5)
AA_5 <- nrow(AA_eqtl_ordered_5)
EUR_specific_5 <- nrow(EUR_eqtl_ordered_specific_5)
AA_specific_5 <- nrow(AA_eqtl_ordered_specific_names_5) 
shared_5 <- EUR_5 - EUR_specific_5




#AA_eqtl_ordered_only_0.00001 <- subset(x = AA_eqtl_ordered_only, subset = pvalue <= 0.00001)
write.csv(AA_eqtl_ordered_specific_names_5, "AA_eqtl_ordered_specific_names_5.csv")

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#EUR eQTLs chromosome 6
#EUR eQTLs chromosome 6 

#------------------------------------------------------------------------

#read table
EUR_eqtl_6 <- read.table("EUR_eQTL_chr6.results.txt")
#add colnames
colnames(x = EUR_eqtl_6) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
EUR_eqtl_ordered_6 <- subset(x = EUR_eqtl_6, subset = pvalue <= 0.0000001)
#order p-value
EUR_eqtl_ordered_6 <- EUR_eqtl_ordered_6[order(EUR_eqtl_ordered_6$pvalue),]
#save SNP column 
EUR_eqtls_SNP_6 <- EUR_eqtl_ordered_6[,6, drop =FALSE]
#write csv for ordered EUR chr
#write.csv(EUR_eqtl_ordered_6, "EUR_eqtl_ordered_6.csv")
#extra commands 
chr6_MP_EUR = plot(-log10(x = EUR_eqtl_6$pvalue) , main='chr6')                 
#------------------------------------------------------------------------
#AA eQTLs chromosome 6 ##threshold for 0.00001 through fastqtl 
#read table
AA_eqtl_6 <- read.table("AA_eQTL_chr6.results.txt")

#rename columns 
colnames(x = AA_eqtl_6) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
AA_eqtl_ordered_6 <- subset(x = AA_eqtl_6, subset = pvalue <= 0.0000001)
#ordered by pvalue 
AA_eqtl_ordered_6 <- AA_eqtl_ordered_6[order(AA_eqtl_ordered_6$pvalue),]
##filter african american for population specific ones 
#subset eqtls not in europeans
AA_eqtl_ordered_specific_6 <- AA_eqtl_ordered_6[!(AA_eqtl_ordered_6$snpid %in% EUR_eqtl_ordered_6$snpid),]
#get names of eqtls ordered pvalue 
AA_eqtl_ordered_specific_names_6 <- AA_eqtl_ordered_specific_6[order(AA_eqtl_ordered_specific_6$pvalue),] 
#subset eqtls not in african americans and only in europeans
EUR_eqtl_ordered_specific_6 <- EUR_eqtl_ordered_6[!(EUR_eqtl_ordered_6$snpid %in% AA_eqtl_ordered_6$snpid),]

#Data for eQTLs number less than 0.00006
EUR_6 <- nrow(EUR_eqtl_ordered_6)
AA_6 <- nrow(AA_eqtl_ordered_6)
EUR_specific_6 <- nrow(EUR_eqtl_ordered_specific_6)
AA_specific_6 <- nrow(AA_eqtl_ordered_specific_names_6) 
shared_6 <- EUR_6 - EUR_specific_6




#AA_eqtl_ordered_only_0.00001 <- subset(x = AA_eqtl_ordered_only, subset = pvalue <= 0.00001)
write.csv(AA_eqtl_ordered_specific_names_6, "AA_eqtl_ordered_specific_names_6.csv")
#------------------------------------------------------------------------
#EUR eQTLs chromosome 7
#EUR eQTLs chromosome 7 

#------------------------------------------------------------------------

#read table
EUR_eqtl_7 <- read.table("EUR_eQTL_chr7.results.txt")
#add colnames
colnames(x = EUR_eqtl_7) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
EUR_eqtl_ordered_7 <- subset(x = EUR_eqtl_7, subset = pvalue <= 0.0000001)
#order p-value
EUR_eqtl_ordered_7 <- EUR_eqtl_ordered_7[order(EUR_eqtl_ordered_7$pvalue),]
#save SNP column 
EUR_eqtls_SNP_7 <- EUR_eqtl_ordered_7[,6, drop =FALSE]
#write csv for ordered EUR chr
#write.csv(EUR_eqtl_ordered_7, "EUR_eqtl_ordered_7.csv")
#extra commands 
chr7_MP_EUR = plot(-log10(x = EUR_eqtl_7$pvalue) , main='chr7')                 
#------------------------------------------------------------------------
#AA eQTLs chromosome 6 ##threshold for 0.00001 through fastqtl 
#read table
AA_eqtl_7 <- read.table("AA_eQTL_chr7.results.txt")

#rename columns 
colnames(x = AA_eqtl_7) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
AA_eqtl_ordered_7 <- subset(x = AA_eqtl_7, subset = pvalue <= 0.0000001)
#ordered by pvalue 
AA_eqtl_ordered_7 <- AA_eqtl_ordered_7[order(AA_eqtl_ordered_7$pvalue),]
##filter african american for population specific ones 
#subset eqtls not in europeans
AA_eqtl_ordered_specific_7 <- AA_eqtl_ordered_7[!(AA_eqtl_ordered_7$snpid %in% EUR_eqtl_ordered_7$snpid),]
#get names of eqtls ordered pvalue 
AA_eqtl_ordered_specific_names_7 <- AA_eqtl_ordered_specific_7[order(AA_eqtl_ordered_specific_7$pvalue),] 
#subset eqtls not in african americans and only in europeans
EUR_eqtl_ordered_specific_7 <- EUR_eqtl_ordered_7[!(EUR_eqtl_ordered_7$snpid %in% AA_eqtl_ordered_7$snpid),]

#Data for eQTLs number less than 0.00006
EUR_7 <- nrow(EUR_eqtl_ordered_7)
AA_7 <- nrow(AA_eqtl_ordered_7)
EUR_specific_7 <- nrow(EUR_eqtl_ordered_specific_7)
AA_specific_7 <- nrow(AA_eqtl_ordered_specific_names_7) 
shared_7 <- EUR_7 - EUR_specific_7




#AA_eqtl_ordered_only_0.00001 <- subset(x = AA_eqtl_ordered_only, subset = pvalue <= 0.00001)
write.csv(AA_eqtl_ordered_specific_names_7, "AA_eqtl_ordered_specific_names_7.csv")
#------------------------------------------------------------------------
#EUR eQTLs chromosome 8 
#EUR eQTLs chromosome 8 

#------------------------------------------------------------------------

#read table
EUR_eqtl_8 <- read.table("EUR_eQTL_chr8.results.txt")
#add colnames
colnames(x = EUR_eqtl_8) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
EUR_eqtl_ordered_8 <- subset(x = EUR_eqtl_8, subset = pvalue <= 0.0000001)
#order p-value
EUR_eqtl_ordered_8 <- EUR_eqtl_ordered_8[order(EUR_eqtl_ordered_8$pvalue),]
#save SNP column 
EUR_eqtls_SNP_8 <- EUR_eqtl_ordered_8[,6, drop =FALSE]
#write csv for ordered EUR chr
#write.csv(EUR_eqtl_ordered_8, "EUR_eqtl_ordered_8.csv")
#extra commands 
chr8_MP_EUR = plot(-log10(x = EUR_eqtl_8$pvalue) , main='chr8')                 
#------------------------------------------------------------------------
#AA eQTLs chromosome 8 ##threshold for 0.00001 through fastqtl 
#read table
AA_eqtl_8 <- read.table("AA_eQTL_chr8.results.txt")

#rename columns 
colnames(x = AA_eqtl_8) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
AA_eqtl_ordered_8 <- subset(x = AA_eqtl_8, subset = pvalue <= 0.0000001)
#ordered by pvalue 
AA_eqtl_ordered_8 <- AA_eqtl_ordered_8[order(AA_eqtl_ordered_8$pvalue),]
##filter african american for population specific ones 
#subset eqtls not in europeans
AA_eqtl_ordered_specific_8 <- AA_eqtl_ordered_8[!(AA_eqtl_ordered_8$snpid %in% EUR_eqtl_ordered_8$snpid),]
#get names of eqtls ordered pvalue 
AA_eqtl_ordered_specific_names_8 <- AA_eqtl_ordered_specific_8[order(AA_eqtl_ordered_specific_8$pvalue),] 
#subset eqtls not in african americans and only in europeans
EUR_eqtl_ordered_specific_8 <- EUR_eqtl_ordered_8[!(EUR_eqtl_ordered_8$snpid %in% AA_eqtl_ordered_8$snpid),]

#Data for eQTLs number less than 0.00006
EUR_8 <- nrow(EUR_eqtl_ordered_8)
AA_8 <- nrow(AA_eqtl_ordered_8)
EUR_specific_8 <- nrow(EUR_eqtl_ordered_specific_8)
AA_specific_8 <- nrow(AA_eqtl_ordered_specific_names_8) 
shared_8 <- EUR_8 - EUR_specific_8




#AA_eqtl_ordered_only_0.00001 <- subset(x = AA_eqtl_ordered_only, subset = pvalue <= 0.00001)
write.csv(AA_eqtl_ordered_specific_names_8, "AA_eqtl_ordered_specific_names_8.csv")

#------------------------------------------------------------------------#------------------------------------------------------------------------
#------------------------------------------------------------------------#------------------------------------------------------------------------
#EUR eQTLs chromosome 9 
#read table
EUR_eqtl_9 <- read.table("EUR_eQTL_chr9.results.txt")
#add colnames
colnames(x = EUR_eqtl_9) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
EUR_eqtl_ordered_9 <- subset(x = EUR_eqtl_9, subset = pvalue <= 0.0000001)
#order p-value
EUR_eqtl_ordered_9 <- EUR_eqtl_ordered_9[order(EUR_eqtl_ordered_9$pvalue),]
#save SNP column 
EUR_eqtls_SNP_9 <- EUR_eqtl_ordered_9[,6, drop =FALSE]
#write csv for ordered EUR chr
#write.csv(EUR_eqtl_ordered_9, "EUR_eqtl_ordered_9.csv")
#extra commands 
chr9_MP_EUR = plot(-log10(x = EUR_eqtl_9$pvalue) , main='chr9') 
#AA eQTLs chromosome 9 ##threshold for 0.00001 through fastqtl 
#read table
AA_eqtl_9 <- read.table("AA_eQTL_chr9.results.txt")

#rename columns 
colnames(x = AA_eqtl_9) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
AA_eqtl_ordered_9 <- subset(x = AA_eqtl_9, subset = pvalue <= 0.0000001)
#ordered by pvalue 
AA_eqtl_ordered_9 <- AA_eqtl_ordered_9[order(AA_eqtl_ordered_9$pvalue),]
##filter african american for population specific ones 
#subset eqtls not in europeans
AA_eqtl_ordered_specific_9 <- AA_eqtl_ordered_9[!(AA_eqtl_ordered_9$snpid %in% EUR_eqtl_ordered_9$snpid),]
#get names of eqtls ordered pvalue 
AA_eqtl_ordered_specific_names_9 <- AA_eqtl_ordered_specific_9[order(AA_eqtl_ordered_specific_9$pvalue),] 
#subset eqtls not in african americans and only in europeans
EUR_eqtl_ordered_specific_9 <- EUR_eqtl_ordered_9[!(EUR_eqtl_ordered_9$snpid %in% AA_eqtl_ordered_9$snpid),]

#Data for eQTLs number less than 0.00006
EUR_9 <- nrow(EUR_eqtl_ordered_9)
AA_9 <- nrow(AA_eqtl_ordered_9)
EUR_specific_9 <- nrow(EUR_eqtl_ordered_specific_9)
AA_specific_9 <- nrow(AA_eqtl_ordered_specific_names_9) 
shared_9 <- EUR_9 - EUR_specific_9




#AA_eqtl_ordered_only_0.00001 <- subset(x = AA_eqtl_ordered_only, subset = pvalue <= 0.00001)
write.csv(AA_eqtl_ordered_specific_names_9, "AA_eqtl_ordered_specific_names_9.csv")

#------------------------------------------------------------------------
#EUR eQTLs chromosome 10 
EUR_eqtl_10 <- read.table("EUR_eQTL_chr10.results.txt")
#add colnames
colnames(x = EUR_eqtl_10) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
EUR_eqtl_ordered_10 <- subset(x = EUR_eqtl_10, subset = pvalue <= 0.0000001)
#order p-value
EUR_eqtl_ordered_10 <- EUR_eqtl_ordered_10[order(EUR_eqtl_ordered_10$pvalue),]
#save SNP column 
EUR_eqtls_SNP_10 <- EUR_eqtl_ordered_10[,6, drop =FALSE]
#write csv for ordered EUR chr
#write.csv(EUR_eqtl_ordered_10, "EUR_eqtl_ordered_10.csv")
#extra commands 
chr10_MP_EUR = plot(-log10(x = EUR_eqtl_10$pvalue) , main='chr10') 
#AA eQTLs chromosome 10 ##threshold for 0.00001 through fastqtl 
#read table
AA_eqtl_10 <- read.table("AA_eQTL_chr10.results.txt")

#rename columns 
colnames(x = AA_eqtl_10) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
AA_eqtl_ordered_10 <- subset(x = AA_eqtl_10, subset = pvalue <= 0.0000001)
#ordered by pvalue 
AA_eqtl_ordered_10 <- AA_eqtl_ordered_10[order(AA_eqtl_ordered_10$pvalue),]
##filter african american for population specific ones 
#subset eqtls not in europeans
AA_eqtl_ordered_specific_10 <- AA_eqtl_ordered_10[!(AA_eqtl_ordered_10$snpid %in% EUR_eqtl_ordered_10$snpid),]
#get names of eqtls ordered pvalue 
AA_eqtl_ordered_specific_names_10 <- AA_eqtl_ordered_specific_10[order(AA_eqtl_ordered_specific_10$pvalue),] 
#subset eqtls not in african americans and only in europeans
EUR_eqtl_ordered_specific_10 <- EUR_eqtl_ordered_10[!(EUR_eqtl_ordered_10$snpid %in% AA_eqtl_ordered_10$snpid),]

#Data for eQTLs number less than 0.00006
EUR_10 <- nrow(EUR_eqtl_ordered_10)
AA_10 <- nrow(AA_eqtl_ordered_10)
EUR_specific_10 <- nrow(EUR_eqtl_ordered_specific_10)
AA_specific_10 <- nrow(AA_eqtl_ordered_specific_names_10) 
shared_10 <- EUR_10 - EUR_specific_10




#AA_eqtl_ordered_only_0.00001 <- subset(x = AA_eqtl_ordered_only, subset = pvalue <= 0.00001)
write.csv(AA_eqtl_ordered_specific_names_10, "AA_eqtl_ordered_specific_names_10.csv")
#------------------------------------------------------------------------
#EUR eQTLs chromosome 11

EUR_eqtl_11 <- read.table("EUR_eQTL_chr11.results.txt")
#add colnames
colnames(x = EUR_eqtl_11) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
EUR_eqtl_ordered_11 <- subset(x = EUR_eqtl_11, subset = pvalue <= 0.0000001)
#order p-value
EUR_eqtl_ordered_11 <- EUR_eqtl_ordered_11[order(EUR_eqtl_ordered_11$pvalue),]
#save SNP column 
EUR_eqtls_SNP_11 <- EUR_eqtl_ordered_11[,6, drop =FALSE]
#write csv for ordered EUR chr
#write.csv(EUR_eqtl_ordered_11, "EUR_eqtl_ordered_11.csv")
#extra commands 
chr11_MP_EUR = plot(-log11(x = EUR_eqtl_11$pvalue) , main='chr11') 
#AA eQTLs chromosome 11 ##threshold for 0.00001 through fastqtl 
#read table
AA_eqtl_11 <- read.table("AA_eQTL_chr11.results.txt")

#rename columns 
colnames(x = AA_eqtl_11) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
AA_eqtl_ordered_11 <- subset(x = AA_eqtl_11, subset = pvalue <= 0.0000001)
#ordered by pvalue 
AA_eqtl_ordered_11 <- AA_eqtl_ordered_11[order(AA_eqtl_ordered_11$pvalue),]
##filter african american for population specific ones 
#subset eqtls not in europeans
AA_eqtl_ordered_specific_11 <- AA_eqtl_ordered_11[!(AA_eqtl_ordered_11$snpid %in% EUR_eqtl_ordered_11$snpid),]
#get names of eqtls ordered pvalue 
AA_eqtl_ordered_specific_names_11 <- AA_eqtl_ordered_specific_11[order(AA_eqtl_ordered_specific_11$pvalue),] 
#subset eqtls not in african americans and only in europeans
EUR_eqtl_ordered_specific_11 <- EUR_eqtl_ordered_11[!(EUR_eqtl_ordered_11$snpid %in% AA_eqtl_ordered_11$snpid),]

#Data for eQTLs number less than 0.00006
EUR_11 <- nrow(EUR_eqtl_ordered_11)
AA_11 <- nrow(AA_eqtl_ordered_11)
EUR_specific_11 <- nrow(EUR_eqtl_ordered_specific_11)
AA_specific_11 <- nrow(AA_eqtl_ordered_specific_names_11) 
shared_11 <- EUR_11 - EUR_specific_11




#AA_eqtl_ordered_only_0.00001 <- subset(x = AA_eqtl_ordered_only, subset = pvalue <= 0.00001)
write.csv(AA_eqtl_ordered_specific_names_11, "AA_eqtl_ordered_specific_names_11.csv")
#------------------------------------------------------------------------
#EUR eQTLs chromosome 12
EUR_eqtl_12 <- read.table("EUR_eQTL_chr12.results.txt")
#add colnames
colnames(x = EUR_eqtl_12) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
EUR_eqtl_ordered_12 <- subset(x = EUR_eqtl_12, subset = pvalue <= 0.0000001)
#order p-value
EUR_eqtl_ordered_12 <- EUR_eqtl_ordered_12[order(EUR_eqtl_ordered_12$pvalue),]
#save SNP column 
EUR_eqtls_SNP_12 <- EUR_eqtl_ordered_12[,6, drop =FALSE]
#write csv for ordered EUR chr
#write.csv(EUR_eqtl_ordered_12, "EUR_eqtl_ordered_12.csv")
#extra commands 
chr12_MP_EUR = plot(-log12(x = EUR_eqtl_12$pvalue) , main='chr12') 
#AA eQTLs chromosome 12 ##threshold for 0.00001 through fastqtl 
#read table
AA_eqtl_12 <- read.table("AA_eQTL_chr12.results.txt")

#rename columns 
colnames(x = AA_eqtl_12) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
AA_eqtl_ordered_12 <- subset(x = AA_eqtl_12, subset = pvalue <= 0.0000001)
#ordered by pvalue 
AA_eqtl_ordered_12 <- AA_eqtl_ordered_12[order(AA_eqtl_ordered_12$pvalue),]
##filter african american for population specific ones 
#subset eqtls not in europeans
AA_eqtl_ordered_specific_12 <- AA_eqtl_ordered_12[!(AA_eqtl_ordered_12$snpid %in% EUR_eqtl_ordered_12$snpid),]
#get names of eqtls ordered pvalue 
AA_eqtl_ordered_specific_names_12 <- AA_eqtl_ordered_specific_12[order(AA_eqtl_ordered_specific_12$pvalue),] 
#subset eqtls not in african americans and only in europeans
EUR_eqtl_ordered_specific_12 <- EUR_eqtl_ordered_12[!(EUR_eqtl_ordered_12$snpid %in% AA_eqtl_ordered_12$snpid),]

#Data for eQTLs number less than 0.00006
EUR_12 <- nrow(EUR_eqtl_ordered_12)
AA_12 <- nrow(AA_eqtl_ordered_12)
EUR_specific_12 <- nrow(EUR_eqtl_ordered_specific_12)
AA_specific_12 <- nrow(AA_eqtl_ordered_specific_names_12) 
shared_12 <- EUR_12 - EUR_specific_12




#AA_eqtl_ordered_only_0.00001 <- subset(x = AA_eqtl_ordered_only, subset = pvalue <= 0.00001)
write.csv(AA_eqtl_ordered_specific_names_12, "AA_eqtl_ordered_specific_names_12.csv")
#------------------------------------------------------------------------
#EUR eQTLs chromosome 13
EUR_eqtl_13 <- read.table("EUR_eQTL_chr13.results.txt")
#add colnames
colnames(x = EUR_eqtl_13) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
EUR_eqtl_ordered_13 <- subset(x = EUR_eqtl_13, subset = pvalue <= 0.0000001)
#order p-value
EUR_eqtl_ordered_13 <- EUR_eqtl_ordered_13[order(EUR_eqtl_ordered_13$pvalue),]
#save SNP column 
EUR_eqtls_SNP_13 <- EUR_eqtl_ordered_13[,6, drop =FALSE]
#write csv for ordered EUR chr
#write.csv(EUR_eqtl_ordered_13, "EUR_eqtl_ordered_13.csv")
#extra commands 
chr13_MP_EUR = plot(-log13(x = EUR_eqtl_13$pvalue) , main='chr13') 
#AA eQTLs chromosome 13 ##threshold for 0.00001 through fastqtl 
#read table
AA_eqtl_13 <- read.table("AA_eQTL_chr13.results.txt")

#rename columns 
colnames(x = AA_eqtl_13) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
AA_eqtl_ordered_13 <- subset(x = AA_eqtl_13, subset = pvalue <= 0.0000001)
#ordered by pvalue 
AA_eqtl_ordered_13 <- AA_eqtl_ordered_13[order(AA_eqtl_ordered_13$pvalue),]
##filter african american for population specific ones 
#subset eqtls not in europeans
AA_eqtl_ordered_specific_13 <- AA_eqtl_ordered_13[!(AA_eqtl_ordered_13$snpid %in% EUR_eqtl_ordered_13$snpid),]
#get names of eqtls ordered pvalue 
AA_eqtl_ordered_specific_names_13 <- AA_eqtl_ordered_specific_13[order(AA_eqtl_ordered_specific_13$pvalue),] 
#subset eqtls not in african americans and only in europeans
EUR_eqtl_ordered_specific_13 <- EUR_eqtl_ordered_13[!(EUR_eqtl_ordered_13$snpid %in% AA_eqtl_ordered_13$snpid),]

#Data for eQTLs number less than 0.00006
EUR_13 <- nrow(EUR_eqtl_ordered_13)
AA_13 <- nrow(AA_eqtl_ordered_13)
EUR_specific_13 <- nrow(EUR_eqtl_ordered_specific_13)
AA_specific_13 <- nrow(AA_eqtl_ordered_specific_names_13) 
shared_13 <- EUR_13 - EUR_specific_13




#AA_eqtl_ordered_only_0.00001 <- subset(x = AA_eqtl_ordered_only, subset = pvalue <= 0.00001)
write.csv(AA_eqtl_ordered_specific_names_13, "AA_eqtl_ordered_specific_names_13.csv")
#------------------------------------------------------------------------
#EUR eQTLs chromosome 14
EUR_eqtl_14 <- read.table("EUR_eQTL_chr14.results.txt")
#add colnames
colnames(x = EUR_eqtl_14) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
EUR_eqtl_ordered_14 <- subset(x = EUR_eqtl_14, subset = pvalue <= 0.0000001)
#order p-value
EUR_eqtl_ordered_14 <- EUR_eqtl_ordered_14[order(EUR_eqtl_ordered_14$pvalue),]
#save SNP column 
EUR_eqtls_SNP_14 <- EUR_eqtl_ordered_14[,6, drop =FALSE]
#write csv for ordered EUR chr
#write.csv(EUR_eqtl_ordered_14, "EUR_eqtl_ordered_14.csv")
#extra commands 
chr14_MP_EUR = plot(-log14(x = EUR_eqtl_14$pvalue) , main='chr14') 
#AA eQTLs chromosome 14 ##threshold for 0.00001 through fastqtl 
#read table
AA_eqtl_14 <- read.table("AA_eQTL_chr14.results.txt")

#rename columns 
colnames(x = AA_eqtl_14) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
AA_eqtl_ordered_14 <- subset(x = AA_eqtl_14, subset = pvalue <= 0.0000001)
#ordered by pvalue 
AA_eqtl_ordered_14 <- AA_eqtl_ordered_14[order(AA_eqtl_ordered_14$pvalue),]
##filter african american for population specific ones 
#subset eqtls not in europeans
AA_eqtl_ordered_specific_14 <- AA_eqtl_ordered_14[!(AA_eqtl_ordered_14$snpid %in% EUR_eqtl_ordered_14$snpid),]
#get names of eqtls ordered pvalue 
AA_eqtl_ordered_specific_names_14 <- AA_eqtl_ordered_specific_14[order(AA_eqtl_ordered_specific_14$pvalue),] 
#subset eqtls not in african americans and only in europeans
EUR_eqtl_ordered_specific_14 <- EUR_eqtl_ordered_14[!(EUR_eqtl_ordered_14$snpid %in% AA_eqtl_ordered_14$snpid),]

#Data for eQTLs number less than 0.00006
EUR_14 <- nrow(EUR_eqtl_ordered_14)
AA_14 <- nrow(AA_eqtl_ordered_14)
EUR_specific_14 <- nrow(EUR_eqtl_ordered_specific_14)
AA_specific_14 <- nrow(AA_eqtl_ordered_specific_names_14) 
shared_14 <- EUR_14 - EUR_specific_14




#AA_eqtl_ordered_only_0.00001 <- subset(x = AA_eqtl_ordered_only, subset = pvalue <= 0.00001)
write.csv(AA_eqtl_ordered_specific_names_14, "AA_eqtl_ordered_specific_names_14.csv")
#------------------------------------------------------------------------
#EUR eQTLs chromosome 15
EUR_eqtl_15 <- read.table("EUR_eQTL_chr15.results.txt")
#add colnames
colnames(x = EUR_eqtl_15) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
EUR_eqtl_ordered_15 <- subset(x = EUR_eqtl_15, subset = pvalue <= 0.0000001)
#order p-value
EUR_eqtl_ordered_15 <- EUR_eqtl_ordered_15[order(EUR_eqtl_ordered_15$pvalue),]
#save SNP column 
EUR_eqtls_SNP_15 <- EUR_eqtl_ordered_15[,6, drop =FALSE]
#write csv for ordered EUR chr
#write.csv(EUR_eqtl_ordered_15, "EUR_eqtl_ordered_15.csv")
#extra commands 
chr15_MP_EUR = plot(-log15(x = EUR_eqtl_15$pvalue) , main='chr15') 
#AA eQTLs chromosome 15 ##threshold for 0.00001 through fastqtl 
#read table
AA_eqtl_15 <- read.table("AA_eQTL_chr15.results.txt")

#rename columns 
colnames(x = AA_eqtl_15) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
AA_eqtl_ordered_15 <- subset(x = AA_eqtl_15, subset = pvalue <= 0.0000001)
#ordered by pvalue 
AA_eqtl_ordered_15 <- AA_eqtl_ordered_15[order(AA_eqtl_ordered_15$pvalue),]
##filter african american for population specific ones 
#subset eqtls not in europeans
AA_eqtl_ordered_specific_15 <- AA_eqtl_ordered_15[!(AA_eqtl_ordered_15$snpid %in% EUR_eqtl_ordered_15$snpid),]
#get names of eqtls ordered pvalue 
AA_eqtl_ordered_specific_names_15 <- AA_eqtl_ordered_specific_15[order(AA_eqtl_ordered_specific_15$pvalue),] 
#subset eqtls not in african americans and only in europeans
EUR_eqtl_ordered_specific_15 <- EUR_eqtl_ordered_15[!(EUR_eqtl_ordered_15$snpid %in% AA_eqtl_ordered_15$snpid),]

#Data for eQTLs number less than 0.00006
EUR_15 <- nrow(EUR_eqtl_ordered_15)
AA_15 <- nrow(AA_eqtl_ordered_15)
EUR_specific_15 <- nrow(EUR_eqtl_ordered_specific_15)
AA_specific_15 <- nrow(AA_eqtl_ordered_specific_names_15) 
shared_15 <- EUR_15 - EUR_specific_15




#AA_eqtl_ordered_only_0.00001 <- subset(x = AA_eqtl_ordered_only, subset = pvalue <= 0.00001)
write.csv(AA_eqtl_ordered_specific_names_15, "AA_eqtl_ordered_specific_names_15.csv")
#------------------------------------------------------------------------
#EUR eQTLs chromosome 16
EUR_eqtl_16 <- read.table("EUR_eQTL_chr16.results.txt")
#add colnames
colnames(x = EUR_eqtl_16) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
EUR_eqtl_ordered_16 <- subset(x = EUR_eqtl_16, subset = pvalue <= 0.0000001)
#order p-value
EUR_eqtl_ordered_16 <- EUR_eqtl_ordered_16[order(EUR_eqtl_ordered_16$pvalue),]
#save SNP column 
EUR_eqtls_SNP_16 <- EUR_eqtl_ordered_16[,6, drop =FALSE]
#write csv for ordered EUR chr
#write.csv(EUR_eqtl_ordered_16, "EUR_eqtl_ordered_16.csv")
#extra commands 
chr16_MP_EUR = plot(-log16(x = EUR_eqtl_16$pvalue) , main='chr16') 
#AA eQTLs chromosome 16 ##threshold for 0.00001 through fastqtl 
#read table
AA_eqtl_16 <- read.table("AA_eQTL_chr16.results.txt")

#rename columns 
colnames(x = AA_eqtl_16) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
AA_eqtl_ordered_16 <- subset(x = AA_eqtl_16, subset = pvalue <= 0.0000001)
#ordered by pvalue 
AA_eqtl_ordered_16 <- AA_eqtl_ordered_16[order(AA_eqtl_ordered_16$pvalue),]
##filter african american for population specific ones 
#subset eqtls not in europeans
AA_eqtl_ordered_specific_16 <- AA_eqtl_ordered_16[!(AA_eqtl_ordered_16$snpid %in% EUR_eqtl_ordered_16$snpid),]
#get names of eqtls ordered pvalue 
AA_eqtl_ordered_specific_names_16 <- AA_eqtl_ordered_specific_16[order(AA_eqtl_ordered_specific_16$pvalue),] 
#subset eqtls not in african americans and only in europeans
EUR_eqtl_ordered_specific_16 <- EUR_eqtl_ordered_16[!(EUR_eqtl_ordered_16$snpid %in% AA_eqtl_ordered_16$snpid),]

#Data for eQTLs number less than 0.00006
EUR_16 <- nrow(EUR_eqtl_ordered_16)
AA_16 <- nrow(AA_eqtl_ordered_16)
EUR_specific_16 <- nrow(EUR_eqtl_ordered_specific_16)
AA_specific_16 <- nrow(AA_eqtl_ordered_specific_names_16) 
shared_16 <- EUR_16 - EUR_specific_16




#AA_eqtl_ordered_only_0.00001 <- subset(x = AA_eqtl_ordered_only, subset = pvalue <= 0.00001)
write.csv(AA_eqtl_ordered_specific_names_16, "AA_eqtl_ordered_specific_names_16.csv")
#------------------------------------------------------------------------
#EUR eQTLs chromosome 17
EUR_eqtl_17 <- read.table("EUR_eQTL_chr17.results.txt")
#add colnames
colnames(x = EUR_eqtl_17) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
EUR_eqtl_ordered_17 <- subset(x = EUR_eqtl_17, subset = pvalue <= 0.0000001)
#order p-value
EUR_eqtl_ordered_17 <- EUR_eqtl_ordered_17[order(EUR_eqtl_ordered_17$pvalue),]
#save SNP column 
EUR_eqtls_SNP_17 <- EUR_eqtl_ordered_17[,6, drop =FALSE]
#write csv for ordered EUR chr
#write.csv(EUR_eqtl_ordered_17, "EUR_eqtl_ordered_17.csv")
#extra commands 
chr17_MP_EUR = plot(-log10(x = EUR_eqtl_17$pvalue) , main='chr17') 
#AA eQTLs chromosome 17 ##threshold for 0.00001 through fastqtl 
#read table
AA_eqtl_17 <- read.table("AA_eQTL_chr17.results.txt")

#rename columns 
colnames(x = AA_eqtl_17) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
AA_eqtl_ordered_17 <- subset(x = AA_eqtl_17, subset = pvalue <= 0.0000001)
#ordered by pvalue 
AA_eqtl_ordered_17 <- AA_eqtl_ordered_17[order(AA_eqtl_ordered_17$pvalue),]
##filter african american for population specific ones 
#subset eqtls not in europeans
AA_eqtl_ordered_specific_17 <- AA_eqtl_ordered_17[!(AA_eqtl_ordered_17$snpid %in% EUR_eqtl_ordered_17$snpid),]
#get names of eqtls ordered pvalue 
AA_eqtl_ordered_specific_names_17 <- AA_eqtl_ordered_specific_17[order(AA_eqtl_ordered_specific_17$pvalue),] 
#subset eqtls not in african americans and only in europeans
EUR_eqtl_ordered_specific_17 <- EUR_eqtl_ordered_17[!(EUR_eqtl_ordered_17$snpid %in% AA_eqtl_ordered_17$snpid),]

#Data for eQTLs number less than 0.00006
EUR_17 <- nrow(EUR_eqtl_ordered_17)
AA_17 <- nrow(AA_eqtl_ordered_17)
EUR_specific_17 <- nrow(EUR_eqtl_ordered_specific_17)
AA_specific_17 <- nrow(AA_eqtl_ordered_specific_names_17) 
shared_17 <- EUR_17 - EUR_specific_17




#AA_eqtl_ordered_only_0.00001 <- subset(x = AA_eqtl_ordered_only, subset = pvalue <= 0.00001)
write.csv(AA_eqtl_ordered_specific_names_17, "AA_eqtl_ordered_specific_names_17.csv")

#------------------------------------------------------------------------
#EUR eQTLs chromosome 18
EUR_eqtl_18 <- read.table("EUR_eQTL_chr18.results.txt")
#add colnames
colnames(x = EUR_eqtl_18) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
EUR_eqtl_ordered_18 <- subset(x = EUR_eqtl_18, subset = pvalue <= 0.0000001)
#order p-value
EUR_eqtl_ordered_18 <- EUR_eqtl_ordered_18[order(EUR_eqtl_ordered_18$pvalue),]
#save SNP column 
EUR_eqtls_SNP_18 <- EUR_eqtl_ordered_18[,6, drop =FALSE]
#write csv for ordered EUR chr
#write.csv(EUR_eqtl_ordered_18, "EUR_eqtl_ordered_18.csv")
#extra commands 
chr18_MP_EUR = plot(-log10(x = EUR_eqtl_18$pvalue) , main='chr18') 
#AA eQTLs chromosome 18 ##threshold for 0.00001 through fastqtl 
#read table
AA_eqtl_18 <- read.table("AA_eQTL_chr18.results.txt")

#rename columns 
colnames(x = AA_eqtl_18) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
AA_eqtl_ordered_18 <- subset(x = AA_eqtl_18, subset = pvalue <= 0.0000001)
#ordered by pvalue 
AA_eqtl_ordered_18 <- AA_eqtl_ordered_18[order(AA_eqtl_ordered_18$pvalue),]
##filter african american for population specific ones 
#subset eqtls not in europeans
AA_eqtl_ordered_specific_18 <- AA_eqtl_ordered_18[!(AA_eqtl_ordered_18$snpid %in% EUR_eqtl_ordered_18$snpid),]
#get names of eqtls ordered pvalue 
AA_eqtl_ordered_specific_names_18 <- AA_eqtl_ordered_specific_18[order(AA_eqtl_ordered_specific_18$pvalue),] 
#subset eqtls not in african americans and only in europeans
EUR_eqtl_ordered_specific_18 <- EUR_eqtl_ordered_18[!(EUR_eqtl_ordered_18$snpid %in% AA_eqtl_ordered_18$snpid),]

#Data for eQTLs number less than 0.00006
EUR_18 <- nrow(EUR_eqtl_ordered_18)
AA_18 <- nrow(AA_eqtl_ordered_18)
EUR_specific_18 <- nrow(EUR_eqtl_ordered_specific_18)
AA_specific_18 <- nrow(AA_eqtl_ordered_specific_names_18) 
shared_18 <- EUR_18 - EUR_specific_18




#AA_eqtl_ordered_only_0.00001 <- subset(x = AA_eqtl_ordered_only, subset = pvalue <= 0.00001)
write.csv(AA_eqtl_ordered_specific_names_18, "AA_eqtl_ordered_specific_names_18.csv")
#------------------------------------------------------------------------
#EUR eQTLs chromosome 19
EUR_eqtl_19 <- read.table("EUR_eQTL_chr19.results.txt")
#add colnames
colnames(x = EUR_eqtl_19) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
EUR_eqtl_ordered_19 <- subset(x = EUR_eqtl_19, subset = pvalue <= 0.0000001)
#order p-value
EUR_eqtl_ordered_19 <- EUR_eqtl_ordered_19[order(EUR_eqtl_ordered_19$pvalue),]
#save SNP column 
EUR_eqtls_SNP_19 <- EUR_eqtl_ordered_19[,6, drop =FALSE]
#write csv for ordered EUR chr
#write.csv(EUR_eqtl_ordered_19, "EUR_eqtl_ordered_19.csv")
#extra commands 
chr19_MP_EUR = plot(-log10(x = EUR_eqtl_19$pvalue) , main='chr19') 
#AA eQTLs chromosome 19 ##threshold for 0.00001 through fastqtl 
#read table
AA_eqtl_19 <- read.table("AA_eQTL_chr19.results.txt")

#rename columns 
colnames(x = AA_eqtl_19) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
AA_eqtl_ordered_19 <- subset(x = AA_eqtl_19, subset = pvalue <= 0.0000001)
#ordered by pvalue 
AA_eqtl_ordered_19 <- AA_eqtl_ordered_19[order(AA_eqtl_ordered_19$pvalue),]
##filter african american for population specific ones 
#subset eqtls not in europeans
AA_eqtl_ordered_specific_19 <- AA_eqtl_ordered_19[!(AA_eqtl_ordered_19$snpid %in% EUR_eqtl_ordered_19$snpid),]
#get names of eqtls ordered pvalue 
AA_eqtl_ordered_specific_names_19 <- AA_eqtl_ordered_specific_19[order(AA_eqtl_ordered_specific_19$pvalue),] 
#subset eqtls not in african americans and only in europeans
EUR_eqtl_ordered_specific_19 <- EUR_eqtl_ordered_19[!(EUR_eqtl_ordered_19$snpid %in% AA_eqtl_ordered_19$snpid),]

#Data for eQTLs number less than 0.00006
EUR_19 <- nrow(EUR_eqtl_ordered_19)
AA_19 <- nrow(AA_eqtl_ordered_19)
EUR_specific_19 <- nrow(EUR_eqtl_ordered_specific_19)
AA_specific_19 <- nrow(AA_eqtl_ordered_specific_names_19) 
shared_19 <- EUR_19 - EUR_specific_19




#AA_eqtl_ordered_only_0.00001 <- subset(x = AA_eqtl_ordered_only, subset = pvalue <= 0.00001)
write.csv(AA_eqtl_ordered_specific_names_19, "AA_eqtl_ordered_specific_names_19.csv")
#------------------------------------------------------------------------
#EUR eQTLs chromosome 20
EUR_eqtl_20 <- read.table("EUR_eQTL_chr20.results.txt")
#add colnames
colnames(x = EUR_eqtl_20) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
EUR_eqtl_ordered_20 <- subset(x = EUR_eqtl_20, subset = pvalue <= 0.0000001)
#order p-value
EUR_eqtl_ordered_20 <- EUR_eqtl_ordered_20[order(EUR_eqtl_ordered_20$pvalue),]
#save SNP column 
EUR_eqtls_SNP_20 <- EUR_eqtl_ordered_20[,6, drop =FALSE]
#write csv for ordered EUR chr
#write.csv(EUR_eqtl_ordered_20, "EUR_eqtl_ordered_20.csv")
#extra commands 
chr20_MP_EUR = plot(-log10(x = EUR_eqtl_20$pvalue) , main='chr20') 
#AA eQTLs chromosome 20 ##threshold for 0.00001 through fastqtl 
#read table
AA_eqtl_20 <- read.table("AA_eQTL_chr20.results.txt")

#rename columns 
colnames(x = AA_eqtl_20) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
AA_eqtl_ordered_20 <- subset(x = AA_eqtl_20, subset = pvalue <= 0.0000001)
#ordered by pvalue 
AA_eqtl_ordered_20 <- AA_eqtl_ordered_20[order(AA_eqtl_ordered_20$pvalue),]
##filter african american for population specific ones 
#subset eqtls not in europeans
AA_eqtl_ordered_specific_20 <- AA_eqtl_ordered_20[!(AA_eqtl_ordered_20$snpid %in% EUR_eqtl_ordered_20$snpid),]
#get names of eqtls ordered pvalue 
AA_eqtl_ordered_specific_names_20 <- AA_eqtl_ordered_specific_20[order(AA_eqtl_ordered_specific_20$pvalue),] 
#subset eqtls not in african americans and only in europeans
EUR_eqtl_ordered_specific_20 <- EUR_eqtl_ordered_20[!(EUR_eqtl_ordered_20$snpid %in% AA_eqtl_ordered_20$snpid),]

#Data for eQTLs number less than 0.00006
EUR_20 <- nrow(EUR_eqtl_ordered_20)
AA_20 <- nrow(AA_eqtl_ordered_20)
EUR_specific_20 <- nrow(EUR_eqtl_ordered_specific_20)
AA_specific_20 <- nrow(AA_eqtl_ordered_specific_names_20) 
shared_20 <- EUR_20 - EUR_specific_20




#AA_eqtl_ordered_only_0.00001 <- subset(x = AA_eqtl_ordered_only, subset = pvalue <= 0.00001)
write.csv(AA_eqtl_ordered_specific_names_20, "AA_eqtl_ordered_specific_names_20.csv")


#------------------------------------------------------------------------
#eQTLs chromosome 21
EUR_eqtl_21 <- read.table("EUR_eQTL_chr21.results.txt")
#add colnames
colnames(x = EUR_eqtl_21) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
EUR_eqtl_ordered_21 <- subset(x = EUR_eqtl_21, subset = pvalue <= 0.0000001)
#order p-value
EUR_eqtl_ordered_21 <- EUR_eqtl_ordered_21[order(EUR_eqtl_ordered_21$pvalue),]
#save SNP column 
EUR_eqtls_SNP_21 <- EUR_eqtl_ordered_21[,6, drop =FALSE]
#write csv for ordered EUR chr
#write.csv(EUR_eqtl_ordered_21, "EUR_eqtl_ordered_21.csv")
#extra commands 
chr21_MP_EUR = plot(-log10(x = EUR_eqtl_21$pvalue) , main='chr21') 
#AA eQTLs chromosome 21 ##threshold for 0.00001 through fastqtl 
#read table
AA_eqtl_21 <- read.table("AA_eQTL_chr21.results.txt")

#rename columns 
colnames(x = AA_eqtl_21) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
AA_eqtl_ordered_21 <- subset(x = AA_eqtl_21, subset = pvalue <= 0.0000001)
#ordered by pvalue 
AA_eqtl_ordered_21 <- AA_eqtl_ordered_21[order(AA_eqtl_ordered_21$pvalue),]
##filter african american for population specific ones 
#subset eqtls not in europeans
AA_eqtl_ordered_specific_21 <- AA_eqtl_ordered_21[!(AA_eqtl_ordered_21$snpid %in% EUR_eqtl_ordered_21$snpid),]
#get names of eqtls ordered pvalue 
AA_eqtl_ordered_specific_names_21 <- AA_eqtl_ordered_specific_21[order(AA_eqtl_ordered_specific_21$pvalue),] 
#subset eqtls not in african americans and only in europeans
EUR_eqtl_ordered_specific_21 <- EUR_eqtl_ordered_21[!(EUR_eqtl_ordered_21$snpid %in% AA_eqtl_ordered_21$snpid),]

#Data for eQTLs number less than 0.00006
EUR_21 <- nrow(EUR_eqtl_ordered_21)
AA_21 <- nrow(AA_eqtl_ordered_21)
EUR_specific_21 <- nrow(EUR_eqtl_ordered_specific_21)
AA_specific_21 <- nrow(AA_eqtl_ordered_specific_names_21) 
shared_21 <- EUR_21 - EUR_specific_21




#AA_eqtl_ordered_only_0.00001 <- subset(x = AA_eqtl_ordered_only, subset = pvalue <= 0.00001)
write.csv(AA_eqtl_ordered_specific_names_21, "AA_eqtl_ordered_specific_names_21.csv")



#------------------------------------------------------------------------
#eQTLs chromosome 22
EUR_eqtl_22 <- read.table("EUR_eQTL_chr22.results.txt")
#add colnames
colnames(x = EUR_eqtl_22) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
EUR_eqtl_ordered_22 <- subset(x = EUR_eqtl_22, subset = pvalue <= 0.0000001)
#order p-value
EUR_eqtl_ordered_22 <- EUR_eqtl_ordered_22[order(EUR_eqtl_ordered_22$pvalue),]
#save SNP column 
EUR_eqtls_SNP_22 <- EUR_eqtl_ordered_22[,6, drop =FALSE]
#write csv for ordered EUR chr
#write.csv(EUR_eqtl_ordered_22, "EUR_eqtl_ordered_22.csv")
#extra commands 
chr22_MP_EUR = plot(-log10(x = EUR_eqtl_22$pvalue) , main='chr22') 
#AA eQTLs chromosome 22 ##threshold for 0.00001 through fastqtl 
#read table
AA_eqtl_22 <- read.table("AA_eQTL_chr22.results.txt")

#rename columns 
colnames(x = AA_eqtl_22) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
AA_eqtl_ordered_22 <- subset(x = AA_eqtl_22, subset = pvalue <= 0.0000001)
#ordered by pvalue 
AA_eqtl_ordered_22 <- AA_eqtl_ordered_22[order(AA_eqtl_ordered_22$pvalue),]
##filter african american for population specific ones 
#subset eqtls not in europeans
AA_eqtl_ordered_specific_22 <- AA_eqtl_ordered_22[!(AA_eqtl_ordered_22$snpid %in% EUR_eqtl_ordered_22$snpid),]
#get names of eqtls ordered pvalue 
AA_eqtl_ordered_specific_names_22 <- AA_eqtl_ordered_specific_22[order(AA_eqtl_ordered_specific_22$pvalue),] 
#subset eqtls not in african americans and only in europeans
EUR_eqtl_ordered_specific_22 <- EUR_eqtl_ordered_22[!(EUR_eqtl_ordered_22$snpid %in% AA_eqtl_ordered_22$snpid),]

#Data for eQTLs number less than 0.00006
EUR_22 <- nrow(EUR_eqtl_ordered_22)
AA_22 <- nrow(AA_eqtl_ordered_22)
EUR_specific_22 <- nrow(EUR_eqtl_ordered_specific_22)
AA_specific_22 <- nrow(AA_eqtl_ordered_specific_names_22) 
shared_22 <- EUR_22 - EUR_specific_22




#AA_eqtl_ordered_only_0.00001 <- subset(x = AA_eqtl_ordered_only, subset = pvalue <= 0.00001)
write.csv(AA_eqtl_ordered_specific_names_22, "AA_eqtl_ordered_specific_names_22.csv")
grid.newpage()
g = draw.pairwise.venn(1234, 1451, 895, category = c("AA", "EUR"), lty = rep("blank", 2), fill = c("snow2", "blue"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))
#add title 
require(gridExtra)
grid.arrange(gTree(children=g), top='Whole Blood eQTLs Chromosome 1')





AA_eqtl <- read.table("AA_eQTL_chr2.results.txt")
colnames(x = AA_eqtl) <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
AA_eqtl_0.00001 <- subset(x = AA_eqtl, subset = pvalue <= 0.00001)


chr1_MP_AA = plot(x= AA_eqtl_snpid$snpid, y = -log10(x = AA_eqtl_ordered$pvalue))

write.csv(AA_eqtl, "AA_eqtl.csv")

AA_eqtl_ordered <- AA_eqtl_ordered_specific_1[order(AA_eqtl_ordered_specific_1$pvalue),]
AA_eqtl_snpid <- AA_eqtl_ordered_specific_1[order(AA_eqtl_ordered_specific_1$snpid),]
write.csv(AA_eqtl_ordered, "AA_eqtl_ordered.csv")

########----------------save snp names 
EUR_eqtls_SNP_1<- EUR_eqtl_ordered_1[,c(2,4), drop =FALSE]
EUR_eqtls_SNP_2<- EUR_eqtl_ordered_2[,c(2,4), drop =FALSE]
EUR_eqtls_SNP_3<- EUR_eqtl_ordered_3[,c(2,4), drop =FALSE]
EUR_eqtls_SNP_4<- EUR_eqtl_ordered_4[,c(2,4), drop =FALSE]
EUR_eqtls_SNP_5<- EUR_eqtl_ordered_5[,c(2,4), drop =FALSE]
EUR_eqtls_SNP_6<- EUR_eqtl_ordered_6[,c(2,4), drop =FALSE]
EUR_eqtls_SNP_7<- EUR_eqtl_ordered_7[,c(2,4), drop =FALSE]
EUR_eqtls_SNP_8<- EUR_eqtl_ordered_8[,c(2,4), drop =FALSE]
EUR_eqtls_SNP_9<- EUR_eqtl_ordered_9[,c(2,4), drop =FALSE]
EUR_eqtls_SNP_10<- EUR_eqtl_ordered_10[,c(2,4), drop =FALSE]
EUR_eqtls_SNP_11<- EUR_eqtl_ordered_11[,c(2,4), drop =FALSE]
EUR_eqtls_SNP_12<- EUR_eqtl_ordered_12[,c(2,4), drop =FALSE]
EUR_eqtls_SNP_13<- EUR_eqtl_ordered_13[,c(2,4), drop =FALSE]
EUR_eqtls_SNP_14<- EUR_eqtl_ordered_14[,c(2,4), drop =FALSE]
EUR_eqtls_SNP_15<- EUR_eqtl_ordered_15[,c(2,4), drop =FALSE]
EUR_eqtls_SNP_16<- EUR_eqtl_ordered_16[,c(2,4), drop =FALSE]
EUR_eqtls_SNP_17<- EUR_eqtl_ordered_17[,c(2,4), drop =FALSE]
EUR_eqtls_SNP_18<- EUR_eqtl_ordered_18[,c(2,4), drop =FALSE]
EUR_eqtls_SNP_19<- EUR_eqtl_ordered_19[,c(2,4), drop =FALSE]
EUR_eqtls_SNP_20<- EUR_eqtl_ordered_20[,c(2,4), drop =FALSE]
EUR_eqtls_SNP_21<- EUR_eqtl_ordered_21[,c(2,4), drop =FALSE]
EUR_eqtls_SNP_22<- EUR_eqtl_ordered_22[,c(2,4), drop =FALSE]
  

#
AA_eqtls_SNP_1<- AA_eqtl_ordered_1[,c(2,4), drop =FALSE]
AA_eqtls_SNP_2<- AA_eqtl_ordered_2[,c(2,4), drop =FALSE]
AA_eqtls_SNP_3<- AA_eqtl_ordered_3[,c(2,4), drop =FALSE]
AA_eqtls_SNP_4<- AA_eqtl_ordered_4[,c(2,4), drop =FALSE]
AA_eqtls_SNP_5<- AA_eqtl_ordered_5[,c(2,4), drop =FALSE]
AA_eqtls_SNP_6<- AA_eqtl_ordered_6[,c(2,4), drop =FALSE]
AA_eqtls_SNP_7<- AA_eqtl_ordered_7[,c(2,4), drop =FALSE]
AA_eqtls_SNP_8<- AA_eqtl_ordered_8[,c(2,4), drop =FALSE]
AA_eqtls_SNP_9<- AA_eqtl_ordered_9[,c(2,4), drop =FALSE]
AA_eqtls_SNP_10<- AA_eqtl_ordered_10[,c(2,4), drop =FALSE]
AA_eqtls_SNP_11<- AA_eqtl_ordered_11[,c(2,4), drop =FALSE]
AA_eqtls_SNP_12<- AA_eqtl_ordered_12[,c(2,4), drop =FALSE]
AA_eqtls_SNP_13<- AA_eqtl_ordered_13[,c(2,4), drop =FALSE]
AA_eqtls_SNP_14<- AA_eqtl_ordered_14[,c(2,4), drop =FALSE]
AA_eqtls_SNP_15<- AA_eqtl_ordered_15[,c(2,4), drop =FALSE]
AA_eqtls_SNP_16<- AA_eqtl_ordered_16[,c(2,4), drop =FALSE]
AA_eqtls_SNP_17<- AA_eqtl_ordered_17[,c(2,4), drop =FALSE]
AA_eqtls_SNP_18<- AA_eqtl_ordered_18[,c(2,4), drop =FALSE]
AA_eqtls_SNP_19<- AA_eqtl_ordered_19[,c(2,4), drop =FALSE]
AA_eqtls_SNP_20<- AA_eqtl_ordered_20[,c(2,4), drop =FALSE]
AA_eqtls_SNP_21<- AA_eqtl_ordered_21[,c(2,4), drop =FALSE]
AA_eqtls_SNP_22<- AA_eqtl_ordered_22[,c(2,4), drop =FALSE]


AA_non <- c(AA_1,AA_2,AA_3,AA_4,AA_5,AA_6,AA_7,AA_8,AA_9,AA_10,AA_11,AA_12,AA_13,AA_14,AA_15,AA_16,AA_17,AA_18,AA_19,AA_20,AA_21,AA_22)
AA_specific_non <- c(AA_specific_1,AA_specific_2,AA_specific_3,AA_specific_4,AA_specific_5,AA_specific_6,AA_specific_7,AA_specific_8,AA_specific_9,AA_specific_10,AA_specific_11,AA_specific_12,AA_specific_13,AA_specific_14,AA_specific_15,AA_specific_16,AA_specific_17,AA_specific_18,AA_specific_19,AA_specific_20,AA_specific_21,AA_specific_22)

EUR_non <- c(EUR_1,EUR_2,EUR_3,EUR_4,EUR_5,EUR_6,EUR_7,EUR_8,EUR_9,EUR_10,EUR_11,EUR_12,EUR_13,EUR_14,EUR_15,EUR_16,EUR_17,EUR_18,EUR_19,EUR_20,EUR_21,EUR_22)
EUR_specific_non <- c(EUR_specific_1,EUR_specific_2,EUR_specific_3,EUR_specific_4,EUR_specific_5,EUR_specific_6,EUR_specific_7,EUR_specific_8,EUR_specific_9,EUR_specific_10,EUR_specific_11,EUR_specific_12,EUR_specific_13,EUR_specific_14,EUR_specific_15,EUR_specific_16,EUR_specific_17,EUR_specific_18,EUR_specific_19,EUR_specific_20,EUR_specific_21,EUR_specific_22)

chr <- c(1:22)

eqtls_data <- matrix(c(chr, EUR_non, EUR_specific_non, AA_non, AA_specific_non), ncol=5, byrow=FALSE)
eqtls_data_chi <- matrix(c( EUR_non, EUR_specific_non, AA_specific_non), ncol=3, byrow=FALSE)
colnames(eqtls_data_chi) <- c('Shared', 'European Population-specific', 'African-American Population-specific')

colnames(eqtls_data) <- c('Chromosome', 'European Non-specific', 'European Population-specific', 'African-American Non-specific', 'African-American Population-specific')

barplot(as.matrix(eqtls_data))

chisq.tests <- chisq.test(eqtls_data_chi)
chisq.tests
chisq.tests$observed
round(chisq.tests$expected,3)

eqtls_data <- matrix(c(chr, EUR_non, EUR_specific_non, AA_non, AA_specific_non), ncol=5, byrow=FALSE)
colnames(eqtls_data) <- c('Chromosome', 'Shared eQTLs', 'European Population-specific', 'African-American Non-specific', 'African-American Population-specific')

eqtls_data <- as.data.frame(eqtls_data)
cols <- c('darkslateblue','darkseagreen','coral');
ylim <- c(0,20000);
par(lwd=4);
barplot(main='Significant eQTLs in Europeans and African-Americans', sub='p-value = 10-7', t(eqtls_data[c('Shared eQTLs', 'European Population-specific', 'African-American Population-specific')]),  beside=T, ylim=ylim, border=cols, col = cols, names.arg = eqtls_data$Chromosome, xlab='Chromosome', ylab ='# of significant eQTLs', legend.text=c('Shared eQTLs', 'European Population-specific', 'African-American Population-specific'),args.legend=list(text.col=cols,col=cols,border=cols)); 
 
eqtls_data <- as.data.frame(eqtls_data)
cols <- c('darkseagreen4');
ylim <- c(0,50);
par(lwd=4);
barplot(main='Significant eQTLs in African-Americans', sub='p-value = 10-7', t(eqtls_data[c('African-American Population-specific')]),  beside=T, ylim=ylim, border=cols, col = cols, names.arg = eqtls_data$Chromosome, xlab='Chromosome', ylab ='# of significant eQTLs', legend.text=c( 'African-American Population-specific'),args.legend=list(text.col=cols,col=cols,border=cols)); 


chr_size <- read.table('chr_size.txt')
colnames(chr_size) <- c("Chromosome", "Size (bp)")
size <- as.matrix(chr_size)
size <- chr_size$`Size (bp)`
barplot(size, main='Chromosome Size Distribution (bp)', xlab = 'Chromosome Number', ylim = c(0,250000000), axes=TRUE , col = 'darkslategray4')

