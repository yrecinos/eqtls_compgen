##Python Script to correctly name files 


for i in range(1,23):
    print ("EUR_eqtls_SNP_",i, "<- EUR_eqtl_ordered_",i, "[,c(2,4), drop =FALSE]", sep='')

#Europeans eQTLs chromosome 1  ##threshold for 0.00001 through fastqtl 
#read table
for i in range(1,23):
    print('EUR_eqtl_',i, '<- read.table("EUR_eQTL_chr"',i,'results.txt"', sep='')
    print('colnames(x = EUR_eqtl_',i,") <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')", sep='')
    print('EUR_eqtl_ordered_',i, '<- subset(x = EUR_eqtl_',i,' subset = pvalue <= 0.000000i)', sep='')
    print('EUR_eqtl_ordered_',i,' <- EUR_eqtl_ordered_',i,'[order(EUR_eqtl_ordered_',i,'$pvalue),]',sep='')

#rename columns 
colnames(x = AA_eqtl_",i,") <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')
#filter for p-value
AA_eqtl_ordered_",i," <- subset(x = AA_eqtl_",i,", subset = pvalue <= 0.000000",i,")
#ordered by pvalue 
AA_eqtl_ordered_",i," <- AA_eqtl_ordered_",i,"[order(AA_eqtl_ordered_",i,"$pvalue),]
##filter african american for population specific ones 
#subset eqtls not in europeans
AA_eqtl_ordered_specific_",i," <- AA_eqtl_ordered_",i,"[!(AA_eqtl_ordered_",i,"$snpid %in% EUR_eqtl_ordered_",i,"$snpid),]
#get names of eqtls ordered pvalue 
AA_eqtl_ordered_specific_names_",i," <- AA_eqtl_ordered_specific_",i,"[order(AA_eqtl_ordered_specific_",i,"$pvalue),] 
#subset eqtls not in african americans and only in europeans
EUR_eqtl_ordered_specific_",i," <- EUR_eqtl_ordered_",i,"[!(EUR_eqtl_ordered_",i,"$snpid %in% AA_eqtl_ordered_",i,"$snpid),]

#Data for eQTLs number less than 0.00005
EUR_",i," <- nrow(EUR_eqtl_ordered_",i,")
AA_",i," <- nrow(AA_eqtl_ordered_",i,")
EUR_specific_",i," <- nrow(EUR_eqtl_ordered_specific_",i,")
AA_specific_",i," <- nrow(AA_eqtl_ordered_specific_names_",i,") 
shared_",i," <- EUR_",i," - EUR_specific_",i,"

  


AA_eqtl_ordered_only_0.0000",i," <- subset(x = AA_eqtl_ordered_only, subset = pvalue <= 0.0000",i,")
write.csv(AA_eqtl_ordered_specific_names_",i,", "AA_eqtl_ordered_specific_names_",i,".csv")     

#manhattan plot for chr ",i,"
chr",i,"_MP_AA = plot(-log",i,"0(x = AA_eqtl$pvalue))


for i in range(1,2):
    print('EUR_eqtl_',i, ' <- read.table("EUR_eQTL_chr',i,'results.txt")', sep='')
    print('colnames(x = EUR_eqtl_',i,") <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')", sep='')
    print('EUR_eqtl_ordered_',i, ' <- subset(x = EUR_eqtl_',i,' subset = pvalue <= 0.0000001)', sep='')
    print('EUR_eqtl_ordered_',i,' <- EUR_eqtl_ordered_',i,'[order(EUR_eqtl_ordered_',i,'$pvalue),]',sep='')
    print("EUR_eqtls_SNP_",i," <- EUR_eqtl_ordered_",i,"[,c(2,4), drop =FALSE]",sep='')
    print("EUR_eqtls_SNP_",i," <- EUR_eqtl_ordered_",i,"[,c(2,4), drop =FALSE]",sep='')
    print("chr",i,"_MP_EUR = plot(-log",i,"0(x = EUR_eqtl_",i,"$pvalue) , main='chr",i,"'",sep='')
    print("AA_eqtl_",i,' <- read.table("AA_eQTL_chr',i,'.results.txt")',sep='')
    print('AA_eqtl_',i, ' <- read.table("AA_eQTL_chr',i,'results.txt"', sep='')
    print('colnames(x = AA_eqtl_',i,") <- c('geneid', 'snpid', 'distance', 'pvalue', 'slope')", sep='')
    print('AA_eqtl_ordered_',i, ' <- subset(x = AA_eqtl_',i,' subset = pvalue <= 0.0000001)', sep='')
    print('AA_eqtl_ordered_',i,' <- AA_eqtl_ordered_',i,'[order(AA_eqtl_ordered_',i,'$pvalue),]',sep='')
    print("AA_eqtls_SNP_",i," <- AA_eqtl_ordered_",i,"[,c(2,4), drop =FALSE]",sep='')
    print("AA_eqtls_SNP_",i," <- AA_eqtl_ordered_",i,"[,c(2,4), drop =FALSE]",sep='')
    print("chr",i,"_MP_AA = plot(-log10(x = AA_eqtl_",i,"$pvalue) , main='chr",i,"')",sep='')

for i in range(0,23):
    EUR_22 <- nrow(EUR_eqtl_ordered_22)
    AA_22 <- nrow(AA_eqtl_ordered_22)
EUR_specific_22 <- nrow(EUR_eqtl_ordered_specific_22)
AA_specific_22 <- nrow(AA_eqtl_ordered_specific_names_22) 
shared_22 <- EUR_22 - EUR_specific_22
    

