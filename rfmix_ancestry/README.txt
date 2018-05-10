The purpose of this pipeline is to map eQTLs and their occurrences in African American genomes and the local ancestry around these eQTLs. 

Requirements:
Most updated version of vcftools
FastQTL
RFMix version 1
Python version 3

Step 1:
Take phased and imputed vcf file file and pare it down using vcftools. Can use any parameters you like as long as the final size of the file is small enough to run on whatever machine you're using. As an example, I used 

vcftools --gzvcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_phased_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --plink --chr 1 --min-alleles 2 --max-alleles 2 --thin 40000 --out chr1_ALL

This must be done one chromosome at a time.

Step 2:
PCA.ipynb was used to generate the PCA plot with input files GTEx_Analysis_v7_eQTL_covariates/Whole_Blood.v7.covariates.txt. 

The IDs of African Americans and Europeans were written to Gtex_AA.txt and Gtex_EUR.txt for downstream eQTL analysis.

Step 3:
Now using File_reformatting.ipynb we reformat the pared down .vcf file to three files to input into RFMix. One is a class file, one is a map file of the SNPs used, the third is a genotype file in binary form. 

Step 4:
Use the files from the step generated above as input into RFMix. Example parameters:

python RunRFMix.py PopPhased ../chr2/binary_chr2_ALL.txt ../classes.txt ../chr2/chr2_snp_locations.txt -G 5 -w 0.01 -o chr2 -x 4

INSERT fastQTL stuff here

Step 6: 
To interpret fastQTL and RFMix results many visualizations were performed in the Matching_eQTLs_to_ancestry.ipynb. 



