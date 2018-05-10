# eqtls_compgen
Final_Compgen_project 
Results are found in the Results folder. 


The purpose of this pipeline is to map eQTLs and their occurrences in African American genomes and the local ancestry around these eQTLs. 

Requirements:
Most updated version of vcftools
FastQTL
RFMix version 1
Python version 3
R Studio
Tabix

Step 1:
Take phased and imputed vcf file and pare it down using vcftools. Can use any parameters you like as long as the final size of the file is small enough to run on whatever machine you're using. As an example, I used 

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

REFER TO PART BELOW FOR fastQTL analysis. 

Step 6: 
To interpret fastQTL and RFMix results many visualizations were performed in the Matching_eQTLs_to_ancestry.ipynb. 




——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
Purpose: identify eQTLs that are specific to African American and Europeans. 

Requirements: 
Tabi
BGZip 
FastQTL 
R Studio 

1. Extracting the selected samples from the genotype and phenotype files was important and necessary in order to perform the eQTL analysis on that select population. The code on how the samples were selected and extracted was done using R. Sample input and output is found in file name: select_samples.sh 

2. File formatting for FastQTL is very important. Once the select samples were grouped, fastQTL was attempted to run. Error messages popped up. Formatting is key in FastQTL. The VCF file and the BED file (containing the expression matrices) did not contain the same name format. Troubleshooting and formatting code for this issue is found in tile_formatting.sh


3.To run FastQTL from command line example below: find the commands in file fastQTL_command.sh In the folder fastQTL_input_clean, all the files necesary for fastQTL (correctly formatted can be found). 
  ~/FastQTL/FastQTL/bin/fastQTL --vcf GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vcf.gz --bed Whole_Blood_Analysis.v6p.normalized.expression.long_names.bed.gz --include-samples Gtex_EUR_long.txt --threshold 0.00001 --region 1 --out EUR_eQTL_chr1.results.gz
  
  
4. Correctly name files in R once the fastQTL results were obtained - code in python can be found in file Python_naming_script.py. 

5. Filtering of population specific eQTLs was done on R and can be found in eQTL_R_analysis.R. 


  
  
