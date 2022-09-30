Title / Post-Imputation Pipeline

Author / Jared Taylor

About / This pipeline combines the pre and post imputation files to add back in the snps lost during imputation. The input for the pipeline is the combined pre-imputation vcf and the directory containing the files collected from the imputation software. If the pre-imputation vcf is hg19, be sure to set the crossmap variable to "yes".  

Modules / plink2, bcftools 1.9, samtools 1.9, vcftools 1.16, SnpSift, and CrossMap

Notes / To run, make sure to update all hard coded paths to the actual path that you have access to. All output files will be in a new directory named "post_imp_output". 

	sample command "./post_imp.sh <pre_imputation_vcf> <post_imputation_directory>"

Help / For any questions, email jtaylor@hudsonalpha.org
