#!/bin/bash

starting_dir=$PWD/
post_imp_dir=$2
pre_imp_vcf=$1
snp_sift=/cluster/home/jtaylor/scripts/Pre-Imputation-QC-Pipeline/dependencies/SnpSift.jar
crossmap="yes"

module load cluster/bcftools/1.9
module load cluster/htslib/1.9
module load cluster/vcftools/0.1.16
module load cluster/python/3.7.6

mkdir processing
mkdir post_imp_output
mkdir logs

cd $post_imp_dir

for i in {1..22} X; do bsub -J "process1_${i}" -o ${starting_dir}logs/process1_${i}.out "java -Xmx64G -Djava.io.tmpdir=${post_imp_dir} -jar ${snp_sift} filter \"((AF>0) & (R2>0.3)) | ((exists TYPED) | (exists TYPED_ONLY))\" chr${i}.dose.vcf.gz > ${starting_dir}processing/chr${i}.dose_hasVar_R2-0.3.vcf && bcftools annotate -x FORMAT ${starting_dir}processing/chr${i}.dose_hasVar_R2-0.3.vcf > ${starting_dir}processing/chr${i}.dose_hasVar_R2-0.3_noFormat.vcf && bgzip -c ${starting_dir}processing/chr${i}.dose_hasVar_R2-0.3_noFormat.vcf > ${starting_dir}processing/chr${i}.dose_hasVar_R2-0.3_noFormat.vcf.gz && tabix -p vcf ${starting_dir}processing/chr${i}.dose_hasVar_R2-0.3_noFormat.vcf.gz"; done

bwait -w "ended(process1_*)"

cd ${starting_dir}processing/

for i in {1..22} X; do bsub -J "process2_${i}" -o ${starting_dir}logs/process2_${i}.out "sed -e \"s/chr//\" chr${i}.dose_hasVar_R2-0.3_noFormat.vcf > chr${i}.dose_hasVar_R2-0.3_noFormat_NOchr.vcf && bgzip -c chr${i}.dose_hasVar_R2-0.3_noFormat_NOchr.vcf > chr${i}.dose_hasVar_R2-0.3_noFormat_NOchr.vcf.gz && tabix -p vcf chr${i}.dose_hasVar_R2-0.3_noFormat_NOchr.vcf.gz"; done

bwait -w "ended(process2_*)"

for i in {1..22} X; do bsub -J "process3_${i}" -o ${starting_dir}logs/process3_${i}.out "bcftools annotate chr${i}.dose_hasVar_R2-0.3_noFormat_NOchr.vcf.gz -a /cluster/home/ncochran/dbSNP/hg38/GCF_000001405.38_ENSEMBL.${i}.vcf.gz -c ID -Ov > chr${i}.dose_hasVar_R2-0.3_noFormat_dbSNP-154.vcf"; done

bwait -w "ended(process3_*)"

java -Xmx64G -Djava.io.tmpdir=${starting_dir}/processing -jar ${snp_sift} split -j chr1.dose_hasVar_R2-0.3_noFormat_dbSNP-154.vcf chr2.dose_hasVar_R2-0.3_noFormat_dbSNP-154.vcf chr3.dose_hasVar_R2-0.3_noFormat_dbSNP-154.vcf chr4.dose_hasVar_R2-0.3_noFormat_dbSNP-154.vcf chr5.dose_hasVar_R2-0.3_noFormat_dbSNP-154.vcf chr6.dose_hasVar_R2-0.3_noFormat_dbSNP-154.vcf chr7.dose_hasVar_R2-0.3_noFormat_dbSNP-154.vcf chr8.dose_hasVar_R2-0.3_noFormat_dbSNP-154.vcf chr9.dose_hasVar_R2-0.3_noFormat_dbSNP-154.vcf chr10.dose_hasVar_R2-0.3_noFormat_dbSNP-154.vcf chr11.dose_hasVar_R2-0.3_noFormat_dbSNP-154.vcf chr12.dose_hasVar_R2-0.3_noFormat_dbSNP-154.vcf chr13.dose_hasVar_R2-0.3_noFormat_dbSNP-154.vcf chr14.dose_hasVar_R2-0.3_noFormat_dbSNP-154.vcf chr15.dose_hasVar_R2-0.3_noFormat_dbSNP-154.vcf chr16.dose_hasVar_R2-0.3_noFormat_dbSNP-154.vcf chr17.dose_hasVar_R2-0.3_noFormat_dbSNP-154.vcf chr18.dose_hasVar_R2-0.3_noFormat_dbSNP-154.vcf chr19.dose_hasVar_R2-0.3_noFormat_dbSNP-154.vcf chr20.dose_hasVar_R2-0.3_noFormat_dbSNP-154.vcf chr21.dose_hasVar_R2-0.3_noFormat_dbSNP-154.vcf chr22.dose_hasVar_R2-0.3_noFormat_dbSNP-154.vcf chrX.dose_hasVar_R2-0.3_noFormat_dbSNP-154.vcf > combined_post_imp_dbSNP-154.vcf && bgzip -@24 -c combined_post_imp_dbSNP-154.vcf > combined_post_imp_dbSNP-154.vcf.gz && tabix -p vcf combined_post_imp_dbSNP-154.vcf.gz && /cluster/home/jtaylor/software/plink2/plink2 --vcf combined_post_imp_dbSNP-154.vcf.gz --make-bed --out post_imp --set-missing-var-ids @:#:\$r:\$a --const-fid --new-id-max-allele-len 50 truncate --vcf-half-call missing

if [ "$crossmap" = "yes" ]
then
	source /cluster/home/jtaylor/software/crossmap/env/bin/activate 
	
	CrossMap.py vcf /cluster/home/jtaylor/reference_files/post_imputation_pipeline_files/hg19ToHg38.over_NOchr.chain $pre_imp_vcf /cluster/home/jtaylor/reference_files/post_imputation_pipeline_files/hg38.fa pre_imp.hg38.vcf
	
	deactivate
else
	cp $pre_imp_vcf pre_imp.hg38.vcf
fi

vcf-sort pre_imp.hg38.vcf > pre_imp_sorted.hg38.vcf

/cluster/home/jtaylor/software/plink2/plink2 --vcf pre_imp_sorted.hg38.vcf --make-bed --out pre_imp --set-all-var-ids @:#:\$r:\$a --const-fid --new-id-max-allele-len 50 truncate --vcf-half-call missing --allow-extra-chr

bgzip -c pre_imp_sorted.hg38.vcf > pre_imp_sorted.hg38.vcf.gz && tabix -p vcf pre_imp_sorted.hg38.vcf.gz

cut -f2 pre_imp.bim > TypedIDs.txt

cut -f2 post_imp.bim > ImputedIDs.txt

cat TypedIDs.txt ImputedIDs.txt | sort | uniq -u > UniqAcrossTypedAndImputed.txt

cat UniqAcrossTypedAndImputed.txt TypedIDs.txt | sort | uniq -d > TypedOnly.txt

/cluster/home/jtaylor/software/plink2/plink2 --bfile pre_imp --extract TypedOnly.txt --make-bed --out TypedOnly --allow-extra-chr

/cluster/home/jtaylor/software/plink2/plink2 --bfile TypedOnly --recode vcf id-paste=iid --allow-extra-chr --out TypedOnly --real-ref-alleles

bgzip -c TypedOnly.vcf > TypedOnly.vcf.gz && tabix -p vcf TypedOnly.vcf.gz

awk '! /\#/' TypedOnly.vcf > TypedOnly_NoHead.txt

cat combined_post_imp_dbSNP-154.vcf TypedOnly_NoHead.txt | vcf-sort -c > post_imp_wTypedOnly.vcf && bgzip -c post_imp_wTypedOnly.vcf > post_imp_wTypedOnly.vcf.gz && tabix -p vcf post_imp_wTypedOnly.vcf.gz

mv post_imp_wTypedOnly.vcf.gz ${starting_dir}post_imp_output/

cd $starting_dir

#rm -r processing/

