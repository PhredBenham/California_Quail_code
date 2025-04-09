#!/usr/bin/env bash

#Script puts together preliminary filtering, calculation of individual het, PCA, admixture
#Specify flags to be set on command line
while getopts v:s:o:d:t: flag
do
    case "${flag}" in
        v) inputVCF=${OPTARG};; #input vcf file to be filtered, etc.
        s) sampleFile=${OPTARG};; #input sample meta-data file that is in same order as the vcf
        o) output=${OPTARG};; #directory for output files.
        d) date=${OPTARG};; #todays date in DayMonYear format, e.g. 15Feb2023
        t) taxa=${OPTARG};; #alpha code for species being analyzed, e.g. SOSP
    esac
done

source activate CCGP_popgen

#input vcf and sample file
#standard filters of vcf file
filter_out="${output}/${taxa}_filtered_${date}"
vcftools --gzvcf $inputVCF --maf 0.01 --max-alleles 2 --min-alleles 2 --minDP 5 --max-missing 1 --remove-indels --out $filter_out --recode
 
#prune filtered vcf
pruned_out="${output}/${taxa}_FilterPrune_${date}.vcf"
bcftools +prune -m 0.25 -w 2500 $inputVCF -Ov -o $pruned_out

#use vcftools to estimate per individual heterozygosity.
het_out="${output}/${taxa}_${date}_heterozygosity"
vcftools --vcf $pruned_out --het --out $het_out

pca_out="${output}/${taxa}_${date}_PCA"
#use custom script to output PCA 
python3 CCGP_PCA.py -i $pruned_out -s $sampleFile -o $pca_out -a Subspecies 

#convert vcf to admixture file .ped
ped_out="${output}/${taxa}_${date}_forAdmix"
plink --vcf $pruned_out --const-fid 0 --allow-extra-chr --recode12 --out $ped_out

admix_out="${output}/${taxa}_${date}_Admixout"
[ -d $admix_out ] || mkdir $admix_out
#run admixture

for K in 1 2 3 4 5 6 7 8 9 10;
	do 
		admixture --cv ${ped_out}.ped $K | tee $admix_out/log${K}.out; 
	done

mv *.Q $admix_out
mv *.P $admix_out

cd $admix_out
grep -h CV log*.out >Admixture_CV_out.txt
