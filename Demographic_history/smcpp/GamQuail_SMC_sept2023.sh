#bash script for running smc++

#bgzip and tabix vcf file

#get list of scaffolds to perform analysis restrict analysis to scaffolds >1Mb
#for loop to convert each scaffold to smc file
while read scaf
do
       for i in MVZCCGP-CaOr114_I-D11 MVZCCGP-CaOr113_I-C11 MVZCCGP-CaOr116_I-E11 MVZCCGP-CaOr119_I-H11 MVZCCGP-CaOr118_I-G11
       do
               smc++ vcf2smc -d $i $i /media/data2/Quail_reanalysis/GamQuail_Autosome_16Sept.recode.vcf.gz Gam_$i.${scaf}.smc.gz "$scaf" --mask ../MissingSites_mask.bed.gz GAM:MVZCCGP-CaOr115_II-D02,MVZCCGP-CaOr123_II-E02,MVZCCGP-CaOr124_I-D12,MVZCCGP-CaOr120_I-A12,MVZCCGP-CaOr110_II-C02,MVZCCGP-CaOr109_I-B11,MVZCCGP-CaOr122_I-C12,MVZCCGP-CaOr117_I-F11,MVZCCGP-CaOr121_I-B12,MVZCCGP-CaOr114_I-D11,MVZCCGP-CaOr113_I-C11,MVZCCGP-CaOr116_I-E11,MVZCCGP-CaOr119_I-H11,MVZCCGP-CaOr118_I-G11
       done
done <"Quail_largeautosome_list.txt"

#run estimate function on all scaffolds
smc++ estimate --timepoints 1e3 1e6 --cores 6 --spline cubic 3.14e-09 *.smc.gz

#plot SMC++ results
smc++ plot --linear GamQuail_SMCplot_final.png -c model.final.json
smc++ plot --linear GamQuail_SMCopt_final.png -c .model.*.json
