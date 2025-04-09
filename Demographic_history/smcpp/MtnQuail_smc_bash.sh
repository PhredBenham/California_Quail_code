#bash script for running smc++

#bgzip and tabix vcf file

#get list of scaffolds to perform analysis restrict analysis to scaffolds >1Mb
#for loop to convert each scaffold to smc file
while read scaf
do
       for i in MVZCCGP-CaOr142_II-B03 MVZCCGP-CaOr162_II-A05 MVZCCGP-CaOr132_I-H12 MVZCCGP-CaOr129_I-E12 MVZCCGP-CaOr131_I-G12 
       do
               smc++ vcf2smc -d $i $i /media/data2/Quail_reanalysis/MtnQuail_Autosomes_11sept.recode.vcf.gz Mtn_$i.${scaf}.smc.gz "$scaf" --mask ../MissingSites_mask.bed.gz MTN:MVZCCGP-CaOr146_II-E03,MVZCCGP-CaOr128_II-G05,MVZCCGP-CaOr142_II-B03,MVZCCGP-CaOr143_II-C03,MVZCCGP-CaOr145_II-D03,MVZCCGP-CaOr147_II-F03,MVZCCGP-CaOr149_II-H03,MVZCCGP-CaOr137_II-G07,MVZCCGP-CaOr166_II-D05,MVZCCGP-CaOr139_II-H07,MVZCCGP-CaOr163_II-B08,MVZCCGP-CaOr162_II-A05,MVZCCGP-CaOr126_II-F07,MVZCCGP-CaOr130_I-F12,MVZCCGP-CaOr127_II-F05,MVZCCGP-CaOr125_II-E07,MVZCCGP-CaOr138_II-H02,MVZCCGP-CaOr165_II-C05,MVZCCGP-CaOr134_II-G02,MVZCCGP-CaOr140_II-A08,MVZCCGP-CaOr141_II-A03,MVZCCGP-CaOr133_II-F02,MVZCCGP-CaOr164_II-B05,MVZCCGP-CaOr131_I-G12,MVZCCGP-CaOr148_II-G03,MVZCCGP-CaOr132_I-H12,MVZCCGP-CaOr167_II-C06,MVZCCGP-CaOr129_I-E12
       done
done <"Quail_largeautosome_list.txt"


#run estimate function on all scaffolds

smc++ estimate --timepoints 1e3 1e6 --cores 6 --spline cubic 3.14e-09 *.smc.gz

#plot SMC++ results
smc++ plot --linear MtnQuail_SMCplot_final.png -c model.final.json
smc++ plot --linear MtnQuail_SMCopt_final.png -c .model.*.json

