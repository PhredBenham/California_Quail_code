from numpy import array
import numpy as np
import moments
import pylab
import os
import pandas as pd

pop_id = ["MTN", "SOU"]
sampleSizes = [23, 5]

data_dict = moments.Misc.make_data_dict_vcf("MtnQuail_Autosomes_11sept.recode.vcf", "Moqu_2pop_samples.txt")
fs = moments.Spectrum.from_data_dict(data_dict, pop_ids=pop_id, projections=sampleSizes, polarized=False)
fs.to_file("Oreortyx_2dsfs.fs")
#fs = moments.Spectrum.from_file("Oreortyx_2dsfs.fs")

seg = fs.S()
print("Number of segregating sites: %d" % seg)

#estimate total length of sequence for point estimates of demographic parameters.
L = int((seg/33021009)*905814793)
print(L)

moments.Misc.bootstrap(data_dict, pop_ids=pop_id, projections = sampleSizes, mask_corners=False, 
					polarized=False, bed_filename=None, num_boots=100,
					save_dir="oreortyx_2pop_bootstraps")

#code for renaming bootstrapped sfs as the suffix .fs is needed for importing into GADMA
cmd = "for f in SFS_*; do mv $f ${f}.fs; done"
os.system(cmd)

#there seems to be a bug of some sort in moments code for generating bootstraps
path="/media/data2/Quail_reanalysis/GADMA/MtnQuail/oreortyx_2pop_bootstraps/*.fs"
print(path)
Bootstraps = glob.glob(path)
#print(Bootstraps)
for fs in Bootstraps:
	#print(fs)
	file_name = fs.split("/")[-1]
	print(file_name)
	boot = moments.Spectrum.from_file(fs)
	folded_fs = boot.fold()
	outfile_boot="/media/data2/Quail_reanalysis/GADMA/MtnQuail/folded_quail_bootstraps/" + file_name
	folded_fs.to_file(outfile_boot)