from numpy import array
import numpy as np
import moments
import pylab
import glob
import os
import pandas as pd

pop_id = ["Gambel", "SoCal", "NorCal"]
sampleSizes = [14, 14, 14]

data_dict = moments.Misc.make_data_dict_vcf("Quail_15Nov.recode.vcf", "Quail_PopMap.txt")
fs = moments.Spectrum.from_data_dict(data_dict, pop_ids=pop_id, projections=sampleSizes, polarized=False)
fs.to_file("Callipepla_3dsfs_15Nov.fs")

seg = fs.S()
print("Number of segregating sites: %d" % seg)

#estimate length of sequence based on segregating sites in sfs.
L = int((seg/37163489)*905814793)
print(L)

moments.Misc.bootstrap(data_dict, pop_ids=pop_id, projections = sampleSizes, mask_corners=False, 
					polarized=False, bed_filename=None, num_boots=100,
					save_dir="quail_bootstraps2")

#code for renaming bootstrapped sfs as the suffix .fs is needed for importing into GADMA
cmd = "for f in SFS_*; do mv $f ${f}.fs; done"
os.system(cmd)

#there seems to be a bug of some sort in moments code for generating bootstraps
path="/media/data2/Quail_reanalysis/GADMA/Callipepla/quail_bootstraps2/*.fs"
print(path)
Bootstraps = glob.glob(path)
#print(Bootstraps)
for fs in Bootstraps:
	#print(fs)
	file_name = fs.split("/")[-1]
	print(file_name)
	boot = moments.Spectrum.from_file(fs)
	folded_fs = boot.fold()
	outfile_boot="/media/data2/Quail_reanalysis/GADMA/Callipepla/folded_quail_bootstraps/" + file_name
	folded_fs.to_file(outfile_boot)
