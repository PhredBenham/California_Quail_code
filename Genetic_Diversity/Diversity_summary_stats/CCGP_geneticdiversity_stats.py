#! /usr/bin/env python3

import numpy as np
import scipy
import pandas
from Bio import SeqIO
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import seaborn as sns
import argparse
sns.set_style('white')
sns.set_style('ticks')
sns.set_context('notebook')
import h5py
import allel; print('scikit-allel', allel.__version__)

##########################################################################################
def SumStat_Scatterplot(pop,yval,data,ylab,filename):
	fig, ax = plt.subplots(figsize=(14, 6))
	sns.despine(ax=ax, offset=5)
	ax = sns.barplot(x=pop, y=yval, data=data, errwidth=0.25, capsize=0.1)
	#x = sns.catplot(x="Population", y='F', hue='SamplePeriod', col="DataType", kind="violin",data=data)hue='DataType',
	#ax = sns.pointplot(x='Population', y='value', hue='Time', data=data,errorbar="ci", errwidth=0.25, capsize=0.1)	
	#ax.legend(bbox_to_anchor=(1, 1), loc='upper left', fontsize='medium')
	#ax.axhline(0, color="black", lw=0.5)
	ax.set_ylabel(ylab)
	fig.tight_layout()
	fig.savefig(filename, dpi=300)
##########################################################################################

Usage="""
Script for estimating genetic diversity in sliding windows using scikit-allel. 
Requires a vcf file and a sample file in .csv format that at minimum includes sample names matching samples in vcf and 
a column with population information.
"""

parser = argparse.ArgumentParser(description = Usage)
parser.add_argument("-i","--input", required=True, help="input vcf or hdf5 file")
parser.add_argument("-o","--output", required=True, help="output file name prefix")
parser.add_argument("-s","--sample_data", required=True, help="csv file with sample info, samples must be in same order as vcf file")
parser.add_argument("-p","--population", nargs='?', const=1, type=str, default=None, help="column name for estimating populations")
parser.add_argument("-w","--window_size", required=True, type=int, default=15000, help="size of windows for estimating pop gen stats")
parser.add_argument("-g", "--genome_size", required=True, type=int, help="Size of genome")
parser.add_argument("-f","--ref_genome", required=True, help="reference genome for creating mask of inaccessible sites")
args = parser.parse_args()

VCFfile=args.input
OutFile=args.output
SampleData=args.sample_data
PopCategory=args.population
WindowSize=args.window_size
GenomeSize=args.genome_size
RefGenome=args.ref_genome


#input reference genome, use bedtools or other tool before hand to mask callable sites or vice versa
fastaSeq = SeqIO.to_dict(SeqIO.parse(RefGenome,'fasta'))

#input data, check if h5 file has already been created. if not create h5 file and import.
VCF_path = VCFfile
h5file = VCF_path + '.h5'
if not os.path.exists(h5file):
	allel.vcf_to_hdf5(VCF_path, h5file, fields='*', overwrite=True)

#load hdf5 file with genotype calls and sample data
capture_callset = h5py.File(h5file, mode='r')

#load csv sample file. Order of samples in CSV must match order in the VCF file!!!
samples = pandas.read_csv(SampleData)
print(samples.head())

#insert check to ensure that vcf file and sample file match one another (kill program if false)
#or maybe if false I can reorder pandas dataframe to match order in vcf.

#create index for each scaffold and position
chrom = capture_callset['variants/CHROM']
#print(chrom[1])
pos = capture_callset['variants/POS']
idx = allel.SortedMultiIndex(chrom,pos)
vtbl = allel.VariantChunkedTable(capture_callset['variants'], names=['CHROM','POS', 'REF', 'ALT', 'QUAL'])
genotypes = allel.GenotypeChunkedArray(capture_callset['calldata/GT'])
scaffolds_uniq = pandas.unique(np.array(chrom))
#print(scaffolds_uniq)
scaf_number = list(range(1,len(scaffolds_uniq)+1))
scaf_dict = dict(zip(scaffolds_uniq,scaf_number))
print(scaf_dict)	

#get array of populations/time period
populations = list(samples[PopCategory].unique())
print(populations)

#output txt file to write diversity stat estimates. 
OutFileName = OutFile + ".txt"
DivOutput = open(OutFileName, "w")

#loop over each population and estimate diversity stats. 
for pop in populations:
	print(pop)
	columns = ['scaffold', 'NCBI_scaf', 'start_POS', "pi", "tajimaD", "wattersonTheta" ]

	windows = int(GenomeSize/5000)
	
	df = pandas.DataFrame(columns=columns, index=range(windows))
	
	print(df)
	
	num_pi = 0
	
	for chr in scaffolds_uniq:
		print(chr)
		#make boolean array of masked positions in the genome
		scaffolds = (str(fastaSeq[chr].seq))
		chr1 = np.array(list(scaffolds))
		chr2 = chr1 == 'N'
		print(len(chr2))
		
		#subset vtbl and genotypes based on chromosome
		loc_chrom = idx.locate_range(chr)
		print(loc_chrom)
		vtbl_chrom = vtbl[loc_chrom]
		vtbl_chrom_pos = vtbl_chrom['POS']
		
		#make genotype table for chromosome
		gt_chrom = genotypes[loc_chrom]
		if len(vtbl_chrom_pos) > 1:
			sample_selection = samples[PopCategory].isin({pop}).values
			#samples_subset = samples[sample_selection]
			#samples_subset.reset_index(drop=True, inplace=True)
			variant_selection = vtbl_chrom.eval('QUAL > 1')[:]
			#ac_pop = genotypes.subset(variant_selection, sample_selection)
			genotypes_subset = gt_chrom.subset(variant_selection, sample_selection)
			ac_pop = genotypes_subset.count_alleles()
	
			#calculate summary stats
			pi, windows, n_bases, counts = allel.windowed_diversity(vtbl_chrom_pos,ac_pop,size=WindowSize, is_accessible=chr2)
			theta,_,_,_ = allel.windowed_watterson_theta(vtbl_chrom_pos,ac_pop,size=WindowSize, is_accessible=chr2)
			tajD,_,_ = allel.windowed_tajima_d(vtbl_chrom_pos,ac_pop,size=WindowSize)
	
			#append to pandas dataframe
			index_range = range(num_pi,num_pi+len(pi))
			window_start = windows[:,0]
			NCBIscaf = [chr]*len(pi)
			scaf_num = [scaf_dict[chr]]*len(pi)
			df.loc[index_range,'scaffold'] = scaf_num
			df.loc[index_range,'NCBI_scaf'] = NCBIscaf
			df.loc[index_range,'start_POS'] = window_start
			df.loc[index_range,"pi"] = pi
			df.loc[index_range,"wattersonTheta"] = theta
			df.loc[index_range,"tajimaD"] = tajD
				
			num_pi += len(pi)
				
	NucDiv = df.dropna(how='all')
	print(NucDiv)
	print(num_pi)
	OutFileCSV = OutFile + pop + ".csv"
	NucDiv.to_csv(OutFileCSV, na_rep='NaN')
	
	#NucDiv1 = pandas.read_csv(OutFile, index_col=0)
	#NucDiv = NucDiv1[NucDiv1['scaffold'] < 120] 
		
	pi_mean=np.nanmean(NucDiv['pi'], dtype='float32')
	pi_max=np.nanmax(NucDiv['pi'])
	pi_min=np.nanmin(NucDiv['pi'])
	pi_std=np.nanstd(NucDiv['pi'], dtype='float32')
	
	
	theta_mean=np.nanmean(NucDiv['wattersonTheta'], dtype='float32')
	theta_max=np.nanmax(NucDiv['wattersonTheta'])
	theta_min=np.nanmin(NucDiv['wattersonTheta'])
	theta_std=np.nanstd(NucDiv['wattersonTheta'],dtype='float32')
	
	tajd_mean=np.nanmean(NucDiv['tajimaD'], dtype='float32')
	tajd_max=np.nanmax(NucDiv['tajimaD'])
	tajd_min=np.nanmin(NucDiv['tajimaD'])
	tajd_std=np.nanstd(NucDiv['tajimaD'], dtype='float32')
	
	#print mean species gendiv
	print(pop, "nuc div:", WindowSize, pi_mean,pi_std, sep='\t', file=DivOutput)
	print(pop, "theta:", WindowSize, theta_mean,theta_std, sep='\t', file=DivOutput)
	print(pop, "tajd:", WindowSize, tajd_mean,tajd_std, sep='\t', file=DivOutput)


DivOutput.close()

'''
# How to plot nuc div by window along chrom.
NucDiv['window'] = range(len(NucDiv))
df_grouped = NucDiv.groupby(('scaffold'))

# manhattan plot
fig = plt.figure(figsize=(16, 6)) # Set the figure size
ax = fig.add_subplot(111)
#colors = ['#FF7B84','#FFBBC0'] #cal_quail 
colors = ['#05BFC4','#69FBFF'] #gam_quail 
#colors = ['#189202',"#71FA5A"] #mtn_quail
x_labels = []
x_labels_pos = []
for num, (name, group) in enumerate(df_grouped):
    group.plot(kind='line', x='window', y='pi',color=colors[num % len(colors)], ax=ax, legend=None)
    x_labels.append(name)
    x_labels_pos.append((group['window'].iloc[-1] - (group['window'].iloc[-1] - group['window'].iloc[0])/2))
ax.set_xticks(x_labels_pos)
ax.set_xticklabels(x_labels)

# set axis limits
ax.set_ylim([0, 0.004])
plt.axhline(y = pi_mean, color ='#189202', linestyle = '--')
plt.axhline(y = 0.000237, color = '#05BFC4', linestyle = '--')
plt.axhline(y = 0.00032101702768886617, color = '#FF7B84', linestyle = '--')


# x axis label
ax.set_xlabel('Scaffold')

# show the graph
PlotName = Species + "_NucDivPlot_June2024.pdf"
fig.savefig(PlotName, dpi=150)


fig, ax = plt.subplots(figsize=(7, 5))
sns.despine(ax=ax, offset=5)
ax = sns.histplot(x="value", hue="variable", data=NucDiv_melt)
fig.savefig("callipepla_hist_nucdiv.pdf", dpi=300)


SumStats = pandas.read_csv(OutFile, header=[0,1], index_col=0)
print(SumStats)
Sumstats_melt = SumStats.melt()
print(Sumstats_melt)
Sumstats_melt.to_csv("Callipepla_sumstats_8Aug2023_melt.csv", na_rep='NaN')

NucDiv = Sumstats_melt[Sumstats_melt['Statistic'].isin(['pi'])]
WattTheta = Sumstats_melt[Sumstats_melt['Statistic'].isin(['wattersonTheta'])]
TajD = Sumstats_melt[Sumstats_melt['Statistic'].isin(['tajimaD'])]

data = [NucDiv, WattTheta, TajD]

ylab = ["Nucleotide diversity", "Watterson's theta", "Tajima's D"]
filename = ["Callipepla_NucDiv_8Aug2023.pdf", "Callipepla_WattTheta_8Aug2023.pdf", "Callipepla_TajD_8Aug2023.pdf"]

for i in range(3):
	SumStat_Scatterplot("Population", 'value', data[i],ylab[i],filename[i])
'''
