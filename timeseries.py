"""
timeseries.py. Zhipeng Lu. 2023-02-01
plot time series for expression of gene sets.
0. get gene sets from input files
1. normalize against median
2. plot violins or other formats, or just lines if smaller gene sets

example command
cd /Users/lu/Documents/lulab/projects/snorna/mES
python ~/Documents/scripts/snorna/timeseries.py 20230318_mES_RNAseq_count.txt \
~/Documents/lulab/projects/mm10/GOBP_OXIDATIVE_PHOSPHORYLATION.txt

~/Documents/lulab/projects/mm10/GOBP_GLYCOLYTIC_PROCESS_THROUGH_GLUCOSE_6_PHOSPHATE.txt
~/Documents/lulab/projects/mm10/mt-mRNAs.txt
~/Documents/lulab/projects/mm10/GOBP_CARDIAC_MUSCLE_TISSUE_MORPHOGENESIS.txt

//www.gsea-msigdb.org/gsea/msigdb/mouse/geneset/
GOBP_FORMATION_OF_PRIMARY_GERM_LAYER

use mouse annotations
~/Documents/lulab/projects/hg38/annotations/MT-mRNAs.txt


"""
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
smpRNA,smpRibo,smpdictRNA,smpdictRibo = {},{},{},{}
RNAtype,RNAdict,Ribodict = {},{},{}




################################################################################
#####1. extract subsets of gene names, e.g., from GO collections. 
genef = open(sys.argv[2], 'r')
genedict = {line.strip() for line in genef}; genef.close()
print(genedict)
################################################################################




################################################################################
#####2. get data for RNAseq into a dictionary and merge replicates. 
dataf = open(sys.argv[1], 'r')    
smps = dataf.readline().split()[2:] #samples
for line in dataf:
    record = line.split()
    RNAdict[record[0]] = [float(i) if float(i)>=1 else 0 for i in record[2:]]
    RNAtype[record[0]] = record[1]
dataf.close()
RNAsRNA = sorted(RNAdict.keys())
for smp in smps: smpRNA[smp] = []
for RNA in RNAsRNA:
    for i in range(len(smps)): smpRNA[smps[i]].append(RNAdict[RNA][i])
print("Sample count:", len(smps))
print("RNA count:", len(RNAsRNA))
print(smps)
################################################################################




################################################################################
#####3. normalize data to median of each sample (med=1)
r=range(len(RNAsRNA))
for smp in smps:
    m=np.median(smpRNA[smp]); smpRNA[smp]=[smpRNA[smp][i]/m for i in r]
    print(smp,sum(smpRNA[smp]))
################################################################################




################################################################################
#####4. plot violins or lines for the subset:
wt,ko = ['mES_WT','EB_WT','CM_WT'],['mES_D133KO5','EB_D133KO5','CM_D133KO5']
subwt,subko = {},{}
for i in range(len(RNAsRNA)):
    if RNAsRNA[i] in genedict:
        subwt[RNAsRNA[i]]=[smpRNA[s][i] for s in wt]
        subko[RNAsRNA[i]]=[smpRNA[s][i] for s in ko]
#print(subwt["mt-Atp6"]); print(subko["mt-Atp6"])

locs = [0,1,2]
medwt = np.median([i[0] for i in list(subwt.values())])
medko = np.median([i[0] for i in list(subko.values())])
fig, ax = plt.subplots()

for i in subwt: 
    #plot lines. geomean is not really correct between two samples. 
    plt.plot(locs,subwt[i]/stats.gmean(subwt[i]),"b-")
    plt.plot(locs,subko[i]/stats.gmean(subko[i]),"r-")
    #plt.plot(locs,subwt[i]/medwt/subwt[i][0],"b-")
    #plt.plot(locs,[j*medko for j in subko[i]]/medwt/medwt/subko[i][0],"r-")
#plt.ylim(0,1)
plt.savefig("d.pdf")
plt.clf()
################################################################################


















