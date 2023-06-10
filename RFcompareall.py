"""
RFcompareall.py.
1. scatter plots for RF coverage
2. start and end distribution plots for RFs. 

Ryvkin 2013. HAMR: high-throughput annotation of modified ribonucleotides

It appears that different fragments have different modifications, e.g., Met_CAT.
I should extract them to quantify the modification levels in RFprofile.py. 
Phe_GAA may have diff on the RT stops.
Val_CAC, seems that the m1A nonmodified ones are cut at prior to this location. 
Gln_CTG pos 59 also seems to be modified differently based on fragment length. 
Need to take care of terminal mismatching problems.

input sample groups:
WT vs. KO, in biological replicates.
siCtrl vs. siFBL and siDKC1 in biological replicates.

1. start by analyzing top ranking tRNAs, then extend to all ~500 tRNA genes. 
2. plot results in a scatter and label by tRNA name, coordinates and (start,end)
3. start analysis with Met_CAT, which showed obvious diff for the 3' half. 
4. normalize within the tRNA gene or against RPM. 

/Users/lu/Documents/lulab/projects/hg38/hg38-tRNAs/\
hg38-tRNAs-confidence-set.bed \
"""

import sys
import pysam
import numpy as np
import matplotlib.pyplot as plt
sys.path.append("/Users/lu/Documents/scripts/snorna/tRF/")
from RF import * #tRNAtable, RFbams, bam2RFcov, siCtrlbams, etc. 
print(bam2RFcov)
from scipy import stats
if len(sys.argv) < 3:
    print("Usage: python RFcompare.py tRNAbed group1bams group2bams")
    print("tRNAbed: format: (chr, start, end, name)")
    print("bams: comma separated list of bam files")
    sys.exit()
bed = sys.argv[1]
bams1 = sys.argv[2].split(",")
bams2 = sys.argv[3].split(",")
countcut = 100 #cutoff for gene expression, based on counts. 



################################################################################
##### 1. extract tRNA coordinates from bed. 
bed = open(sys.argv[1], 'r')
#example: chr1 16861921 16861991 Gly_CCC
anticodons = {}
for line in bed:
    record = line.split()
    if record[3] not in anticodons: anticodons[record[3]] = []
    anticodons[record[3]].append((record[0],int(record[1]),int(record[2])))
bed.close()
coordslist = sum(list(anticodons.values()),[])
print("Number of genes:", len(coordslist))
#print(coordslist[0])
################################################################################





################################################################################
##### 2. perform pairwise comparisons of all RFs using (start,end):[counts]
RFcovbams = []
RFnamesall = []
for bam in bams1+bams2:
    RFcovgenes = [] #one list for all fragments in a bam file. 
    for anticodon in anticodons:
        for coord in anticodons[anticodon]:
            RFcovdict =bam2RFcov(bam,anticodon,coord)#(RNA,coord,(start,end)):n
            covsum = sum(RFcovdict.values())
            RFnames = sorted(list(RFcovdict.keys()))
            RFcovlist = []
            if covsum>=countcut:RFcovlist=[RFcovdict[k]/covsum for k in RFnames]
            else: RFcovlist=[0 for k in RFnames]
            RFcovgenes.append(RFcovlist)
            RFnamesall.append(RFnames)
    RFcovbams.append(RFcovgenes)
RFcovbams = [sum(i,[]) for i in RFcovbams]
RFnamesall = sum(RFnamesall,[])
#len(RFcovbams) is number of bam files in bams1+bams2.
#each item in the list is a list of counts for all combinations of (start,end).
#len(RFcovbams[0]) == len(RFnamesall)
################################################################################






################################################################################
##### 5. plot RFs from all tRNAs (highest) in a scatter

RFcovbams1, RFcovbams2 = RFcovbams[:len(bams1)], RFcovbams[len(bams1):]
RFlist1 = list(np.mean(np.array(RFcovbams1),axis=0))
RFlist2 = list(np.mean(np.array(RFcovbams2),axis=0))
RFstd1 = list(np.std(np.array(RFcovbams1),axis=0))
RFstd2 = list(np.std(np.array(RFcovbams2),axis=0))

RFnameslt0 = []
RFlist1lt0 = []; RFlist2lt0 = []; RFstd1lt0 = []; RFstd2lt0 = []
for i in range(len(RFlist1)):
    if RFlist1[i]>=0.005 and RFlist2[i]>=0.005:
       RFnameslt0.append(RFnamesall[i])
       RFlist1lt0.append(RFlist1[i]); RFlist2lt0.append(RFlist2[i])
       RFstd1lt0.append(RFstd1[i]); RFstd2lt0.append(RFstd2[i])
print("Before and after filtering out zeros:", len(RFlist1), len(RFlist1lt0))


"""
##### plot all dots with the same color between samples 
plt.scatter(RFlist1lt0,RFlist2lt0)
plt.errorbar(RFlist1lt0,RFlist2lt0,xerr=RFstd1lt0,yerr=RFstd2lt0,fmt="o")
plt.gca().set_xlim(0,max(RFlist1lt0+RFlist2lt0)*1.2)
plt.gca().set_ylim(0,max(RFlist1lt0+RFlist2lt0)*1.2)
plt.gca().set_aspect('equal')
plt.savefig("a.pdf")
"""


##### plot short and long fragments as different colored dots
fig, ax = plt.subplots()
RFnamesshort = []
RFlist1short = []; RFlist2short = []; RFstd1short = []; RFstd2short = []
RFnameslong = []
RFlist1long = []; RFlist2long = []; RFstd1long = []; RFstd2long = []
for i in range(len(RFlist1lt0)):
    k = RFnameslt0[i]
    if k[2][1]-k[2][0]<=0.6*(k[1][2]-k[1][1]):
       RFnamesshort.append(k)
       RFlist1short.append(RFlist1lt0[i]); RFlist2short.append(RFlist2lt0[i])
       RFstd1short.append(RFstd1lt0[i]); RFstd2short.append(RFstd2lt0[i])
    else:
       RFnameslong.append(k)
       RFlist1long.append(RFlist1lt0[i]); RFlist2long.append(RFlist2lt0[i])
       RFstd1long.append(RFstd1lt0[i]); RFstd2long.append(RFstd2lt0[i])

plt.rcParams['pdf.use14corefonts'] = True
#plt.scatter(RFlist1long,RFlist2long,color="blue",alpha=0.5)
#plt.scatter(RFlist1short,RFlist2short,color="red",alpha=0.5)
plt.errorbar(RFlist1short,RFlist2short,xerr=RFstd1short,yerr=RFstd2short,
             fmt=".",color="red",alpha=0.5)
plt.errorbar(RFlist1long,RFlist2long,xerr=RFstd1long,yerr=RFstd2long,
             fmt=".",color="blue",alpha=0.5)
regshort = stats.linregress(RFlist1short, RFlist2short)
reglong = stats.linregress(RFlist1long, RFlist2long)
print(f"R-squared for short fragments: {regshort.rvalue**2:.6f}")
print(regshort.intercept,regshort.slope)
print(f"R-squared for long fragments: {reglong.rvalue**2:.6f}")
print(reglong.intercept,reglong.slope)

plt.plot(RFlist1short,regshort.intercept+regshort.slope*np.array(RFlist1short),
         'red')
plt.plot(RFlist1long,reglong.intercept+reglong.slope*np.array(RFlist1long),
         'blue')


[plt.text(x=RFlist1short[i],
          y=RFlist2short[i],
          s='  '+RFnamesshort[i][0]+' '+str(RFnamesshort[i][2]),
          size=5, rotation="vertical")
          for i in range(len(RFnamesshort)) if RFlist2short[i]>=0.15 or
          RFlist2short[i]-RFlist1short[i]>0.05]
ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)
plt.gca().set_xlim(2**-8, 2**0) #0,max(RFlist1lt0+RFlist2lt0)*1.2
plt.gca().set_ylim(2**-8, 2**0)
plt.gca().set_aspect('equal')
plt.savefig("siDKC1.pdf")

################################################################################





"""
example command for scatterplots:
cd /Users/lu/Documents/lulab/projects/snorna/tRFs
python ~/Documents/scripts/snorna/RFcompareall.py \
20220901_tRFs_highest_tRNAs.bed \
HEK293_WT_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_WT_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam \
HEK293_D101KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D101KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam
Before and after filtering out zeros: 667721 1010
R-squared for short fragments: 0.863564
0.007044212049574693 0.8409709768765936
R-squared for long fragments: 0.532527
0.007776055374384857 0.5808066716612933


D103
python ~/Documents/scripts/snorna/RFcompareall.py \
20220901_tRFs_highest_tRNAs.bed \
HEK293_WT_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_WT_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam \
HEK293_D103KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D103KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam
Before and after filtering out zeros: 667721 723
R-squared for short fragments: 0.691059
0.019660901831301262 0.9129857230093641
R-squared for long fragments: 0.338937
0.007165767629881668 0.5423375762249477


D133
python ~/Documents/scripts/snorna/tRF/RFcompareall.py \
20220901_tRFs_highest_tRNAs.bed \
HEK293_WT_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_WT_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam \
HEK293_D133KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D133KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam
Before and after filtering out zeros: 667721 723
R-squared for short fragments: 0.691059
0.019660901831301262 0.9129857230093641
R-squared for long fragments: 0.338937
0.007165767629881668 0.5423375762249477


D32A
python ~/Documents/scripts/snorna/RFcompareall.py \
20220901_tRFs_highest_tRNAs.bed \
HEK293_WT_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_WT_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam \
HEK293_D32AKO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D32AKO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam
Before and after filtering out zeros: 667721 1102
R-squared for short fragments: 0.922772
0.003333654938330568 0.8034708917511524
R-squared for long fragments: 0.888906
0.0013819198054959222 0.9577227549791689


D32A_D33
python ~/Documents/scripts/snorna/RFcompareall.py \
20220901_tRFs_highest_tRNAs.bed \
HEK293_WT_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_WT_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam \
HEK293_D32A_33KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D32A_33KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam
Before and after filtering out zeros: 667721 1095
R-squared for short fragments: 0.931079
0.003698852862180791 0.8034703328417855
R-squared for long fragments: 0.891904
0.0038983843120464212 0.8042931444999942


D33
python ~/Documents/scripts/snorna/RFcompareall.py \
20220901_tRFs_highest_tRNAs.bed \
HEK293_WT_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_WT_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam \
HEK293_D33KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D33KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam
Before and after filtering out zeros: 667721 1119
R-squared for short fragments: 0.954233
0.0031393230766620477 0.8801732182726864
R-squared for long fragments: 0.956970
0.002820232728726141 0.884182413131563


D33_D35A
python ~/Documents/scripts/snorna/RFcompareall.py \
20220901_tRFs_highest_tRNAs.bed \
HEK293_WT_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_WT_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam \
HEK293_D33_D35AKO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D33_D35AKO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam
Before and after filtering out zeros: 667721 1104
R-squared for short fragments: 0.943714
0.004790243525662676 0.8174951215123251
R-squared for long fragments: 0.850416
0.005649927153178063 0.6024354346882177


D35A
python ~/Documents/scripts/snorna/RFcompareall.py \
20220901_tRFs_highest_tRNAs.bed \
HEK293_WT_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_WT_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam \
HEK293_D35AKO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D35AKO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam
Before and after filtering out zeros: 667721 1063
R-squared for short fragments: 0.897711
0.005218663110863317 0.8061092983467194
R-squared for long fragments: 0.817166
0.003362958617106185 0.754077964970291


D97
python ~/Documents/scripts/snorna/RFcompareall.py \
20220901_tRFs_highest_tRNAs.bed \
HEK293_WT_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_WT_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam \
HEK293_D97KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D97KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam
Before and after filtering out zeros: 667721 1082
R-squared for short fragments: 0.938881
0.004487697523915257 0.9676541021285786
R-squared for long fragments: 0.641829
0.005669543595903312 0.6452653446752019


D97_D133
python ~/Documents/scripts/snorna/RFcompareall.py \
20220901_tRFs_highest_tRNAs.bed \
HEK293_WT_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_WT_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam \
HEK293_D97_133KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D97_133KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam
#49 top tRNA genes, 667721 and 1097 RFs before and after filtering RFs <0.005.
Before and after filtering out zeros: 667721 1097
R-squared for short fragments: 0.943048
0.0026216392217316004 0.9701548556662352
R-squared for long fragments: 0.887392
0.0029069937535827275 0.8138779390328409


siDKC1
python ~/Documents/scripts/snorna/RFcompareall.py \
20220901_tRFs_highest_tRNAs.bed \
HEK293_siCtrl_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_siCtrl_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam \
HEK293_siDKC1_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_siDKC1_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam,\
HEK293_siDKC2_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_siDKC2_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam
#49 top tRNA genes, 667721 and 890 RFs before and after filtering RFs <0.005.


siFBL
python ~/Documents/scripts/snorna/tRF/RFcompareall.py \
20220901_tRFs_highest_tRNAs.bed \
HEK293_siCtrl_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_siCtrl_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam \
HEK293_siFBL1_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_siFBL1_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam,\
HEK293_siFBL2_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_siFBL2_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam
#49 top tRNA genes, 667721 and 1035 RFs before and after filtering RFs <0.005.
R-squared for short fragments: 0.910277
0.0009630090957723457 0.9484256829535367
R-squared for long fragments: 0.904071
0.00011304477219443998 0.6560848496303587


siFBL
python ~/Documents/scripts/snorna/RFcompareall.py Met_CAT_top.bed \
HEK293_siCtrl_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_siCtrl_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam \
HEK293_siFBL1_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_siFBL1_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam,\
HEK293_siFBL2_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_siFBL2_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam


siFBL
python ~/Documents/scripts/snorna/RFcompareall.py Arg_TCT_top.bed \
HEK293_siCtrl_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_siCtrl_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam \
HEK293_siFBL1_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_siFBL1_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam,\
HEK293_siFBL2_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_siFBL2_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam

"""



"""
################################################################################
##### 6. plot starts and ends distribution.
#first convert the RFcovbams to a new dictionary of starts and ends coverage.
RFcovbams = [] #normalized against sum of each gene
RFendsbams = [] #newcode
RFnamesall = [] 
for bam in bams1+bams2:
    RFcovgenes = [] #one list for all fragments in a bam file.
    RFendsgenes = [] #one list for all fragments in a bam file. #newcode
    for anticodon in anticodons:
        for coord in anticodons[anticodon]:
            RFcovdict =bam2RFcov(bam,anticodon,coord)#(RNA,coord,(start,end)):n
            covsum = sum(RFcovdict.values())
            RFnames = sorted(list(RFcovdict.keys()))
            RFcovlist = []
            if covsum>=countcut:RFcovlist=[RFcovdict[k]/covsum for k in RFnames]
            else: RFcovlist=[0 for k in RFnames]
            RFcovgenes.append(RFcovlist)
            RFnamesall.append(RFnames)
            RFendsgenes.append([RFcovdict[k] for k in RFnames]) #newcode
    RFcovbams.append(RFcovgenes)
    RFendsbams.append(RFendsgenes) #newcode
RFcovbams = [sum(i,[]) for i in RFcovbams]
RFnamesall = sum(RFnamesall,[])
RFendsbams = [sum(i,[]) for i in RFendsbams] #newcode
#len(RFcovbams) is number of bam files in bams1+bams2.
#each item in the list is a list of counts for all combinations of (start,end).
#len(RFcovbams[0]) == len(RFnamesall)
################################################################################
"""




"""
################################################################################
##### 6. plot starts and ends distribution.
#input files:
RFnamesall #ranked list, correlate with the data lists RFcovbams and RFendsbams
RFcovbams  #normalized to sum of a tRNA gene
RFendsbams #raw numbers

#step1. convert the (RNA,coord,(start,end)) to standard coordinates.
leader = list(range(-19,1))
region1 = list(range(1,21))+["20a","20b"]+list(range(21,38))
region2 = ["intron"]+list(range(38,46))
vloop = ["v"+str(i) for i in range(1,16)]
region3 = list(range(46,94))
poslist = [str(i) for i in leader+region1+region2+vloop+region3]
startdict = {pos:0 for pos in poslist}
enddict = {pos:0 for pos in poslist}
print(poslist)

print(len(RFnamesall), len(RFcovbams[0]), len(RFendsbams[0]))
for i in range(len(RFendsbams[0])):
    RFinfo = RFnamesall[i]
    #count = RFendsbams[2][i] #raw read counts
    count = RFcovbams[2][i]  #normalized for each tRNA. 
    RNA,coord,(start,end) = RFinfo
    coord = coord[0]+":"+str(coord[1])+"-"+str(coord[2])
    if start==-20 or end==-20: continue #ignore this location due to 0/1-base
    start = tRNAtable[("tRNA_"+RNA,coord,start)]
    end = tRNAtable[("tRNA_"+RNA,coord,end)]
    startdict[start] += count; enddict[end] += count
for k in poslist: print("\t".join([k,str(startdict[k]),str(enddict[k])]))
################################################################################
"""
















