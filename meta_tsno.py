"""
meta_tsno.py. Zhipeng Lu, 2023-01-01. To do: consolidate tRNA genes into groups.
meta analysis of tRNA-snoRNA interactions. Input file contains these fields: 

1. snoRNA_name  e.g., SNORD97
2. tRNA_name    e.g., tRNA_Met_CAT-chr6:28944575:28944647
3. snoRNA_seq   
4. tRNA_seq
5. gudie_box    D, D'
6. guide_seq    GACGCGTT, variable length
7. tRNA_snoRNA_duplex     e.g., CAGCGCGTCA&AGACGCGTTA
8. tRNA_snoRNA_structure  e.g., .((((((((.&.)))))))).
9. Nm_pos      e.g., 22
10.mini_MFE    e.g., -13.6
11.shift_MFE   e.g., "-9.1,-7.5,-4.9,...", could be negative or positive
12.PARIS_reads e.g., 84
13.CLIP_reads  e.g., 3

#extract information from the table
snoRNA_name,tRNA_name,snoRNA_seq,tRNA_seq,gudie_box,guide_seq,\
tRNA_snoRNA_duplex,tRNA_snoRNA_structure,Nm_pos,mini_MEF,shift_MEF,\
PARIS_reads,CLIP_reads =tuple(record)

Example command:
cd ~/Documents/lulab/projects/snorna/snoRNA_tRNA_intrxn
python ~/Documents/scripts/snorna/meta_tsno.py \
20230301_snoRNA_tRNA_interaction_table_shuffled_exponly.txt

#this file contains all the shuffled duplexes
20230301_snoRNA_tRNA_interaction_table_shuffled.txt
#this file contains all the shuffled duplexes, for the experimental pairs. 
20230301_snoRNA_tRNA_interaction_table_shuffled_exponly.txt
#this file contains all interactions, predicted and experimentally determined. 
snoRNA_tRNA_interaction_table.txt
#this file contains only the experimentally determined ones, why only 490. 
snoRNA_tRNA_interaction_table_exp.txt
"""


import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

tsnofile = open(sys.argv[1], 'r')
miniall = [] #all interactions, minimal MFE
miniexp = [] #experiment derived interactions
metaall = [[] for i in range(60)] #60 lists
metaexp = [[] for i in range(60)] #60 lists

################################################################################
shuffled = True
################################################################################

#1. extract MFEs from experimental and shuffled data. 
for line in tsnofile:
    if line[:11] == "snoRNA_name": continue
    record = line.split()
    mini_MFE,shift_MFE = 0,record[6] #for shuffled data. 
    if not shuffled: mini_MFE,shift_MFE = tuple(record[9:11]) #for original data
    MFEs = shift_MFE.strip("\"").split(',')
    miniall.append(float(mini_MFE))
    for i in range(60):
        if MFEs[i] != "NA": metaall[i].append(float(MFEs[i]))
        else: metaall[i].append(0)
    if shuffled or (int(record[11]) or int(record[12])): #support by PARIS/CLIP
        miniexp.append(float(mini_MFE))
        for i in range(60):
            if MFEs[i] != "NA": metaexp[i].append(float(MFEs[i]))
            else: metaexp[i].append(0)
#for lis in metaall: print(sum(lis)/len(lis))
#for lis in metaexp: print(sum(lis)/len(lis))
print("mini MFEs, all:", sum(miniall)/len(miniall))
print("mini MFEs, exp:", sum(miniexp)/len(miniexp))
print("No. total predicted interactions", len(miniall))
print("No. experimental interactions", len(miniexp))

"""
plt.plot(sorted(miniexp))
plt.savefig("tsno_miniexp.pdf")
plt.close()
plt.plot(sorted(sum(metaexp,[])))
plt.savefig("tsno_metaexp.pdf")
plt.close()
plt.plot(sorted(miniall))
plt.savefig("tsno_miniall.pdf")
plt.close()
plt.plot(sorted(sum(metaall,[])))
plt.savefig("tsno_metaall.pdf")
plt.close()
"""

#2. plot all interactions supported by PARIS and CLIP, medians, or shuffled.
MFEexpmatrix = metaexp[:30]+[miniexp]+metaexp[30:]
plt.plot([np.mean(i) for i in MFEexpmatrix])
plt.ylim(-12,0)
plt.savefig("tsno_MFEexpmatrix_median_shuffled.pdf")
plt.close()
plt.clf()

#3. plot matrix for all interactions from PARIS and CLIP.
#ranked by interaction read number, plot heat map of all interactions.
matrix = MFEexpmatrix
fig, ax = plt.subplots(figsize=(3, 2))
psm = plt.pcolor(matrix, cmap='bwr', vmin=-20, vmax=20)
fig.colorbar(psm, ax=ax)
plt.savefig("tsno_MFEexpmatrix_shuffled.pdf")



"""
#plot median of MFE from the matrix. calculating confidence intervals?
plt.plot(list(range(0,61)))
ax = sns.lineplot(list(range(0,61)),MFEexpmatrix)
plt.savefig("tsno_MFEexpmatrix_median_shuffled.pdf")
"""












