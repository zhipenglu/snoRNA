"""
paris_snoatlas.py
compare PARIS and snoAtlas for snoRNA-rRNA and snoRNA-snRNA interactions

input files:
snoatlas_1117.txt
Target1: 54, with reported targets, 43 rRNAs and 11 snRNAs. 
Target2: 123, with reported targets, 110 rRNAs and 13 snRNAs.
tables1_snoRNA_rRNA_DGbedpe_Atlas_Plexy.txt
2663 records, including only SNORD-rRNA interactions.

##### example command: 
cd ~/Documents/scripts/snorna/
python paris_snoatlas.py ./sno/snoatlas_1117.txt \
~/Dropbox/_2023_sno/sourcedata/tables1_snoRNA_rRNA_DGbedpe_Atlas_Plexy.txt
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
if len(sys.argv) < 3:
    print("Usage: python paris_snoatlas.py snoatlas parisdata")
snoatlas = open(sys.argv[1], 'r')
parisdata = open(sys.argv[2], 'r')


################################################################################
#####1. extract snoRNA information from snoatlas_1117.txt.
#copied from snoatlas.py, modified to include up to 4 targets. 
"""
1 SNORD38A snoID_0002 RF00212 STK FA CD chr1 + 45243514 45243584 45243515
45243582 184203407 Profile Vertebrates Interaction2 NA RCH:28S-1858
target annotations: R - reported, C - conserved, H - hypothetical/predicted

0 No.
1 Name
2 ID
3 Rfam ID
4 STK
5 FA
6 Type: CD
15 Conservation
16-17 Target 1/2

"""
snoRfamdict = {} #key = Rfam if Rfam record exists
snoIDdict = {} #key = snoID for all records, including ones with Rfam records.  
for line in snoatlas:
    record = line.split("\t")
    if record[0] == "No.": continue
    Name,ID,Rfam,_,_,Type=record[1:7]; cons=record[15]; t1,t2=record[17:19]
    info = {"Name":Name,"Type":Type,"cons":cons,"t1":t1,"t2":t2}
    if Rfam != "NA": snoRfamdict[Rfam]={"ID":ID, **info}
    snoIDdict[ID]={"Rfam":Rfam, **info}
    #print(snoIDdict); break
snoatlas.close()
print("snoRNAs with Rfam and snoID records:", len(snoRfamdict),len(snoIDdict))
################################################################################



################################################################################
#####2. extract info from tables1_snoRNA_rRNA_DGbedpe_Atlas_Plexy.txt
#contains info for snoAtlas, plot the top 200 after ignoring U3/U8/U13. 
parisNA = [] #not in snoAtlas
parissnoAtlas = [] # in snoAtlas
nlines = 0
for line in parisdata:
    record = line.split()
    if record[1] == "start1": continue
    if "U3" in record or "SNORD13" in record or "SNORD118" in record: continue
    nlines+=1
    if nlines>50: break
    if record[10] == "NA": parisNA.append([nlines, int(record[7])])
    else: parissnoAtlas.append([nlines, int(record[7])])
parisdata.close()
print("SNORD-rRNA interactions in and out of snoAtlas:",
      len(parissnoAtlas), len(parisNA))
################################################################################




################################################################################
#####3. plot comparisons between snoAtlas and PARIS data.
plt.plot([i[0] for i in parisNA], [i[1] for i in parisNA], 'bo')
plt.plot([i[0] for i in parissnoAtlas], [i[1] for i in parissnoAtlas], 'ro')
plt.yscale("log"); plt.ylim(0.8,3000)
plt.savefig("a.pdf")
################################################################################













