"""
RFprofile.py
#Extract most abundant tRNA for each anticodon from 20220901_tRFs_tRNAs.txt
#Extract tRNA sequences from hg38-mature-tRNAs.fa
#We do need to include ones with introns, which should be less than 100nt. 
#After plotting individual input sam files, compare between conditions.
#I expect some tRNAs to have strong differences in the tRF distribution.

1. get coordinates for a tRNA gene, get all reads overlapping it.
2. get tRNA sequence from a file. 
3. take reads, count each nucleotide, starts, ends, and compare to reference.
4. plot base composition in a bar graph
5. plot starts and ends in rectangles.


##### plot individual tRNA genes one by one. 

example command:
cd /Users/lu/Documents/lulab/projects/snorna/tRFs
python ~/Documents/scripts/snorna/RFprofile.py \
HEK293_WT_ROS_tRFs_rep1Aligned_chr15_65869062_65869133_Gln_CTG.sam \
~/Documents/lulab/projects/hg38/hg38-tRNAs/hg38-mature-tRNAs.fa 0.1

zsh; for i j k in $(cat 20220901_tRFs_highest_tRNAs.txt | tr '\t' '\n'); \
do (l=${i/tRNA}; m=${j/:/_}; n=${m/-/_}; \
samtools view HEK293_D97KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam $j > \
HEK293_D97KO_ROS_tRFs_rep1Aligned_$n$l.sam); done

for sam in *D97KO*rep1*sam; do (python ~/Documents/scripts/snorna/RFprofile.py $sam \
~/Documents/lulab/projects/hg38/hg38-tRNAs/hg38-mature-tRNAs.fa 0.1); done

"""


import sys
import matplotlib.pyplot as plt
plt.rcParams['pdf.use14corefonts'] = True
if len(sys.argv) < 4:
    print("Usage: python RFprofile.py insam tRNA_fasta diffcut")
    sys.exit()
insam = open(sys.argv[1], 'r')
tRNAfile = open(sys.argv[2], 'r')
diffcut = float(sys.argv[3]) #diff cutoff between reference and fragments, 0.1
#file name: HEK293_WT_ROS_tRFs_rep1Aligned_chr15_65869062_65869133_Gln_CTG.sam
info = sys.argv[1][:-4].split("_")[-5:]
coord = info[0]+":"+info[1]+"-"+info[2]
coords = [int(info[1]),int(info[2])]
RNAlen = coords[1]-coords[0]+1
print(coord, coords, RNAlen)

################################################################################
##### 1. extract coverage and base composition from sam
nreads = {} #store start end coverage as a dictionary, (start,end):n
compos = {} #store base composition for each position, pos:[NA,NC,NG,NT]
for i in range(-20, RNAlen+20):
    for j in range(-20, RNAlen+20): nreads[i,j] = 0
    compos[i] = [0,0,0,0]
FLAG = 0
for line in insam:
    if line[0] == "@": continue
    FLAG,RNAME,POS,_,CIGAR = line.split()[1:6]
    POS = int(POS)
    cigars = [1 if c in "IDNSHPX=" else 0 for c in CIGAR]
    if sum(cigars): continue
    if POS <= coords[0]-20 or POS+int(CIGAR[:-1]) >= coords[1]+20: continue
    nreads[(POS-coords[0],POS-coords[0]+int(CIGAR[:-1]))]+=1
    SEQ = line.split()[9]
    for i in range(len(SEQ)):
        if POS+i<coords[0]-20 or POS+i>coords[1]+20: continue
        if SEQ[i] == "A":   compos[POS+i-coords[0]][0]+=1
        elif SEQ[i] == "C": compos[POS+i-coords[0]][1]+=1
        elif SEQ[i] == "G": compos[POS+i-coords[0]][2]+=1
        elif SEQ[i] == "T": compos[POS+i-coords[0]][3]+=1
STRAND = '-' if '{0:012b}'.format(int(FLAG))[-5]=='1' else '+'
insam.close()
#for key in sorted(nreads.keys()):
#    if nreads[key]: print(key, nreads[key])
################################################################################




################################################################################
##### 2. extract tRNA references.
tRNAfa = {} #coord:[name,strand,seq]. coord:"chr1:1-100". name:"tRNA_Ala_AGC".
name,strand,seq = '','',''
coord1 = 0
for line in tRNAfile:
    record = line.split()
    if line[0] == ">":
        if seq: tRNAfa[coord1] = [name,strand,seq]
        name,coord1,strand=record[0],record[12],record[13][1:2]
        name = "_".join(name.split("_")[-1].split("-")[:3])
        seq = ''
    else: seq += line.split()[0]
tRNAfa[coord1] = [name,strand,seq]
#print("Number of tRNA genes:", len(tRNAfa))
#print(tRNAfa["chr5:181173650-181173722"])
#print(tRNAfa["chr6:27777885-27777956"])
tRNAfile.close()
################################################################################






################################################################################
##### 3. calculate sequence variations.
name,strand,seq = tRNAfa[coord]
stacked = {} #only ones that differ from reference by diffcut, e.g., 0.1
freqs = {}
for i in range(RNAlen):
    freqs[i] = [j/sum(compos[i]) for j in compos[i]]
for i in range(RNAlen): 
    diff = 0; pos = i; nt = [0,1,2,3]
    if strand == "-": pos = RNAlen-i-1; nt = [3,2,1,0]
    #print(len(seq))
    if seq[i] == "A": diff = 1-freqs[pos][nt[0]]
    elif seq[i] == "C": diff = 1-freqs[pos][nt[1]]
    elif seq[i] == "G": diff = 1-freqs[pos][nt[2]]
    elif seq[i] in "UT": diff = 1-freqs[pos][nt[3]]
    if diff >= diffcut: stacked[i] = freqs[pos]
if strand == "-":
    for key in stacked: stacked[key] = stacked[key][::-1]
#for key in stacked: print(key,stacked[key])
################################################################################






################################################################################
##### 4. Fragment profile plots 
fig, axs = plt.subplots(3,sharex=True,figsize=(6,3))

##### A. total coverage in bars
xcoords = sorted(compos.keys())
bars = []
for key in xcoords: bars.append(sum(compos[key])) 
if STRAND == "-": xcoords = xcoords[::-1]
axs[0].bar(xcoords,bars,color="b",width=1)
axs[0].set_ylim(0,max(bars)*1.2)

##### B. coverage of individual tRFs in rectangle format.
lines = []
for key in nreads.keys():
    if not nreads[key]: continue
    start,end = key; n = nreads[key]
    if STRAND == "-":
        start,end = RNAlen-end, RNAlen-start
    lines.append((start,start)); lines.append((0,n)); lines.append('b')
    lines.append((end,end)); lines.append((0,n)); lines.append('b')
    lines.append((start,end)); lines.append((n,n)); lines.append('b')
axs[1].plot(*lines, linestyle="-")
axs[1].set_ylim(0,max(nreads.values())*1.2)

##### C. sequence variation in stacked bars.
varpos = sorted(stacked.keys())
NA = [stacked[i][0] for i in varpos]
NC = [stacked[i][1] for i in varpos]
NG = [stacked[i][2] for i in varpos]
NT = [stacked[i][3] for i in varpos]
bG = [NA[i]+NC[i] for i in range(len(NA))]
bT = [NA[i]+NC[i]+NG[i] for i in range(len(NA))]
axs[2].bar(varpos, NA, width=1, label='A', color="#00FF00")
axs[2].bar(varpos, NC, width=1, label='C', color="#0000FF", bottom=NA)
axs[2].bar(varpos, NG, width=1, label='G', color="#000000", bottom=bG)
axs[2].bar(varpos, NT, width=1, label='T', color="#FF0000", bottom=bT)
axs[2].set_ylim(0,1.2)
plt.xlim(-20,110)
plt.savefig(sys.argv[1][:-4]+".pdf")
print(sys.argv[1],coord,name,strand,seq)
################################################################################

























