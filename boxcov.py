"""
boxcov.py. Zhipengluchina@gmail.com. 2023-03-26. 
calculates enrichment of the arm on the snoRNA side.


################################################################################
Input sam file: 
Total_onearm2snoRNA.sam                 269625 reads
Total_onearm2snoRNA.cliques.t_o0.2.sam  43369 reads
SAM record: 23 fields. Useful: 
ID-HEK293snoRNA FLAG chr1 start 255 20M * 0 0 seq qual NH:i:1 HI:i:1 AS:i:16
nM:i:0 NM:i:0 MD:Z:20 jM:B:c,-1 jI:B:i,-1 ch:A:1
SA:Z:chr7,98881709,+,23S27M,255,1;
DG:Z:MIR3609_exon1,PUSL1_exon1,0 NG:i:0


################################################################################
snoatlas_1117.txt, 17 fields. Useful: coordinates of D and D' boxes. 
1 SNORD38A snoID_0002 RF00212 STK FA CD chr1 + 45243514 45243584 45243515
45243582 184203407 Profile Vertebrates Interaction2 NA RCH:28S-1858		


################################################################################
snoatlas_fa folder: 521 files. 
some of the files contain multiple human genes. To get all the records (n=1107):
CD: 627. HACA: 480. 
cd /Users/lu/Documents/scripts/snorna/sno/snoatlas_fa
cat *fa | grep H.sapiens > H.sapiens
Need to lift the hg19 coordinates to hg38 for comparison with SAM. Format: 
>CD_1968-1_H.sapiens_(chr12)_50850354,50850569_+_hg19_GTGATGA_8_CTGA_207
CAGCCCAGTGATGATCACTATTCCTACTTAGAGAGAATAGGACTAACTTTCAGAAATCCAGGCATTTTTCTACCTTTCATACTATCTTTCTTTCACTTTACTTCTCTTTTCTGTCTTTTATCACTTCTTTCTTTCTTCATTCTTTCTCTCTTTTGCCTGGATCGAGATTGTTAAGTCCCTCTCAGTGAAGGGTAAGATTATGAGATCTGAGGGCTG
>CD_142-1_H.sapiens_(chr11)_17097325,17097415_-_hg19_GTGTTGG_57_CAGA_40_GTGATGA_6_CTGA_83
TCACTGTGATGATGGTTTTCCAACATTCGCAGTTTCCACCAGAAAGGTTTTCCTTAGTGTTGGGTAAACCTTCCTTGGATGTCTGAGTGAG


################################################################################
Some of the annotations are incorrect:
>CD_169-2_H.sapiens_(chr17)_16342823,16342870_+_
hg19_CTGACGA_42_CAGA_30_CTGATGA_1_CTGA_42
Even for the CD box annotations, we do not have a consistent format: 
grep ">CD" H.sapiens | sed 's/[^_]//g' | awk '{ print length }'
 153 10
  62 13
 373 14
  28 17
  10 18
   1 9

Annotation is incorrect for SNORD97. 
>CD_807-1_H.sapiens_(chr11)_10823014,10823155_-_
hg19_hsa_CD_807-97_99.3-100.71-96.8_NNNNNNN_0_CTGG_46_ATGATTA_8_ATGA_133
TTGCCCGATGATTATAAAAAGACGCGTTATTAAGAGGACTTTATGCTGGAGTTCTTGACGTTTTTCTCTCTTTTCTATA
CTTCTTTTTCTTTCTTTGAATGTCCAGCGTCCTGTGAGCGAAGATTATGAGATATGAGGGCAA

################################################################################
Example command:
cd /Users/lu/Documents/lulab/projects/snorna/snoRNA_chimbam
python ~/Documents/scripts/snorna/boxcov.py \
/Users/lu/Documents/scripts/snorna/sno/H.sapiens \
Total_onearm2snoRNA.cliques.t_o0.2_SNORD97.sam

Total_onearm2snoRNA.cliques.t_o0.2.sam
SNORD97
hg38: chr11:10801467-10801608
hg19: chr11:10823014-10823155
"""


import sys
import matplotlib.pyplot as plt
from liftover import get_lifter
if len(sys.argv)<3:
    print("Usage: python boxcov.py snoRef intrxn.sam")
    exit()
w=100 #window size on each side



################################################################################
#####1. extract snoRNA annotations from H.sapiens. 
snoref = open(sys.argv[1], 'r')
CDmotifs = []; HACAmotifs = [] #information for D/D' and H/ACA motifs
CDlen1,CDlen2 = [],[]; HACAlen1,HACAlen2 = [],[] #full length sno distributions
#CDlen and HACAlen are used to show that chim reads do not cover full length sno
c = get_lifter('hg19', 'hg38')
for line in snoref:
    record = line.strip().split("_")
    chrom = record[3][1:-1]; start,end = [int(i) for i in record[4].split(",")]
    try: start,end = c[chrom][start][0][1],c[chrom][end][0][1]
    except: continue #only deal with alignments on the normal chromosomes. 
    if record[0] == ">CD":
        try: motifs = sorted([int(record[i]) for i in [-1,-3,-5,-7]])
        except: continue        
        CDmotifs.append([chrom,start,end,record[5]]+motifs)
        for i in range(-motifs[1],end-start-motifs[1]): CDlen1.append(i)
        for i in range(-motifs[3],end-start-motifs[3]): CDlen2.append(i)
    elif record[0] == ">HACA":
        motifs = sorted([int(record[i]) for i in [-1,-3]])
        HACAmotifs.append([chrom,start,end,record[5]]+motifs)
        for i in range(-motifs[0],end-start-motifs[0]): HACAlen1.append(i)
        for i in range(-motifs[1],end-start-motifs[1]): HACAlen2.append(i)
snoref.close()
CDlen1pos,CDlen1count = sorted(set(CDlen1)),[]
CDlen2pos,CDlen2count = sorted(set(CDlen2)),[]
for i in CDlen1pos: CDlen1count.append(CDlen1.count(i))
for i in CDlen2pos: CDlen2count.append(CDlen2.count(i))
HACAlen1pos,HACAlen1count = sorted(set(HACAlen1)),[]
HACAlen2pos,HACAlen2count = sorted(set(HACAlen2)),[]
for i in HACAlen1pos: HACAlen1count.append(HACAlen1.count(i))
for i in HACAlen2pos: HACAlen2count.append(HACAlen2.count(i))
print("CD and HACA snoRNAs analyzed:", len(CDmotifs), len(HACAmotifs))
#plt.plot(CDlen1pos,CDlen1count)
#plt.plot(CDlen2pos,CDlen2count)
#plt.xlim(-100,100)
#plt.savefig("CDlendist.pdf")
#plt.clf()
#print(CDmotifs[0],HACAmotifs[0]) #Examples as follows:
#['chr11', 17075778, 17075868, '-', 6, 40, 57, 83]
#['chr12', 26008145, 26008256, '+', 'AGAGAA', '36', 'ACA', '107']
################################################################################





################################################################################
#####2. convert CDmotifs to dictionaries
#key=(RNAME,POS),value=[windowstart,windowend]
#key spans the the snoRNA length, covers each nt.
#value is the window start and end, e.g., D'/D box or H/ACA pos +/- 100nt. 
CDplusdict1,CDplusdict2,CDminusdict1,CDminusdict2 = {},{},{},{}
HACAplusdict1,HACAplusdict2,HACAminusdict1,HACAminusdict2 = {},{},{},{}
for i in CDmotifs: #[chrom,start,end,strand,C,D',C',D]
    if i[3] == "+":
        for j in range(i[1],i[2]):
            CDplusdict1[(i[0],j)] = [i[1]+i[5]-w,i[1]+i[5]+w]+i[:3]
            CDplusdict2[(i[0],j)] = [i[1]+i[7]-w,i[1]+i[7]+w]+i[:3]
    elif i[3] == "-":
        for j in range(i[1],i[2]):
            CDminusdict1[(i[0],j)] = [i[2]-i[5]-w,i[2]-i[5]+w]+i[:3]
            CDminusdict2[(i[0],j)] = [i[2]-i[7]-w,i[2]-i[7]+w]+i[:3]
for i in HACAmotifs: #[chrom,start,end,strand,H,ACA]
    if i[3] == "+":
        for j in range(i[1],i[2]):
            HACAplusdict1[(i[0],j)] = [i[1]+i[4]-w,i[1]+i[4]+w]
            HACAplusdict2[(i[0],j)] = [i[1]+i[5]-w,i[1]+i[5]+w]
    elif i[3] == "-": 
        for j in range(i[1],i[2]):
            HACAminusdict1[(i[0],j)] = [i[2]-i[4]-w,i[2]-i[4]+w]
            HACAminusdict2[(i[0],j)] = [i[2]-i[5]-w,i[2]-i[5]+w]
print("CD snoRNA nt pos dictionary:", len(CDplusdict1),len(CDminusdict1)) 
################################################################################





################################################################################
#####3. extract sam records aligned to refs.
n=[0,1,2,3]
CDpluscov1,CDpluscov2,CDminuscov1,CDminuscov2 = ([0]*w*2 for i in n)
HACApluscov1,HACApluscov2,HACAminuscov1,HACAminuscov2 = ([0]*w*2 for i in n)
intrxn = open(sys.argv[2], 'r') #input sam file, after DG assembly 
for line in intrxn: #input sam file
    if line[0] == "@": continue
    record = line.split() #print(record[21])
    r = [] #range of the window, len=w*2, [start,end]
    try: RNAME,POS,LEN = record[2],int(record[3]),int(record[5][:-1])
    except: continue #ignore complex CIGAR strings for now.

    ############################################################################
    #####Note: choose types of reads. 
    if "U3" in record[21]or "U8" in record[21]or "U13" in record[21]: continue
    if "18S" not in record[21] and "28S" not in record[21] and \
       "5S" not in record[21] and "5.8S" not in record[21]: continue
    #if "tRNA" not in record[21]: continue
    ############################################################################

    #A. plot coverage relative to D' boxes
    if (RNAME,POS) in CDplusdict1 or (RNAME,POS+LEN) in CDplusdict1:
        try: r = CDplusdict1[(RNAME,POS)]
        except: r = CDplusdict1[(RNAME,POS+LEN)]
        for i in range(max(0,POS-r[0]),min(w*2,POS+LEN-r[0])): CDpluscov1[i]+=1
    if (RNAME,POS) in CDminusdict1 or (RNAME,POS+LEN) in CDminusdict1:
        try: r = CDminusdict1[(RNAME,POS)]
        except: r = CDminusdict1[(RNAME,POS+LEN)]
        for i in range(max(0,POS-r[0]),min(w*2,POS+LEN-r[0])): CDminuscov1[i]+=1
    """
    #B. plot coverage relative to D boxes
    if (RNAME,POS) in CDplusdict2 or (RNAME,POS+LEN) in CDplusdict2:
        try: r = CDplusdict2[(RNAME,POS)]
        except: r = CDplusdict2[(RNAME,POS+LEN)]
        for i in range(max(0,POS-r[0]),min(w*2,POS+LEN-r[0])): CDpluscov2[i]+=1
    if (RNAME,POS) in CDminusdict2 or (RNAME,POS+LEN) in CDminusdict2:
        try: r = CDminusdict2[(RNAME,POS)]
        except: r = CDminusdict2[(RNAME,POS+LEN)]
        for i in range(max(0,POS-r[0]),min(w*2,POS+LEN-r[0])): CDminuscov2[i]+=1
    """
    """
    #C. plot coverage relative to D' box, but separate left/right-side reads
    if (RNAME,POS) in CDplusdict1 or (RNAME,POS+LEN) in CDplusdict1:
        try: r = CDplusdict1[(RNAME,POS)] #[wstart,wend,chrom,snostart,snoend]
        except: r = CDplusdict1[(RNAME,POS+LEN)]
        for i in range(max(0,POS-r[0]),min(w*2,POS+LEN-r[0])):
            if POS+LEN/2<(r[3]+r[4])/2: CDpluscov1[i]+=1
            else: CDpluscov2[i]+=1
    if (RNAME,POS) in CDminusdict1 or (RNAME,POS+LEN) in CDminusdict1:
        print(RNAME,POS,LEN)
        try: r = CDminusdict1[(RNAME,POS)] #[wstart,wend,chrom,snostart,snoend]
        except: r = CDminusdict1[(RNAME,POS+LEN)]
        print(r)
        r = CDminusdict2[(RNAME,POS)]
        print(r)
        print(max(0,POS-r[0]),min(w*2,POS+LEN-r[0]))
        print(POS+LEN/2,(r[3]+r[4])/2)
        for i in range(max(0,POS-r[0]),min(w*2,POS+LEN-r[0])):
            if POS+LEN/2<(r[3]+r[4])/2: CDminuscov2[i]+=1
            else: CDminuscov1[i]+=1
    """
    #D. plot coverage relative to H boxes
    if (RNAME,POS) in HACAplusdict1 or (RNAME,POS+LEN) in HACAplusdict1:
        try: r = HACAplusdict1[(RNAME,POS)]
        except: r = HACAplusdict1[(RNAME,POS+LEN)]
        for i in range(max(0,POS-r[0]),min(w*2,POS+LEN-r[0])):
            HACApluscov1[i]+=1
    if (RNAME,POS) in HACAminusdict1 or (RNAME,POS+LEN) in HACAminusdict1:
        try: r = HACAminusdict1[(RNAME,POS)]
        except: r = HACAminusdict1[(RNAME,POS+LEN)]
        for i in range(max(0,POS-r[0]),min(w*2,POS+LEN-r[0])):
            HACAminuscov1[i]+=1
intrxn.close()
################################################################################






################################################################################
#####4. plot reads coverage.
D1 = [CDpluscov1[i]+CDminuscov1[::-1][i] for i in range(w*2)]
D2 = [CDpluscov2[1]+CDminuscov2[::-1][i] for i in range(w*2)]
H1 = [HACApluscov1[i]+HACAminuscov1[::-1][i] for i in range(w*2)]
H2 = [HACApluscov2[i]+HACAminuscov2[::-1][i] for i in range(w*2)]

c1,c2 = HACAlen1count,HACAlen2count
plt.plot([i+w for i in HACAlen1pos],[i/max(HACAlen1count)for i in c1])
plt.plot(list(range(w*2)),[i/max(H1) for i in H1])
plt.xlim(0,w*2); plt.savefig("HACA_box1.pdf"); plt.clf()
#plt.plot([i+w for i in HACAlen2pos],[i/max(HACAlen2count)for i in c2])
#plt.plot(list(range(w*2)),[i/max(H2) for i in H2])
#plt.xlim(0,w*2); plt.savefig("HACA_box2.pdf"); plt.clf()

"""
plt.plot([i+w for i in CDlen1pos],[i/max(CDlen1count) for i in CDlen1count])
plt.plot(list(range(w*2)),[i/max(D1) for i in D1])
plt.xlim(0,w*2); plt.savefig("CD_box1.pdf"); plt.clf()
plt.plot([i+w for i in CDlen2pos],[i/max(CDlen2count) for i in CDlen2count])
plt.plot(list(range(w*2)),[i/max(D2) for i in D2])
plt.xlim(0,w*2); plt.savefig("CD_box2.pdf"); plt.clf()
"""
################################################################################






################################################################################
#####5. plot reads on the two halves of each C/D snoRNA separately:
#r = list(range(w*2))
#plt.plot([i+w for i in CDlen1pos],[i/max(CDlen1count) for i in CDlen1count])
#plt.plot(r,[i/max(CDpluscov1+CDpluscov2) for i in CDpluscov1])
#plt.plot(r,[i/max(CDpluscov1+CDpluscov2) for i in CDpluscov2])
#plt.plot(r,[i/max(CDminuscov1+CDminuscov2) for i in CDminuscov1][::-1])
#plt.plot(r,[i/max(CDminuscov1+CDminuscov2) for i in CDminuscov2][::-1])
#plt.xlim(0,w*2)
#plt.savefig("CD_box1.pdf")
#plt.clf()
################################################################################









