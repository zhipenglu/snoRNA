"""
gscu.py
Gene Specific Codon Usage, originally described by Begley et al. 2007 Mol Cell.
Only analyze the nuclear proteome. The mito-proteome is too small (n=13) and
does not interact with the snoRNAs.
The following two values are calculated: gscu_gene and gscu_group
In addition to the 64 codons, iMet-AUG and SeC-UGA are calculated separately. 

1. human APPRIS non-redundant transcripts: hg38_appris_principal.txt
36144 records, 20345 genes, 18459 MANE_Select, 23281 PRINCIPAL:1, 
Homo_sapiens.GRCh38.cds.all.fa  120712 records
MANE: Matched Annotation from NCBI and EMBL-EBI (MANE)
https://appris.bioinfo.cnio.es/#/downloads

2. Whole genome codon usage properties should be pre-calculated for each gene. 
2. example one gene (e.g., HIF1A, Rapino 2018).
4. cluster and plot the GSCU in a heatmap (Begley 2007). 
5. Analyze specific codons (e.g., AUG for Met tRNA, targeted by D97/D133).

hg38_appris_cds_codonusage.txt can be used to calculate GSEA for each codon.
cut -f34 hg38_appris_cds_codonusage.txt > 
cut -f1 hg38_appris_cds_codonusage.txt | tr '_' '\t' | cut -f1 \
> hg38_appris_cds_names.txt

6. use a RNA-seq or ribo-seq dataset to calculate weighted codon usage. 
multiply gscu by raw read number should work reasonably well, however, relative
CDS and UTR length may skew this analysis. A better way is to use CDS raw reads. 

For human CDS principal isoforms: 
cd /Users/lu/Documents/scripts/snorna/codon
python gscu.py hg38_appris_principal.txt Homo_sapiens.GRCh38.cds.all.fa \
hg38_appris_cds_codonusage.txt \
/Users/lu/Documents/lulab/projects/snorna/riboseq/remap/\
20230216_RNAseq_ReadCounts.txt

For mouse CDS principal isoforms:
python gscu.py 
mm10_appris_principal_Ensembl105.txt
20230216_RNAseq_ReadCounts.txt

for file in /Users/lu/Documents/lulab/projects/snorna/riboseq/TE/\
*protein_coding.txt; do (python gscu.py hg38_appris_principal.txt \
Homo_sapiens.GRCh38.cds.all.fa \
hg38_appris_cds_codonusage.txt $file); done

"""

import sys, statistics
import matplotlib.pyplot as plt
import numpy as np

if len(sys.argv) < 4:
    print("python gscu.py appris hg38cds cutable [RNAseq|Riboseq]")
    sys.exit()
appris = open(sys.argv[1], 'r')
cds = open(sys.argv[2], 'r')
cu = open(sys.argv[3], 'w')
plt.rcParams['pdf.use14corefonts'] = True




################################################################################
#####0. set up codon information. 
codons = {
"UUU":"F","CUU":"L","AUU":"I","GUU":"V","UUC":"F","CUC":"L","AUC":"I","GUC":"V",
"UUA":"L","CUA":"L","AUA":"I","GUA":"V","UUG":"L","CUG":"L","AUG":"M","GUG":"V",
"UCU":"S","CCU":"P","ACU":"T","GCU":"A","UCC":"S","CCC":"P","ACC":"T","GCC":"A",
"UCA":"S","CCA":"P","ACA":"T","GCA":"A","UCG":"S","CCG":"P","ACG":"T","GCG":"A",
"UAU":"Y","CAU":"H","AAU":"N","GAU":"D","UAC":"Y","CAC":"H","AAC":"N","GAC":"D",
"UAA":"*","CAA":"Q","AAA":"K","GAA":"E","UAG":"*","CAG":"Q","AAG":"K","GAG":"E",
"UGU":"C","CGU":"R","AGU":"S","GGU":"G","UGC":"C","CGC":"R","AGC":"S","GGC":"G",
"UGA":"*","CGA":"R","AGA":"R","GGA":"G","UGG":"W","CGG":"R","AGG":"R","GGG":"G",
"aug":"m","uga":"U"}
#special cases: #aas += ["U","m"]; aas.sort()
#"AUG":"m", initiator Met, only the first codon
#"UGA":"U", selenocysteine, only in 26 proteins.
#AtoI editing is ignored in this analysis
#non-standard codon table is not considered for the human genome
aas = sorted(list(codons.values())) #amino acids
tri = sorted(list(codons.keys())) #trinucleotide
aa_codons = sorted([codons[codon]+"_"+codon for codon in codons]) #e.g. "F_UUU"
seleno = ["DIO1","DIO2","DIO3","GPX1","GPX2","GPX3","GPX4","GPX6","SEPHS2",
          "SPS2","TXNRD1","TXNRD2","TXNRD3","MSRB1","SELENOP","SELENOM",
          "SELENOO","SELENOW","SELENOH","SELENOT","SELENOK","SELENON",
          "SELENOV","SELENOS","SELENOI","SELENOF"] #Is it complete?
#https://www.selenodb.org/
################################################################################







################################################################################
#####1. extract a primary CDS list based on APPRIS annotation.
#genedict = {} #name:[ENST,APPRIS,MANE]
MANEdict = {} #name:ENST
MANErevdict = {} #ENST:name
for line in appris:
    record = line.split()
    if record[-1] == "MANE_Select":
        MANEdict[record[0]] = record[2]
        MANErevdict[record[2]] = record[0]
#print(MANEdict)
appris.close()
################################################################################





################################################################################
#####2. find the transcripts from Homo_sapiens.GRCh38.cds.all.fa
ENSTdict = {} #ENST: [seqs]
ENSTMANE = {} #name: [ENST, seq]
ENST = ""
for line in cds:
    if line[0] == ">": ENST = line.split()[0][1:-2]; ENSTdict[ENST] = []
    else: ENSTdict[ENST].append(line.strip('\n'))
cds.close()
for ENST in ENSTdict:
    if ENST in MANErevdict:
        ENSTMANE[MANErevdict[ENST]]=[ENST,''.join(ENSTdict[ENST])]
#print(ENSTMANE["ACTB"])
################################################################################





################################################################################
#####3. calculate codon frequency for each CDS.
ENSTcu = {} #name:[ENST,cu]; cu = {codon:freq}
for name in ENSTMANE: #ENSTs based on MANE selection.
    seq = ENSTMANE[name][1].upper()
    seq = seq.replace("T","U")
    ENSTcu[name] = [ENSTMANE[name][0],{}]
    if len(seq)%3: print("Sequence length not 3x.")
    triplets = [seq[i*3:i*3+3] for i in range(int(len(seq)/3))]
    if triplets[0] == "AUG": triplets[0] = "aug" #initiator Met
    for i in range(len(triplets)-1):
        if triplets[i] == "UGA": triplets[i] = "uga" #selenocysteine
    #print(seq, triplets); break
    for codon in codons: ENSTcu[name][1][codon]=triplets.count(codon)*3/len(seq)
#print(ENSTcu["ACTB"])
cutable = []
for i in range(66):
    list1 = []
    for name in ENSTcu: list1.append(ENSTcu[name][1][aa_codons[i][2:]])
    cutable.append(list1)
"""
#plot the cu values in a violin and box plot.
for i in cutable: print(statistics.mean(i))
plt.violinplot(cutable, showmedians=False, showextrema=False)
plt.ylim(-0.02,0.10)
plt.boxplot(cutable, widths=0.1, sym='')
plt.savefig(sys.argv[3]+"cutable.pdf")
"""
################################################################################





################################################################################
#####4. calculate standard deviations for all codons.
codoninfo = {} #codon: sd
for codon in codons: codoninfo[codon] = []
for codon in codons:
    for name in ENSTcu: codoninfo[codon].append(ENSTcu[name][1][codon])
    codoninfo[codon] = [statistics.mean(codoninfo[codon]),
                        statistics.stdev(codoninfo[codon])]
#print(codoninfo)
################################################################################






################################################################################
#####5. calculate z values for all genes:
ENSTz = {} #name:[ENST,z]; z = {codon:zvalue}
for name in ENSTcu:
    ENSTz[name] = [ENSTcu[name][0],{}]
    for codon in codons:
        ENSTz[name][1][codon] = \
        (ENSTcu[name][1][codon]-codoninfo[codon][0])/codoninfo[codon][1]
#print(ENSTz["ACTB"][1])
#plot the z values in a violin plot, and box plot.
ztable = []
for i in range(64):
    list1 = []
    for name in ENSTz: list1.append(ENSTz[name][1][aa_codons[i][2:]])
    ztable.append(list1)
"""
plt.violinplot(ztable, showmedians=False, showextrema=False)
plt.ylim(-5,5)
plt.boxplot(ztable, widths=0.1, sym='')
plt.savefig(sys.argv[3]+"ztable.pdf")
################################################################################
"""





################################################################################
#####6. export the genes and codons and z values
#labels as follows: 
#name_ENST: gene name and ENST
#aa_codon: amino acid name and codon

cu.write('\t'.join(["codons"]+aa_codons)+'\n')
for name in ENSTz:
    record = [name+"_"+ENSTz[name][0]]
    for aa_codon in aa_codons: record.append(str(ENSTz[name][1][aa_codon[2:]]))
    cu.write('\t'.join(record)+'\n')
cu.close()

"""
example selenoprotein information: 
SELENOO_ENST00000380903	4.998734908686786
TXNRD1_ENST00000525566	5.153397356543974
TXNRD3_ENST00000524230	5.201669331729214
TXNRD2_ENST00000400521	6.387014500166781
SEPHS2_ENST00000478753	7.472817384709616
SELENOV_ENST00000335426	9.677608241196376
SELENON_ENST00000361547	11.369074395934826
DIO1_ENST00000361921	13.443304561689235
GPX3_ENST00000388825	14.808217697338153
GPX6_ENST00000361902	15.142362119500538
GPX1_ENST00000419783	16.480905363574568
SELENOT_ENST00000471696	17.154729853788705
GPX2_ENST00000389614	17.604534160083485
SELENOS_ENST00000526049	17.697335890645356
SELENOM_ENST00000400299	23.039156052576555
DIO2_ENST00000438257	24.554503178804552
SELENOH_ENST00000534355	27.35248829715146
MSRB1_ENST00000361871	28.75660537119368
SELENOK_ENST00000495461	35.422466427962526
SELENOW_ENST00000601048	38.24237355889935
SELENOP_ENST00000514985	88.13384938710469
all other proteins has a z value of -0.027794646671814076, indicating 0 Sec.

iMet is absent from some proteins that start with Leu, e.g. MYC:
BAG1_ENST00000634734    -1.1720419857149387
C16orf82_ENST00000505035        -1.1720419857149387
CACNG8_ENST00000270458  -1.1720419857149387
DDX17_ENST00000403230   -1.1720419857149387
FNDC5_ENST00000373471   -1.1720419857149387
HMHB1_ENST00000289448   -1.1720419857149387
KCTD11_ENST00000333751  -1.1720419857149387
MORN2_ENST00000644631   -1.1720419857149387
MYC_ENST00000621592     -1.1720419857149387
NELFB_ENST00000343053   -1.1720419857149387
NPTXR_ENST00000333039   -1.1720419857149387
NPW_ENST00000566435     -1.1720419857149387
NR1I2_ENST00000393716   -1.1720419857149387
PANK2_ENST00000610179   -1.1720419857149387
PIGX_ENST00000392391    -1.1720419857149387
PRPS1L1_ENST00000506618 -1.1720419857149387
RNF187_ENST00000305943  -1.1720419857149387
STIM2_ENST00000467087   -1.1720419857149387
SWI5_ENST00000418976    -1.1720419857149387
TEAD1_ENST00000527636   -1.1720419857149387
TEAD4_ENST00000359864   -1.1720419857149387
TRPV6_ENST00000359396   -1.1720419857149387
TXNRD3_ENST00000524230  -1.1720419857149387
WDR26_ENST00000414423   -1.1720419857149387

TTN_ENST00000589042     -1.161776348572715
OBSCN_ENST00000680850   -1.130648213361935
NEB_ENST00000397345     -1.1287062113749307
MACF1_ENST00000564288   -1.1231429899469518
TTN is the biggest protein, and therefore, the lowest iMet. 
"""
################################################################################






################################################################################
#####7. cluster and plot codon usage z values using CLUSTER and TreeView. 
"""
cluster genes, uncentered
cluster arrays, uncentered
single linkage.
complete and average linkages are the most commonly used, but much slower than
single. Centroid is not commonly used, and it often results in artifacts.
sox2 is highly enriched in AUG codons, and may be significantly affected by
SNORD97/133. 
"""
################################################################################







################################################################################
#####8. organize RNA expression weights for codon usage analysis.
"""
for RNA-seq, use the following samples
3 WT1
5,6 D32A_1 D32A_2
8 D33_2
12 D32A_33_2
15,16 D97_1 D97_2
17,18 D133_1 D133_2
14 D97_133_2
9,10 D35A_1 D35A_2
19,20 D101_1 D101_2
21,22 D103_1 D103_2

#20 D33_35_2

cd /Users/lu/Documents/lulab/projects/snorna/riboseq/remap
awk '$2=="mRNA"' 20230216_RNAseq_ReadCounts.txt | awk \
'{print $1 "\t" $3}' > 20230216_RNAseq_ReadCounts_WT1_mRNA.txt
awk '$2=="mRNA"' 20230216_RNAseq_ReadCounts.txt | awk \
'{print $1 "\t" $5+$6}' > 20230216_RNAseq_ReadCounts_D32A_mRNA.txt
awk '$2=="mRNA"' 20230216_RNAseq_ReadCounts.txt | awk \
'{print $1 "\t" $8}' > 20230216_RNAseq_ReadCounts_D33_mRNA.txt
awk '$2=="mRNA"' 20230216_RNAseq_ReadCounts.txt | awk \
'{print $1 "\t" $12}' > 20230216_RNAseq_ReadCounts_D32AD33_mRNA.txt
awk '$2=="mRNA"' 20230216_RNAseq_ReadCounts.txt | awk \
'{print $1 "\t" $15+$16}' > 20230216_RNAseq_ReadCounts_D97_mRNA.txt
awk '$2=="mRNA"' 20230216_RNAseq_ReadCounts.txt | awk \
'{print $1 "\t" $17+$18}' > 20230216_RNAseq_ReadCounts_D133_mRNA.txt
awk '$2=="mRNA"' 20230216_RNAseq_ReadCounts.txt | awk \
'{print $1 "\t" $14}' > 20230216_RNAseq_ReadCounts_D97D133_mRNA.txt
awk '$2=="mRNA"' 20230216_RNAseq_ReadCounts.txt | awk \
'{print $1 "\t" $9+$10}' > 20230216_RNAseq_ReadCounts_D35A_mRNA.txt
awk '$2=="mRNA"' 20230216_RNAseq_ReadCounts.txt | awk \
'{print $1 "\t" $19+$20}' > 20230216_RNAseq_ReadCounts_D101_mRNA.txt
awk '$2=="mRNA"' 20230216_RNAseq_ReadCounts.txt | awk \
'{print $1 "\t" $21+$22}' > 20230216_RNAseq_ReadCounts_D103_mRNA.txt

awk '$2=="protein_coding"' RNAseq_gene_matrix_counts.txt | awk \
'{print $1 "\t" $20}' > RNAseq_gene_matrix_counts_D33D35A_protein_coding.txt

awk '$2=="protein_coding"' Riboseq_gene_matrix_counts.txt | awk \
'{print $1 "\t" $3+$4}' > Riboseq_gene_matrix_counts_WT1_protein_coding.txt
awk '$2=="protein_coding"' Riboseq_gene_matrix_counts.txt | awk \
'{print $1 "\t" $5+$6}' > Riboseq_gene_matrix_counts_D32A_protein_coding.txt
awk '$2=="protein_coding"' Riboseq_gene_matrix_counts.txt | awk \
'{print $1 "\t" $7+$8}' > Riboseq_gene_matrix_counts_D33_protein_coding.txt
awk '$2=="protein_coding"' Riboseq_gene_matrix_counts.txt | awk \
'{print $1 "\t" $9+$10}' > Riboseq_gene_matrix_counts_D32AD33_protein_coding.txt
awk '$2=="protein_coding"' Riboseq_gene_matrix_counts.txt | awk \
'{print $1 "\t" $11+$12}' > Riboseq_gene_matrix_counts_D97_protein_coding.txt
awk '$2=="protein_coding"' Riboseq_gene_matrix_counts.txt | awk \
'{print $1 "\t" $13+$14}' > Riboseq_gene_matrix_counts_D133_protein_coding.txt
awk '$2=="protein_coding"' Riboseq_gene_matrix_counts.txt | awk \
'{print $1 "\t" $15+$16}' > Riboseq_gene_matrix_counts_D97D133_protein_coding.txt
awk '$2=="protein_coding"' Riboseq_gene_matrix_counts.txt | awk \
'{print $1 "\t" $17+$18}' > Riboseq_gene_matrix_counts_D35A_protein_coding.txt
awk '$2=="protein_coding"' Riboseq_gene_matrix_counts.txt | awk \
'{print $1 "\t" $19}' > Riboseq_gene_matrix_counts_D33D35A_protein_coding.txt
awk '$2=="protein_coding"' Riboseq_gene_matrix_counts.txt | awk \
'{print $1 "\t" $21+$22}' > Riboseq_gene_matrix_counts_D101_protein_coding.txt
awk '$2=="protein_coding"' Riboseq_gene_matrix_counts.txt | awk \
'{print $1 "\t" $23+$24}' > Riboseq_gene_matrix_counts_D103_protein_coding.txt
"""
################################################################################





"""
################################################################################
#####9. calculate weighted codon usage, roughly equivalent to global codon usage
#expression file: table of gene name and value, from RNA-seq or Ribo-seq
#obviously, we need to use only protein coding genes.
#Dittmar, Goodenbour and Pan 2006, Plos Genetics
#Tissue-Specific Differences in Human Transfer RNA Expression
#export all globalcu tables, plot with parse_gscu.py and tRNA_augc_bias.py
#ENSTcu = {} #name:[ENST,cu]; cu = {codon:freq}
expdict = {}
globalcu = {}
for codon in codons: globalcu[codon] = 0
if len(sys.argv) == 5: 
    expfile = open(sys.argv[4], 'r')
    for line in expfile: record=line.split();expdict[record[0]]=float(record[1])
    total = sum(expdict.values()); print("Total reads:", total)
    expfile.close()
    for gene in expdict: expdict[gene] = expdict[gene]/total
    for gene in ENSTcu:
        for codon in codons:
            if gene not in expdict: continue
            globalcu[codon] += ENSTcu[gene][1][codon]*expdict[gene]

    globalfile = open(sys.argv[4][:-4]+"_globalcu.txt", 'w')
    globalfile.write("codons\t"+sys.argv[4].split("_")[4]+"\n")
    for aa_codon in aa_codons:
        globalfile.write(aa_codon+"\t"+str(globalcu[aa_codon[2:]])+'\n')
    globalfile.close()
################################################################################
"""







################################################################################
#####1. build dictionaries for samples. copied from RNAseq_compare.py
infile = open(sys.argv[4], 'r')
sampledict = {}; RNAdict = {}
samples = infile.readline().split()[2:]
for line in infile:
    record = line.split()
    if record[1] != "mRNA": continue  #only calculate mRNAs. 
    values = [float(i) if float(i)>=1 else 0 for i in record[2:]]
    if np.prod(values): RNAdict[record[0]] = values
infile.close()
RNAs = sorted(RNAdict.keys())
print("Number of samples:", len(samples))
print("Number of RNAs:", len(RNAs))
for sample in samples: sampledict[sample] = []
for RNA in RNAs:     
    for i in range(len(samples)): sampledict[samples[i]].append(RNAdict[RNA][i])
################################################################################







################################################################################
#####10. check codon bias in up and down regulated protein coding genes
#calculate global codon and aa usage in top/bottom 10% of diff genes in KO/WT
for i in range(1,int(len(samples)/2)):
    #merge data for WT and KO, only use data with >=10 reads in both.
    rl = range(len(RNAs))
    s1=[sampledict[samples[0]][j]+sampledict[samples[1]][j] for j in rl]
    s2=[sampledict[samples[2*i]][j]+sampledict[samples[2*i+1]][j] for j in rl]
    t1 = [s1[j] for j in rl if s1[j]>=10 and s2[j]>=10]
    t2 = [s2[j] for j in rl if s1[j]>=10 and s2[j]>=10]
    RNAsf = [RNAs[j] for j in rl if s1[j]>=10 and s2[j]>=10]
    ratiosdict = {RNAsf[j]:t2[j]/t1[j] for j in range(len(RNAsf))}
    t2dict = {RNAsf[j]:t2[j] for j in range(len(RNAsf))}
    print("Number of RNAs after filtering:", len(t1))

    #calculate total codon bias and aa bias in up and down regulated 10%
    dn10pct = {codon:0 for codon in tri}
    for r in sorted(ratiosdict.items(),key=lambda x:x[1])[:int(len(RNAsf)/10)]:
        if r[0] not in ENSTcu: continue
        for codon in tri: dn10pct[codon]+=ENSTcu[r[0]][1][codon]*t2dict[r[0]]    
    for codon in tri: dn10pct[codon] = dn10pct[codon]/sum(dn10pct.values())
    up10pct = {codon:0 for codon in tri}
    for r in sorted(ratiosdict.items(),key=lambda x:x[1])[-int(len(RNAsf)/10):]:
        if r[0] not in ENSTcu: continue
        for codon in tri: up10pct[codon]+=ENSTcu[r[0]][1][codon]*t2dict[r[0]]
    for codon in tri: up10pct[codon] = up10pct[codon]/sum(up10pct.values())
    tribias = {x:up10pct[x[2:]]/dn10pct[x[2:]] if dn10pct[x[2:]] else 1
               for x in aa_codons} #trinucleotide bias
    sorttribias = dict(sorted(tribias.items(), key=lambda x:x[1]))
    sorttri = list(sorttribias.keys())
    dn10pctaa = {aa:0 for aa in aas}; up10pctaa = {aa:0 for aa in aas}
    for codon in tri:
        dn10pctaa[codons[codon]] += dn10pct[codon]
        up10pctaa[codons[codon]] += up10pct[codon]
    aabias = {x:up10pctaa[x]/dn10pctaa[x] if dn10pctaa[x] else 1 for x in aas}
    sortaabias = dict(sorted(aabias.items(), key=lambda x:x[1]))
    sortaa = list(sortaabias.keys())

    plt.plot(sorttri, [np.log2(x) for x in list(sorttribias.values())])
    plt.ylim(-2.5,2.5)
    plt.savefig("20230216_Riboseq_topbottom10pct_codon_"+samples[2*i][:-2]+".pdf")
    plt.clf()
    plt.plot(sortaa, [np.log2(x) for x in list(sortaabias.values())])
    plt.ylim(-1.2,1.2)
    plt.savefig("20230216_Riboseq_topbottom10pct_aa_"+samples[2*i][:-2]+".pdf")
    plt.clf()
################################################################################
#sys.exit()





################################################################################
#####11. calculate codon usage in specific gene ontologies.
#e.g., mito-mRNAs, oxphos mRNAs, glycolysis
"""
~/Documents/lulab/projects/hg38/annotations/\
HALLMARK_glycolysis.txt
HALLMARK_oxidative_phosphorylation.txt
GOCC_cytoskeleton.txt
GOCC_nuclear_chromosome.txt
"""
folder = "/Users/lu/Documents/lulab/projects/hg38/annotations/"
cytosk = folder + "GOCC_cytoskeleton.txt"
chrom = folder + "GOCC_nuclear_chromosome.txt"
glyco = folder + "HALLMARK_glycolysis.txt"
oxphos = folder + "HALLMARK_oxidative_phosphorylation.txt"

GOfile = open(glyco, 'r')
GOfile.readline() #skip first 2 lines
GOfile.readline()
GOcu = {i:[] for i in tri}
g=0 #number of genes analyzed
for line in GOfile:
    if line.strip() in ENSTcu: 
        g+=1
        for i in tri: GOcu[i].append(ENSTcu[line.strip()][1][i])
GOfile.close()
for i in GOcu: print(i, '\t', sum(GOcu[i])/len(GOcu[i]))
print(g)
#plot violins for each GO category for 64 codons. 
#this difference on Met usage is not very obvious.
################################################################################















