"""
tRNA_codonbias.py

1. plot the results from gscu.py: 22 files
Similar to Fig. 2 in Bornelov 2019 Genome Biol. But here we first quantify
global occupancy, weighted by expression, and then analyzed the KO/WT ratios.
2. plot the correlation between aa usage in mRNA levels vs. tRNA levels
3. plot variations in aa usage:
4. pseudo time series to show codon and aa usage from single and double KOs

example command:
cd /Users/lu/Documents/lulab/projects/snorna/riboseq/TE
python ~/Documents/scripts/snorna/codon/tRNA_codonbias.py 
"""

import sys
import numpy as np
import scipy
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 6})
plt.rcParams['pdf.use14corefonts'] = True
from sklearn.linear_model import LinearRegression
from decimal import Decimal
SMALL_SIZE, MEDIUM_SIZE, BIGGER_SIZE = 6, 10, 12
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=SMALL_SIZE)  # fontsize of the figure title


################################################################################
def rc(seq):
    rule = {"A":"T", "a":"T", "T":"A", "t":"A", "U":"A", "u":"A", \
            "G":"C", "g":"C", "C":"G", "c":"G", "Y":"R", "y": "R", \
            "R":"Y", "r":"Y", "N":"N", "n":"N", "-":"-", ".":"."}
    return "".join([rule[base] for base in reversed(seq)])
################################################################################



################################################################################
#####0. input files and amino acid codons
infiles = [
"RNAseq_gene_matrix_counts_D101_protein_coding_globalcu.txt",
"RNAseq_gene_matrix_counts_D103_protein_coding_globalcu.txt",
"RNAseq_gene_matrix_counts_D133_protein_coding_globalcu.txt",
"RNAseq_gene_matrix_counts_D32AD33_protein_coding_globalcu.txt",
"RNAseq_gene_matrix_counts_D32A_protein_coding_globalcu.txt",
"RNAseq_gene_matrix_counts_D33D35A_protein_coding_globalcu.txt",
"RNAseq_gene_matrix_counts_D33_protein_coding_globalcu.txt",
"RNAseq_gene_matrix_counts_D35A_protein_coding_globalcu.txt",
"RNAseq_gene_matrix_counts_D97D133_protein_coding_globalcu.txt",
"RNAseq_gene_matrix_counts_D97_protein_coding_globalcu.txt",
"RNAseq_gene_matrix_counts_WT1_protein_coding_globalcu.txt",
"Riboseq_gene_matrix_counts_D101_protein_coding_globalcu.txt",
"Riboseq_gene_matrix_counts_D103_protein_coding_globalcu.txt",
"Riboseq_gene_matrix_counts_D133_protein_coding_globalcu.txt",
"Riboseq_gene_matrix_counts_D32AD33_protein_coding_globalcu.txt",
"Riboseq_gene_matrix_counts_D32A_protein_coding_globalcu.txt",
"Riboseq_gene_matrix_counts_D33D35A_protein_coding_globalcu.txt",
"Riboseq_gene_matrix_counts_D33_protein_coding_globalcu.txt",
"Riboseq_gene_matrix_counts_D35A_protein_coding_globalcu.txt",
"Riboseq_gene_matrix_counts_D97D133_protein_coding_globalcu.txt",
"Riboseq_gene_matrix_counts_D97_protein_coding_globalcu.txt",
"Riboseq_gene_matrix_counts_WT1_protein_coding_globalcu.txt"]
codons = {
"UUU":"F","CUU":"L","AUU":"I","GUU":"V","UUC":"F","CUC":"L","AUC":"I","GUC":"V",
"UUA":"L","CUA":"L","AUA":"I","GUA":"V","UUG":"L","CUG":"L","AUG":"M","GUG":"V",
"UCU":"S","CCU":"P","ACU":"T","GCU":"A","UCC":"S","CCC":"P","ACC":"T","GCC":"A",
"UCA":"S","CCA":"P","ACA":"T","GCA":"A","UCG":"S","CCG":"P","ACG":"T","GCG":"A",
"UAU":"Y","CAU":"H","AAU":"N","GAU":"D","UAC":"Y","CAC":"H","AAC":"N","GAC":"D",
"UAA":"*","CAA":"Q","AAA":"K","GAA":"E","UAG":"*","CAG":"Q","AAG":"K","GAG":"E",
"UGU":"C","CGU":"R","AGU":"S","GGU":"G","UGC":"C","CGC":"R","AGC":"S","GGC":"G",
"UGA":"*","CGA":"R","AGA":"R","GGA":"G","UGG":"W","CGG":"R","AGG":"R","GGG":"G"}
#"aug":"m","uga":"U"}
aa1to3 = {"A":"Ala","C":"Cys","D":"Asp","E":"Glu","F":"Phe","G":"Gly","H":"His",
          "I":"Ile","K":"Lys","L":"Leu","M":"Met","N":"Asn","P":"Pro","Q":"Gln",
          "R":"Arg","S":"Ser","T":"Thr","V":"Val","W":"Trp","Y":"Tyr","U":"SeC",
          "m":"iMet"}
#match tRNAs to mRNAs using Fig. 5A from Bornelov 2019.
anticodons = {} #human tRNA anticodons.
"""

"""
wobbles = {"AGC":["GCU","GCC"], "ACG":["CGU","CGC"], "GCC":["GGU","GGC"],
           "GUG":["CAU","CAC"], "AAC":["GUU","GUC"], "AGA":["UCU","UCC"],
           "GAA":["UUU","UUC"], "GUU":["AAU","AAC"], "AAG":["CUU","CUC"],
           "AGG":["CCU","CCC"], "GUC":["GAU","GAC"], "GUA":["UAU","UAC"],
           "AGU":["ACU","ACC"], "AAU":["AUU","AUC"], "GCU":["AGU","AGC"],
           "GCA":["UGU","UGC"]}
wobblecodons = sum(list(wobbles.values()),[])
others = {rc(i).replace("T","U"):[i] for i in codons if i not in
          sum(list(wobbles.values()),[])}
decoding = {**wobbles, **others}
#print(len(decoding))
#for i in decoding: print(i)

aa3to1 = {i[1]:i[0] for i in aa1to3.items()} #print(aa3to1)
aas = list(sorted(set(list(codons.values())))) #amino acids
tri = sorted(list(codons.keys())) #trinucleotide
aa_codons = sorted([codons[codon]+"_"+codon for codon in codons])
################################################################################





################################################################################
#####1. convert input globalcu codon usage files to dictionaries.
globalcus = {}
for file in infiles:
    fileh = open(file, 'r')
    key = file.split("_")[0]+"_"+file.split("_")[4]; globalcus[key] = {}
    for line in fileh:
        record = line.split()
        if record[0] == "codons": continue
        globalcus[key][record[0]] = float(record[1])
    fileh.close()
################################################################################







################################################################################
#####2. calculate amino acid usage frequency. add up synonymous codons
globalaas = {} #print(globalcus["RNAseq_D35A"])
for k in globalcus: #k is sample KO name, e.g., "RNAseq_D35A"
    globalaas[k] = {}
    for aa_codon in aa_codons:
        aa = aa_codon[0]
        if aa not in globalaas[k]: globalaas[k][aa] = 0
        globalaas[k][aa] += globalcus[k][aa_codon]
for k in globalaas: #sum each file to 1
    total = sum(globalaas[k].values())
    globalaas[k] = {i:globalaas[k][i]/total for i in globalaas[k]}
for k in globalaas: #normalize against WT. k, e.g.: RNAseq_WT1
    intype = k.split("_")[0] #intype: RNAseq or Riboseq
    for i in globalaas[k]: globalaas[k][i] /= globalaas[intype+"_WT1"][i]
aamarkers = ["'$" + aa + "$'" for aa in aas]
#print(globalaas["RNAseq_D35A"])
################################################################################







################################################################################
#####3. normalize codon usage for each RNAseq and Riboseq dataset, vs. WT. 
for k in globalcus: #sum each file to 1
    total = sum(globalcus[k].values())
    globalcus[k] = {i:globalcus[k][i]/total for i in globalcus[k]}
    #for i in globalcus[k]: print(globalcus[k][i])
for k in globalcus: #normalize against WT. k, e.g.: RNAseq_WT1
    intype = k.split("_")[0] #intype: RNAseq or Riboseq
    for i in globalcus[k]: globalcus[k][i] /= globalcus[intype+"_WT1"][i]
################################################################################






################################################################################
#####4. extract anticodon and amino acid frequencies for tRNA levels from: 
#~/Documents/lulab/projects/snorna/riboseq/TE/20230125_RNAseq_tRNA_count.txt
tRNAfile = open("20230125_RNAseq_count.txt", 'r')
samples = tRNAfile.readline().split()[2:]
tRNAcounts = {i:{} for i in samples}
tRNAaacounts = {i:{} for i in samples} #key = samplename, WT, D32A, ...
tRNAac = {}
for line in tRNAfile:
    record = line.split()
    if record[1] != "tRNA": continue
    gene = record[0].split("-")[0][5:].split("_"); tRNA1 = aa3to1[gene[0]]
    for i in range(len(samples)):
        ac = gene[1].replace("T","U")
        if ac not in tRNAac: tRNAac[ac] = 0
        if ac not in tRNAcounts[samples[i]]: tRNAcounts[samples[i]][ac]=0
        tRNAcounts[samples[i]][ac] += float(record[i+2])
        if tRNA1 not in tRNAaacounts[samples[i]]:
            tRNAaacounts[samples[i]][tRNA1]=0
        tRNAaacounts[samples[i]][tRNA1] += float(record[i+2])
tRNAfile.close() #print(tRNAcounts); print(tRNAaacounts)
tRNAratios, tRNAaaratios = {}, {}
tRNAac = [i for i in tRNAac]
for s in samples:
    sW = sum(tRNAcounts["WT"].values()); sS = sum(tRNAcounts[s].values())
    tRNAratios[s]=[tRNAcounts[s][i]*sW/tRNAcounts["WT"][i]/sS if
                   tRNAcounts["WT"][i] else 0 for i in tRNAcounts["WT"]]
    aas22 = aas[1:] #ignore stop codon
    aW = sum(tRNAaacounts["WT"].values()); aS = sum(tRNAaacounts[s].values())
    tRNAaaratios[s]=[tRNAaacounts[s][i]*aW/tRNAaacounts["WT"][i]/aS
                     for i in aas22]
################################################################################ 






################################################################################
#####5. plot amino acid and codon scatters. 
endAU = [c for c in aa_codons if c[-1] in "AU"]
endGC = [c for c in aa_codons if c[-1] in "GC"]
endA = [c for c in aa_codons if c[-1] == "A"]
endU = [c for c in aa_codons if c[-1] == "U"]
endG = [c for c in aa_codons if c[-1] == "G"]
endC = [c for c in aa_codons if c[-1] == "C"]

aavarRNA  = {} #standard deviations of aa usage in RNAseq mRNAs. 
aavarRibo = {} #standard deviations of aa usage in Riboseq mRNAs.
#aavartRNA = {} #standard deviations of aa usage in RNAseq tRNAs.
for KO in globalcus:
    if "RNA" in KO:
        RNAseq = [globalcus[KO][i] for i in aa_codons] 
        RNAseqA = [globalcus[KO][i] for i in endA]
        RNAseqU = [globalcus[KO][i] for i in endU]
        RNAseqG = [globalcus[KO][i] for i in endG]
        RNAseqC = [globalcus[KO][i] for i in endC]
        RNAseqAU = [globalcus[KO][i] for i in endAU]
        RNAseqGC = [globalcus[KO][i] for i in endGC]
        RNAseqaa = [globalaas[KO][i] for i in aas]
        Riboseq = [globalcus["Ribo"+KO[3:]][i] for i in aa_codons]
        RiboseqA = [globalcus["Ribo"+KO[3:]][i] for i in endA]
        RiboseqU = [globalcus["Ribo"+KO[3:]][i] for i in endU]
        RiboseqG = [globalcus["Ribo"+KO[3:]][i] for i in endG]
        RiboseqC = [globalcus["Ribo"+KO[3:]][i] for i in endC]
        RiboseqAU = [globalcus["Ribo"+KO[3:]][i] for i in endAU]
        RiboseqGC = [globalcus["Ribo"+KO[3:]][i] for i in endGC]
        Riboseqaa = [globalaas["Ribo"+KO[3:]][i] for i in aas]
        aavarRNA[KO] = np.std(RNAseqaa)
        aavarRibo["Ribo"+KO[3:]] = np.std(Riboseqaa)
        
        """
        ##### plot codon frequencies for mRNAs vs. tRNAs. 
        #0. merge codon usage according to anticodon decoding rules.
        #1. Some synonymous codon tRNAs are too similar. We need unique reads. 
        #2. only the absolutely required wobble codons are included
        #3. iMet is ignored as all Met codons in mRNAs are counted together.
        #4. Sec is ignored for now.
        if "WT1" in KO: continue
        decodeRNAseq = [0 for i in range(len(tRNAac))]
        decodeRiboseq = [0 for i in range(len(tRNAac))]
        for i in range(len(tRNAac)):
            if tRNAac[i] in ["GAU","AUA"]: continue
            #ignore two intron containing tRNAs for Ile and Tyr
            for j in decoding[tRNAac[i]]:
                for aa_codon in globalcus[KO]:
                    if j in aa_codon:
                        decodeRNAseq[i]+=globalcus[KO][aa_codon]
                        decodeRiboseq[i]+=globalcus["Ribo"+KO[3:]][aa_codon]
        for i in range(len(tRNAac)):
            if tRNAac[i] in ["GAU","AUA"]: continue
            decodeRNAseq[i]/=len(decoding[tRNAac[i]])
            decodeRiboseq[i]/=len(decoding[tRNAac[i]])
        
        for i in list(range(len(tRNAratios[KO.split("_")[1]])))[::-1]:
            if tRNAratios[KO.split("_")[1]][i]<=0 or decodeRNAseq[i]<=0 or \
               decodeRiboseq[i]<=0:
                del tRNAratios[KO.split("_")[1]][i]
                del decodeRNAseq[i]; del decodeRiboseq[i]
        
        tRNAratioslog = [np.log2(i) for i in tRNAratios[KO.split("_")[1]]]
        decodeRNAseqlog = [np.log2(i) for i in decodeRNAseq]
        decodeRiboseqlog = [np.log2(i) for i in decodeRiboseq]
        fig, ax = plt.subplots()
        plt.plot(tRNAratioslog, decodeRNAseqlog, "bo")
        plt.xlim(-5,4)
        plt.ylim(-0.4,0.4)
        plt.vlines(x=0,ymin=-0.4, ymax=0.4, colors='k')
        plt.hlines(y=0,xmin=-5, xmax=4, colors='k')
        
        a, b = np.polyfit(tRNAratioslog, decodeRNAseqlog, 1)
        plt.plot(tRNAratioslog, np.array(tRNAratioslog)*a+b, "k-")
        plt.savefig("tRNAanticodons_vs_RNAseq"+KO+".pdf")
        pr = scipy.stats.pearsonr(tRNAratioslog,decodeRNAseqlog)
        print(KO, "r="+"{:.3f},".format(pr[0]),"p="+'%.2E' % Decimal(pr[1]))
        print("y="+"{:.3f}".format(a)+"x + "+"{:.3f}".format(b))
        plt.clf()
        """


        

        """
        ##### plot codon frequencies for A,U,G,C,AU,GC-ending ones
        fig, ax = plt.subplots(figsize=(1,1))
        a, b = np.polyfit(RNAseq, Riboseq, 1)
        plt.plot(RNAseq, np.array(RNAseq)*a+b, "k-")
        print(KO)
        print("y="+"{:.3f}".format(a)+"x + "+"{:.3f}".format(b))
        print("r="+"{:.3f},".format(scipy.stats.pearsonr(RNAseq,Riboseq)[0]),
              "p="+'%.2E' % Decimal(scipy.stats.pearsonr(RNAseq,Riboseq)[1]))
        plt.plot(RNAseqAU, RiboseqAU, "bo", markersize=0.5, alpha=0.75)
        plt.plot(RNAseqGC, RiboseqGC, "ro", markersize=0.5, alpha=0.75)
        #plt.plot(RNAseqA, RiboseqA, "bo", markersize=0.5, alpha=0.75)
        #plt.plot(RNAseqU, RiboseqU, "b+", markersize=0.5, alpha=0.75)
        #plt.plot(RNAseqG, RiboseqG, "ro", markersize=0.5, alpha=0.75)        
        #plt.plot(RNAseqC, RiboseqC, "r+", markersize=0.5, alpha=0.75)
        
        ##### format the plot, set individual parameters. 
        plt.xticks([0.7,0.8,0.9,1.0,1.1,1.2,1.3])
        plt.yticks([0.7,0.8,0.9,1.0,1.1,1.2,1.3])
        plt.xlim(0.7,1.3) #0.85-1.15 for aa frequencies, 0.7-1.3 for codons
        plt.ylim(0.7,1.3)
        ax.set_aspect('equal')
        plt.savefig(KO+"_AUGCbias.pdf")
        plt.clf()
        RNAratio = np.mean(RNAseqGC)/np.mean(RNAseqAU)
        RNAttest = scipy.stats.ttest_ind(RNAseqGC,RNAseqAU)[1]#2-side, same var
        Riboratio = np.mean(RiboseqGC)/np.mean(RiboseqAU)
        Ribottest = scipy.stats.ttest_ind(RiboseqGC,RiboseqAU)[1]          
        print("Diff in RNA-seq and Ribo-seq:") 
        print("RNA-seq\nratio = "+"{:.3f},".format(RNAratio),
              "p = "+'%.2E' % Decimal(RNAttest))
        print("Ribo-seq\nratio = "+"{:.3f},".format(Riboratio),
              "p = "+'%.2E' % Decimal(Ribottest))
        """

        
        
        ##### plot amino acid usage in RNAseq vs. Riboseq
        fig, ax = plt.subplots(figsize=(1.5,1.5))
        a, b = np.polyfit(RNAseqaa, Riboseqaa, 1)
        plt.plot(RNAseqaa, np.array(RNAseqaa)*a+b, "k-")
        for i in range(len(aas)):
            plt.plot(RNAseqaa[i], Riboseqaa[i], color="black", ls="None",
                     marker = aamarkers[i][1:-1], markersize=3, alpha=1)
        print(KO)
        print("y="+"{:.3f}".format(a)+"x + "+"{:.3f}".format(b))
        print("r="+"{:.3f},".format(scipy.stats.pearsonr\
                                      (RNAseqaa,Riboseqaa)[0]),
              "p="+'%.2E' % Decimal(scipy.stats.pearsonr\
                                      (RNAseqaa,Riboseqaa)[1]))
        ##### format the plot, set individual parameters. 
        #plt.xticks([0.85,0.90,0.95,1.00,1.05,1.10,1.15])
        #plt.yticks([0.85,0.90,0.95,1.00,1.05,1.10,1.15])
        plt.xlim(0.9,1.1) #0.85-1.15 for aa frequencies, 0.7-1.3 for codons
        plt.ylim(0.9,1.1)
        ax.set_aspect('equal')
        plt.savefig(KO+"_aabias_0.1range.pdf")
        plt.clf()
        


        """
        ##### plot amino acid usage in tRNAs vs. RNAseq or Riboseq
        print(KO)
        if KO == "Riboseq_WT1": continue
        x = [np.log2(i) for i in tRNAaaratios[KO.split("_")[1]]]
        y = [np.log2(i) for i in Riboseqaa[1:]] #aa usage from RNAseq/Riboseq
        fig, ax = plt.subplots(figsize=(1.5,1.5))
        a, b = np.polyfit(x, y, 1)
        plt.plot(x, np.array(x)*a+b, "k-")
        plt.vlines(x=0,ymin=-0.2, ymax=0.2, colors='k')
        plt.hlines(y=0,xmin=-6, xmax=3, colors='k')
        print("y="+"{:.3f}".format(a)+"x+"+"{:.3f}".format(b))
        print("r="+"{:.3f},".format(scipy.stats.pearsonr(x,y)[0]),
              "p="+'%.2E' % Decimal(scipy.stats.pearsonr(x,y)[1]))
        for i in range(len(x)): 
            plt.plot(x[i],y[i],color="k", ls="None",
                     marker = aamarkers[i+1][1:-1], markersize=3, alpha=1)
        #for i in range(len(aas)-1):print("\t".join([aas[i+1],str(x[i]),str(y[i])]))
        plt.xlim(-6,3); plt.ylim(-0.2,0.2)
        plt.savefig("tRNA_vs_"+"Ribo"+KO[3:]+"_aausage.pdf"); plt.clf()
        """
################################################################################








################################################################################
#####5. plot variations in aa usage:
#KOlist=['D32A','D33','D32AD33','D35A','D33D35A','D97','D133','D97D133',
#        'D101','D103']
#print(globalcus.keys())
#for KO in aavarRNA: print(KO, aavarRNA["RNA_"+KO], aavarRibo["Ribo_"+KO])
################################################################################





"""
################################################################################
#####5. pseudo time series to show codon and aa usage from single and double KOs
fig, (ax1,ax2) = plt.subplots(1,2,figsize=(2,1))
cindex = list(range(len(codons)))
cseries1 = [] #codon, series1 for WT-D97-D97D133
cseries2 = [] #codon, series2 for WT-D133-D97D133
for i in cindex:
    intype = "RNAseq_"
    #intype = "Riboseq_"
    WTcu = globalcus[intype+"WT1"][aa_codons[i]]
    D97cu = globalcus[intype+"D97"][aa_codons[i]]
    D133cu = globalcus[intype+"D133"][aa_codons[i]]
    D97D133cu = globalcus[intype+"D97D133"][aa_codons[i]]
    cseries1.append([0,np.log2(D97cu/WTcu),np.log2(D97D133cu/WTcu)])
    cseries2.append([0,np.log2(D133cu/WTcu),np.log2(D97D133cu/WTcu)])
for i in cindex: ax1.plot([0,1,2],cseries1[i]); ax2.plot([0,1,2],cseries2[i])
ax1.set_ylim(-0.4,0.4); ax2.set_ylim(-0.4,0.4)
plt.savefig("20230125_Riboseq_ReadCounts_D97D133_tRNA_codon_pseudoseries.pdf")
plt.clf()

fig, (ax1,ax2) = plt.subplots(1,2,figsize=(2,1))
aindex = list(range(len(aas)))
aseries1 = [] #amino acid, series1 for WT-D97-D97D133
aseries2 = [] #amino acid, series2 for WT-D133-D97D133
for i in aindex:
    intype = "RNAseq_"
    #intype = "Riboseq_"
    WTau = globalaas[intype+"WT1"][aas[i]]
    D97au = globalaas[intype+"D97"][aas[i]]
    D133au = globalaas[intype+"D133"][aas[i]]
    D97D133au = globalaas[intype+"D97D133"][aas[i]]
    aseries1.append([0,np.log2(D97au/WTau),np.log2(D97D133au/WTau)])
    aseries2.append([0,np.log2(D133au/WTau),np.log2(D97D133au/WTau)])
for i in aindex: ax1.plot([0,1,2],aseries1[i]); ax2.plot([0,1,2],aseries2[i])
ax1.set_ylim(-0.2,0.2); ax2.set_ylim(-0.2,0.2)
plt.savefig("20230125_Riboseq_ReadCounts_D97D133_tRNA_aa_pseudoseries.pdf")
plt.clf()

fig, ax = plt.subplots(figsize=(1,1))
D97std = np.std([i[1] for i in aseries1])
D133std = np.std([i[1] for i in aseries2])
D97D133std = np.std([i[2] for i in aseries1])
plt.bar([0,1,2],[D97std,D133std,D97D133std], width=0.5)
plt.savefig("a.pdf")
################################################################################
"""


