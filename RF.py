"""
RF.py, a collection of common functions and data for RF analysis. 
working directary: /Users/lu/Documents/lulab/projects/snorna/tRFs/
associated file: "/projects/hg38/hg38-tRNAs/hstRNA_pos_adjust.txt"

1. tRNAtable: a dictionary converting genomic to standard coordinates
2. tRNAadjusttable, a function to make the above table
3. bam2RFcov: a function to calculate coverage
4. RFbams: a dictionary of bam files for RNA fragments.
"""

import pysam



################################################################################
##### 1. build a conversion table (dictionary) for all tRNAs.
af = "/Users/lu/Documents/lulab/projects/hg38/hg38-tRNAs/hstRNA_pos_adjust.txt"
humantRNAintrons = {
    ("tRNA_Arg_TCT","chr1:93847573-93847657"): [38,-37,12],
    ("tRNA_Arg_TCT","chr17:8120925-8121012"): [38,-37,15],
    ("tRNA_Arg_TCT","chr9:128340076-128340166"): [38,-37,18],
    ("tRNA_Arg_TCT","chr11:59551294-59551379"): [38,-37,13],
    ("tRNA_Arg_TCT","chr6:27562184-27562270"): [38,-37,14],
    ("tRNA_Leu_CAA","chr6:28896223-28896328"): [39,-46,23],
    ("tRNA_Leu_CAA","chr6:28941053-28941157"): [39,-46,22],
    ("tRNA_Leu_CAA","chr6:27605638-27605745"): [39,-47,24],
    ("tRNA_Leu_CAA","chr6:27602569-27602675"): [39,-47,23],
    ("tRNA_Leu_CAA","chr1:248873855-248873960"): [39,-47,22],
    ("tRNA_Ile_TAT","chr19:39412168-39412260"): [39,-37,19],
    ("tRNA_Ile_TAT","chr2:42810536-42810628"): [39,-37,19],
    ("tRNA_Ile_TAT","chr6:27020346-27020439"): [39,-37,20],
    ("tRNA_Ile_TAT","chr6:27631421-27631514"): [39,-37,20],
    ("tRNA_Ile_TAT","chr6:28537590-28537683"): [39,-37,20],
    ("tRNA_Tyr_ATA","chr2:218245826-218245918"): [38,-37,20],
    ("tRNA_Tyr_GTA","chr6:26568858-26568948"): [38,-37,18],
    ("tRNA_Tyr_GTA","chr2:27050782-27050870"): [38,-37,16],
    ("tRNA_Tyr_GTA","chr6:26577104-26577192"): [38,-37,16],
    ("tRNA_Tyr_GTA","chr14:20657464-20657557"): [38,-37,21],
    ("tRNA_Tyr_GTA","chr8:66113367-66113459"): [38,-37,20],
    ("tRNA_Tyr_GTA","chr8:66113988-66114076"): [38,-37,16],
    ("tRNA_Tyr_GTA","chr14:20653099-20653192"): [38,-37,21],
    ("tRNA_Tyr_GTA","chr14:20663192-20663285"): [38,-37,21],
    ("tRNA_Tyr_GTA","chr14:20683273-20683361"): [38,-37,16],
    ("tRNA_Tyr_GTA","chr6:26594874-26594962"): [38,-37,16],
    ("tRNA_Tyr_GTA","chr14:20659958-20660051"): [38,-37,21],
    ("tRNA_Tyr_GTA","chr6:26575570-26575659"): [38,-37,17],
    ("tRNA_Tyr_GTA","chr8:65697297-65697384"): [38,-36,15]}
def tRNAadjusttable(adjustfile,tRNAintrons): 
    """
    Input:
    adjustfile: name of file containing liftover relationship for coordinates.
    format: tRNA_Ala_AGC-chr6:28795964:28796035	1 1
    tRNAintrons: Intron positions for intron-containing genes
    format: (tRNA,genomic_coordinate):[start,end,length] #1-based intron coords
    Output:
    adjusttable: (tRNA,genomic_coordinate,pos),adjustpos
    """    
    tRNAtable = {}
    infile = open(adjustfile, 'r')
    for line in infile:
        tRNA_coord,pos,adjustpos = line.split()
        tRNA,coord = tRNA_coord.split('-'); pos = int(pos)
        coord = coord.split(":")
        RNAlen = int(coord[2])-int(coord[1])+1
        coord = coord[0]+":"+coord[1]+"-"+coord[2]
        if (tRNA,coord) in tRNAintrons:
            introninfo = tRNAintrons[(tRNA,coord)]
            for i in range(introninfo[0],introninfo[0]+introninfo[2]):
                tRNAtable[(tRNA,coord,i)] = "intron"
            if pos>=introninfo[0]: pos+=introninfo[2]
        tRNAtable[(tRNA,coord,pos)] = adjustpos
    infile.close()
    
    #process -20 and +20 extra coordinates.
    tRNAgenes = list(set([k[:2] for k in tRNAtable.keys()]))
    #print(tRNAgenes)
    for tRNAgene in tRNAgenes:
        tRNA,coord = tRNAgene
        coord = [int(i) for i in coord.split(":")[1].split("-")]
        RNAlen = int(coord[1])-int(coord[0])+1
        for i in range(-19,1): tRNAtable[tRNAgene+(i,)] = str(i)
        for i in range(1,21): tRNAtable[tRNAgene+(i+RNAlen,)] = str(i+73)
    return tRNAtable

tRNAtable = tRNAadjusttable(af, humantRNAintrons)
#print("tRNA conversion table length: ", len(tRNAtable))
################################################################################





################################################################################
##### 2. A function to convert bam alignments into an RF coverage dictionary
#copied from RFcompareall.py
def bam2RFcov(bam, RNA, coord):
    """
    input:  indexed bam file name and RNA coordinates ("chr1", 100, 200)
    output: dictionary of (name, coord, (start,end)): count
    """
    RNAlen = coord[2]-coord[1]+1
    RFcovdict = {}
    for i in range(-20, RNAlen+20):
        for j in range(-20, RNAlen+20): RFcovdict[RNA,coord,(i,j)] = 0
    FLAG = 0
    samfile = pysam.AlignmentFile(bam, "rb")
    for line in samfile.fetch(*coord):
        FLAG,RNAME,POS,_,CIGAR = str(line).split()[1:6]
        POS = int(POS)
        cigars = [1 if c in "IDNSHPX=" else 0 for c in CIGAR]
        if sum(cigars): continue
        if POS<=coord[1]-20 or POS+int(CIGAR[:-1])>=coord[2]+20: continue
        RFcovdict[RNA,coord,(POS-coord[1],POS-coord[1]+int(CIGAR[:-1]))]+=1
    samfile.close()
    STRAND = '-' if '{0:012b}'.format(int(FLAG))[-5]=='1' else '+'
    if STRAND == "-":
        t = {}
        for i in range(-20, RNAlen+20):
            for j in range(-20, RNAlen+20):
                t[RNA,coord,(i,j)]=RFcovdict[RNA,coord,(RNAlen-j-1,RNAlen-i-1)]
        RFcovdict = t
    return RFcovdict
################################################################################




################################################################################
#30 files, excluding the file WT_3 and the D133KO
# D133 samples missed the longer fragments, leading to enrichment artifacts. 
#"HEK293_D133KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam":"D133KO",
#"HEK293_D133KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam":"D133KO",
"""
RFbams = {
"HEK293_D101KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam":"D101KO",
"HEK293_D101KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam":"D101KO",
"HEK293_D103KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam":"D103KO",
"HEK293_D103KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam":"D103KO",
"HEK293_D32AKO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam":"D32AKO",
"HEK293_D32AKO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam":"D32AKO",
"HEK293_D32A_33KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam":"D32AD33KO",
"HEK293_D32A_33KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam":"D32AD33KO",
"HEK293_D33KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam":"D33KO",
"HEK293_D33KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam":"D33KO",
"HEK293_D33_D35AKO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam":"D33D35AKO",
"HEK293_D33_D35AKO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam":"D33D35AKO",
"HEK293_D35AKO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam":"D35AKO",
"HEK293_D35AKO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam":"D35AKO",
"HEK293_D97KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam":"D97KO",
"HEK293_D97KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam":"D97KO",
"HEK293_D97_133KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam":"D97D133KO",
"HEK293_D97_133KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam":"D97D133KO",
"HEK293_WT_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam":"WT",
"HEK293_WT_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam":"WT",
"HEK293_siCtrl_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam":"siCtrl",
"HEK293_siCtrl_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam":"siCtrl",
"HEK293_siDKC1_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam":"siDKC1",
"HEK293_siDKC1_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam":"siDKC1",
"HEK293_siDKC2_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam":"siDKC1",
"HEK293_siDKC2_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam":"siDKC1",
"HEK293_siFBL1_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam":"siFBL",
"HEK293_siFBL1_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam":"siFBL",
"HEK293_siFBL2_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam":"siFBL",
"HEK293_siFBL2_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam":"siFBL"}
"""
RFbams = "HEK293_D101KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D101KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam,\
HEK293_D103KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D103KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam,\
HEK293_D32AKO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D32AKO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam,\
HEK293_D32A_33KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D32A_33KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam,\
HEK293_D33KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D33KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam,\
HEK293_D33_D35AKO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D33_D35AKO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam,\
HEK293_D35AKO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D35AKO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam,\
HEK293_D97KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D97KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam,\
HEK293_D97_133KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D97_133KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam,\
HEK293_WT_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_WT_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam,\
HEK293_siCtrl_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_siCtrl_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam,\
HEK293_siDKC1_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_siDKC1_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam,\
HEK293_siDKC2_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_siDKC2_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam,\
HEK293_siFBL1_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_siFBL1_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam,\
HEK293_siFBL2_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_siFBL2_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam"
################################################################################





################################################################################
##### 1. input bam files, comma-separated lists
siCtrlbams = "\
HEK293_siCtrl_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_siCtrl_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam"
siFBLbams = "\
HEK293_siFBL1_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_siFBL1_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam,\
HEK293_siFBL2_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_siFBL2_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam"
siDKC1bams = "\
HEK293_siDKC1_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_siDKC1_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam,\
HEK293_siDKC2_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_siDKC2_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam"
WTbams = "\
HEK293_WT_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_WT_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam"
D101KObams = "\
HEK293_D101KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D101KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam"
D103KObams = "\
HEK293_D103KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D103KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam"
D32AKObams = "\
HEK293_D32AKO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D32AKO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam"
D32AD33KObams = "\
HEK293_D32A_33KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D32A_33KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam"
D33KObams = "\
HEK293_D33KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D33KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam"
D33D35AKObams = "\
HEK293_D33_D35AKO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D33_D35AKO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam"
D35AKObams = "\
HEK293_D35AKO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D35AKO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam"
D97KObams = "\
HEK293_D97KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D97KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam"
D97D133KObams = "\
HEK293_D97_133KO_ROS_tRFs_rep1Aligned.sortedByCoord.out.bam,\
HEK293_D97_133KO_ROS_tRFs_rep2Aligned.sortedByCoord.out.bam"
################################################################################


