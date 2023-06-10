"""
goterm2goid.py. From GSEA analysis output, find GO_ID from the term names. 

example command. 
cd /Users/lu/gsea_home/output/20230125_TE
python ~/Documents/scripts/snorna/goterm2goid.py \
20230125_TE_D133_c5.all/gsea_report*neg*.tsv 0.05 \
20230125_TE_D133_c5.all_neg0.05.txt

cd /Users/lu/gsea_home/output/20230318
python ~/Documents/scripts/snorna/goterm2goid.py \
EB_D133KO5_vs_mES_D133KO5/gsea_report*neg*.tsv 0.05 \
EB_D133KO5_vs_mES_D133KO5_neg0.05.txt

cd /Users/lu/gsea_home/output/20230511_dKO
python ~/Documents/scripts/snorna/goterm2goid.py \
mESC_D97D133_rnaseq_norm_logratio/gsea_report*neg*.tsv 0.05 \
mESC_D97D133_rnaseq_norm_logratio_neg0.05.txt

In R:
setwd("/Users/lu/gsea_home/output/20230511_dKO/")
library(simplifyEnrichment)
set.seed(0)
infile <- "./mESC_D97D133_rnaseq_norm_logratio_neg0.05.txt"
goids = unlist(read.table(infile))
mat = GO_similarity(goids, ont="BP")
simplifyGO(mat)


"""
import sys
pcut = float(sys.argv[2])
gofile = open("/Users/lu/Documents/scripts/snorna/goterms.txt", 'r')
godict = {}
for line in gofile:
    goid, goterm = line.split("\t")
    godict[goterm[1:-2].lower()] = goid.strip("\"")
gofile.close()
print("GO dictionary size:", len(godict))
#print(godict["nephron epithelium development"])

goidlist = []
gseafile = open(sys.argv[1], 'r')
for line in gseafile:
    record = line.split()
    if record[0] == "NAME" or record[-7] == "---": continue
    nompval = float(record[-7])
    if nompval <= pcut:
        name = " ".join(record[1].lower().split("_")[1:])
        if name in godict: goidlist.append(godict[name])
gseafile.close()
for goid in goidlist: print(goid)
outfile = open(sys.argv[3], 'w')
outfile.write("\n".join(goidlist)+'\n')
outfile.close()
