"""
compare_rRNA_snRNA_tRNA.py. 2022-06-21. Zhipeng Lu
plot a comparison of snoRNA targets across rRNAs, tRNAs and snRNAs. 

Input files:
snoRNA_rRNA_DGbedpe_Atlas_Plexy.txt  102 snoRNAs
sno_t_intrxn.txt  269 snoRNAs

cd /Users/lu/Dropbox/_2022_sno/chimera/
python ~/Documents/scripts/snorna/compare_rRNA_snRNA_tRNA.py \
snoRNA_rRNA_DGbedpe_Atlas_Plexy.txt sno_t_intrxn.txt

"""
import sys
import matplotlib.pyplot as plt

rRNAf = open(sys.argv[1], 'r')
tRNAf = open(sys.argv[2], 'r')

sno_r = {} #snoRNAs and trans reads connecting to rRNAs
for line in rRNAf:
    record = line.split()
    if record[6]=="DG": continue
    r, sno, _, _ = record[6].split(",")
    n = int(record[7])
    if sno not in sno_r: sno_r[sno] = n
    else: sno_r[sno] += n
rRNAf.close()

sno_t = {} #snoRNAs and trans reads connecting to tRNAs
for line in tRNAf:
    sno, t, n = line.split()
    n = int(n)
    if sno not in sno_t: sno_t[sno] = n
    else: sno_t[sno] += n
tRNAf.close()

snodict = {}
snolist = []
for sno in sno_r: snodict[sno] = [sno_r[sno],0]
for sno in sno_t:
    if sno in snodict: snodict[sno][1] = sno_t[sno]
    else: snodict[sno] = [0, sno_t[sno]]

#print(len(sno_r), len(sno_t))
snolist = list(snodict.items())
snolist.sort(key=lambda x: sum(x[1]), reverse=True)
snolist = [i for i in snolist if sum(i[1])>=1]
maxn = max([sum(i[1]) for i in snolist])

fig = plt.figure(figsize=(10,3))
ax = fig.add_subplot(2, 1, 1)
plt.plot([i[1][0] for i in snolist],'bo', markersize=1)
plt.plot([i[1][1] for i in snolist],'ro', markersize=1)
for i in range(len(snolist)): #add a vertical line to connect them. 
    plt.vlines(x=i,ymin=0,ymax=max(snolist[i][1]),colors='gray',ls=':',lw=0.2)
    if snolist[i][1][0]*snolist[i][1][1]:
        plt.vlines(x=i,ymin=min(snolist[i][1]), ymax=max(snolist[i][1]),
                   colors='gray', ls='-', lw=0.5) #label=''

plt.ylim([0.7,maxn*1.5])
ax.set_yscale('log')
plt.savefig("compare_rRNA_snRNA_tRNA.pdf")
for i in snolist: print(i[1][1])







