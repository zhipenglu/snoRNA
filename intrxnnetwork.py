"""
intrxnnetwork.py
plotting the connections between two groups of RNAs based on experimental data,
like PARIS. Examples include the snoRNAs and tRNAs. The two groups of RNAs are
first collected into dictionaries and sorted by numbers of reads connecting
each RNA. The connections are scaled based on square root of read numbers. 

Example command
cd /Users/lu/Dropbox/_2022_sno/chimera
python ~/Documents/scripts/snorna/intrxnnetwork.py sno_t_intrxn.txt \
snoRNA tRNA 2

"""

################################################################################
import sys
import matplotlib.pyplot as plt
if len(sys.argv) < 4:
    print("Usage: python intrxnnetwork.py infile RNA1 RNA2 cxlimit")
    sys.exit()
infile = open(sys.argv[1], 'r')
RNA1file = open(sys.argv[2], 'w')
RNA2file = open(sys.argv[3], 'w')
cxlimit = int(sys.argv[4]) #minimal No. of reads in each connection
################################################################################





################################################################################
#####build 3 dictionaries: RNA1, RNA2, intrxndict
RNA1dict = {}
RNA2dict = {}
intrxndict = {}
for line in infile:
    RNA1, RNA2, connections = tuple(line.split())
    if RNA1 not in RNA1dict: RNA1dict[RNA1] = 0
    if RNA2 not in RNA2dict: RNA2dict[RNA2] = 0
    RNA1dict[RNA1] += int(connections)
    RNA2dict[RNA2] += int(connections)
    intrxndict[(RNA1,RNA2)] = int(connections)
RNA1num = len(RNA1dict)
RNA2num = len(RNA2dict)
print("Numbers of RNA1 and RNA2:", RNA1num, RNA2num)
################################################################################





################################################################################
#####filter the dictionaries by No. of reads
tempdict = {}
for key in intrxndict:
    # ("SNORD97" in key[0] or "SNORD133" in key[0])
    # ("SNORD32" in key[0] or "SNORD33" in key[0] or "SNORD34" in key[0] or "SNORD35" in key[0]) 
    if ("His_GTG" in key[1]) and intrxndict[key] >=cxlimit: #"
        tempdict[key] = intrxndict[key] 
intrxndict = tempdict
RNA1filter = [i[0] for i in intrxndict]
RNA2filter = [i[1] for i in intrxndict]
RNA1dicttemp = {}; RNA2dicttemp = {}
for RNA in RNA1dict:
    if RNA in RNA1filter: RNA1dicttemp[RNA] = RNA1dict[RNA]
for RNA in RNA2dict: 
    if RNA in RNA2filter: RNA2dicttemp[RNA] = RNA2dict[RNA]
RNA1dict = RNA1dicttemp
RNA2dict = RNA2dicttemp
RNA1num = len(RNA1dict)
RNA2num = len(RNA2dict)
print("Numbers of RNA1 and RNA2 passing filter:", RNA1num, RNA2num)
print("Number of interactions passing filter:", len(intrxndict))
################################################################################





################################################################################
#####calculate the coordinates of each RNA
#if RNA1 has fewer elements, we need to switch RNA1 and RNA2. ???
#RNA1 on top row, RNA2 on bottom row
#recalculate numbers of reads supporting each RNA.
RNA1dict = {}; RNA2dict = {}
for key in intrxndict:
    if key[0] not in RNA1dict: RNA1dict[key[0]] = 0
    RNA1dict[key[0]] += intrxndict[key]
    if key[1] not in RNA2dict: RNA2dict[key[1]] = 0
    RNA2dict[key[1]] += intrxndict[key]
RNA1coords = {} 
RNA2coords = {}
i = 1; j = 1
RNA1list = list(sorted(RNA1dict, key=RNA1dict.get, reverse=True))
for RNA in RNA1list: RNA1coords[RNA] = (i, 2); i+=1
RNA2list = list(sorted(RNA2dict, key=RNA2dict.get, reverse=True))
for RNA in RNA2list: RNA2coords[RNA] = (j*RNA1num/RNA2num, 1); j+=1
#print(RNA1list, RNA2list)
################################################################################





################################################################################
#####plot the lines connecting all the RNAs. 
fig, ax = plt.subplots(figsize=(10,1))
for key in intrxndict:
    dot1 = RNA1coords[key[0]]
    dot2 = RNA2coords[key[1]]
    xs = [dot1[0], dot2[0]]
    ys = [dot1[1], dot2[1]]
    plt.plot(xs, ys, linewidth= intrxndict[key]**0.5/3) #c='b'
plt.savefig(sys.argv[1]+".pdf")

RNA1out = []; RNA2out = []
for RNA in RNA1list: RNA1out.append(RNA+', '+str(RNA1dict[RNA]))
for RNA in RNA2list: RNA2out.append(RNA+', '+str(RNA2dict[RNA]))
RNA1file.write('\n'.join(RNA1out))
RNA1file.close()
RNA2file.write('\n'.join(RNA2out))
RNA2file.close()
################################################################################










