"""
aligntrim.py
After CLUSTAL alignment of sequences, e.g. from Rfam collection, trim alignments
so that only positions repreented by >= e.g. 30% seqs are kept.
cd /Users/lu/Documents/lulab/projects/snorna/Rfam
python ~/Documents/scripts/snorna/aligntrim.py RF01290.fa RF01290trim.fa
Note: In the seed alignments for RF01290, the last CUGA was not aligned

"""
import sys
if len(sys.argv)<4:
    print("Usage: python aligntrim.py infa outfa covlimit")
    print("infa: input aligned fasta sequences")
    print("covlimit: minimal coverage of each position, default 0.3")
    sys.exit()

infa = open(sys.argv[1], 'r')
outfa = open(sys.argv[2], 'w')
covlimit = float(sys.argv[3])

seqdict = {} #dictionary to store input sequences. 
name = ''
for line in infa:
    if line[0]==">": name=line; seqdict[name]=[]
    else: seqdict[name].append(line.strip('\n'))
infa.close()
for name in seqdict: seqdict[name]=''.join(seqdict[name])

covered = [] #list of positions with >=covlimit coverage
seqs = list(seqdict.values())
seqlen = len(seqs[0])
nseqs = len(seqs)
for i in range(seqlen):
    if list(zip(*seqs))[i].count('-')<(1-covlimit)*nseqs: covered.append(i)
seqlentrim = 0
for name in seqdict:
    seqtrim=''.join([seqdict[name][i] for i in covered])
    outfa.write(name+seqtrim+'\n')
    seqlentrim = len(seqtrim)
outfa.close()
print("Length of trimmed sequences:", seqlentrim)
