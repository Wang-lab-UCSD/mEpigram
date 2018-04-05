import os 
from decimal import Decimal
from sys import argv

"""Usage:
motifpvalfifler.py memefile enrichementfile pvalcutoff outfile
"""

#wd = "mepigram/CEBPB_0.1_1000.faa.mepi.results.results/"

#memefile = wd + 'motifs.mepigram.meme'
#enrichmentfile = wd + 'enrichments.tsv'
#pvalcutoff = 0.001
#outfile = wd + "motifs.mepigram.final.meme"

memefile = argv[1]
enrichmentfile = argv[2]
pvalcutoff = float(argv[3])
outfile = argv[4]

# read enrichement file:
lines = open(enrichmentfile).readlines()

#print lines[0]
#for line in lines[:10]:
#    print line.strip()
lines = lines[1:]
pvalues = []
motifs = []
for line in lines:
    tmp = line.strip().split()
    pvalues += [float(tmp[1])]
    motifs += tmp[0]

sorted_lines =[x for _,x in sorted(zip(pvalues,lines))]
#for line in sorted_lines[:5]:
#    print line.strip()
    
#print "-"*30
    
# rename each motif by pval rank and enrichement score:
oldnames = []
newnames = []
for i in range(len(sorted_lines)):
    line = sorted_lines[i]
    tmp = line.strip().split()
    oldname = tmp[0]
    if float(tmp[1]) > pvalcutoff:
        continue
    pval = '%.2E' % Decimal(tmp[1])
    newname = str(i) + '_' + str(round(float(tmp[-1]),4)) + '_' + pval
    print newname
    oldnames += [oldname]
    newnames += [newname]
# get get only the motifs in gotten
#print "-"*40
pwms = open(memefile).read().split("MOTIF")
header = pwms[0]
newpwms = []
for i in range(len(oldnames)):
    oldname = oldnames[i]
    #print motif
    for pwm in pwms:
        if oldname in pwm:
            newpwms += ["MOTIF\t" + pwm.replace(oldname,newnames[i]).strip()+'\n']
#print header + '\n'.join(newpwms)

out = open(outfile,'w')
out.write(header + '\n'.join(newpwms))
out.close()

#rename the meme files accordingly
#os.system("mv %s %s" %(memefile,memefile+'.tmp'))
#os.system("mv %s %s" %(outfile,memefile)