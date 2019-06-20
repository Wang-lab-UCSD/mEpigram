from numpy.random import choice
import numpy as np
from sys import argv


#use 1st-order Markov chain to make the background

def makeModel(seqs,alphabet):
    """read through the sequences and make a MC model"""
    #use the alphabet to make all of the possible combination
    dimer_counts = {}
    for a in alphabet:
        dimer_counts[a] = {}
        for b in alphabet:
            dimer_counts[a][b] = 0 
    
    #read all of the di-mers
    for seq in seqs:
        seq = seq.upper()
        for i in range(len(seq) - 2 + 1):
            try:
                dimer_counts[seq[i]][seq[i + 1]] += 1
            except KeyError:
                continue

            
            
    # normalize the dict
    for a in dimer_counts:
        total = 0.0
        for b in dimer_counts[a]:
            total += dimer_counts[a][b]
            
        for b in alphabet:
            try:
                dimer_counts[a][b] = dimer_counts[a][b]/total #normalized
            except ZeroDivisionError: #total is zero
                dimer_counts[a][b] = 0.0
                
    model = {}
    for a in alphabet:
        model[a] = []
        for b in alphabet:
            model[a] += [dimer_counts[a][b]]

    return model



def makeRandomSeqs(mc_model, alphabet, seed, seqlengths): #takes a seed to make reproducible numbers
    """takes MC transitional matrix and produces a random sequence """
    np.random.seed(0) ; np.random.rand(seed)

    outseqs = []
    for l in seqlengths:
        
        newseq = [choice(["A","C","G","T"])] #initialize the seq, minus E and F
        for i in range(l - 1):
            newseq += [choice(alphabet, p = mc_model[newseq[-1]])]
        outseqs += [''.join(newseq)]
        
    return outseqs
    
    
    
def main():
    """ 
    IN: takes the set of fasta sequences, a seed, a number of copies to make
    OUT: a fasta file containing the sequences generated randomly from the mc_model made using the input fasta file
    """
    seed = int(argv[1]) #seed for random process
    fastafile = argv[2] #
    copynum = int(argv[3]) # number of copies to make
    outfile = argv[4]
    
    alphabet = ["A","G","C","T","E","F"]
    
    
    #read the sequences and determine their lengths
    lines = open(fastafile).read().split(">")[1:]
    seqs = []
    for l in lines:
        seqs += [''.join(l.split()[1:])]
    
    seqlengths = []
    for s in seqs:
        seqlengths += [len(s)]
    seqlengths = seqlengths*copynum
    
    print "making markov model"
    mc_model = makeModel(seqs,alphabet)
    print "generating random sequences"
    newseqs = makeRandomSeqs(mc_model, alphabet, seed, seqlengths)
    
    outfile = open(outfile, 'w')
    c = 0
    for s in newseqs:
        line = ">seq"+str(c)+'\n'+s+'\n'
        outfile.write(line)
        c += 1
    outfile.close()
    

if __name__ == "__main__":
    main()
    
    
    
    
    
    
   
