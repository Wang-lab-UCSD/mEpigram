#!/usr/bin/env python
'''

Takes fasta file, motif file, background, seed if required, and outputfilename
calculate the number of of sequences in fasta file and length
dinucleotide shuffle them if they are longer than a threshold
output into a tmp.DS.faa file for motifscanner part B to use
run motifscanner part B with fasta file, tmp.DS.faa, motif meme file, background, outputfile
'''

from sys import argv
import random as rd
import os
import sys



def main():
    #parsing argument
    faafile = None
    memefile = None
    outfile = None
    backgroundfile = None
    backgroundnotprovided = True
    seed = None
    pvalue = None
    typeEF = False
    background = {}

    curdir=os.getcwd()
    print "Current directory is",curdir
    #change working directory to the script's 
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)

    
    # get command line arguments
    #
    
    usage = """USAGE: 
    python %s [options] -f fastafile -m motiffile -o output file

        -f <filename>       input file, FASTA format 
        -m <filename>       input motif file, MEME format
        --typeEF            signal that this is type EF       
        -s <number>         seed for the dinucleotide shuffling, if not specified, random seed will be assigned (optional)
        -p <number>         P-value cutoff for calling matches, default is 0.0001 (optional)
        -b <filename>       background distribution of the bases, default is equal probability of each base (optional)
        -o <filename>       name of the output scan file to be created (optional)
        -h                  print this usage message
    """ % (sys.argv[0])

    # parse command line
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if (arg == "-f"):
            i += 1
            try: faafile = sys.argv[i]
            except:
            	print " Error in -f" 
            	sys.exit(1)
        elif (arg == "-m"):
            i += 1
            try: memefile = sys.argv[i]
            except:
            	print "Error in -m" 
            	sys.exit(1)
        elif (arg == "--typeEF"):
            print "RUNNING TYPE EF ..."
            #i += 1
            typeEF = True
        elif (arg == "-b"):
            i += 1
            try: 
                backgroundfile = sys.argv[i]
                backgroundnotprovided = False
            except:
            	print "Error in -b" 
            	sys.exit(1)
        elif (arg == "-p"):
            i += 1
            try: pvalue = float(sys.argv[i])
            except:
            	print "Error in -p" 
            	sys.exit(1)
        elif (arg == "-s"):
            i += 1
            try: seed = int(sys.argv[i])
            except:
            	print "Error in -s, please check the command", sys.argv
            	sys.exit(1)
        elif (arg == "-o"):
            i += 1
            try: outfile = sys.argv[i]
            except:
            	print "Error in -o" 
            	sys.exit(1)
        elif (arg == "-h"):
        	print >> sys.stderr, usage; sys.exit(1)
        else:
            print >> sys.stderr, "Unknown command line argument: " + arg + " \nPlease check the usage with -h"
            sys.exit(1)
        i += 1
    if len(argv)==1:
    	print "No input, exit."
    	print >> sys.stderr, usage; sys.exit(1)
    	sys.exit(1)
    # check that required arguments given
    if (faafile == None):
    	print "No input fasta file, exit."
    	print >> sys.stderr, usage; sys.exit(1)
    	sys.exit(1)
    elif (memefile == None):
    	print "No input motif file, exit."
        print >> sys.stderr, usage; sys.exit(1)
    	sys.exit(1)
    elif (seed == None):
    	print "No seed specified, random seed is chosen."
        seed=rd.randint(0, 10000)
    if (pvalue == None):
    	print "No pvalue cutoof specified, using default 0.0001"
        pvalue=0.0001
    if (outfile == None):
        outfile = faafile+".scanned.txt"
    	print "No output file specified, will use the default",outfile
    if (backgroundfile == None):

    	if not (typeEF):
            print "No background file specified, will use the default: equal probabilities for A, C, G, T, E"
            background={'A':0.2,'C':0.2,'G':0.2,'T':0.2,'E':0.2}
        else:
            print "No background file specified, will use the default: equal probabilities for A, C, G, T, E, F"
            background={'A':0.16666666666,'C':0.16666666666,'G':0.16666666666,'T':0.16666666666,'E':0.16666666666, 'F':0.16666666666 }
        backgroundfile="./baseComposition.tmp.txt"
        outbackground=open(backgroundfile,'w')
        print backgroundfile
        for char in background:
            outbackground.write(char+'\t'+str(background[char])+'\n')
        outbackground.close()




    
    
    #loading data:
    print "Reading fasta file"
    faaname=faafile
    '''if os.path.isfile(faafile): #if file in current dir
        faaname=faafile
    elif os.path.isfile(curdir+'/'+faafile):
        faaname=curdir+'/'+faafile
    else:
        print "File does not exist:",faaname'''

    memename=memefile
    '''if os.path.isfile(memefile): #if file in current dir
        memename=memefile
    elif os.path.isfile(curdir+'/'+memefile):
        memename=curdir+'/'+memefile
    else:
        print "File does not exist:",memename'''
    outname=outfile


    faafile=open(faaname).read().strip().split(">")
    faafile=faafile[1:]
    faaseqs={}
    totalbasenum=0
    for line in faafile:
        tmp=line.strip().split("\n")
        faaseqs[tmp[0]]=tmp[1]
        totalbasenum+=len(tmp[1])
    print "Number of sequences:",len(faaseqs),". Number of total bases:",totalbasenum
    
    estimatekmernum=totalbasenum-len(faaseqs)*10 #assuming kmer length is 10
    timestoshuffle=1
    if estimatekmernum<1000000:
        timestoshuffle=1000000/estimatekmernum+1
        print "Less than 10^6 kmers, run Di-nuc shuffling",timestoshuffle,"times"
    else:
    	#print "More than 1000000"
    	print "Di-nuc shuffling 1 times"
    #run di-nuc shuffling
    shufflefile=faaname+".tmp.DS.faa"
    command=''
    if typeEF:
        command= "python fasta-dinucleotide-shuffle_typeEF.py -s "+str(seed)+" -f "+faaname+" -c "+str(timestoshuffle)+" > "+shufflefile
    else:
        command= "python fasta-dinucleotide-shuffle_typeE.py -s "+str(seed)+" -f "+faaname+" -c "+str(timestoshuffle)+" > "+shufflefile
    os.system(command)
	#run motif scanning
	#usage: julia motifscanner_B.jl p-valuecutoff memefile positivefile DSnegativefile output
    print "Scanning motifs"
    if typeEF:
        command= "julia motifscannerB_typeEF.jl "+str(pvalue)+" "+memename+" " +faaname+" "+shufflefile+" "+outname+" "+backgroundfile
    else:
        command= "julia motifscannerB_typeE.jl "+str(pvalue)+" "+memename+" " +faaname+" "+shufflefile+" "+outname+" "+backgroundfile
    print command
    os.system(command)	
    os.system("rm "+shufflefile)
    if backgroundnotprovided:
        os.system("rm "+backgroundfile)

if __name__ == "__main__":
    main()
