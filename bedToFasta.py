#!/usr/bin/env python
import time
import sys
from sys import argv
import os

'''Usage:
    python BedtoFasta.py genome bedfile output

    The genome is a single FASTA file
    '''
def main():
    """metgenome=argv[1]
    bedfile=argv[2]
    output=argv[3]"""

    metgenome = None
    bedfile = None
    output = None
    ref_is_file = True #Default input reference genome is a FASTA file

    # get command line arguments
    #
    usage = """USAGE:
    %s [options]

        -f <filename>       input file, in BED format
        -r <ref_file>       the reference genome (default to be in the form of a single FASTA file)
        -d <ref_directory>  indicates that the reference genome in the form of a directory containing files for chromosome in FASTA format (optional)
        -o <output file>    name of the output FASTA file to be created
        -h                  print this usage message
    """ % (sys.argv[0])

    # parse command line
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if (arg == "-f"):
            i += 1
            try: bedfile = sys.argv[i]
            except:
                print " Error in -f"
                sys.exit(1)
        elif (arg == "-d"):
            ref_is_file = False
        elif (arg == "-r"):
            i += 1
            try: metgenome = sys.argv[i]
            except:
                print "Error in -r, please check the command", sys.argv
                sys.exit(1)
        elif (arg == "-o"):
            i += 1
            try: output = sys.argv[i]
            except:
                print "Error in -o"
                sys.exit(1)
        elif (arg == "-h"):
            print >> sys.stderr, usage; sys.exit(1)
        else:
            print >> sys.stderr, "Unknown command line argument: " + arg + " \nPlease check the usage with -h"
            sys.exit(1)
        i += 1
    # check that required arguments given
    if (bedfile == None):
        print "No input, exit."
        #sys.exit(1)
        print >> sys.stderr, usage; sys.exit(1)
    elif (output == None):
        print "No output, exit."
        print >> sys.stderr, usage; sys.exit(1)
    elif (metgenome== None):
        print "No reference genome, exit. See usage with -h"
        sys.exit(1)

    print "Working with",bedfile, metgenome, output

    #load the peaks bed file
    #example:
    print "Reading BED file..."
    file=open(bedfile)
    met_reads=[]
    for line in file:
        tmp=line.strip().split("\t")
        chrom=tmp[0]
        start=int(tmp[1])
        end=int(tmp[2])
        met_reads+=[[chrom,start,end]]
    print "Number of regions",len(met_reads)
    print "First region in BED file",met_reads[0]
    #sys.exit(main(0))

    #load metgenome
    print "Reading the reference genome..."
    #check if the ref is in the file or directory format
    genome={}
    if ref_is_file:
    	#example: file=open("/Users/vungo/Downloads/Work/genomes/mm10.fa")
    	file = open(metgenome)
    	seq = file.read()
    	seq = seq.split(">")
    	seq = seq[1:]
    	#print len(seq)
    	for chro in range(len(seq)):
        	temp = seq[chro].strip().split('\n')
        	chromosome = temp[0] #name of the chromosome
        	print chromosome
        	dna = "".join(temp[1:])
        	genome[chromosome] = dna
    	print genome.keys()
    else:
    	print "Reference is a directory. Reading the files..."
    	chromfiles=os.listdir(metgenome)
    	#check if the names of the chr files are chrXX or chrXX.faa, convert the latter to the former
    	chromnames=[]
    	for filename in chromfiles:
    		print filename
    		chro = filename.strip().split(".")[0]
    		chromnames += [chro]
    		seq = open(metgenome+"/"+filename).read().strip().split("\n")
    		seq = seq[1:] #the first element is the chromosome name, discard it
    		genome[chro] = "".join(seq)
    	print genome.keys()

    print "Outputing the FASTA file"
    target=open(output,"w")
    #DONT CARE ABOUT DIRECTION! because of the methylation ... its complicated
    num=0
    time_start=time.time()
    for row in range(len(met_reads)):
        if num%1000==0:
            print "number of regions converted",num
        num+=1
        chrom=met_reads[row][0]
        start=met_reads[row][1]
        end=met_reads[row][2]
        #direction=klf4.loc[row,0]
        name=">%s%s%s\n" %(chrom,"_",str(start))
        #name=">"+chrom+str(start)+"\n"
        seq=genome[chrom][start:end]
        target.write(name)
        target.write(seq+"\n")
    print "Running time is",time.time()-time_start
    target.close()

if __name__=="__main__":
    main()
