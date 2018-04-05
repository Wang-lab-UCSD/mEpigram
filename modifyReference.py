#!/usr/bin/env python

"""
usage: Takes a wig file, a genome_dir , a genome type, it makes a genome_met.faa directory

TODO:
-add a parameter to choose whether the G should be converted or not. 
-or, make the user deal with it before inserting E into the genome.

Check input file format.

"""
from sys import argv
import os
import sys
def readFASTA(faafile):
    lines=open(faafile).read().split()
    seq=[]
    for line in lines:
        if ">" not in line:
            seq+=[line.strip()]
    seq=''.join(seq)

    return seq

def insertEtochrom(chrom,hashdata,chromseq,threshold):
    #use this because it's faster
    #use loc -1 because the coordinates are 1-based while here things are 0-based
    if chrom not in hashdata:
        print chrom,"not in methylation data"
        return chromseq
    print "converting chromosome",chrom,"to list"
    chromseq=list(chromseq)
    print 'done converting string to list'
    num=0
    nummet=0
    numerr=0
    for loc in hashdata[chrom]:
        if num%100000==0:print num,nummet,numerr
        if hashdata[chrom][loc]>=threshold and (chromseq[loc-1].upper()=="C" or chromseq[loc-1].upper()=="E"):
        	chromseq[loc-1]='E'
        	nummet+=1
        else:
        	numerr+=1
        num+=1
    print "done inserting E to list.",nummet,"positions converted to E.",numerr,"positions were not converted."
    chromreturn=''.join(chromseq)
    return chromreturn

def insertEFtochrom(chrom,hashdata,chromseq,threshold):
    #use this because it's faster
    #use loc -1 because the coordinates are 1-based while here things are 0-based
    if chrom not in hashdata:
        print chrom,"not in methylation data"
        return chromseq
    print "converting chromosome",chrom,"to list"
    chromseq=list(chromseq)
    print 'done converting string to list'
    num=0
    nummetE=0
    nummetF=0
    numerr=0
    for loc in hashdata[chrom]:
        if num%100000==0:print num,nummetE,nummetF,numerr
        if hashdata[chrom][loc]>=threshold and (chromseq[loc-1].upper()=="C" or chromseq[loc-1].upper()=="E"):
            chromseq[loc-1]='E'
            nummetE+=1
        elif hashdata[chrom][loc]>=threshold and (chromseq[loc-1].upper()=="G" or chromseq[loc-1].upper()=="F"):
            chromseq[loc-1]='F'
            nummetF+=1
        else:
            numerr+=1
        num+=1
    print "done inserting E,F to list.",nummetE,nummetF,"positions converted to E,F.",numerr,"positions were not converted."
    chromreturn=''.join(chromseq)
    return chromreturn


#main 


def main():
    #parsing argument
    wigfile = None
    outdir = None
    genome_dir = None
    threshold = 0.5
    isWig = False
    typeEF = False
    
    # get command line arguments
    #
    usage = """USAGE: 
    %s [options]

        -f <filename>  		input file, default is wig format
        --wig                   indicates that the input file is in wig format
        --typeEF                convert the modified base to E and the complementary base to F, default is False (optional)
        -t <number>             methylation threshold above which C will be considered E, default is 0.5 (optional)       
        -r <ref_genome> 	directory of the reference genome (best put in a separate directory)
        -o <output dir> 	name of the output dir to be created (eg. hg19_methylated)
        -h              	print this usage message
    """ % (sys.argv[0])

    # parse command line
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if (arg == "-f"):
            i += 1
            try: wigfile = sys.argv[i]
            except:
            	print " Error in -f" 
            	sys.exit(1)
        elif (arg == "-t"):
            i += 1
            try: threshold = float(sys.argv[i])
            except:
            	print "Error in -t" 
            	sys.exit(1)
        elif (arg == "--wig"):
            #i += 1
            #print "Is Bed"
            try: isWig = True
            except:
                print "Error in --bed" 
                sys.exit(1)
        elif (arg == "--typeEF"):
            #i += 1
            #print "Is Bed"
            try: typeEF = True
            except:
                print "Error in --typeEF" 
                sys.exit(1)
        elif (arg == "-r"):
            i += 1
            try: genome_dir = sys.argv[i]
            except:
            	print "Error in -r, please check the command", sys.argv
            	sys.exit(1)
        elif (arg == "-o"):
            i += 1
            try: outdir = sys.argv[i]
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
    if (wigfile == None):
    	print "No input, exit."
    	print >> sys.stderr, usage; sys.exit(1)
    	sys.exit(1)
    elif (outdir == None):
    	print "No output, exit."
    	sys.exit(1)
    elif (genome_dir == None):
    	print "No reference genome selected, exit."
    	sys.exit(1)

    
    print "Working with",wigfile,genome_dir,outdir

    

    metdata={}
    if isWig == True:
        print "Loading the wig input file..."
        seq=open(wigfile).read()
        seq=seq.strip().split('chrom')
        seq=seq[1:]
        for i in range(len(seq)):
            seq[i]=seq[i].strip().split('\n')
        for i in range(len(seq)):
            chrom=seq[i][0][1:]
            print chrom
            metdata[chrom]={}
            print len(seq[i]),"positions"
            for j in range(1,len(seq[i])):
                tmp=seq[i][j].strip().split("\t")
                try:
                    coor=int(tmp[0])
                    beta=float(tmp[1])
                    metdata[chrom][coor]=beta
                except:
                    #print seq[i][j] #this is at the end of each chromosome
                    continue
    if isWig == False: #this is 0-based, while the metdata assumes 1-based coordinates
        bedfile = wigfile
        print "Loading the BEDGRAPH input file..."
        seqs=open(bedfile).read().strip().split('\n')
        for line in seqs:
            tmp=line.split('\t')
            chrom=tmp[0]
            loc=int(tmp[1])
            value=float(tmp[3])
            try:
                metdata[chrom][loc+1]=value
            except KeyError:
                metdata[chrom]={}
                metdata[chrom][loc+1]=value
        for chrom in metdata:
            print chrom, len(metdata[chrom]),'positions'
    #print metdata

    #load the genome
    print "Loading the reference genome..."
    genome={}
    filelist=os.listdir(genome_dir)
    for f in filelist:
        #if "chr" == f[:3]:
        seq=readFASTA(genome_dir+'/'+f)
        chrom = f.strip()
        genome[chrom]=seq
        print "loaded",chrom
    chromnames=sorted(genome.keys())
    #print chromnames
    #inserting E to each chromosome
    if typeEF:
        for chrom in chromnames:
            if chrom not in metdata:
                print "No methylation information for chromosome",chrom
                #continue
            else:
                temp = insertEFtochrom(chrom,metdata,genome[chrom],threshold)
                genome[chrom]=temp
    else:
        for chrom in chromnames:
            if chrom not in metdata:
                print "No methylation information for chromosome",chrom
                #continue
            else:
                temp = insertEtochrom(chrom,metdata,genome[chrom],threshold)
                genome[chrom]=temp
        #print genome
	
	#error checking step
    num=0
    numerr=0
    for chrom in metdata.keys():
        if chrom not in genome:
            print chrom, "not in the reference genome. Skipping..."
            continue
        else:
            for loc in metdata[chrom]:
                #if genome[chrom][loc-1].upper()=='E' and genome[chrom][loc-1].upper()=='G' and metdata[chrom][loc]<0.5: 
                if genome[chrom][loc-1]=='E' or genome[chrom][loc-1]=="F" and  metdata[chrom][loc]<threshold:
                    numerr+=1
                    #print "error", chrom, loc, genome[chrom][loc-1]
                num+=1
                if num%1000000==0: print num,numerr
                #if num>10000000:break

    #print the file:
    #mkdir 
    print "Outputing to",outdir
    os.system("mkdir "+outdir )


    for chrom in genome:
    	print chrom
        target=open(outdir+"/"+chrom,"w")
        target.write(">"+chrom+"\n")
        target.write(genome[chrom])
        target.close()


if __name__ == "__main__":
    main()
