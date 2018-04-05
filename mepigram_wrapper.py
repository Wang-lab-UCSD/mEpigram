#!/usr/bin/env python

"""
TODO: 


Workflow:

	- Takes a set of sequences
	- Di-nuc shuffles it .
	- Takes the background model
	- Takes the number of of maximum motifs

"""
from sys import argv
import random as rd
import os
import sys
import os.path
import argparse


def load_motifs_typeE(filename):
	file=open(filename)
	seq=file.read().split("MOTIF")
	seq=seq[1:]
	motifs={}
	infos={}
	for s in range(len(seq)):
		t=seq[s].strip().split("\n")
		name=int(t[0].split('_')[0])
		motifs[name]=t[2:]
		infos[name]=t[0:2]
	for m in motifs:
		tdict={'A':[],'C':[],'G':[],'T':[],'E':[],}
		for pos in range(len(motifs[m])):
			tmp=motifs[m][pos].strip().split("\t")
			tdict['A']+=[float(tmp[0])]
			tdict['C']+=[float(tmp[1])]
			tdict['G']+=[float(tmp[2])]
			tdict['T']+=[float(tmp[3])]
			tdict['E']+=[float(tmp[4])]
		motifs[m]=tdict
	
	return motifs,infos

def load_motifs_typeEF(filename):
	file=open(filename)
	seq=file.read().split("MOTIF")
	seq=seq[1:]
	motifs={}
	infos={}
	for s in range(len(seq)):
		t=seq[s].strip().split("\n")
		name=int(t[0].split('_')[0])
		motifs[name]=t[2:]
		infos[name]=t[0:2]
	for m in motifs:
		tdict={'A':[],'C':[],'G':[],'T':[],'E':[],'F':[]}
		for pos in range(len(motifs[m])):
			tmp=motifs[m][pos].strip().split("\t")
			tdict['A']+=[float(tmp[0])]
			tdict['C']+=[float(tmp[1])]
			tdict['G']+=[float(tmp[2])]
			tdict['T']+=[float(tmp[3])]
			tdict['E']+=[float(tmp[4])]
			tdict['F']+=[float(tmp[5])]
		motifs[m]=tdict
	
	return motifs,infos

def calcBaseComp(backgroundfile,alphabet,basecompfile):
	#backgroundfile=argv[1]
	#alphabetfile = argv[2] # alphabet, a squence of alphabets that are used, separated by commas. Example: A,G,C,T
	output = basecompfile


	seqs=open(backgroundfile).read().split('\n')
	total=float(seqs[0].strip().split('\t')[1])
	seqs=seqs[1:]
	#seqs[:10]
	totalalphacounts={}
	for a in alphabet:
		totalalphacounts[a] = 0
	kmerlen=len(seqs[0].strip().split('\t')[0])
	for line in seqs:
		tmp=line.strip().split("\t")
		if len(tmp)<2:
			#print tmp
			continue
		kmercount=int(tmp[1])
		alphacounts={}
		for char in tmp[0]:
			if char not in alphacounts:
				alphacounts[char]=1
			else:
				alphacounts[char]+=1
		for char in alphacounts:
			if char not in totalalphacounts:
				totalalphacounts[char]=alphacounts[char]*kmercount
			else:
				totalalphacounts[char]+=alphacounts[char]*kmercount

	outfile=open(output,'w')
	for char in totalalphacounts:
		line=char+'\t'+str(float(totalalphacounts[char])/total/kmerlen)
		outfile.write(line+'\n')
	return

def taggingmotifs(filename,outfile):
	print "tagging m-motifs in",filename
	'''This function adds a tag to motifs with P(E) >= 0.5'''
	#filename="./test_mepigram_pipeline_complete.meme"
	#outfile=filename.replace(".meme",".tagged.meme")
	motifs,infos=load_motifs_typeE(filename)
	header='''MEME version 4.5 - modififed
ALPHABET= ACGTE 
strands: + 
Background letter frequencies
A 0.295 C 0.205 G 0.205 T 0.295 E 0.0076

'''
	#print infos
	#print the motifs with the E or non-E or EF motif tag
	target=open(outfile,'w')
	target.write(header)
	threshold=0.5
	#print header
	for m in motifs:
		modified=False
		for pos in motifs[m]['E']:
			if pos > threshold:
				modified=True
				break
		firstline=''
		if modified==True:
			firstline="MOTIF"+'\t'+infos[m][0].strip()+"_m-motif"
		else:
			#print m,"motif"
			firstline="MOTIF"+'\t'+infos[m][0].strip()
		target.write(firstline+'\n')
		target.write(infos[m][1]+'\n')
		for i in range(len(motifs[m]['A'])):
			line=[]
			for j in ['A','C','G','T','E']:
				line+=[motifs[m][j][i]]
			line='\t'.join([str(x) for x in line])
			target.write(line+'\n')
	target.close()
	print "Finished tagging"
	return

def main():
	#parsing argument
	faafile = None
	#memefile = None
	outfile = None #this contains the meme file, the enrichment file .
	mode = None # typeEF or typeE
	backgroundfile = None
	graphdir = None 
	filter_boolean = False # Not used right now
	maxmotifnum = None
	enrichmentmode = "none"  #choose whether to calculate the enrichment of all motifs found or just the m-motifs, or just non-m motifs, or none
	#enrichmentmode is not used right now, default is all

	####Testing using hardcodes
	maxmotifnum=200
	mode=None
	seed=rd.randint(0, 10000000)


	# parse command line
	# Required arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("faafile", help = "input file, FASTA format")
	parser.add_argument("backgroundfile", help = "the background model, which contains k-mer counts")
	parser.add_argument("graphdir", help = "the directory that contains the graph")
	parser.add_argument("mode", help = "typeE|typeEF: specify whether using typeE or typeEF")
	#Optional arguments
	parser.add_argument("-o", "--output", default = None, help = "output directory to be created")
	parser.add_argument("-n", "--maxmotifnum",default = 200, type=int, help="integer, maximum number of motifs to find; default is 200")
	parser.add_argument("-ml", "--makelogo", default = "n", help = "y|n: whether to create motif logos; default is No")

	args = parser.parse_args()
	faafile = args.faafile
	backgroundfile = args.backgroundfile
	graphdir = args.graphdir 
	mode = args.mode # typeEF or typeE

	outfile = args.output #this contains the motif file, the enrichment file .
	maxmotifnum = args.maxmotifnum
	makelogo = False
	if args.makelogo != "y" and args.makelogo != "n":
		print "ERROR: parameter not recognized in --makelogo:", args.makelogo
		sys.exit(1)
		if args.makelogo == "y":
			makelogo = True
	if mode != "typeE" and mode != "typeEF":
		print "ERROR: paramether not recognized in --mode:", args.mode

	#check if files are there
	if os.path.isfile(faafile) == False: 
		print "ERROR: FASTA input file doesn't exist at",faafile
		sys.exit(1)
	if os.path.isfile(backgroundfile) == False:
		print "ERROR: Background model file doesn't exist at", backgroundfile
		sys.exit(1)
	if os.path.isdir(graphdir) == False:
		print "ERROR: graph dir doesn't exist at",graphdir
		sys.exist(1)
	if (outfile == None):
		outfile = faafile+".mepigram"
		print "No output file specified, will use the default",outfile


	# Now load the data:
	# shuffle the data
	print "Reading fasta file"
	faaname=faafile
	faafile=open(faafile).read().strip().split(">")
	faafile=faafile[1:]
	faaseqs={}
	totalbasenum=0
	#remove lines with less than a certain length 
	badcount=0
	minlen=20
	for line in faafile:
		tmp=line.strip().split("\n")
		seq=''.join(tmp[1:])
		if len(seq)>=minlen:
			faaseqs[tmp[0]]=seq
			totalbasenum+=len(seq)
		else:
			badcount+=1
	if badcount>0:
		print badcount,"sequences shorter than "+str(minlen)+" base pairs. They are skipped."
	print "Number of sequences:",len(faaseqs),". Number of total bases:",totalbasenum
	
	estimatekmernum=totalbasenum-len(faaseqs)*10 #assuming kmer length is 10
	timestoshuffle=1
	minkmernumber=1 # Maybe use 20 millions here later, testing this out
	if estimatekmernum<minkmernumber:
		timestoshuffle=minkmernumber/estimatekmernum+1
		print "Less than "+str(minkmernumber)+" kmers, run Di-nuc shuffling",timestoshuffle,"times"
	else:
		#print "More than 1000000"
		print "Di-nuc shuffling 1 times"
	
	#run di-nuc shuffling
	shufflefile=faaname+'.'+str(seed)+".DS.tmp.faa"
	command=''
	if mode == 'typeE':
		command="python fasta-dinucleotide-shuffle_typeE.py -s "+str(seed)+" -f "+faaname+" -c "+str(timestoshuffle)+" > "+shufflefile
	else:
		command="python fasta-dinucleotide-shuffle_typeEF.py -s "+str(seed)+" -f "+faaname+" -c "+str(timestoshuffle)+" > "+shufflefile

	print command
	os.system(command)

	print "Running mEpigram... "
	if mode == 'typeE':
		command="python mepigram_typeE.py "+faaname+" "+shufflefile+" "+backgroundfile+" "+graphdir +" "+outfile+" "+str(maxmotifnum)
	else:
		command="python mepigram_typeEF.py "+faaname+" "+shufflefile+" "+backgroundfile+" "+graphdir +" "+outfile+" "+str(maxmotifnum)
	print command
	os.system(command)

	'''This part renames the motifs into meth and unmeth motifs
    it is not used in this version'''
	memefile=outfile#+'.tagged'
	#taggingmotifs(outfile,memefile)
	#memefile="../ENCFF002CQR.2motifs.meme"

	resultdir='/'.join(memefile.split('/')[:-1])
	resultdir=memefile+".results"
	os.system("rm -r "+resultdir)
	os.system("mkdir "+resultdir)
	

	"""
	Filtering step will go in here
	Combine all of the motifs together into a list, then calculate pairwise similarity distance of all of them. 
	If two motifs have similarity higher than a certain threshold, keep one, discard one. 
	"""
	#print "Filtering results..."


	"""
	Calculate enrichment of motifs: 2 methods
	
		* Fisher p-value: <Use this one for now>
			- Scan through to get the max score of each region
			- Determine a cutoff such that the p-value is maximized
		* Straight enrichment (Vu's method)
			- Get a score distribution of k-mers for each motif by scanning through randomly shuffled regions
			- Determine a cutoff accordingly to a p-value
			- Use the cutoff to call for matches.
	"""
	print "Calculating enrichment by scanning..."
	# use method 1,

	print "Estimating background base compostion..."

	baseCompositionFile=backgroundfile+".basecomposition.tmp"
	alphabet = []
	if mode == "typeEF":
		alphabet = ['A','C','G','T','E','F']
	elif mode == "typeE":
		alphabet = ['A','C','G','T','E']
	calcBaseComp(backgroundfile,alphabet,baseCompositionFile)
	
	"""
	command="python baseComposition.py "+backgroundfile+" "+baseCompositionFile
	print command
	os.system(command)
	"""

	
	if mode == 'typeE':
		print "Scanning on positive sequences..."
		command="julia quickPssmScanBestMatchLiteTypeE.jl "+memefile+" "+faaname+" "+"quickscan.positive.tmp"+" "+resultdir+" "+baseCompositionFile
		print command
		os.system(command)	
		print "Scanning on negative sequences..."
		command="julia quickPssmScanBestMatchLiteTypeE.jl "+memefile+" " +shufflefile+" "+"quickscan.negative.tmp"+" "+resultdir+" "+baseCompositionFile
		print command
		os.system(command)

	else:
		print "Scanning on positive sequences..."
		command="julia quickPssmScanBestMatchLiteTypeEF.jl "+memefile+" "+faaname+" "+"quickscan.positive.tmp"+" "+resultdir+" "+baseCompositionFile
		print command
		os.system(command)  
		print "Scanning on negative sequences..."
		command="julia quickPssmScanBestMatchLiteTypeEF.jl "+memefile+" " +shufflefile+" "+"quickscan.negative.tmp"+" "+resultdir+" "+baseCompositionFile
		print command
		os.system(command)

	#sdos.system("rm "+shufflefile)
	print "Calculating Fisher P-values..."
	#get the result files together
	scannedresults={}
	for file in os.listdir(resultdir):
		if "quickscan.positive" in file:
			motif=file.split('.quickscan')[0]
			if motif not in scannedresults:
				scannedresults[motif]={}
				scannedresults[motif]['P']=file
			else:
				if "P" in scannedresults[motif]:
					"ERROR!!!: Already found positive scan file for motif",motif
				else:
					scannedresults[motif]['P']=file
			
		elif "quickscan.negative" in file:
			motif=file.split('.quickscan')[0]
			if motif not in scannedresults:
				scannedresults[motif]={}
				scannedresults[motif]['N']=file
			else:
				if "N" in scannedresults[motif]:
					print "ERROR!!!: Already found negative scan file for motif",motif
				else:
					scannedresults[motif]['N']=file
		else:
			continue
	print "Found scanned results for",len(scannedresults),"motifs"

	fisherresults=[]
	for motif in scannedresults:
		print motif
		if len(scannedresults[motif])!=2:
			print "ERROR!!!: number of scanned files for motif",motif,"is not 2"
		command="python2.7 fisher_P-value.py "+resultdir+'/'+scannedresults[motif]['P']+" "+resultdir+'/'+scannedresults[motif]['N']+" "+resultdir+'/'+motif+".fisher.tmp" +" "+motif
		#print command
		os.system(command)
		fisherresults+=[open(resultdir+'/'+motif+".fisher.tmp").read().strip()]

	print "Calculating enrichments..."
	out=open(resultdir+'/enrichments.tsv','w')
	out.write("MOTIF\tp-value\tscoreCutoff\tPosMatches\tPosNonMatches\tNegMatches\tNegNonMatches\tEnrichment"+'\n')
	for line in fisherresults:
		out.write(line+'\n')
	out.close()
	#concat the fisher results

	if makelogo:
		logodir = memefile + '.LOGOS'
		if mode == "typeE":
			cmd = "python makeLOGO.py -m %s -o %s" %(memefile, logodir)
		else:
			cmd = "python makeLOGO.py --typeEF -m %s -o %s" %(memefile, logodir)

		print "Making motif LOGOs..."
		os.system(cmd)
		os.system("mv %s %s"%(logodir, resultdir))


	os.system("mv %s %s" %(memefile, resultdir+"/motifs.mepigram.meme"))
	os.system("mv %s %s" %(outfile, resultdir))



	#print "Cleaning up temporary files..."
	os.system("rm "+shufflefile)
	#ftoremove=[]
	for file in os.listdir(resultdir):
		if "tmp" in file:
			command="rm "+resultdir+'/'+file
			#print command
			os.system(command)

	
if __name__ == "__main__":
	main()



