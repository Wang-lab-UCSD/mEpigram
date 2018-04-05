import argparse
import os 
from collections import Counter
import sys

def revcomplE(word):
	revComp=''
	transdict={'A':'T','C':'G','G':'C','T':'A'}
	index=0
	while index <len(word):
		if word[index] == 'E' :#if is methylated C
			if index <len(word)-1 and word[index+1]=='G': #if G is next to it
				revComp += 'GE'
				index+=2   
			else:
				revComp += 'G'
				index+=1
		else:
			try: 
				revComp += transdict[word[index]]
			except KeyError:
				revComp += 'N'
			index+=1		 
	return revComp[::-1]

def revcomplEF(word):
	revComp=''
	transdict={'A':'T','C':'G','G':'C','T':'A','E':'F','F':'E'}
	index=0
	while index <len(word):
		try:
			revComp += transdict[word[index]]
		except KeyError:
			revComp += 'N'
		index += 1		 
	return revComp[::-1]

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-gd", "--genomedir", help="a directory that contains genome fasta sequence files, ideally one chromosome per file")
	parser.add_argument("-k", "--klength",type =int, help="integer, the length of k-mers")
	parser.add_argument("-m", "--mode",type =str, help="specify whether using typeE or typeEF")
	#parser.add_argument('--typeEF', dest='typeEF',action='store_true', help="using typeEF, default is typeE")
	#parser.set_defaults(typeEF=False)

	args = parser.parse_args()

	genomedir = args.genomedir
	if genomedir == None:
		print "GENOMEDIR is missing, please specify"
		return
	try:
		k = int(args.klength)
	except:
		print "KLENGTH must be specified"
		return

	mode = args.mode
	print mode

	typeEF = False
	if mode == None:
		print "Please choose a mode, typeE or typeEF"
		sys.exit(1)
	elif mode != "typeE" and mode != "typeEF":
		print "ERROR: unrecognized mode", mode
		sys.exit(1)
	else:
		if mode == typeEF:
			typeEF = True

	print "GENOMEDIR is",genomedir
	print "klength is", k
	print "Using ",mode



	kmercounts = {} # 
	total = 0
	for file in os.listdir(genomedir):
		print "Reading ...",file
		for line in open(genomedir+'/'+file):
			if line[0] == ">": # this is the name, skip
				continue
			seq = line.strip().upper()
			print len(seq)
			
			#count the kmers
			for i in range(len(seq) - k):
				kmer = seq[i:i+k]
				if "N" in kmer:
					continue
				try: 
					kmercounts[kmer] += 1
				except KeyError:
					kmercounts[kmer] = 1
				total += 1
			
			if typeEF: #check if type E or EF
				seq = revcomplEF(seq)
			else:
				seq = revcomplE(seq)
			
			for i in range(len(seq) - k):
				kmer = seq[i:i+k]
				if "N" in kmer:
					continue
				try: 
					kmercounts[kmer] += 1
				except KeyError:
					kmercounts[kmer] = 1
				total += 1


	#output
	outfile = ""
	if typeEF:
		outfile = "background_typeEF-"+str(k) +".tsv"
	else:
		outfile = "background_typeE-"+str(k) +".tsv"
	target = open(outfile,'w')
	target.write("TOTAL\t"+str(total)+'\n')
	for kmer in kmercounts:
		target.write(kmer+'\t'+str(kmercounts[kmer])+'\n')
	target.close()
		


	return 


if __name__ == "__main__":
	main()