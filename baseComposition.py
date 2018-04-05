
#!/usr/bin/env python

"""
Calculate the base composition from the background file

"""

from sys import argv

def main():
	backgroundfile=argv[1]
	alphabetfile = argv[2] # alphabet, a squence of alphabets that are used, separated by commas. Example: A,G,C,T
	output = argv[3]


	seqs=open(backgroundfile).read().split('\n')
	total=float(seqs[0].strip().split('\t')[1])
	seqs=seqs[1:]
	#seqs[:10]
	totalalphacounts={}
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


if __name__=="__main__":
	main()