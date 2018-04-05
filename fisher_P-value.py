#!/usr/bin/env python


"""
This determines a score cutoff based on the scan results of a motif and calculates the fisher-pvalue at that cutoff
Usage: python fisher_p-value.py positivescannedfile negativescannedfile outputfile
"""
from sys import argv
from scipy.stats import stats
import numpy as np
def loadscanresult(filename):
	"""
	Each file contains the sequence names and scores of a motif scan result
	Example: chr1_1910\t1.01013
				.....
	"""
	file=open(filename)
	seqnames=[]
	scores=[]
	for line in file:
		tmp=line.strip().split('\t')
		try:
			sn=tmp[0]
			s=float(tmp[1])
			seqnames+=[sn]
			scores+=[s]
		except:
			pass

	return seqnames,scores

def FisherScoreCutoff(posScores,negScores):
	"""Define the upper and lower 5 percent bounds of the distribution of scores, 
		then increase the score cutoff by 1/100 the difference between the two scores at a time
		compute the Fishter p-value for each score cutoff and record the best one """
	#posScores=sorted(posScores,reverse=True)
	negScores=sorted(negScores,reverse=True)
	
	contingencyTable=[[0,0],[0,0]]
	scoreCutoff=negScores[len(negScores)/100]
	posCount=0
	negCount=0
	for i in range(len(posScores)):
		if posScores[i]>=scoreCutoff:
			posCount+=1
	for j in range(len(negScores)):
		if negScores[j]>=scoreCutoff:
			negCount+=1
	contingencyTable[0][0]=posCount
	contingencyTable[0][1]=len(posScores)-posCount
	contingencyTable[1][0]=negCount
	contingencyTable[1][1]=len(negScores)-negCount
	ob,p_value=stats.fisher_exact(contingencyTable)
	return p_value,scoreCutoff,contingencyTable[0]+contingencyTable[1]
	'''upperBound=negScores[len(negScores)/10000]
	lowerBound=negScores[len(negScores)/100*50]
	print "Lower index",len(negScores)/100*50
	print len(posScores),len(negScores)
	print "UPPER LOWER",upperBound,lowerBound
	intervalnum=100
	interval=(upperBound-lowerBound)/intervalnum
	print "INTERVAL",interval
	bestcuttoff=-100.0
	bestpvalue=1.0
	contingencyTable=[[0,0],[0,0]]
	besttable=contingencyTable
	for i in range(intervalnum):
		scoreCutoff=lowerBound+interval*i
		if scoreCutoff >upperBound:
			break
		posCount=0
		negCount=0
		for j in range(len(posScores)):
			if posScores[j]>=scoreCutoff:
				posCount+=1
		for j in range(len(negScores)):
			if negScores[j]>=scoreCutoff:
				negCount+=1
		contingencyTable[0][0]=posCount
		contingencyTable[0][1]=len(posScores)-posCount
		contingencyTable[1][0]=negCount
		contingencyTable[1][1]=len(negScores)-negCount
		ob,p_value=stats.fisher_exact(contingencyTable)
		#print p_value,scoreCutoff,contingencyTable
		if p_value<bestpvalue:
			bestpvalue=p_value
			bestcuttoff=scoreCutoff
			besttable=contingencyTable[0]+contingencyTable[1]
	return bestpvalue,bestcuttoff,besttable'''


	return

def main():
	PosResultfile=argv[1]
	NegResultfile=argv[2]
	outfile=argv[3]
	motifname=argv[4]
	#outfile="test.fisher.txt"
	seqnames,posScores=loadscanresult(PosResultfile)
	seqnames,negScores=loadscanresult(NegResultfile)
	#posScores=np.random.normal(9.0,1.0,70000/100*99)
	#posScores=list(posScores)+list(np.random.normal(11.0,1.0,70000/100*1))
	#negScores=np.random.normal(9.1,1.0,70000)
	#print "Min Pos Neg",min(posScores),min(negScores)
	#print "Max Pos Neg",max(posScores),max(negScores)
	#print "Calculating Fisher p-value..."
	pvalue,scoreCutoff,table=FisherScoreCutoff(posScores,negScores)
	#print table
	out=open(outfile,'w')
	#header="MOTIF\tp-value\tscoreCutoff\tPosMatches\tPosNonMatches\tNegMatches\tNegNonMatches\tEnrichment"
	#out.write(header+'\n')
	totalPos=table[0]+table[1]
	totalNeg=table[2]+table[3]
	enrichment=table[0]/float(totalPos)/(float(table[2])/float(totalNeg))
	out.write(motifname+'\t'+str(pvalue)+'\t'+str(scoreCutoff)+'\t'+str(table[0])+'\t'+str(table[1])+'\t'+str(table[2])+'\t'+str(table[3])+'\t'+str(enrichment)+'\n')
	out.close()
	#print pvalue,scoreCutoff
	#print table
	return

if __name__=="__main__":
	main()