#!/usr/bin/env python

### This is the primary program, it works with CpG methylation only.
'''
The function is in like this 
python mepigram_typeEF.py POSset.faa NEGset.faa backgroundmodel graph output.meme maxNo.motif
'''
from sys import argv
import sys
import os
import math
import time



def sortlex(L):
    L=map(str,L)
    return map(int,sorted(L))
#### functions:
#revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','E':'F','}[B] for B in x][::-1]) #make revcompl #this this for non-CpG algorithm
#revcompl_num = lambda x: ''.join([{'0':'3','1':'2','2':'1','3':'0','4':'5','5':'4'}[B] for B in x][::-1]) #make revcompl #this this for non-CpG algorithm
basetonum = lambda x: ''.join([{'A':'0','C':'1','G':'2','T':'3','E':'4','F':'5'}[B] for B in x])
def revcompl(word):
    revComp=''
    transdict={'A':'T','C':'G','G':'C','T':'A','E':'F','F':'E'}
    index=0
    while index <len(word):
        revComp+=transdict[word[index]]
        index+=1         
    return revComp[::-1]

def revcompl_num(word):
    revComp=''
    transdict={'0':'3','1':'2','2':'1','3':'0','4':'5','5':'4'}
    index=0
    while index <len(word):
        revComp+=transdict[word[index]]
        index+=1
    return revComp[::-1]

def nonRedWord(word):
    #this function makes a reverse compliment of a word and return the 
    #first one lexicographically out of the two (original, compliment)
    revComp=revcompl(word)
    words = sorted([word, revComp]) #sort to see which one is first
    return words[0];

def getOtherEightMer(word):
    revComp =  revcompl_num(word)
    return revComp

def getMotifScore(ref2Countmatrix,offsets,wc,expansion):
#ref2Countmatrix is a dict of dict, I think
    scores = [] 
    bestScore = None
    #print "ref2countmatrix is",ref2Countmatrix
    #print "offsets is",offsets
    #print "wc is",wc
    #print "expansion is",expansion
    for offset in offsets:
        pwm = {}
        length = 0
        pcBoo = 0
        lowestScore = None 
        for i in range(offset,WORDLENGTH+offset):
            #print "i is", i
            for j in range(0,6): #use 4 because of 4 bases
                #print "j is",j
                if i not in ref2Countmatrix or j not in ref2Countmatrix[i]:
                    pcBoo = 1 
                if i in ref2Countmatrix and j in ref2Countmatrix[i] and (lowestScore == None or ref2Countmatrix[i][j] < lowestScore) and ref2Countmatrix[i][j] >0:
                    #print "refcount is",ref2Countmatrix[i][j]
                    
                    lowestScore = ref2Countmatrix[i][j]
                    #print "set the lowestscore to",lowestScore
    
        pseudocount = 0
        if pcBoo == 1:
            pseudocount = lowestScore*0.1
        #print "pseudocount is",pseudocount
        for i in range(offset, WORDLENGTH + offset):
            colTotal = 0
            colLowest = 0
            for j in range(0,6): #use 4 because of 4 bases
                if (i in ref2Countmatrix and j in ref2Countmatrix[i]):
                    colTotal +=ref2Countmatrix[i][j]
                    
            if colTotal > 0:
                colTotal += pseudocount*6 #multiply by the number of base types
                #print "coltotal is", colTotal
                for j in range(0,6): #use 4 because of 4 bases
                    if (i not in ref2Countmatrix or j not in ref2Countmatrix[i]) or ref2Countmatrix[i][j] ==0:
                        if length not in pwm:
                            pwm[length]={}
                        pwm[length][j]=pseudocount/colTotal
                        #print pseudocount/colTotal
                    else:
                        if length not in pwm:
                            pwm[length]={}
                        pwm[length][j] = (ref2Countmatrix[i][j]+pseudocount)/colTotal
                        #print pwm[length][j]
                length +=1
                #print length
        motifScore = 0
        #print "pwm is",pwm
        for i in range(EIGTHMERCOUNT):
            motifMatchScore = getMotifMatchScore(pwm,eightMer2array[i],eightMer2revArray[i])
            #print motifMatchScore
            score = motifMatchScore*preCalculatedErichmentScore[i]
            #print score
            motifScore +=score
            #print motifScore
        if bestScore == None or motifScore >bestScore:
            bestScore = motifScore
        scores+=[motifScore]

	
    return ((bestScore * math.log(wc)), scores)

def getMotifMatchScore(motif,word,wordRev):
  
    score = {}
    for i in range(WORDLENGTH):
        if (i not in motif or word[i] not in motif[i]) or (i not in motif or wordRev[i] not in motif[i]):
            #print "exit at sub getMotifMatchScore!";
            #print "i = ",i,"\nmotif =",motif,"\nword =",word,"\nrevword =",wordRev
            sys.exit(0)
        if i==0:
            score[0] = motif[i][word[i]]
            score[1] = motif[i][wordRev[i]]
        else:
            score[0] = score[0]*motif[i][word[i]]
            score[1] = score[1]*motif[i][wordRev[i]]
    if score[0]>score[1]:
        return score[0]
    else:
        return score[1]

def makeGraph():
    graph = {} #
    slideGraph = {}
    directions = {}  #note that I use a hash here whereas the orginal epigram uses a list 
    slides = {}
    
    if(os.path.exists(GRAPHDIR)==True): #if graphdir exists
	
        #print "Pre-existing graph :)\n";
    
        n2inc = {}
        try:
            print "reading from",GRAPHDIR+"/word2id.tsv"
            file = open (GRAPHDIR+'/word2id.tsv')
            #print "opened word2id.tsv"
        except:
            print "cannot open the graph"
            sys.exit(0)
        for line in file:
		
            tmp = line.strip().split('\t')
            number = int(tmp[0].strip())
            word = tmp[1].strip()
            #print number, word
            if word in wordsForGraph:
                #print number, word
                n2inc[number] = wordsForGraph[word]  #convert the number from the graph to the number in the program
           
        try:
            print "reading from ", GRAPHDIR+"/direction.tsv"
            file = open(GRAPHDIR+'/direction.tsv')
        except:
            print "can't read direction.tsv \n"
            sys.exit(0)
        for line in file:
            tmp = line.strip().split('\t')
            i = int(tmp[0])
            j = int(tmp[1])
            dir = int(tmp[2])
            #print i, j , dir 
            if i in n2inc and j in n2inc:
                if n2inc[i] not in directions:
                    directions[n2inc[i]]={}
                    #print n2inc[i]
                directions[n2inc[i]][n2inc[j]]=dir
                
            
        #print directions
        try:
            print "reading from ", GRAPHDIR+"/graph.tsv"
            file = open(GRAPHDIR+'/graph.tsv')
        except:
            print "can't read graph.tsv \n"
            sys.exit()
        #count = 0
        for line in file:
            tmp = line.strip().split('\t')
            i = int(tmp[0].strip())
            j = int(tmp[1].strip())
            
            if i in n2inc and j in n2inc:
                
                if n2inc[i] not in graph:
                    graph[n2inc[i]] = []
                    #count +=1
                
                graph[n2inc[i]] += [n2inc[j]]
                #print n2inc[i], n2inc[j], i ,j
              
        #print count
            
        try:
            print "reading from ", GRAPHDIR+"/slideGraph.tsv"
            file = open(GRAPHDIR+'/slideGraph.tsv')
        except:
            print "can't read slideGraph.tsv"
            sys.exit()
        #count = 0
        for line in file:
            tmp = line.strip().split("\t")
            i = int(tmp[0])
            j = int(tmp[1])
           
            try:
                slideGraph[n2inc[i]].append(n2inc[j])
                #count +=1
                
            except:
                if i in n2inc and j in n2inc:
                    slideGraph[n2inc[i]] = [n2inc[j]]
                    
                    #count +=1
        #print "count is", count
		

        try: 
            file = open(GRAPHDIR+'/slide.tsv')
        except:
            print "can't read slide.tsv"
            sys.exit()
        for line in file:
            tmp = line.strip().split('\t')
            i = int(tmp[0])
            j = int(tmp[1])
            slide = int(tmp[2])
            try:
                slides[n2inc[i]][n2inc[j]] = slide
            except KeyError:
                if i in n2inc and j in n2inc:
                    slides[n2inc[i]] = {}
                    slides[n2inc[i]][n2inc[j]] = slide
	

    else:
	
        print "no pre-existing graph, please create graph before running\n";
        sys.exit()
    '''
	do this later
	}'''
    #print len(n2inc)
    print len(directions),len(graph),len(slideGraph),len(slides)
    return (graph, slideGraph, directions, slides)


    
    

#Main program
if len(argv)<=1: #no arguments:
	print "Usage: python mepigram_typeEF.py foreground dinucleotide-shuffled-foreground hg19_meth/background_met-8.tsv graphdir resultMEMEfile maxmotiftofind"
	sys.exit(1)


POSITIVESEQUENCEFILE = argv[1]
NEGATIVESEQUNCEFILE = argv[2]
#generated from bgModel.pl or bgModel_met.pl #it takes the genome sequence and count the frenquency of each word = no.word +no.revword
WHOLEGENOME8MERCOUNTS = argv[3]
GRAPHDIR = argv[4]
resultFile= argv[5]
MAXMOTIFTOFIND = int(argv[6])
MODE = 'n';
DENSITYORPEAKFILE = '';
PADDINGROUNDS = 4;
EFCUTFOFF = 1.5;
#GRAPHDIR = './metgraph-8mer/';
#GRAPHDIR = './metgraph-6mer/';


# this part checks whether the arguments are good:


#WORDLENGTH=8 #hardcode, change it later
if "background_typeEF-" in WHOLEGENOME8MERCOUNTS and ".tsv" in WHOLEGENOME8MERCOUNTS:
    #WHOLEGENOME8MERCOUNTS =~ /background-(\d+)\.tsv/ ,is the variable matches the form specified?
    temp = WHOLEGENOME8MERCOUNTS.split("background_typeEF-") 
    WORDLENGTH = int(temp[1].split('.')[0])  #word length is from the name of the WHOLEGENOME8MERCOUNTS
    print "Word length is "+ str(WORDLENGTH)
elif "background_typeEF-" in WHOLEGENOME8MERCOUNTS and ".tsv" in WHOLEGENOME8MERCOUNTS:
    #WHOLEGENOME8MERCOUNTS =~ /background-(\d+)\.tsv/ ,is the variable matches the form specified?
    temp = WHOLEGENOME8MERCOUNTS.split("background_typeEF-") 
    WORDLENGTH = int(temp[1].split('.')[0])  #word length is from the name of the WHOLEGENOME8MERCOUNTS
    print "Word length is "+ str(WORDLENGTH) 
else:
    print "The format of the background file name is incorrect. Can't extract the 'word length' from it using regular expression\n\nFile name format should be \"background_typeEF-WORDLENGTH.tsv\" where WORDLENGTH = the length of the word used. Instead I found %s.\n" % WHOLEGENOME8MERCOUNTS
    sys.exit()

if resultFile=="":  #if the resultFile is not defined
    resultFile = POSITIVESEQUENCEFILE + '.epigram.out.meme'
    print "No result file given therefore writing to %s\n" %resultFile
	
if MODE=="": #if MODE is not defined
	
	MODE = 'n'
	
elif MODE == 'p' or MODE == 'dnn':
	if DENSITYORPEAKFILE=="": #DENSITYORPEAKFILE not defined
		
		print "Running in mode 'n' as no density file was specified\n"
		MODE = 'n'
		
	
if PADDINGROUNDS<=0:
	
	print "No number of widening rounds given. Going with default (4).\n";
	PADDINGROUNDS = 4
	
if EFCUTFOFF <=0.0:
	
	print "No enrichment threshold set (ET). Going with the default (1.5).\n";
	EFCUTFOFF = 1.5;
	
if GRAPHDIR=="":
	
	GRAPHDIR  = './graph-' + str(WORDLENGTH) + 'mer/'
	print "- Trying $GRAPHDIR\n"
	if os.path.exists(GRAPHDIR)!=True:  #checking the graphdir directory if there is anything
		
		print "- No graph dir for specified WORDLENGTH - therefore making my own graph. Warning: this might take a while for long words!\n";
		
	else:
		
		print "- exists :)\n";
		
	
if MAXMOTIFTOFIND==0:
	
	MAXMOTIFTOFIND = 200
#print WORDLENGTH


WORDLENGTHMINUS1 = WORDLENGTH - 1;
#below are other options that can be varied
WORDUSELIMIT = 1; #this allows the same K-mer to be used in more than 1 motif, default is 1
STARTMOTIFWORDLIMIT = 5; #this alters the number of k-mers that are needed to start a motif, default is 5
EXPANDMOTIFWORDLIMIT = 5; #this alters the number of k-mers that are needed to expand a motif, default is 5
FINALWORDCOUNTLIMIT = 20; #this alters the final number of k-mers that need to be included for a motif to be reported at the end, default is 20
FINALPOSITIONCOVERAGELIMIT = 15; #this alters the number of k-mers that need to cover each position in a motif for that position to be reported, default is 15

wholeGenome8merCounts = {}
wholeGenomeTotalCount = 0.0

file= open(WHOLEGENOME8MERCOUNTS)
for line in file:
    tmp=line.strip().split("\t") #this is to remove the \n from the string.
    eightMer = tmp[0]
    count = float(tmp[1])
    if line=="": 
        break
    if eightMer == 'TOTAL':
		wholeGenomeTotalCount = count
    else:
        eightMer = nonRedWord(eightMer)
        try: 
            wholeGenome8merCounts[eightMer] += count
        except KeyError:
            wholeGenome8merCounts[eightMer] = count
            

chr2peak2score = {}
if MODE == 'p':
	
	print "reading in %s \n" %DENSITYORPEAKFILE
	file= open (DENSITYORPEAKFILE) 
	for line in file:
		a=line.strip().split("\t")
		
        chr = a[0]
        start = a[1]
        end = a[2]
        score = a[3]
        try:
            chr2peak2score[chr][start] = score #it's a hash within a hash,not taking the $end
        except KeyError:
            chr2peak2score[chr]={}
        
elif MODE == 'dnn':
	
	print "reading in %s \n" %DENSITYORPEAKFILE
	file = open (DENSITYORPEAKFILE) 
	for line in file:
		a = line.strip().split("\t")
		chr = a[0]
		start = a[1]
        try:
            chr2peak2score[chr][start] = a[2:]  #add the reference of this array to this hash
        except KeyError:
            chr2peak2score[chr]={}
            

set2totalWords = {}
word2set2count = {}
seq2methylation = {}
totalSeq = [0,0]
proportion = {}

SEQFILES = [NEGATIVESEQUNCEFILE, POSITIVESEQUENCEFILE]
for id in range(2):
    inputFile = SEQFILES[id]
    print "Reading %s \n" %inputFile
    wordCounts = {}
    seq = 'filler'
    file = open(inputFile)
    for line in file:
        tmp=line.strip()
        if '>' == line[0]:#if the line starts with >, it's the header --> the sequence name
            seq = tmp[1:] #add the seq name to this variable
        else: #line is the sequence
            tmp = tmp.upper() #convert the line to upper case
            totalSeq[id]+=1  #increase the number of sequences in this file by one
            a = list(tmp) #splits the sequence into array of characters
            chr = ''
            pos = 0
            if((MODE == 'dnn' or MODE == 'p' or MODE == 'd') and  'chr' == seq[0:3]):
                b = seq.split("_")
                chr = b[0]
                pos = b[1] #put the chromosome and the position of the read into these
            length = len(a);   #get the length of the read
            #print "length of read is ", length
            limit = length - WORDLENGTH+1 #get the limit of how many words you can get per read
            #print limit
            thisSeqWordCounts = {} #find wordcounts for just this read
            for i in range(limit):
                end = i + WORDLENGTH #index of the end of the word
                word = ''.join(a[i:end])
                #for j in range(i,end):
                #    word = word + a[j] #put the word into $word but concating each character  ###make sure to optimize this 
                if 'N' not in word:   #if word doesnt contain N (unknown base)
                    score = 1.0; #set score =1
                    if (id == 1 and MODE != 'n'): #if id is Negativesequencefile and MODE not equal to 'n'
                        #print "mode is n"
                        if( MODE == 'p' and id == 1):
                            score =  chr2peak2score[chr][pos]; #use this score instead, because the dict is already made
                        elif((MODE == 'd') or (MODE == 'dnn' and id == 1)):
                            bin = int ((i + (WORDLENGTH/2)) / 100); #the array is in bins, so this calculate the bin of the 
                            #read and get the appropriate score. there will be about 100 scores like each other each 
                            score = chr2peak2score[chr][pos][bin];
                    usedWord = nonRedWord(word);
                    if(score < 0):
                        #print "less than 0"
                        score = 0.0
                    try:
                        thisSeqWordCounts[usedWord]+=[score]
                    except KeyError:
                         #creating a list of scores in hashtable{key} 
                        thisSeqWordCounts[usedWord]=[score]
            for word in thisSeqWordCounts.keys(): #for each word in this
                #print word
                count = 0.0
                for item in thisSeqWordCounts[word]: #for each item in this (which is 0 or 1), add them all up for each word
                    count += item
                try:
                    wordCounts[word] += count
                except KeyError:
                   #then add the number of the word's occurrence in this seq to the hash wordCounts, 
                    #which contain all words in the file
                    wordCounts[word] = count
                try:
                    proportion[word][id] += 1 #add 1 for each word in the seq
                except KeyError:
                    if word not in proportion:
                        proportion[word]=[0.0,0.0]
                        proportion[word][id] += 1.0
    
    for word in wordCounts.keys():		#for each word in the file
        try:
            word2set2count[word][id] += wordCounts[word]; #put it in another hash
        except KeyError:
            if word not in word2set2count:
                #print "word not here"
                word2set2count[word]=[0.0,0.0]
                word2set2count[word][id] += wordCounts[word]
            
        try:
            set2totalWords[id] += wordCounts[word]; #add to total words in each set
        except KeyError:
            set2totalWords[id] = wordCounts[word]
            

for id in range(2):# for each id, 
    set2totalWords[id] = float(int(set2totalWords[id] + 0.5)); #round up the values
    if set2totalWords[id] == 0:
        set2totalWords[id] = 1.0    #base case is 1
        

proportionOfPeaksWithWord = {}
'''this part is odd, fix it somehow'''
for word in word2set2count: #for word in word2setcount
    #print word
    # totalSeq[1] is the number of sequences of the positive
    #print proportion[word][1]
    try:
        if proportion[word][1]!=0:
            pro1 = float(proportion[word][1]) / totalSeq[1]
            #print proportion[word][1],word,pro1,totalSeq[1]
        else:
            pro1 = 0.1 / totalSeq[1]
        proportionOfPeaksWithWord[word] = pro1
    except KeyError:
        print "keyerror"
        pro1 = 0.1 / totalSeq[1] #if the word is not in there, use this pseudocount
    #proportionOfPeaksWithWord[word] = pro1
    #print pro1

#testing the program up to this part
sum =0
for i in proportion:
    sum+=proportion[i][1]
print sum, totalSeq[1],len(proportionOfPeaksWithWord)

'''
#testing: setting the things in proportionOfPeaksWithWord = pro1 = 0.1 / totalSeq[1]
for word in proportionOfPeaksWithWord:
    proportionOfPeaksWithWord[word]=0.1/totalSeq[1]
### this makes p-values right
'''

words = []
revWords = []
pvalues = []
efs = []
totalWordCount = 0.0
eightMer2array = []
eightMer2revArray = [] #stores 
preCalculatedErichmentScore = []
wordsAsStrings = []
for word in word2set2count.keys():
    #print word,
    count0 = 0
    count1 = 1
    #for the neg file
    if word2set2count[word][0] !=None and word2set2count[word][0] > 0:
        count0 = int (word2set2count[word][0] + 0.5)
        if count0 == 0:
            count0 = 1 # adding pseudo count 
    else:
        count0 = 1
    #do the same for the pos file
    if word2set2count[word][1] !=None and word2set2count[word][1] > 0:
        count1 = int (word2set2count[word][1] + 0.5)
        if count1 == 0:
            count1 = 1 #adding pseudo count
    else:
        count1 = 1
    try:
        enrichmentFactorWholeGenome = count1/set2totalWords[1] / (float(wholeGenome8merCounts[word]) /float(wholeGenomeTotalCount)) 
        #try this: (to exclude regions of the foreground from the background )
        #enrichmentFactorWholeGenome = count1/set2totalWords[1] / (float(wholeGenome8merCounts[word]-count1)/float(wholeGenomeTotalCount-set2totalWords[1]))
    except:
        #this means the word is not seen in the wholeGenome, set pseudocount for that =1 
        #print "error division by zero",float(wholeGenome8merCounts[word]-count1),word
        enrichmentFactorWholeGenome = count1/set2totalWords[1] / (1.0 /float(wholeGenomeTotalCount)) 

    enrichmentFactorShuffledPeaks = count1/set2totalWords[1] / (count0/set2totalWords[0]);
    
    try:
        dualEnrichmentFactor = (count1/set2totalWords[1]) / (((count0/set2totalWords[0]) + (float(wholeGenome8merCounts[word]) /wholeGenomeTotalCount)) / 2);
    except KeyError:
        dualEnrichmentFactor = (count1/set2totalWords[1]) / (((count0/set2totalWords[0]) + (1.0 /wholeGenomeTotalCount)) / 2);
    
    if(enrichmentFactorWholeGenome > 1 and enrichmentFactorShuffledPeaks > 1):
        wordsAsStrings.append(word) #add the word to word as string list
        
        pvalue =  proportionOfPeaksWithWord[word] * (math.log(enrichmentFactorWholeGenome) + math.log(enrichmentFactorShuffledPeaks))
        #the p-value actcally refers to 'W' - the weight
        word = basetonum(word) #translate the word to 0123
        #print word
        revCompWord = getOtherEightMer(word)
        #print revCompWord   this is same as perl
        word = list(word)
        revCompWord = list(revCompWord)
        words.append(word)
        revWords.append(revCompWord)
        pvalues.append(pvalue)
        efs.append(dualEnrichmentFactor)
        totalWordCount+=1
        #print pvalue
    if(dualEnrichmentFactor > EFCUTFOFF or dualEnrichmentFactor < (1/EFCUTFOFF)):
        
        try: 
            word = basetonum(word)
            w = list(word)
        except KeyError:
            w = list(word)
        eightMer2array.append(w)
        otherEightMer = getOtherEightMer(word);
        #print otherEightMer
        o = list(otherEightMer)
        eightMer2revArray.append(o)
        score = (math.log(enrichmentFactorWholeGenome) + math.log(enrichmentFactorShuffledPeaks));
        #print round(score,2)
        preCalculatedErichmentScore.append(score)
#print len(wordsAsStrings)

###testing
sum=0.0
for i in preCalculatedErichmentScore:
    sum+=i
print sum
print len(preCalculatedErichmentScore)
#wholeGenome8merCounts['AAAAAA']
#done

#adding preCalculated enrichemnt score normalisation
postiveSum = 0.0
negativeSum = 0.0
for score in preCalculatedErichmentScore:
	
	if(score > 0):
		
		postiveSum += score;
		#print "pos",postiveSum
	else:
		
		negativeSum += score;
		#print "neg", negativeSum
if negativeSum !=0.0:
    mulitplier = postiveSum / abs(negativeSum)
    print mulitplier

    for i in range(len(preCalculatedErichmentScore)):
        if(preCalculatedErichmentScore[i] < 0):
       
            preCalculatedErichmentScore[i] *= mulitplier;
		
print negativeSum

####### the program starts here
EIGTHMERCOUNT = len(preCalculatedErichmentScore)
print "Words read as input = ",totalWordCount
print "Total words to be used in scoring = ",EIGTHMERCOUNT
print "Sorting words by weights (W)."
print "EF = ",EFCUTFOFF,"and PR = ",PADDINGROUNDS," The result file will be", resultFile
e = totalWordCount - 1;
si = sorted(range(len(pvalues)), key=lambda k: pvalues[k],reverse=True)  #sort by pvalues, get the index
#then reorder the other arrays based on that index
words = [words[i] for i in si]
revWords = [revWords[i] for i in si]
pvalues = [pvalues[i] for i in si]
#print pvalues
efs = [efs[i] for i in si]
#for p in efs:
#    print p
wordsAsStrings = [wordsAsStrings[i] for i in si]

wordsForGraph ={} 
for wordID in range(int(totalWordCount)):
	
	w = wordsAsStrings[wordID]
	wordsForGraph[w]= wordID
	
#wordsAsStrings = []; #empty it

print "Getting the graph\n";
(graph,slideGraph,directions,slides) = makeGraph(); #make the graphs

motifs = [] 
usedWords = {} #hash
motifScores = [] 
motifsWordsUsed = []
motif2length = [] 
colCounts = []
countMatrices = []
wordsUsed = []
motifCount = 0.0
print "lengraph, lendirections, lenslidegraph, lenslides are "

print "Making the motifs\n";
#print map(int,words)
#should make words and revwords as list of numbers
#also the eightmer2array and eightmer2revarray
for i in range(len(words)):
    words[i]=map(int,words[i])
for i in range(len(revWords)):
    revWords[i]=map(int,revWords[i])
for i in range(len(eightMer2array)):
    eightMer2array[i]=map(int,eightMer2array[i])
for i in range(len(eightMer2revArray)):
    eightMer2revArray[i]=map(int,eightMer2revArray[i])
print totalWordCount

#ids=[str(i) for i in range(int(totalWordCount))]
#print ids
time_start = time.time()
for wordID in range(int(totalWordCount)): # the loop 'WORD'
    #print "WORDID is",wordID,"--------------------------------"
    if((wordID not in usedWords or usedWords[wordID] < WORDUSELIMIT) and efs[wordID] > EFCUTFOFF):
        #the most enriched remaining word is taken and a motif is made
        
        newMotif = {} 
        colCount = {}
        for j in range(WORDLENGTH):
            i = j + PADDINGROUNDS
            if i not in newMotif:
                newMotif[i]={}
            #if words[wordID][j]==4: print "YAHOOOOOOOO"
            if words[wordID][j] not in newMotif[i]:
                newMotif[i][words[wordID][j]] = 0
            newMotif[i][words[wordID][j]] += pvalues[wordID]
            try: 
                colCount[i]+=1
            except KeyError:
                colCount[i] = 1
            #print "colCount",colCount[i]
        #print "NEWMOTIF is ",newMotif
        """checked"""
        words2include = {}
        words2slide = {}
        slideWord2direction = {}
        oneMisMatchWords = {}
        #print "graph is ",graph[wordID]
        if wordID not in graph:
            #print "error 1311 wordID is not in graph", wordID
            continue
        for matchWordID in sortlex(graph[wordID]): 
            #print "matchword ID ",matchWordID
            '''checked'''
            if (matchWordID not in usedWords or usedWords[matchWordID]< WORDUSELIMIT):
                words2include[matchWordID] = 0
                oneMisMatchWords[matchWordID] = 0
                for match2matchWordID in sortlex(graph[matchWordID]):
                    if (match2matchWordID not in usedWords or usedWords[match2matchWordID] < WORDUSELIMIT) and match2matchWordID != wordID:
                        words2include[match2matchWordID] = 0
                        #print "match2matchWordID",match2matchWordID
        #print words2include
        #print oneMisMatchWords
        '''checked'''
        #this part identifies words which are one slide away from the alignemnt
        if PADDINGROUNDS >= 1:
            if wordID not in slideGraph:
                #print "error 1312 wordID not in slideGraph", wordID
                continue
            for slideWordID in sortlex(slideGraph[wordID]): #first we just get slide words from the main word
                if (slideWordID not in usedWords or usedWords[slideWordID] < WORDUSELIMIT) and slideWordID not in words2include and slideWordID not in words2slide and slideWordID != wordID:
                    #print "checking",slideWordID
                    slideWord2direction[slideWordID] = directions[wordID][slideWordID]
                    words2slide[slideWordID] = slides[wordID][slideWordID]
                    #print slideWordID,"slidewordID"
                    if slideWordID not in graph:
                        #print "not defined"
                        continue
                    for match2theSlideWordID in sortlex(graph[slideWordID]):
                        #print "match2theSlideWordID",match2theSlideWordID
                        if (match2theSlideWordID not in usedWords or usedWords[match2theSlideWordID] < WORDUSELIMIT) and match2theSlideWordID not in words2include and match2theSlideWordID not in words2slide and match2theSlideWordID != wordID:
                            direction1 = slideWord2direction[slideWordID]
                            direction2 = directions[slideWordID][match2theSlideWordID]
                            direction = 0 
                            if direction1 != direction2:
                                direction = 1 
                            words2slide[match2theSlideWordID] = words2slide[slideWordID]
                            slideWord2direction[match2theSlideWordID] = direction
            #print "words2slide",words2slide
            #print "slideWord2direction",slideWord2direction
            word2checkForSlide = sortlex(oneMisMatchWords.keys())
           
            #hardcode this to test:
            #word2checkForSlide =[124,60]
            #print "word2checkForSlide",word2checkForSlide
            for matchWordID in sortlex(word2checkForSlide):
                direction1 = directions[wordID][matchWordID]
                #print "Loop with matchwordID as-------------", matchWordID
                #print "direction1",direction1
                if (matchWordID not in slideGraph): 
                    #print "error 1313 thsis number is not in slideGraph", matchWordID
                    continue
                for slideWordID in sortlex(slideGraph[matchWordID]):
                    if (slideWordID not in usedWords or usedWords[slideWordID] < WORDUSELIMIT) and slideWordID not in words2include and slideWordID not in words2slide and slideWordID != wordID:
                        direction2 = directions[matchWordID][slideWordID]
                        #print "direction2",direction2
                        if direction1 == direction2:
                            slideWord2direction[slideWordID] = 0
                            #print "EQUAL, Adding",slideWordID,"to slideword2direction"
                        else:
                            slideWord2direction[slideWordID] = 1
                            #print "NOT EQUAL, adding",slideWordID,"to slideword2direction"
                        if direction1 ==1:
                            s = -1
                            if slides[matchWordID][slideWordID] == -1:
                                s = 1
                            #print "s is",s
                            words2slide[slideWordID] = s
                        else:
                            words2slide[slideWordID] = slides[matchWordID][slideWordID]
                            #print "word2slide adding",words2slide[slideWordID], matchWordID
                        
                        #print "gets things one mismatch from the slide"
                        if slideWordID not in graph:
                            #print "error 1314 slideWordID not in graph",slideWordID
                            continue
                        for slide2match in sortlex(graph[slideWordID]):
                            if (slide2match not in usedWords or usedWords[slide2match]<WORDUSELIMIT) and slide2match not in words2include and slide2match not in words2slide and slide2match!=wordID:
                                words2slide[slide2match] = words2slide[slideWordID]
                                #print "word2slide slide2match",words2slide[slide2match],slide2match
                                #print slideWord2direction[slideWordID], 'compare', directions[slideWordID][slide2match]
                                if (slideWord2direction[slideWordID]) == int(directions[slideWordID][slide2match]):
                                    slideWord2direction[slide2match] = 0
                                    #print "adding", slide2match,"to slideword2direction at 0"
                                else:
                                    slideWord2direction[slide2match] = 1
                                    #print "adding",slide2match,"to slideword2direction at 1"
                        #print "slideword2direction",slideWord2direction
                        for slideWordID2 in sortlex(slideGraph[slideWordID]):
                            if PADDINGROUNDS >= 2 and (slideWordID2 not in usedWords or usedWords[slideWordID2] < WORDUSELIMIT) and slideWordID2 not in words2include and slideWordID2 not in words2slide and slideWordID2 != wordID:
                                #sorts out the direction
                                #print "PADDINGROUNDS is", PADDINGROUNDS
                                if directions[slideWordID][slideWordID2] == 0:
                                    slideWord2direction[slideWordID2] = slideWord2direction[slideWordID]
                                    #print "adding21",slideWordID2,"to slideWord2direction at",slideWord2direction[slideWordID]
                                elif slideWord2direction[slideWordID] == 1:
                                    slideWord2direction[slideWordID2] = 0
                                    #print "adding22",slideWordID2,"to slideWord2direction at 0"
                                else:
                                    slideWord2direction[slideWordID2] = 1
                                    #print "adding23",slideWordID2,"to slideWord2direction at 1"
                                #sort out the slide:
                                if words2slide[slideWordID] == -1:
                                    words2slide[slideWordID2] = -2
                                else:
                                    words2slide[slideWordID2] = 2
                                if slideWordID2 not in graph:
                                    #print "error 1315 slideWordID2 not in graph",slideWordID2
                                    continue
                                for slide2match in sortlex(graph[slideWordID2]):
                                    #print "slide2match is",slide2match
                                    if (slide2match not in usedWords or usedWords[slide2match] < WORDUSELIMIT) and slide2match not in words2include and slide2match not in words2slide and slide2match != wordID:
                                        words2slide[slide2match] = words2slide[slideWordID2]
                                        #print "MARKER"
                                        if slideWord2direction[slideWordID2] == directions[slideWordID2][slide2match]:
                                            slideWord2direction[slide2match] = 0
                                            #print "adding24",slide2match,"to slideWord2direction at 0"
                                        else:
                                            slideWord2direction[slide2match] = 1
                                            #print "adding25",slide2match,"to slideWord2direction at 1"
                                for slideWordID3 in sortlex(slideGraph[slideWordID2]):
                                    if PADDINGROUNDS >= 3 and (slideWordID3 not in usedWords or usedWords[slideWordID3] < WORDUSELIMIT) and slideWordID3 not in words2include and slideWordID3 not in words2slide and slideWordID3 != wordID:
                                        #sorts out the direction:
                                        #print "PADDINGROUNDS is", PADDINGROUNDS
                                        #print "directions[slideWordID2][slideWordID3]",slideWordID2, slideWordID3, directions[slideWordID2][slideWordID3]
                                        if directions[slideWordID2][slideWordID3] == 0:
                                            slideWord2direction[slideWordID3] = slideWord2direction[slideWordID2]
                                            #print "adding31",slideWordID3,"to slideWord2direction at",slideWord2direction[slideWordID2]
                                        elif slideWord2direction[slideWordID2] == 1:
                                            slideWord2direction[slideWordID3] = 0
                                            #print "adding32",slideWordID3,"to slideWord2direction at 0"
                                        else:
                                            slideWord2direction[slideWordID3] = 1
                                            #print "adding33",slideWordID3,"to slideWord2direction at 1"
                                        #sorts out the slide:
                                        if words2slide[slideWordID] == -1:
                                            words2slide[slideWordID3] = -3
                                        else:
                                            words2slide[slideWordID3] = 3
                                        if slideWordID3 not in graph:
                                            #print "error 1316 slideWordID3 not in graph",slideWordID3
                                            continue
                                        for slide2match in sortlex(graph[slideWordID3]):
                                            if (slide2match not in usedWords or usedWords[slide2match] < WORDUSELIMIT) and slide2match not in words2include and slide2match not in words2slide and slide2match != wordID:
                                                words2slide[slide2match] = words2slide[slideWordID3]
                                                if slideWord2direction[slideWordID3] == directions[slideWordID3][slide2match]:
                                                    slideWord2direction[slide2match] = 0
                                                    #print "adding34",slide2match,"to slideWord2direction at 0"
                                                else:
                                                    slideWord2direction[slide2match] = 1
                                                    #print "adding35",slide2match,"to slideWord2direction at 1"
                                        
                                    
                                            
                                        for slideWordID4 in sortlex(slideGraph[slideWordID3]):
                                            if PADDINGROUNDS >=4 and (slideWordID4 not in usedWords or usedWords[slideWordID4] < WORDUSELIMIT) and slideWordID4 not in words2include and slideWordID4 not in words2slide and slideWordID4 != wordID:
                                                #sorts out the direction
                                                #print "PADDINGROUNDS is", PADDINGROUNDS
                                                if directions[slideWordID3][slideWordID4] == 0:
                                                    slideWord2direction[slideWordID4] = slideWord2direction[slideWordID3]
                                                    #print "adding41",slideWordID4,"to slideWord2direction at",slideWord2direction[slideWordID3]
                                                elif slideWord2direction[slideWordID3] == 1:
                                                    slideWord2direction[slideWordID4] = 0
                                                    #print "adding42",slideWordID4,"to slideWord2direction at 0"
                                                else:
                                                    slideWord2direction[slideWordID4] = 1
                                                    #print "adding43",slideWordID4,"to slideWord2direction at 1"
                                                #sort out the slide
                                                if words2slide[slideWordID] == -1:
                                                    words2slide[slideWordID4] = -4
                                                else:
                                                    words2slide[slideWordID4] = 4
                                                if slideWordID4 not in graph:
                                                    #print "error 1317 slideWordID4 not in graph", slideWordID4
                                                    continue
                                                for slide2match in sortlex(graph[slideWordID4]):
                                                    if (slide2match not in usedWords or usedWords[slide2match] < WORDUSELIMIT) and slide2match not in words2include and slide2match not in words2slide and slide2match!=wordID:
                                                        words2slide[slide2match] = words2slide[slideWordID4]
                                                        if slideWord2direction[slideWordID4] == directions[slideWordID4][slide2match]:
                                                            slideWord2direction[slide2match] = 0
                                                            #print "adding44",slide2match,"to slideWord2direction at 0"
                                                        else:
                                                            slideWord2direction[slide2match] = 1
                                                            #print "adding45",slide2match,"to slideWord2direction at 1"
        
    
        #if enough words are identified they are added to the motif
        wordsIncluded = len(words2include.keys()) + len(words2slide) +1 #it's +1 because there is the word which instilaized it as well
        wordsForTheFirstPart = len(words2include) + 1
        #print "wordsIncluded",wordsIncluded
        #print "MOTIFs before the if",motifs 
        if wordsForTheFirstPart > STARTMOTIFWORDLIMIT:
            #here words one or two mistmatches from the main word are added
            for word2include in sortlex(words2include):
                if word2include not in directions[wordID] or directions[wordID][word2include] == 0:
                    for j in range(WORDLENGTH):
                        i = j + PADDINGROUNDS
                        if i not in newMotif:
                            newMotif[i] = {}
                        if words[word2include][j] not in newMotif[i]:
                            newMotif[i][words[word2include][j]] = 0
                        newMotif[i][words[word2include][j]] +=pvalues[word2include]
                        try: 
                            colCount[i]+= 1
                        except KeyError:
                            colCount[i] = 1
                else:
                    for j in range(WORDLENGTH):
                        i = j + PADDINGROUNDS
                        if i not in newMotif:
                            newMotif[i] = {}
                        if revWords[word2include][j] not in newMotif[i]:
                            newMotif[i][revWords[word2include][j]] = 0
                        newMotif[i][revWords[word2include][j]] += pvalues[word2include]
                        try: 
                            colCount[i]+= 1
                        except KeyError:
                            colCount[i] = 1
                    
            #print "NEWMOTIF BEFORE THE IF",newMotif
            '''checked that newmotif is correct'''
            slides2 = [PADDINGROUNDS]   #not sure about this  my $slides = [$PADDINGROUNDS];
            #print "slide is", slides2, ", padding is", PADDINGROUNDS
    
            (intialMotifScore,partScores) = getMotifScore(newMotif,slides2,wordsForTheFirstPart,1)
            #print "Initial Motif score =", intialMotifScore
            #print "partScores is",partScores
            if intialMotifScore>0:
                wordsIncludedInThisMotif = []
                wordCount = wordsForTheFirstPart
                #print "wordCount is", wordCount
                for word2include in sortlex(words2include):
                    wordsIncludedInThisMotif+=[word2include]
                    #print "pushing",word2include
                bestMotifScore = intialMotifScore
                bestPartScores = partScores
                bestMotif = newMotif
                #here the slide words are added
                goRight = [1,2,3,4]
                goLeft = [-1,-2,-3,-4]
                slidesToTry = []
                slidesToTry +=[goRight]
                slidesToTry +=[goLeft]
                #print slidesToTry
                slidesBeingUsed = [PADDINGROUNDS] #not sure my @slidesBeingUsed = ($PADDINGROUNDS);
                expansionCount = 0;
                for direction in range(2):
                    for slideToTry in sortlex(slidesToTry[direction]): #EXPAND loop
                        #print "slideToTry is --------------",slideToTry
                        motif = bestMotif.copy()
                        roundColCount = colCount.copy()
                        addingThisRoundCount = 0
                        thisRoundsWords = []
                        for  word2slide in sortlex(words2slide):
                            if slideToTry != words2slide[word2slide]:
                                #print "slideToTry, words2slide[word2slide]",slideToTry, words2slide[word2slide]
                                continue #skip this iteration, use the next one
                            addingThisRoundCount += 1
                            slide = words2slide[word2slide] + PADDINGROUNDS
                            
                            line = {}
                            thisRoundsWords+=[word2slide]
                            if slideWord2direction[word2slide] == 0:
                                for j in range(WORDLENGTH):
                                    i = j + slide 
                                    #print "i j slide wordlength are",i, j, slide,WORDLENGTH
                                    try:
                                        motif[i][words[word2slide][j]] += pvalues[word2slide]
                                    except KeyError:
                                        try:
                                            motif[i][words[word2slide][j]] =pvalues[word2slide]
                                        except KeyError:
                                            motif[i] = {}
                                            motif[i][words[word2slide][j]] =pvalues[word2slide]
                                            #print "adding new values to motif",i,words[word2slide][j],pvalues[word2slide]
                                    line[i] = words[word2slide][j]
                                    try:
                                        roundColCount[i]+=1
                                    except KeyError:
                                        #print i, "is not in roundcolcount1, adding it!!"
                                        roundColCount[i] = 1
                                    #print "roundColCount1 i",roundColCount[i],i
                            else:
                                for j in range(WORDLENGTH):
                                    i = j + slide
                                    #print "revwordsword2slidej", revWords[word2slide][j]
                                    #print "motif is ",motif[i]
                                    #if the key is not in the motif, add it to it?
                                    try:
                                        motif[i][revWords[word2slide][j]] +=pvalues[word2slide]
                                    except KeyError:
                                        try:
                                            motif[i][revWords[word2slide][j]] =pvalues[word2slide]
                                        except KeyError:
                                            #motif[i]={}
                                            motif[i]={}
                                            motif[i][revWords[word2slide][j]] =pvalues[word2slide]
                                    #print "i is", i, "word2slide is ",revWords[word2slide][j],"j is",j
                                    line[i] = revWords[word2slide][j]
                                    try: 
                                        roundColCount[i]+=1
                                    except KeyError:
                                        #print i,"is not in roundcolcount2, adding it!!"
                                        roundColCount[i] = 1
                                    #print "roundColCount2 i",roundColCount[i],i
                        #print "MOTIF IS ",motif
                        #print "WORDCOUNT IS",wordCount
                        #print "addingThisRoundCount is",addingThisRoundCount
                        if addingThisRoundCount > EXPANDMOTIFWORDLIMIT:
                            slidesBeingUsed += [PADDINGROUNDS +slideToTry]
                            #print "PADDINGROUNDS +slideToTry",PADDINGROUNDS +slideToTry
                            (updatedMotifScore, updatedPartScores) = getMotifScore(motif, slidesBeingUsed,(wordCount + addingThisRoundCount), (expansionCount +2))
                            #print bestMotifScore,'->',updatedMotifScore,"after adding",addingThisRoundCount
                            parts = len(bestPartScores)
                            #it checks that the expansion hasn't decreased the current score of any of the parts
                            boo = 0
                            for i in range(parts):
                                if updatedPartScores[i] < bestPartScores[i]:
                                    boo = 1
                                    #print "boo is 1"
                            if boo == 1 or updatedPartScores[parts-1] <0:
                                slidesBeingUsed.pop()
                                #print "endinng this EXPAND round ROUNDCOLCOUND 5 IS",roundColCount[5],colCount[5]
                                break #this breaks out of the EXPAND loop
                            else:
                                for w in sortlex(thisRoundsWords):
                                    wordsIncludedInThisMotif+=[w]
                                bestMotifScore = updatedMotifScore
                                bestPartScores = updatedPartScores
                                bestMotif = motif
                                colCount = roundColCount
                                #print "ROUNDCOLCOUNT 5 IS",roundColCount[5]
                                wordCount += addingThisRoundCount
                                expansionCount += 1
                #print "-Final motif score =",bestMotifScore,"with",wordCount,"words"
                if wordCount <FINALWORDCOUNTLIMIT:
                    continue #skip the rest, go to the next WORD iteration
                wordsIncludedInThisMotif+=[wordID]
                for w in sortlex(wordsIncludedInThisMotif):
                    try:
                        usedWords[w]+=1
                    except KeyError:
                        #print "error 1317", w,"not in usedwords"
                        usedWords[w]=1
                
                newPwm = {}
                colTotals = []
                finalColMatrix = {}
                length = 0
                colMaxPtotal = 0
                #print "bestmotif is", bestMotif
                for i in range(WORDLENGTH +2*PADDINGROUNDS):
                    colTotal = 0
                    #print "i is",i
                    for j in range(6): #use 5 instead of 4
                        if i in bestMotif and j in bestMotif[i]:
                            #print "bestmotif",i,j,bestMotif[i][j]
                            colTotal += bestMotif[i][j]
                        else:
                            try:
                                bestMotif[i][j] = 0
                            except KeyError:   #check this 
                                #print "key error"
                                bestMotif[i]={}
                                bestMotif[i][j] = 0
                    colTotals+=[colTotal]
                    if colMaxPtotal<colTotal:
                        colMaxPtotal = colTotal
                    #print "colTotal is", colTotal
                    if colTotal >0:
                        #print "length1 is", length
                        for j in range(6):#use 5 instead of 4
                            
                            if i not in bestMotif or j not in bestMotif[i]:
                                try: 
                                    newPwm[length][j] = 0
                                    #print "newPwm1 i j", length,j,0
                                except KeyError:
                                    newPwm[length]={}
                                    newPwm[length][j] = 0
                                    #print "newPwm2 i j", length,j,0
                                try: 
                                    finalColMatrix[length][j] = 0
                                except KeyError:
                                    finalColMatrix[length] = {}
                                    finalColMatrix[length][j] = 0
                            else:
                                try: 
                                    newPwm[length][j] = bestMotif[i][j] / colTotal;
                                    #print "newPwm3 i j", length,j,bestMotif[i][j],  colTotal
                                except KeyError:
                                    newPwm[length]={}
                                    newPwm[length][j] = bestMotif[i][j] / colTotal;
                                    #print "newPwm4 i j", length,j,bestMotif[i][j],  colTotal
                                try: 
                                    finalColMatrix[length][j] = bestMotif[i][j]
                                except KeyError:
                                    finalColMatrix[length] = {}
                                    finalColMatrix[length][j] = bestMotif[i][j]
                            newPwm[length][j] = round(newPwm[length][j],6)
                        length+=1
                        #print "length is", length
                        
    
                #print "NEWPWM IS", newPwm
        
                
                col2remove = {}
                for i in range(WORDLENGTH + (PADDINGROUNDS*2)):
                    if i not in colCount:
                        colCount[i] = 0
                #find positions which arent covered by many word
                pos = 0
                for i in range(WORDLENGTH + (PADDINGROUNDS*2)):
                    #print "colCount",i,"FINALPOSITIONCOVERAGELIMIT",FINALPOSITIONCOVERAGELIMIT,colCount[i]
                    if colCount[i] < FINALPOSITIONCOVERAGELIMIT and colCount[i]>0:
                        #we don't want to removed the 0 count positions becasue they were never included in the pwm
                        print "!!removing", i,"(", pos,") as col count =",colCount[i],"!!"
                        col2remove[pos] = 0
                    if colCount[i]>0:
                        pos+=1
                countOfCols2remove = len(col2remove)
                colsDeleted = 0
                for col in sortlex(col2remove):
                    col = col - colsDeleted
                    #print "deleting column",col
                    #del newPwm[col]
                    #del finalColMatrix[col]
                    #del colCount[col]
                    colsDeleted+=1
                #print "length countOfCols2remove",length,countOfCols2remove
                length = length - countOfCols2remove
                #print colCount
                tmp  = [str(colCount[i]) for i in colCount]
                #print tmp
                
                colCount = '-'.join(tmp)
                colCounts+=[colCount]
                motifs+=[newPwm]
                motifScores+=[bestMotifScore]
                motifsWordsUsed+=[wordCount]
                motif2length+=[length]
                #print motif2length
                countMatrices+=[finalColMatrix]
                wordsUsed+=[wordsIncludedInThisMotif]
                motifCount+=1
                if motifCount == MAXMOTIFTOFIND:
                    break #stop looking for more motifs 
                
                    
print "Motif count =",motifCount;
print "Time elapsed",time.time()-time_start


#output the results:
e = motifCount -1 

si = sorted(range(len(motifScores)), key=lambda k: motifScores[k],reverse=True)  
#sort by motifScores, get the index
#then reorder the other arrays based on that index
motifs = [motifs[i] for i in si]

motifScores = [motifScores[i] for i in si]
motifsWordsUsed = [motifsWordsUsed[i] for i in si];
motif2length = [motif2length[i] for i in si];
colCounts = [colCounts[i] for i in si]
countMatrices = [countMatrices[i] for i in si];
wordsUsed = [wordsUsed[i] for i in si]

'''Attention!!!: check whether to use motif2length or length(motifs[motifID])'''

for p in range(len(motifScores)):
	
	motifScores[p] = round(motifScores[p],3) #add 3 decimal places to the end of each p
	
HEADER="mEpigram Type EF\nALPHABET= ACGTEF\n"
outfile=open(resultFile,'w')
print resultFile
outfile.write(HEADER+"\n")
#print HEADER
for motifID in range(len(motifScores)): 
    motifname="MOTIF"+" "+str(motifID)+"_"+str(motifScores[motifID])#+"_"+str(colCounts[motifID])+"_"+str(motif2length[motifID])+"_"+str(motifsWordsUsed[motifID])
    #print colCounts[motifID]
    #print motifname
    head="letter-probability matrix: alength= "+str(6)+" w= "+str(len(motifs[motifID]))+" nsites= "+str(motifsWordsUsed[motifID])+" score= "+str(motifScores[motifID])
    #print head
    outfile.write(motifname+"\n")
    outfile.write(head+"\n")
    for i in range(len(motifs[motifID])):
        line =""
        for j in range(6):
            line+=str(motifs[motifID][i][j])+"\t"
        line =line.strip()
        #print line
        outfile.write(line+"\n")
    #print "\n"
    outfile.write("\n")
outfile.close()








