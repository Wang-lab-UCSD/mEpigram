#!/usr/bin/env python
"""
INPUT: typeE or typeEF MEME files. 

OUTPUT: motif logo for each PWM
"""

import os
import random as rd
from sys import argv
import os
import random as rd
#from sys import argv
import sys
def load_motifs_typeEF(filename):
    file=open(filename)
    seq=file.read().split("MOTIF")
    seq=seq[1:]
    motifs={}
    infos={}
    for s in range(len(seq)):
        t=seq[s].strip().split("\n")
        #name=int(t[0].split('_')[0]) #load only the number
        name=t[0] #load whole name
        motifs[name]=t[2:]
        infos[name]=t[0:2]
    for m in motifs:
        tdict={'A':[],'C':[],'G':[],'T':[],'E':[],'F':[]}
        for pos in range(len(motifs[m])):
            tmp=motifs[m][pos].strip().split("\t")
            if len(tmp)!=6:
                print "ERROR: Motif matrix does not contain 6 collumn. It is not typeEF." 
                sys.exit(1)
            tdict['A']+=[float(tmp[0])]
            tdict['C']+=[float(tmp[1])]
            tdict['G']+=[float(tmp[2])]
            tdict['T']+=[float(tmp[3])]
            tdict['E']+=[float(tmp[4])]
            tdict['F']+=[float(tmp[5])]
        motifs[m]=tdict
    
    return motifs,infos

def randomkmer_typeEF(PWM):
    #take a PWM and randomly create a Kmer based on it
    # chance goes A->C->G->T->E
    kmer=""
    for pos in range(len(PWM['A'])):
        die=rd.randint(0,100000)
        if die<=PWM['A'][pos]*100000:
            kmer+='A'
        elif die<=(PWM['A'][pos]+PWM['C'][pos])*100000:
            kmer+='C'
        elif die<=(PWM['A'][pos]+PWM['C'][pos]+PWM['G'][pos])*100000:
            kmer+='G'
        elif die<=(PWM['A'][pos]+PWM['C'][pos]+PWM['G'][pos]+PWM['T'][pos])*100000:
            kmer+='T'
        elif die<=(PWM['A'][pos]+PWM['C'][pos]+PWM['G'][pos]+PWM['T'][pos]+PWM['E'][pos])*100000:
            kmer+='E'
        else:
            kmer+='F'
    return kmer

def load_motifs_typeE(filename):
    file=open(filename)
    seq=file.read().split("MOTIF")
    seq=seq[1:]
    motifs={}
    infos={}
    for s in range(len(seq)):
        t=seq[s].strip().split("\n")
        #name=int(t[0].split('_')[0]) #load only number
        name = t[0] #load whole name
        motifs[name]=t[2:]
        infos[name]=t[0:2]
    for m in motifs:
        tdict={'A':[],'C':[],'G':[],'T':[],'E':[],}
        for pos in range(len(motifs[m])):
            tmp=motifs[m][pos].strip().split("\t")
            if len(tmp)!=5:
                print "ERROR: Motif matrix does not contain 5 collumn. It is not typeE." 
                sys.exit(1)
            tdict['A']+=[float(tmp[0])]
            tdict['C']+=[float(tmp[1])]
            tdict['G']+=[float(tmp[2])]
            tdict['T']+=[float(tmp[3])]
            tdict['E']+=[float(tmp[4])]
        motifs[m]=tdict
    
    return motifs,infos
    
def randomkmer_typeE(PWM):
    #take a PWM and randomly create a Kmer based on it
    # chance goes A->C->G->T->E
    kmer=""
    for pos in range(len(PWM['A'])):
        die=rd.randint(0,100000)
        if die<=PWM['A'][pos]*100000:
            kmer+='A'
        elif die<=(PWM['A'][pos]+PWM['C'][pos])*100000:
            kmer+='C'
        elif die<=(PWM['A'][pos]+PWM['C'][pos]+PWM['G'][pos])*100000:
            kmer+='G'
        elif die<=(PWM['A'][pos]+PWM['C'][pos]+PWM['G'][pos]+PWM['T'][pos])*100000:
            kmer+='T'
        else:
            kmer+='E'
    return kmer

def main():
    #memefile=argv[1]
    #resultdir=argv[2]
    #nametag=argv[3] #use this to add to the front of the motif
    typeEF=False
    memefile=None
    resultdir=None
    nametag=""


    usage = """USAGE: 
    %s [options]

        -m <filename>       input motif file
        --typeEF            convert the modified base to E and the complementary base to F, default is False (optional)
        -t                  name tag to add to the motif files (optional)
        -o <output dir>     name of the output dir to be created (eg. hg19_methylated)
        -h                  print this usage message
    """ % (sys.argv[0])

    # parse command line
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if (arg == "-m"):
            i += 1
            try: memefile = sys.argv[i]
            except:
                print " Error in -m" 
                sys.exit(1)
        elif (arg == "-t"):
            i += 1
            try: nametag = float(sys.argv[i])
            except:
                print "Error in -t" 
                sys.exit(1)
        elif (arg == "--typeEF"):
            #i += 1
            #print "Is Bed"
            try: typeEF = True
            except:
                print "Error in --typeEF" 
                sys.exit(1)
        elif (arg == "-o"):
            i += 1
            try: resultdir = sys.argv[i]
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
    if (memefile == None):
        print "No input, exit."
        print >> sys.stderr, usage; sys.exit(1)
        sys.exit(1)
    elif (resultdir == None):
        print "No output, exit."
        sys.exit(1)
    



    if typeEF:
        print "type EF"
        os.system("mkdir "+resultdir)
        motifs,infos=load_motifs_typeEF(memefile)
        #for m in range(len(motifs)): #only print the top 25 
        for m in motifs:
            print "motif",m
            PWM=motifs[m]
            name=infos[m][0] 
            name=name.replace("MOTIF",'').strip()

            if "m-motif" in name:
                name=nametag+"_"+name.split("_")[0]+"_m-motif"
            else:
                name=nametag+"_"+name.split("_")[0]
                
            name = m #use motif name, NO TAG

            target=open("tmp.weblogo.faa",'w')
            for i in range(3000):
                kmer=randomkmer_typeEF(PWM)
                target.write(">"+str(i)+'\n')
                target.write(kmer+"\n")
            target.close()
            #os.system("weblogo --format png --size large --title 'Motif "+str(m)+"' -a ACGTE  --color '#CC0000' A 'a' --color '#008000' T 'Py' --color '#0000CC' C 'c' -C '#CC00CC' E 'e' -C '#FFB300' G 'g' <tmp.faa > " +resultdir+"/Motif"+str(m)+".png")
            os.system("weblogo --format eps --size large -a ACGTEF  --color '#CC0000' A 'a' --color '#008000' T 'Py' --color '#0000CC' C 'c' -C '#CC00CC' E 'e' -C '#FFB300' G 'g' -C '#66ffff' F 'f' <tmp.weblogo.faa> " +resultdir+"/"+name+".eps")
            #os.system("weblogo --format eps --size large -a ACGTEF  --color '#CC0000' A 'a' --color '#008000' T 'Py' --color '#0000CC' C 'c' -C '#CC00CC' E 'e' -C '#FFB300' G 'g' -C '#66ffff' F 'f' -D transfac <tmp.weblogo.transfact.dat> " +resultdir+"/"+name+".eps")

        os.system('rm tmp.weblogo.faa')

    else:

        os.system("mkdir "+resultdir)
        motifs,infos=load_motifs_typeE(memefile)
        #for m in range(len(motifs)): #only print the top 25 
        for m in motifs:
            print "motif",m
            PWM=motifs[m]
            name=infos[m][0]
            name=name.replace("MOTIF",'').strip()

            if "m-motif" in name:
                name=nametag+"_"+name.split("_")[0]+"_m-motif"
            else:
                name=nametag+"_"+name.split("_")[0]

            name = m #use motif name, NO TAG

            target=open("tmp.weblogo.faa",'w')
            for i in range(3000):
                kmer=randomkmer_typeE(PWM)
                target.write(">"+str(i)+'\n')
                target.write(kmer+"\n")
            target.close()
            #os.system("weblogo --format png --size large --title 'Motif "+str(m)+"' -a ACGTE  --color '#CC0000' A 'a' --color '#008000' T 'Py' --color '#0000CC' C 'c' -C '#CC00CC' E 'e' -C '#FFB300' G 'g' <tmp.faa > " +resultdir+"/Motif"+str(m)+".png")
            os.system("weblogo --format eps --size large -a ACGTE  --color '#CC0000' A 'a' --color '#008000' T 'Py' --color '#0000CC' C 'c' -C '#CC00CC' E 'e' -C '#FFB300' G 'g' <tmp.weblogo.faa> " +resultdir+"/"+name+".eps")
            #os.system("weblogo --format eps --size large -a ACGTE  --color '#CC0000' A 'a' --color '#008000' T 'Py' --color '#0000CC' C 'c' -C '#CC00CC' E 'e' -C '#FFB300' G 'g' -D transfac <tmp.weblogo.transfact.dat> " +resultdir+"/"+name+".eps")
        os.system('rm tmp.weblogo.faa')
    return
if __name__=="__main__":
	main()

