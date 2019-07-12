# mEpigram  

This tool allows users to find methylated motifs in CpG context and also motifs that contain modified bases. 

Version 0.07

Contact: vqngo@ucsd.edu

## Installation

- The programs were written in [Python2.7](https://www.python.org/downloads/) and [Julia](https://julialang.org/downloads/) (tested with Julia v0.6.2 )

- mEpigram requires a graph of possible k-mer interactions to function. Download the graphs here: 
	* CpG mode: http://wanglab.ucsd.edu/star/mEpigram/graphE-8mer.tar.gz
	* non-CpG mode: http://wanglab.ucsd.edu/star/mEpigram/graphEF-7mer.tar.gz

- If you want to generate motif logos (Optional), please install [WebLOGO](https://pypi.python.org/pypi/weblogo) on your computer. 



## Usage
 
To run mEpigram: Use mepigram_wrapper.py script.<br /> Also, the package includes several tools to help preprocess your data before using the program.

#### mEpigram Pipeline: 
You can use the included pipeline to run mepigram, it will output motifs and their enrichement scores.

To very quickly test the pipeline, go to the program main directory and execute:

TypeE:	

`python mepigram_wrapper.py testfiles/data_typeE/test.typeE.faa testfiles/data_typeE/background_typeE-5.tsv testfiles/data_typeE/graphE-5mer/ typeE -o test.typeE`

TypeEF:

`python mepigram_wrapper.py testfiles/data_typeEF/test.typeEF.faa testfiles/data_typeEF/background_typeEF-5.tsv testfiles/data_typeEF/graphEF-5mer/ typeEF -o test.typeEF`


If you use k=8 by inputting background_typeE-8.tsv (you need to generate this), graphE-8mer (download it from our website) instead, you should be able to find several highly enriched m-motifs. 

*Note: This pipeline must be executed in the mepigram main directory. For more information, execute: *python mepigram_wrapper.py -h*


#### Work flow: preparing your own data for the program

1. Insert methylation information into the genome, the input is assumed to be in BED format by default. WIG format can be used with --wig. In BED format, each line contains chromosome name, start location (0-based index), start location +1. An output directory will be created to contain the new genome with methylation information. The reference genome should be in a directory format, with each chromosomal sequence contained in a separate file, labeled by its chromosome name. 
	
	`python modifyReference.py -f testfiles/input_modified_base/test.mCInput.bedgraph -r testfiles/samplegenome/ -t 0.5 -o meth_genomE`
	
	OR
	
	`python modifyReference.py --typeEF -f testfiles/input_modified_base/test.mCInput.bedgraph -r testfiles/samplegenome/ -t 0.5 -o meth_genomeEF`

2. Make methylated sequences from bed files and the genome above:
	
	`python bedToFasta.py -f input.bed -r methyl_ref_genomeA -o output.faa`

3. Make background model: Count the number of k-mers in the genome. This might take a while (a few hours on human whole genome) but you only need to do this once per reference genome. Example: Count the number of 8-mers in the sample genome (it's recommended to use at least k=7 for typeEF mode and k=8 for typeE mode:
	
	`python bgModel.py -gd testfiles/samplegenome/ -k 8 -m typeE`
	
	OR 
	
	`python bgModel.py -gd testfiles/samplegenome/ -k 8 -m typeEF`
	
4. Running the pipeline: After getting all the above steps done, you can now use mepigram_wrapper.py on your own data. Enjoy!


## Optional steps:
#### Motif scanning: 
After discovering motifs from the pipeline, you might want to find where they are in a set of sequences. To identify locations of matches with your motifs, you can use the motif scanning tool. The program takes a FASTA file, a motif PWM file, and a background file (optional) that specifies background base composition. The program will use the nucleotide composition from the input sequences if a background is not provided.
	
`python motifscannerA.py -f testfiles/data_typeE/test.typeE.faa -m testfiles/data_typeE/test.typeE.meme -o test.scanned.txt`

Example of the content of a background nucleotide distribution file:

A	0.295484488819<br/>
C	0.199430583447<br/>
E	0.00508492759333<br/>
T	0.295484488959<br/>
G	0.204515511181<br/>
	
*Note: The motifscannerA.py file should be executed in the main mEpigram directory

#### Making motif LOGOs: 
You can generate both types of logo using the `makeLOGO.py` script. Use the flag --typeEF to generate typeEF motifs. Once you've installed Weblogo, you can run this command to test the program:

`python makeLOGO.py -m testfiles/data_typeE/test.typeE.meme -o testLogo`

#### Data used in the manuscript:
[H1_hg19_methylated_genome_typeE](http://tabit.ucsd.edu/mepigram/hg19_data/hg19_H1_typeE.tar.gz) <br />
[GM12878_hg19_methylated_genome_typeE](http://tabit.ucsd.edu/mepigram/hg19_data/hg19_GM12878_typeE.tar.gz)


