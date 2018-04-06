# mEpigram  

This tool allows users to find methylated motifs in CpG context and also motifs that contain modified bases. 

Version 0.07

Contact: vqngo@ucsd.edu

## Installation

- The programs were written in [Python2.7](https://www.python.org/downloads/) and [Julia](https://julialang.org/downloads/) (tested with Julia v0.6.2 )

- mEpigram requires a graph of possible k-mer interactions to function. Download the graphs here: 
	* CpG mode: http://wanglab.ucsd.edu/star/mepigram/graphE-8mer.tar.gz
	* non-CpG mode: http://wanglab.ucsd.edu/star/mepigram/graphEF-7mer.tar.gz

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

*Note: This pipeline must be executed in the mepigram main directory.


#### mEpigram preprocessing scripts:

1. Insert methylation information into the genome, the input is assumed to be in BED format by default. WIG format can be used with --wig. In BED format, each line contains chromosome name, start location (0-based index), start location +1. An output directory will be created to contain the new genome with methylation information. The reference genome should be in a directory format, with each chromosomal sequence contained in a separate file, labeled by its chromosome name. 
	
	`python modifyReference.py -f input.bed -r reference_genome_directory -o methyl_ref_genomeA`
	
	OR
	
	`python modifyReference.py --typeEF -f input.bed -r reference_genome_directory -o methyl_ref_genomeA`

2. Make methylated sequences from bed files and the genome above:
	
	`python bedToFasta.py -f input.bed -r methyl_ref_genomeA -o output.faa`

3. Make background model: Count the number of k-mers in the genome. This might take a while (a few hours on human whole genome) but you only need to do this once per reference genome. Example: Count the number of 5-mers in the sample genome:
	
	`python bgModel.py -gd testfiles/samplegenome/ -k 5 -m typeE`
	
	OR 
	
	`python bgModel.py -gd testfiles/samplegenome/ -k 5 -m typeEF`

#### Motif scanning: 
To identify locations of matches using your motifs, you can use the motif scanning tool. The program takes a FASTA file, a motif PWM file, and a background file that states background base composition. Although the background is optional, it is recommended that you use the appropriate background as the program will assume equal nucleotide distribution if the background is not provided.
	
`python motifscannerA.py [options] -f fastafile -m motiffile -o output_file -b backgroundBaseComposition`
	
*Note: The motifscannerA.py file should be executed in the main mEpigram directory

#### Making motif LOGOs: 
You can generate both types of logo using the `makeLOGO.py` script. Use the flag --typeEF to generate typeEF motifs.

