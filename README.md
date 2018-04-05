# mEpigram  

This tool allows users to find methylated motifs in CpG context and also motifs that contain modified bases. 

Version 0.05

Contact: vqngo@ucsd.edu

## Installation

- The programs were written in Python2.7 and Julia (tested with Julia v0.6.2 )

- mEpigram requires a graph of possible k-mer interactions to function. Download the graphs here: 
	* CpG mode: http://wanglab.ucsd.edu/star/mepigram/graphE-8mer.tar.gz
	* non-CpG mode: http://wanglab.ucsd.edu/star/mepigram/graphEF-7mer.tar.gz

- Although the Julia scripts are not required for motif discovery, it's recommended that you installed Julia to use motif scanning and enrichment calculation.

- If you want to generate motif logos, please install WebLOGO on your computer https://pypi.python.org/pypi/weblogo  



## Usage
 
To run mEpigram: It's recommended to run mEpigram using the mepigram_wrapper.py script. However, you can run each modules separately as well.

#### mEpigram Pipeline: 
You can use the included pipeline to run mepigram, it will perform dinucleotide-shulffing, mepigram_typeE, and enrichment calculation (which uses the motif scanning module).

To very quickly test the pipeline: 

TypeE:	

`python mepigram_wrapper.py -f testfiles/test_data_typeE/ENCFF002CQR.narrowPeak.1000.faa -m typeE -b testfiles/test_data_typeE/background_typeE-5.tsv -g testfiles/test_data_typeE/graphE-5mer/`

TypeEF:

`python mepigram_wrapper.py -f testfiles/test_data_typeEF/ENCFF002CQR.narrowPeak.1000.faa -m typeEF -b testfiles/test_data_typeEF/background_typeEF-5.tsv -g testfiles/test_data_typeEF/graphEF-5mer/`


If you use k=8 by inputting background_typeE-8.tsv, graphE-8mer (download it from our website) instead, you should be able to find several highly enriched m-motifs. 

*Note: This pipeline must be executed in the mepigram main directory. It also requires Julia installed.


#### Running each module separately:

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

4. Shuffle the sequences: this step will di-nucleotide shuffle your FASTA input sequences to be used. 
	
	`python fasta-dinucleotide-shuffle_typeE.py -f sequences.faa > sequences.DS.faa`
	
	OR
	
	`python fasta-dinucleotide-shuffle_typeEF.py -f sequences.faa > sequences.DS.faa`


5. Run mEpigram (all parameters must be provided):
	
	`python mepigram_typeE.py testInput.faa testInputCTCF.DS.faa your_specific_genome/background_met-8.tsv metgraph-8mer/ resultfile.meme max_No_motifs`
	
	OR
	
	`python mepigram_typeEF.py testInput.faa testInputCTCF.DS.faa your_specific_genome/background_met-8.tsv metgraph-8mer/ resultfile.meme max_No_motifs`
	
* max_No_motifs: integer, the maximum number of motifs to be produced.

## Other tools

#### Motif scanning: 
To identify locations of matches using your motifs, you can use the motif scanning tool. The motif-scanning program is written in julia so it's necessary to have the Julia language installed. The program takes a FASTA file, a motif PWM file, and a background file that states background base composition. Although the background is optional, it is recommended that you use the appropriate background as the program will assume equal nucleotide distribution if the background is not provided.
	
`python motifscannerA.py [options] -f fastafile -m motiffile -o output_file -b backgroundBaseComposition`
	
*Note: The motifscannerA.py file should be executed in the main mEpigram directory

#### Making motif LOGOs: 
You can generate both types of logo using the `makeLOGO.py` script. Use the flag --typeEF to generate typeEF motifs.

