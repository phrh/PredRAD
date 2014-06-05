### Genome-wide predictability of restriction sites across the eukaryotic tree of life

The software developed is available in this page. This software will help guide the design of any study using RAD sequencing and related methods.

----------------
#### Requirements

- Python 2.7 and above
- [Biopython](http://biopython.org/wiki/Main_Page)
- [Bowtie](https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.0.1)

----------------
#### Install

Download python and shell scritps 


----------------
#### Usage




- **pattern_for_genome.sh**.  This shell script will search all the restriction sites from the input file (patternfilename) in every genome from the input file (genomefilename). As a result this scripts creates  

	The input arguments are are: 
	- genomefilename: name of file with table with two columns (1) species code and (2) link to whole genome fasta file 
    (see example/test_pattern_for_genome/example_genome_table )
    
	- patternfilename - name of file with table with tow columns (1) restriction site regular expression and (2) restriction site name 
    (see example/test_pattern_for_genome/patternfilename)

	To run, just write on shell

	_patternforgenome.sh genomefilename patternfilename_

- **genome_nucleotide_distrib.py**. This python module will 


	 The input files are:

	- genomefilename 
	- resultsfile 

	To run, just write on shell

	_python genomenucleotidedistrib.py genomefilename resultsfile_

- **sequence_probability.py**

	To run, just write on shell
    
	_python sequenceprobability.py distributionfile patternsfile_

