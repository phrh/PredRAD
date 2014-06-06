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

For the shell script (change execute permissions)

----------------
#### Usage




- **pattern_for_genome_paper.sh**.  This shell script will search all the restriction sites from the input file (patternfilename) in every genome from the input file (genomefilename). As a result the script provides the following files:
 
	- ALL.aligned.txt, ALL.failed.txt, ALL.processed.txt,  ALL.suppressed.txt - each file with a table summarizing bowtie output(reads aligned, failed, processed and suppressed) for each genome.
	- ALL.count.txt - contains a table with the number of restriciton sites found in each genome
	- ALL.size.txt - contains a table with the size of each genome

	The input arguments are: 
	- genomefilename: name of file with table with two columns (1) species code and (2) link to whole genome fasta file 
    (see  [test/genomeFileExample.txt](https://github.com/phrh/Genome-wide-predictability-of-restriction-sites-across-the-eukaryotic-tree-of-life/blob/master/test/genomeFileExample.txt))
	- patternfilename - name of file with table with tow columns (1) restriction site regular expression and (2) restriction site name 
    (see [test/Patterns_list.txt](https://github.com/phrh/Genome-wide-predictability-of-restriction-sites-across-the-eukaryotic-tree-of-life/blob/master/test/Patterns_list.txt))

	To run, just write on shell

	_./patternforgenome.sh genomefilename patternfilename_

----------------
- **genome_nucleotide_distrib_paper.py**. This python script obtains the nucleotides, dinucleotide and trinucleotides distribution for each genome from the input file (genomefilename)


	 The input arguments are:

	- genomefilename: name of file with table with two columns (1) species code and (2) link to whole genome fasta file 
    (see  [test/genomeFileExample.txt](https://github.com/phrh/Genome-wide-predictability-of-restriction-sites-across-the-eukaryotic-tree-of-life/blob/master/test/genomeFileExample.txt))
	- resultsfile : name of the outputfile 

	To run, just write on shell

	_python genome_nucleotide_distrib_paper.py genomefilename resultsfile_

----------------

- **sequence_probability_paper.py**. This python script obtains the probability for each restriction site from the input file (patternfilename) in every genome considering nt, dint and trint frequencies (distributionfile). As a result the script provides the following files:

	- $distributionfile$_nt    - contains a table with the sequences probabilities (based on nucleotide probabilities)
	- $distributionfile$_dint  - contains a table with the sequences probabilities (based on dinucleotides probabilities)
	- $distributionfile$_trint - contains a table with the sequences probabilities (based on trinucleotides probabilities)

	The input arguments are:
	- distributionfile - output from genome_nucleotide_distrib_paper (see [test/DistributionFile.txt](https://github.com/phrh/Genome-wide-predictability-of-restriction-sites-across-the-eukaryotic-tree-of-life/blob/master/test/DistributionFile.txt))
	- patternfilename - name of file with table with tow columns (1) restriction site regular expression and (2) restriction site name 
    (see [test/Patterns_list.txt](https://github.com/phrh/Genome-wide-predictability-of-restriction-sites-across-the-eukaryotic-tree-of-life/blob/master/test/Patterns_list.txt))

	To run, just write on shell
    
	_python sequence_probability_paper.py distributionfile patternsfile_

