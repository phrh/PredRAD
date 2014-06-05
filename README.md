Genome-wide predictability of restriction sites across the eukaryotic tree of life
-----------------------------------------------------------------------------------
The software developed is availabe in this page.  will help guide the design of any study using RAD sequencing and related methods.


Requirements
------------
Python 2.7 and above


Install
-------
Download python scritps 



Usage
-----

The first script is pattern_for_genome.sh. 
The input files are: 


* pattern_for_genome.sh genomefilename patternfilename

  1. genomefilename  - table with species code and link to whole genome fasta file (see example/test_pattern_for_genome/example_genome_table )
  2. patternfilename - table with restriction site regular expression and restriction site name (see example/test_pattern_for_genome/patternfilename)


* python genome_nucleotide_distrib.py genomefilename resultsfile 



* python sequence_probability.py distributionfile patternsfile 

