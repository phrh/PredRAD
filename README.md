Genome-wide predictability of restriction sites across the eukaryotic tree of life
-----------------------------------------------------------------------------------
The software developed is available in this page. This software will help guide the design of any study using RAD sequencing and related methods.


Requirements
------------
Python 2.7 and above


Install
-------
Download python scritps 



Usage
-----

1.  The first script is **pattern_for_genome.sh**.  The input files are: 


* pattern_for_genome.sh genomefilename patternfilename

  1. genomefilename  - table with two columns (1) species code and (2) link to whole genome fasta file 
    (see example/test_pattern_for_genome/example_genome_table )
  2. patternfilename - table with tow columns (1) restriction site regular expression and (2) restriction site name 
    (see example/test_pattern_for_genome/patternfilename)


* python genome_nucleotide_distrib.py genomefilename resultsfile 

  1. genomefilename - 
  2. results file

* python sequence_probability.py distributionfile patternsfile 

