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

* **pattern_for_genome.sh**.  This script will search all the restriction sites from the input file (patternfilename) in every genome from the input file (genomefilename) 

> The input arguments are are: 

  .1. genomefilename  - name of file with table with two columns (1) species code and (2) link to whole genome fasta file 
    (see example/test_pattern_for_genome/example_genome_table )
  .2. patternfilename - name of file with table with tow columns (1) restriction site regular expression and (2) restriction site name 
    (see example/test_pattern_for_genome/patternfilename)

> To run, just write on shell

  pattern_for_genome.sh genomefilename patternfilename

 

2. **genome_nucleotide_distrib.py**.


> The input files are:

  1. genomefilename - 
  2. resultsfile 

* python genome_nucleotide_distrib.py genomefilename resultsfile 

3. **sequence_probability.py**

* python sequence_probability.py distributionfile patternsfile 

