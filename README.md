## PredRAD

High-throughput sequencing of reduced representation libraries obtained through digestion with restriction enzymes–generally known as restriction-site associated DNA sequencing (RAD-seq)–is now one most commonly used strategies to generate single nucleotide polymorphism data in eukaryotes. The choice of restriction enzyme is critical for the design of any RAD-seq study as it determines the number of genetic markers that can be obtained for a given species, and ultimately the success of a project. 

For the design of a study using RAD-seq, or a related methodology, there are two general fundamental questions that researchers face: i) what is the best restriction enzyme to use to obtain a desired number of RAD tags in the organism of interest? And ii) how many markers can be obtained with a particular enzyme in the organism of interest? Thi software pipeline will allow any researcher to obtain an approximate answer to these questions and this will help guide the design of any study using RAD sequencing and related methods.

This Git contains the software code and output results from [Herrera S., P.H. Reyes-Herrera & T.M. Shank (2014) Genome-wide predictability of restriction sites across the eukaryotic tree of life. bioRxiv preprint doi: http://dx.doi.org/10.1101/007781](biorxiv.org/content/early/2014/08/08/007781)

----------------
#### Requirements

- Python 2.7 and above
- [Biopython](http://biopython.org/wiki/Main_Page)
- [Bowtie](https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.0.1)

----------------
#### Install

Download python and shell scritps 

For the shell script (change execute permissions using chmod u+x)

----------------
#### Usage




- **restriction_site_search.sh**.  This shell script will search all the restriction sites from the input file (patternfilename) in every genome from the input file (genomefilename). As a result the script provides the following files:
 
	- ALL.aligned.txt, ALL.failed.txt, ALL.processed.txt,  ALL.suppressed.txt - each file with a table summarizing bowtie output(reads aligned, failed, processed and suppressed) for each genome.
	- ALL.count.txt - contains a table with the number of restriciton sites found in each genome
	- ALL.size.txt - contains a table with the size of each genome

	The input arguments are: 
	- genomefilename: name of file with table with two columns (1) species code and (2) link to whole genome fasta file 
    (see  [test/genomeFileExample.txt](https://github.com/phrh/Genome-wide-predictability-of-restriction-sites-across-the-eukaryotic-tree-of-life/blob/master/test/genomeFileExample.txt))
	- patternfilename - name of file with table with tow columns (1) restriction site regular expression and (2) restriction site name 
    (see [test/Patterns_list.txt](https://github.com/phrh/Genome-wide-predictability-of-restriction-sites-across-the-eukaryotic-tree-of-life/blob/master/test/Patterns_list.txt))

	To run, just write on shell

	_./restriction_site_search.sh genomefilename patternfilename_

----------------
- **obtain_nucleotides_model.py**. This python script obtains the nucleotides, dinucleotide and trinucleotides distribution for each genome from the input file (genomefilename)


	 The input arguments are:

	- genomefilename: name of file with table with two columns (1) species code and (2) link to whole genome fasta file 
    (see  [test/genomeFileExample.txt](https://github.com/phrh/Genome-wide-predictability-of-restriction-sites-across-the-eukaryotic-tree-of-life/blob/master/test/genomeFileExample.txt))
	- resultsfile : name of the outputfile 

	To run, just write on shell

	_python obtain_nucleotides_model.py genomefilename resultsfile_

----------------

- **sequence_probability.py**. This python script obtains the probability for each restriction site from the input file (patternfilename) in every genome considering nt, dint and trint frequencies (distributionfile). As a result the script provides the following files:

	- $distributionfile$_nt    - contains a table with the sequences probabilities (based on nucleotide probabilities)
	- $distributionfile$_dint  - contains a table with the sequences probabilities (based on dinucleotides probabilities)
	- $distributionfile$_trint - contains a table with the sequences probabilities (based on trinucleotides probabilities)

	The input arguments are:
	- distributionfile - output from genome_nucleotide_distrib_paper (see [test/DistributionFile.txt](https://github.com/phrh/Genome-wide-predictability-of-restriction-sites-across-the-eukaryotic-tree-of-life/blob/master/test/DistributionFile.txt))
	- patternfilename - name of file with table with tow columns (1) restriction site regular expression and (2) restriction site name 
    (see [test/Patterns_list.txt](https://github.com/phrh/Genome-wide-predictability-of-restriction-sites-across-the-eukaryotic-tree-of-life/blob/master/test/Patterns_list.txt))

	To run, just write on shell
    
	_python sequence_probability.py distributionfile patternsfile_


#### License

Created by Santiago Herrera and Paula H. Reyes-Herrera on 11 June 2014
Copyright (c) 2014 Santiago Herrera and Paula H. Reyes-Herrera. All rights reserved.
 
PredRAD is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 2.
