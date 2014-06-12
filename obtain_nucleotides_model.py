# PredRAD - guides the design of any study using RAD sequencing and related methods.
#
# Created by Santiago Herrera and Paula H. Reyes-Herrera on 11 June 2014
# Copyright (c) 2014 Santiago Herrera and Paula H. Reyes-Herrera. All rights reserved.
# 
# This file is part of PredRAD.
#
# PredRAD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2.
#-------------------------------------------------------------------------------------
#!/usr/bin/env python
# Usage: python obtain_nucleotides_model.py genomefilename resultsfile

#-------------------------------FUNCTIONS-LIST----------------------------------------
# writelog(filename,text)
# process_inputfile(url, scode, logfile)
# count_nt_freq(scode,filenameall)
# int_main(genomefile,resultfile)

#---------------------------------FUNCTIONS-------------------------------------------

def writelog(filename,text):
	import os
	from datetime import date
	import getpass
	f = open(filename, 'a+')
	now = date.today().ctime()
	myhost = os.uname()[1]
	user = getpass.getuser()
	f.write(str(now)+' '+user+'@'+myhost+'\t'+text+'\n')
	f.close()

def process_inputfile(url, scode, logfile):
#Importing libraries
	import os, gzip, StringIO, urllib2	
	i=1
#Obtain files from web site and concatenate if multiple
	if url.count('wgs')>0 :
		commandwget='wget '+url+'.'+str(i)+'.fsa_nt.gz -O '+scode+'.fsa.'+str(i)+'.gz'
		print commandwget
		os.system(commandwget)
		size = os.path.getsize(scode+'.fsa.'+str(i)+'.gz')
		if size < 100 :
			writelog(logfile,'fasta '+scode+' size '+str(size))
			commandrm='rm '+scode+'.fsa.'+str(i)+'.gz'
			os.system(commandrm)
		if size > 100 :
			while size > 100 :
				commandgunzip='gunzip '+scode+'.fsa.'+str(i)+'.gz'
				os.system(commandgunzip)
				commandrm='rm '+scode+'.fsa.'+str(i)+'.gz'
				os.system(commandrm)
				i=i+1
				commandwget='wget '+url+'.'+str(i)+'.gz -O '+scode+'.fsa.'+str(i)+'.gz'
				os.system(commandwget)
				size = os.path.getsize(scode+'.fsa.'+str(i)+'.gz')
			commandcat='cat '+scode+'.fsa.* > '+scode+'.fasta'
			os.system(commandcat)
			commandrm='rm '+scode+'.fsa.*'
			os.system(commandrm)
	
	elif url.count('assembled_chromosomes')>0 :
		commandwget='wget '+url+str(i)+'.fa.gz -O '+scode+'.fsa.'+str(i)+'.gz'
		print commandwget
		os.system(commandwget)
		size = os.path.getsize(scode+'.fsa.'+str(i)+'.gz')
		if size < 100 :
			writelog(logfile,'fasta '+scode+' size '+str(size))
			commandrm='rm '+scode+'.fsa.'+str(i)+'.gz'
			os.system(commandrm)
		if size > 100 :
			while size > 100 :
				commandgunzip='gunzip '+scode+'.fsa.'+str(i)+'.gz'
				os.system(commandgunzip)
				commandrm='rm '+scode+'.fsa.'+str(i)+'.gz'
				os.system(commandrm)
				i=i+1
				commandwget='wget '+url+str(i)+'.fa.gz -O '+scode+'.fsa.'+str(i)+'.gz'
				os.system(commandwget)
				commandgunzip='gunzip '+scode+'.fsa.'+str(i)+'.gz'
				os.system(commandgunzip)
				commandrm='rm '+scode+'.fsa.'+str(i)+'.gz'
				os.system(commandrm)
				i=i+1
				commandwget='wget '+url+str(i)+'.fa.gz -O '+scode+'.fsa.'+str(i)+'.gz'
				os.system(commandwget)
				size = os.path.getsize(scode+'.fsa.'+str(i)+'.gz')
			commandcat='cat '+scode+'.fsa.* > '+scode+'.fasta'
			os.system(commandcat)
			commandrm='rm '+scode+'.fsa.*'
			os.system(commandrm)
		
	else :
		if url.count('gz')>0 :		
			commandwget='wget '+url+' -O '+scode+'.fasta.gz'
			os.system(commandwget)
			commandgunzip='gunzip '+scode+'.fasta.gz'
			os.system(commandgunzip)
		else :
			commandwget='wget '+url+' -O '+scode+'.fasta'
			os.system(commandwget)
	

def count_nt_freq(scode,filenameall):
#Importing libraries	
	from Bio import SeqIO
	import re
	import os
#handling sequences from fasta file
	handle = open(scode+'.fasta','rU')
	countA=0
	countT=0
	countG=0
	countC=0
	countN=0
	countAA=0
	countAG=0
	countAC=0
	countAT=0
	countGA=0
	countGG=0
	countGC=0
	countGT=0
	countCA=0
	countCG=0
	countCC=0
	countCT=0
	countTA=0
	countTG=0
	countTC=0
	countTT=0
	countAAA=0
	countAAG=0
	countAAC=0
	countAAT=0
	countAGA=0
	countAGG=0
	countAGC=0
	countAGT=0
	countACA=0
	countACG=0
	countACC=0
	countACT=0
	countATA=0
	countATG=0
	countATC=0
	countATT=0
	countGAA=0
	countGAG=0
	countGAC=0
	countGAT=0
	countGGA=0
	countGGG=0
	countGGC=0
	countGGT=0
	countGCA=0
	countGCG=0
	countGCC=0
	countGCT=0
	countGTA=0
	countGTG=0
	countGTC=0
	countGTT=0
	countCAA=0
	countCAG=0
	countCAC=0
	countCAT=0
	countCGA=0
	countCGG=0
	countCGC=0
	countCGT=0
	countCCA=0
	countCCG=0
	countCCC=0
	countCCT=0
	countCTA=0
	countCTG=0
	countCTC=0
	countCTT=0
	countTAA=0
	countTAG=0
	countTAC=0
	countTAT=0
	countTGA=0
	countTGG=0
	countTGC=0
	countTGT=0
	countTCA=0
	countTCG=0
	countTCC=0
	countTCT=0
	countTTA=0
	countTTG=0
	countTTC=0
	countTTT=0
	countATCG=0
	countnt=0
#counting for each record
	for record in SeqIO.parse(handle,'fasta') :
	#Regular expressions
		countnt=countnt+len(re.findall('[AGCTKMRYSWBVHDN]',str(record.seq.upper())))
		countATCG=countATCG+len(re.findall('[AGCT]',str(record.seq.upper())))		
		countA=countA+record.seq.upper().count('A')	
		countG=countG+record.seq.upper().count('G')
		countC=countC+record.seq.upper().count('C')
		countT=countT+record.seq.upper().count('T')
		countN=countN+record.seq.upper().count('N')
		countAA=countAA+record.seq.upper().count('AA')
		countAG=countAG+record.seq.upper().count('AG')
		countAC=countAC+record.seq.upper().count('AC')
		countAT=countAT+record.seq.upper().count('AT')
		countGA=countGA+record.seq.upper().count('GA')
		countGG=countGG+record.seq.upper().count('GG')
		countGC=countGC+record.seq.upper().count('GC')
		countGT=countGT+record.seq.upper().count('GT')
		countCA=countCA+record.seq.upper().count('CA')
		countCG=countCG+record.seq.upper().count('CG')
		countCC=countCC+record.seq.upper().count('CC')
		countCT=countCT+record.seq.upper().count('CT')
		countTA=countTA+record.seq.upper().count('TA')
		countTG=countTG+record.seq.upper().count('TG')
		countTC=countTC+record.seq.upper().count('TC')
		countTT=countTT+record.seq.upper().count('TT')
		countAAA=countAAA+record.seq.upper().count('AAA')
		countAAG=countAAG+record.seq.upper().count('AAG')
		countAAC=countAAC+record.seq.upper().count('AAC')
		countAAT=countAAT+record.seq.upper().count('AAT')
		countAGA=countAGA+record.seq.upper().count('AGA')
		countAGG=countAGG+record.seq.upper().count('AGG')
		countAGC=countAGC+record.seq.upper().count('AGC')
		countAGT=countAGT+record.seq.upper().count('AGT')
		countACA=countACA+record.seq.upper().count('ACA')
		countACG=countACG+record.seq.upper().count('ACG')
		countACC=countACC+record.seq.upper().count('ACC')
		countACT=countACT+record.seq.upper().count('ACT')
		countATA=countATA+record.seq.upper().count('ATA')
		countATG=countATG+record.seq.upper().count('ATG')
		countATC=countATC+record.seq.upper().count('ATC')
		countATT=countATT+record.seq.upper().count('ATT')
		countGAA=countGAA+record.seq.upper().count('GAA')
		countGAG=countGAG+record.seq.upper().count('GAG')
		countGAC=countGAC+record.seq.upper().count('GAC')
		countGAT=countGAT+record.seq.upper().count('GAT')
		countGGA=countGGA+record.seq.upper().count('GGA')
		countGGG=countGGG+record.seq.upper().count('GGG')
		countGGC=countGGC+record.seq.upper().count('GGC')
		countGGT=countGGT+record.seq.upper().count('GGT')
		countGCA=countGCA+record.seq.upper().count('GCA')
		countGCG=countGCG+record.seq.upper().count('GCG')
		countGCC=countGCC+record.seq.upper().count('GCC')
		countGCT=countGCT+record.seq.upper().count('GCT')
		countGTA=countGTA+record.seq.upper().count('GTA')
		countGTG=countGTG+record.seq.upper().count('GTG')
		countGTC=countGTC+record.seq.upper().count('GTC')
		countGTT=countGTT+record.seq.upper().count('GTT')
		countCAA=countCAA+record.seq.upper().count('CAA')
		countCAG=countCAG+record.seq.upper().count('CAG')
		countCAC=countCAC+record.seq.upper().count('CAC')
		countCAT=countCAT+record.seq.upper().count('CAT')
		countCGA=countCGA+record.seq.upper().count('CGA')
		countCGG=countCGG+record.seq.upper().count('CGG')
		countCGC=countCGC+record.seq.upper().count('CGC')
		countCGT=countCGT+record.seq.upper().count('CGT')
		countCCA=countCCA+record.seq.upper().count('CCA')
		countCCG=countCCG+record.seq.upper().count('CCG')
		countCCC=countCCC+record.seq.upper().count('CCC')
		countCCT=countCCT+record.seq.upper().count('CCT')
		countCTA=countCTA+record.seq.upper().count('CTA')
		countCTG=countCTG+record.seq.upper().count('CTG')
		countCTC=countCTC+record.seq.upper().count('CTC')
		countCTT=countCTT+record.seq.upper().count('CTT')
		countTAA=countTAA+record.seq.upper().count('TAA')
		countTAG=countTAG+record.seq.upper().count('TAG')
		countTAC=countTAC+record.seq.upper().count('TAC')
		countTAT=countTAT+record.seq.upper().count('TAT')
		countTGA=countTGA+record.seq.upper().count('TGA')
		countTGG=countTGG+record.seq.upper().count('TGG')
		countTGC=countTGC+record.seq.upper().count('TGC')
		countTGT=countTGT+record.seq.upper().count('TGT')
		countTCA=countTCA+record.seq.upper().count('TCA')
		countTCG=countTCG+record.seq.upper().count('TCG')
		countTCC=countTCC+record.seq.upper().count('TCC')
		countTCT=countTCT+record.seq.upper().count('TCT')
		countTTA=countTTA+record.seq.upper().count('TTA')
		countTTG=countTTG+record.seq.upper().count('TTG')
		countTTC=countTTC+record.seq.upper().count('TTC')
		countTTT=countTTT+record.seq.upper().count('TTT')
	handle.close()
	os.remove(scode+'.fasta')
#printing into file
	f = open(filenameall,'a')
	text=scode+'\t'+str(countnt)+'\t'+str(countATCG)
	text=text+'\t'+str(countA)+'\t'+str(countG)+'\t'+str(countC)+'\t'+str(countT)+'\t'+ \
	str(countN)+'\t'+str(countAA)+'\t'+str(countAG)+'\t'+str(countAC)+'\t'+str(countAT)+ \
	'\t'+str(countGA)+'\t'+str(countGG)+'\t'+str(countGC)+'\t'+str(countGT)+'\t'+str(countCA)+ \
	'\t'+str(countCG)+'\t'+str(countCC)+'\t'+str(countCT)+'\t'+str(countTA)+'\t'+str(countTG)+ \
	'\t'+str(countTC)+'\t'+str(countTT)+'\t'+str(countAAA)+'\t'+str(countAAG)+'\t'+str(countAAC)+ \
	'\t'+str(countAAT)+'\t'+str(countAGA)+'\t'+str(countAGG)+'\t'+str(countAGC)+'\t'+str(countAGT)+ \
	'\t'+str(countACA)+'\t'+str(countACG)+'\t'+str(countACC)+'\t'+str(countACT)+'\t'+str(countATA)+ \
	'\t'+str(countATG)+'\t'+str(countATC)+'\t'+str(countATT)+'\t'+str(countGAA)+'\t'+str(countGAG)+ \
	'\t'+str(countGAC)+'\t'+str(countGAT)+'\t'+str(countGGA)+'\t'+str(countGGG)+'\t'+str(countGGC)+ \
	'\t'+str(countGGT)+'\t'+str(countGCA)+'\t'+str(countGCG)+'\t'+str(countGCC)+'\t'+str(countGCT)+ \
	'\t'+str(countGTA)+'\t'+str(countGTG)+'\t'+str(countGTC)+'\t'+str(countGTT)+'\t'+str(countCAA)+ \
	'\t'+str(countCAG)+'\t'+str(countCAC)+'\t'+str(countCAT)+'\t'+str(countCGA)+'\t'+str(countCGG)+ \
	'\t'+str(countCGC)+'\t'+str(countCGT)+'\t'+str(countCCA)+'\t'+str(countCCG)+'\t'+str(countCCC)+ \
	'\t'+str(countCCT)+'\t'+str(countCTA)+'\t'+str(countCTG)+'\t'+str(countCTC)+'\t'+str(countCTT)+ \
	'\t'+str(countTAA)+'\t'+str(countTAG)+'\t'+str(countTAC)+'\t'+str(countTAT)+'\t'+str(countTGA)+ \
	'\t'+str(countTGG)+'\t'+str(countTGC)+'\t'+str(countTGT)+'\t'+str(countTCA)+'\t'+str(countTCG)+ \
	'\t'+str(countTCC)+'\t'+str(countTCT)+'\t'+str(countTTA)+'\t'+str(countTTG)+'\t'+str(countTTC)+ \
	'\t'+str(countTTT)+'\n'
	f.write(text)


#MAIN PROGRAM
def int_main(genomefile,resultfile):
	import os
	from datetime import date
	now = date.today()	
	logfile='genome_nucleotide_distribution_'+str(now)
	writelog(logfile,"--------------------")
	writelog(logfile,"Input file:"+genomefile+' '+resultfile) 
	writelog(resultfile,'Species\tnt\tAGCT\tA\tG\tC\tT\tN\tAA\tAG\tAC\tAT\tGA\tGG\
	GC\tGT\tCA\tCG\tCC\tCT\tTA\tTG\tTC\tTT\tAAA\tAAG\tAAC\tAAT\tAGA\tAGG\tAGC\tAGT\tACA\
	ACG\tACC\tACT\tATA\tATG\tATC\tATT\tGAA\tGAG\tGAC\tGAT\tGGA\tGGG\tGGC\tGGT\tGCA\tGCG\
	GCC\tGCT\tGTA\tGTG\tGTC\tGTT\tCAA\tCAG\tCAC\tCAT\tCGA\tCGG\tCGC\tCGT\tCCA\tCCG\tCCC\
	CCT\tCTA\tCTG\tCTC\tCTT\tTAA\tTAG\tTAC\tTAT\tTGA\tTGG\tTGC\tTGT\tTCA\tTCG\tTCC\tTCT\
	TTA\tTTG\tTTC\tTTT')
	g = open(genomefile,'r')
	for line in g :
		scode, url = line.split("\t")
		url=url.rstrip("\n")
		print 'scode '+scode
		print 'url' +url
		writelog(logfile,scode+" process starts")
		process_inputfile(url,scode,logfile)
		count_nt_freq(scode,resultfile)

	
#----------------------------------MAIN-------------------------------------------

import sys
int_main(sys.argv[1],sys.argv[2])


