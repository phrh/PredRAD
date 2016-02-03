#!/bin/bash
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
# Please cite :      Santiago Herrera, Paula H. Reyes-Herrera, and Timothy M. Shank.  Predicting RAD-seq Marker Numbers across the Eukaryotic Tree of Life. Genome Biol Evol (2015) Vol. 7 3207-3225 
#-------------------------------------------------------------------------------------
# Usage:   ./restriction_site_search.sh genomefilename patternfilename localfile
# Arguments: 
# 	$1 genomefile
# 	$2 patternfile
#       $3 localfile flag (Yes - if the genome file contains paths instead of url, No - otherwise) default value NO
#	$4 bowtieflag (Yes to use bowtie to align and get bowtie statistics, No - otherwise) default value YES 
#-------------------------------FUNCTIONS-LIST----------------------------------------
# pattern_upstream_downstream
# process_inputfile
# count_nt_sites
# extract_bowtie_stats
# extract_bowtie_stats_by_genome
# extract_bowtie_stats_for_all
# check_inputfile
# writelog
#---------------------------------FUNCTIONS-------------------------------------------

function writelog(){
	echo "[$(date)] ($USER@$(hostname)) - $1" >> $2
}

function check_inputfile(){
	#FUNCTION check_inputfile
	#Arguments:
	# $1 file name 
	# $2 temporary folder name
	# This function checks for line break at end of file and add if missing and Remove blank lines
	filename=$(basename $1);
	od -c $1 | tail -2 > $2/last_characters.txt;
	cat $2/last_characters.txt |while read number character; 
		do 
		if [ "$character" != "\n" ]
		then 
			while read line || [ -n "$line" ]; 
			do echo $line; 
			done < $1 > $2/$filename.corrected;
			cp $2/$filename.corrected $1;
		fi; 
	done;
	sed -i '/^[[:space:]]*$/d' $1;
}

function process_inputfile(){
	#FUNCTION process_inputfile
	#Arguments: 
	# $1 genome file url
	# $2 species code
	# $3 Tempfolder
#if the genome does not come from ncbi site, the gbff file is not there ... write it in the log file__
	i=1
	stringa="ncbi"
	echo $url | grep $stringa >> /dev/null
	if [ $? = 0 ]
	then
	#ncbi site
		wget $1.$i.fsa_nt.gz -O $2.fsa.$i.gz
		#wget $1.$i.gz -O $2.fsa.$i.gz
		size=$(stat -c %s $2.fsa.$i.gz)
		echo "SIZE $size"
		while [ $size -gt 100 ]
		do	
			echo "Inside while"			
			if [ $size -lt 100 ]
			then 
				echo "Inside if - then"
				writelog "fasta $i size $size" 	 restriction_site_search.log
				echo "size < 100 in $i"
				rm $2.fsa.$i.gz
			else
				echo "Inside else"
				echo "gunzip $2.fsa.$i.gz"
				gunzip $2.fsa.$i.gz
				#echo "rm $2.fsa.$i.gz"
				#rm $2.fsa.$i.gz
			fi;
			echo " Increasing i -> $i"					
			i=$((i+1))		
			echo "wget $1.$i.gz -O $2.fsa.$i.gz"
			wget $1.$i.gz -O $2.fsa.$i.gz
			size=$(stat -c %s $2.fsa.$i.gz)

		done;
	writelog "fasta found: $2"  restriction_site_search.log
	cat $2.fsa.* > $2.fasta
	rm $2.fsa.*
	tr '[:lower:]' '[:upper:]' < $2.fasta > $2.UP
	rm $2.fasta
	else
			#other site
	writelog "fasta not found: $2"  restriction_site_search.log
	fi;
 
}

function count_nt_sites(){
	#FUNCTION count_nt_sites
	#Arguments: 
	# $1 species code
	# $2 patterns filename
	# $3 temp foldername
	#Count the number of nucleotides in the genome (exclude fasta sequence names), including ambiguities
	echo -e -n "$1\t" >> ALL.size.txt
	grep -v '[>]' $1.UP | grep -o '[AGCTKMRYSWBVHDN]' | wc -w > $1.txt
	#Count the number of nucleotides in the genome (exclude fasta sequence names), excluding ambiguities
	grep -v '[>]' $1.UP | grep -o '[AGCT]' | wc -w >> $1.txt
	#Count the number of GC nucleotides in the genome, #Count the number of nucleotides in the genome (exclude fasta sequence names), excluding ambiguities
	grep -v '[>]' $1.UP | grep -o '[GC]' | wc -w >> $1.txt
	while read size; do echo -n -e "$size \t" >> ALL.size.txt; done < $1.txt
	echo -e -n '\n' >> ALL.size.txt
	echo -n -e 'Counts\t'>> $1.txt
	#Count the number of cut sites per enzyme (exclude fasta sequence names)	
	while read pattern patternname; do eval "grep -o '$pattern' $1.UP | grep -v '[>|.,0123456789]' | wc -w >> $3/$1.count.txt" ; done < $2
	while read count; do echo -n -e "$count \t" >> $1.txt; done < $3/$1.count.txt 
}

function pattern_upstream_downstream(){
	#FUNCTION pattern_upstream_downstream
	#Arguments: 
	# $1 speciescode
	# $2 pattern 
	# $3 patternname
	# $4 temp foldername
	# This function search for recognition sequence patterns in the genome and "sequences" 100bp up- and down-stream of each restriction site
	input="$2";
	reverse="";
	len=${#input};
	for (( i=$len-1; i>=0; i-- ))
	do	
		if [ "${input:$i:1}" == "]" ]
		then 
			reverse="$reverse["
		else
			if [ "${input:$i:1}" == "[" ]
			then reverse="$reverse]"
			else
				reverse="$reverse${input:$i:1}"
			fi;
		fi; 
	done
	
	eval "
	tr -d '\n' < $1.UP > $4/$1_nb
	grep -o '$2[AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT]' $4/$1_nb | grep -v '[>|.,0123456789]' > $4/$1_$3_up.txt
	rev $4/$1_nb > $4/$1_rev
	grep -o '$reverse[AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT][AGCT]' $4/$1_rev | grep -v '[>|.,0123456789]' > $4/$1_$3_down_rev.txt
	rev $4/$1_$3_down_rev.txt > $4/$1_$3_down.txt
	cat $4/$1_$3_up.txt $4/$1_$3_down.txt > $4/$1_$3_all
	"
	awk '{print ">Seq"rand()"\n"$0}' $4/$1_$3_all > FASTA/$1.$3.fasta
}

function extract_bowtie_stats(){
	#FUNCTION extract_bowtie_stats
	#Arguments: 
	scode=$1
	pattern=$2
	patternname=$3
	tempfolder=$4
	# This function parses the output of bowtie alignments
	tr -d '[a-z,A-z,#:]' < $tempfolder/$scode.$patternname.stats.txt > $tempfolder/$scode.$patternname.stats2.txt;                        
	awk 'BEGIN{FS=" "}NR==1{print $1}' $tempfolder/$scode.$patternname.stats2.txt >> $tempfolder/$scode.processed.txt
	awk 'BEGIN{FS=" "}NR==2{print $1}' $tempfolder/$scode.$patternname.stats2.txt >> $tempfolder/$scode.aligned.txt
	awk 'BEGIN{FS=" "}NR==3{print $1}' $tempfolder/$scode.$patternname.stats2.txt >> $tempfolder/$scode.failed.txt
	awk 'BEGIN{FS=" "}NR==4{print $2}' $tempfolder/$scode.$patternname.stats2.txt >> $tempfolder/$scode.suppressed.txt
}

function extract_bowtie_stats_by_genome(){
	#FUNCTION extract_bowtie_stats_by_genome
	#Arguments: 
	scode=$1
	tempfolder=$2
	# This function concatenates the outputs of bowtie alignments for all enzymes per genome
	echo -e -n '\nProcessed\t' >> $scode.txt
	while read data; do echo -n -e "$data\t" >> $scode.txt; done < $tempfolder/$scode.processed.txt	
	echo -n -e "\nAligned\t" >> $scode.txt
	while read data; do echo -n -e "$data\t" >> $scode.txt; done < $tempfolder/$scode.aligned.txt
	echo -n -e "\nFailed\t" >> $scode.txt	
	while read data; do echo -n -e "$data\t" >> $scode.txt; done < $tempfolder/$scode.failed.txt
	echo -n -e "\nSuppressed\t"  >> $scode.txt
	while read data; do echo -n -e "$data\t" >> $scode.txt; done < $tempfolder/$scode.suppressed.txt
}		

function extract_bowtie_stats_for_all(){
	#FUNCTION extract_bowtie_stats_for_all
	#Arguments: 
	scode=$1
	tempfolder=$2
	bowtieflag=$3
	# This function concatenates the outputs of bowtie alignments for all enzymes for all genomes - and $scode.count.txt
	if [ "$bowtieflag" == "YES" ]
		then
		echo -e -n "\n$scode\t" >> ALL.processed.txt
		while read data; do echo -n -e "$data\t" >> ALL.processed.txt; done < $tempfolder/$scode.processed.txt
		echo -e -n "\n$scode\t" >> ALL.aligned.txt
		while read data; do echo -n -e "$data\t" >> ALL.aligned.txt; done < $tempfolder/$scode.aligned.txt
		echo -e -n "\n$scode\t" >> ALL.failed.txt
		while read data; do echo -n -e "$data\t" >> ALL.failed.txt; done < $tempfolder/$scode.failed.txt
		echo -e -n "\n$scode\t" >> ALL.suppressed.txt
		while read data; do echo -n -e "$data\t" >> ALL.suppressed.txt; done < $tempfolder/$scode.suppressed.txt
	fi;	
	echo -e -n "\n$scode\t" >> ALL.count.txt
	while read data; do echo -n -e "$data\t" >> ALL.count.txt; done < $tempfolder/$scode.count.txt	
}

#---------------------------------MAIN-------------------------------------------
writelog "--------------------"  restriction_site_search.log
# Create directories to store the bowtie databases and the aligned tags against the genome
# default values flags
localfile="NO"
bowtieflag="YES"
mkdir bowtie_db aligned FASTA
tempfolder="TEMP_$(date +"%m_%d_%y_%H_%M")"
mkdir $tempfolder
paramsfile="$1"
. ./$paramsfile
writelog "Input arguments: patternfile $patternfile genomefile $genomefile localfile $localfile bowtieflag $bowtieflag"  restriction_site_search.log
# check if files are in other directory and copy to temp  directory
# get filename (not the path)
patternfilename=$(basename $patternfile)
#cp $patternfile $tempfolder/$patternfilename
# Sort patterns file
sort -k 2 $patternfile > $tempfolder/$patternfilename.sort
cp $tempfolder/$patternfilename.sort $patternfile	
# Check format of input files
check_inputfile $patternfile $tempfolder
check_inputfile $genomefile  $tempfolder
# flags
localfile="$(echo -e "${localfile}" | tr -d '[[:space:]]')"
bowtieflag="$(echo -e "${bowtieflag}" | tr -d '[[:space:]]')"
localfile=$(echo $localfile | tr '[a-z]' '[A-Z]')
bowtieflag=$(echo $bowtieflag | tr '[a-z]' '[A-Z]')
# For each assembly file
while read scode url; 
	do
	writelog "$scode process starts"  restriction_site_search.log
	# Obtain and process each assembly file
	# check if there is a local fasta file or we need to download the file from ncbi
	if [ "$localfile" == "YES" ]
	then 
		tr '[:lower:]' '[:upper:]' < $url > $scode.UP
		# obtain .up file from local file
		writelog "local file $url converted to $scode.UP"  restriction_site_search.log
	else
		process_inputfile $url $scode;
	fi;
	#check scode.UP content
	size=$(stat -c %s $scode.UP)
	if [ "$?" != "0" ]; 
		then
		#echo "Problem with $scode.UP"
		writelog "Error with $scode.UP size"  restriction_site_search.log
		exit 1
	fi;	
	# Perform initial counts of patterns
	count_nt_sites $scode $patternfile $tempfolder ;
	# Build bowtie index database for the genome
	if [ "$bowtieflag" == "YES" ]
		then
		cd bowtie_db
		bowtie-build -f ../$scode.UP $scode.Genome
		if [ "$?" != "0" ]; 
			then
			echo "bowtie-build cannot find input file $url or $scode.UP"
			writelog "bowtie-build cannot find input file $url or $scode.UP"  restriction_site_search.log
			exit 1
		fi;
		cd ..
		writelog "$scode processed files and bowtie index done"  restriction_site_search.log
	fi;	
	# Perform in silico RAD sequencing and map reads back to the genome assembly
	while read pattern patternname; 
		do 
		pattern_upstream_downstream $scode $pattern $patternname $tempfolder; 	
		if [ "$bowtieflag" == "YES" ]
			then
			writelog "bowtie $scode $patternname"  restriction_site_search.log
			eval "bowtie -v 3 --best --strata -m 1 -p 8 -f ./bowtie_db/$scode.Genome ./FASTA/$scode.$patternname.fasta ./aligned/$scode.$patternname.bowtie >> $tempfolder/$scode.$patternname.stats.txt 2>&1";
# Parse alignment outputs
		extract_bowtie_stats $scode $pattern $patternname $tempfolder;
		fi;
	done < $patternfile;
	if [ "$bowtieflag" == "YES" ]
		then
			extract_bowtie_stats_by_genome $scode $tempfolder
			writelog "$scode process done"	 restriction_site_search.log
	fi;	
	rm -rf aligned/* bowtie_db/* FASTA/*
#	rm -rf aligned/ bowtie_db/ FASTA/
	rm $scode.UP
done < $genomefile
while read scode url;
	do 
	extract_bowtie_stats_for_all $scode $tempfolder $bowtieflag
	writelog "stats for all done"	 restriction_site_search.log
done < $genomefile
rm -rf aligned/ bowtie_db/ FASTA/
rm -R $tempfolder


