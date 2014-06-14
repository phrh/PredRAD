# Created by Santiago Herrera on 11 June 2014
# Copyright (c) 2014 Santiago Herrera. All rights reserved.
#
# This file contains the analyses performed for the manuscript Herrera S, Reyes-Herrera, PH, Shank TM (in prep) Genome-wide predictability of restriction sites across the eukaryotic tree of life.
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2.

#Change working directory
#setwd("/Users/tiagohe/Documents/RADtag/Genomes/paper/scripts/")
# Load library
library(gplots)
library(plyr)
library(RColorBrewer)
library(Hmisc)
#Import output files
pre_counts<-read.table("ALL.count.txt", header=FALSE, sep="\t", blank.lines.skip=TRUE)
pre_size<-read.table("ALL.size.txt", header=FALSE, sep="\t", blank.lines.skip=TRUE)
pre_suppressed<-read.table("ALL.suppressed.txt", header=FALSE, sep="\t", blank.lines.skip=TRUE)
pre_aligned<-read.table("ALL.aligned.txt", header=FALSE, sep="\t", blank.lines.skip=TRUE)
ids<-read.table("phylo_ordered_names.txt", header=FALSE, sep="\t", blank.lines.skip=TRUE)
pre_gcditrinuc<-read.table("DistributionFile.txt", header=TRUE, sep="\t", blank.lines.skip=TRUE)
pre_monoprobs<-read.table("DistributionFile.txt_nt", header=TRUE, sep="\t", blank.lines.skip=TRUE)
pre_diprobs<-read.table("DistributionFile.txt_dint", header=TRUE, sep="\t", blank.lines.skip=TRUE)
pre_triprobs<-read.table("DistributionFile.txt_trint", header=TRUE, sep="\t", blank.lines.skip=TRUE)


# Add column names
colnames(pre_size) <- c("id","gen_size_amb","gen_size_unamb","gc")
colnames(ids) <- c("id")
colnames(pre_counts) <- c("id","AgeI","ApoI","BsrFI","EcoRI","FatI","KpnI","MluCI","MseI","MspI","NcoI","NgoMIV","NotI","NsiI","NspI","PciI","PstI","SbfI","SgrAI")
colnames(pre_suppressed) <- c("id","AgeI","ApoI","BsrFI","EcoRI","FatI","KpnI","MluCI","MseI","MspI","NcoI","NgoMIV","NotI","NsiI","NspI","PciI","PstI","SbfI","SgrAI")
colnames(pre_aligned) <- c("id","AgeI","ApoI","BsrFI","EcoRI","FatI","KpnI","MluCI","MseI","MspI","NcoI","NgoMIV","NotI","NsiI","NspI","PciI","PstI","SbfI","SgrAI")
colnames(pre_monoprobs) <- c("id","AgeI","ApoI","BsrFI","EcoRI","FatI","KpnI","MluCI","MseI","MspI","NcoI","NgoMIV","NotI","NsiI","NspI","PciI","PstI","SbfI","SgrAI")
colnames(pre_diprobs) <- c("id","AgeI","ApoI","BsrFI","EcoRI","FatI","KpnI","MluCI","MseI","MspI","NcoI","NgoMIV","NotI","NsiI","NspI","PciI","PstI","SbfI","SgrAI")
colnames(pre_triprobs) <- c("id","AgeI","ApoI","BsrFI","EcoRI","FatI","KpnI","MluCI","MseI","MspI","NcoI","NgoMIV","NotI","NsiI","NspI","PciI","PstI","SbfI","SgrAI")
size<-join(ids, pre_size, by= "id")
counts<-join(ids, pre_counts, by= "id")
suppressed<-join(ids, pre_suppressed, by= "id")
aligned<-join(ids, pre_aligned, by= "id")
gcditrinuc<-join(ids, pre_gcditrinuc, by= "id")
monoprobs<-join(ids, pre_monoprobs[,1:19], by= "id")
diprobs<-join(ids, pre_diprobs[,1:19], by= "id")
triprobs<-join(ids, pre_triprobs[,1:19], by= "id")

#Extract the species codes for each genome
scode<-size[,1]
scode_list<-as.vector(size[,1])
#Extract the unambiguous genome size
genome_size<-gcditrinuc[,3] #from [AGCT] counts in biophyton
#Extract the number of bases the genome and calculate frequencies
A_freq<-(gcditrinuc[,4])/genome_size
G_freq<-(gcditrinuc[,5])/genome_size
C_freq<-(gcditrinuc[,6])/genome_size
T_freq<-(gcditrinuc[,7])/genome_size
N_freq<-(gcditrinuc[,8])/genome_size
#Calculate the GC content of the genome
gc_content<-(G_freq+C_freq)
at_content<-1-(G_freq+C_freq)
gc_matrix<-cbind(ids,gc_content)

#Calculate observed dinucleotide and trinucleotide frequencies
di_total=0
for (i in 9:24) {di_total=di_total+gcditrinuc[,i] }
AA_freq<-(gcditrinuc[,9])/di_total
AG_freq<-(gcditrinuc[,10])/di_total
AC_freq<-(gcditrinuc[,11])/di_total
AT_freq<-(gcditrinuc[,12])/di_total
GA_freq<-(gcditrinuc[,13])/di_total
GG_freq<-(gcditrinuc[,14])/di_total
GC_freq<-(gcditrinuc[,15])/di_total
GT_freq<-(gcditrinuc[,16])/di_total
CA_freq<-(gcditrinuc[,17])/di_total
CG_freq<-(gcditrinuc[,18])/di_total
CC_freq<-(gcditrinuc[,19])/di_total
CT_freq<-(gcditrinuc[,20])/di_total
TA_freq<-(gcditrinuc[,21])/di_total
TG_freq<-(gcditrinuc[,22])/di_total
TC_freq<-(gcditrinuc[,23])/di_total
TT_freq<-(gcditrinuc[,24])/di_total
tri_total=0
for (i in 25:88) {tri_total=tri_total+gcditrinuc[,i] }
AAA_freq<-(gcditrinuc[,25])/tri_total
AAG_freq<-(gcditrinuc[,26])/tri_total
AAC_freq<-(gcditrinuc[,27])/tri_total
AAT_freq<-(gcditrinuc[,28])/tri_total
AGA_freq<-(gcditrinuc[,29])/tri_total
AGG_freq<-(gcditrinuc[,30])/tri_total
AGC_freq<-(gcditrinuc[,31])/tri_total
AGT_freq<-(gcditrinuc[,32])/tri_total
ACA_freq<-(gcditrinuc[,33])/tri_total
ACG_freq<-(gcditrinuc[,34])/tri_total
ACC_freq<-(gcditrinuc[,35])/tri_total
ACT_freq<-(gcditrinuc[,36])/tri_total
ATA_freq<-(gcditrinuc[,37])/tri_total
ATG_freq<-(gcditrinuc[,38])/tri_total
ATC_freq<-(gcditrinuc[,39])/tri_total
ATT_freq<-(gcditrinuc[,40])/tri_total
GAA_freq<-(gcditrinuc[,41])/tri_total
GAG_freq<-(gcditrinuc[,42])/tri_total
GAC_freq<-(gcditrinuc[,43])/tri_total
GAT_freq<-(gcditrinuc[,44])/tri_total
GGA_freq<-(gcditrinuc[,45])/tri_total
GGG_freq<-(gcditrinuc[,46])/tri_total
GGC_freq<-(gcditrinuc[,47])/tri_total
GGT_freq<-(gcditrinuc[,48])/tri_total
GCA_freq<-(gcditrinuc[,49])/tri_total
GCG_freq<-(gcditrinuc[,50])/tri_total
GCC_freq<-(gcditrinuc[,51])/tri_total
GCT_freq<-(gcditrinuc[,52])/tri_total
GTA_freq<-(gcditrinuc[,53])/tri_total
GTG_freq<-(gcditrinuc[,54])/tri_total
GTC_freq<-(gcditrinuc[,55])/tri_total
GTT_freq<-(gcditrinuc[,56])/tri_total
CAA_freq<-(gcditrinuc[,57])/tri_total
CAG_freq<-(gcditrinuc[,58])/tri_total
CAC_freq<-(gcditrinuc[,59])/tri_total
CAT_freq<-(gcditrinuc[,60])/tri_total
CGA_freq<-(gcditrinuc[,61])/tri_total
CGG_freq<-(gcditrinuc[,62])/tri_total
CGC_freq<-(gcditrinuc[,63])/tri_total
CGT_freq<-(gcditrinuc[,64])/tri_total
CCA_freq<-(gcditrinuc[,65])/tri_total
CCG_freq<-(gcditrinuc[,66])/tri_total
CCC_freq<-(gcditrinuc[,67])/tri_total
CCT_freq<-(gcditrinuc[,68])/tri_total
CTA_freq<-(gcditrinuc[,69])/tri_total
CTG_freq<-(gcditrinuc[,70])/tri_total
CTC_freq<-(gcditrinuc[,71])/tri_total
CTT_freq<-(gcditrinuc[,72])/tri_total
TAA_freq<-(gcditrinuc[,73])/tri_total
TAG_freq<-(gcditrinuc[,74])/tri_total
TAC_freq<-(gcditrinuc[,75])/tri_total
TAT_freq<-(gcditrinuc[,76])/tri_total
TGA_freq<-(gcditrinuc[,77])/tri_total
TGG_freq<-(gcditrinuc[,78])/tri_total
TGC_freq<-(gcditrinuc[,79])/tri_total
TGT_freq<-(gcditrinuc[,80])/tri_total
TCA_freq<-(gcditrinuc[,81])/tri_total
TCG_freq<-(gcditrinuc[,82])/tri_total
TCC_freq<-(gcditrinuc[,83])/tri_total
TCT_freq<-(gcditrinuc[,84])/tri_total
TTA_freq<-(gcditrinuc[,85])/tri_total
TTG_freq<-(gcditrinuc[,86])/tri_total
TTC_freq<-(gcditrinuc[,87])/tri_total
TTT_freq<-(gcditrinuc[,88])/tri_total


#Calculate frequency of each mono, di and trinucleotide taking int account the antiparallel
#structure of double stranded DNA

A_freqstar<-(A_freq+T_freq)/2
T_freqstar<-(T_freq+A_freq)/2
G_freqstar<-(G_freq+C_freq)/2
C_freqstar<-(C_freq+G_freq)/2
AA_freqstar<-(AA_freq+TT_freq)/2
AG_freqstar<-(AG_freq+CT_freq)/2
AC_freqstar<-(AC_freq+GT_freq)/2
AT_freqstar<-(AT_freq+AT_freq)/2
GA_freqstar<-(GA_freq+TC_freq)/2
GG_freqstar<-(GG_freq+CC_freq)/2
GC_freqstar<-(GC_freq+GC_freq)/2
GT_freqstar<-(GT_freq+AC_freq)/2
CA_freqstar<-(CA_freq+TG_freq)/2
CG_freqstar<-(CG_freq+CG_freq)/2
CC_freqstar<-(CC_freq+GG_freq)/2
CT_freqstar<-(CT_freq+AG_freq)/2
TA_freqstar<-(TA_freq+TA_freq)/2
TG_freqstar<-(TG_freq+CA_freq)/2
TC_freqstar<-(TC_freq+GA_freq)/2
TT_freqstar<-(TT_freq+AA_freq)/2
AAA_freqstar<-(AAA_freq+TTT_freq)/2
AAG_freqstar<-(AAG_freq+CTT_freq)/2
AAC_freqstar<-(AAC_freq+GTT_freq)/2
AAT_freqstar<-(AAT_freq+ATT_freq)/2
AGA_freqstar<-(AGA_freq+TCT_freq)/2
AGG_freqstar<-(AGG_freq+CCT_freq)/2
AGC_freqstar<-(AGC_freq+GCT_freq)/2
AGT_freqstar<-(AGT_freq+ACT_freq)/2
ACA_freqstar<-(ACA_freq+TGT_freq)/2
ACG_freqstar<-(ACG_freq+CGT_freq)/2
ACC_freqstar<-(ACC_freq+GGT_freq)/2
ACT_freqstar<-(ACT_freq+AGT_freq)/2
ATA_freqstar<-(ATA_freq+TAT_freq)/2
ATG_freqstar<-(ATG_freq+CAT_freq)/2
ATC_freqstar<-(ATC_freq+GAT_freq)/2
ATT_freqstar<-(ATT_freq+AAT_freq)/2
GAA_freqstar<-(GAA_freq+TTC_freq)/2
GAG_freqstar<-(GAG_freq+CTC_freq)/2
GAC_freqstar<-(GAC_freq+GTC_freq)/2
GAT_freqstar<-(GAT_freq+ATC_freq)/2
GGA_freqstar<-(GGA_freq+TCC_freq)/2
GGG_freqstar<-(GGG_freq+CCC_freq)/2
GGC_freqstar<-(GGC_freq+GCC_freq)/2
GGT_freqstar<-(GGT_freq+ACC_freq)/2
GCA_freqstar<-(GCA_freq+TGC_freq)/2
GCG_freqstar<-(GCG_freq+CGC_freq)/2
GCC_freqstar<-(GCC_freq+GGC_freq)/2
GCT_freqstar<-(GCT_freq+AGC_freq)/2
GTA_freqstar<-(GTA_freq+TAC_freq)/2
GTG_freqstar<-(GTG_freq+CAC_freq)/2
GTC_freqstar<-(GTC_freq+GAC_freq)/2
GTT_freqstar<-(GTT_freq+AAC_freq)/2
CAA_freqstar<-(CAA_freq+TTG_freq)/2
CAG_freqstar<-(CAG_freq+CTG_freq)/2
CAC_freqstar<-(CAC_freq+GTG_freq)/2
CAT_freqstar<-(CAT_freq+ATG_freq)/2
CGA_freqstar<-(CGA_freq+TCG_freq)/2
CGG_freqstar<-(CGG_freq+CCG_freq)/2
CGC_freqstar<-(CGC_freq+GCG_freq)/2
CGT_freqstar<-(CGT_freq+ACG_freq)/2
CCA_freqstar<-(CCA_freq+TGG_freq)/2
CCG_freqstar<-(CCG_freq+CGG_freq)/2
CCC_freqstar<-(CCC_freq+GGG_freq)/2
CCT_freqstar<-(CCT_freq+AGG_freq)/2
CTA_freqstar<-(CTA_freq+TAG_freq)/2
CTG_freqstar<-(CTG_freq+CAG_freq)/2
CTC_freqstar<-(CTC_freq+GAG_freq)/2
CTT_freqstar<-(CTT_freq+AAG_freq)/2
TAA_freqstar<-(TAA_freq+TTA_freq)/2
TAG_freqstar<-(TAG_freq+CTA_freq)/2
TAC_freqstar<-(TAC_freq+GTA_freq)/2
TAT_freqstar<-(TAT_freq+ATA_freq)/2
TGA_freqstar<-(TGA_freq+TCA_freq)/2
TGG_freqstar<-(TGG_freq+CCA_freq)/2
TGC_freqstar<-(TGC_freq+GCA_freq)/2
TGT_freqstar<-(TGT_freq+ACA_freq)/2
TCA_freqstar<-(TCA_freq+TGA_freq)/2
TCG_freqstar<-(TCG_freq+CGA_freq)/2
TCC_freqstar<-(TCC_freq+GGA_freq)/2
TCT_freqstar<-(TCT_freq+AGA_freq)/2
TTA_freqstar<-(TTA_freq+TAA_freq)/2
TTG_freqstar<-(TTG_freq+CAA_freq)/2
TTC_freqstar<-(TTC_freq+GAA_freq)/2
TTT_freqstar<-(TTT_freq+AAA_freq)/2

ANA_freqstar<-(AAA_freqstar+AGA_freqstar+ACA_freqstar+ATA_freqstar)
ANG_freqstar<-(AGG_freqstar+ACG_freqstar+ATG_freqstar+AAG_freqstar)
ANC_freqstar<-(ACC_freqstar+ATC_freqstar+AAC_freqstar+AGC_freqstar)
ANT_freqstar<-(ATT_freqstar+AAT_freqstar+AGT_freqstar+ACT_freqstar)
GNA_freqstar<-(GAA_freqstar+GGA_freqstar+GCA_freqstar+GTA_freqstar)
GNG_freqstar<-(GGG_freqstar+GCG_freqstar+GTG_freqstar+GAG_freqstar)
GNC_freqstar<-(GCC_freqstar+GTC_freqstar+GAC_freqstar+GGC_freqstar)
GNT_freqstar<-(GTT_freqstar+GAT_freqstar+GGT_freqstar+GCT_freqstar)
CNA_freqstar<-(CAA_freqstar+CGA_freqstar+CCA_freqstar+CTA_freqstar)
CNG_freqstar<-(CGG_freqstar+CCG_freqstar+CTG_freqstar+CAG_freqstar)
CNC_freqstar<-(CCC_freqstar+CTC_freqstar+CAC_freqstar+CGC_freqstar)
CNT_freqstar<-(CTT_freqstar+CAT_freqstar+CGT_freqstar+CCT_freqstar)
TNA_freqstar<-(TAA_freqstar+TGA_freqstar+TCA_freqstar+TTA_freqstar)
TNG_freqstar<-(TGG_freqstar+TCG_freqstar+TTG_freqstar+TAG_freqstar)
TNC_freqstar<-(TCC_freqstar+TTC_freqstar+TAC_freqstar+TGC_freqstar)
TNT_freqstar<-(TTT_freqstar+TAT_freqstar+TGT_freqstar+TCT_freqstar)

#Calculate dinucleotide odds ratio (rho star) that accounts for the complementary antiparallel
#structure of double stranded DNA
CG_rhostar<-(CG_freqstar/(C_freqstar*G_freqstar))
GC_rhostar<-(GC_freqstar/(G_freqstar*C_freqstar))
AT_rhostar<-(AT_freqstar/(A_freqstar*T_freqstar))
TA_rhostar<-(TA_freqstar/(T_freqstar*A_freqstar))
GG_rhostar<-(GG_freqstar/(G_freqstar*G_freqstar))
CC_rhostar<-(CC_freqstar/(C_freqstar*C_freqstar))
AA_rhostar<-(AA_freqstar/(A_freqstar*A_freqstar))
TT_rhostar<-(TT_freqstar/(T_freqstar*T_freqstar))
AC_rhostar<-(AC_freqstar/(A_freqstar*C_freqstar))
GT_rhostar<-(GT_freqstar/(G_freqstar*T_freqstar))
CA_rhostar<-(CA_freqstar/(C_freqstar*A_freqstar))
TG_rhostar<-(TG_freqstar/(T_freqstar*G_freqstar))
AG_rhostar<-(AG_freqstar/(A_freqstar*G_freqstar))
CT_rhostar<-(CT_freqstar/(C_freqstar*T_freqstar))
GA_rhostar<-(GA_freqstar/(G_freqstar*A_freqstar))
TC_rhostar<-(TC_freqstar/(T_freqstar*C_freqstar))
dinuc_rhostar<-cbind(ids,CG_rhostar,GC_rhostar,AT_rhostar,TA_rhostar,GG_rhostar,CC_rhostar,AA_rhostar,TT_rhostar,AC_rhostar,GT_rhostar,CA_rhostar,TG_rhostar,AG_rhostar,CT_rhostar,GA_rhostar,TC_rhostar)
colnames(dinuc_rhostar) <- c("id","CG","GC","AT","TA","GG","CC","AA","TT","AC","GT","CA","TG","AG","CT","GA","TC")
row.names(dinuc_rhostar) <- dinuc_rhostar$id
dinuc_rhostar_frame<- dinuc_rhostar[,2:17]
dinuc_rhostar_matrix <- data.matrix(dinuc_rhostar_frame)

#Calculate trinucleotide odds ratio (gamma star) that accounts for the complementary antiparallel
#structure of double stranded DNA
AAA_gammastar<-((AAA_freqstar*A_freqstar*A_freqstar*A_freqstar)/(AA_freqstar*AA_freqstar*ANA_freqstar))
AAG_gammastar<-((AAG_freqstar*A_freqstar*A_freqstar*G_freqstar)/(AA_freqstar*AG_freqstar*ANG_freqstar))
AAC_gammastar<-((AAC_freqstar*A_freqstar*A_freqstar*C_freqstar)/(AA_freqstar*AC_freqstar*ANC_freqstar))
AAT_gammastar<-((AAT_freqstar*A_freqstar*A_freqstar*T_freqstar)/(AA_freqstar*AT_freqstar*ANT_freqstar))
AGA_gammastar<-((AGA_freqstar*A_freqstar*G_freqstar*A_freqstar)/(AG_freqstar*GA_freqstar*ANA_freqstar))
AGG_gammastar<-((AGG_freqstar*A_freqstar*G_freqstar*G_freqstar)/(AG_freqstar*GG_freqstar*ANG_freqstar))
AGC_gammastar<-((AGC_freqstar*A_freqstar*G_freqstar*C_freqstar)/(AG_freqstar*GC_freqstar*ANC_freqstar))
AGT_gammastar<-((AGT_freqstar*A_freqstar*G_freqstar*T_freqstar)/(AG_freqstar*GT_freqstar*ANT_freqstar))
ACA_gammastar<-((ACA_freqstar*A_freqstar*C_freqstar*A_freqstar)/(AC_freqstar*CA_freqstar*ANA_freqstar))
ACG_gammastar<-((ACG_freqstar*A_freqstar*C_freqstar*G_freqstar)/(AC_freqstar*CG_freqstar*ANG_freqstar))
ACC_gammastar<-((ACC_freqstar*A_freqstar*C_freqstar*C_freqstar)/(AC_freqstar*CC_freqstar*ANC_freqstar))
ACT_gammastar<-((ACT_freqstar*A_freqstar*C_freqstar*T_freqstar)/(AC_freqstar*CT_freqstar*ANT_freqstar))
ATA_gammastar<-((ATA_freqstar*A_freqstar*T_freqstar*A_freqstar)/(AT_freqstar*TA_freqstar*ANA_freqstar))
ATG_gammastar<-((ATG_freqstar*A_freqstar*T_freqstar*G_freqstar)/(AT_freqstar*TG_freqstar*ANG_freqstar))
ATC_gammastar<-((ATC_freqstar*A_freqstar*T_freqstar*C_freqstar)/(AT_freqstar*TC_freqstar*ANC_freqstar))
ATT_gammastar<-((ATT_freqstar*A_freqstar*T_freqstar*T_freqstar)/(AT_freqstar*TT_freqstar*ANT_freqstar))
GAA_gammastar<-((GAA_freqstar*G_freqstar*A_freqstar*A_freqstar)/(GA_freqstar*AA_freqstar*GNA_freqstar))
GAG_gammastar<-((GAG_freqstar*G_freqstar*A_freqstar*G_freqstar)/(GA_freqstar*AG_freqstar*GNG_freqstar))
GAC_gammastar<-((GAC_freqstar*G_freqstar*A_freqstar*C_freqstar)/(GA_freqstar*AC_freqstar*GNC_freqstar))
GAT_gammastar<-((GAT_freqstar*G_freqstar*A_freqstar*T_freqstar)/(GA_freqstar*AT_freqstar*GNT_freqstar))
GGA_gammastar<-((GGA_freqstar*G_freqstar*G_freqstar*A_freqstar)/(GG_freqstar*GA_freqstar*GNA_freqstar))
GGG_gammastar<-((GGG_freqstar*G_freqstar*G_freqstar*G_freqstar)/(GG_freqstar*GG_freqstar*GNG_freqstar))
GGC_gammastar<-((GGC_freqstar*G_freqstar*G_freqstar*C_freqstar)/(GG_freqstar*GC_freqstar*GNC_freqstar))
GGT_gammastar<-((GGT_freqstar*G_freqstar*G_freqstar*T_freqstar)/(GG_freqstar*GT_freqstar*GNT_freqstar))
GCA_gammastar<-((GCA_freqstar*G_freqstar*C_freqstar*A_freqstar)/(GC_freqstar*CA_freqstar*GNA_freqstar))
GCG_gammastar<-((GCG_freqstar*G_freqstar*C_freqstar*G_freqstar)/(GC_freqstar*CG_freqstar*GNG_freqstar))
GCC_gammastar<-((GCC_freqstar*G_freqstar*C_freqstar*C_freqstar)/(GC_freqstar*CC_freqstar*GNC_freqstar))
GCT_gammastar<-((GCT_freqstar*G_freqstar*C_freqstar*T_freqstar)/(GC_freqstar*CT_freqstar*GNT_freqstar))
GTA_gammastar<-((GTA_freqstar*G_freqstar*T_freqstar*A_freqstar)/(GT_freqstar*TA_freqstar*GNA_freqstar))
GTG_gammastar<-((GTG_freqstar*G_freqstar*T_freqstar*G_freqstar)/(GT_freqstar*TG_freqstar*GNG_freqstar))
GTC_gammastar<-((GTC_freqstar*G_freqstar*T_freqstar*C_freqstar)/(GT_freqstar*TC_freqstar*GNC_freqstar))
GTT_gammastar<-((GTT_freqstar*G_freqstar*T_freqstar*T_freqstar)/(GT_freqstar*TT_freqstar*GNT_freqstar))
CAA_gammastar<-((CAA_freqstar*C_freqstar*A_freqstar*A_freqstar)/(CA_freqstar*AA_freqstar*CNA_freqstar))
CAG_gammastar<-((CAG_freqstar*C_freqstar*A_freqstar*G_freqstar)/(CA_freqstar*AG_freqstar*CNG_freqstar))
CAC_gammastar<-((CAC_freqstar*C_freqstar*A_freqstar*C_freqstar)/(CA_freqstar*AC_freqstar*CNC_freqstar))
CAT_gammastar<-((CAT_freqstar*C_freqstar*A_freqstar*T_freqstar)/(CA_freqstar*AT_freqstar*CNT_freqstar))
CGA_gammastar<-((CGA_freqstar*C_freqstar*G_freqstar*A_freqstar)/(CG_freqstar*GA_freqstar*CNA_freqstar))
CGG_gammastar<-((CGG_freqstar*C_freqstar*G_freqstar*G_freqstar)/(CG_freqstar*GG_freqstar*CNG_freqstar))
CGC_gammastar<-((CGC_freqstar*C_freqstar*G_freqstar*C_freqstar)/(CG_freqstar*GC_freqstar*CNC_freqstar))
CGT_gammastar<-((CGT_freqstar*C_freqstar*G_freqstar*T_freqstar)/(CG_freqstar*GT_freqstar*CNT_freqstar))
CCA_gammastar<-((CCA_freqstar*C_freqstar*C_freqstar*A_freqstar)/(CC_freqstar*CA_freqstar*CNA_freqstar))
CCG_gammastar<-((CCG_freqstar*C_freqstar*C_freqstar*G_freqstar)/(CC_freqstar*CG_freqstar*CNG_freqstar))
CCC_gammastar<-((CCC_freqstar*C_freqstar*C_freqstar*C_freqstar)/(CC_freqstar*CC_freqstar*CNC_freqstar))
CCT_gammastar<-((CCT_freqstar*C_freqstar*C_freqstar*T_freqstar)/(CC_freqstar*CT_freqstar*CNT_freqstar))
CTA_gammastar<-((CTA_freqstar*C_freqstar*T_freqstar*A_freqstar)/(CT_freqstar*TA_freqstar*CNA_freqstar))
CTG_gammastar<-((CTG_freqstar*C_freqstar*T_freqstar*G_freqstar)/(CT_freqstar*TG_freqstar*CNG_freqstar))
CTC_gammastar<-((CTC_freqstar*C_freqstar*T_freqstar*C_freqstar)/(CT_freqstar*TC_freqstar*CNC_freqstar))
CTT_gammastar<-((CTT_freqstar*C_freqstar*T_freqstar*T_freqstar)/(CT_freqstar*TT_freqstar*CNT_freqstar))
TAA_gammastar<-((TAA_freqstar*T_freqstar*A_freqstar*A_freqstar)/(TA_freqstar*AA_freqstar*TNA_freqstar))
TAG_gammastar<-((TAG_freqstar*T_freqstar*A_freqstar*G_freqstar)/(TA_freqstar*AG_freqstar*TNG_freqstar))
TAC_gammastar<-((TAC_freqstar*T_freqstar*A_freqstar*C_freqstar)/(TA_freqstar*AC_freqstar*TNC_freqstar))
TAT_gammastar<-((TAT_freqstar*T_freqstar*A_freqstar*T_freqstar)/(TA_freqstar*AT_freqstar*TNT_freqstar))
TGA_gammastar<-((TGA_freqstar*T_freqstar*G_freqstar*A_freqstar)/(TG_freqstar*GA_freqstar*TNA_freqstar))
TGG_gammastar<-((TGG_freqstar*T_freqstar*G_freqstar*G_freqstar)/(TG_freqstar*GG_freqstar*TNG_freqstar))
TGC_gammastar<-((TGC_freqstar*T_freqstar*G_freqstar*C_freqstar)/(TG_freqstar*GC_freqstar*TNC_freqstar))
TGT_gammastar<-((TGT_freqstar*T_freqstar*G_freqstar*T_freqstar)/(TG_freqstar*GT_freqstar*TNT_freqstar))
TCA_gammastar<-((TCA_freqstar*T_freqstar*C_freqstar*A_freqstar)/(TC_freqstar*CA_freqstar*TNA_freqstar))
TCG_gammastar<-((TCG_freqstar*T_freqstar*C_freqstar*G_freqstar)/(TC_freqstar*CG_freqstar*TNG_freqstar))
TCC_gammastar<-((TCC_freqstar*T_freqstar*C_freqstar*C_freqstar)/(TC_freqstar*CC_freqstar*TNC_freqstar))
TCT_gammastar<-((TCT_freqstar*T_freqstar*C_freqstar*T_freqstar)/(TC_freqstar*CT_freqstar*TNT_freqstar))
TTA_gammastar<-((TTA_freqstar*T_freqstar*T_freqstar*A_freqstar)/(TT_freqstar*TA_freqstar*TNA_freqstar))
TTG_gammastar<-((TTG_freqstar*T_freqstar*T_freqstar*G_freqstar)/(TT_freqstar*TG_freqstar*TNG_freqstar))
TTC_gammastar<-((TTC_freqstar*T_freqstar*T_freqstar*C_freqstar)/(TT_freqstar*TC_freqstar*TNC_freqstar))
TTT_gammastar<-((TTT_freqstar*T_freqstar*T_freqstar*T_freqstar)/(TT_freqstar*TT_freqstar*TNT_freqstar))
trinuc_gammastar<-cbind(ids,AAA_gammastar,AAG_gammastar,AAC_gammastar,AAT_gammastar,AGA_gammastar,AGG_gammastar,AGC_gammastar,AGT_gammastar,ACA_gammastar,ACG_gammastar,ACC_gammastar,ACT_gammastar,ATA_gammastar,ATG_gammastar,ATC_gammastar,ATT_gammastar,GAA_gammastar,GAG_gammastar,GAC_gammastar,GAT_gammastar,GGA_gammastar,GGG_gammastar,GGC_gammastar,GGT_gammastar,GCA_gammastar,GCG_gammastar,GCC_gammastar,GCT_gammastar,GTA_gammastar,GTG_gammastar,GTC_gammastar,GTT_gammastar,CAA_gammastar,CAG_gammastar,CAC_gammastar,CAT_gammastar,CGA_gammastar,CGG_gammastar,CGC_gammastar,CGT_gammastar,CCA_gammastar,CCG_gammastar,CCC_gammastar,CCT_gammastar,CTA_gammastar,CTG_gammastar,CTC_gammastar,CTT_gammastar,TAA_gammastar,TAG_gammastar,TAC_gammastar,TAT_gammastar,TGA_gammastar,TGG_gammastar,TGC_gammastar,TGT_gammastar,TCA_gammastar,TCG_gammastar,TCC_gammastar,TCT_gammastar,TTA_gammastar,TTG_gammastar,TTC_gammastar,TTT_gammastar)
colnames(trinuc_gammastar) <- c("id","AAA","AAG","AAC","AAT","AGA","AGG","AGC","AGT","ACA","ACG","ACC","ACT","ATA","ATG","ATC","ATT","GAA","GAG","GAC","GAT","GGA","GGG","GGC","GGT","GCA","GCG","GCC","GCT","GTA","GTG","GTC","GTT","CAA","CAG","CAC","CAT","CGA","CGG","CGC","CGT","CCA","CCG","CCC","CCT","CTA","CTG","CTC","CTT","TAA","TAG","TAC","TAT","TGA","TGG","TGC","TGT","TCA","TCG","TCC","TCT","TTA","TTG","TTC","TTT")
row.names(trinuc_gammastar) <- trinuc_gammastar$id
trinuc_gammastar_frame<- trinuc_gammastar[,2:65]
trinuc_gammastar_matrix <- data.matrix(trinuc_gammastar_frame)

#Calculate expected frequency of cut sites with each enzyme for gc
AgeI_expgcfreq<-(at_content/2)*(gc_content/2)*(gc_content/2)*(gc_content/2)*(gc_content/2)*(at_content/2)
ApoI_expgcfreq<-((at_content/2)+(gc_content/2))*(at_content/2)*(at_content/2)*(at_content/2)*(at_content/2)*((at_content/2)+(gc_content/2))
BsrFI_expgcfreq<-((at_content/2)+(gc_content/2))*(gc_content/2)*(gc_content/2)*(gc_content/2)*(gc_content/2)*((at_content/2)+(gc_content/2))
EcoRI_expgcfreq<-(gc_content/2)*(at_content/2)*(at_content/2)*(at_content/2)*(at_content/2)*(gc_content/2)
FatI_expgcfreq<-(gc_content/2)*(at_content/2)*(at_content/2)*(gc_content/2)
KpnI_expgcfreq<-(gc_content/2)*(gc_content/2)*(at_content/2)*(at_content/2)*(gc_content/2)*(gc_content/2)
MluCI_expgcfreq<-(at_content/2)*(at_content/2)*(at_content/2)*(at_content/2)
MseI_expgcfreq<-(at_content/2)*(at_content/2)*(at_content/2)*(at_content/2)
MspI_expgcfreq<-(gc_content/2)*(gc_content/2)*(gc_content/2)*(gc_content/2)
NcoI_expgcfreq<-(gc_content/2)*(gc_content/2)*(at_content/2)*(at_content/2)*(gc_content/2)*(gc_content/2)
NgoMIV_expgcfreq<-(gc_content/2)*(gc_content/2)*(gc_content/2)*(gc_content/2)*(gc_content/2)*(gc_content/2)
NotI_expgcfreq<-(gc_content/2)*(gc_content/2)*(gc_content/2)*(gc_content/2)*(gc_content/2)*(gc_content/2)*(gc_content/2)*(gc_content/2)
NsiI_expgcfreq<-(at_content/2)*(at_content/2)*(gc_content/2)*(gc_content/2)*(at_content/2)*(at_content/2)
NspI_expgcfreq<-((at_content/2)+(gc_content/2))*(gc_content/2)*(at_content/2)*(at_content/2)*(gc_content/2)*((at_content/2)+(gc_content/2))
PciI_expgcfreq<-(at_content/2)*(gc_content/2)*(at_content/2)*(at_content/2)*(gc_content/2)*(at_content/2)
PstI_expgcfreq<-(gc_content/2)*(at_content/2)*(gc_content/2)*(gc_content/2)*(at_content/2)*(gc_content/2)
SbfI_expgcfreq<-(gc_content/2)*(gc_content/2)*(at_content/2)*(gc_content/2)*(gc_content/2)*(at_content/2)*(gc_content/2)*(gc_content/2)
SgrAI_expgcfreq<-(gc_content/2)*((at_content/2)+(gc_content/2))*(gc_content/2)*(gc_content/2)*(gc_content/2)*(gc_content/2)*((at_content/2)+(gc_content/2))*(gc_content/2)

expgcfreq<-cbind(ids,AgeI_expgcfreq,ApoI_expgcfreq,BsrFI_expgcfreq,EcoRI_expgcfreq,FatI_expgcfreq,KpnI_expgcfreq,MluCI_expgcfreq,MseI_expgcfreq,MspI_expgcfreq,NcoI_expgcfreq,NgoMIV_expgcfreq,NotI_expgcfreq,NsiI_expgcfreq,NspI_expgcfreq,PciI_expgcfreq,PstI_expgcfreq,SbfI_expgcfreq,SgrAI_expgcfreq)
colnames(expgcfreq) <- c("id","AgeI","ApoI","BsrFI","EcoRI","FatI","KpnI","MluCI","MseI","MspI","NcoI","NgoMIV","NotI","NsiI","NspI","PciI","PstI","SbfI","SgrAI")
row.names(expgcfreq) <- expgcfreq$id
expgcfreq_frame<- expgcfreq[,2:19]
expgcfreq_matrix <- data.matrix(expgcfreq_frame)

#Calculate expected number of cut sites with each enzyme for gc
AgeI_expgc<-AgeI_expgcfreq*genome_size
ApoI_expgc<-ApoI_expgcfreq*genome_size
BsrFI_expgc<-BsrFI_expgcfreq*genome_size
EcoRI_expgc<-EcoRI_expgcfreq*genome_size
FatI_expgc<-FatI_expgcfreq*genome_size
KpnI_expgc<-KpnI_expgcfreq*genome_size
MluCI_expgc<-MluCI_expgcfreq*genome_size
MseI_expgc<-MseI_expgcfreq*genome_size
MspI_expgc<-MspI_expgcfreq*genome_size
NcoI_expgc<-NcoI_expgcfreq*genome_size
NgoMIV_expgc<-NgoMIV_expgcfreq*genome_size
NotI_expgc<-NotI_expgcfreq*genome_size
NsiI_expgc<-NsiI_expgcfreq*genome_size
NspI_expgc<-NspI_expgcfreq*genome_size
PciI_expgc<-PciI_expgcfreq*genome_size
PstI_expgc<-PstI_expgcfreq*genome_size
SbfI_expgc<-SbfI_expgcfreq*genome_size
SgrAI_expgc<-SgrAI_expgcfreq*genome_size
expgcected<-cbind(ids,AgeI_expgc,ApoI_expgc,BsrFI_expgc,EcoRI_expgc,FatI_expgc,KpnI_expgc,MluCI_expgc,MseI_expgc,MspI_expgc,NcoI_expgc,NgoMIV_expgc,NotI_expgc,NsiI_expgc,NspI_expgc,PciI_expgc,PstI_expgc,SbfI_expgc,SgrAI_expgc)
colnames(expgcected) <- c("id","AgeI","ApoI","BsrFI","EcoRI","FatI","KpnI","MluCI","MseI","MspI","NcoI","NgoMIV","NotI","NsiI","NspI","PciI","PstI","SbfI","SgrAI")
row.names(expgcected) <- expgcected$id
expgc_frame<- expgcected[,2:19]
expgc_matrix <- data.matrix(expgc_frame)

#Expected frequency of cut sites with each enzyme for mono
AgeI_expmonofreq<-monoprobs[,2]
ApoI_expmonofreq<-monoprobs[,3]
BsrFI_expmonofreq<-monoprobs[,4]
EcoRI_expmonofreq<-monoprobs[,5]
FatI_expmonofreq<-monoprobs[,6]
KpnI_expmonofreq<-monoprobs[,7]
MluCI_expmonofreq<-monoprobs[,8]
MseI_expmonofreq<-monoprobs[,9]
MspI_expmonofreq<-monoprobs[,10]
NcoI_expmonofreq<-monoprobs[,11]
NgoMIV_expmonofreq<-monoprobs[,12]
NotI_expmonofreq<-monoprobs[,13]
NsiI_expmonofreq<-monoprobs[,14]
NspI_expmonofreq<-monoprobs[,15]
PciI_expmonofreq<-monoprobs[,16]
PstI_expmonofreq<-monoprobs[,17]
SbfI_expmonofreq<-monoprobs[,18]
SgrAI_expmonofreq<-monoprobs[,19]
expmonofreq<-cbind(ids,AgeI_expmonofreq,ApoI_expmonofreq,BsrFI_expmonofreq,EcoRI_expmonofreq,FatI_expmonofreq,KpnI_expmonofreq,MluCI_expmonofreq,MseI_expmonofreq,MspI_expmonofreq,NcoI_expmonofreq,NgoMIV_expmonofreq,NotI_expmonofreq,NsiI_expmonofreq,NspI_expmonofreq,PciI_expmonofreq,PstI_expmonofreq,SbfI_expmonofreq,SgrAI_expmonofreq)
colnames(expmonofreq) <- c("id","AgeI","ApoI","BsrFI","EcoRI","FatI","KpnI","MluCI","MseI","MspI","NcoI","NgoMIV","NotI","NsiI","NspI","PciI","PstI","SbfI","SgrAI")
row.names(expmonofreq) <- expmonofreq$id
expmonofreq_frame<- expmonofreq[,2:19]
expmonofreq_matrix <- data.matrix(expmonofreq_frame)


#Calculate expected number of cut sites with each enzyme for mono
AgeI_expmono<-AgeI_expmonofreq*genome_size
ApoI_expmono<-ApoI_expmonofreq*genome_size
BsrFI_expmono<-BsrFI_expmonofreq*genome_size
EcoRI_expmono<-EcoRI_expmonofreq*genome_size
FatI_expmono<-FatI_expmonofreq*genome_size
KpnI_expmono<-KpnI_expmonofreq*genome_size
MluCI_expmono<-MluCI_expmonofreq*genome_size
MseI_expmono<-MseI_expmonofreq*genome_size
MspI_expmono<-MspI_expmonofreq*genome_size
NcoI_expmono<-NcoI_expmonofreq*genome_size
NgoMIV_expmono<-NgoMIV_expmonofreq*genome_size
NotI_expmono<-NotI_expmonofreq*genome_size
NsiI_expmono<-NsiI_expmonofreq*genome_size
NspI_expmono<-NspI_expmonofreq*genome_size
PciI_expmono<-PciI_expmonofreq*genome_size
PstI_expmono<-PstI_expmonofreq*genome_size
SbfI_expmono<-SbfI_expmonofreq*genome_size
SgrAI_expmono<-SgrAI_expmonofreq*genome_size
expected<-cbind(ids,AgeI_expmono,ApoI_expmono,BsrFI_expmono,EcoRI_expmono,FatI_expmono,KpnI_expmono,MluCI_expmono,MseI_expmono,MspI_expmono,NcoI_expmono,NgoMIV_expmono,NotI_expmono,NsiI_expmono,NspI_expmono,PciI_expmono,PstI_expmono,SbfI_expmono,SgrAI_expmono)
colnames(expected) <- c("id","AgeI","ApoI","BsrFI","EcoRI","FatI","KpnI","MluCI","MseI","MspI","NcoI","NgoMIV","NotI","NsiI","NspI","PciI","PstI","SbfI","SgrAI")
row.names(expected) <- expected$id
expmono_frame<- expected[,2:19]
expmono_matrix <- data.matrix(expmono_frame)

#Expected frequency of cut sites with each enzyme for di
AgeI_expdifreq<-diprobs[,2]
ApoI_expdifreq<-diprobs[,3]
BsrFI_expdifreq<-diprobs[,4]
EcoRI_expdifreq<-diprobs[,5]
FatI_expdifreq<-diprobs[,6]
KpnI_expdifreq<-diprobs[,7]
MluCI_expdifreq<-diprobs[,8]
MseI_expdifreq<-diprobs[,9]
MspI_expdifreq<-diprobs[,10]
NcoI_expdifreq<-diprobs[,11]
NgoMIV_expdifreq<-diprobs[,12]
NotI_expdifreq<-diprobs[,13]
NsiI_expdifreq<-diprobs[,14]
NspI_expdifreq<-diprobs[,15]
PciI_expdifreq<-diprobs[,16]
PstI_expdifreq<-diprobs[,17]
SbfI_expdifreq<-diprobs[,18]
SgrAI_expdifreq<-diprobs[,19]
expdifreq<-cbind(ids,AgeI_expdifreq,ApoI_expdifreq,BsrFI_expdifreq,EcoRI_expdifreq,FatI_expdifreq,KpnI_expdifreq,MluCI_expdifreq,MseI_expdifreq,MspI_expdifreq,NcoI_expdifreq,NgoMIV_expdifreq,NotI_expdifreq,NsiI_expdifreq,NspI_expdifreq,PciI_expdifreq,PstI_expdifreq,SbfI_expdifreq,SgrAI_expdifreq)
colnames(expdifreq) <- c("id","AgeI","ApoI","BsrFI","EcoRI","FatI","KpnI","MluCI","MseI","MspI","NcoI","NgoMIV","NotI","NsiI","NspI","PciI","PstI","SbfI","SgrAI")
row.names(expdifreq) <- expdifreq$id
expdifreq_frame<- expdifreq[,2:19]
expdifreq_matrix <- data.matrix(expdifreq_frame)

#Calculate expected number of cut sites with each enzyme for di
AgeI_expdi<-AgeI_expdifreq*genome_size
ApoI_expdi<-ApoI_expdifreq*genome_size
BsrFI_expdi<-BsrFI_expdifreq*genome_size
EcoRI_expdi<-EcoRI_expdifreq*genome_size
FatI_expdi<-FatI_expdifreq*genome_size
KpnI_expdi<-KpnI_expdifreq*genome_size
MluCI_expdi<-MluCI_expdifreq*genome_size
MseI_expdi<-MseI_expdifreq*genome_size
MspI_expdi<-MspI_expdifreq*genome_size
NcoI_expdi<-NcoI_expdifreq*genome_size
NgoMIV_expdi<-NgoMIV_expdifreq*genome_size
NotI_expdi<-NotI_expdifreq*genome_size
NsiI_expdi<-NsiI_expdifreq*genome_size
NspI_expdi<-NspI_expdifreq*genome_size
PciI_expdi<-PciI_expdifreq*genome_size
PstI_expdi<-PstI_expdifreq*genome_size
SbfI_expdi<-SbfI_expdifreq*genome_size
SgrAI_expdi<-SgrAI_expdifreq*genome_size
expected<-cbind(ids,AgeI_expdi,ApoI_expdi,BsrFI_expdi,EcoRI_expdi,FatI_expdi,KpnI_expdi,MluCI_expdi,MseI_expdi,MspI_expdi,NcoI_expdi,NgoMIV_expdi,NotI_expdi,NsiI_expdi,NspI_expdi,PciI_expdi,PstI_expdi,SbfI_expdi,SgrAI_expdi)
colnames(expected) <- c("id","AgeI","ApoI","BsrFI","EcoRI","FatI","KpnI","MluCI","MseI","MspI","NcoI","NgoMIV","NotI","NsiI","NspI","PciI","PstI","SbfI","SgrAI")
row.names(expected) <- expected$id
expdi_frame<- expected[,2:19]
expdi_matrix <- data.matrix(expdi_frame)

#Expected frequency of cut sites with each enzyme for tri
AgeI_exptrifreq<-triprobs[,2]
ApoI_exptrifreq<-triprobs[,3]
BsrFI_exptrifreq<-triprobs[,4]
EcoRI_exptrifreq<-triprobs[,5]
FatI_exptrifreq<-triprobs[,6]
KpnI_exptrifreq<-triprobs[,7]
MluCI_exptrifreq<-triprobs[,8]
MseI_exptrifreq<-triprobs[,9]
MspI_exptrifreq<-triprobs[,10]
NcoI_exptrifreq<-triprobs[,11]
NgoMIV_exptrifreq<-triprobs[,12]
NotI_exptrifreq<-triprobs[,13]
NsiI_exptrifreq<-triprobs[,14]
NspI_exptrifreq<-triprobs[,15]
PciI_exptrifreq<-triprobs[,16]
PstI_exptrifreq<-triprobs[,17]
SbfI_exptrifreq<-triprobs[,18]
SgrAI_exptrifreq<-triprobs[,19]
exptrifreq<-cbind(ids,AgeI_exptrifreq,ApoI_exptrifreq,BsrFI_exptrifreq,EcoRI_exptrifreq,FatI_exptrifreq,KpnI_exptrifreq,MluCI_exptrifreq,MseI_exptrifreq,MspI_exptrifreq,NcoI_exptrifreq,NgoMIV_exptrifreq,NotI_exptrifreq,NsiI_exptrifreq,NspI_exptrifreq,PciI_exptrifreq,PstI_exptrifreq,SbfI_exptrifreq,SgrAI_exptrifreq)
colnames(exptrifreq) <- c("id","AgeI","ApoI","BsrFI","EcoRI","FatI","KpnI","MluCI","MseI","MspI","NcoI","NgoMIV","NotI","NsiI","NspI","PciI","PstI","SbfI","SgrAI")
row.names(exptrifreq) <- exptrifreq$id
exptrifreq_frame<- exptrifreq[,2:19]
exptrifreq_matrix <- data.matrix(exptrifreq_frame)

#Calculate expected number of cut sites with each enzyme for tri
AgeI_exptri<-AgeI_exptrifreq*genome_size
ApoI_exptri<-ApoI_exptrifreq*genome_size
BsrFI_exptri<-BsrFI_exptrifreq*genome_size
EcoRI_exptri<-EcoRI_exptrifreq*genome_size
FatI_exptri<-FatI_exptrifreq*genome_size
KpnI_exptri<-KpnI_exptrifreq*genome_size
MluCI_exptri<-MluCI_exptrifreq*genome_size
MseI_exptri<-MseI_exptrifreq*genome_size
MspI_exptri<-MspI_exptrifreq*genome_size
NcoI_exptri<-NcoI_exptrifreq*genome_size
NgoMIV_exptri<-NgoMIV_exptrifreq*genome_size
NotI_exptri<-NotI_exptrifreq*genome_size
NsiI_exptri<-NsiI_exptrifreq*genome_size
NspI_exptri<-NspI_exptrifreq*genome_size
PciI_exptri<-PciI_exptrifreq*genome_size
PstI_exptri<-PstI_exptrifreq*genome_size
SbfI_exptri<-SbfI_exptrifreq*genome_size
SgrAI_exptri<-SgrAI_exptrifreq*genome_size
exptriected<-cbind(ids,AgeI_exptri,ApoI_exptri,BsrFI_exptri,EcoRI_exptri,FatI_exptri,KpnI_exptri,MluCI_exptri,MseI_exptri,MspI_exptri,NcoI_exptri,NgoMIV_exptri,NotI_exptri,NsiI_exptri,NspI_exptri,PciI_exptri,PstI_exptri,SbfI_exptri,SgrAI_exptri)
colnames(exptriected) <- c("id","AgeI","ApoI","BsrFI","EcoRI","FatI","KpnI","MluCI","MseI","MspI","NcoI","NgoMIV","NotI","NsiI","NspI","PciI","PstI","SbfI","SgrAI")
row.names(exptriected) <- exptriected$id
exptri_frame<- exptriected[,2:19]
exptri_matrix <- data.matrix(exptri_frame)

#Extract the observed number of cut sites with each enzyme
AgeI_obs<-counts[,2]
ApoI_obs<-counts[,3]
BsrFI_obs<-counts[,4]
EcoRI_obs<-counts[,5]
FatI_obs<-counts[,6]
KpnI_obs<-counts[,7]
MluCI_obs<-counts[,8]
MseI_obs<-counts[,9]
MspI_obs<-counts[,10]
NcoI_obs<-counts[,11]
NgoMIV_obs<-counts[,12]
NotI_obs<-counts[,13]
NsiI_obs<-counts[,14]
NspI_obs<-counts[,15]
PciI_obs<-counts[,16]
PstI_obs<-counts[,17]
SbfI_obs<-counts[,18]
SgrAI_obs<-counts[,19]
observed<-cbind(ids,AgeI_obs,ApoI_obs,BsrFI_obs,EcoRI_obs,FatI_obs,KpnI_obs,MluCI_obs,MseI_obs,MspI_obs,NcoI_obs,NgoMIV_obs,NotI_obs,NsiI_obs,NspI_obs,PciI_obs,PstI_obs,SbfI_obs,SgrAI_obs)
colnames(observed) <- c("id","AgeI","ApoI","BsrFI","EcoRI","FatI","KpnI","MluCI","MseI","MspI","NcoI","NgoMIV","NotI","NsiI","NspI","PciI","PstI","SbfI","SgrAI")
row.names(observed) <- observed$id
obs_frame<- observed[,2:19]
obs_matrix <- data.matrix(obs_frame)

#Calculate observed frequency of cut sites with each enzyme
AgeI_obsfreq<-AgeI_obs/genome_size
ApoI_obsfreq<-ApoI_obs/genome_size
BsrFI_obsfreq<-BsrFI_obs/genome_size
EcoRI_obsfreq<-EcoRI_obs/genome_size
FatI_obsfreq<-FatI_obs/genome_size
KpnI_obsfreq<-KpnI_obs/genome_size
MluCI_obsfreq<-MluCI_obs/genome_size
MseI_obsfreq<-MseI_obs/genome_size
MspI_obsfreq<-MspI_obs/genome_size
NcoI_obsfreq<-NcoI_obs/genome_size
NgoMIV_obsfreq<-NgoMIV_obs/genome_size
NotI_obsfreq<-NotI_obs/genome_size
NsiI_obsfreq<-NsiI_obs/genome_size
NspI_obsfreq<-NspI_obs/genome_size
PciI_obsfreq<-PciI_obs/genome_size
PstI_obsfreq<-PstI_obs/genome_size
SbfI_obsfreq<-SbfI_obs/genome_size
SgrAI_obsfreq<-SgrAI_obs/genome_size
obsfreq<-cbind(ids,AgeI_obsfreq,ApoI_obsfreq,BsrFI_obsfreq,EcoRI_obsfreq,FatI_obsfreq,KpnI_obsfreq,MluCI_obsfreq,MseI_obsfreq,MspI_obsfreq,NcoI_obsfreq,NgoMIV_obsfreq,NotI_obsfreq,NsiI_obsfreq,NspI_obsfreq,PciI_obsfreq,PstI_obsfreq,SbfI_obsfreq,SgrAI_obsfreq)
colnames(obsfreq) <- c("id","AgeI","ApoI","BsrFI","EcoRI","FatI","KpnI","MluCI","MseI","MspI","NcoI","NgoMIV","NotI","NsiI","NspI","PciI","PstI","SbfI","SgrAI")
row.names(obsfreq) <- obsfreq$id
obsfreq_frame<- obsfreq[,2:19]
obsfreq_matrix <- data.matrix(obsfreq_frame)
obsfreq_matrix_perMb <- obsfreq_matrix*10^6
log_obsfreq_matrix_perMb <- log10(obsfreq_matrix_perMb)
## Replace -Inf with 0.0
log_obsfreq_matrix_perMb[273,12]<-0.0
log_obsfreq_matrix_perMb[274,12]<-0.0
log_obsfreq_matrix_perMb[433,12]<-0.0


#Calculate similarity index for each enzyme ((observed-expgcected)/expgcected) for gc
AgeI_percgc<-(AgeI_obsfreq-AgeI_expgcfreq)/AgeI_expgcfreq
ApoI_percgc<-(ApoI_obsfreq-ApoI_expgcfreq)/ApoI_expgcfreq
BsrFI_percgc<-(BsrFI_obsfreq-BsrFI_expgcfreq)/BsrFI_expgcfreq
EcoRI_percgc<-(EcoRI_obsfreq-EcoRI_expgcfreq)/EcoRI_expgcfreq
FatI_percgc<-(FatI_obsfreq-FatI_expgcfreq)/FatI_expgcfreq
KpnI_percgc<-(KpnI_obsfreq-KpnI_expgcfreq)/KpnI_expgcfreq
MluCI_percgc<-(MluCI_obsfreq-MluCI_expgcfreq)/MluCI_expgcfreq
MseI_percgc<-(MseI_obsfreq-MseI_expgcfreq)/MseI_expgcfreq
MspI_percgc<-(MspI_obsfreq-MspI_expgcfreq)/MspI_expgcfreq
NcoI_percgc<-(NcoI_obsfreq-NcoI_expgcfreq)/NcoI_expgcfreq
NgoMIV_percgc<-(NgoMIV_obsfreq-NgoMIV_expgcfreq)/NgoMIV_expgcfreq
NotI_percgc<-(NotI_obsfreq-NotI_expgcfreq)/NotI_expgcfreq
NsiI_percgc<-(NsiI_obsfreq-NsiI_expgcfreq)/NsiI_expgcfreq
NspI_percgc<-(NspI_obsfreq-NspI_expgcfreq)/NspI_expgcfreq
PciI_percgc<-(PciI_obsfreq-PciI_expgcfreq)/PciI_expgcfreq
PstI_percgc<-(PstI_obsfreq-PstI_expgcfreq)/PstI_expgcfreq
SbfI_percgc<-(SbfI_obsfreq-SbfI_expgcfreq)/SbfI_expgcfreq
SgrAI_percgc<-(SgrAI_obsfreq-SgrAI_expgcfreq)/SgrAI_expgcfreq
percgc<-cbind(ids,AgeI_percgc,ApoI_percgc,BsrFI_percgc,EcoRI_percgc,FatI_percgc,KpnI_percgc,MluCI_percgc,MseI_percgc,MspI_percgc,NcoI_percgc,NgoMIV_percgc,NotI_percgc,NsiI_percgc,NspI_percgc,PciI_percgc,PstI_percgc,SbfI_percgc,SgrAI_percgc)
colnames(percgc) <- c("id","AgeI","ApoI","BsrFI","EcoRI","FatI","KpnI","MluCI","MseI","MspI","NcoI","NgoMIV","NotI","NsiI","NspI","PciI","PstI","SbfI","SgrAI")
row.names(percgc) <- percgc$id
percgc_frame<- percgc[,2:19]
percgc_matrix <- data.matrix(percgc_frame)

#Calculate similarity index for each enzyme ((observed-expected)/expected) for mono
AgeI_percmono<-(AgeI_obsfreq-AgeI_expmonofreq)/AgeI_expmonofreq
ApoI_percmono<-(ApoI_obsfreq-ApoI_expmonofreq)/ApoI_expmonofreq
BsrFI_percmono<-(BsrFI_obsfreq-BsrFI_expmonofreq)/BsrFI_expmonofreq
EcoRI_percmono<-(EcoRI_obsfreq-EcoRI_expmonofreq)/EcoRI_expmonofreq
FatI_percmono<-(FatI_obsfreq-FatI_expmonofreq)/FatI_expmonofreq
KpnI_percmono<-(KpnI_obsfreq-KpnI_expmonofreq)/KpnI_expmonofreq
MluCI_percmono<-(MluCI_obsfreq-MluCI_expmonofreq)/MluCI_expmonofreq
MseI_percmono<-(MseI_obsfreq-MseI_expmonofreq)/MseI_expmonofreq
MspI_percmono<-(MspI_obsfreq-MspI_expmonofreq)/MspI_expmonofreq
NcoI_percmono<-(NcoI_obsfreq-NcoI_expmonofreq)/NcoI_expmonofreq
NgoMIV_percmono<-(NgoMIV_obsfreq-NgoMIV_expmonofreq)/NgoMIV_expmonofreq
NotI_percmono<-(NotI_obsfreq-NotI_expmonofreq)/NotI_expmonofreq
NsiI_percmono<-(NsiI_obsfreq-NsiI_expmonofreq)/NsiI_expmonofreq
NspI_percmono<-(NspI_obsfreq-NspI_expmonofreq)/NspI_expmonofreq
PciI_percmono<-(PciI_obsfreq-PciI_expmonofreq)/PciI_expmonofreq
PstI_percmono<-(PstI_obsfreq-PstI_expmonofreq)/PstI_expmonofreq
SbfI_percmono<-(SbfI_obsfreq-SbfI_expmonofreq)/SbfI_expmonofreq
SgrAI_percmono<-(SgrAI_obsfreq-SgrAI_expmonofreq)/SgrAI_expmonofreq
percmono<-cbind(ids,AgeI_percmono,ApoI_percmono,BsrFI_percmono,EcoRI_percmono,FatI_percmono,KpnI_percmono,MluCI_percmono,MseI_percmono,MspI_percmono,NcoI_percmono,NgoMIV_percmono,NotI_percmono,NsiI_percmono,NspI_percmono,PciI_percmono,PstI_percmono,SbfI_percmono,SgrAI_percmono)
colnames(percmono) <- c("id","AgeI","ApoI","BsrFI","EcoRI","FatI","KpnI","MluCI","MseI","MspI","NcoI","NgoMIV","NotI","NsiI","NspI","PciI","PstI","SbfI","SgrAI")
row.names(percmono) <- percmono$id
percmono_frame<- percmono[,2:19]
percmono_matrix <- data.matrix(percmono_frame)

#Calculate similarity index for each enzyme ((observed-expected)/expected) for di
AgeI_percdi<-(AgeI_obsfreq-AgeI_expdifreq)/AgeI_expdifreq
ApoI_percdi<-(ApoI_obsfreq-ApoI_expdifreq)/ApoI_expdifreq
BsrFI_percdi<-(BsrFI_obsfreq-BsrFI_expdifreq)/BsrFI_expdifreq
EcoRI_percdi<-(EcoRI_obsfreq-EcoRI_expdifreq)/EcoRI_expdifreq
FatI_percdi<-(FatI_obsfreq-FatI_expdifreq)/FatI_expdifreq
KpnI_percdi<-(KpnI_obsfreq-KpnI_expdifreq)/KpnI_expdifreq
MluCI_percdi<-(MluCI_obsfreq-MluCI_expdifreq)/MluCI_expdifreq
MseI_percdi<-(MseI_obsfreq-MseI_expdifreq)/MseI_expdifreq
MspI_percdi<-(MspI_obsfreq-MspI_expdifreq)/MspI_expdifreq
NcoI_percdi<-(NcoI_obsfreq-NcoI_expdifreq)/NcoI_expdifreq
NgoMIV_percdi<-(NgoMIV_obsfreq-NgoMIV_expdifreq)/NgoMIV_expdifreq
NotI_percdi<-(NotI_obsfreq-NotI_expdifreq)/NotI_expdifreq
NsiI_percdi<-(NsiI_obsfreq-NsiI_expdifreq)/NsiI_expdifreq
NspI_percdi<-(NspI_obsfreq-NspI_expdifreq)/NspI_expdifreq
PciI_percdi<-(PciI_obsfreq-PciI_expdifreq)/PciI_expdifreq
PstI_percdi<-(PstI_obsfreq-PstI_expdifreq)/PstI_expdifreq
SbfI_percdi<-(SbfI_obsfreq-SbfI_expdifreq)/SbfI_expdifreq
SgrAI_percdi<-(SgrAI_obsfreq-SgrAI_expdifreq)/SgrAI_expdifreq
percdi<-cbind(ids,AgeI_percdi,ApoI_percdi,BsrFI_percdi,EcoRI_percdi,FatI_percdi,KpnI_percdi,MluCI_percdi,MseI_percdi,MspI_percdi,NcoI_percdi,NgoMIV_percdi,NotI_percdi,NsiI_percdi,NspI_percdi,PciI_percdi,PstI_percdi,SbfI_percdi,SgrAI_percdi)
colnames(percdi) <- c("id","AgeI","ApoI","BsrFI","EcoRI","FatI","KpnI","MluCI","MseI","MspI","NcoI","NgoMIV","NotI","NsiI","NspI","PciI","PstI","SbfI","SgrAI")
row.names(percdi) <- percdi$id
percdi_frame<- percdi[,2:19]
percdi_matrix <- data.matrix(percdi_frame)

#Calculate similarity index for each enzyme ((observed-exptriected)/exptriected) for tri
AgeI_perctri<-(AgeI_obsfreq-AgeI_exptrifreq)/AgeI_exptrifreq
ApoI_perctri<-(ApoI_obsfreq-ApoI_exptrifreq)/ApoI_exptrifreq
BsrFI_perctri<-(BsrFI_obsfreq-BsrFI_exptrifreq)/BsrFI_exptrifreq
EcoRI_perctri<-(EcoRI_obsfreq-EcoRI_exptrifreq)/EcoRI_exptrifreq
FatI_perctri<-(FatI_obsfreq-FatI_exptrifreq)/FatI_exptrifreq
KpnI_perctri<-(KpnI_obsfreq-KpnI_exptrifreq)/KpnI_exptrifreq
MluCI_perctri<-(MluCI_obsfreq-MluCI_exptrifreq)/MluCI_exptrifreq
MseI_perctri<-(MseI_obsfreq-MseI_exptrifreq)/MseI_exptrifreq
MspI_perctri<-(MspI_obsfreq-MspI_exptrifreq)/MspI_exptrifreq
NcoI_perctri<-(NcoI_obsfreq-NcoI_exptrifreq)/NcoI_exptrifreq
NgoMIV_perctri<-(NgoMIV_obsfreq-NgoMIV_exptrifreq)/NgoMIV_exptrifreq
NotI_perctri<-(NotI_obsfreq-NotI_exptrifreq)/NotI_exptrifreq
NsiI_perctri<-(NsiI_obsfreq-NsiI_exptrifreq)/NsiI_exptrifreq
NspI_perctri<-(NspI_obsfreq-NspI_exptrifreq)/NspI_exptrifreq
PciI_perctri<-(PciI_obsfreq-PciI_exptrifreq)/PciI_exptrifreq
PstI_perctri<-(PstI_obsfreq-PstI_exptrifreq)/PstI_exptrifreq
SbfI_perctri<-(SbfI_obsfreq-SbfI_exptrifreq)/SbfI_exptrifreq
SgrAI_perctri<-(SgrAI_obsfreq-SgrAI_exptrifreq)/SgrAI_exptrifreq
perctri<-cbind(ids,AgeI_perctri,ApoI_perctri,BsrFI_perctri,EcoRI_perctri,FatI_perctri,KpnI_perctri,MluCI_perctri,MseI_perctri,MspI_perctri,NcoI_perctri,NgoMIV_perctri,NotI_perctri,NsiI_perctri,NspI_perctri,PciI_perctri,PstI_perctri,SbfI_perctri,SgrAI_perctri)
colnames(perctri) <- c("id","AgeI","ApoI","BsrFI","EcoRI","FatI","KpnI","MluCI","MseI","MspI","NcoI","NgoMIV","NotI","NsiI","NspI","PciI","PstI","SbfI","SgrAI")
row.names(perctri) <- perctri$id
perctri_frame<- perctri[,2:19]
perctri_matrix <- data.matrix(perctri_frame)

###############examine genome mapping results
## extract the number of processed tags per genome,enzyme
AgeI_proc<-aligned[,2]+suppressed[,2]
ApoI_proc<-aligned[,3]+suppressed[,3]
BsrFI_proc<-aligned[,4]+suppressed[,4]
EcoRI_proc<-aligned[,5]+suppressed[,5]
FatI_proc<-aligned[,6]+suppressed[,6]
KpnI_proc<-aligned[,7]+suppressed[,7]
MluCI_proc<-aligned[,8]+suppressed[,8]
MseI_proc<-aligned[,9]+suppressed[,9]
MspI_proc<-aligned[,10]+suppressed[,10]
NcoI_proc<-aligned[,11]+suppressed[,11]
NgoMIV_proc<-aligned[,12]+suppressed[,12]
NotI_proc<-aligned[,13]+suppressed[,13]
NsiI_proc<-aligned[,14]+suppressed[,14]
NspI_proc<-aligned[,15]+suppressed[,15]
PciI_proc<-aligned[,16]+suppressed[,16]
PstI_proc<-aligned[,17]+suppressed[,17]
SbfI_proc<-aligned[,18]+suppressed[,18]
SgrAI_proc<-aligned[,19]+suppressed[,19]
processed<-cbind(ids,AgeI_proc,ApoI_proc,BsrFI_proc,EcoRI_proc,FatI_proc,KpnI_proc,MluCI_proc,MseI_proc,MspI_proc,NcoI_proc,NgoMIV_proc,NotI_proc,NsiI_proc,NspI_proc,PciI_proc,PstI_proc,SbfI_proc,SgrAI_proc)
colnames(processed) <- c("id","AgeI","ApoI","BsrFI","EcoRI","FatI","KpnI","MluCI","MseI","MspI","NcoI","NgoMIV","NotI","NsiI","NspI","PciI","PstI","SbfI","SgrAI")
row.names(processed) <- processed$id
processed_frame<- processed[,2:19]
processed_matrix <- data.matrix(processed_frame)

## extract the number of suppresed alingments (of tags) per genome,enzyme
AgeI_supp<-suppressed[,2]
ApoI_supp<-suppressed[,3]
BsrFI_supp<-suppressed[,4]
EcoRI_supp<-suppressed[,5]
FatI_supp<-suppressed[,6]
KpnI_supp<-suppressed[,7]
MluCI_supp<-suppressed[,8]
MseI_supp<-suppressed[,9]
MspI_supp<-suppressed[,10]
NcoI_supp<-suppressed[,11]
NgoMIV_supp<-suppressed[,12]
NotI_supp<-suppressed[,13]
NsiI_supp<-suppressed[,14]
NspI_supp<-suppressed[,15]
PciI_supp<-suppressed[,16]
PstI_supp<-suppressed[,17]
SbfI_supp<-suppressed[,18]
SgrAI_supp<-suppressed[,19]

### calculate the percentage of suppresed alignments (of tags) per genome,enzyme
AgeI_supp_perc<-AgeI_supp/AgeI_proc*100
ApoI_supp_perc<-ApoI_supp/ApoI_proc*100
BsrFI_supp_perc<-BsrFI_supp/BsrFI_proc*100
EcoRI_supp_perc<-EcoRI_supp/EcoRI_proc*100
FatI_supp_perc<-FatI_supp/FatI_proc*100
KpnI_supp_perc<-KpnI_supp/KpnI_proc*100
MluCI_supp_perc<-MluCI_supp/MluCI_proc*100
MseI_supp_perc<-MseI_supp/MseI_proc*100
MspI_supp_perc<-MspI_supp/MspI_proc*100
NcoI_supp_perc<-NcoI_supp/NcoI_proc*100
NgoMIV_supp_perc<-NgoMIV_supp/NgoMIV_proc*100
NotI_supp_perc<-NotI_supp/NotI_proc*100
NsiI_supp_perc<-NsiI_supp/NsiI_proc*100
NspI_supp_perc<-NspI_supp/NspI_proc*100
PciI_supp_perc<-PciI_supp/PciI_proc*100
PstI_supp_perc<-PstI_supp/PstI_proc*100
SbfI_supp_perc<-SbfI_supp/SbfI_proc*100
SgrAI_supp_perc<-SgrAI_supp/SgrAI_proc*100
supp_perc<-cbind(ids,AgeI_supp_perc,ApoI_supp_perc,BsrFI_supp_perc,EcoRI_supp_perc,FatI_supp_perc,KpnI_supp_perc,MluCI_supp_perc,MseI_supp_perc,MspI_supp_perc,NcoI_supp_perc,NgoMIV_supp_perc,NotI_supp_perc,NsiI_supp_perc,NspI_supp_perc,PciI_supp_perc,PstI_supp_perc,SbfI_supp_perc,SgrAI_supp_perc)
colnames(supp_perc) <- c("id","AgeI","ApoI","BsrFI","EcoRI","FatI","KpnI","MluCI","MseI","MspI","NcoI","NgoMIV","NotI","NsiI","NspI","PciI","PstI","SbfI","SgrAI")
row.names(supp_perc) <- supp_perc$id
supp_perc_frame<- supp_perc[,2:19]
supp_perc_matrix <- data.matrix(supp_perc_frame)
supp_perc_matrix[is.nan(supp_perc_matrix)] <- 0
supp_perc_matrix[is.infinite(supp_perc_matrix)] <- 0


############################################################################################


#### output files

write.table(gc_matrix, "~/Documents/RADtag/Genomes/paper/scripts/gc_matrix.txt", sep="\t")
write.table(percgc_matrix, "~/Documents/RADtag/Genomes/paper/scripts/percgc_matrix.txt", sep="\t")
write.table(percmono_matrix, "~/Documents/RADtag/Genomes/paper/scripts/percmono_matrix.txt", sep="\t")
write.table(percdi_matrix, "~/Documents/RADtag/Genomes/paper/scripts/percdi_matrix.txt", sep="\t")
write.table(perctri_matrix, "~/Documents/RADtag/Genomes/paper/scripts/perctri_matrix.txt", sep="\t")
write.table(processed_matrix, "~/Documents/RADtag/Genomes/paper/scripts/processed_matrix.txt", sep="\t")
write.table(supp_perc_matrix, "~/Documents/RADtag/Genomes/paper/scripts/supp_perc_matrix.txt", sep="\t")
write.table(supp_perc_matrix, "~/Documents/RADtag/Genomes/paper/scripts/supp_perc_matrix.txt", sep="\t")
write.table(expgc_matrix, "~/Documents/RADtag/Genomes/paper/scripts/expgc_matrix.txt", sep="\t")
write.table(expmono_matrix, "~/Documents/RADtag/Genomes/paper/scripts/expmono_matrix.txt", sep="\t")
write.table(expdi_matrix, "~/Documents/RADtag/Genomes/paper/scripts/expdi_matrix.txt", sep="\t")
write.table(exptri_matrix, "~/Documents/RADtag/Genomes/paper/scripts/exptri_matrix.txt", sep="\t")
write.table(expgcfreq_matrix, "~/Documents/RADtag/Genomes/paper/scripts/expgcfreq_matrix.txt", sep="\t")
write.table(expmonofreq_matrix, "~/Documents/RADtag/Genomes/paper/scripts/expmonofreq_matrix.txt", sep="\t")
write.table(expdifreq_matrix, "~/Documents/RADtag/Genomes/paper/scripts/expdifreq_matrix.txt", sep="\t")
write.table(exptrifreq_matrix, "~/Documents/RADtag/Genomes/paper/scripts/exptrifreq_matrix.txt", sep="\t")
write.table(obs_matrix, "~/Documents/RADtag/Genomes/paper/scripts/obs_matrix.txt", sep="\t")
write.table(obsfreq_matrix_perMb, "~/Documents/RADtag/Genomes/paper/scripts/obsfreq_matrix_perMb.txt", sep="\t")
write.table(log_obsfreq_matrix_perMb, "~/Documents/RADtag/Genomes/paper/scripts/log_obsfreq_matrix_perMb.txt", sep="\t")
write.table(trinuc_gammastar_matrix, "~/Documents/RADtag/Genomes/paper/scripts/trinuc_gammastar_matrix.txt", sep="\t")
write.table(dinuc_rhostar_matrix, "~/Documents/RADtag/Genomes/paper/scripts/dinuc_rhostar_matrix.txt", sep="\t")

############PLOTTING
##### FOR PAPER #####

#Figure 1.  Observed restriction site frequencies. Heatmap of the observed frequency of restriction sites. 
log_obsfreq_matrix_perMb_reord<-cbind(log_obsfreq_matrix_perMb[,12],log_obsfreq_matrix_perMb[,18],log_obsfreq_matrix_perMb[,3],log_obsfreq_matrix_perMb[,11],log_obsfreq_matrix_perMb[,1],log_obsfreq_matrix_perMb[,9],log_obsfreq_matrix_perMb[,17],log_obsfreq_matrix_perMb[,16],log_obsfreq_matrix_perMb[,13],log_obsfreq_matrix_perMb[,2],log_obsfreq_matrix_perMb[,4],log_obsfreq_matrix_perMb[,7],log_obsfreq_matrix_perMb[,8],log_obsfreq_matrix_perMb[,14],log_obsfreq_matrix_perMb[,10],log_obsfreq_matrix_perMb[,15],log_obsfreq_matrix_perMb[,5],log_obsfreq_matrix_perMb[,6])
colnames(log_obsfreq_matrix_perMb_reord) <- c("NotI","SgrAI","BsrFI","NgoMIV","AgeI","MspI","SbfI","PstI","NsiI","ApoI","EcoRI","MluCI","MseI","NspI","NcoI","PciI","FatI","KpnI")

hv_log_obsfreq_perMb<-heatmap.2((log_obsfreq_matrix_perMb_reord), Rowv=NA, Colv=NA,scale="none", col= colorpanel(41,"blue","yellow","red"), main="Log (observed # restriction sites per Mb)", trace="none", denscol="grey")


#Figure 2. Center: heatmap of the ρ_XY^* odds ratio values.
hv<-heatmap.2((dinuc_rhostar_matrix), Rowv=NA, Colv=NA, scale="none", col= colorpanel(42,"green","#000000","red"), main="dinucleotide rho*", trace="none")


#Figure 2. Right: heatmap of the ρ_XY^* odds ratio significant values ρ_XY^*<0.78 and ρ_XY^*>1.23.
ditri_matrix_breaks<-c(0,0.78,1.23,2)
hv<-heatmap.2((dinuc_rhostar_matrix), Rowv=NA, Colv=NA, scale="none", breaks=ditri_matrix_breaks, col= colorpanel(3,"green","#000000","red"), main="dinucleotide rho* significant", trace="none")


#Figure 3. Right: heatmap of the γ_XYZ^* odds ratio values. 
gamma_breaks<-c(0.1430325,0.1854307,0.2278288,0.2702270,0.3126251,0.3550232,0.3974214,0.4398195,0.4822177,0.5246158,0.5670139,0.6094121,0.6518102,0.6942084,0.7366065,0.7790046,0.8214028,0.8638009,0.9061990,0.9485972,0.9909953,1.0333935,1.0757916,1.1181897,1.1605879,1.2029860,1.2453842,1.2877823,1.3301804,1.3725786,1.4149767,1.4573749,1.4997730,1.5421711,1.5845693,1.6269674,1.6693655,1.7117637,1.7541618,1.7965600,1.8389581,1.8813562,1.9237544)
hv<-heatmap.2((trinuc_gammastar_matrix), Rowv=NA, scale="none", col= colorpanel(42,"green","#000000","red"), breaks=gamma_breaks,main="trinucleotide gamma*", trace="none")


#Figure 4. Right: heatmap of the γ_XYZ^* odds ratio significant values ρ_XY^*<0.78 and ρ_XY^*>1.23. 
ditri_matrix_breaks<-c(0,0.78,1.23,2)
hv<-heatmap.2((trinuc_gammastar_matrix), Rowv=NA, scale="none", breaks=ditri_matrix_breaks, col= colorpanel(3,"green","#000000","red"), main="trinucleotide gamma* significant", trace="none")

#Figure 5. Box and whisker plots of the similarity index (SI) for each species per enzyme. 
AgeI_SI<-cbind(percgc[,2],percmono[,2],percdi[,2],perctri[,2])
colnames(AgeI_SI)<-c("AgeI_gc","AgeI_mono","AgeI_di","AgeI_tri")
pdf(file='~/Documents/RADtag/Genomes/paper/Figures/AgeI_SI.pdf')
boxplot(AgeI_SI,col="orange",ylab="SI",pars=list(pch=20))
dev.off()
ApoI_SI<-cbind(percgc[,3],percmono[,3],percdi[,3],perctri[,3])
colnames(ApoI_SI)<-c("ApoI_gc","ApoI_mono","ApoI_di","ApoI_tri")
pdf(file='~/Documents/RADtag/Genomes/paper/Figures/ApoI_SI.pdf')
boxplot(ApoI_SI,col="orange",ylab="SI",pars=list(pch=20))
dev.off()
BsrFI_SI<-cbind(percgc[,4],percmono[,4],percdi[,4],perctri[,4])
colnames(BsrFI_SI)<-c("BsrFI_gc","BsrFI_mono","BsrFI_di","BsrFI_tri")
pdf(file='~/Documents/RADtag/Genomes/paper/Figures/BsrFI_SI.pdf')
boxplot(BsrFI_SI,col="orange",ylab="SI",pars=list(pch=20))
dev.off()
EcoRI_SI<-cbind(percgc[,5],percmono[,5],percdi[,5],perctri[,5])
colnames(EcoRI_SI)<-c("EcoRI_gc","EcoRI_mono","EcoRI_di","EcoRI_tri")
pdf(file='~/Documents/RADtag/Genomes/paper/Figures/EcoRI_SI.pdf')
boxplot(EcoRI_SI,col="orange",ylab="SI",pars=list(pch=20))
dev.off()
FatI_SI<-cbind(percgc[,6],percmono[,6],percdi[,6],perctri[,6])
colnames(FatI_SI)<-c("FatI_gc","FatI_mono","FatI_di","FatI_tri")
pdf(file='~/Documents/RADtag/Genomes/paper/Figures/FatI_SI.pdf')
boxplot(FatI_SI,col="orange",ylab="SI",pars=list(pch=20))
dev.off()
KpnI_SI<-cbind(percgc[,7],percmono[,7],percdi[,7],perctri[,7])
colnames(KpnI_SI)<-c("KpnI_gc","KpnI_mono","KpnI_di","KpnI_tri")
pdf(file='~/Documents/RADtag/Genomes/paper/Figures/KpnI_SI.pdf')
boxplot(KpnI_SI,col="orange",ylab="SI",pars=list(pch=20))
dev.off()
MluCI_SI<-cbind(percgc[,8],percmono[,8],percdi[,8],perctri[,8])
colnames(MluCI_SI)<-c("MluCI_gc","MluCI_mono","MluCI_di","MluCI_tri")
pdf(file='~/Documents/RADtag/Genomes/paper/Figures/MluCI_SI.pdf')
boxplot(MluCI_SI,col="orange",ylab="SI",pars=list(pch=20))
dev.off()
MseI_SI<-cbind(percgc[,9],percmono[,9],percdi[,9],perctri[,9])
colnames(MseI_SI)<-c("MseI_gc","MseI_mono","MseI_di","MseI_tri")
pdf(file='~/Documents/RADtag/Genomes/paper/Figures/MseI_SI.pdf')
boxplot(MseI_SI,col="orange",ylab="SI",pars=list(pch=20))
dev.off()
MspI_SI<-cbind(percgc[,10],percmono[,10],percdi[,10],perctri[,10])
colnames(MspI_SI)<-c("MspI_gc","MspI_mono","MspI_di","MspI_tri")
pdf(file='~/Documents/RADtag/Genomes/paper/Figures/MspI_SI.pdf')
boxplot(MspI_SI,col="orange",ylab="SI",pars=list(pch=20))
dev.off()
NcoI_SI<-cbind(percgc[,11],percmono[,11],percdi[,11],perctri[,11])
colnames(NcoI_SI)<-c("NcoI_gc","NcoI_mono","NcoI_di","NcoI_tri")
pdf(file='~/Documents/RADtag/Genomes/paper/Figures/NcoI_SI.pdf')
boxplot(NcoI_SI,col="orange",ylab="SI",pars=list(pch=20))
dev.off()
NgoMIV_SI<-cbind(percgc[,12],percmono[,12],percdi[,12],perctri[,12])
colnames(NgoMIV_SI)<-c("NgoMIV_gc","NgoMIV_mono","NgoMIV_di","NgoMIV_tri")
pdf(file='~/Documents/RADtag/Genomes/paper/Figures/NgoMIV_SI.pdf')
boxplot(NgoMIV_SI,col="orange",ylab="SI",pars=list(pch=20))
dev.off()
NotI_SI<-cbind(percgc[,13],percmono[,13],percdi[,13],perctri[,13])
#Exclude outlier Enthis
NotI_SI<-NotI_SI[-37,]
colnames(NotI_SI)<-c("NotI_gc","NotI_mono","NotI_di","NotI_tri")
pdf(file='~/Documents/RADtag/Genomes/paper/Figures/NotI_SI.pdf')
boxplot(NotI_SI,col="orange",ylab="SI",pars=list(pch=20))
dev.off()
NsiI_SI<-cbind(percgc[,14],percmono[,14],percdi[,14],perctri[,14])
colnames(NsiI_SI)<-c("NsiI_gc","NsiI_mono","NsiI_di","NsiI_tri")
pdf(file='~/Documents/RADtag/Genomes/paper/Figures/NsiI_SI.pdf')
boxplot(NsiI_SI,col="orange",ylab="SI",pars=list(pch=20))
dev.off()
NspI_SI<-cbind(percgc[,15],percmono[,15],percdi[,15],perctri[,15])
colnames(NspI_SI)<-c("NspI_gc","NspI_mono","NspI_di","NspI_tri")
pdf(file='~/Documents/RADtag/Genomes/paper/Figures/NspI_SI.pdf')
boxplot(NspI_SI,col="orange",ylab="SI",pars=list(pch=20))
dev.off()
PciI_SI<-cbind(percgc[,16],percmono[,16],percdi[,16],perctri[,16])
colnames(PciI_SI)<-c("PciI_gc","PciI_mono","PciI_di","PciI_tri")
pdf(file='~/Documents/RADtag/Genomes/paper/Figures/PciI_SI.pdf')
boxplot(PciI_SI,col="orange",ylab="SI",pars=list(pch=20))
dev.off()
PstI_SI<-cbind(percgc[,17],percmono[,17],percdi[,17],perctri[,17])
colnames(PstI_SI)<-c("PstI_gc","PstI_mono","PstI_di","PstI_tri")
pdf(file='~/Documents/RADtag/Genomes/paper/Figures/PstI_SI.pdf')
boxplot(PstI_SI,col="orange",ylab="SI",pars=list(pch=20))
dev.off()
SbfI_SI<-cbind(percgc[,18],percmono[,18],percdi[,18],perctri[,18])
colnames(SbfI_SI)<-c("SbfI_gc","SbfI_mono","SbfI_di","SbfI_tri")
pdf(file='~/Documents/RADtag/Genomes/paper/Figures/SbfI_SI.pdf')
boxplot(SbfI_SI,col="orange",ylab="SI",pars=list(pch=20))
dev.off()
SgrAI_SI<-cbind(percgc[,19],percmono[,19],percdi[,19],perctri[,19])
colnames(SgrAI_SI)<-c("SgrAI_gc","SgrAI_mono","SgrAI_di","SgrAI_tri")
pdf(file='~/Documents/RADtag/Genomes/paper/Figures/SgrAI_SI.pdf')
boxplot(SgrAI_SI,col="orange",ylab="SI",pars=list(pch=20))
dev.off()


#Figure 6. Center: heatmap of the similarity indexes for the dinucleotide model. 
matrix_breaks<-c(-1.01,-1,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1)

percdi_matrix_reord<-cbind(percdi_matrix[,12],percdi_matrix[,18],percdi_matrix[,3],percdi_matrix[,11],percdi_matrix[,1],percdi_matrix[,9],percdi_matrix[,17],percdi_matrix[,16],percdi_matrix[,13],percdi_matrix[,2],percdi_matrix[,4],percdi_matrix[,7],percdi_matrix[,8],percdi_matrix[,14],percdi_matrix[,10],percdi_matrix[,15],percdi_matrix[,5],percdi_matrix[,6])
colnames(percdi_matrix_reord) <- c("NotI","SgrAI","BsrFI","NgoMIV","AgeI","MspI","SbfI","PstI","NsiI","ApoI","EcoRI","MluCI","MseI","NspI","NcoI","PciI","FatI","KpnI")

hv_di<-heatmap.2((percdi_matrix_reord), Rowv=NA, Colv=NA, scale="none", breaks=matrix_breaks, col= colorpanel(42,"#00FFFF","#000000","#FFE600"), main="Similarity Index  Dinucleotide", trace="none", denscol="red")


#Figure 6. Right: heatmap of the similarity indexes for the dinucleotide model. 
matrix_breaks<-c(-1.01,-1,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1)

perctri_matrix_reord<-cbind(perctri_matrix[,12],perctri_matrix[,18],perctri_matrix[,3],perctri_matrix[,11],perctri_matrix[,1],perctri_matrix[,9],perctri_matrix[,17],perctri_matrix[,16],perctri_matrix[,13],perctri_matrix[,2],perctri_matrix[,4],perctri_matrix[,7],perctri_matrix[,8],perctri_matrix[,14],perctri_matrix[,10],perctri_matrix[,15],perctri_matrix[,5],perctri_matrix[,6])
colnames(perctri_matrix_reord) <- c("NotI","SgrAI","BsrFI","NgoMIV","AgeI","MspI","SbfI","PstI","NsiI","ApoI","EcoRI","MluCI","MseI","NspI","NcoI","PciI","FatI","KpnI")

hv_tri<-heatmap.2((perctri_matrix_reord), Rowv=NA, Colv=NA, scale="none", breaks=matrix_breaks, col= colorpanel(42,"#00FFFF","#000000","#FFE600"), main="Heatmap of Bias Index Trinucleotides", trace="none", denscol="red")


#Figure 7. Center: Heatmap of the similarity indexes for the GC content model 
matrix_breaks<-c(-1.01,-1,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1)

percgc_matrix_reord<-cbind(percgc_matrix[,12],percgc_matrix[,18],percgc_matrix[,3],percgc_matrix[,11],percgc_matrix[,1],percgc_matrix[,9],percgc_matrix[,17],percgc_matrix[,16],percgc_matrix[,13],percgc_matrix[,2],percgc_matrix[,4],percgc_matrix[,7],percgc_matrix[,8],percgc_matrix[,14],percgc_matrix[,10],percgc_matrix[,15],percgc_matrix[,5],percgc_matrix[,6])
colnames(percgc_matrix_reord) <- c("NotI","SgrAI","BsrFI","NgoMIV","AgeI","MspI","SbfI","PstI","NsiI","ApoI","EcoRI","MluCI","MseI","NspI","NcoI","PciI","FatI","KpnI")

hv_gc<-heatmap.2((percgc_matrix_reord), Rowv=NA, Colv=NA,scale="none", breaks=matrix_breaks, col= colorpanel(42,"#00FFFF","#000000","#FFE600"), main="Similarity Index GC content", trace="none", denscol="red")


#Figure 7. Right: heatmap of the similarity indexes for the mononucleotide model. 
matrix_breaks<-c(-1.01,-1,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1)

percmono_matrix_reord<-cbind(percmono_matrix[,12],percmono_matrix[,18],percmono_matrix[,3],percmono_matrix[,11],percmono_matrix[,1],percmono_matrix[,9],percmono_matrix[,17],percmono_matrix[,16],percmono_matrix[,13],percmono_matrix[,2],percmono_matrix[,4],percmono_matrix[,7],percmono_matrix[,8],percmono_matrix[,14],percmono_matrix[,10],percmono_matrix[,15],percmono_matrix[,5],percmono_matrix[,6])
colnames(percmono_matrix_reord) <- c("NotI","SgrAI","BsrFI","NgoMIV","AgeI","MspI","SbfI","PstI","NsiI","ApoI","EcoRI","MluCI","MseI","NspI","NcoI","PciI","FatI","KpnI")

hv_mono<-heatmap.2((percmono_matrix_reord), Rowv=NA, Colv=NA, scale="none", breaks=matrix_breaks, col= colorpanel(42,"#00FFFF","#000000","#FFE600"), main="Similarity Index Mononucleotide", trace="none", denscol="red")


#Figure 8. Recovery of RAD-tags after in silico genome digestion and sequencing. Heatmap of the percentage of RAD-tags that produced more than one unique alignment to their reference genome. 
supp_perc_matrix_reord<-cbind(supp_perc_matrix[,12],supp_perc_matrix[,18],supp_perc_matrix[,3],supp_perc_matrix[,11],supp_perc_matrix[,1],supp_perc_matrix[,9],supp_perc_matrix[,17],supp_perc_matrix[,16],supp_perc_matrix[,13],supp_perc_matrix[,2],supp_perc_matrix[,4],supp_perc_matrix[,7],supp_perc_matrix[,8],supp_perc_matrix[,14],supp_perc_matrix[,10],supp_perc_matrix[,15],supp_perc_matrix[,5],supp_perc_matrix[,6])
colnames(supp_perc_matrix_reord)<-c("NotI","SgrAI","BsrFI","NgoMIV","AgeI","MspI","SbfI","PstI","NsiI","ApoI","EcoRI","MluCI","MseI","NspI","NcoI","PciI","FatI","KpnI")
hv<-heatmap.2((supp_perc_matrix_reord), Rowv=NA, Colv=NA, scale="none", col= colorpanel(17, "#000000","#FF0000"), main="Percentage of non-unique alignments", trace="none", denscol="green")

