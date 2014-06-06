#!/usr/bin/env python
# Usage: python sequence_probability_paper.py distributionfile patternsfile

#-------------------------------FUNCTIONS-LIST----------------------------------------
# writelog(filename,text)
# scalar_multiply(matrix, scalar,sizerow,sizecol,r)
# scalar_multiply_vec(vector,scalar,sizecol,r)
# obtain_model(fields,nt,dint,freqdint,trint)
# sequence2numbers(sequence)
# prob_nt(sequence,nt)
# prob_dint(sequence,dint,nt)
# prob_trint(sequence,trint,freqdint,nt)
# obtain_paths(pattern)
# find_paths(sequences,path,pospath)
# int_main(distributionfile, patternsfile)

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

def scalar_multiply(matrix, scalar,sizerow,sizecol,r):          
        for row in range(sizerow):
            for column in range(sizecol):
                r[row][column]=matrix[row][column]*scalar

def scalar_multiply_vec(vector,scalar,sizecol,r):           
            for column in range(sizecol):
                r[column]=vector[column]*scalar

def obtain_model(fields,nt,dint,freqdint,trint):
    totalnt=float(fields[1])
    totalacgt=float(fields[2])
    #nucleotides frequencies
    countnt=0
    #AGCT
    for i in range(4):      
        nt[0+i]=float(fields[3+i])/totalacgt
    #AGCT
    #dinucleotides frequencies
    totaldint=0
    for i in range(4):
        for j in range(4):
            dint[i][0+j]=float(fields[8+j+i*4])
            totaldint=totaldint+dint[i][0+j]
            freqdint[i][0+j]=dint[i][0+j]
            
    #obtain frequencies dividing counters by totaldint
    scalar_multiply(dint,float(1/totaldint),4,4,dint)    
    scalar_multiply(freqdint,float(1/totaldint),4,4,freqdint)    
    #obtain conditional probabilities
    #divide frequencies by the sum of each row of frecuencies
    for i in range(4):
        scalar_multiply_vec(dint[i],float(1/sum(dint[i])),4,dint[i])
    #AGCT
    #trinucleotides frequencies     
    totaltrint=0
    for i in range(16):     
        for j in range(4):
            trint[i][0+j]=float(fields[24+j+i*4])
            totaltrint=totaltrint+trint[i][0+j]
    scalar_multiply(trint,float(1/totaltrint),16,4,trint)           
    for i in range(16):
        scalar_multiply_vec(trint[i],float(1/sum(trint[i])),4,trint[i])

def sequence2numbers(sequence):
    #Maps sequence2numbers
    #A->0, G->1, C->2, T->3
    options={'A':0,'G':1,'C':2,'T':3}
    size=len(sequence)
    result=[-1 for x in xrange(size)]
    for i in range(size):         
        result[i]=options[sequence[i]]        
    return result

def prob_nt(sequence,nt):
        numbers=sequence2numbers(sequence)
        prob=1
        for j in range(len(numbers)):
                prob=prob*nt[numbers[j]]
        return prob        
    
def prob_dint(sequence,dint,nt):
        numbers=sequence2numbers(sequence)
        prob=nt[numbers[0]]
        for j in range(len(numbers)-1):
                prob=prob*dint[numbers[j]][numbers[j+1]]
        return prob        
                
def prob_trint(sequence,trint,freqdint,nt):
        numbers=sequence2numbers(sequence)
        prob=freqdint[numbers[0]][numbers[1]]
        for j in range(len(numbers)-2):
                row=numbers[j]*4+numbers[j+1]
                prob=prob*trint[row][numbers[j+2]]
        return prob 
        
            
def obtain_paths(pattern):    
        #obtain all possible paths for a motif
        #Replacing IUPAC dictionary
        iupacdict = {'M':'[AC]','R':'[AG]','W':'[AT]','S':'[CG]','Y':'[CT]','K':'[GT]','V':'[ACG]','H':'[ACT]','D':'[AGT]','B':'[CGT]','X':'[ACGT]','N':'[ACGT]'}    
        for key, value in iupacdict.iteritems():
                pattern=pattern.replace(key,value)

        paths=[]
        npaths=[]
        pospaths=[]
        tempstring=[]
        len1=len(pattern)
        flag=0
        match=''
        count=0
        #Obtaining template(match), paths and positions
        for i in range(len1):
            if flag==1 :
                if pattern[i]==']':
                    paths.append(tempstring)
                    npaths.append(len(tempstring))
                    tempstring=[]
                    flag=0                
                else:
                    tempstring.append(pattern[i])
            elif pattern[i]=='[':
                 count=count+1
                 flag=1
                 match=match+'-'
                 pospaths.append(count-1)
            else:
                match=match+pattern[i]
                count=count+1
        #Finding paths and saving everything in the variable sequences
        sequences=[]
        sequences.append(match)    
        for j in range(len(paths)):
            sequences=find_paths(sequences,paths[j],pospaths[j])
        return sequences     

                    
                
            
def find_paths(sequences,path,pospath):
    templates=sequences
    sequences=[]
    for j in range(len(templates)):
        temp=templates[j]
        for i in range(len(path)):
            seq=temp[:pospath]+path[i]+temp[pospath+1:]
            sequences.append(seq)
    return  sequences      

    

def int_main(distributionfile, patternsfile):
        import os
        from datetime import date
        now = date.today()  
        logfile='sequence_probability_'+str(now)
        writelog(logfile,"--------------------")
        writelog(logfile,"Input files: distributions -"+distributionfile+', patterns '+patternsfile) 
        d = open(distributionfile,'r')
        p=open(patternsfile,'r')
        text='\t'
        for line in p:
                 pattern,patternname=line.split()
                 text=text+patternname+'\t'        
        p.close()
        lines=d.read().split('\n')
        #fOR GENOMES
        filent = open(distributionfile+'_nt','a')
        filedint = open(distributionfile+'_dint','a')
        filetrint = open(distributionfile+'_trint','a')
        filent.write(text+'\n')
        filedint.write(text+'\n')
        filetrint.write(text+'\n')
	# remove first line (header) and last line(empty)
        for line in lines[1:-1] :                
                #mononucleotide frequencies
                nt = [0 for x in xrange(4)]
                #dinucleotide transition prob 
                dint = [[0 for x in xrange(4)] for x in xrange(4)]
                #dinucleotide frequencies
                fdint = [[0 for x in xrange(4)] for x in xrange(4)]
                #trinucleotide transition prob.
                trint = [[0 for x in xrange(4)] for x in xrange(16)]
                fields = line.split("\t")                
                scode=fields[0]                
                textnt=scode+'\t'
                textdint=scode+'\t'
                texttrint=scode+'\t'
                #obtain model (frequencies, conditional probabilities)
                if(len(fields[1])>0):
                        obtain_model(fields,nt,dint,fdint,trint)
                p=open(patternsfile,'r')
                for line in p:
                            pattern,patternname=line.split()
                            #obtain all sequences from a pattern
                            sequences=obtain_paths(pattern)
                            probnt=0
                            probdint=0
                            probtrint=0
                            #adding probabilities from all possible paths in a pattern
                            for nseq in range(len(sequences)):
                                    probnt=probnt+prob_nt(sequences[nseq],nt)
                                    probdint=probdint+ prob_dint(sequences[nseq],dint,nt)
                                    probtrint=probtrint+prob_trint(sequences[nseq],trint,fdint,nt)
                            textnt=textnt+str(probnt)+'\t'
                            textdint=textdint+str(probdint)+'\t'
                            texttrint=texttrint+str(probtrint)+'\t'
               
                textnt=textnt+'\n'
                textdint=textdint+'\n'
                texttrint=texttrint+'\n'
                filent.write(textnt)
                filedint.write(textdint)
                filetrint.write(texttrint)
                p.close()
        filent.close()
        filedint.close()
        filetrint.close()
               
    
#----------------------------------------------------------------------------------        

import sys
#Pass distribution file and patternsfile from shell
int_main(sys.argv[1],sys.argv[2])


        

