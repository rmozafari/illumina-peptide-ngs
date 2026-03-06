nameof_inputfile="20170111_All_DNA.txt"
nameof_outputfile="20170111_All_DNA_Protease.txt"
forward_primer="GCCGGCGGAGGT"
#reverse_primer="GGCTCCAGCCAA"

#------------------------------------------------------------
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna #These two are necessary for operations like translation and reverse complement. 

inputfile=open(nameof_inputfile)
outputfile=open(nameof_outputfile, "w")
library_composition = {}

total_seq_num=0 #just to know how many sequences existed
translatable_seq_num=0 #and how many of them are translatable

for line in inputfile: 
    seq_is_correct=0 #I restart this variable everytime in the loop, since later, i will use it to figure out if the 5' primer is in the sequence or not
    total_seq_num +=1
    
    if line[len(line)-1:len(line)]=='\n':#for some unknown reason, sometimes line has a \n on it, which needs to be removed
        line = line[0:len(line)-1]

    #read_dnaseq=Seq(line, generic_dna) #necessary so that sequence can be translated by Biopython (now sequence is not just a string, but a sequence class that is recognized by biopython)
    #read_dnaseq_revcomp=read_dnaseq.reverse_complement() #takes reverse complement of readstart_read=line.find("(")+1+startingresidue
    #read_dnaseq_revcomp_str=str(read_dnaseq_revcomp)#I make it into a string so I can search it
    #primer_read=read_dnaseq_revcomp_str.find(forward_primer) #basically I look upstream of ATG, and find things that have the right primers

    #if primer_read>-1:#If the primers are in the right place
    #    new_frame=read_dnaseq_revcomp_str[primer_read:len(read_dnaseq_revcomp_str)]
    #    seq_is_correct=1
    #else:
    primer_read=line.find(forward_primer)
    if primer_read>-1:
        #new_frame=line[primer_read:len(line)]
        seq_is_correct=1

            
    if seq_is_correct==1:
        #start_read=new_frame.find("ATG")
        #end_read=new_frame.find(reverse_primer)
        #if end_read==-1:
        #    end_read=start_read-9

        #read=new_frame[start_read:end_read+9]
        #read_dnaseq=Seq(read, generic_dna)
        #read_aaseq=read_dnaseq.translate() 
        outputfile.write(str(line))
        outputfile.write("\n")
        
##        if end_read in library_composition.keys():
##            library_composition[end_read]=library_composition[end_read]+1
##        else:
##            library_composition[end_read]=1    

        translatable_seq_num +=1


outputfile.close()
