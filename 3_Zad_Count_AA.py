#Enter input file name below.

nameof_inputfile="TLN_Tryp.txt"
nameof_outputfile="TLN_Tryp_Count_AA.txt"
#-----------------------------------------------------------------------

outputfile=open(nameof_outputfile, "w")
aminoacids_str = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','*','X']

total_Sequence_Count=0
Sequence_Coll={} #Name of the dictionary containing everything

#we need to first make an entry for every possibility, before filling them up
for i in range(0,len(aminoacids_str)):
    Sequence_Coll[aminoacids_str[i]]=0

#Now we fill them up
for line in open(nameof_inputfile):
    #since we need to go one amino acid at a time:
    for i in range(0,len(line.strip())):
        Sequence_Coll[line.strip()[i]]+=1 

    #this is just to make sure we know about the things that are going on.
    total_Sequence_Count +=1
    if total_Sequence_Count%1000000==0:
        print(total_Sequence_Count)
    
#I have to sort it before I can write it line by line
sequencecounts=sorted(Sequence_Coll.items())
for i in range (0, len(sequencecounts)):
    outputfile.write(str(sequencecounts[i]))
    outputfile.write("\n")
    
outputfile.close()


