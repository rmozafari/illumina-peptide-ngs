# Prints the contents of a file on screen line by line. Useful for seeing what's inside large files.




#Enter the name of file and how many lines you want to read below
nameof_file="Sample_1_NoIndex_L001_R1_001.fastq"
numberof_linestoread=100

#-------------------------------------------
i=0
for line in open(nameof_file):
    print (line)
    i=i+1
    if i==numberof_linestoread:
        break
    


    
