import sys
import math
import os

#input arguments:
#1. name of input file
#2. the size of the window wanted
#3. the prefix to append to the output files
input_file = sys.argv[1]
windowsize = int(sys.argv[2])
missingthresh = float(sys.argv[3]) #any window with an individual with this proportion of missing data or higher will be removed
prefix = sys.argv[4]



#generate info about alignment and number of cycles needed
with open(input_file, 'r') as fasta_input:
    infile = fasta_input.read().splitlines()
    alignmentlength = len(infile[1])
    cycles = int(math.ceil(alignmentlength/windowsize))


#for loop, which cycles through each line, checks if the line is a header, and splits the read into the desired length
#also a missing data threshold, where 'N' is missing data
#fasta files should not be 50bp breaks, and there should not be spaces between new scaffolds or samples
for i in range(1,cycles+1):
    with open(input_file, 'r') as fasta_input:
        windowstart = -windowsize+windowsize*i
        windowend = windowsize*i
        has_lines = True
        name_line = ''
        seq_line = ''
        if windowend < alignmentlength:
            filename = prefix+"_bp_"+str(windowstart+1)+"-"+str(windowend)+".fa"
        else:
            filename = prefix+"_bp_"+str(windowstart+1)+"-"+str(alignmentlength)+".fa"
        with open(filename, 'w') as fasta_output:
            deletefile = False
            while has_lines:
                fasta_line = fasta_input.readline()
                if len(fasta_line) == 0:
                    has_lines = False
                else:
                    fasta_line = fasta_line.strip()
                    if fasta_line.startswith('>'):
                        name_line = fasta_line
                    elif len(fasta_line) > 0:
                        seq_line = fasta_line
                        if name_line.startswith('>REF '):#remove reference sequence
                            missingcount = 0
                        else:
                            fasta_output.write(name_line+'\n')
                            fasta_output.write(seq_line[windowstart:windowend]+'\n')
                            missingcount = seq_line[windowstart:windowend].count('N')
                        if missingcount/windowsize >= missingthresh:
                            deletefile = True
        if deletefile == True:
            os.remove(filename)

