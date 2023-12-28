import sys
import math
import os
import argparse


##Define variables
######
parser = argparse.ArgumentParser()

parser.add_argument('--input_file', '-i', help="Input file, in fasta format. Fasta must be single-line fasta format. Mandatory parameter.")

parser.add_argument('--outdir', '-o', help="Output directory in which to write the output files. Mandatory parameter.")

parser.add_argument('--missingthresh', '-m', nargs='?', default=1, type=float, help="The missing data threshold, per individual. Argument should be between 0-1. Includes a sample in the fasta if the percentage of missing data ('N's) is less than or equal to this value. 0 = no missing data allowed, 1 = all missing data allowed. Default = 1. Optional parameter.")

args = parser.parse_args()
######


#Main code
#Some of the name manipulation could be altered if desired, but should function fine as-is
with open(args.input_file, 'r') as fasta_input:
    has_lines = True
    name_line = ''
    seq_line = ''
    while has_lines:
        fasta_line = fasta_input.readline()
        if len(fasta_line) == 0:
            has_lines = False
        else:
            fasta_line = fasta_line.strip()
            if fasta_line.startswith('>'):
                name_line = fasta_line
                newname_line = name_line.replace(">", ">"+os.path.basename(args.input_file)+ " ")
                outfile_name = os.path.join(args.outdir, name_line.replace(">", "").replace(":", "_")+".fa")
            elif len(fasta_line) > 0:
                seq_line = fasta_line
                missingcount = seq_line.count('N')
                if missingcount/len(seq_line) <= args.missingthresh:
                    with open(outfile_name, 'a') as fasta_output:
                        fasta_output.write(newname_line+'\n')
                        fasta_output.write(seq_line+'\n')

