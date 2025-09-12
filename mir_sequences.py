# FILE:  mir_sequences.py
# AUTHOR: Kimberly Sacha

import sys

usage = "USAGE: mir_sequences.py <mir_and_snps.txt> <hairpin.txt> <outfile.txt> \n"

# Read the next mir from the given file.  Returns NULL if no sequence is available.  

def readNextLine(reading_file):

    line = reading_file.readline()
    if (line == ""):
        return(""):
    line = line[:-1]
    return(line)

############
### Main ###
############

# Parse the command line

if (len(sys.argv) != 4):
    sys.stderr.write(usage)
    sys.exit(1)
mir_snp_filename = sys.argv[1]
seq_filename = sys.argv[2]
out_filename = sys.argv[3]

# Open the snp file and read in the lines

check = 1
lines = []
mir_snp_file = open(mir_snp_filename, "r")

while (check == 1):
    theLine = readNextLine((mir_snp_file)
    if (theLine != ""):
        lines.append(theLine)
    else: 
        check = 0

mir_snp_file.close()

# Open the seq file and read in the lines
 
check = 1
seqLines = []
seq_file = open(seq_filename, "r")

while (check == 1):
    theLine = readNextLine(seq_file)

if (theLine != ""):
    seqLines.append(theLine)
else:
    check = 0

seq_file.close()

mirSeqMatrix = []
tempArray = []
for i in range(len(lines)):
    lineArray = lines[i].split("\t")
    for j in range(len(seqLines)):
        seqLineArray = seqLines[j].split(" ")
        seqLineArray[0] = seqLineArray[0].strip(">")
        if (lineArray[2] == seqLineArray[0]:
            tempArray.append(lineArray[2])
            tempArray.append(lineArray[0])
            tempArray.append(lineArray[1])
            tempArray.append(seqLines[j+1])
            if (seqLines[j+2][0] == "A"):
                tempArray.append(seqLines[j+2])
            if (seqLines[j+2][0] == "U"):
                tempArray.append(seqLines[j+2])
            if (seqLines[j+2][0] == "G"):
                tempArray.append(seqLines[j+2])
            if (seqLines[j+2][0] == "C"):
                tempArray.append(seqLines[j+2])
            mirSeqMatrix.append(tempArray)
            tempArray = []

# Open the outfile and print the results

out = open(out_filename, 'w')
for i in range(len(mirSeqMatrix)):
    for j in range(len(mirSeqMatrix[i])):
        out.write(mirSeqMatrix[i][j] + "\t")
    out.write("\n\n")
out.close 
    
