# FILE: mir_2_sequences.py
# AUTHOR: Kimberly Sacha

import sys

usage = "USAGE: mir_2_sequences.py <mir_and_snps.txt> <hairpin.txt> <mir_locations.txt> <alleles.txt> <snp_loc.txt> <outfile.txt> \n"


# Read the next mir from the given file.  Returns NULL if no sequence is available.
def readNextLine(reading_file):

  line = reading_file.readline()
  if (line == ""):
    return("")
  line = line[:-1]
  return(line)

# Return the complementary sequency string
def complement(s):
  basecomplement = {'A':'U', 'C':'G', 'G':'C','U':'A'}
  letters = list(s)
  letters = [basecomplement[base] for base in letters]
  return ''.join(letters)

 
############
### Main ###
############

# Parse the command line.
if (len(sys.argv) != 7):
  sys.stderr.write(usage)
  sys.exit(1)
mir_snp_filename = sys.argv[1]    #usually mir_names_and_rs_nums.txt
seq_filename = sys.argv[2]        #usually hairpin.txt
loc_filename = sys.argv[3]        #usually mir_locations.txt which is hsa.gff from miRBase
allele_filename = sys.argv[4]     #usualy allele_out.txt
snp_loc_filename = sys.argv[5]    #usually snomir_with_snps_070607.txt
out_filename = sys.argv[6]        #usually seqout_deluxe.txt


#open the snp file and read in the lines
check = 1
lines = []
mir_snp_file = open(mir_snp_filename, "r")

while (check == 1):
  theLine = readNextLine(mir_snp_file)
  if (theLine != ""):
      lines.append(theLine)
  else:
    check = 0

mir_snp_file.close()

#open the seq file and read in the lines
check = 1
seqLines = []                                       #seqLines is an array of hairpin.txt lines
seq_file = open(seq_filename, "r")

while (check == 1):
  theLine = readNextLine(seq_file)
  if (theLine != ""):
      seqLines.append(theLine)
  else:
    check = 0

seq_file.close()

#open the location file and read in the lines
check = 1
locLines = []                                       #locLines is an array of mir_locations.txt lines
loc_file = open(loc_filename, "r")

while (check == 1):
  locLine = readNextLine(loc_file)
  if (locLine != ""):
      locLines.append(locLine)
  else:
    check = 0

loc_file.close()


#open the allele file and read in the lines
check = 1
alleleLines = []                                     #alleleLines is an array of allele.txt lines
allele_file = open(allele_filename, "r")

while (check == 1):
  alleleLine = readNextLine(allele_file)
  if (alleleLine != ""):
      alleleLines.append(alleleLine)
  else:
    check = 0

allele_file.close()


#open the SNP location file and read in the lines
check = 1
snpLocLines = []                                     #snpLocLines is an array of snp_loc.txt lines
snp_loc_file = open(snp_loc_filename, "r")

while (check == 1):
  snpLocLine = readNextLine(snp_loc_file)
  if (snpLocLine != ""):
      snpLocLines.append(snpLocLine)
  else:
    check = 0

snp_loc_file.close()


#############################
####  mirSeqMatrix  #########
#############################

#use seqLines make a matrix called mirSeqMatrix with name, rs num, chrom, dna
mirSeqMatrix = []
tempArray = []
for i in range(len(lines)):
    lineArray = lines[i].split("\t")
    for j in range(len(seqLines)):
        seqLineArray = seqLines[j].split(" ")         #seqLineArray is an array of hairpin.txt information
        seqLineArray[0] = seqLineArray[0].strip(">")
        if (lineArray[2] == seqLineArray[0]):         #if miRNA names match
            tempArray.append(lineArray[2])            #add mimRNA name
            tempArray.append(lineArray[0])            #add rs number
            tempArray.append(lineArray[1])            #add chrom
            tempArray.append(seqLines[j+1])           #add a line of dna
            if (seqLines[j+2][0] == "A"): 
                tempArray.append(seqLines[j+2])       #add a line of dna starting with 'A'
            if (seqLines[j+2][0] == "U"): 
                tempArray.append(seqLines[j+2])
            if (seqLines[j+2][0] == "G"): 
                tempArray.append(seqLines[j+2])
            if (seqLines[j+2][0] == "C"): 
                tempArray.append(seqLines[j+2])
            mirSeqMatrix.append(tempArray)
            tempArray = []

for i in range(len(mirSeqMatrix)):
  mirSeqMatrix[i][3] = mirSeqMatrix[i][3] + mirSeqMatrix[i][4]
  mirSeqMatrix[i][4] = ""


#############################
####  alleleMatrix  #########
#############################
            
#use mirSeqMatrix and locations to make a matrix called alleleMatrix with name, rs num, allele, loc start, loc stop
name = ""
alleleLinesArray = []
alleleLinesArray1 = []
alleleMatrix = []
tempArray2 = []
locStart = ""
locEnd = ""

for i in range(len(mirSeqMatrix)):
  for j in range(len(locLines)):
    locLineArray = locLines[j].split("\t")
    if locLineArray[0][0] != '#':
      name = locLineArray[8][21:]
      name = name.strip('";')
      locStart = locLineArray[3]
      locEnd = locLineArray[4]
      read = locLineArray[6]
    if mirSeqMatrix[i][0] == name:
      for k in range(len(alleleLines)):
        alleleLinesArray = alleleLines[k].split("\t")
        alleleLinesArray1 = alleleLinesArray[0].split(" ")
        if alleleLinesArray1[0] == mirSeqMatrix[i][1]:      #if ... match
          tempArray2.append(name)                           #add name
          tempArray2.append(mirSeqMatrix[i][1])             #add rs num
          tempArray2.append(alleleLinesArray[1][9:12])      #add alleles
          tempArray2.append(locStart)                       #add start location
          tempArray2.append(locEnd)                         #add end location
          tempArray2.append(read)                           #add reading direction
          alleleMatrix.append(tempArray2)
          tempArray2 = []


#############################
####  snpLocMatrix  #########
#############################

snpLocMatrix = []
tempArray3 = []
locChromArray = []
thelocArray = []

#use snpLocLines to make a matrix called snpLocMatrix with SNP rs num, chrom, locations
for i in range(len(snpLocLines)):
  snpLocLinesArray = snpLocLines[i].split(" ")
  for j in range(len(mirSeqMatrix)):
    if snpLocLinesArray[0] == mirSeqMatrix[j][1]:       #if the rs nums match
      locChromArray = snpLocLinesArray[2].split(":")
      theLocArray = locChromArray[1].split("-") 
      
      #theLoc = theLoc.strip("
      tempArray3.append(snpLocLinesArray[0])            #add rs num
      tempArray3.append(locChromArray[0])               #add chrom
      tempArray3.append(theLocArray[0])                 #add location
      snpLocMatrix.append(tempArray3)
      tempArray3 = []

################################
#####outputting DNA#############

snp = 0
start = 0
place = 0
dna = ""
dnaReverse = ""
newdna1 = ""
newdna2 = ""
newReverse1 = ""
newReverse2 = ""
alleles = ""
newallele = ""

for i in range(len(snpLocMatrix)):
  for j in range(len(alleleMatrix)):
    if (snpLocMatrix[i][0] == alleleMatrix[j][1]):
      for k in range(len(mirSeqMatrix)):
        if (alleleMatrix[j][1] == mirSeqMatrix[k][1]):
          alleles = alleleMatrix[j][2]
          alleles = alleles.replace("T", "U")
          snp = int(snpLocMatrix[i][2])
          start = int(alleleMatrix[j][3])
          place = snp - start
          dna = mirSeqMatrix[k][3]
          newdna1 = dna[:place] + " " + alleles[0] + " " + dna[place+1:]
          newdna2 = dna[:place] + " " + alleles[2] + " " + dna[place+1:]
          dnaReverse = dna[::-1]
          dnaReverse = complement(dnaReverse)   
          newReverse1 = dnaReverse[:place] + " " + alleles[0] + " " + dnaReverse[place+1:]
          newReverse2 = dnaReverse[:place] + " " + alleles[2] + " " + dnaReverse[place+1:]
          print "\n" + mirSeqMatrix[k][0] + "\t" + alleleMatrix[j][5] 
          print "Forward - replaced " + dna[place] + " with " + alleles[0] + "\t" + newdna1
          print "Forward - replaced " + dna[place] + " with " + alleles[2] + "\t" + newdna2
          print "Reverse - replaced " + dnaReverse[place] + " with " + alleles[0] + "\t" + newReverse1
          print "Reverse - replaced " + dnaReverse[place] + " with " + alleles[2] + "\t" + newReverse2
        

#open the outfile and print the results
out = open(out_filename, 'w')
for i  in range(len(mirSeqMatrix)):
  for j in range(len(mirSeqMatrix[i])):
    out.write(mirSeqMatrix[i][j] + "\t")
  out.write("\n\n")

out.close()
