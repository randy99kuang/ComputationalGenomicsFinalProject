import sys
sys.path.insert(1, '../../')
from src.data_structures.set_abc import Set
import string

#Parse an HIV file into the individual sequences
def parse_file(txt):
    genomeList = []
    currentString = ""
    #keep going until you reach the end of the file
    while True:
        #read each line
        line = txt.readline()
        if len(line) == 0:
            # End of file
            break

        if line[0] == ">":  # append whenever we see the mark for a new strain's genome
            #get list of genomes
            if len(currentString) > 0:
                currentString = currentString.translate({ord(c): None for c in string.whitespace})  # strip whitespaces
                genomeList.append(currentString)
                currentString = ""
            continue
        else:
            line = line.replace("-", "")    # remove any dashes in the line
            currentString = currentString + line

    # append list string to list
    currentString = currentString.translate({ord(c): None for c in string.whitespace})  # strip whitespaces
    genomeList.append(currentString)

    return genomeList

#same as parse_file, except for E. coli
def parse_file_ecoli(txt):
    currentString = ""
    txt.readline()
    genomeList = []
    while True:
        line = txt.readline()
        if len(line) == 0:
            # End of file
            break
        if line[0] == ">":  # append whenever we see the mark for a new strain's genome
            continue
        else:
            line = line.replace("-", "")  # remove any dashes in the line
            currentString = currentString + line

    # append list string to list
    currentString = currentString.translate({ord(c): None for c in string.whitespace})  # strip whitespaces
    genomeList.append(currentString)
    return genomeList

#Break a genome into kmer_size kmers and insert into a data structure given
def break_kmers(genome, ds, kmer_size):
    for i in range(len(genome) - kmer_size + 1):  # for each k-mer
        kmer = genome[i:i+kmer_size]
        ds.insert(kmer)
    return ds
