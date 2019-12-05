from src.data_structures.set_abc import Set
import string


def parse_file(txt):
    genomeList = []
    currentString = ""
    while True:
        line = txt.readline()
        if len(line) == 0:
            # End of file
            break

        if line[0] == ">":  # append whenever we see the mark for a new strain's genome
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


def break_kmers(genome, ds, kmer_size):
    for i in range(len(genome) - kmer_size + 1):  # for each k-mer
        kmer = genome[i:i+kmer_size]
        ds.insert(kmer)
    return ds
