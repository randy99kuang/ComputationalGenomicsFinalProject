from src.data_structures.set_abc import Set
import sys


def parse_file(txt, ds):
    count = -1
    genomeList = []
    while True:
        line = txt.readline()
        if len(line) == 0:
            # End of file
            break
        if line[0] == ">":
            count += 1
            genomeList.append("")
            continue
        # line = line.replace("-", "")
        genomeList[count] = genomeList[count] + line
    
def break_kmers(genome, ds, kmer_size):
    for i in range(len(genome) - kmer_size + 1):  # for each k-mer
            kmer = genome[i:i+kmer_size]
            ds.insert(kmer)
    return ds