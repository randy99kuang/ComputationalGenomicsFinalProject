import sys
from parse import *
from src.data_structures.hashset import HashSet

kmer_size = 20

g = open(sys.argv[1], "r")
genomeList = parse_file(g)

ds_list = []
for i in range(len(genomeList)):
    ds_list.append(HashSet)

for i in range(len(genomeList)):
    ds_list[i] = break_kmers(genomeList[i], ds_list[i], kmer_size)