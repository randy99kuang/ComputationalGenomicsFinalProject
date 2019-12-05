import sys
sys.path.insert(1, '../../')

from src.analysis.parse import *
from src.data_structures.hashset import HashSet
from src.data_structures.bloom_filter import BloomFilter

HIV_GENOME_AVERAGE_SIZE = 9700  # approximate average length of HIV genomes used for testing
FPR = 0.005                     # the false positive rate we will use


def getDataStructure(ds):
    if ds == "HashSet":
        return HashSet()
    else:
        return BloomFilter(HIV_GENOME_AVERAGE_SIZE, FPR)


def readHIV(kmer_size, ds):
    g1 = open("hiv16.fasta", "r")
    list1 = parse_file(g1)

    g2 = open("hiv32.fasta", "r")
    list2 = parse_file(g2)

    genomeList = list1 + list2

    ds_list = []
    for i in range(len(genomeList)):
        ds_list.append(getDataStructure(ds))

    for i in range(len(genomeList)):
        # print(len(genomeList[i]))
        ds_list[i] = break_kmers(genomeList[i], ds_list[i], kmer_size)
        # print(sum(ds_list[i].getDictionary().values()))

    return ds_list


def merge(numIntersections, ds_list):
    for i in range(numIntersections):   # start by intersecting all sets numIntersection times
        new_list = []
        currentSize = int(len(ds_list) / 2)
        for j in range(currentSize):
            a = ds_list[0]
            b = ds_list[1]
            new_list.append(a.intersection(b))
            del ds_list[:2]             # remove the first two elements
        ds_list = new_list

    while len(ds_list) > 1:
        new_list = []
        currentSize = int(len(ds_list) / 2)
        for j in range(currentSize):
            a = ds_list[0]
            b = ds_list[1]
            new_list.append(a.union(b))
            del ds_list[:2]             # remove the first two elements
        ds_list = new_list

    return ds_list[0]                   # return the only element left in the list
