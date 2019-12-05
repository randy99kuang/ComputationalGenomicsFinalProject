import sys
from src.analysis.parse import *
from src.data_structures.hashset import HashSet
from src.data_structures.bloom_filter import BloomFilter


def getDataStructure(ds):
    if ds == "HashSet":
        return HashSet()
    else:
        return BloomFilter()


def readHIV(kmer_size, ds):
    g1 = open("..\\data\\HIV\\hiv16.fasta", "r")
    list1 = parse_file(g1)

    g2 = open("..\\data\\HIV\\hiv32.fasta", "r")
    list2 = parse_file(g2)

    genomeList = list1 + list2

    ds_list = []
    for i in range(len(genomeList)):
        ds_list.append(getDataStructure(ds))

    for i in range(len(genomeList)):
        ds_list[i] = break_kmers(genomeList[i], ds_list[i], kmer_size)

    return ds_list

def merge(numIntersections, ds_list):
    for i in range(numIntersections):   # start by intersecting all sets numIntersection times

    return ds_list[0]   # return the only element left in the list