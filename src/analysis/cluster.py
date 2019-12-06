from src.analysis.parse import *
from src.data_structures.hashset import HashSet
from src.data_structures.bloom_filter import BloomFilter
from src.data_structures.counting_filter import CountingFilter
from pathlib import Path

import sys
sys.path.insert(1, '../../')

HIV_GENOME_AVERAGE_SIZE = 9700  # approximate average length of HIV genomes used for testing
FPR = 0.005                     # the false positive rate we will use
ECOLI_GENOME_AVERAGE_SIZE = 4700000

def getDataStructure(ds):
    if ds == "HashSet":
        return HashSet()
    elif ds == "BloomFilter":
        return BloomFilter(HIV_GENOME_AVERAGE_SIZE, FPR)
    else:
        return CountingFilter(HIV_GENOME_AVERAGE_SIZE, FPR)

def getDataStructureEcoli(ds):
    if ds == "HashSet":
        return HashSet()
    elif ds == "BloomFilter":
        return BloomFilter(ECOLI_GENOME_AVERAGE_SIZE, FPR)
    else:
        return CountingFilter(ECOLI_GENOME_AVERAGE_SIZE, FPR)


def readHIV(kmer_size, ds):
    hiv16 = Path("../data/HIV/hiv16.fasta")
    g1 = open(hiv16, "r")
    list1 = parse_file(g1)

    hiv32 = Path("../data/HIV/hiv32.fasta")
    g2 = open(hiv32, "r")
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


def readECOLI(kmer_size, ds):
    ecoli1 = Path("../data/ECOLI/ecoli1.fas")
    g1 = open(ecoli1, "r")
    list1 = parse_file_ecoli(g1)

    ecoli2 = Path("../data/ECOLI/ecoli2.fas")
    g2 = open(ecoli1, "r")
    list2 = parse_file_ecoli(g2)

    ecoli3 = Path("../data/ECOLI/ecoli3.fas")
    g3 = open(ecoli1, "r")
    list3 = parse_file_ecoli(g3)

    ecoli4 = Path("../data/ECOLI/ecoli4.fas")
    g4 = open(ecoli1, "r")
    list4 = parse_file_ecoli(g4)

    genomeList = list1 + list2 + list3 + list4

    ds_list = []
    for i in range(len(genomeList)):
        ds_list.append(getDataStructureEcoli(ds))

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
