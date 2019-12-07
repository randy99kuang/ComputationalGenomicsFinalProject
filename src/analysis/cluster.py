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

#Return a data structure with values that fit our specific data sets
def getDataStructure(ds):
    if ds == "HashSet":
        return HashSet()
    elif ds == "BloomFilter":
        return BloomFilter(HIV_GENOME_AVERAGE_SIZE, FPR)
    else:
        return CountingFilter(HIV_GENOME_AVERAGE_SIZE, FPR)

#Return a data structure with values that fit our specific data sets
def getDataStructureEcoli(ds):
    if ds == "HashSet":
        return HashSet()
    elif ds == "BloomFilter":
        return BloomFilter(ECOLI_GENOME_AVERAGE_SIZE, FPR)
    else:
        return CountingFilter(ECOLI_GENOME_AVERAGE_SIZE, FPR)

#Read the HIV file, parse data, and put into data structure with correct kmer_size
def readHIV(kmer_size, ds):
    hiv16 = Path("../data/HIV/hiv16.fasta")
    g1 = open(hiv16, "r")
    list1 = parse_file(g1)

    hiv32 = Path("../data/HIV/hiv32.fasta")
    g2 = open(hiv32, "r")
    list2 = parse_file(g2)

    genomeList = list1 + list2

    #Put all strains into a list of data structures,  one DS for each strain
    ds_list = []
    for i in range(len(genomeList)):
        ds_list.append(getDataStructure(ds))

    #put kmers into the data structures
    for i in range(len(genomeList)):
        # print(len(genomeList[i]))
        ds_list[i] = break_kmers(genomeList[i], ds_list[i], kmer_size)
        # print(sum(ds_list[i].getDictionary().values()))

    return ds_list

#Same thing as readHIV above, except for E. coli
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

#Merge numIntersections amount of times between the filters in ds_list.
#This runs numIntersections on pairs in ds_list, halving the number of data
#structures outputted as inputted for each numIntersection, and unions the rest
#until there is one output data structure
def merge(numIntersections, ds_list):
    for i in range(numIntersections):   # start by intersecting all sets numIntersection times
        new_list = []
        currentSize = int(len(ds_list) / 2)
        #Intersect half of  them with one another
        for j in range(currentSize):
            a = ds_list[0]
            b = ds_list[1]
            new_list.append(a.intersection(b))
            del ds_list[:2]             # remove the first two elements
        ds_list = new_list

    #Keep unioning until the list is length one length, which is the final data
    #structure that contains all the values from all the nonintersected sets
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
