import sys
sys.path.insert(1, '../../')
from src.analysis.cluster import *
from pympler import asizeof

hs_final = None
bf_final = None
hs_test_list = []
bf_test_list = []


def preprocessHIV(kmer_length):
    """
    Method that reads the 32 HIV strains into both a hashset and a bloom filter, and then reads the 5 test strains
    into both as well. This method will be called before many of the data-analysis methods
    """

    global hs_final, bf_final, hs_test_list, bf_test_list   # mark these as global variables so we can edit them

    hs_list = readHIV(kmer_length, "HashSet")
    hs_final = merge(0, hs_list)
    print("Final length of all union hash set:", hs_final.getSize())
    print("true size in bytes of hash set:", asizeof.asizeof(hs_final))

    bf_list = readHIV(kmer_length, "BloomFilter")
    bf_final = merge(0, bf_list)
    print("bit size of bloom filter:", bf_final.getBitSize())
    print("true size in bytes of bloom filter:", asizeof.asizeof(bf_final))

    f = open("hiv5Test.fasta", "r")
    genome_test_list = parse_file(f)

    # Fill in bf_test_list
    for i in range(len(genome_test_list)):
        bf_test_list.append(getDataStructure("BloomFilter"))

    for i in range(len(genome_test_list)):
        bf_test_list[i] = break_kmers(genome_test_list[i], bf_test_list[i], kmer_length)

    # Fill in hs_test_list
    for i in range(len(genome_test_list)):
        hs_test_list.append(getDataStructure("HashSet"))

    for i in range(len(genome_test_list)):
        hs_test_list[i] = break_kmers(genome_test_list[i], hs_test_list[i], kmer_length)


def kmerLength_vs_size():
    return 0


def spaceAnalysis():
    """
    This will use Pympler's asizeof.asizeof() method to find the deep size of bloom filters and hash sets
    """

    preprocessHIV(100)
    print("Final length of all union hash set:", hs_final.getSize())
    print("true size in bytes of hash set:", asizeof.asizeof(hs_final))
    print("bit size of bloom filter:", bf_final.getBitSize())
    print("true size in bytes of bloom filter:", asizeof.asizeof(bf_final))


def timeAnalysis():
    """
    This will use Python's timeit method to find how long it takes to build the bloom filters and hash sets
    """

    return 0
