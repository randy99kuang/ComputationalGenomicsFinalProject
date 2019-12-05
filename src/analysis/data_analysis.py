from src.analysis.cluster import *
from pympler import asizeof

hs_final = None
bf_final = None
bf_test_list = []


def preprocess():
    hs_list = readHIV(100, "HashSet")
    hs_final = merge(0, hs_list)
    print("Final size of all union hash set:", hs_final.getSize())
    print("true size in bytes of hash set:", asizeof.asizeof(hs_final))

    bf_list = readHIV(100, "BloomFilter")
    bf_final = merge(3, bf_list)
    print("bit size of bloom filter:", bf_final.getBitSize())
    print("true size in bytes of bloom filter:", asizeof.asizeof(bf_final))

    f = open("..\\data\\HIV\\hiv5Test.fasta", "r")
    genome_test_list = parse_file(f)

    for i in range(len(genome_test_list)):
        bf_test_list.append(getDataStructure("BloomFilter"))

    for i in range(len(genome_test_list)):
        bf_test_list[i] = break_kmers(genome_test_list[i], bf_test_list[i], 100)
        print(bf_test_list[i].getBitSize())
