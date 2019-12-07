import matplotlib.pyplot as plt
import timeit
from src.analysis.cluster import *
from pympler import asizeof

sys.path.insert(1, '../../')
preprocess_bloom_filter = False
hs_final = None
bf_final = None
cf_final = None
hs_test_list = []

"""
These functions within this file share the exact same architecture and logic
as the ones within data_analysis_hiv.py. This means that the only things that
have changed within this file compared to that are that the name HIV in
variables and functions are changed to ecoli (or ecoli is added). This file also
references data in locations where we expect to find e. coli data.
"""

def preprocessAllECOLI(kmer_length, numIntersections):
    """
    Method that reads the 32 ECOLI strains into both a hashset and a bloom filter, and then reads the 5 test strains
    into both as well. This method will be called before many of the data-analysis methods
    """
    #preprocessHashSetECOLI(kmer_length)
    # print("test1")
    preprocessBloomFilterECOLI(kmer_length, numIntersections)
    # # print("test2")
    preprocessCountingFilterECOLI(kmer_length, numIntersections)
    # # print("test3")
    preprocessTestDataECOLI(kmer_length)
    # # print("test4")


def preprocessHashSetECOLI(kmer_length):
    global hs_final                     # mark this as global variables so we can edit them

    hs_list = readECOLI(kmer_length, "HashSet")
    hs_final = merge(0, hs_list)


def preprocessBloomFilterECOLI(kmer_length, numIntersections):
    global bf_final                    # mark this as global variables so we can edit them
    bf_list = readECOLI(kmer_length, "BloomFilter")

    bf_final = merge(numIntersections, bf_list)


def preprocessCountingFilterECOLI(kmer_length, numIntersections):
    global cf_final                    # mark this as global variables so we can edit them

    cf_list = readECOLI(kmer_length, "CountingFilter")
    cf_final = merge(numIntersections, cf_list)


def preprocessTestDataECOLI(kmer_length):
    global hs_test_list, bf_test_list  # mark these as global variables so we can edit them

    ECOLI1 = Path("../data/ECOLI/testecoli1.fas")
    ECOLI2 = Path("../data/ECOLI/testecoli2.fas")
    # ECOLI3 = Path("../data/ECOLI/testecoli3.fas")
    # ECOLI4 = Path("../data/ECOLI/testecoli4.fas")
    f1 = open(ECOLI1, "r")
    f2 = open(ECOLI2, "r")
    # f3 = open(ECOLI3, "r")
    # f4 = open(ECOLI4, "r")
    genome_test_list = [parse_file_ecoli(f1)[0], parse_file_ecoli(f2)[0]]  # , parse_file_ecoli(f3), parse_file_ecoli(f4)]

    # Fill in hs_test_list
    for i in range(len(genome_test_list)):
        hs_test_list.append(getDataStructure("HashSet"))

    for i in range(len(genome_test_list)):
        hs_test_list[i] = break_kmers(genome_test_list[i], hs_test_list[i], kmer_length)


def kmerLength_vs_hashset_size_ecoli():
    """
    This will use Pympler's asizeof.asizeof() method to find the deep size of hash sets, and use the internal bit array
    size to approximate the deep size of bloom filters
    """

    x = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    y_hs_bytes = []
    y_bf_bytes = []
    y_cf_bytes = []
    for i in x:
        preprocessBloomFilterECOLI(i, 0)
        preprocessHashSetECOLI(i)
        preprocessCountingFilterECOLI(i, 0)
        y_bf_bytes.append(bf_final.getBitSize() / 8)
        y_hs_bytes.append(asizeof.asizeof(hs_final))
        y_cf_bytes.append(cf_final.getBitSize())
    fig = plt.figure()
    plt.plot(x, y_hs_bytes, label='HashSet')
    plt.plot(x, y_bf_bytes, label='BloomFilter')
    plt.plot(x, y_cf_bytes, label='CountingFilter')
    plt.xlabel('k-mer length (nucleotides)')
    plt.ylabel('size of data structure (bytes)')
    plt.title('Impact of k-mer length on the size of HashSets, BloomFilters, and CountingFilters')
    plt.legend()
    plt.ticklabel_format(style='plain')
    plt.show()


def hashsetTimeAnalysisECOLI():
    """
    This will use Python's timeit method to find how long it takes to build the bloom filters and hash sets
    """

    SETUP_CODE = '''
from __main__ import preprocessHashSetECOLI'''

    TEST_CODE = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    TEST_CODE[0] = '''
preprocessHashSetECOLI(50)'''
    TEST_CODE[1] = '''
preprocessHashSetECOLI(60)'''
    TEST_CODE[2] = '''
preprocessHashSetECOLI(70)'''
    TEST_CODE[3] = '''
preprocessHashSetECOLI(80)'''
    TEST_CODE[4] = '''
preprocessHashSetECOLI(90)'''
    TEST_CODE[5] = '''
preprocessHashSetECOLI(100)'''
    TEST_CODE[6] = '''
preprocessHashSetECOLI(110)'''
    TEST_CODE[7] = '''
preprocessHashSetECOLI(120)'''
    TEST_CODE[8] = '''
preprocessHashSetECOLI(130)'''

    hashSetTimes = [0, 1, 2, 3, 4, 5, 6, 7, 8]

    for i in range(9):
        times = timeit.repeat(setup=SETUP_CODE, stmt=TEST_CODE[i], number=3)
        hashSetTimes[i] = min(times)
    return hashSetTimes


def bloomfilterTimeAnalysisECOLI():
    """
    This will use Python's timeit method to find how long it takes to build the bloom filters and hash sets
    """

    SETUP_CODE = '''
from __main__ import preprocessBloomFilterECOLI'''

    TEST_CODE = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    TEST_CODE[0] = '''
preprocessBloomFilterECOLI(50, 0)'''
    TEST_CODE[1] = '''
preprocessBloomFilterECOLI(60, 0)'''
    TEST_CODE[2] = '''
preprocessBloomFilterECOLI(70, 0)'''
    TEST_CODE[3] = '''
preprocessBloomFilterECOLI(80, 0)'''
    TEST_CODE[4] = '''
preprocessBloomFilterECOLI(90, 0)'''
    TEST_CODE[5] = '''
preprocessBloomFilterECOLI(100, 0)'''
    TEST_CODE[6] = '''
preprocessBloomFilterECOLI(110, 0)'''
    TEST_CODE[7] = '''
preprocessBloomFilterECOLI(120, 0)'''
    TEST_CODE[8] = '''
preprocessBloomFilterECOLI(130, 0)'''

    bloomFilterTimes = [0, 1, 2, 3, 4, 5, 6, 7, 8]

    for i in range(9):
        times = timeit.repeat(setup=SETUP_CODE, stmt=TEST_CODE[i], number=3)
        bloomFilterTimes[i] = min(times)
    return bloomFilterTimes


def countingfilterTimeAnalysisECOLI():
    """
    This will use Python's timeit method to find how long it takes to build the bloom filters and hash sets
    """
    SETUP_CODE = '''
from __main__ import preprocessCountingFilterECOLI'''

    TEST_CODE = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    TEST_CODE[0] = '''
preprocessCountingFilterECOLI(50, 0)'''
    TEST_CODE[1] = '''
preprocessCountingFilterECOLI(60, 0)'''
    TEST_CODE[2] = '''
preprocessCountingFilterECOLI(70, 0)'''
    TEST_CODE[3] = '''
preprocessCountingFilterECOLI(80, 0)'''
    TEST_CODE[4] = '''
preprocessCountingFilterECOLI(90, 0)'''
    TEST_CODE[5] = '''
preprocessCountingFilterECOLI(100, 0)'''
    TEST_CODE[6] = '''
preprocessCountingFilterECOLI(110, 0)'''
    TEST_CODE[7] = '''
preprocessCountingFilterECOLI(120, 0)'''
    TEST_CODE[8] = '''
preprocessCountingFilterECOLI(130, 0)'''


    countingFilterTimes = [0, 1, 2, 3, 4, 5, 6, 7, 8]

    for i in range(9):
        times = timeit.repeat(setup=SETUP_CODE, stmt=TEST_CODE[i], number=3)
        countingFilterTimes[i] = min(times)
    return countingFilterTimes


def compareTimeAnalysesEcoli():
    """
        This function calls the TimeAnalysis() functions for all the data
        structures, collects timing data vs kmer size, and plots them together
    """
    kmerSizeList = [50, 60, 70, 80, 90, 100, 110, 120, 130]

    fig = plt.figure()
    bftaec = bloomfilterTimeAnalysisECOLI()
    print(bftaec)
    cftaec = countingfilterTimeAnalysisECOLI()
    print(cftaec)
    #plt.plot(kmerSizeList, hashsetTimeAnalysisECOLI(), label='HashSet')
    plt.plot(kmerSizeList, bftaec, label='BloomFilter')
    plt.plot(kmerSizeList, cftaec, label='CountingFilter')
    plt.xlabel('k-mer length (nucleotides)')
    plt.ylabel('time to construct (seconds)')
    plt.title('Impact of k-mer length on construction time of HashSets, BloomFilters, and CountingFilters')
    plt.legend()
    plt.ticklabel_format(style='plain')
    plt.show()


def accuracyAnalysisECOLI():
    for i in range(1, 3, 1):
        preprocessTestDataECOLI(130)
        print("number of intersections: ", i)
        preprocessBloomFilterECOLI(130, i)
        for j in range(len(hs_test_list)):
            kmers_contained_in_bf = 0
            for key, value in hs_test_list[j].getDictionary().items():
                if bf_final.contains(key):
                    kmers_contained_in_bf += 1
            print("Test Strain Similarity to Bloom Filter/Counting Filter, Strain", j + 1, ":", kmers_contained_in_bf/hs_test_list[j].getSize())
        hs_test_list.clear()

#Finds whether a k-mer is common across strains of E. coli
def mark_test_strain_kmers_ecoli(kmer, intersections):
    global preprocess_bloom_filter
    if not preprocess_bloom_filter:
        preprocessBloomFilterECOLI(len(kmer), intersections)
        preprocess_bloom_filter = True
    if bf_final.contains(kmer):
        return True
    else:
        return False
