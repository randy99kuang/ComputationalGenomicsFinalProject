import matplotlib.pyplot as plt
import timeit
from src.analysis.cluster import *
from pympler import asizeof

sys.path.insert(1, '../../')

hs_final = None
bf_final = None
hs_test_list = []
bf_test_list = []


def preprocessAllHIV(kmer_length, numIntersections):
    """
    Method that reads the 32 HIV strains into both a hashset and a bloom filter, and then reads the 5 test strains
    into both as well. This method will be called before many of the data-analysis methods
    """

    preprocessHashSetHIV(kmer_length)
    preprocessBloomFilterHIV(kmer_length, numIntersections)
    preprocessTestDataHIV(kmer_length)


def preprocessHashSetHIV(kmer_length):
    global hs_final                     # mark this as global variables so we can edit them

    hs_list = readHIV(kmer_length, "HashSet")
    hs_final = merge(0, hs_list)


def preprocessBloomFilterHIV(kmer_length, numIntersections):
    global bf_final                    # mark this as global variables so we can edit them

    bf_list = readHIV(kmer_length, "BloomFilter")
    bf_final = merge(numIntersections, bf_list)


def preprocessTestDataHIV(kmer_length):
    global hs_test_list, bf_test_list  # mark these as global variables so we can edit them

    hiv5 = Path("../data/HIV/hiv5Test.fasta")
    f = open(hiv5, "r")
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


def kmerLength_vs_hashset_size():
    """
    This will use Pympler's asizeof.asizeof() method to find the deep size of hash sets, and use the internal bit array
    size to approximate the deep size of bloom filters
    """

    x = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    y_hs_bytes = []
    y_bf_bytes = []
    for i in x:
        preprocessBloomFilterHIV(i, 0)
        preprocessHashSetHIV(i)
        y_bf_bytes.append(bf_final.getBitSize() / 8)
        y_hs_bytes.append(asizeof.asizeof(hs_final))
    fig = plt.figure()
    plt.plot(x, y_hs_bytes, label='HashSet')
    plt.plot(x, y_bf_bytes, label='BloomFilter')
    plt.xlabel('k-mer length (nucleotides)')
    plt.ylabel('size of data structure (bytes)')
    plt.title('Impact of k-mer length on the size of HashSets and BloomFilters')
    plt.legend()
    plt.ticklabel_format(style='plain')
    plt.show()


def hashsetTimeAnalysis():
    """
    This will use Python's timeit method to find how long it takes to build the bloom filters and hash sets
    """

    SETUP_CODE = ''' 
from __main__ import preprocessHashSetHIV'''

    TEST_CODE = ''' 
preprocessHashSetHIV(100)'''

    times = timeit.repeat(setup=SETUP_CODE, stmt=TEST_CODE, number=10)
    print('Hash set creation and merge time: {}'.format(min(times)))


def bloomfilterTimeAnalysis():
    """
    This will use Python's timeit method to find how long it takes to build the bloom filters and hash sets
    """

    SETUP_CODE = ''' 
from __main__ import preprocessBloomFilterHIV'''

    TEST_CODE = ''' 
preprocessBloomFilterHIV(100, 0)'''

    times = timeit.repeat(setup=SETUP_CODE, stmt=TEST_CODE, number=10)
    print('Bloom filter creation and merge time: {}'.format(min(times)))





