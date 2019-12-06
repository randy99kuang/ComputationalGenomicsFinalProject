<<<<<<< HEAD
import matplotlib.pyplot as plt
import timeit
from src.analysis.cluster import *
from pympler import asizeof

sys.path.insert(1, '../../')

hs_final = None
bf_final = None
cf_final = None
hs_test_list = []
bf_test_list = []
cf_test_list = []


def preprocessAllHIV(kmer_length, numIntersections):
    """
    Method that reads the 32 HIV strains into both a hashset and a bloom filter, and then reads the 5 test strains
    into both as well. This method will be called before many of the data-analysis methods
    """

    preprocessHashSetHIV(kmer_length)
    preprocessBloomFilterHIV(kmer_length, numIntersections)
    preprocessCountingFilterHIV(kmer_length, numIntersections)
    preprocessTestDataHIV(kmer_length)


def preprocessHashSetHIV(kmer_length):
    global hs_final                     # mark this as global variables so we can edit them

    hs_list = readHIV(kmer_length, "HashSet")
    hs_final = merge(0, hs_list)


def preprocessBloomFilterHIV(kmer_length, numIntersections):
    global bf_final                    # mark this as global variables so we can edit them

    bf_list = readHIV(kmer_length, "BloomFilter")
    bf_final = merge(numIntersections, bf_list)


def preprocessCountingFilterHIV(kmer_length, numIntersections):
    global bf_final                    # mark this as global variables so we can edit them

    bf_list = readHIV(kmer_length, "CountingFilter")
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

    TEST_CODE = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    TEST_CODE[0] = '''
preprocessHashSetHIV(5)'''
    TEST_CODE[1] = '''
preprocessHashSetHIV(10)'''
    TEST_CODE[2] = '''
preprocessHashSetHIV(30)'''
    TEST_CODE[3] = '''
preprocessHashSetHIV(50)'''
    TEST_CODE[4] = '''
preprocessHashSetHIV(70)'''
    TEST_CODE[5] = '''
preprocessHashSetHIV(100)'''
    TEST_CODE[6] = '''
preprocessHashSetHIV(120)'''
    TEST_CODE[7] = '''
preprocessHashSetHIV(150)'''
    TEST_CODE[8] = '''
preprocessHashSetHIV(190)'''
    TEST_CODE[9] = '''
preprocessHashSetHIV(250)'''
    TEST_CODE[10] = '''
preprocessHashSetHIV(350)'''
    TEST_CODE[11] = '''
preprocessHashSetHIV(500)'''

    hashSetTimes = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

    for i in range(12):
        times = timeit.repeat(setup=SETUP_CODE, stmt=TEST_CODE[i], number=3)
        hashSetTimes[i] = min(times)
    return hashSetTimes


def bloomfilterTimeAnalysis():
    """
    This will use Python's timeit method to find how long it takes to build the bloom filters and hash sets
    """

    SETUP_CODE = '''
from __main__ import preprocessBloomFilterHIV'''

    TEST_CODE = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    TEST_CODE[0] = '''
preprocessBloomFilterHIV(5, 0)'''
    TEST_CODE[1] = '''
preprocessBloomFilterHIV(10, 0)'''
    TEST_CODE[2] = '''
preprocessBloomFilterHIV(30, 0)'''
    TEST_CODE[3] = '''
preprocessBloomFilterHIV(50, 0)'''
    TEST_CODE[4] = '''
preprocessBloomFilterHIV(70, 0)'''
    TEST_CODE[5] = '''
preprocessBloomFilterHIV(100, 0)'''
    TEST_CODE[6] = '''
preprocessBloomFilterHIV(120, 0)'''
    TEST_CODE[7] = '''
preprocessBloomFilterHIV(150, 0)'''
    TEST_CODE[8] = '''
preprocessBloomFilterHIV(190, 0)'''
    TEST_CODE[9] = '''
preprocessBloomFilterHIV(250, 0)'''
    TEST_CODE[10] = '''
preprocessBloomFilterHIV(350, 0)'''
    TEST_CODE[11] = '''
preprocessBloomFilterHIV(500, 0)'''

    bloomFilterTimes = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

    for i in range(12):
        times = timeit.repeat(setup=SETUP_CODE, stmt=TEST_CODE[i], number=3)
        bloomFilterTimes[i] = min(times)
    return bloomFilterTimes

def countingfilterTimeAnalysis():
    """
    This will use Python's timeit method to find how long it takes to build the bloom filters and hash sets
    """

    SETUP_CODE = '''
from __main__ import preprocessCountingFilterHIV'''

    TEST_CODE = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    TEST_CODE[0] = '''
preprocessCountingFilterHIV(5, 0)'''
    TEST_CODE[1] = '''
preprocessCountingFilterHIV(10, 0)'''
    TEST_CODE[2] = '''
preprocessCountingFilterHIV(30, 0)'''
    TEST_CODE[3] = '''
preprocessCountingFilterHIV(50, 0)'''
    TEST_CODE[4] = '''
preprocessCountingFilterHIV(70, 0)'''
    TEST_CODE[5] = '''
preprocessCountingFilterHIV(100, 0)'''
    TEST_CODE[6] = '''
preprocessCountingFilterHIV(120, 0)'''
    TEST_CODE[7] = '''
preprocessCountingFilterHIV(150, 0)'''
    TEST_CODE[8] = '''
preprocessCountingFilterHIV(190, 0)'''
    TEST_CODE[9] = '''
preprocessCountingFilterHIV(250, 0)'''
    TEST_CODE[10] = '''
preprocessCountingFilterHIV(350, 0)'''
    TEST_CODE[11] = '''
preprocessCountingFilterHIV(500, 0)'''

    countingFilterTimes = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

    for i in range(12):
        times = timeit.repeat(setup=SETUP_CODE, stmt=TEST_CODE[i], number=3)
        countingFilterTimes[i] = min(times)
    return countingFilterTimes


def compareTimeAnalyses():
    """
        This function calls the TimeAnalysis() functions for all the data
        structures, collects timing data vs kmer size, and plots them together
    """
    kmerSizeList = [5, 10, 30,  50, 70, 100, 120, 150, 190, 250, 350, 500]

    fig = plt.figure()
    plt.plot(kmerSizeList, hashsetTimeAnalysis(), label='HashSet')
    plt.plot(kmerSizeList, bloomfilterTimeAnalysis(), label='BloomFilter')
    plt.plot(kmerSizeList, countingfilterTimeAnalysis(), label='CountingFilter')
    plt.xlabel('k-mer length (nucleotides)')
    plt.ylabel('time to construct (seconds)')
    plt.title('Impact of k-mer length on construction time of HashSets, BloomFilters, and CountingFilters')
    plt.legend()
    plt.ticklabel_format(style='plain')
    plt.show()


def accuracyAnalysisHIV():
    for i in range(1, 6, 1):
        print("number of intersections: ", i)
        preprocessAllHIV(100, i)
        hashDictionary = hs_final.getDictionary()
        totalCorrect = 0
        total = len(hashDictionary)
        falseNegativeCount = 0
        falseNegativeSum = 0
        for key, value in hashDictionary.items():
            if value >= 2 ** i:
                if bf_final.contains(key):
                    totalCorrect += 1
                else:
                    # print("Error: bloom filter does not contain key when it should, of value:", value)
                    falseNegativeCount += 1
                    falseNegativeSum += value
            else:
                if not bf_final.contains(key):
                    totalCorrect += 1
                # else:
                #   print("Error: bloom filter does contain key when it shouldn't, of value:", value)

        print(totalCorrect / total)
        print("false negative average:", falseNegativeSum / falseNegativeCount)
=======
import matplotlib.pyplot as plt
import timeit
from src.analysis.cluster import *
from pympler import asizeof

sys.path.insert(1, '../../')

hs_final = None
bf_final = None
hs_test_list = []
bf_test_list = []


def preprocessAllECOLI(kmer_length, numIntersections):
    """
    Method that reads the 32 ECOLI strains into both a hashset and a bloom filter, and then reads the 5 test strains
    into both as well. This method will be called before many of the data-analysis methods
    """

    preprocessHashSetECOLI(kmer_length)
    preprocessBloomFilterECOLI(kmer_length, numIntersections)
    preprocessTestDataECOLI(kmer_length)


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

    bf_list = readECOLI(kmer_length, "CountingFilter")
    bf_final = merge(numIntersections, cf_list)


def preprocessTestDataECOLI(kmer_length):
    global hs_test_list, bf_test_list  # mark these as global variables so we can edit them

    ECOLI5 = Path("../data/ECOLI/ECOLI5Test.fasta")
    f = open(ECOLI5, "r")
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
        preprocessBloomFilterECOLI(i, 0)
        preprocessHashSetECOLI(i)
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
from __main__ import preprocessHashSetECOLI'''

    TEST_CODE = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    TEST_CODE[0] = '''
preprocessHashSetECOLI(5)'''
    TEST_CODE[1] = '''
preprocessHashSetECOLI(10)'''
    TEST_CODE[2] = '''
preprocessHashSetECOLI(30)'''
    TEST_CODE[3] = '''
preprocessHashSetECOLI(50)'''
    TEST_CODE[4] = '''
preprocessHashSetECOLI(70)'''
    TEST_CODE[5] = '''
preprocessHashSetECOLI(100)'''
    TEST_CODE[6] = '''
preprocessHashSetECOLI(120)'''
    TEST_CODE[7] = '''
preprocessHashSetECOLI(150)'''
    TEST_CODE[8] = '''
preprocessHashSetECOLI(190)'''
    TEST_CODE[9] = '''
preprocessHashSetECOLI(250)'''
    TEST_CODE[10] = '''
preprocessHashSetECOLI(350)'''
    TEST_CODE[11] = '''
preprocessHashSetECOLI(500)'''

    hashSetTimes = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

    for i in range(12):
        times = timeit.repeat(setup=SETUP_CODE, stmt=TEST_CODE[i], number=3)
        hashSetTimes[i] = min(times)
    return hashSetTimes


def bloomfilterTimeAnalysis():
    """
    This will use Python's timeit method to find how long it takes to build the bloom filters and hash sets
    """

    SETUP_CODE = '''
from __main__ import preprocessBloomFilterECOLI'''

    TEST_CODE = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    TEST_CODE[0] = '''
preprocessBloomFilterECOLI(5, 0)'''
    TEST_CODE[1] = '''
preprocessBloomFilterECOLI(10, 0)'''
    TEST_CODE[2] = '''
preprocessBloomFilterECOLI(30, 0)'''
    TEST_CODE[3] = '''
preprocessBloomFilterECOLI(50, 0)'''
    TEST_CODE[4] = '''
preprocessBloomFilterECOLI(70, 0)'''
    TEST_CODE[5] = '''
preprocessBloomFilterECOLI(100, 0)'''
    TEST_CODE[6] = '''
preprocessBloomFilterECOLI(120, 0)'''
    TEST_CODE[7] = '''
preprocessBloomFilterECOLI(150, 0)'''
    TEST_CODE[8] = '''
preprocessBloomFilterECOLI(190, 0)'''
    TEST_CODE[9] = '''
preprocessBloomFilterECOLI(250, 0)'''
    TEST_CODE[10] = '''
preprocessBloomFilterECOLI(350, 0)'''
    TEST_CODE[11] = '''
preprocessBloomFilterECOLI(500, 0)'''

    bloomFilterTimes = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

    for i in range(12):
        times = timeit.repeat(setup=SETUP_CODE, stmt=TEST_CODE[i], number=3)
        bloomFilterTimes[i] = min(times)
    return bloomFilterTimes

def countingfilterTimeAnalysis():
    """
    This will use Python's timeit method to find how long it takes to build the bloom filters and hash sets
    """

    SETUP_CODE = '''
from __main__ import preprocessCountingFilterECOLI'''

    TEST_CODE = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    TEST_CODE[0] = '''
preprocessCountingFilterECOLI(5, 0)'''
    TEST_CODE[1] = '''
preprocessCountingFilterECOLI(10, 0)'''
    TEST_CODE[2] = '''
preprocessCountingFilterECOLI(30, 0)'''
    TEST_CODE[3] = '''
preprocessCountingFilterECOLI(50, 0)'''
    TEST_CODE[4] = '''
preprocessCountingFilterECOLI(70, 0)'''
    TEST_CODE[5] = '''
preprocessCountingFilterECOLI(100, 0)'''
    TEST_CODE[6] = '''
preprocessCountingFilterECOLI(120, 0)'''
    TEST_CODE[7] = '''
preprocessCountingFilterECOLI(150, 0)'''
    TEST_CODE[8] = '''
preprocessCountingFilterECOLI(190, 0)'''
    TEST_CODE[9] = '''
preprocessCountingFilterECOLI(250, 0)'''
    TEST_CODE[10] = '''
preprocessCountingFilterECOLI(350, 0)'''
    TEST_CODE[11] = '''
preprocessCountingFilterECOLI(500, 0)'''

    countingFilterTimes = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

    for i in range(12):
        times = timeit.repeat(setup=SETUP_CODE, stmt=TEST_CODE[i], number=3)
        countingFilterTimes[i] = min(times)
    return countingFilterTimes


def compareTimeAnalyses():
    """
        This function calls the TimeAnalysis() functions for all the data
        structures, collects timing data vs kmer size, and plots them together
    """
    kmerSizeList = [5, 10, 30,  50, 70, 100, 120, 150, 190, 250, 350, 500]

    fig = plt.figure()
    plt.plot(kmerSizeList, hashsetTimeAnalysis(), label='HashSet')
    plt.plot(kmerSizeList, bloomfilterTimeAnalysis(), label='BloomFilter')
    plt.plot(kmerSizeList, countingfilterTimeAnalysis(), label='CountingFilter')
    plt.xlabel('k-mer length (nucleotides)')
    plt.ylabel('time to construct (seconds)')
    plt.title('Impact of k-mer length on construction time of HashSets, BloomFilters, and CountingFilters')
    plt.legend()
    plt.ticklabel_format(style='plain')
    plt.show()


def accuracyAnalysisECOLI():
    for i in range(1, 6, 1):
        print("number of intersections: ", i)
        preprocessAllECOLI(100, i)
        hashDictionary = hs_final.getDictionary()
        totalCorrect = 0
        total = len(hashDictionary)
        falseNegativeCount = 0
        falseNegativeSum = 0
        for key, value in hashDictionary.items():
            if value >= 2 ** i:
                if bf_final.contains(key):
                    totalCorrect += 1
                else:
                    # print("Error: bloom filter does not contain key when it should, of value:", value)
                    falseNegativeCount += 1
                    falseNegativeSum += value
            else:
                if not bf_final.contains(key):
                    totalCorrect += 1
                # else:
                #   print("Error: bloom filter does contain key when it shouldn't, of value:", value)

        print(totalCorrect / total)
        print("false negative average:", falseNegativeSum / falseNegativeCount)
>>>>>>> cebd255d41d5afa87686720c71fdc98798316f85
