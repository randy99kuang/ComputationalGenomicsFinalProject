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
bf_test_list = []


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
    global cf_final                    # mark this as global variables so we can edit them

    cf_list = readHIV(kmer_length, "CountingFilter")
    cf_final = merge(numIntersections, cf_list)


def preprocessTestDataHIV(kmer_length):
    global hs_test_list, bf_test_list  # mark these as global variables so we can edit them

    hiv5 = Path("../data/HIV/hiv5Test.fasta")
    f = open(hiv5, "r")
    genome_test_list = parse_file(f)

    # for i in range(len(genome_test_list)):
    #     bf_test_list[i] = break_kmers(genome_test_list[i], bf_test_list[i], kmer_length)

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
    y_cf_bytes = []
    for i in x:
        preprocessBloomFilterHIV(i, 0)
        preprocessCountingFilterHIV(i, 0)
        preprocessHashSetHIV(i)
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


def preprocessHashSetHIVStrainTime(kmer_length, numStrains):
    global hs_final                     # mark this as global variables so we can edit them

    hs_list = readHIV(kmer_length, "HashSet")
    hs_final = merge(0, hs_list[0:numStrains])


def preprocessBloomFilterHIVStrainTime(kmer_length, numIntersections, numStrains):
    global bf_final                    # mark this as global variables so we can edit them

    bf_list = readHIV(kmer_length, "BloomFilter")
    bf_final = merge(numIntersections, bf_list[0:numStrains])


def preprocessCountingFilterHIVStrainTime(kmer_length, numIntersections, numStrains):
    global cf_final                    # mark this as global variables so we can edit them

    cf_list = readHIV(kmer_length, "CountingFilter")
    cf_final = merge(numIntersections, cf_list[0:numStrains])

def compareStrainTimeAnalysis():
    SETUP_CODE_HS = '''
from __main__ import preprocessHashSetHIVStrainTime'''
    SETUP_CODE_BF = '''
from __main__ import preprocessBloomFilterHIVStrainTime'''
    SETUP_CODE_CF = '''
from __main__ import preprocessCountingFilterHIVStrainTime'''

    #k-mers of length 100
    TEST_CODE_HS = [0, 1, 2, 3]
    TEST_CODE_HS[0] = '''
preprocessHashSetHIVStrainTime(100, 32)'''
    TEST_CODE_HS[1] = '''
preprocessHashSetHIVStrainTime(100, 16)'''
    TEST_CODE_HS[2] = '''
preprocessHashSetHIVStrainTime(100, 8)'''
    TEST_CODE_HS[3] = '''
preprocessHashSetHIVStrainTime(100, 4)'''

    TEST_CODE_BF = [0, 1, 2, 3]
    TEST_CODE_BF[0] = '''
preprocessBloomFilterHIVStrainTime(100, 0, 32)'''
    TEST_CODE_BF[1] = '''
preprocessBloomFilterHIVStrainTime(100, 0, 16)'''
    TEST_CODE_BF[2] = '''
preprocessBloomFilterHIVStrainTime(100, 0, 8)'''
    TEST_CODE_BF[3] = '''
preprocessBloomFilterHIVStrainTime(100, 0, 4)'''

    TEST_CODE_CF = [0, 1, 2, 3]
    TEST_CODE_CF[0] = '''
preprocessCountingFilterHIVStrainTime(100, 0, 32)'''
    TEST_CODE_CF[1] = '''
preprocessCountingFilterHIVStrainTime(100, 0, 16)'''
    TEST_CODE_CF[2] = '''
preprocessCountingFilterHIVStrainTime(100, 0, 8)'''
    TEST_CODE_CF[3] = '''
preprocessCountingFilterHIVStrainTime(100, 0, 4)'''


    numStrains = [32, 16, 8, 4]
    hashSetTimes = [0, 1, 2, 3]
    bloomFilterTimes = [0, 1, 2, 3]
    countingFilterTimes = [0, 1, 2, 3]

    for i in range(4):
        timesHs = timeit.repeat(setup=SETUP_CODE_HS, stmt=TEST_CODE_HS[i], number=3)
        hashSetTimes[i] = min(timesHs)
        timesBf = timeit.repeat(setup=SETUP_CODE_BF, stmt=TEST_CODE_BF[i], number=3)
        bloomFilterTimes[i] = min(timesBf)
        timesCf = timeit.repeat(setup=SETUP_CODE_CF, stmt=TEST_CODE_CF[i], number=3)
        countingFilterTimes[i] = min(timesCf)

    fig = plt.figure()
    plt.plot(numStrains, hashSetTimes, label='HashSet', marker='o')
    plt.plot(numStrains, bloomFilterTimes, label='BloomFilter', marker='o')
    plt.plot(numStrains, countingFilterTimes, label='CountingFilter', marker='o')
    plt.xlabel('number of strains of HIV')
    plt.ylabel('time to construct (seconds)')
    plt.title('Impact of strain count on merge time of HashSets, BloomFilters, and CountingFilters with length 100 kmers')
    plt.legend()
    plt.ticklabel_format(style='plain')
    plt.show()


def false_positive_accuracy():
    count = 0
    for i in range(10000):
        if bf_final.contains(str(i)):
            count += 1
    return count/10000


def accuracyAnalysisHIV():
    bloomFilterStrainSimilarity = [1, 0, 1, 2, 3, 4]
    hashSetStrainSimilarity = [1, 0, 1, 2, 3, 4]
    bloomFilterAccuracy = [1, 0, 1, 2, 3, 4]
    bloomFilterFalseNegative = [1, 0, 1, 2, 3, 4]
    bloomFilterFalsePositive = [0, 1, 2, 3, 4, 5]
    for i in range(1, 6, 1):
    # print("number of intersections: ", i)
        preprocessAllHIV(100, i)
        hashDictionary = hs_final.getDictionary()
        totalCorrect = 0
        total = len(hashDictionary)
        falseNegativeCount = 0
        falseNegativeSum = 0

        bloomFilterFalsePositive[i] = false_positive_accuracy()

        for key, value in hashDictionary.items():
            if value >= 2:
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
        bloomFilterAccuracy[i] = totalCorrect / total
        bloomFilterFalseNegative[i] = falseNegativeSum / falseNegativeCount

        #print("Bloom filter accuracy:", totalCorrect / total)
        #print("Bloom filter false negative average count:", falseNegativeSum / falseNegativeCount)
        bloomSum = 0
        hashSum = 0
        for j in range(len(hs_test_list)):
            kmers_contained_in_bf = 0
            kmers_contained_in_hash = 0
            for key, value in hs_test_list[j].getDictionary().items():
                # if value >= 2:
                if bf_final.contains(key):
                    kmers_contained_in_bf += 1
                if hs_final.contains(key):
                    kmers_contained_in_hash += 1
            bloomSum += kmers_contained_in_bf/hs_test_list[j].getSize()
            hashSum += kmers_contained_in_hash / hs_test_list[j].getSize()
            #print("Test Strain Similarity to Bloom Filter, Strain", j + 1, ":", kmers_contained_in_bf/hs_test_list[j].getSize())
            #print("Test Strain Similarity to Hash, Strain", j + 1, ":", kmers_contained_in_hash / hs_test_list[j].getSize())
        bloomFilterStrainSimilarity[i] = bloomSum / 5
        hashSetStrainSimilarity[i] = hashSum / 5
        hs_test_list.clear()
    bloomSum = 0
    hashSum = 0
    preprocessAllHIV(100, 0)
    bloomFilterFalsePositive[0] = false_positive_accuracy()
    for j in range(len(hs_test_list)):
        kmers_contained_in_bf = 0
        kmers_contained_in_hash = 0
        for key, value in hs_test_list[j].getDictionary().items():
            # if value >= 2:
            if bf_final.contains(key):
                kmers_contained_in_bf += 1
            if hs_final.contains(key):
                kmers_contained_in_hash += 1
        bloomSum += kmers_contained_in_bf/hs_test_list[j].getSize()
        hashSum += kmers_contained_in_hash / hs_test_list[j].getSize()
    bloomFilterStrainSimilarity[0] = bloomSum / 5
    hashSetStrainSimilarity[0] = hashSum / 5
    hs_test_list.clear()


    intersectionX = [0, 1, 2, 3, 4, 5]

    fig1 = plt.figure()
    plt.plot(intersectionX, bloomFilterStrainSimilarity, marker='o')
    plt.xlabel('number of intersections')
    plt.ylabel('bloom filter average strain similarity')
    plt.title('Average strain similarity for bloom filter versus number of intersections')
    plt.legend()
    plt.ticklabel_format(style='plain')
    plt.show()

    fig2 = plt.figure()
    plt.plot(intersectionX, hashSetStrainSimilarity, marker='o')
    plt.xlabel('number of intersections')
    plt.ylabel('hash set average strain similarity')
    plt.title('Average strain similarity for hash set versus number of intersections')
    plt.legend()
    plt.ticklabel_format(style='plain')
    plt.show()

    fig3 = plt.figure()
    plt.plot(intersectionX, bloomFilterAccuracy, marker='o')
    plt.xlabel('number of intersections')
    plt.ylabel('bloom filter accuracy')
    plt.title('Bloom filter accuracy versus number of intersections')
    plt.legend()
    plt.ticklabel_format(style='plain')
    plt.show()

    fig4 = plt.figure()
    plt.plot(intersectionX, bloomFilterFalseNegative, marker='o')
    plt.xlabel('number of intersections')
    plt.ylabel('bloom filter false negatives')
    plt.title('Bloom filter false negatives versus number of intersections')
    plt.legend()
    plt.ticklabel_format(style='plain')
    plt.show()

    fig5 = plt.figure()
    plt.plot(intersectionX, bloomFilterFalsePositive, marker='o')
    plt.xlabel('number of intersections')
    plt.ylabel('bloom filter false negatives')
    plt.title('Bloom filter false positives versus number of intersections')
    plt.legend()
    plt.ticklabel_format(style='plain')
    plt.show()


def mark_test_strain_kmers(kmer):
    global preprocess_bloom_filter
    if not preprocess_bloom_filter:
        preprocessBloomFilterHIV(len(kmer))
        preprocess_bloom_filter = True
    if bf_final.contains(kmer):
        return True
    else:
        return False
