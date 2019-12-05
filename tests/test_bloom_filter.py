import sys
sys.path.insert(1, '../')
import unittest
import random
import math
from src.data_structures.bloom_filter import BloomFilter

n = 1000        # Number of items we will put in our Bloom Filter
p = 0.05        # Desired false positive rate
buffer = 0.05   # Small buffer for false positive rate when unit testing


class TestHashSet(unittest.TestCase):
    
    def test_test(self):
        self.assertTrue(True)
        self.assertFalse(False)

    @classmethod
    def get_false_positive_rate(self, numElements, numBits, numHash):
        """
        From class, we are given that for false positive rate follows the formula:
        p = (1 - e^(-hn/m))^h
        """
        return (1 - math.exp(-numHash * numElements / numBits))**numHash

    def test_insert_contains(self):
        s = BloomFilter(n, p)
        tempList = []
        for i in range(1000):
            tempList.append(i)

        random.shuffle(tempList)
        for i in tempList:
            s.insert(str(i))

        # If the Bloom Filter actually contains an element, it should always say it does
        for i in range(1000):
            self.assertTrue(s.contains(str(i)))

        # If a Bloom Filter does not contain an element, there is a ~5% chance of a false positive
        errors = 0
        for i in range(-1, -1001, -1):
            if s.contains(str(i)):
                errors = errors + 1
        print(errors)
        self.assertTrue(errors < n * (p + buffer))

    def test_union1(self):
        a = BloomFilter(n, p)
        b = BloomFilter(n, p)
        for i in range(1000):
            a.insert(str(i))

        for i in range(1000, 2000):
            b.insert(str(i))

        a.union(b)

        # If the union Bloom Filter actually contains an element, it should always say it does
        for i in range(2000):
            self.assertTrue(a.contains(str(i)))

        # If a Bloom Filter does not contain an element, there is a chance of false positives, which will be calculated
        # by the get_false_positive_rate method
        errors = 0
        for i in range(-1, -1001, -1):
            if a.contains(str(i)):
                errors = errors + 1
        print("number of errors is" ,errors)
        fpr = self.get_false_positive_rate(n * 2, a.getBitSize(), a.getHashCount())
        print("false positive rate is ", fpr)
        self.assertTrue(errors < n * (fpr + buffer))

    def test_union2(self):
        a = BloomFilter(n, p)
        b = BloomFilter(n, p)
        for i in range(1000):
            a.insert(str(i))

        for i in range(500, 1000):
            b.insert(str(i))

        a.union(b)

        # If the union Bloom Filter actually contains an element, it should always say it does
        for i in range(1000):
            self.assertTrue(a.contains(str(i)))

        # If a Bloom Filter does not contain an element, there is a chance of false positives
        errors = 0
        for i in range(-1, -1001, -1):
            if a.contains(str(i)):
                errors = errors + 1
        print("number of errors is", errors)
        fpr = self.get_false_positive_rate(n, a.getBitSize(), a.getHashCount())
        print("false positive rate is", fpr)
        self.assertTrue(errors < n * (fpr + buffer))

    def test_intersection1(self):
        a = BloomFilter(n, p / 50)  # we need to use a slightly smaller fpr here so we don't have too many intersects
        b = BloomFilter(n, p / 50)
        for i in range(1000):
            a.insert(str(i))

        for i in range(1000, 2000):
            b.insert(str(i))

        a.intersection(b)

        # There are no elements in the intersection, so this bloom filter should contain nothing except for a few
        # elements caused by hash collisions
        errors = 0
        for i in range(2000):
            if a.contains(str(i)):
                errors = errors + 1

        print("number of intersections is", errors)
        self.assertTrue(errors < n * 0.01)

    def test_intersection2(self):
        a = BloomFilter(n, p)
        b = BloomFilter(n, p)
        for i in range(1000):
            a.insert(str(i))

        for i in range(500, 1000):
            b.insert(str(i))

        a.intersection(b)

        # There are no elements from 0 - 499, so in this range the bloom filter should only contain false postives
        errors = 0
        for i in range(500):
            if a.contains(str(i)):
                errors = errors + 1
        print("intersect2 errors", errors)
        fpr = self.get_false_positive_rate(500, a.getBitSize(), a.getHashCount())  # only 500 elements in intersect
        print("intersect2 fpr", fpr)
        self.assertTrue(errors < n * (fpr + buffer))

        # These elements all should be in the intersection, so they must be in the intersect bloom filter
        for i in range(500, 1000, 1):
            self.assertTrue(a.contains(str(i)))

    def test_intersection3(self):
        a = BloomFilter(n, p)
        b = BloomFilter(n, p)
        c = BloomFilter(n, p)
        for i in range(1000):
            a.insert(str(i))

        for i in range(500, 1000):
            b.insert(str(i))

        for i in range(600, 700):
            c.insert(str(i))

        b.intersection(c)
        a.intersection(b)

        fpr = self.get_false_positive_rate(100, a.getBitSize(), a.getHashCount())  # only 100 elements in intersect

        # False positive errors should be expected in this range
        errors = 0
        for i in range(600):
            if a.contains(str(i)):
                errors = errors + 1
        self.assertTrue(errors < n * (fpr + buffer))

        # 100% accuracy is expected in this range
        for i in range(600, 700):
            self.assertTrue(a.contains(str(i)))

        # False positive errors should be expected in this range
        errors = 0
        for i in range(700, 1000):
            if a.contains(str(i)):
                errors = errors + 1
        self.assertTrue(errors < n * (fpr + buffer))







