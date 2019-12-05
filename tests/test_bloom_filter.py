import unittest
import random
from src.data_structures.bloom_filter import BloomFilter

n = 1000    # Number of items we will put in our Bloom Filter
p = 0.05    # Desired false positive rate


class TestHashSet(unittest.TestCase):
    
    def test_test(self):
        self.assertTrue(True)
        self.assertFalse(False)

    def test_insert_contains(self):
        s = BloomFilter(n, p)
        tempList = []
        for i in range(1000):
            tempList.append(i)

        random.shuffle(tempList)
        for i in tempList:
            s.insert(str(i))

        errors = 0
        for i in range(1000):
            if not s.contains(str(i)):
                errors = errors + 1
        print(errors)
        self.assertTrue(errors < 55)

        errors = 0
        for i in range(-1, -1000, -1):
            if s.contains(str(i)):
                errors = errors + 1
        print(errors)
        self.assertTrue(errors < 55)

    def test_union1(self):
        a = BloomFilter(n, p)
        b = BloomFilter(n, p)
        for i in range(1000):
            a.insert(i)

        for i in range(1000, 2000):
            b.insert(i)

        a.union(b)
        for i in range(2000):
            self.assertEqual(a.getCount(i), 1)

    def test_union2(self):
        a = BloomFilter(n, p)
        b = BloomFilter(n, p)
        for i in range(1000):
            a.insert(i)

        for i in range(500, 1000):
            b.insert(i)

        a.union(b)
        for i in range(500):
            self.assertEqual(a.getCount(i), 1)
        for i in range(500, 1000, 1):
            self.assertEqual(a.getCount(i), 2)

    def test_intersection1(self):
        a = BloomFilter(n, p)
        b = BloomFilter(n, p)
        for i in range(1000):
            a.insert(i)

        for i in range(1000, 2000):
            b.insert(i)

        a.intersection(b)
        for i in range(2000):
            self.assertEqual(a.getCount(i), 0)

    def test_intersection2(self):
        a = BloomFilter(n, p)
        b = BloomFilter(n, p)
        for i in range(1000):
            a.insert(i)

        for i in range(500, 1000):
            b.insert(i)

        a.intersection(b)
        for i in range(500):
            self.assertEqual(a.getCount(i), 0)
        for i in range(500, 1000, 1):
            self.assertEqual(a.getCount(i), 1)

    def test_intersection2(self):
        a = BloomFilter(n, p)
        b = BloomFilter(n, p)
        c = BloomFilter(n, p)
        for i in range(1000):
            a.insert(i)
            a.insert(i)

        for i in range(500, 1000):
            b.insert(i)
            b.insert(i)
            b.insert(i)

        for i in range(600, 700):
            c.insert(i)
            c.insert(i)
            c.insert(i)
            c.insert(i)

        b.intersection(c)
        a.intersection(b)
        for i in range(600):
            self.assertEqual(a.getCount(i), 0)
        for i in range(600, 700):
            self.assertEqual(a.getCount(i), 2)
        for i in range(700, 1000):
            self.assertEqual(a.getCount(i), 0)







