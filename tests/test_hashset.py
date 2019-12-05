import sys
sys.path.insert(1, '../')
import unittest
import random
from src.data_structures.hashset import HashSet


class TestHashSet(unittest.TestCase):

    def test_test(self):
        self.assertTrue(True)
        self.assertFalse(False)

    def test_insert_contains(self):
        s = HashSet()
        tempList = []
        for i in range(1000):
            tempList.append(i)

        random.shuffle(tempList)
        for i in tempList:
            s.insert(i)

        for i in range(1000):
            self.assertTrue(s.contains(i))

        for i in range(-1, -1000, -1):
            self.assertFalse(s.contains(i))

    def test_get_count(self):
        s = HashSet()
        for i in range(100):
            for j in range(i):
                s.insert(i)

        for i in range(100):
            self.assertEqual(s.getCount(i), i)

    def test_union1(self):
        a = HashSet()
        b = HashSet()
        for i in range(1000):
            a.insert(i)

        for i in range(1000, 2000):
            b.insert(i)

        a.union(b)
        for i in range(2000):
            self.assertEqual(a.getCount(i), 1)

    def test_union2(self):
        a = HashSet()
        b = HashSet()
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
        a = HashSet()
        b = HashSet()
        for i in range(1000):
            a.insert(i)

        for i in range(1000, 2000):
            b.insert(i)

        a.intersection(b)
        for i in range(2000):
            self.assertEqual(a.getCount(i), 0)

    def test_intersection2(self):
        a = HashSet()
        b = HashSet()
        for i in range(1000):
            a.insert(i)

        for i in range(500, 1000):
            b.insert(i)

        a.intersection(b)
        for i in range(500):
            self.assertEqual(a.getCount(i), 0)
        for i in range(500, 1000, 1):
            self.assertEqual(a.getCount(i), 1)

    def test_intersection3(self):
        a = HashSet()
        b = HashSet()
        c = HashSet()
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







