import unittest, random
from src.data_structures.hashset import HashSet


class TestHashSet(unittest.TestCase):

    def test_test(self):
        self.assertTrue(True)
        self.assertFalse(False)

    def test_hashset(self):
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





