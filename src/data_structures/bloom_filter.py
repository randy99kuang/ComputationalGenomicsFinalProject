import math
import mmh3
from bitarray import bitarray
from src.data_structures.set_abc import Set


class BloomFilter(Set):
    """
    Bloom Filter implementation that uses bitarray internally as well as hash functions from mmh3
    """

    def __init__(self, items_count, fp_prob):
        """
        items_count : int
            Number of items expected to be stored in bloom filter
        fp_prob : float
            False Positive probability in decimal
        """

        # False positive probability in decimal
        self.fp_prob = fp_prob

        # Size of bit array to use
        self.size = self.get_size(items_count, fp_prob)

        # number of hash functions to use
        self.hash_count = self.get_hash_count(self.size, items_count)

        # Bit array of given size
        self.bit_array = bitarray(self.size)

        # initialize all bits as 0
        self.bit_array.setall(0)

    def getBitSize(self):
        return self.size

    def getHashCount(self):
        return self.hash_count

    def insert(self, item):
        """
        Add an item in the filter
        """
        digests = []
        for i in range(self.hash_count):
            # create digest for given item.
            # i work as seed to mmh3.hash() function
            # With different seed, digest created is different
            digest = mmh3.hash(item, i) % self.size
            digests.append(digest)

            # set the bit True in bit_array
            self.bit_array[digest] = True

    def contains(self, item):
        """
        Check for existence of an item in filter
        """
        for i in range(self.hash_count):
            digest = mmh3.hash(item, i) % self.size
            if not self.bit_array[digest]:
                # if any of bit is False then,its not present
                # in filter
                # else there is probability that it exist
                return False
        return True

    def getBitArray(self):
        return self.bit_array

    def union(self, b):
        """
        b  :  bloom_filter to union with

        unions self and b, and stores the union in self
        b will be automatically discarded by garbage collector
        """

        bBitArray = b.getBitArray()
        for i in range(self.size):
            if bBitArray[i]:
                self.bit_array[i] = True

        return self

    def intersection(self, b):
        """
        b  :  bloom_filter to intersect with

        intersects self and b, and stores the intersection in self
        b will be automatically discarded by garbage collector
        """

        bBitArray = b.getBitArray()
        newBitArray = bitarray(self.size)
        newBitArray.setall(0)
        for i in range(self.size):
            if bBitArray[i] and self.bit_array[i]:
                newBitArray[i] = 1

        self.bit_array = newBitArray

        return self

    @classmethod
    def get_size(self, n, p):
        """
        n : int
            number of items expected to be stored in filter
        p : float
            False Positive probability in decimal

        Return the optimal size of the bit array given the expected number of elements and the false positive rate.
        From class, we are given that for false positive rate follows the formula:
        p = (1 - e^(-hn/m))^h

        From class, we also know that the optimal number of hash functions is given by the formula:
        h = (m/n) * ln(2)

        If we substitute this value of h into the previous equation and then isolate m, we get the equation below. Note
        that a formal proof of this will be given in the paper we turn in.
        m = -(n * ln(p)) / (ln(2)^2)
        """

        m = -(n * math.log(p)) / (math.log(2) ** 2)
        return int(m)

    @classmethod
    def get_hash_count(self, m, n):
        """
        m : int
            size of bit array
        n : int
            number of items expected to be stored in filter

        Return the number of hash functions needed to minimize the false positive rate.
        From class, we are given that the formula for this is
        h = (m/n) * ln(2)
        """

        h = (m / n) * math.log(2)
        return int(h)
