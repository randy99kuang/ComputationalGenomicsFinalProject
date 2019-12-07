import math
import mmh3
from array import array
from src.data_structures.set_abc import Set


class CountingFilter(Set):
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

        # Counting array of given size
        self.counting_array = [0] * self.size


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
            self.counting_array[digest] = self.counting_array[digest] + 1
            # Add to each element in bit_array
        # self.bit_array = [x+1 for x in self.bit_array[digests]]

    def delete(self, item):
        for i in range(self.hash_count):
            digest = mmh3.hash(item, i) % self.size
            self.counting_array[digest] = self.counting_array[digest] - 1

    def contains(self, item):
        """
        Check for existence of an item in filter
        """
        for i in range(self.hash_count):
            digest = mmh3.hash(item, i) % self.size
            if self.counting_array[digest] == 0:
                # if any of bit is False then,its not present
                # in filter
                # else there is probability that it exist
                return False
        return True

    def howMany(self, item):
        low = float("inf")
        for i in range(self.hash_count):
            digest = mmh3.hash(item, i) % self.size
            if self.counting_array[digest] < low:
                low = self.counting_array[digest]
        return low

    def getCountingArray(self):
        return self.counting_array

    def union(self, b):
        """
        b  :  counting_filter to union with

        unions self and b, and stores the union in self
        b will be automatically discarded by garbage collector
        """
        mp = b.getCountingArray()
        for index in range(self.size):
            self.counting_array[index] += mp[index]

        return self

    def intersection(self, b):
        """
        b  :  counting_filter to intersect with

        intersects self and b, and stores the intersection in self
        b will be automatically discarded by garbage collector
        """
        mp = b.getCountingArray()
        for index in range(self.size):
            if self.counting_array[index] < mp[index]:
                continue
            elif mp[index] < self.counting_array[index]:
                self.counting_array[index] = mp[index]
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
