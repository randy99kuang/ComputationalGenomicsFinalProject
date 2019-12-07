from src.data_structures.set_abc import Set


"""
HashSet implementation that uses an underlying dictionary which maps a k-mer
to the number of times it appears in the genome.
"""
class HashSet(Set):

    #Constructor
    def __init__(self):
        self.elements = {}

    def getDictionary(self):
        return self.elements

    def getSize(self):
        return len(self.elements)

    """
    Tests whether or not a kmer is present inside the hashset.
    """
    def contains(self, x):
        return x in self.elements

    """
    Since this is a set, we will check for repeats - if repeats, we will
    increment the number of k-mers that are present to have an actual count
    """
    def insert(self, x):
        #Check if contained
        if not self.contains(x):
            self.elements[x] = 1
        else:
            self.elements[x] += 1

    #Get the number of times a particular kmer x is in the set
    def getCount(self, x):
        if self.contains(x):
            return self.elements[x]
        else:
            return 0

    #Take the union of two HashSets
    def union(self, b):
        bDict = b.getDictionary()
        #We take all the keys, union the keys, and then put them all in a set
        keyIntersect = set(self.elements.keys()).union(set(bDict.keys()))
        aNew = {}
        
        for i in keyIntersect:
            aNew[i] = self.getCount(i) + b.getCount(i)

        self.elements.clear()
        bDict.clear()
        self.elements = aNew
        return self

    def intersection(self, b):
        bDict = b.getDictionary()
        keyIntersection = set(self.elements.keys()).intersection(set(bDict.keys()))
        aNew = {}
        for i in keyIntersection:
            aNew[i] = min(self.getCount(i), b.getCount(i))

        self.elements.clear()
        bDict.clear()
        self.elements = aNew
        return self
