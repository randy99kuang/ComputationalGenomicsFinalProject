from src.data_structures.set_abc import Set


class HashSet(Set):

    def __init__(self):
        self.elements = {}

    def getDictionary(self):
        return self.elements

    def getSize(self):
        return len(self.elements)

    def contains(self, x):
        return x in self.elements

    def insert(self, x):
        if not self.contains(x):
            self.elements[x] = 1
        else:
            self.elements[x] += 1

    def getCount(self, x):
        if self.contains(x):
            return self.elements[x]
        else:
            return 0

    def union(self, b):
        bDict = b.getDictionary()
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
