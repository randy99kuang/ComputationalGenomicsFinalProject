from src.data_structures.set_abc import Set


class HashSet(Set):

    def __init__(self):
        self.elements = {}

    def contains(self, x):
        return x in self.elements

    def insert(self, x):
        if not self.contains(x):
            self.elements[x] = 1
        else:
            self.elements[x] += 1

    def intersect(self, other):
        keyIntersect = set(self.keys()) & set(other.keys())
        ret = {}
        for i in keyIntersect:
            ret[i] = self[i] + other[i]
        return ret

    def union(self, other):
        keyUnion = self.keys().union(other.keys())
        ret = {}
        for i in keyUnion:
            #Get value if exists, add 0 if it doesn't exist
            ret[i] = self.get(i, 0) + other.get(i, 0)
        return ret
