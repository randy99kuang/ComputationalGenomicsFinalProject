from data_structures.set_abc import Set


class HashSet(Set):

    def __init__(self):
        self.elements = set()

    def contains(self, x):
        return x in self.elements

    def insert(self, x):
        self.elements.add(x)

