from abc import ABC, abstractmethod


class Set(ABC):
    def contains(self):
        pass

    def insert(self):
        pass

    def union(self, b):
        pass

    def intersection(self, b):
        pass
