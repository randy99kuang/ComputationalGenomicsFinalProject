def hashfn(item):
    h = hash(item)
    return (1 << (h%64)) | (1 << (h/64%64))


def mask(val):
    return bin(hashfn(val))[2:]


class counting_filter(object):
    def __init__(self, genome):
        self.items = [0] * len(genome)

    def insert(self, item):
        bits = mask(item)
        for index, bit in enumerate(bits):
            if bit == '1':
                self.items[index] += 1

    def contains(self, item):
        bits = mask(item)
        for index, bit in enumerate(bits):
            if bit == '1' and self.items[index] == 0:
                return False
        return True

    def getCount(self, x):
        b = contains(x)
        if x == False:
            return 0
        bits = mask(x)
        low = float("inf")
        for index, bit in enumerate(bits):
            if bit == '1' and self.items[index] < low:
                low = self.items[index]
        return low

    def intersect(self, other):
        for index in range(len(self.items)):
            if self.items[index] < other[index]:
                continue
            elif other[index] < self.items[index]:
                self.items[index] = other[index]
        return self.items

    def union(self, other):
        for index in range(len(self.items)):
            self.items[index] += other[index]
        return self.items

    def remove(self, item):
        bits = mask(item)
        for index, bit in enumerate(bits):
            if bit == '1' and self.items[index]:
                self.items[index] -= 1
