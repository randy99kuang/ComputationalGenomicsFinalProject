def hashfn(item):
    h = hash(item)
    return (1 << (h%64)) | (1 << (h/64%64))


def mask(val):
    return bin(hashfn(val))[2:]


class counting_filter(object):
    def __init__(self, genome):
        self.items = [0] * len(genome)

    def add(self, item):
        bits = mask(item)
        for index, bit in enumerate(bits):
            if bit == '1':
                self.items[index] += 1

    def query(self, item):
        bits = mask(item)
        for index, bit in enumerate(bits):
            if bit == '1' and self.items[index] == 0:
                return False
        return True

    def remove(self, item):
        bits = mask(item)
        for index, bit in enumerate(bits):
            if bit == '1' and self.items[index]:
                self.items[index] -= 1
    
