#from src.data_structures.set_abc import Set
import sys


def parse(txt):
    count = -1
    genomeList = []
    while True:
        line = txt.readline()
        if len(line) == 0:
            # End of file
            break
        if line[0] == ">":
            count += 1
            genomeList.append("")
            continue
        # line = line.replace("-", "")
        genomeList[count] = genomeList[count] + line
    return genomeList


g = open(sys.argv[1], "r")
readList = parse(g)
test = {}
for i in range(len(readList)):
    if readList[i] in test:
        print("same string")
        print(i)
    else:
        test[readList[i]] = 1