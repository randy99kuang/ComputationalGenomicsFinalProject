from src.analysis.cluster import *

ds_list = readHIV(10, "HashSet")
final = merge(0, ds_list)
print(final.getSize())

for key, value in final.getDictionary().items():
    if value > 10:
        print(key, value)
