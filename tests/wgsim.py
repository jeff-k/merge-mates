import sys
import pysam

mismatches = 0
total = 0

for record in pysam.FastxFile(sys.argv[1]):
    _, start, end, *k = record.name.split('_')
    length = len(record.sequence)
    gap = (int(end) - int(start)) + 1
    if length != gap and length != (gap - 1):
        mismatches += 1
        print("wrong: {} ({} (-1))".format(length, gap))
        print(record.name)
    total += 1

print("{}/{} incorrect overlaps".format(mismatches, total))
