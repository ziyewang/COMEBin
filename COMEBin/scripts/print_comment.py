#!/usr/bin/env python

# This script prints any comment in a structured and pretty way.

import sys

comm = sys.argv[1]
delim = str(sys.argv[2])

print('\n' + delim * 120)

max_len = 90

cut = comm.split(" ")
line = ""
for word in cut:
    if (len(line) + 1 + len(word)) > max_len:
        edge1 = int((120 - len(line)) / 2 - 5)
        edge2 = int(120 - edge1 - len(line) - 10)
        print(delim * 5, end='')
        print(" " * edge1, end='')
        print(line, end='')
        print(" " * edge2, end='')
        print(delim * 5)
        line = word
    else:
        line = line + " " + word
edge1 = int((120 - len(line)) / 2 - 5)
edge2 = int(120 - edge1 - len(line) - 10)
print(delim * 5, end='')
print(" " * edge1, end='')
print(line, end='')
print(" " * edge2, end='')
print(delim * 5)

print(delim * 120, '\n')