#!/usr/bin/env python
import sys

annotationFile = sys.argv[1]
with open(annotationFile, 'r') as f:
    first_line = f.readline()
    names = set()
    for l in f:
        lst = l.split()[1]
        names.add(lst)
    with open('Sample_names.txt', 'a') as g:
        for name in names:
            g.write(name + "\n")
