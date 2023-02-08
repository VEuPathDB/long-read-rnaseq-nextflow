#!/usr/bin/env python
import sys
samID = sys.argv[1].replace('[', ' ').replace(']', ' ').replace(',', ' ').split()
with open('Sample_names.txt', 'a') as f:
    for i in range(len(samID)):
        f.write(samID[i] + "\n")