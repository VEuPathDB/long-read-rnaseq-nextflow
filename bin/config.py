#!/usr/bin/env python
import sys

samID = sys.argv[1].replace('[', ' ').replace(']', ' ').replace(',', ' ').split()
build = sys.argv[2].strip()
Platform = sys.argv[3].strip()
SampleName = sys.argv[4].replace('[', ' ').replace(']', ' ').replace(',', ' ').split()

with open('config.txt', 'a') as f:
    for i in range(len(samID)):
        f.write(samID[i] + "," + str(build) + "," + str(Platform) + "," + SampleName[i] + "\n")


with open('Sample_names.txt', 'a') as f:
    for i in range(len(samID)):
        f.write(samID[i] + "\n")
