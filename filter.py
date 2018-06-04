#!/usr/bin/env python

file1 = "TEST_DEXSeq_IDs_ALL.txt"
with open(file1) as code:
	ID = code.readlines()
	ID_str = ''.join(ID)
code.close()

for i in ID:
    if '+' in i:
        print i
