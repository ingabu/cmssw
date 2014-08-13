#!/usr/bin/python

from sys import argv
import os

if not os.path.exists('dataCard_original.txt'): os.system('mv dataCard.txt dataCard_original.txt')

mu = argv[1]

def modifyfile(filename):
    f = open(filename, 'rU')
    lines = f.readlines()
    f.close()
    with open('dataCard.txt', "w") as fout:
		for line in lines:
			if line.startswith("rate"):
				data = line.strip().split()
				i=0
				for item in data:
					if (i==2 or i==12 or i==22 or i==32):
						itemtemp = float(mu)*float(item)
						fout.write(str(itemtemp))
						fout.write("\t")
					else:
						fout.write(item)
						fout.write("\t")
					i+=1
				fout.write("\n")
			else:
				fout.write(line)
				
modifyfile('dataCard_original.txt')
