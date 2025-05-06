import sys
import os
import logging
import argparse
import subprocess

from itertools import zip_longest
from lib import allparsers



def filter_aberrfasta(args):

	n = args.number_entries

	infilename = args.inputfastafile
#	path = args.inputfastafile.split("/")[-1]

	infile = open(infilename, "r")

	outfastaname = infilename.replace(".fasta","_filtered.fasta")
	filteredFasta = open(outfastaname, "w")

	lines = infile.readlines()

	for line in lines:
		line = line.strip()
		#count += 1

		if line.startswith(">"):
			header = line
		else:
			seq = line
			if header.split("|")[1] != "TTN": # or any other filters
				#print(record.id, len(record.seq))
				filteredFasta.write(header+'\n')
				filteredFasta.write(seq+'\n')

	return outfastaname

def grouper(n, iterable, padvalue=None):
	return zip_longest(*[iter(iterable)]*n, fillvalue=padvalue)

def write_batch_fastas(outfastaname, args, cwd):

	nr_batchfiles = 0 

	with open(outfastaname,"r") as filteredFasta:
		for i, g in enumerate(grouper( args.number_entries*2, filteredFasta, padvalue=""), 1):
			with open(cwd + "/dom/batch" + str(i) + ".fa", "w") as batchFasta:
				batchFasta.writelines(g)

			nr_batchfiles += 1

	return nr_batchfiles

def main():

	cwd = os.getcwd()

	logging.basicConfig(level = logging.INFO , format='{asctime} - {levelname} - {message}',filemode='w', style = "{", datefmt = "[%d-%m-%Y] [%H:%M] ")

	parser = argparse.ArgumentParser()
	allparsers.parser_separate_fasta(parser)
	args = parser.parse_args()
    
	print(args)

	outfastaname = filter_aberrfasta(args)

	nr_batchfiles = write_batch_fastas(outfastaname, args, cwd)

	logstring = "separated fasta into " + str(nr_batchfiles) + " files" 
	logging.info(logstring)

####some logging .. number of batch files i guess.




main()