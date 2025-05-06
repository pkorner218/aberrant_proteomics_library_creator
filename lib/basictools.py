import tempfile
import os
import itertools

def dna_to_codon(dna_seq): # turn a dna sequence from str format into a codon triplet sequence in list format

	codons = []

	for i in range(0, len(dna_seq), 3):
		codons.append(dna_seq[i:i + 3])
	return codons


def codon_AA_dictionary(): # returns a dictionary of codon-AminoAcid pairings

	codon2AA = {
                'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
                'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
                'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
                'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
                'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
                'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
                'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
                'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
                'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
                'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
                'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
                'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
                'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
                'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
                'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
                'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W' }

	return codon2AA

def AA_to_codon(codon2AA, AA): # for a given AminoAcid all codons are returned by reversing the codon2AA dictionary for this pair

	AA2codon = {}

	[AA2codon.setdefault(k,[]).append(v) for k,v in itertools.chain.from_iterable([itertools.product(vals,[key]) for key,vals in codon2AA.items()])]

	return AA2codon[AA]


def codon_to_AA(codon2AA, codon_seq, mode): # takes a codon sequence in list format and translates it into an amino acid string # if non codon triples are found a warning is returned additionally

	warning = False

	AAseq = []

	if not mode:
		for codon in codon_seq:
			if codon not in codon2AA.keys(): break

			if codon2AA[codon] == "*" : break

			AAseq.append(codon2AA[codon])
	else:
		for codon in codon_seq:
			if codon not in codon2AA.keys():
				AAseq.append("-")
				warning = True
			else:
				AAseq.append(codon2AA[codon])

	AAseq = "".join(AAseq)

	return AAseq, warning



def check_inframe(pos): # check if a given nucleotide pos is start of a codon that is inframe

	inFrame = False

	if pos % 3 == 0:
		inFrame = True

	return inFrame


def returnAAseq(mRNA_codon_aminoacid, seq):

	if mRNA_codon_aminoacid == "mRNA":
		codonseq = dna_to_codon(seq)
		codon2AA = codon_AA_dictionary()
		AAseq = codon_to_AA(codon2AA, codonseq, False)[0]

	elif mRNA_codon_aminoacid == "codon":
		codon2AA = codon_AA_dictionary()
		AAseq = codon_to_AA(codon2AA, seq, False)[0]
	else:
		AAseq = seq

	return AAseq

#####

#def str2bool(str): # turn basic str into a bool
#
#	if str == "True":
#		return True
#	if str == "False":
#		return False

def get_fastadict(infilename): # take a file in fasta format and return it as a dictionary # assumes unique fasta headers !

	main_dict = {}

	Mainfasta = open(infilename, "r")

	lines = Mainfasta.readlines()

	for line in lines:
		line = line.strip()

		if line.startswith(">"):
			header = line
			seqs = []

		else:
			seq = line
			seqs.append(seq)

		seqs2 = "".join(seqs)
		main_dict[header] = seqs2

	return main_dict


def filter_fasta_dict(filterlist, fastadict): # filter a fasta dictionary based on a list of transcript IDs

	filtereddict = {}

	for header, seq in fastadict.items():
		ID = header.split("|")[0][1:]

		if ID in filterlist:
			filtereddict[header] = seq

	return filtereddict


def write_fasta_to_temp(temp, writeoutdict, outlevel): # write a fasta dictionary to a temporrary file

	for header,seq in writeoutdict.items():

		if outlevel != "":
			AAseq = returnAAseq(outlevel, seq)
		else:
			AAseq = "".join(seq)

		temp.write(header.encode("utf-8"))
		temp.write(b"\n")
		temp.write(AAseq.encode("utf-8"))
		temp.write(b"\n")


def write_meta_to_temp(temp2, WTheader, metapos, abpos, aberrantdict): # write metainformation to a temporrary file

	for abkey in aberrantdict.keys():

		outstring = "\t".join([WTheader,abkey, str(abpos),str(metapos)])

		temp2.write(outstring.encode("utf-8"))
		temp2.write(b"\n")

def write_fasta_to_file(outfilename, writeoutdict, outlevel): # write a fasta dictionary to an outfilename

	outfile = open(outfilename, "w")

	for header,seq in writeoutdict.items():

		if outlevel != "":
			AAseq = "".join(returnAAseq(outlevel, seq))
		else:
			AAseq = "".join(seq)

		outfile.write(header)
		outfile.write("\n")
		outfile.write(AAseq)
		outfile.write("\n")

	outfile.close()

def duplicated_seqs(l): # take a list and return all found duplicates

	seen = set()
	seen_add = seen.add 	# adds all elements it doesn't know yet to seen and all other to seen_twice
	seen_twice = set(x for x in l if x in seen or seen_add(x)) 	# turn the set into a list (as requested)

	return list(seen_twice)


def obtain_basedicts(fastadict):

	WTdict ={}

	for header, sequence in fastadict.items():

		if header.split("|")[-1] == "WT":
			if header.split("|")[-2] != "WT":

				base = "|".join(header.split("|")[:-1])
				WTdict[base] = [sequence]

	return WTdict

