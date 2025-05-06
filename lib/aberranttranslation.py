import re
from lib import basictools
import itertools

def exec_p1_m1_function(input_key, *params):
	func_dict = {"p1":p1_FS, "m1":m1_FS}
	func = func_dict.get(input_key)
	return func(*params)

def m1_FS(FSseq):
	m1seq = FSseq
	m1seq = m1seq[:3] + m1seq[2] + m1seq[3 :]
	return m1seq

def p1_FS(FSseq):
	p1seq = FSseq
	p1seq = p1seq[:3] + p1seq[4:]
	return p1seq

def frameshift(mainseq, abpos, args, codonseq):

	frameshifted = {}

	if args.mRNA_codon_aminoacid != "mRNA":
		lengthmainseq = len(codonseq)
	else:
		lengthmainseq = len(mainseq)
        
	if args.frameshift_direction != "both":
		FSdirections = [args.frameshift_direction]
	else:
		FSdirections = ["p1","m1"]

	#print("######", "bef",abpos, lengthmainseq, len(codonseq)) 
	#print("######", "bef",abpos,abpos[0], lengthmainseq, len(codonseq), abpos[1]) 
    
	if (int(abpos[0]) > 2) and (lengthmainseq > int(abpos[1])+1):

		#print("after",abpos[0]/3, lengthmainseq/3, abpos[1]/3,len(codonseq)/3) 

		for FSdirect in FSdirections:
			if args.mRNA_codon_aminoacid == "codon":
				seqrange = mainseq[abpos[0]-1:]
				before = mainseq[:abpos[0]-1]
				FSseqrange = "".join(seqrange)

				frameshiftedseq = exec_p1_m1_function(FSdirect, FSseqrange)
				FScodonseqrange = basictools.dna_to_codon(frameshiftedseq) ###
				frameshifted[FSdirect] = before + FScodonseqrange

			elif args.mRNA_codon_aminoacid == "mRNA":

					FSseqrange = mainseq[abpos[0]-3:]
					before = mainseq[:abpos[0]-3]

					frameshifted[FSdirect] = before + exec_p1_m1_function(FSdirect, FSseqrange)
			else:
				seqrange = codonseq[abpos[0]-1:]
				before = codonseq[:abpos[0]-1]
				FSseqrange = "".join(seqrange)

				frameshiftedseq = exec_p1_m1_function(FSdirect, FSseqrange)
				FScodonseqrange = basictools.dna_to_codon(frameshiftedseq)
				#frameshifted[FSdirect] = before + FScodonseqrange

				if args.codon_aware:
					FSdirect = FSdirect + ":" + str(codonseq[abpos[0]])

				totalseq = before + FScodonseqrange

				frameshifted[FSdirect] = basictools.returnAAseq("codon", totalseq)

	return frameshifted

# returns a dict {+1:p1seq, -1:m1seq}

###################################################

def repeat(mainseq, abpos, args, codonseq):

	matches = []

	before = mainseq[:int(abpos[0])]
	after = mainseq[int(abpos[1]):]

	if args.mRNA_codon_aminoacid == "mRNA":
		abposseq = mainseq[int(abpos[0]):int(abpos[1])]

		for ab_pos in re.finditer(re.escape(args.aberrantsequence), abposseq):
			matches.append([ab_pos.start(),ab_pos.end()])

		lastmatch = matches[-1] # consider last occurrence of the aberrant string to be the one to repeat
		repnumber = int(args.number)

		if lastmatch[1]-lastmatch[0] == 1: #for singe nucleotides repeat once less, they are already part of the seq
			repnumber = int(args.number) - 1

		repetitive = abposseq[int(lastmatch[0]):int(lastmatch[1])] * repnumber

		befab = abposseq[:int(lastmatch[0])] # position string sequence before repeat location
		aftab = abposseq[int(lastmatch[1]):] # position string sequence after repeat location
		newseq = "".join([before,befab,repetitive,aftab,after])

	if args.mRNA_codon_aminoacid == "codon":
		abposseq = mainseq[int(abpos[0])]
		repnumber = int(args.number) -1

		repetitive = basictools.dna_to_codon(abposseq * repnumber)

		newseq = before + repetitive + after


	if args.mRNA_codon_aminoacid == "aminoacid":
		mainseq = list(mainseq)
		abposseq = mainseq[int(abpos[0])]
		repnumber = int(args.number)

		repetitive = abposseq * int(repnumber)
		newseq = "".join(before + repetitive + after)

	return newseq

# returns a sequence or in codonseq a list

#####################################################


def skip(mainseq, abpos, args, codonseq):

	if args.notincludelast:

		if args.mRNA_codon_aminoacid != "codon":
			before = mainseq[:abpos[1]-1]
			after = mainseq[abpos[1] + int(args.length):]
		else:
			before = mainseq[:abpos[1]]
			after = mainseq[abpos[1] + int(args.length):]
	else:

		if args.mRNA_codon_aminoacid != "codon":
			before = mainseq[:abpos[1]]
			after = mainseq[abpos[1] + int(args.length):]
		else: # codon
			before = mainseq[:abpos[1]+1]
			after = mainseq[abpos[1]+1 + int(args.length):]

	newseq = before + after

	return newseq


####################################################

def truncate(mainseq, abpos, args, codonseq):

	if args.mRNA_codon_aminoacid == "codon":
		abpos[1] = abpos[1]+1


	newseq = mainseq[:int(abpos[1])]

	if args.notincludelast:
		newseq = mainseq[:int(abpos[0])]

	return newseq


###################################################

def altstart(mainseq, abpos, args, codonseq):

	if args.mRNA_codon_aminoacid == "codon":
		abpos[1] = abpos[1]+1

	newseq = mainseq[int(abpos[0]):]

	if args.notincludelast:
		newseq = mainseq[int(abpos[1]):]

	return newseq


#################################################

def insertion(mainseq, abpos, args, codonseq):

	if args.notincludelast:
		if args.mRNA_codon_aminoacid != "codon":

			if args.mRNA_codon_aminoacid == "mRNA":
				before = mainseq[:abpos[0]]
				after = mainseq[abpos[0]:]
			else: # amino acid
				before = mainseq[:abpos[1]-1]
				after = mainseq[abpos[0]:]
		else: # codon
			before = mainseq[:abpos[1]]
			after = mainseq[abpos[1]:]
	else:
		if args.mRNA_codon_aminoacid != "codon":
			before = mainseq[:abpos[1]]
			after = mainseq[abpos[1]:]
		else: # codon
			before = mainseq[:abpos[1]+1]
			after = mainseq[abpos[1]+1:]

	if args.mRNA_codon_aminoacid != "codon":
		inserted_afterseq = args.aberrantsequence + after
		newseq = before + inserted_afterseq
	else:
		afterseq = "".join(after)

		inserted_afterseq = basictools.dna_to_codon(args.aberrantsequence + afterseq)
		newseq = before + inserted_afterseq

	return newseq


###########################################


def reversion(mainseq, abpos, args, codonseq):

#def skip(mainseq, abpos, args):

	if args.notincludelast:

		if args.mRNA_codon_aminoacid != "codon":
			before = mainseq[:abpos[1]-1]
			after = mainseq[abpos[1] + int(args.length):]

			middle = mainseq[abpos[1]-1:abpos[1] + int(args.length)]
		else:
			before = mainseq[:abpos[1]]
			after = mainseq[abpos[1] + int(args.length):]

			middle = "".join(mainseq[abpos[1]:abpos[1] + int(args.length)])

	else:

		if args.mRNA_codon_aminoacid != "codon":
			before = mainseq[:abpos[1]]
			after = mainseq[abpos[1] + int(args.length):]

			middle = mainseq[abpos[1]:abpos[1] + int(args.length)]

		else: # codon
			before = mainseq[:abpos[1]+1]
			after = mainseq[abpos[1]+1 + int(args.length):]

			middle = "".join(mainseq[abpos[1]+1:abpos[1]+1 + int(args.length)])

	revmiddle = middle[::-1]

	if args.mRNA_codon_aminoacid == "codon":
		newseq = before + basictools.dna_to_codon(revmiddle) + after
	else:
		newseq = before + revmiddle + after

	return newseq

###########################################


def substitution(mainseq, abpos, args, codonseq):

	sub2 = args.substitutant.split("-")[-1]
	codon2AA = basictools.codon_AA_dictionary()

	AA2codon = basictools.AA_to_codon(codon2AA, sub2)
	sub2codon = AA2codon[0]

#	print(AA2codon, sub2)
#	print(mainseq, "---mainseq")
#	print(codonseq, "---codonseq")


	newseq = codonseq
	newseq[abpos[0]] = sub2codon 
#	print(codonseq, "---newseq")



	return newseq