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
