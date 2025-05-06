import sys
import time
import argparse
import itertools
import os
from lib import basictools
from lib import allparsers
import logging

### $transcriptREF $proteinREF $CDSFramefasta
### R shiny inputs
### input files as radio buttons | gencode [19] [38] or NCBI [...] or own
### output name
# get the start time

def change_header(original_header, prot_trans):

	NCBI_list = [">lcl",">NP_",">XP_",">YP_",">NM_",">NR_",">XM_", ">XR_"]

	if (original_header[1:6] == "ENST0") or (original_header[1:6] == "ENSTR"): # ensembl gencode transcript # or ensembl X Y chromosome annotation # or ensembl db

		if "|" in original_header: # gencode

			ENSTID = (original_header.split("|")[0].split(".")[0])[1:]

			Gene = original_header.split("|")[5]

			if prot_trans == "prot":
				return "ENST", ENSTID
			else:
				if "CDS:" in original_header:
					CDS = [x for x in original_header.split("|") if "CDS:" in x][0]
					return "ENST", ENSTID, Gene, CDS
				else:
					return "unknown"

		else:
			ENSTID = (original_header.split(" ")[0].split(".")[0])[1:]
			listline = original_header.split(" ")
			#print(listline)
			if len(listline) > 6:

				if "gene_symbol" in listline[6]:
					Gene = listline[6].split("gene_symbol:")[1]
					return "ENST", ENSTID, Gene
				else:
					return "unknown"
			else:
				return "unknown"
                
	elif (original_header[1:6] == "ENSP0"):
		if "|" in original_header: # gencode
			ENSTID = (original_header.split("|")[1].split(".")[0])
			return "ENST", ENSTID
		else: # ensembl db
			listline = original_header.split(" ")
			#print(listline)
			ENSTID = listline[4].split("transcript:")[1]
			return "ENST", ENSTID
            
	elif original_header[:4] in NCBI_list: # NCBI refseq transcript # cds_from_genomic.fna

		#print("")
		#print(prot_trans, original_header, "-------------")
		if prot_trans == "prot":
			ID = original_header.split(" ")[0][1:]

			#print("NCBI", ID)
			return "NCBI", ID
		else:
			if ("protein_id=" in original_header) and ("gene=" in original_header):

				ID = original_header.split("|")[1].split("_cds")[0]

				NCBIparts = [part for part in original_header.split(" ") if part.startswith("[")]

				Gene = [x for x in NCBIparts if "gene=" in x][0].split("gene=")[1][:-1]

				proteinID = [x for x in NCBIparts if "protein_id=" in x][0].split("protein_id=")[1][:-1]

				#print("NCBI",proteinID, Gene, ID)
				return "NCBI", proteinID, Gene, ID

			elif ("protein_id=" not in original_header) and ("gene=" not in original_header): 
				ID = original_header[1:].split(" ")[0]
				Gene = original_header.split("(")[1].split(")")[0]
				return "NCBI",ID, Gene

			else:
				return "unknown"


	else:
	#	print(header)
		ID = original_header.split("|")[0][1:]
		Gene = original_header.split("|")[1]

		if prot_trans == "trans":

			if "CDS:" in original_header:
				CDS = original_header.split("|")[2]
				return "unknown", ID, Gene, CDS
			else:
				return "unknown"
		else:
			return "unknown", ID, Gene


def open_fasta(infilename, prot_trans):

	infile = open(infilename, "r")
	lines = infile.readlines()
	linedict = {}

	for line in lines:
		line = line.strip()

		if line != "":
			if line.startswith(">"):
#				print(line, prot_trans)
				header = list(change_header(line,prot_trans))
				header = prot_trans + "_" + "|".join(header)
				#print("1.2", header)

				seqs = []
			else:
				seq = line

				if seq not in seqs:
					seqs.append(seq)
					linedict[header] = ["".join(seqs)]
	return linedict


def validate(transcriptREFname,peptideREFname):

	validated_dict= {}
	not_found_keys = []

	transcriptdict = open_fasta(transcriptREFname, "trans")  ### shiny input or default
	#print("1.2")
	#print(transcriptdict)

	AAdict = open_fasta(peptideREFname, "prot")	   ### shiny input or default

	for header, seq in transcriptdict.items():
		seq = seq[0]

		if "NCBI" in header.split("|")[0]:
			codonseq = basictools.dna_to_codon(seq)
			outkey = ">" + header.split("|")[1] + "|" + header.split("|")[2] + "|" + header.split("|")[-1]
			#print(codonseq, outkey, "!!!!!!!")

#>NP_001357096.1|SYTL4|NM_001370167.1

		else:
			if "CDS:" in header:
				CDSstart = int(header.split("|")[-1].split("CDS:")[1].split("-")[0])
				CDSend = int(header.split("|")[-1].split("CDS:")[1].split("-")[1]) 

				CDSseq = seq[CDSstart-1:CDSend]
				codonseq = basictools.dna_to_codon(CDSseq)

				outkey = ">" + "|".join(header.split("|")[1:])
				seq = seq[CDSstart-1:]

			else:
				codonseq = basictools.dna_to_codon(seq)
				outkey = ">" + "|".join(header.split("|")[1:])
				#seq = seq[CDSstart-1:]

		codon2AA = basictools.codon_AA_dictionary()

		if codonseq != []:
			if (basictools.codon_to_AA(codon2AA,[codonseq[0]], True)[0] == "M") and (basictools.codon_to_AA(codon2AA,[codonseq[-1]],True)[0] == "*"):

				Codon_AAseq = basictools.codon_to_AA(codon2AA,codonseq,False)[0]
				#print(Codon_AAseq)
				prot_key = "|".join(header.split("|")[:2]).replace("trans","prot")
				#print(prot_key)
				#print(AAdict.keys())
                
                
				if prot_key in AAdict.keys():
					Protein_AAseq = AAdict[prot_key][0]
					#print("herrerer")

					if Codon_AAseq == Protein_AAseq:
						validated_dict[outkey] = seq
				else:
					not_found_keys.append(prot_key)

#	print(validated_dict)

	#print("not found keys", len(not_found_keys))

	validated_dict2 = {}

	for key, v in validated_dict.items():
		if key[:-1] == "|":
			key2 = key[:-1]
			validated_dict2[key2] = v
		else:
			validated_dict2[key] = v

	return validated_dict2


def get_canonicals(valdict):

	fastadict = valdict #NCBget_fastadict(infilename)
	canonical_dict = {}
	canonicals = []

	for header in fastadict.keys():

		NCBI_list = [">lcl",">NP_",">XP_",">YP_"]

		if header.startswith(">ENST0"):

			if "CDS:" in header:
				CDS = header.split("CDS:")[1]
				length = int(CDS.split("-")[1]) - int(CDS.split("-")[0]) + 1
			else:
				length = len(fastadict[header])

		else:
			 length = len(fastadict[header])

		Gene = header.split("|")[1]
		ID = header.split("|")[0][1:]

		if Gene not in canonical_dict.keys():
			canonical_dict[Gene] = {}
			canonical_dict[Gene][length] = [ID]
		else:
			if length not in canonical_dict[Gene].keys():
				canonical_dict[Gene][length] = [ID]
			else:
				canonical_dict[Gene][length].append(ID)

	#print(canonical_dict)

	for Gene in canonical_dict.keys():
		lengths = canonical_dict[Gene].keys()
		highest =  sorted(lengths)[-1]
		canonical = canonical_dict[Gene][highest][0]

		canonicals.append(canonical)

	return canonicals # list of longest transcript ID per Gene



def main():

	starttime = time.time()

	logging.basicConfig(level = logging.INFO , format='{asctime} - {levelname} - {message}',filemode='w', style = "{", datefmt = "[%d-%m-%Y] [%H:%M] ")

	parser = argparse.ArgumentParser()
	allparsers.parser_get_validated(parser)
	args = parser.parse_args()

	logging.info("run started")

	validateddict = validate(args.transcriptREF, args.proteinREF) # compare transcript file to protein file # ensure both are the same AA sequence, ensuring they start with 'ATG' end with stop codon. no prior frameshift occured etc

	if args.aminoacid_outlevel:
		outlevel = "mRNA"
	else:
		outlevel = ""

	if args.CanonicalFilter:
    
		infotext = "Number of IDs in file before being canonical " +  str(len(validateddict))
		logging.info(infotext)

		canonicals = get_canonicals(validateddict) # get list of isoforms of same gene and their lengths, longest is considered canonical
		fastadict = basictools.filter_fasta_dict(canonicals,validateddict ) # filter fasta file for only canonicals

		infotext = "Number of IDs in file after being canonical " +  str(len(fastadict))
		logging.info(infotext)

	else:

		infotext = "Number of IDs in file " +  str(len(validateddict))
		logging.info(infotext)

		fastadict = validateddict # keep all isoforms etc

	duplicated_headers = basictools.duplicated_seqs(fastadict.keys()) # find duplicated headers !!!
	#print(duplicated_headers)
	if duplicated_headers != []:
		warningtext = "Warning !!! headers " + duplicated_headers + " were found more than once in your file"
		logging.warning(warningtext)
		warningtext = "duplicated headers " +  str(len(duplicated_headers)) + " out of all header " + str(len(fastadict.keys()))
		logging.warning(warningtext)

	duplicated_sequences = basictools.duplicated_seqs(fastadict.values()) # how many duplicates sequences 
	if duplicated_sequences != []:
		warningtext =  "duplicated sequences " + str(len(duplicated_sequences)) + " out of all sequences " + str(len(fastadict.values()))
		logging.warning(warningtext)

	outfilename = os.getcwd() + "/" + args.outfilename.split("/")[-1]

	basictools.write_fasta_to_file(outfilename, fastadict, outlevel) # write dict to outfile

	endtime = time.time()
	elapsed_time = endtime - starttime
	infotext = 'Execution time: ' + str( elapsed_time/60)  +  ' minutes'
	logging.info(infotext)

main()

