import sys
from lib import basictools
import os

def validate_argparse_findaberrant_inputs(args):
#	print(args)
	valid_pst = ["codon","aminoacid","sequence", "startATGpos", "STOPpos", "FILEpos","regex"]
	fileendings = (".fa",".fasta",".fas",".tsv",".vcf")

	insequence = False

###

	if args.translation_error != "frameshift":
		args.remove_truncation = False

	if os.path.isfile(args.input): # input is a valid file

		if args.translation_error != 'mutation':
			if not (args.input.endswith(fileendings)):
				print("Invalid input file :", args.input, " file detected that is not in fasta, tsv, or vcf format")
				sys.exit(1)
                

	else: # input is supposedly a sequence # its not a valid file

		for char in [".","/",",","*","~"]:
			if char in args.input:
				print("Invalid input file :", args.input, " given filename does not seem to exist")
				sys.exit(1)
			else:
				insequence = True # input sequence is valid

				if args.translation_error == "frameshift":
					if args.position_string_type == "aminoacid":
						print("Please Note that frameshifting can not be performed without a given mRNA sequence. Please provide mRNA sequence instead.")
						sys.exit(1)



	if os.path.isfile(args.outfilename):
		print("Invalid output filename :", args.outfilename, " given outfilename does already exist")
		sys.exit(1)

	outendings = (".fa",".fasta",".fas")
	if not (args.outfilename.endswith(outendings)):
		print("Invalid output file name:", args.outfilename, " filename given that is not in .fasta .fa or .fas format")
		sys.exit(1)

	if "position_string_type" in args and args.position_string_type in valid_pst:

		if args.position_string_type == "regex":
			if args.mRNA_codon_aminoacid == "codon":
				print("regular expression can only be used on mRNA or AminoAcid level")
				sys.exit(1)


		if args.position_string_type == "codon":
			if not args.position_string in basictools.codon_AA_dictionary().keys():
				print("Invalid codon :", args.position_string)
				sys.exit(1)

			if args.mRNA_codon_aminoacid != "codon":
				print("Please be aware that position_string_type [codon] choice requires RNA_codon_aminoacid to also be [codon].\n Alternatively it is suggested to use position_string [sequence] with the codon as three letter nucleotide and the mRNA_codon_aminoacid [mRNA] option.")
				sys.exit(1)

		if args.position_string_type == "aminoacid":
			if not args.position_string in basictools.codon_AA_dictionary().values():
				print("Invalid amino acid :", args.position_string)
				sys.exit(1)

		if args.position_string_type == "sequence":
			for letter in args.position_string:
				if not letter in  basictools.codon_AA_dictionary().values():
					print("Invalid nucleoptide or amino acid found :", letter )
					sys.exit(1)

		if args.position_string_type == "FILEpos":
			if args.Reference is not None:
				if os.path.isfile(args.Reference):
					if not (args.Reference.endswith(fileendings)): 
						print("The given reference file", args.ref ,"is not a valid fasta file")
						sys.exit()
				else:
					if args.Reference not in ["hg19","hg38"]:
						print("The given reference is not valid", args.ref, "options are hg19/hg38")
						sys.exit()
			else:
				print("A Reference file -ref is required if using a position based fileinput")
				sys.exit(1)

			if args.position_string != "FILE":
				print("You chose to provide a file with positions (in -psl option). This requires the -ps option to be FILE")
				sys.exit(1)



	else:
		if args.translation_error != 'mutation' and args.translation_error != "vcf_file":
			print("Invalid input for position_string_type")
			sys.exit(1)

############


	if "mRNA_codon_aminoacid" in args and  args.mRNA_codon_aminoacid not in ["mRNA","codon","aminoacid"]:

		print("Invalid sequence level :", args.mRNA_aminoacid)
		sys.exit(1)


###################

	if insequence: # input is a sequence. check if it has proper nucleotides or amino acid letters
		if args.mRNA_codon_aminoacid == "mRNA":
			for letter in args.input:
				if not letter in ["A","C","G","T"]:
					print("Invalid nucleotide in the input sequence :", letter )
					sys.exit(1)

		if args.mRNA_codon_aminoacid == "aminoacid":
			for letter in args.input:
				if not letter in basictools.codon_AA_dictionary().values():
					print("Invalid amino acid in the input sequence :", letter )
					sys.exit(1)

####################


####################



def parse_all_findaberrant(parser):

	subparsers = parser.add_subparsers(required = True, dest="translation_error")

	p = subparsers.add_parser('frameshift', help='mRNA or Amino Acid sequence should be frameshifted')
	parser_Frameshift(p)

	p = subparsers.add_parser('substitution', help='')
	parser_substitution(p)

def parser_findaberrant_common(parser):
	parser.add_argument("-i","--input", type = str, help = "Name of input. Either a fasta file [file.fasta], tab separated pos file (ENSTID/GeneID\tpos or chr:1-2) [file.tsv], a mutation file [file.vcf] or a single input sequence (mRNA or amino acid level) [ACGTGGAGACGAT]", required = True)
	parser.add_argument("-pst","--position_string_type", type = str, help = "The type of position where aberrant translation starts. options are ([codon], [aminoacid], [sequence], [startATGpos] (from start codon), [STOPpos] (from stop codon), [FILEpos] (from given tsv file or vcf file), [regex] regular expression in quotation marks ''", required = True)
	parser.add_argument("-ff","--force_frame", help= "In case you gave a sequence as -pst argument, this will only consider occurrences were the sequence was found to be in the normal reading frame", action='store_true')
	parser.add_argument("-ca","--codon_aware", help = "If you chose the -spt aminoacid option but still want the output sequences to be aware and include the different codons", action='store_true') 
	parser.add_argument("-ps","--position_string",  type = str, help = "The actual position where aberrance starts? (e.g [TTA] for codon, aminoacid [L], sequence [CTGGTGATTGATGCG], a position [9], based on input file [FILE])", required = True)
	parser.add_argument("-mca","--mRNA_codon_aminoacid", type = str, help= "is the event to be on mRNA [mRNA], codon [codon] or amino acid [aminoacid] level. This determines also the distances and lengths of other options to be counted either in nucleotides, or codons/amino acids! If you provide any vcf or chr based input mRNA level is required.", required = True)

	parser.add_argument("-trp","--trypsin", help = "If this flag is used output amino acid sequences will be trypsinized into peptides at K/R. Peptide length is 5-55 amino acid", action='store_true') 

	parser.add_argument("-o","--outfilename",type=str, help = "name of the outfile", required = True)

	parser.add_argument("-aao", "--aminoacid_outlevel", help= "Level of your output file sequence default is aminoacid. Option changes fasta output format to mRNA sequence. Please be aware that this level can not be used in further pipeline",action='store_false')



def parser_Frameshift(parser):
	parser.add_argument("-fsd","--frameshift_direction", type=str, help = "In which direction should the frameshift be considered? [p1],[m1],[both])", required = True, default =None)
	parser.add_argument("-rt","--remove_truncation",  help = "Flag whether frameshifts resulting in immediate stop codons should be kept (default) or be removed (if -rt).", action="store_true")
	parser_findaberrant_common(parser)



def parser_get_validated(parser):
	parser.add_argument("-t","--transcriptREF", type=str, help="transcript reference file (fasta format)", required=True)
	parser.add_argument("-p","--proteinREF", type=str, help="protein reference file (fasta format)", required=True)
	parser.add_argument("-o","--outfilename", type=str, help="name of the outfile to create (ending .fa or .fasta)", required=True)
	parser.add_argument("-c", "--CanonicalFilter",help="Only use longest transript per gene", action='store_true')
	parser.add_argument("-aao", "--aminoacid_outlevel", help= "Level of your output file sequence default is mRNA. Option changes fasta output format to aminoacid/peptide sequence. Please be aware that this level can not be used for frameshifting. Other aberrant translation options are possible though.",action='store_true' )

def parser_separate_fasta(parser):
	parser.add_argument("-n", "--number_entries", type = int, help = "number of entries to split fasta file into. Number of lines / 2 (>header\sequence)", required = False, default = 60000)
	parser.add_argument("-i", "--inputfastafile",  type=str,help = "inputfile aberrant fasta", required = True)

def parser_substitution(parser):
	parser.add_argument("-sub", "--substitutant", type=str, help = "W-F", required = True)
	parser_findaberrant_common(parser)
