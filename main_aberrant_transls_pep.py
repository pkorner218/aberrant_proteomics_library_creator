import time
import re
import argparse
import itertools
import sys
import os
import math
import random
import tempfile
import logging

from lib import basictools
from lib import allparsers
from lib import aberranttranslation



def exec_function(input_key, *params):

	func_dict = {
	"frameshift":aberranttranslation.frameshift,
	"substitution":aberranttranslation.substitution,
	}

	func = func_dict.get(input_key)

	return func(*params)



def find_all_Regex_matches(pattern, string):

	pat = re.compile(pattern)
	pos = 0
	out2 = []

	while (match := pat.search(string, pos)) is not None:
		pos = match.start() + 1
		out2.append(match)
	#print(out2)

	return out2

def find_aberrant_positions(ID, mainseq, args ):

	position_string_type = args.position_string_type
	position_string = args.position_string

	ab_start_poss = []
	ab_end_poss = []

	finalposs = []

	metapos = {}

	wtlength = len(basictools.returnAAseq(args.mRNA_codon_aminoacid, mainseq)) ## ab pos can not be after CDS

	if args.mRNA_codon_aminoacid == "mRNA":
		wtlength = wtlength * 3

	if position_string_type == "codon":
		ab_start_poss = [i for i, x in enumerate(mainseq) if x == position_string]
		ab_end_poss = ab_start_poss

	elif position_string_type == "aminoacid":

		for ab_pos in re.finditer(re.escape(position_string), mainseq):
			if args.force_frame:
				if basictools.check_inframe(ab_pos.start()):
					ab_start_poss.append(ab_pos.start())
					ab_end_poss.append(ab_pos.start())
			else:
				ab_start_poss.append(ab_pos.start())
				ab_end_poss.append(ab_pos.start())

	elif position_string_type  == "regex":

		for ab_pos in find_all_Regex_matches(position_string, mainseq):
			if args.force_frame:
				if basictools.check_inframe(ab_pos.start()):
					ab_start_poss.append(ab_pos.start())
					ab_end_poss.append(ab_pos.end())
			else:
				ab_start_poss.append(ab_pos.start())
				ab_end_poss.append(ab_pos.end())

	elif position_string_type  == "sequence":
        
		for ab_pos in re.finditer(re.escape(position_string), mainseq):
			if args.force_frame:
				if basictools.check_inframe(ab_pos.start()):
					ab_start_poss.append(ab_pos.start())
					ab_end_poss.append(ab_pos.end())
			else:
				print(ab_pos)
				ab_start_poss.append(ab_pos.start())
				ab_end_poss.append(ab_pos.end())

	else: ##  position_string_type is rigid position based

		if position_string_type == "startATGpos": # start codon based
			ab_start_poss = [int(args.position_string)]
			ab_end_poss = ab_start_poss
			#print(ab_start_poss,ab_end_poss, (len(mainseq)))

		if position_string_type == "STOPpos": # stop codon based
			ab_start_poss = [wtlength - int(args.position_string)]
			ab_end_poss = ab_start_poss

	for i in range(len(ab_start_poss)):
		if ab_end_poss[i] < wtlength:
			finalposs.append([ab_start_poss[i],ab_end_poss[i]])
			metapos[ab_start_poss[i]] = True
                
	return finalposs, metapos#, regexheaderseq


###########################

def spacearound(potentialposstart,potentialposend, NRbefore, NRafter, length):

	if (int(potentialposstart) > NRbefore) and (int(potentialposend) < (length - NRafter)):
		return True
	else:
		return False

############

def aberrant_header_names(header, abpos, mainseq, newseq, args):
	position = abpos[0]

	if (args.position_string_type == "sequence") and (len(args.position_string) > 10):
			positionstr = "pos" + args.position_string[:3] + "." + str(len(args.position_string)-6) + "." + args.position_string[len(args.position_string)-3:] + "_" + args.mRNA_codon_aminoacid

	elif args.position_string_type == "regex":
		positionstr = "pos" + str(mainseq[abpos[0]:abpos[1]]) + "_" + str(position) + "_" + args.mRNA_codon_aminoacid

	elif "pos" not in args.position_string_type:
		positionstr = "pos" + args.position_string + "_" + str(position) + "_" + args.mRNA_codon_aminoacid
	else:
		positionstr = args.position_string + "_" + args.mRNA_codon_aminoacid

	if "aberrantsequence" in args:
		if len(args.aberrantsequence) > 10:
			aberstr = args.aberrantsequence[:3] + "." + str(len(args.aberrantsequence)-6) + "." + args.aberrantsequence[len(args.aberrantsequence)-3:]
		else:
			aberstr = args.aberrantsequence
##############

	newheaders = {}

	if args.codon_aware: ## if you wanted codon awareness it adds the codon to the header
		codon = header.split("|")[-1]
		header = "|".join(header.split("|")[:-1])
		positionstr = "".join(["pos",args.position_string,"codon",codon,"_",str(position),"_",args.mRNA_codon_aminoacid])
        
	if args.translation_error == "frameshift":
		fsrelabels = {"+1":"p1","-1":"m1"}

		for fsdir, fsseq in newseq.items():

			translation_name = args.translation_error

			if args.mRNA_codon_aminoacid == "mRNA":
				if len(basictools.returnAAseq(args.mRNA_codon_aminoacid, fsseq))*3 == int(position):
					translation_name = translation_name + "Truncation"
			else:
				if (len(basictools.returnAAseq(args.mRNA_codon_aminoacid, fsseq)) == int(position)):
					translation_name = translation_name + "Truncation"

			fslabel = fsdir.split(":")[0]

			headindex = "_".join([translation_name,fslabel,args.position_string_type,positionstr])
			newheader = "|".join([header,headindex])
			newheaders[newheader] = newseq[fsdir]

			if not args.position_string_type == "FILEpos":
				headindex = "_".join([translation_name,fslabel,args.position_string_type,positionstr])
			else:
				posstr = args.position_string_type + str(abpos[0])
				headindex = "_".join([translation_name,fslabel,args.position_string_type,args.position_string_type, str(abpos[0]), args.mRNA_codon_aminoacid ])

			newheader = "|".join([header,headindex])
			newheaders[newheader] = newseq[fsdir]


	if args.translation_error == "substitution":
		headindex = "_".join([args.translation_error, args.substitutant ,args.position_string_type,positionstr])
		newheaders["|".join([header,headindex])] = newseq

	if args.translation_error != "frameshift": # not frameshift
		newheaders["|".join([header,headindex])] = newseq

	if args.remove_truncation:
		TRremovs = []
		for key in newheaders.keys():
			if "frameshiftTruncation" in key:
				TRremovs.append(key)
		for keyname in TRremovs:
			del newheaders[keyname]

	if args.position_string_type == "FILEpos":
		removs = []
		for key in newheaders.keys():
			if "FILEpos_FILE_" in key:
				removs.append(key)
		for keyname in removs:
			del newheaders[keyname]

	return newheaders


#################################

def aberrant_translation_file(temp, temp2, args, infastadict, posdict):

	if posdict != {}: # the positions were already given, but the fastafile is not in the right format yet
		fastadict = {}

		for ID, subdict in infastadict.items():
			seq = list(infastadict[ID].values())[0]
			WTheader = list(infastadict[ID].keys())[0]
			fastadict[WTheader] = seq
	else:
		fastadict = infastadict

	for header, rnaseq in fastadict.items():

		wt_already = False

		codonseq = basictools.dna_to_codon(rnaseq)
		codon2AA = basictools.codon_AA_dictionary()
		AAseq = basictools.codon_to_AA(codon2AA, codonseq, False)[0]
		mRNA_codon_aminoacid_dict = {"mRNA":rnaseq, "codon":codonseq, "aminoacid":AAseq}
		mainseq = mRNA_codon_aminoacid_dict[args.mRNA_codon_aminoacid]
		originallength = len(mainseq)

		ID = header.split("|")[0][1:]
		Gene = header.split("|")[1]

		if posdict == {}:
			ab_poss, metapos = find_aberrant_positions(ID,mainseq, args) # positions to be found based on input
		else:

			ab_poss = posdict[header] # positions given with input file
			metapos = {}
            
			for abpos in ab_poss:
				metapos[abpos[0]] = True

		if not args.aminoacid_outlevel:
			outlevel = ""
		else:
			outlevel = args.mRNA_codon_aminoacid

		if not wt_already:
			WTdict = {"|".join([header,"WT"]):mainseq}
			wt_already = True
			WTout = WTdict

			if len(WTout.keys()) > 1:
				basictools.write_fasta_to_temp(temp, WTout, outlevel)
			else:
				basictools.write_fasta_to_temp(temp, WTdict, outlevel)

		for abpos in ab_poss: # for each pair (startpos,endpos)

			if int(abpos[0]) > 1: # and int(abpos[1] < len(codonseq)): # pos has to be between start and stop

				if args.codon_aware: ## if you wanted codon awareness. get codonseq abpos position and at it to the header
					ca_codonpos = codonseq[abpos[0]]
					cheader = header + "|" + str(ca_codonpos)
					newseq = exec_function(args.translation_error, mainseq, abpos, args, codonseq) # for FS is a dict {+1:seq, -1:seq} # for all others its the new sequence 
					aberrantdict = aberrant_header_names(cheader, abpos, mainseq, newseq, args)
				else:
					newseq = exec_function(args.translation_error, mainseq, abpos, args, codonseq)
					aberrantdict = aberrant_header_names(header, abpos, mainseq, newseq, args)

				basictools.write_fasta_to_temp(temp, aberrantdict, outlevel)
				basictools.write_meta_to_temp(temp2, "|".join([header,"WT"]), metapos[abpos[0]], abpos[0], aberrantdict)


def prepwriteout(args, temp, controldict):

	outdict = {}


	premaindict = basictools.get_fastadict(temp.name)

	maindict = premaindict

	if args.trypsin != False:
		trp_outdict = TRYPSIN(maindict)
		outdict = Trypsin_header(trp_outdict)

	if outdict == {}:
		outdict = maindict


	return outdict



#############################


################################


def TRYPSIN(maindict):

	cutAA = ["R","K"]

	outdict = {}

	for header, seq in maindict.items():
		cut_sites = [0]

		if header.split("|")[-1] == "WT":
			for i in range(len(seq)-1):
				if seq[i] in cutAA and seq[i+1] != "P":
					cut_sites.append(i)

			if cut_sites[-1] != len(seq):
				cut_sites.append(len(seq))

			if len(cut_sites) > 2:
				for j in range(0, len(cut_sites) - 1):
					selected = seq[cut_sites[j]: cut_sites[j + 1]]

					if 5 <= len(selected) <= 55:
						outheader = header + "|" + str(cut_sites[j]) + "_" + str(cut_sites[j + 1])
						outdict[outheader] = selected
			else:
				if 5 <= len(seq) <= 55:
					outheader = header + "|fulllength"
					outdict[outheader] = seq

		else: ######### FS ones

			for i in range(len(seq)-1):
				if seq[i] in cutAA and seq[i+1] != "P":
					cut_sites.append(i)
			if cut_sites[-1] != len(seq):
				cut_sites.append(len(seq))

			if len(cut_sites) > 2:
				for j in range(0, len(cut_sites) - 1):
					selected = seq[cut_sites[j]: cut_sites[j + 1]]

					if 5 <= len(selected) <= 55:
						if selected not in outdict.values():
							Abpart = header.split("|")[-1]
							Abpos = int(Abpart.split("_")[4])

							if cut_sites[j] < Abpos < cut_sites[j+1]: 
								outheader = header + "|" + str(cut_sites[j]) + "_" + str(cut_sites[j+1]) + "_chimera"
								outdict[outheader] = selected
							else:
								outheader = header + "|" + str(cut_sites[j]) + "_" + str(cut_sites[j+1])
								outdict[outheader] = selected
			else:
				if 5 <= len(seq) <= 55:
					if selected not in outdict.values():
						outheader = header + "|fulllength_chimera"
						outdict[outheader] = selected

	return outdict


def Trypsin_header(trp_outdict):
	counter = 1
	outdict = {}

	for k,v in trp_outdict.items():
		k2 = k.split("|")

		if not "WT" in k2[-2]:
			ID = k2[0]
			Gene = k2[1]
			cond = k2[3]
			cuts = k2[-1]
			cond2 = cond.split("_")
			cond3 = cond2[0] + "_" + cond2[1] + "_" +  cond2[4] + "_" + cuts
			outs = "_".join([ID,Gene,cond3, str(counter)]) 
            
			outdict[outs] = v
            
		else:
			ID = k2[0]
			Gene = k2[1]
			cond = k2[3]
			cuts = k2[-1]
			outs = "_".join([ID,Gene,cond,cuts, str(counter)])
			outdict[outs] = v

		counter = counter + 1
  
	return outdict


###############################

def main():

	logging.basicConfig(level = logging.INFO , format='{asctime} - {levelname} - {message}',filemode='w', style = "{", datefmt = "[%d-%m-%Y] [%H:%M] ")

	starttime = time.time()

	parser = argparse.ArgumentParser()
	allparsers.parse_all_findaberrant(parser)
	args = parser.parse_args()


#	print(args) #
#	print("")   #

	allparsers.validate_argparse_findaberrant_inputs(args)

	temp = tempfile.NamedTemporaryFile(delete=False)
	temp2 = tempfile.NamedTemporaryFile(delete=False)
#####################

	if os.path.isfile(args.input): 	# input is a file [fasta , tsv, vcf, or mt.txt]

		fastadict = basictools.get_fastadict(args.input) # fastafile of all transcripts in dictionary { ENST : RNAseq }
		posdict = {}
		logging.info("Fasta dict created")

#####################

	aberrant_translation_file(temp, temp2, args, fastadict, posdict)

########################################################################################

	aberranttime = time.time() #

	elapsed_time = aberranttime - starttime #
	logging.info("Aberrant sequences finished")
	temp.seek(0)
	temp2.seek(0)

	outdict = prepwriteout(args, temp, {})

############################

	cwd = os.getcwd()
	outfilename = cwd + "/"+ args.outfilename.split("/")[-1]

	infotext = "Results written to " + outfilename
	logging.info(infotext)

	basictools.write_fasta_to_file(outfilename, outdict, args.aminoacid_outlevel)

	########################################################################

	endtime = time.time() #
	endtime = endtime - starttime #
	endtime = round((endtime/60),2)
	xstring =  "Pipe finished in " + str(endtime) + " m"
	logging.info(xstring )

	temp.close()
	temp2.close()

main()

