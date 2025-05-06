****
# aberrant proteomics library creator
Singled out scripts from the *aberrant pipeline*. 


# separate_fasta.py
splits a fasta file in given number of subfiles, while making sure that ">headers" and sequences can not be separated.


# 1 get_validated_fasta.py
Since many translation fasta files show problematic inconsistencies to their transcript counterparts, this script ensures that both are the same.
Further problems were found with regards to annotated coding sequences not starting with ATG or ending with stop codons. The script removes transcripts with wrongly annotated CDS.
An optional filter possibility is `-C ` which keeps only the longest transcript per protein (here taken to be canonical)

# 1 Usage
`python3 get_validated_fasta.py -t [transcript referencefile] -p [protein referencefile] -o [outfilename] -c [canonical] -aao [amino acid as outlevel]`


# 2 main_aberrant_transls_pep.py
Creates frameshifted or substituted peptide libraries for mass spectronomy. 
This is a simplified version of the pipeline version that requires the aberrant translation event to happen on codon level.
Therefore parsed options position_string_type *-pst* and mRNA_codon_AminoAcid *-mca* should remain at codon.

The positionstring *-ps* should be a codon such as e.g. *TGG*.  

# 2 Usage for frameshifts
The frameshift direction option *-fsd* can be either *m1*, *p1* or *both*. 

`python3 main_aberrant_transls_pep.py frameshift -fsd [both|m1|p1] -i [input fastafile] -pst [codon] -ps [TGG] -mca [codon] -o [outfilename] -trp [trypsin]` 


# 2 Usage for substitutions
The substitutions option *-sub* should be a codon substituting to an Amino Acid. 

`python3 main_aberrant_transls_pep.py substitution -sub [TGG-F] -i [input fastafile] -pst [codon] -ps [TGG] -mca [codon] -o [outfilename] -trp [trypsin]` 

