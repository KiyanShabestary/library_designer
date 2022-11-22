# coding=utf8

import glob
import re
import os.path
import sys
#import os

# Kiyan Shabestary 13.06.2020

## Converting ORF fasta file(s) to a dictionary. Make sure fasta file has unique identifier
#
# INPUT
# path_to_ORFs			path to fasta file (.txt format) with ORFs and sequence. Make sure that gene identifier is unique.)
#
# OUTPUT
# ORFs					dictionary of ORFs. Key is gene name.
#
def get_ORFs(path_to_ORFs):

	ORFs= {}

	fh=open(path_to_ORFs,'r')

	for line in fh.readlines():
		if line[0] == '>':
			ID = line.split('|')[0]
			ORFs[ID.split('>')[1]] = ''
		else:
			ORFs[ID.split('>')[1]] += line.strip()
	
	fh.close()

	return ORFs


# Reads and stores fasta file(s) in complete_sequence folder for off-target search
#
# INPUT
# path_to_complete_seq	path to fasta file (.txt format) with sequence. 
#
# OUTPUT
# genome				string of the complete genomic sequence
#
def get_complete_seq(path_to_complete_seq):

	genome= ''

	for filename in glob.glob(os.path.join(path_to_complete_seq, '*.txt')):
		with open(os.path.join(os.getcwd(), filename), 'r') as fh: 
			for line in fh.readlines():
				if line[0] != '>': genome += line.strip()
		fh.close()

	return genome


## Finding gRNAs respecting hard constraints (GC_content, bad_seed) using regular expressions. 
#
# INPUT
# ORFs 					dictionary of ORFs. Key is gene name.
# complete_seqs			dictionary of complete sequences. Key is file name.
# N_min 				minimal gRNA length
# N_max 				maximal gRNA length
# bad_seed				list of strings defining bad_seed sequence(s) not desired in gRNA (see Cui et al., 2018 & Yao et al., 2020)
# forbidden_site		list of strings defining sequences not wanted in gRNA for downstream processing (restriction sites)
# GC_min				float defining minimal GC content in gRNA
# GC_max				float defining maximal GC content in gRNA
# 
# OUTPUT
# gRNAs					two-dimensional dictionary of all possible photospacer regions. Key #1 is gene name. Key #2 is gRNA sequence. Item is gRNA position within ORF.
#
def find_gRNAs(ORFs, l_min, l_max, bad_seeds, forbidden_sites, GC_min, GC_max):
	gRNAs={}

	for ORF in ORFs.keys():
		search={}

		for N in range(l_min,l_max+1):
			for m in re.finditer(r'(?=(CC).{%s})'%(str(N)),ORFs[ORF]):
				#print '%02d-%02d: %s' % (m.start(), m.end(), m.group(1))
				#print '%02d-%02d: %s' % (m.start(), m.start()+N+2, ORFs[ORF][int(m.start()):int(m.start()+N+2)])
				if is_good(ORFs[ORF][int(m.start()):int(m.start()+N+2)], bad_seeds, forbidden_sites, GC_min, GC_max): search[ORFs[ORF][int(m.start()):int(m.start()+N+2)]]=m.start()
		
		#if search: gRNAs[ORF]=search
		gRNAs[ORF]=search

	return gRNAs


## Finding offtargets for a given gRNA. Offtargets are screened on forward and reverse strands and for NGG (default PAM) and NAG (alternative PAM) sites.
#
# INPUT
# gRNA 					string with gRNA sequence to be tested
# genome 				string of genome to test gRNA against
# n 					integer representing the length of the region of the gRNA that is checked for off-targets
# m 					integer representing the number of mismatches to consider a site an offtarget location
# 
# OUTPUT
# mismatches			interger representing the number of mismatch allowed to consider a site to be offtarget 		
#
def offtargets(gRNA, genome, reverse_genome, n, m):
	offtargets=0
	

	for start in range(len(genome)-(n+1)):#Checking only n-first bases for CCN PAM (default library PAM)
		position_CCN=0
		mismatch_CCN=0
		position_NGG=0
		mismatch_NGG=0
		position_CTN=0
		mismatch_CTN=0
		position_NAG=0
		mismatch_NAG=0
		gRNA_CTN = 'CTN'+gRNA[3:len(gRNA)]

		##CCN forward search
		while (position_CCN<n and mismatch_CCN<m):
			if (position_CCN == 0 or position_CCN == 1) and (genome[start+position_CCN] != gRNA[position_CCN]): mismatch_CCN = 20#mismatch in PAM is critical for offtarget binding. Mismatch set to an arbitrary high value
			if ((genome[start+position_CCN] != gRNA[position_CCN]) and position_CCN != 2): mismatch_CCN += 1 #Position 3 of the sgRNA is not taken into account for offtarget screen
			position_CCN += 1 #Update position
		if mismatch_CCN<m: offtargets +=1 

		##NGG reverse search
		while (position_NGG<n and mismatch_NGG<m):
			if (position_NGG == 0 or position_NGG == 1) and (reverse_genome[start+position_NGG] != gRNA[position_NGG]): mismatch_NGG = 20
			if ((reverse_genome[start+position_NGG] != gRNA[position_NGG]) and position_NGG != 2): mismatch_NGG += 1
			position_NGG += 1 #Update position
		if mismatch_NGG<m: offtargets +=1 

		##CTN forward search
		while (position_CTN<n and mismatch_CTN<m):
			if (position_CTN == 0 or position_CTN == 1) and (genome[start+position_CTN] != gRNA[position_CTN]): mismatch_CTN = 20
			if ((genome[start+position_CTN] != gRNA_CTN[position_CTN]) and position_CTN != 2): mismatch_CTN += 1
			position_CTN += 1 #Update position
		if mismatch_CTN<m: offtargets +=1 

		##NAG reverse search
		while (position_NAG<n and mismatch_NAG<m):
			if (position_NAG == 0 or position_NAG == 1) and (reverse_genome[start+position_NAG] != gRNA[position_NAG]): mismatch_NAG = 20
			if ((reverse_genome[start+position_NAG] != gRNA_CTN[position_NAG]) and position_NAG != 2): mismatch_NAG += 1
			position_NAG += 1 #Update position
		if mismatch_NAG<m: offtargets +=1 


	return offtargets


## Select gRNAs based on quality, position, and off-targets. Calls offtargets function to get offtargets. 
#
# INPUT
# library				library at previous step. If first step, keys are ORFs but empty
# gRNAs					dictionary of gRNAs
# spacing				integer defining the spacing between sgRNAs
# N 					targeted number of gRNA per ORF
# genome 				complete sequence to check offtarget against
# reverse_genome		complete sequence to check offtarget against
#
# OUTPUT
# library				dictionary of filtered gRNAs
#
def select_gRNAs(library, gRNAs, spacing, N, genome, reverse_genome, n, m):

	for ORF in gRNAs.keys():
		library[ORF]={}
		for gRNA in sorted(gRNAs[ORF], key=gRNAs[ORF].get, reverse=False): # sorts gRNA from binding distance

			if len(library[ORF])==0: # First sgRNA, double if to reduce computations
				if offtargets(gRNA, genome, reverse_genome, n, m)==1: 
					library[ORF][gRNA]=gRNAs[ORF][gRNA]
					#print library[ORF][gRNA]

			elif is_not_within(gRNAs[ORF][gRNA],library[ORF],spacing):
				if offtargets(gRNA, genome, reverse_genome, n, m)==1: 
					library[ORF][gRNA]=gRNAs[ORF][gRNA]
					#print library[ORF][gRNA]

			if len(library[ORF])==N: break
				
	return library


# Makes reverse strand
def reverse_strand(fwd_strand):
	rev_strand = ''
	# find complementary 3'-5'
	for base in fwd_strand:
		if base == 'A':
			rev_strand = rev_strand + 'T'
		if base == 'G':
			rev_strand = rev_strand + 'C'
		if base == 'C':
			rev_strand = rev_strand + 'G'
		if base == 'T':
			rev_strand = rev_strand + 'A'
		if base == 'N':
			rev_strand = rev_strand + 'N'
		if base == 'X':
			rev_strand = rev_strand + 'X'
	# translate to 5'-3'
	return rev_strand[::-1]

def is_not_within(position_to_test, entries, minimal_distance):
	for entry in entries.keys():
		if abs(entries[entry]-position_to_test)<=minimal_distance: return False 
	
	return True

def get_GC(sequence):
	GC = 0
	position = 0

	for base in sequence:
		if position >=3:
			if str(base) == 'G' or str(base) == 'C': GC += 1
		position += 1

	return (float(GC)/(len(sequence)-3))

def is_good(gRNA,bad_seeds,forbidden_sites,GC_min,GC_max):

	if ((get_GC(gRNA) < GC_min) or (get_GC(gRNA) > GC_max)): return False

	for bad_seed in bad_seeds:
		if bad_seed in gRNA: return False

	for forbidden_site in forbidden_sites:
		if forbidden_site in gRNA: return False

	return True

def make_file(library,path_to_library):

	fh=open(path_to_library,'w')
	
	for ORF in library.keys():
		for gDNA in library[ORF].keys():
			fh.write('>%s|%s \nXXXXXX%sYYYYYY\n' % (ORF,str(library[ORF][gDNA]),gDNA))

	fh.close()

	return 0

def main():

	path_to_ORFs = 'ORF_sequences/ORF_'+str(sys.argv[1])+'.txt'
	path_to_library = 'output/library_ORF_'+str(sys.argv[1])+'.txt'
	path_to_complete_seq = 'complete_sequences/'

	l_min=18 #gRNA minimal length (without PAM)
	l_max=23 #gRNA maximal length (without PAM)
	bad_seeds=['ACCCA','ATACT','TGGAA','GGGGGG','TTTT']
	forbidden_sites=[]#Add your forbidden sequence here as above
	GC_min=0.4
	GC_max=0.8
	n=15 #bp checked for offtargets on sgRNA
	m=1 #mismatch allowed to consider a target an offtarget
	N=5 # number of gRNA per ORF
	spacing=5 # minimal spacing between gRNAs within an ORF

	## get ORFs and their sequences
	ORFs=get_ORFs(path_to_ORFs)

	## get complete sequences (ORF and intergenic region) for off-target screen
	genome=get_complete_seq(path_to_complete_seq)
	reverse_genome=reverse_strand(genome)# For search on the reverse strand

	## find all possible gRNAs for each ORF
	gRNAs=find_gRNAs(ORFs,l_min,l_max,bad_seeds,forbidden_sites,GC_min,GC_max)

 	## select gRNAs for the library (1st round)
 	library={} #initialising library
	library=select_gRNAs(library, gRNAs, spacing, N, genome, reverse_genome, n, m)

	print library

	## Writes library in a fasta format
	make_file(library,path_to_library)

	return 0

main()


# def select_gRNAs(library, gRNAs, spacing, N, perc_in_ORF, dist_in_ORF, ORFs):

# 	for ORF in gRNAs.keys():

# 		#Put if loop checking off-target and removing if not good
# 		while(offtargets(library[ORF][min(gRNAs[ORF], key=gRNAs[ORF].get)])>0)


		
# 		library[ORF]={}
# 		library[ORF][min(gRNAs[ORF], key=gRNAs[ORF].get)]=gRNAs[ORF][min(gRNAs[ORF], key=gRNAs[ORF].get)]


# 		for gRNA in gRNAs[ORF].keys():
			
# 			if len(library[ORF])==N: break #Get out of the loop if enough gRNAs for this ORF
# 			if is_not_within(gRNAs[ORF][gRNA],library[ORF],spacing): 
# 				if offtargets(gRNA) == 0: library[ORF][gRNA]=gRNAs[ORF][gRNA] #keep good sgRNA

# 			#else: #Case no sgRNA yet -> assign first from the list

# 			#print gRNAs[ORF][gRNA]

				
# 	return library

## Testing the current library to check if we have number of gRNAs above a minimal treshold for each ORF. Prints results in a log
#	
# INPUT
# library 				current library to be tested

# def test(library):

# 	" write log file 'a' mode"

# 	return acceptance