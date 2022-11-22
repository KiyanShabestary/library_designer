
import random

def get_GC(sequence):
	GC = 0
	position = 0

	for base in sequence:
		if position >=3:
			if str(base) == 'G' or str(base) == 'C': GC += 1
		position += 1

	return (float(GC)/(len(sequence)-3))

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

def is_good(gRNA,bad_seeds,forbidden_sites,GC_min,GC_max):

	if ((get_GC(gRNA) < GC_min) or (get_GC(gRNA) > GC_max)): return False

	for bad_seed in bad_seeds:
		if (bad_seed in gRNA) or (reverse_strand(bad_seed) in gRNA): return False

	for forbidden_site in forbidden_sites:
		if forbidden_site in gRNA: return False

	return True

def main():
	N=10 #numbers of random gRNAs to generate
	l_min=18 #gRNA minimal length
	l_max=23 #gRNA maximal length 
	GC_min=0.4
	GC_max=0.8
	bad_seeds=['ACCCA','ATACT','TGGAA','GGGGGG','TTTT']
	forbidden_sites=[]

	fh=open('output/random_gRNAs.txt','w')
	rdm_gRNAs=[]


	while(len(rdm_gRNAs)<N):
		rdm_gRNA=''

		n=random.randint(l_min,l_max)

		for i in range(0,n):
			base=random.randint(1,4)
			if base == 1: rdm_gRNA+='A'
			if base == 2: rdm_gRNA+='T'
			if base == 3: rdm_gRNA+='G'
			if base == 4: rdm_gRNA+='C'

		if is_good(rdm_gRNA,bad_seeds,forbidden_sites,GC_min,GC_max): 
			rdm_gRNAs.append(rdm_gRNA)
			fh.write('>Ctrl%s\n%s\n' % (str(len(rdm_gRNAs)),str(rdm_gRNA)))

	fh.close()

	return 0


main()