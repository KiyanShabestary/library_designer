# coding=utf8



def divide(path_to_ORFs,n):

	fh=open(path_to_ORFs,'r')

	i=0
	sub=1

	for line in fh.readlines():
		if line[0] == '>':
			if (i % n) == 0:
				path_to_sub_ORF= 'ORF_sequences/ORF_'+str(sub)+'.txt'
				gh=open(path_to_sub_ORF,'w')
				sub+=1
			i+=1

		gh.write(line)

		

	fh.close()

	return 0

def main():

	
	#-----------
	# CHANGE HERE FOR YOUR APPLICATION
	n=250 #ORFs_per_file 
	path_to_ORFs="ORF_sequences/Synechocystis_ORFs.txt" #Path_to_ORF_file 
	#-----------

	divide(path_to_ORFs,n)

	return 0


main()