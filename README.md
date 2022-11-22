# library_designer
This script generates a CRISPRi library for any organism of interest


## Inputs

1. Genomic DNA for off-target search must be deposited in a fasta format as a .txt file in the folder `complete_sequences`
2. ORFs for sgRNA design must be deposited in a fasta format as a .txt file in the folder `ORF_sequences`

## Output

1. The final library will be available as output in the folder `output`


## Execution

For efficiency purposes, few instances of the script can be run in parallel on subsets of the ORF list and then reassemble to form the final library as following:

**1. ORF list subset division.**

The list of ORF needs to be divided in multiple files. This can be accomplished using the helper script: `divide_ORFs.py`

Path to your target ORF list must be specified in the script as shown below.

	#-----------
	# CHANGE HERE FOR YOUR APPLICATION
	n=250 #ORFs_per_file 
	path_to_ORFs="ORF_sequences/Synechocystis_ORFs.txt" #Path_to_ORF_file 
	#-----------

This code will divide your ORF lists in multiple files 'ORF_sequences/ORF_N', where N is the N-th subset of the ORF list. 

**2. Designing guiding RNA.**

The script is then run on each file independently. Running the script for a given subset is executed in the terminal as shown below:
```
python library_designer_parallel N
```

N is the integer to run the script on the Nth subset of the library. If you divided your library in three subsets, you will need to run the code three times:

```
python library_designer_parallel 1
python library_designer_parallel 2
python library_designer_parallel 3
```

Output sublibraries will end up in the output folder and can then be further combined to form the final library. 

**3. Generating random control sgRNA.** ***(Optional)***

A helper script can be used to generate random control sgRNA. However, these have not been tested for potential off-target. 

KS 22.11.22






