from Bio.Restriction import *
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

# First, user inputs RE
# 	We will then modify this string and place it in a list, and then a restrictionbatch
# Second, user inputs nucleotide sequence 1 (without insert)
# Third, user inputs nucleotide sequence 2 (with insert)
# Fourth, run an analysis on both sequences with the restrictionbatch and then we will have 2 dictionaries with cut information
# Fifth, we will start narrowing this down
#		1-3 enzymes, 1-2 is optimal
#		Fragments not too small (>300)
#		Fragments not too large (<10000)
#		Differentiation needs to be possible (spacing between the cuts)

def userInputs():
	global restriction_enzymes = raw_input('Please input all the restriction enzymes you have available separated by spaces. E.g., >>EcoRI SpoI BamHI')
	global plasmid_WithoutInsert = Seq(raw_input('Please enter your second nucleotide sequence that compiles with the IUPAC alphabet. This sequence should include the plasmid you wish you insert'), IUPAC)
	global plasmid_WithInsert = Seq(raw_input('Please enter a nucleotide sequence that compiles with the IUPAC alphabet. This sequence should be your original chromosome (no plasmid)'), IUPAC)
	global int smallestSize = raw_input('What is the smallest band you can measure?')
	global int largestSize = raw_input('What is the large band you can measure')
	global IUPAC = IUPACAmbiguousDNA()

	# Then we are going to manipulate the restriction_enzyme string and place all the enzmyes in a list, which then goes into a restrictionbatch
	
	# Then we will place these into a restrictionbatch
	global restrictionbatch = enzyme_list[1]
	for x in xrange(1,enzyme_list.length()):
		restrictionbatch.add(enzyme_list[x])

def initalAnalysis(restrictionbatch, plasmid_WithInsert, plasmid_WithoutInsert):

	# We will now build a restriction analysis

	plasmid_WithoutInsert_analysis = Analysis(restrictionbatch, plasmid_WithoutInsert, linear = false)
	plasmid_WithInsert_analysis = Analysis(restrictionbatch, plasmid_WithInsert, linear = false)

	# First, we will do a full restriction analysis
	
	full_plasmid_without_insert_analysis = plasmid_WithoutInsert_analysis.full(self, linear = false) #Two dictionaries	
	full_plasmid_with_insert_analysis = plasmid_WithInsert_analysis.full(self, linear = false) 

	# Next, we will eliminate any enzyme that cuts the same number of times in both plasmids
	# To do this, first, we will get a dictionary that contains all the enzymes that don't cut. If an enzyme appears in both dictionaries, we 
	# remove it from the restrictionbatch
	
	zero_cut_list_without_insert = full_plasmid_without_insert_analysis.without_site(self, dct = None) #Two dictionaries with enzymes that do not cut
	zero_cut_list_with_insert = full_plasmid_with_insert_analysis.without_site(self, dct = None)

	for key in zero_cut_list_with_insert:
		if (zero_cut_list_with_insert[key] == zero_cut_list_without_insert[key])
			restrictionbatch.remove(key)

	# Then we repeat this process with with_N_sites(self, N, dct=None) (possibly with a for-loop) and keep narrowing down the RE choices
	# Then we will brute force the combinations and/or add in other factors to generate a final answer	
