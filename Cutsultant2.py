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

#List of all the class variables I think we'll need
IUPAC = IUPACAmbiguousDNA()	
restriction_enzymes = raw_input('Please input all the restriction enzymes you have available separated by spaces. E.g., >>EcoRI SpeI BamHI: ')
plasmid_WithoutInsert = Seq(raw_input('Please enter your second nucleotide sequence that compiles with the IUPAC alphabet. This sequence should include the plasmid you wish you insert: '), IUPAC)
plasmid_WithInsert = Seq(raw_input('Please enter a nucleotide sequence that compiles with the IUPAC alphabet. This sequence should be your original chromosome (no plasmid): '), IUPAC)
smallestSize = raw_input('What is the smallest band you can measure? ')
largestSize = raw_input('What is the large band you can measure? ')
enzyme_list = restriction_enzymes.split() #Creates a list of restriction enzymes, rofl
restrictionbatch = RestrictionBatch(enzyme_list)

def enzymeEliminator(dict1, dict2, restrictionbatch):
	"""Takes in two dictionaries which map restriction enzymes to number of times it cuts and determines whether they cut the same number or not"""
	for key in dict1:
		if (dict1.has_key(key) and dict2.has_key(key)):
			if (dict1[key] == dict2[key]): 
				restrictionbatch.remove(key)

def plasmidAnalysis(restrictionbatch, plasmid_WithInsert, plasmid_WithoutInsert):
	"""Performs the analysis of everything"""

	# We will now build a restriction analysis
	plasmid_WithoutInsert_analysis = Analysis(restrictionbatch, plasmid_WithoutInsert)
	plasmid_WithInsert_analysis = Analysis(restrictionbatch, plasmid_WithInsert)

	# First, we will do a full restriction analysis
	# After reviewing the code I am not sure why we are doing this. Keeping it in for now, I'll remove it later
	full_without_analysis = plasmid_WithoutInsert_analysis.full() #Two dictionaries	
	full_with_analysis = plasmid_WithInsert_analysis.full() 

	# Next, we will eliminate any enzyme that cuts the same number of times in both plasmids
	# To do this, first, we will get a dictionary that contains all the enzymes that don't cut. If an enzyme appears in both dictionaries, we 
	# remove it from the restrictionbatch
	# We repeat the process until we reach 3 (this is arbitrary, I am not sure if there's a better way of doing this)
	
	#These dictionaries will give the enzymes mapped to their cut site, and none if they don't cut
	#zero_cut_dict_without_insert = plasmid_WithoutInsert_analysis.without_site(dct = None) #Two dictionaries with enzymes that do not cut
	#zero_cut_dict_with_insert = plasmid_WithInsert_analysis.without_site(dct = None)

	#one_cut_dict_without_insert = plasmid_WithoutInsert_analysis.with_N_sites(1, dct = None)  #Self may not be needed 
	#one_cut_dict_with_insert = plasmid_WithInsert_analysis.with_N_sites(1, dct = None) #Note: The value is the location of the cut

	# SpeI EcoRI BamHI
	# GAATTC
	# GGAATTC

	for n in range(0,30):
		temp_dict1 = plasmid_WithoutInsert_analysis.with_N_sites(n, dct = None)
		temp_dict2 = plasmid_WithInsert_analysis.with_N_sites(n, dct = None)
		enzymeEliminator(temp_dict1, temp_dict2, restrictionbatch)

	#two_cut_dict_without_insert = plasmid_WithoutInsert_analysis.with_N_sites(2, dct = None)
	#two_cut_dict_with_insert = plasmid_WithInsert_analysis.with_N_sites(2, dct = None)

	#three_cut_dict_without_insert = plasmid_WithoutInsert_analysis.with_N_sites(3, dct = None)
	#three_cut_dict_with_insert = plasmid_WithInsert_analysis.with_N_sites(3, dct = None) 

	# At this point we should have dictionaries that contain all the enzymes that have specific cuts. 
	# enzymeEliminator should then eliminate the enzymes from the restrictionbatch that cut the same number of times
	
	#enzymeEliminator(zero_cut_dict_without_insert, zero_cut_dict_with_insert, restrictionbatch) 
	#enzymeEliminator(one_cut_dict_without_insert, one_cut_dict_with_insert, restrictionbatch)
	#enzymeEliminator(two_cut_dict_without_insert, two_cut_dict_with_insert, restrictionbatch)
	#enzymeEliminator(three_cut_dict_without_insert, three_cut_dict_with_insert, restrictionbatch)


	# Then we repeat this process with with_N_sites(self, N, dct=None) (possibly with a for-loop) and keep narrowing down the RE choices
	# Then we will brute force the combinations and/or add in other factors to generate a final answer	
	# For size, when the cuts are made, assign a value to each of hte segments that is the difference between the two opposite cut sites 
	# 	This way we can get the size of the band and we can see how close they are with respect to a pre-defined tolerance value	

plasmidAnalysis(restrictionbatch, plasmid_WithInsert, plasmid_WithoutInsert)
print restrictionbatch
