import sys
sys.path.insert(1, '../../')
from src.analysis.data_analysis_hiv import *
from src.analysis.data_analysis_ecoli import *

"""
This file is used to run methods which have been written in the
data_analysis_hiv.py and the data_analysis_ecoli.py file for data collection
and sometimes graphing figures of merit of our project's use on our data.

Usage:
The only functions for replicating our data collection are commented out, line
by line, in this file. One method at a time may be uncommented and run.
"""



"""
This part contains the analytic methods for HIV strains
"""
"""
This function compares the size of each data structure (hashset, bloom filter,
and counting filter) on the same data, as k-mer size is varied.
"""
# kmerLength_vs_hashset_size()
"""
This function compares the building time of each data structure as k-mer size
is varied.
"""
# compareTimeAnalyses()
"""
This function compares the amount of time it takes to merge each method as
the number of different strains is increased.
"""
# compareStrainTimeAnalysis()
"""
This function compares the accuracy of the bloom filter (false negatives,
false positives, and correctness) to the hashset (which is 100% accurate). It
also compares the similarity of an external string to the strains that we put
into the filters
"""
#accuracyAnalysisHIV()
"""
This is the function for future researchers that tells you whether a specific
k-mer is common across many strains of HIV.
The input argument 'kmer' should be replaced with the actual k-mer.
Intersections should be replaced with the number of intersections the user wants.
"""
#mark_test_strain_kmers(kmer, intersections)


"""
This function does what the above function does for HIV but with our own test strain.
Returns a dictionary with k-mers that are common (true) and k-mers that are not common (false)
"""
#mark_test_strain_kmers_default()




"""
This part contains the analytic methods for E. coli strains (significantly
longer than HIV)
"""
"""
This function compares the time it takes to build hashsets*, bloom filters, and
counting filters as the size of k-mers are varied.

*The functionality to build hashsets is written but commented out in the
function declaration, as our computers do not have enough memory to run that
that on data sets that big.
"""
#compareTimeAnalysesEcoli()

"""
This function outputs the similarity of an external strain to the strains that
we put in the filters
"""
#accuracyAnalysisECOLI()

"""
This is the function for future researchers that tells you whether a specific
k-mer is common across many strains of E. coli.
The input argument 'kmer' should be replaced with the actual k-mer.
Intersections should be replaced with the number of intersections the user wants.
"""
#mark_test_strain_kmers_ecoli(kmer, intersections)
