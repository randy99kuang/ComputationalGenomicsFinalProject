# ComputationalGenomicsFinalProject

Team Members: Eric Rong, Amritpal Singh, Dikshith Kasimahanthi, Randy Kuang 

## Research Goals 
### What is the method you want to develop?
Our plan is to first find many different strains of an organism such as HIV, and split each of the organism strain’s genomes into k-mers. We will keep track of each k-mer in every strain by storing the k-mers in a bloom filter, one for each strain. After we do this for several different strains of our organism, we will start intersecting or unioning the bloom filters to only account for the most common kmers within all of the strains. We plan to figure out which formulation of intersect and union across the sets will give us a more accurate picture of the most common kmers in these strains. Once we do this to each strain’s k-mer set, we will have a single bloom filter that contains the most common kmers. 

From here, we can take a new strain of our organism and compare it to this intersection/union bloom filter. We will use this filter to identify which k-mers in this new strain are common across other strains, which will allow us to identify functionally distinct regions of the new strain. 

Finally, we plan on duplicating the methods listed above with other set implementations, such as quotient filters or counting filters.  

### How will you evaluate your method?
1. First, we want to evaluate the space usage of our various set implementations. Since data structures like bloom filters trade off accuracy for less space, the space usage will always be evaluated in comparison to a regular HashMap, which we know is error-free but has comparatively high space usage.
2. Moreover, we can measure the time usage of our various implementations versus a regular HashMap by using Python’s built-in timing features. 
3. Next, we want to analyze how the set of intersection/union k-mers grows in comparison to our input data. This can be done in several ways. For example, we can measure average edit distance of strains versus the size of their intersection/union k-mers, or the number of total strains versus the size of their intersection/union k-mers. 
4. Finally, we want to measure the accuracy of our implementations versus the error-free HashMap. Specifically, in the final intersection/union bloom or quotient filter, we can measure how many false positives k-mers exist versus k-mers that actually should be in the intersection/union.

### What input data do you need?
For this project, we need genomes of organisms that have several closely related biological variants, or strains. As a result, we are planning on starting with the genome of the virus HIV (see reference #11 below). This website includes consensus and ancestral sequences, subtype reference alignments, and complete alignments. While virus genomes are a starting point for data, they are oftentimes too small in size to reflect differences in the space/accuracy of our various implementations. As a result, we are also planning on using E. Coli genomes (reference #12 below), which also exist in multiple variants and are far longer than the genome of HIV. 

### List of at least 3 milestones: tasks to accomplish on the way to the research goal
1. Implement a functional bloom filter and quotient filter. 
2. Split the kmers from each strain into their own bloom filter.
3. Write an algorithm to find the most common k-mers among the many different strains by trying various orders of intersect and unions.
4. Evaluate our algorithm using the evaluation steps listed above. 

### List of at least 2 stretch goals: things to do if you finish milestones sooner than expected 
1. Try set implementations that include count statistics (ex. Counting Bloom Filters, Count Min-Sketch, and Counting Quotient Filters), and adjust our intersection/union algorithm to reflect the counts of specific k-mers in these data structures. For example, if one specific k-mer appears in one counting bloom filter 3 times and in another counting bloom filters 5 times, the intersection of the two will have that k-mer appear 3 times. 
2. Test our implementation on even larger data sets. For example, we can use eukaryotic genomes and manually simulate variants by adding in random edits. 
3. Create a graph for better visualization of the differences and similarities of these viral strains. Either do this initially using a modified De Brujin graph or a modified Overlap Graph. If we have enough time, find a better way to graph and visualize this data.
