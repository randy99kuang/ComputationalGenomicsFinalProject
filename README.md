#### Installation Instructions
Before running our project, a few non-default Python packages need to be installed. The instructions below summarized installation instructions.
1. The package ```bitarray``` needs to be installed in order to utilize bit arrays in the Bloom Filter
2. The package ```Pympler``` needs to be installed to run deep recursive space tests
3. The package ```Matplotlib``` needs to be installed to create plots of data
4. The package ```mmh3``` needs to be installed to get random hash functions. In our experience mmh3 still doesn't work even after installation on Windows Computers. To fix this problem we needed to download VisualStudio BuildTools, which mmh3 relies on. To download BuildTools, go to this link https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2017 and search for "Build Tools for Visual Studio 2019". Then, run the executable file that is downloaded, and when prompted, make sure to select "C++ Build Tools" in the download options. For the Windows users in our group, after downloading build tools, mmh3 ran properly. Note that Mac users in our group did not need to download build tools.

#### Project Description
- The ```src``` folder contains all code related to data structures, data analysis, and union/intersection algorithms.
- The ```Figures``` folder contains plots/graphs that we created using methods in the ```analysis``` folder. More instructions about running analysis methods are in the section below.
- The ```tests``` folder contains unit tests for our Set data structures

#### Running the Project
The only Python file that should actually be run is ```main.py```, which is found in ```src/analysis```. This file calls functions that are in ```data_analysis_ecoli.py``` and ```data_analysis_hiv.py``` to populate data structures and create figures/graphs. All methods that should be called are already contained in ```main.py``` and commented out. Thus, to run one of these tests, simply uncomment its method. Each of the functions and what they test are described in ```main.py```, but more details about their implementation can be found in the method descriptions in ```data_analysis_ecoli.py``` and ```data_analysis_hiv.py```.
For those interested, unit tests of our data structures may be run on our data structures (bloom filters, counting filters, and hashsets). These can be found in the ```tests``` folder. To run the unit tests, run this from the command line: ```python -m unittest tests/test_hashset.py```.

#### Determining if a test strains k-mer is common across many strains
One of the goals in our project was to identify what regions of a test strain are common across many different strains. The function mark_test_strain_kmers can do this for us. The function takes in an individual k-mer and the number of intersections. To correctly call and use this function, one would need to split a test strain into its k-mers and keep track of which k-mers returned common (true) and which k-mers returned not-common (false). The user would call the function like this: mark_test_strain_kmers(k-mer, intersections). The number of intersections is something the user will need to specify. They should keep in mind that as the intersections go up, the similarity of the test strain to the database of k-mers goes down. After looping through all the k-mers of a test-strain, they will be able to tell which k-mers account for regions of the genome that are likely regions of high mutations, i.e. those k-mers that return false and so are not common k-mers across strains.
