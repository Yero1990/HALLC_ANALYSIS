
# Hall C Base Analyzer

Author: C. Yero <br>
email: cyero@jlab.org | cyero002@fiu.edu

Brief: This README file gives a brief description and  overview
of the generic Hall C data analyzer.

Overview: The analyzer is based on ROOT C++ and consists of a single
class which is made up of various methods or functions with specified
tasks.


## Directory Content 

`main.cpp`:
The 'main.cpp' reads in a 'run_list.txt' file with the list of runs
to analyze, and calls an instance of the baseAnalyzer class for each run.

`baseAnalyzer.h (.cpp)`: The main analyzer is written as class called **baseAnalyzer**
and consists of a header file (.h) and a .cpp file. The header
file has the methods/function prototypes as well as most variables
used in the .cpp file.  The actual methods/functions are defined in the
.cpp file. The last method in the .cpp file, **run\_data\_analysis()** calls
all the necessary methods to carry out the analysis. In addition, the
code reads in two inputs files described below.

### Contents of the`inp/` directory:

1. `main_controls.inp ` (symbolic link):
This input file contains the main controls and
analysis cuts which the user can modify. For example, the user can set:
input/output fileName patterns, BCM type, beam current threshold cut,
which trigger type to analyze which will be used in dead-time corrections, etc.

2. `set_basic_cuts.inp` (symbolic link) :
This input file contains the main controls and
analysis cuts which the user can modify. For example, the user can set
s a variety of data analysis cuts that can be turned ON or OFF, as well
as modify the ranges. 

3. `set_basic_histos.inp` (symbolic link) :
This input file controls the histogram
binning. This is an extremely important file, as it is intended to
be used by both DATA and SIMULATION to ensure that the exact same
binning is used for both analysis types. I have added a wide variety
of histograms, some of which might not even be needed in some analysis.

4. `runlist.txt` (symbolic link):
This file
