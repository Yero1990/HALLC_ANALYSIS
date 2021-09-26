=======================
Hall C Generic Analyzer
=======================
Author: C. Yero

email: cyero@jlab.org | cyero002@fiu.edu

Brief: This README file gives a brief description and  overview
of the generic Hall C data analyzer.

Overview: The analyzer is based on ROOT C++ and consists of a single
class which is made up of various methods or functions with specified
tasks.

The analyzer consists of:

DIRECTORY: ./
-----------
main.cpp  :
-----------
The 'main.cpp' reads in a 'run_list.txt' file with the list of runs
to analyze, and calls an instance of the baseAnalyzer class for each run.

----------------------
baseAnalyzer.h (.cpp) :
----------------------
The main analyzer is written as class called 'baseAnalyzer'
and consists of a header file (.h) and a .cpp file. The header
file has the methods/function prototypes as well as most variables
used in the .cpp file.  The actual methods/functions are defined in the
.cpp file. The last method in the .cpp file, "run_data_analysis()" calls
all the necessary methods to carry out the analysis. In addition, the
code reads in two inputs files described below.

DIRECTORY: /inp
--------------------
main_controls.inp :
--------------------
The 'main_controls.inp' input file contains the main controls and
analysis cuts which the user can modify. For example, the user can set:
input/output fileName patterns, BCM type, beam current threshold cut,
which trigger type to analyze which will be used in dead-time corrections, etc.

--------------------
set_basic_cuts.inp :
--------------------
The 'set_basic_cuts.inp' input file contains the main controls and
analysis cuts which the user can modify. For example, the user can set
s a variety of data analysis cuts that can be turned ON or OFF, as well
as modify the ranges. 

----------------------
set_basic_histos.inp :
----------------------
The 'set_basic_histos.inp' input file controls the histogram
binning. This is an extremely important file, as it is intended to
be used by both DATA and SIMULATION to ensure that the exact same
binning is used for both analysis types. I have added a wide variety
of histograms, some of which might not even be needed in some analysis.



For information on how to run the code, see HOWTO in this directory.


-------------------------------------

NOTES FOR THE FUTURE:

This generic analyzer will be the starting point to
develop other analysis specific codes such as the beam asymmetry
analysis.

My intention is to keep adding additional methods specific to the
analysis. I also need to start adding SIMC methods, already created
in the DEUTERON analysis. 
