# HALL C ANALYSIS
This repository stores a general as well as experiment-specific upper-level analysis codes. Each experiment-specific code inherits its methods from the base analyzer code. 

I will also (in the future) add a directory for generic 
low-level analysis codes (for example, codes to set reference times, time windows, detector calibrations, etc.)

## Upper-Level Analysis Overview
Upper-Level Analysis refers to data-analysis post reference-time, detector time windows, and detector calibrations as well as any spectrometer offset determinations / optics optimizations using H(e,e'p) elastics or Carbon data.

The generic Hall C Analyzer code in the `BaseAnalyzer/analysis_scripts/`directory, the Kaon LT Beam-Asymmetry Analyzer in the `BeamPolAsymmetry/analysis_scripts/`directory, and any other future analyzers that inherit from the base analyzer are structured in a similar manner. 

The relevant scripts/directories within the `analysis_scripts/` necessary to run the analysis are: <br>

* `setup.sh` : shell script that checks if the required directories have been created (`ROOTfiles/`, `REPORT_OUTPUT/`, `OUTPUT/`) as well creates an input file (set\_basic\_filenames.inp) with the paths and file name patterns.

*  `main.cpp` : steering script that needs to be executed to run the analysis

*  `XXXXAnalyzer.cpp (.h)` : C++ class/header files where the relevant methods and functions to carry out the analysis are defined.

* `main_controls.inp` : main controls input file where parameters ara initialized, and additional input file paths are defined. These additional files are located in the `inp/` directory and look something like: <br>
	* `inp/set_basic_filenames.inp` 
	* `inp/set_basic_cuts.inp`
	* `inp/set_basic_histos.inp`
	* `inp/runlist.txt`	

These input file are actually symbolically linked to the directory relevant to the study being done. This allows for a more flexible and easier way to carry out multiple studies requiring multiple filename structures, cuts, histogram ranges, etc. while at the same time maintaining a generic filename structure via the symbolic link.

## Upper-Level Analysis HOW-TOs:

The `main.cpp` code reads the `main_controls.inp` input file containing the initialization parameters/paths and then creates an instance of the analysis .cpp class with the initialization parameters as input. For every run in the `inp/runlist.txt`, the instantiated object calls an analysis method like `run_data_analysis()` which calls all the methods used in the analysis. The `ROOTfiles/` and `REPORT_OUTPUT/` must contain the relevant ROOTfiles and report files corresponding to the runs in the `inp/runlist.txt`, as these will be read by during the analysis. For each run analyzed, the output ROOTfile and report file will be saved to `OUTPUT/` directory. In the case of multiple runs, there is an option in the `main_controls.inp` to combine the runs (simply adding histograms using the TH1F::Add() method)

Assuming all the necessary input files are in place, then from relevant `analysis_scripts/` directory do: <br>

`root -l main.cpp`

**NOTE**: Usually in experiment, there are multiplle kinematics measured, and for each kinematic, multiple runs are taken to gain the desired statistical precision. Currently, our analysis code is structured to take multiple runs of the same kinematics. So if  multiple kinematics settings are being analyzed, it is suggested that the user makes additional directories with the kinematics identifier (KIN-I, KIN-II, etc.), to store these files, and keep the directorty structure organized and clean. There are no specific guides or direction as to how to manage multiple kinematics files, so the user must get creative.

## To-Do List
Currently, analysis class is ONLY set up to analyze coincidence mode DAQ
and only experimental data. The code needs to be updated to analyze SIMC
as well.
