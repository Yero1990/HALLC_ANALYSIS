
## Hall C Kaon LT (2018) Analysis
This directory contains Helicity Analyzer C++ class, (`helicityAnalyzer.cpp (.h)`) derived from the more generic base Hall C Base Analyzer, `baseAnalyzer.cpp(.h)`, in order to carry out the beam-asymmetry analysis part for the Kaon LT (2018) Hall C Experiment.

**The main nuclear reaction of interest is <br>**
H(*e, e'h)*X: where the incident beam electron *e* interacts with a hydrogen target atom ( H ) and scatters ( *e'* ) into the HMS where it is detected in coincidence with the knocked-out hadron ( *h* ) in the SHMS, and the missing mass continuum, X, is reconstructed from momentum conservation laws.  

The knocked-out hadron can be selected to be either a Proton, Pion or Kaon,
depending on the analysis to be carried out. The *Missing Mass* spectrum is a continuum of multiple particle, such as: Lambda, Sigma, Eta, neutron etc. 

![Reaction Kinematics](./presentations/reaction_kin.png)


Figure 1: General coincidence reaction kinematics summarizing the detected and "Missing" particles during the KaonLT 2018. The helicity states are only relevant for the beam-spin asymmetry analysis to be done in parallel with the main KaonLT analysis.

**Data Acquisition (DAQ) Mode** <br>
During production data-taking, the DAQ was set to COIN MODE (coincidence mode) and the main coincidence trigger used was PTRIG5 = PTRIG1 X PTRIG3 ---> PTRIG5 = ( SHMS 3/4 ) X HMS EL-REAL

**Particle Identification (PID)** <br>

* HMS (electron) identification is done via a Calorimeter and Gas Cherenkov cuts
* SHMS (hadron) identification is primarily done with a coincidence time cut on the relevant hadron in question
	* **Pion / Kaon Separation**: 	Provided by the Heavy Gas Cherenkov for central momenta P >  3.4 GeV/c. The type of gas/ pressure are configured such that the detector will **NOT** emit Cherenkov radiation (signal) for a Kaon, but it will for a Pion. Usually, a cut on the number of photoelectrons is done ( hgcer_npeSum < 1.5 npe ), for example.
	* **Proton / Kaon Separation**: Provided by the Aerogel Cherenkov for central momenta P > 3 GeV/c. The Aerogel Chereknov is configured such that a signal will be emitted if the particle is a Kaon. A cut on the number of photoelectrons is done, (paero_npeSum > 1.5 npe), for example.

## How-To Guide
The `main.cpp` reads the relevant parameters from the `main_controls.inp` file and uses these parameters to initialize an instance of the `helicityAnalyzer` class.  The instantiated object then calls a method of the class: `run_helicity_analysis()`. This method calls all the other necessary methods (in a particular order) to carry out the full helicity analysis.

To successfully run the analyzer, here are some helpful hints: <br>

* check and execute the `setup.sh` shell script, as this creates the input file where the paths and name patterns of the input/output `ROOTfile` and `REPORT` files are defined <br> 
<b>NOTE:</b> 1. make sure to define the filename patterns of the files to be read in as well as the files to be created 2. make sure the symbolic link of the `inp/set_basic_filenames.inp` points to the correct filename created by the shell script.

* check the parameters set in the `main_controls.inp` file and make sure they make sense for the specific analysis to be carried out
	* pay specific attention to setting the proper input filenames <br> 
	  <b>NOTE:</b> the files in the `./inp` directory point to symbolic links which
	  should be checkeed to make sure they are pointing to the correct filenames 
	  (execute the `ls -ltrh` command to reveal the link path) 
	  
* once the input filenames are properly setup and read in the `main_controls.inp` file, run the analysis code: `root -l main.cpp`	
