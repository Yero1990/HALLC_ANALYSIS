#!/bin/sh

#Run this script to set relevant file paths and patterns

echo "Setting Filename PATH Variables . . ."


# set path to wtite the filenames input file (make sure this coincides with
# the symlink filename path in the main_controls.inp file
ofname="./inp/KaonLT/set_basic_filenames_KaonLT_ifarm.inp"


### Check if the required directories do not exist ###
if [ ! -d ${PWD}"/ROOTfiles" ] 
then
    echo "Directory ./ROOTfiles DOES NOT exist."
    echo "You MUST create (or symbolically link) ./ROOTfiles directory. Exiting Now . . .  "
    exit 2
fi

if [ ! -d ${PWD}"/REPORT_OUTPUT" ] 
then
    echo "Directory ./REPORT_OUTPUT DOES NOT exist. "
    echo "You MUST create (or symbolically link) ./REPORT_OUTPUT directory. Exiting Now . . .  "
    exit 2
fi

if [ ! -d ${PWD}"/OUTPUT" ] 
then
    echo "Directory ./OUTPUT DOES NOT exists. Creating now . . ." 
    mkdir ${PWD}"/OUTPUT"
fi


#define file paths and filename patterns (where %d is a placeholder for the run number)
input_ROOTfilePattern=${PWD}"/ROOTfiles/coin_replay_Full_%d_-1.root"            #input analysis ROOTfile naming pattern
input_REPORTPattern=${PWD}"/REPORT_OUTPUT/replay_coin_production_%d_-1.report"      #REPORT OUTPUT (usually created by hallc_replay) naming pattersn

output_ROOTfilePattern=${PWD}"/OUTPUT/basic_histos_%d.root"                         # output ROOTfile per run results
output_ROOTfilePattern_final=${PWD}"/OUTPUT/basic_histos_final.root"                # output ROOTfile (combined runs)
output_REPORTPattern=${PWD}"/OUTPUT/basic_report.txt"                               # output REPORT (this file summarizes trk eff, live time, etc. per each run analyzed)



# Write file paths / name patterns to external input file to be read in by the code
# (the single '>' is to write output to new .txt file, the '>>' is to write a new line to an existing .txt file)
echo "input_ROOTfilePattern="$input_ROOTfilePattern > $ofname
echo "input_REPORTPattern="$input_REPORTPattern >> $ofname

echo "output_ROOTfilePattern="$output_ROOTfilePattern >> $ofname
echo "output_ROOTfilePattern_final="$output_ROOTfilePattern_final >> $ofname
echo "output_REPORTPattern="$output_REPORTPattern >> $ofname
