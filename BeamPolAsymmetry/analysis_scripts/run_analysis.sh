#!/bin/bash
 
# shell script to execute analysis code

CMD1="cd /u/group/c-kaonlt/USERS/cyero/HALLC_ANALYSIS/BeamPolAsymmetry/analysis_scripts"
CMD2="root -l -q -b \"main.cpp\""  

eval ${CMD1}

#source /site/12gev_phys/softenv.sh 2.3
    
eval ${CMD2}

