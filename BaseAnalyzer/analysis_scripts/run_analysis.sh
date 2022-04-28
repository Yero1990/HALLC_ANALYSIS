#!/bin/bash

# shell script to run CaFe analysis

runNum=3289
daq_mode="coin"
e_arm="SHMS"
analysis_type="data"
hel_flag=0
target="LD2"
bcm_type="BCM4A"
bcm_thrs=5
trig_type="trig6"
combine_runs=0

CMD="root -l -q -b \"main_analysis.cpp(${runNum},    \\\"${daq_mode}\\\", 
                                   \\\"${e_arm}\\\", \\\"${analysis_type}\\\",
                                      ${hel_flag},   \\\"${target}\\\",
                                   \\\"${bcm_type}\\\", ${bcm_thrs},
                                   \\\"${trig_type}\\\", ${combine_runs}
                     )\""

echo $CMD
eval $CMD
