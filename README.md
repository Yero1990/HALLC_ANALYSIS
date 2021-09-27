# HALL C ANALYSIS
This repository is to store a general as well as experiment-specific upper-level analysis codes. Each experiment-specific code inherits its methods from the base analyzer code. 

I will also (in the future) add a directory for genric 
low-level analysis codes (for example, codes to set reference times, time windows, detector calibrations, etc.)

## HOW-TO: Upper-Level Analysis
The generic Hall C Analyzer code in the `BaseAnalyzer`, the Kaon LT Beam-Asymmetry Analyzer in the `BeamPolAsymmetry`, and any other future analyzers that inherit from the base analyzer are structured and executed in a similar manner.