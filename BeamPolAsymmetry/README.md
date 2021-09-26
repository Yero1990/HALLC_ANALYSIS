# Beam Polarization/Helicity Asymmetry

This is a project Steve Wood started to work on by analyzing beam-polarized data
from the Kaon-LT and part of SIDIS experiment (2018-2019) at Hall C. At the time,
these experiments had polarized beam, but did not really need it. However, it may be
possible to extract beam asymmetry data for 'free'. 


Link: https://redmine.jlab.org/projects/kltexp/wiki/Single_Spin_Asymmetry_Analysis


----Meeting with Steve Wood (July 31, 2020)----
1)In the beam-asymmetry term in the cross section, the sin(phi) dependence gets
integrated out when extracting the unpolarized cross section, so this term
disappears. That is, the cross section is NOT binned in sin(phi)

2)All the sin(theta) dependencies in the general cross-section expression are inside
the structure functons, and NOT explicitly written.

----ANALYSIS STEPS----
** 1st step: Try to reproduce a beam asymmetry using a
single small Kaon-LT data set to compare with Steve's result.

COMMENTS:
0) Ask R. Trotta / Tanja Horn if there are presently other students working on the Kaon-LT beam asymmetry 
1) Ask R. Trotta / Tanja Horn if they have included the helicity leafs inside the ROOTfiles
2) Ask R. Trotta about status of Kaon-LT python scripts, and consider possibility of using python for analysis
3) Beam Asymmetry uses MAID parametrization. Ask M. Jones whether this MAID or other models for Kaon-LT / SIDIS are
   used in SIMC.
4) Planned meeting with P. Bosted, T. Horn and R. Trotta on beam-spin asymmetry of Kaon-LT data