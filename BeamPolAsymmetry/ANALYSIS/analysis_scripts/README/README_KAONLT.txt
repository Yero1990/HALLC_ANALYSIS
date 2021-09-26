==============
README KAONLT
=============

This directory contains the generic Hall C analyzer code
as well as the derived helicity class to carry out the beam asymmetry
analysis for Kaon LT. The derived helicity class needs to be modified
to the first part of the Hall C Kaon LT experiment that ran on Fall 2018.

The main reactions were:

H(e, e'h)X : electron interacts with hydrogen target (H) and scatters (e')
in coincidence with a hadron (h), and the missing mass (X) is reconstructed

The hadron can be selected to be either a Proton, Pion or Kaon, depending on the
analysis. The Missing mass is a continuum of multiple particles, such a
Lambda, Sigma, Eta, etc.

The Kinematics:

DAQ: coincidence mode |  coincidence trigger: ptrig5, Ps5_factor = 1 (for the runs we currently have locally: run4868 and run4917)

electron arm: HMS                (neg. polarity)
hadron arm: SHMS (pions, kaons, protons)  (pos. polarity)


Particle Identification:

* Hadron identification in the SHMS is primarily done with a coincidence time cut.


Pi+ / K+ Separation:
----------------------
Provided by the heavy gas cherenkov for P_shms > 3.4 GeV/c. The gas cherenkov/pressure
are configured such that the cherenkov will *NOT* emit a signal for a K+ particle (Usually, a cut
on the number of photoelectrons is done, hgcer_npeSum < 1.5 npe, for example).


Proton / K+ Separation:
------------------------
Provided by the aerogel cherenkov for P_shms > 3 GeV/c. The aerogel chereknov is configured such that
a signal will be emitted if the particle is a Kaon.
(A cut on the number of photoelectrons is done, paero_npeSum > 1.5, for example)

* electron identification in the HMS is primarily done with a cut on the calorimeter normalized energy
   in combination with the HMS gas chereknov number of photo electrons cut.


