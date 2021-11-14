#include "baseAnalyzer.h"
#include "helicityAnalyzer.h"
#include <iostream>
#include <stdio.h>
using namespace std;

//_______________________________________________________________________________
helicityAnalyzer::helicityAnalyzer(int irun, string mode, string earm, string ana_type, Bool_t hel_flag=0) :
  baseAnalyzer(irun, mode, earm, ana_type, hel_flag)  //initialize helicity constructor with base class (min. requirements, we can add more variables later)

{
  cout << "Calling Derived Constructor  . . ." << endl;
  
  //NOTE: The base constructor initialize the general histograms as well as BCM related variables
  // among others. The base constructor is automatically called before the derived (this) constructor

  //Below, we can add additional variables we may want to initialize that are specific to the helicity analysis
  //-----Initialize Helicity (Derived) Histogram Pointers-----


  //1D Coincidence Times: Reals/Background
  H_eK_ctime_real  = NULL;
  H_ePi_ctime_real = NULL;
  H_ep_ctime_real  = NULL;
  
  H_eK_ctime_rand  = NULL;
  H_ePi_ctime_rand = NULL;
  H_ep_ctime_rand  = NULL;
  
  //1D Missing Mass Plots
  H_MM_K_PID       = NULL;
  H_MM_Pi_PID      = NULL;
  H_MM_P_PID       = NULL;
  
  H_MM_K_real  = NULL;
  H_MM_Pi_real = NULL;
  H_MM_P_real  = NULL;

  //1D Missing Mass Plots for random coincidences selection
  H_MM_K_rand  = NULL;
  H_MM_Pi_rand = NULL;
  H_MM_P_rand  = NULL;

  //1D Missing Mass Plots after randoms subtraction
  H_MM_K_rand_sub  = NULL;
  H_MM_Pi_rand_sub = NULL;
  H_MM_P_rand_sub  = NULL;

  // Particle Beta (no PID cuts)
  H_Beta_noPID   = NULL;

  // Particle Beta (PID cuts)
  H_Beta_K_PID     = NULL;
  H_Beta_Pi_PID    = NULL;
  H_Beta_P_PID     = NULL;

  // Particle Beta (real coincidence selection)
  H_Beta_K_real     = NULL;
  H_Beta_Pi_real    = NULL;
  H_Beta_P_real     = NULL;

  // Particle Beta (random coincidence selection)
  H_Beta_K_rand     = NULL;
  H_Beta_Pi_rand    = NULL;
  H_Beta_P_rand     = NULL;

  // Particle Beta (random coincidence subtraction)
  H_Beta_K_rand_sub     = NULL;
  H_Beta_Pi_rand_sub    = NULL;
  H_Beta_P_rand_sub     = NULL;
  
  //2D Histograms

  // Beta vs. Coincidence Time (no PID selection cuts)
  H_Beta_vs_ctime_K_noPID  =  NULL;
  H_Beta_vs_ctime_Pi_noPID =  NULL;
  H_Beta_vs_ctime_P_noPID  =  NULL;

  // Beta vs. Coincidence Time (PID selection cuts)
  H_Beta_vs_ctime_K_PID    =  NULL;
  H_Beta_vs_ctime_Pi_PID   =  NULL;
  H_Beta_vs_ctime_P_PID    =  NULL;


  //Kaon LT Kinematic Histos
  H_Q2_vs_W   = NULL;
  H_t_vs_phxq = NULL;   

  //NOTE: The advantage of using the 2D histos below is that one can extract asymmetry information from ph_xq for each th_xq_cm bin
  //easily via a loop over each of the 2D bins
  
  //2D: In-plane angles of X about q, th_xq_cm (in CM frame) "VS" Out-of-plane angles of X about q, ph_xq  (Assume detected X are Kaons)
  H_thxqCM_vs_phxq_K_pos = NULL;    //+ helicity
  H_thxqCM_vs_phxq_K_neg = NULL;   //- helicity

  //Same as above, but for randoms (events out side the main coin. peak)
  H_thxqCM_vs_phxq_K_pos_rand = NULL;  //+ helicity
  H_thxqCM_vs_phxq_K_neg_rand = NULL;  //- helicity

  //Same as above, but for randoms subtracted (Will be used to extract beam asymmetry)
  H_thxqCM_vs_phxq_K_pos_rand_sub = NULL;  //+ helicity
  H_thxqCM_vs_phxq_K_neg_rand_sub = NULL;  //- helicity

  
  //2D: In-plane angles of X about q, th_xq_cm (in CM frame) "VS" Out-of-plane angles of X about q, ph_xq  (Assume detected X are Pions)
  H_thxqCM_vs_phxq_Pi_pos = NULL;  //+ helicity
  H_thxqCM_vs_phxq_Pi_neg = NULL;  //- helicity

  //Same as above, but for randoms (events outside the main coin. peak)
  H_thxqCM_vs_phxq_Pi_pos_rand = NULL;  //+ helicity
  H_thxqCM_vs_phxq_Pi_neg_rand = NULL;  //- helicity

  //Same as above, but for randoms subtracted (Will be used to extract beam asymmetry)
  H_thxqCM_vs_phxq_Pi_pos_rand_sub = NULL;  //+ helicity
  H_thxqCM_vs_phxq_Pi_neg_rand_sub = NULL;  //- helicity

  //2D: In-plane angles of X about q, th_xq_cm (in CM frame) "VS" Out-of-plane angles of X about q, ph_xq  (Assume detected X are Protons)
  H_thxqCM_vs_phxq_P_pos = NULL;  //+ helicity
  H_thxqCM_vs_phxq_P_neg = NULL;  //- helicity

  //Same as above, but for randoms (events outside the main coin. peak)
  H_thxqCM_vs_phxq_P_pos_rand = NULL;  //+ helicity
  H_thxqCM_vs_phxq_P_neg_rand = NULL;  //- helicity

  //Same as above, but for randoms subtracted (Will be used to extract beam asymmetry)
  H_thxqCM_vs_phxq_P_pos_rand_sub = NULL;  //+ helicity
  H_thxqCM_vs_phxq_P_neg_rand_sub = NULL;  //- helicity

}

//_______________________________________________________________________________
helicityAnalyzer::~helicityAnalyzer()
{
  cout << "Calling Derived Destructor  . . ." << endl;

  //The Base Destructor is automatically called LAST, so it deletes and frees the memory
  //of the pointers initialized in the base constructor

  //The derived destructor (this method) is called before the base, and here were delete the
  //variables that were initialized in the derived consstructor.
  
   //-----Delete Helicity (Derived) Histogram Pointers-----

  //1D Coincidence Times: Reals/Background                                             
  delete H_eK_ctime_real;   H_eK_ctime_real  = NULL;
  delete H_ePi_ctime_real;  H_ePi_ctime_real = NULL;
  delete H_ep_ctime_real;   H_ep_ctime_real  = NULL;
  
  delete H_eK_ctime_rand;   H_eK_ctime_rand  = NULL;
  delete H_ePi_ctime_rand;  H_ePi_ctime_rand = NULL;
  delete H_ep_ctime_rand;   H_ep_ctime_rand  = NULL;

  //1D Missing Mass Plots (only PID cuts, not coin cut)
  delete H_MM_K_PID;        H_MM_K_PID       = NULL;  
  delete H_MM_Pi_PID;       H_MM_Pi_PID      = NULL;
  delete H_MM_P_PID;        H_MM_P_PID       = NULL;

  //1D Missing Mass Plots for real coincidences
  delete H_MM_K_real;  H_MM_K_real  = NULL;  
  delete H_MM_Pi_real; H_MM_Pi_real = NULL;
  delete H_MM_P_real;  H_MM_P_real  = NULL;

  //1D Missing Mass Plots for random coincidences selection
  delete H_MM_K_rand;  H_MM_K_rand  = NULL;
  delete H_MM_Pi_rand; H_MM_Pi_rand = NULL;
  delete H_MM_P_rand;  H_MM_P_rand  = NULL;

  //1D Missing Mass Plots after randoms subtraction
  delete H_MM_K_rand_sub;  H_MM_K_rand_sub  = NULL;
  delete H_MM_Pi_rand_sub; H_MM_Pi_rand_sub = NULL;
  delete H_MM_P_rand_sub;  H_MM_P_rand_sub  = NULL;

  // Particle Beta (no PID cuts)
  delete H_Beta_noPID;   H_Beta_noPID   = NULL;

  // Particle Beta (PID cuts)
  delete H_Beta_K_PID;   H_Beta_K_PID     = NULL;
  delete H_Beta_Pi_PID;  H_Beta_Pi_PID    = NULL;
  delete H_Beta_P_PID;   H_Beta_P_PID     = NULL;

  // Particle Beta (real coincidence selection)
  delete H_Beta_K_real;   H_Beta_K_real     = NULL;
  delete H_Beta_Pi_real;  H_Beta_Pi_real    = NULL;
  delete H_Beta_P_real;   H_Beta_P_real     = NULL;

  // Particle Beta (random coincidence selection)
  delete H_Beta_K_rand;   H_Beta_K_rand     = NULL;
  delete H_Beta_Pi_rand;  H_Beta_Pi_rand    = NULL;
  delete H_Beta_P_rand;   H_Beta_P_rand     = NULL;

  // Particle Beta (random coincidence subtracted)
  delete H_Beta_K_rand_sub;   H_Beta_K_rand_sub     = NULL;
  delete H_Beta_Pi_rand_sub;  H_Beta_Pi_rand_sub    = NULL;
  delete H_Beta_P_rand_sub;   H_Beta_P_rand_sub     = NULL;
  
  //2D Histograms

  // Beta vs. Coincidence Time (no PID selection cuts)
  delete H_Beta_vs_ctime_K_noPID;   H_Beta_vs_ctime_K_noPID  =  NULL;
  delete H_Beta_vs_ctime_Pi_noPID;  H_Beta_vs_ctime_Pi_noPID =  NULL;
  delete H_Beta_vs_ctime_P_noPID;   H_Beta_vs_ctime_P_noPID  =  NULL;

  // Beta vs. Coincidence Time (PID selection cuts)
  delete H_Beta_vs_ctime_K_PID;   H_Beta_vs_ctime_K_PID    =  NULL;
  delete H_Beta_vs_ctime_Pi_PID;  H_Beta_vs_ctime_Pi_PID   =  NULL;
  delete H_Beta_vs_ctime_P_PID;   H_Beta_vs_ctime_P_PID    =  NULL;
  
  delete H_Q2_vs_W;   H_Q2_vs_W   = NULL;
  delete H_t_vs_phxq; H_t_vs_phxq = NULL;

  //2D: In-plane angles of X about q, th_xq_cm (in CM frame) "VS" Out-of-plane angles of X about q, ph_xq  (Assume detected X are Kaons)
  delete H_thxqCM_vs_phxq_K_pos; H_thxqCM_vs_phxq_K_pos = NULL;    //+ helicity
  delete H_thxqCM_vs_phxq_K_neg; H_thxqCM_vs_phxq_K_neg = NULL;   //- helicity

  //Same as above, but for randoms (events out side the main coin. peak)
  delete H_thxqCM_vs_phxq_K_pos_rand; H_thxqCM_vs_phxq_K_pos_rand = NULL;  //+ helicity
  delete H_thxqCM_vs_phxq_K_neg_rand; H_thxqCM_vs_phxq_K_neg_rand = NULL;  //- helicity

  //Same as above, but for randoms subtracted
  delete H_thxqCM_vs_phxq_K_pos_rand_sub; H_thxqCM_vs_phxq_K_pos_rand_sub = NULL;  //+ helicity
  delete H_thxqCM_vs_phxq_K_neg_rand_sub; H_thxqCM_vs_phxq_K_neg_rand_sub = NULL;  //- helicity

  
  //2D: In-plane angles of X about q, th_xq_cm (in CM frame) "VS" Out-of-plane angles of X about q, ph_xq (Assume detected X are Pions)
  delete H_thxqCM_vs_phxq_Pi_pos; H_thxqCM_vs_phxq_Pi_pos = NULL;  //+ helicity
  delete H_thxqCM_vs_phxq_Pi_neg; H_thxqCM_vs_phxq_Pi_neg = NULL;  //- helicity

  //Same as above, but for randoms (events outside the main coin. peak)
  delete H_thxqCM_vs_phxq_Pi_pos_rand; H_thxqCM_vs_phxq_Pi_pos_rand = NULL;  //+ helicity
  delete H_thxqCM_vs_phxq_Pi_neg_rand; H_thxqCM_vs_phxq_Pi_neg_rand = NULL;  //- helicity

  //Same as above, but for randoms subtracted
  delete H_thxqCM_vs_phxq_Pi_pos_rand_sub; H_thxqCM_vs_phxq_Pi_pos_rand_sub = NULL;  //+ helicity
  delete H_thxqCM_vs_phxq_Pi_neg_rand_sub; H_thxqCM_vs_phxq_Pi_neg_rand_sub = NULL;  //- helicity
  

  //2D: In-plane angles of X about q, th_xq_cm (in CM frame) "VS" Out-of-plane angles of X about q, ph_xq (Assume detected X are Protons)
  delete H_thxqCM_vs_phxq_P_pos; H_thxqCM_vs_phxq_P_pos = NULL;  //+ helicity
  delete H_thxqCM_vs_phxq_P_neg; H_thxqCM_vs_phxq_P_neg = NULL;  //- helicity

  //Same as above, but for randoms (events outside the main coin. peak)
  delete H_thxqCM_vs_phxq_P_pos_rand; H_thxqCM_vs_phxq_P_pos_rand = NULL;  //+ helicity
  delete H_thxqCM_vs_phxq_P_neg_rand; H_thxqCM_vs_phxq_P_neg_rand = NULL;  //- helicity

  //Same as above, but for randoms subtracted
  delete H_thxqCM_vs_phxq_P_pos_rand_sub; H_thxqCM_vs_phxq_P_pos_rand_sub = NULL;  //+ helicity
  delete H_thxqCM_vs_phxq_P_neg_rand_sub; H_thxqCM_vs_phxq_P_neg_rand_sub = NULL;  //- helicity
  
}


//_______________________________________________________________________________
void helicityAnalyzer::SetFileNames()
{
  cout << "Calling Derived SetFileNames() " << endl;

  //As of now, there is nothing to add to the existing methods from baseAnalyzer.
  
  //Call the existing method, as is is needed to set file names (filenames are set in main_controls.inp)
  baseAnalyzer::SetFileNames();

  //HERE: Read in additional filenames if needed
  
}

void helicityAnalyzer::ReadInputFile(string ftype="")
{
  cout << "Calling Derived ReadInputFile() " << endl;
  
  //Call the existing method, as is is needed to read in parameters from
  //main controls 
  if(ftype=="main_controls")
    {
      cout << "reading main_controls.inp . . . " << endl;
      baseAnalyzer::ReadInputFile("main_controls");
    }

  //Call the existing method, as is is needed to read in parameters for
  //tracking efficiency cuts input file
  if(ftype=="trk_eff_cuts")
    {
      cout << "reading trk_eff_cuts in set_basic_cuts.inp . . . " << endl;  
      baseAnalyzer::ReadInputFile("trk_eff_cuts");
    }

  //HERE: Read in data analysis cuts (RE-DEFINE FORMAT FOR ANALYSIS CUTS) ! ! !
  //NOTE: Since the trk eff cuts and data analysis cuts are in the same cuts file,
  //leave the trk eff. cuts alone, and modify ONLY the format for data analysis cuts
  if(ftype=="analysis_cuts")
    {

      //Read data analysis parameters here (in accordance with the analysis requirements)
      //Example of how to read in parameters is in the baseAnalyzer::ReadInputFile()
      
      cout << "reading analysis_cuts in set_basic_cuts.inp . . . " << endl;    

      //Read IN acceptance cuts
      edelta_cut_flag = stoi(split(FindString("edelta_cut_flag", input_CutFileName.Data())[0], '=')[1]);
      edelta_min = stod(split(FindString("edelta_min", input_CutFileName.Data())[0], '=')[1]);
      edelta_max = stod(split(FindString("edelta_max", input_CutFileName.Data())[0], '=')[1]);

      hdelta_cut_flag = stoi(split(FindString("hdelta_cut_flag", input_CutFileName.Data())[0], '=')[1]);
      hdelta_min = stod(split(FindString("hdelta_min", input_CutFileName.Data())[0], '=')[1]);
      hdelta_max = stod(split(FindString("hdelta_max", input_CutFileName.Data())[0], '=')[1]);

      // Z-Reaction Vertex Difference Cut
      ztarDiff_cut_flag = stoi(split(FindString("ztarDiff_cut_flag", input_CutFileName.Data())[0], '=')[1]);
      ztarDiff_min = stod(split(FindString("ztarDiff_min", input_CutFileName.Data())[0], '=')[1]);
      ztarDiff_max = stod(split(FindString("ztarDiff_max", input_CutFileName.Data())[0], '=')[1]);
    
      //Read IN PID cuts

      //-electron-
      hcer_pidCut_flag = stoi(split(FindString("hcer_pidCut_flag", input_CutFileName.Data())[0], '=')[1]);
      elec_hcer_npe_thrs = stod(split(FindString("elec_hcer_npe_thrs", input_CutFileName.Data())[0], '=')[1]);
      
      hetot_trkNorm_pidCut_flag = stoi(split(FindString("hetot_trkNorm_pidCut_flag", input_CutFileName.Data())[0], '=')[1]);
      elec_hcal_thrs     = stod(split(FindString("elec_hcal_thrs", input_CutFileName.Data())[0], '=')[1]);

      //-Kaon-
      K_paero_npe_flag  = stoi(split(FindString("K_paero_npe_flag", input_CutFileName.Data())[0], '=')[1]);
      K_paero_npe_thrs  = stod(split(FindString("K_paero_npe_thrs", input_CutFileName.Data())[0], '=')[1]);

      K_phgcer_npe_flag  = stoi(split(FindString("K_phgcer_npe_flag", input_CutFileName.Data())[0], '=')[1]);
      K_phgcer_npe_thrs = stod(split(FindString("K_phgcer_npe_thrs", input_CutFileName.Data())[0], '=')[1]);

      K_beta_flag       = stoi(split(FindString("K_beta_flag", input_CutFileName.Data())[0], '=')[1]);
      K_beta_min        = stod(split(FindString("K_beta_min", input_CutFileName.Data())[0], '=')[1]);
      K_beta_max        = stod(split(FindString("K_beta_max", input_CutFileName.Data())[0], '=')[1]);

      eK_ctime_flag     = stoi(split(FindString("eK_ctime_flag", input_CutFileName.Data())[0], '=')[1]);
      eK_ctime_thrs     = stod(split(FindString("eK_ctime_thrs", input_CutFileName.Data())[0], '=')[1]); 

      //-Pion-
      Pi_phgcer_npe_flag  = stoi(split(FindString("Pi_phgcer_npe_flag", input_CutFileName.Data())[0], '=')[1]);
      Pi_phgcer_npe_thrs = stod(split(FindString("Pi_phgcer_npe_thrs", input_CutFileName.Data())[0], '=')[1]);

      Pi_beta_flag       = stoi(split(FindString("Pi_beta_flag", input_CutFileName.Data())[0], '=')[1]);
      Pi_beta_min        = stod(split(FindString("Pi_beta_min", input_CutFileName.Data())[0], '=')[1]);
      Pi_beta_max        = stod(split(FindString("Pi_beta_max", input_CutFileName.Data())[0], '=')[1]);

      ePi_ctime_flag     = stoi(split(FindString("ePi_ctime_flag", input_CutFileName.Data())[0], '=')[1]);
      ePi_ctime_thrs     = stod(split(FindString("ePi_ctime_thrs", input_CutFileName.Data())[0], '=')[1]); 

      //-Proton-
      P_paero_npe_flag  = stoi(split(FindString("P_paero_npe_flag", input_CutFileName.Data())[0], '=')[1]);
      P_paero_npe_thrs  = stod(split(FindString("P_paero_npe_thrs", input_CutFileName.Data())[0], '=')[1]);

      P_phgcer_npe_flag  = stoi(split(FindString("P_phgcer_npe_flag", input_CutFileName.Data())[0], '=')[1]);
      P_phgcer_npe_thrs = stod(split(FindString("P_phgcer_npe_thrs", input_CutFileName.Data())[0], '=')[1]);

      P_beta_flag       = stoi(split(FindString("P_beta_flag", input_CutFileName.Data())[0], '=')[1]);
      P_beta_min        = stod(split(FindString("P_beta_min", input_CutFileName.Data())[0], '=')[1]);
      P_beta_max        = stod(split(FindString("P_beta_max", input_CutFileName.Data())[0], '=')[1]);

      eP_ctime_flag     = stoi(split(FindString("eP_ctime_flag", input_CutFileName.Data())[0], '=')[1]);
      eP_ctime_thrs     = stod(split(FindString("eP_ctime_thrs", input_CutFileName.Data())[0], '=')[1]); 

      
      //----Read IN Kinematic cuts----

      //-Kaon Missing Mass-
      MM_K_cut_flag = stoi(split(FindString("MM_K_cut_flag", input_CutFileName.Data())[0], '=')[1]);
      MM_K_min= stod(split(FindString("MM_K_min", input_CutFileName.Data())[0], '=')[1]);
      MM_K_max= stod(split(FindString("MM_K_max", input_CutFileName.Data())[0], '=')[1]);

      //-Pion Missing Mass-
      MM_Pi_cut_flag = stoi(split(FindString("MM_Pi_cut_flag", input_CutFileName.Data())[0], '=')[1]);
      MM_Pi_min= stod(split(FindString("MM_Pi_min", input_CutFileName.Data())[0], '=')[1]);
      MM_Pi_max= stod(split(FindString("MM_Pi_max", input_CutFileName.Data())[0], '=')[1]);

      //-Proton Missing Mass-
      MM_P_cut_flag = stoi(split(FindString("MM_P_cut_flag", input_CutFileName.Data())[0], '=')[1]);
      MM_P_min= stod(split(FindString("MM_P_min", input_CutFileName.Data())[0], '=')[1]);
      MM_P_max= stod(split(FindString("MM_P_max", input_CutFileName.Data())[0], '=')[1]);

      
      //Read IN  offset parameters to center coincidence time around 0 ns
      K_ctime_offset = stod(split(FindString("K_ctime_offset", input_CutFileName.Data())[0], '=')[1]);
      Pi_ctime_offset = stod(split(FindString("Pi_ctime_offset", input_CutFileName.Data())[0], '=')[1]);
      P_ctime_offset = stod(split(FindString("P_ctime_offset", input_CutFileName.Data())[0], '=')[1]);
      
      //Read IN "NEW" parameters to set the upper limit of coin. time to be a multiple of the minimum thrs cut 
      eK_mult = stod(split(FindString("eK_mult", input_CutFileName.Data())[0], '=')[1]);
      ePi_mult = stod(split(FindString("ePi_mult", input_CutFileName.Data())[0], '=')[1]);
      eP_mult = stod(split(FindString("eP_mult", input_CutFileName.Data())[0], '=')[1]);
      
    }
  
}


//_______________________________________________________________________________
void helicityAnalyzer::SetCuts()
{
  cout << "Calling Derived SetCuts() " << endl;

  //As of now, there is nothing to add to the existing method from baseAnalyzer.

  //Call the ReadInputFile, as it is needed to set analysis cuts
  
  ReadInputFile("main_controls");
  ReadInputFile("trk_eff_cuts");
  ReadInputFile("analysis_cuts");

}

//_______________________________________________________________________________
void helicityAnalyzer::ReadReport()
{
  cout << "Calling Derived ReadReport() " << endl;

  //As of now, there is nothing to add to the existing method from baseAnalyzer.

  //Call the existing method, as it is needed to read report file
  baseAnalyzer::ReadReport();

  //HERE: Read in additional variables to read from report if needed
  
  
}

//_______________________________________________________________________________
void helicityAnalyzer::SetHistBins()
{
  cout << "Calling Derived SetHistBins() " << endl;

  //As of now, there is nothing to add to the existing method from baseAnalyzer.

  //Call the existing method, as it is needed to set histogram binning range
  baseAnalyzer::SetHistBins();

  //HERE: Read in additional histogram binning if needed
  
  
}


//_______________________________________________________________________________
void helicityAnalyzer::CreateHist()
{
  cout << "Calling Derived CreateHist() " << endl;

  //Call the existing method, as it is needed to create generic histograms and add them to lists
  baseAnalyzer::CreateHist();


  
  //Create TLists to store helicity histograms
  hel_HList = new TList();
  
  //HERE: Add additional histograms for helicity analysis if needed

  //----------------------------------------------------------------------------
  //---------HISTOGRAM CATEGORY: Kaon LT (2018) Helicity Analysis Histograms----
  //----------------------------------------------------------------------------

  //1D Coincidence time histos (after selecting the main coincidence time peak, which contains reals+bkg)  -- redefined histos from baseAnalyzer
  H_eK_ctime_real       = new TH1F("H_eK_ctime_real", "eK Coincidence Time (real); eK Coincidence Time [ns]; Counts / mC", coin_nbins, coin_xmin, coin_xmax);
  H_ePi_ctime_real      = new TH1F("H_ePi_ctime_real", "e#pi Coincidence Time (real) ; e#pi Coincidence Time [ns]; Counts / mC", coin_nbins, coin_xmin, coin_xmax);
  H_ep_ctime_real       = new TH1F("H_ep_ctime_real", "ep Coincidence Time (real) ; ep Coincidence Time [ns]; Counts / mC", coin_nbins, coin_xmin, coin_xmax);
  //1D Coincidence time histos (after selecting randoms outside the main coin. peak)
  H_eK_ctime_rand       = new TH1F("H_eK_ctime_rand", "eK Coincidence Time (bkg); eK Coincidence Time [ns]; Counts / mC", coin_nbins, coin_xmin, coin_xmax);
  H_ePi_ctime_rand      = new TH1F("H_ePi_ctime_rand", "e#pi Coincidence Time (bkg); e#pi Coincidence Time [ns]; Counts / mC", coin_nbins, coin_xmin, coin_xmax);
  H_ep_ctime_rand       = new TH1F("H_ep_ctime_rand", "ep Coincidence Time (bkg); ep Coincidence Time [ns]; Counts / mC", coin_nbins, coin_xmin, coin_xmax);

  //1D Missing Mass Histograms (after selecting PID cuts, but not yet coincidenece cut)
  H_MM_K_PID    = new TH1F("H_MM_K_PID",  "Kaon Missing Mass   (PID); M_{miss, K^{+}} [GeV/c^{2}]; Counts / mC ", MM_nbins, MM_xmin, MM_xmax);
  H_MM_Pi_PID   = new TH1F("H_MM_Pi_PID", "Pion Missing Mass   (PID); M_{miss, #pi^{+}} [GeV/c^{2}]; Counts / mC", MM_nbins, MM_xmin, MM_xmax);
  H_MM_P_PID    = new TH1F("H_MM_P_PID",  "Proton Missing Mass (PID); M_{miss, p} [GeV/c^{2}]; Counts / mC ", MM_nbins, MM_xmin, MM_xmax);
  //1D Missing Mass Histograms (after selecting real coincidences, within the main peak)
  H_MM_K_real  = new TH1F("H_MM_K_real", "Kaon Missing Mass (real+bkg); M_{miss, K^{+}} [GeV/c^{2}]; Counts / mC ", MM_nbins, MM_xmin, MM_xmax);
  H_MM_Pi_real = new TH1F("H_MM_Pi_real", "Pion Missing Mass (real+bkg); M_{miss, #pi^{+}} [GeV/c^{2}]; Counts / mC", MM_nbins, MM_xmin, MM_xmax);
  H_MM_P_real  = new TH1F("H_MM_P_real", "Proton Missing Mass (real+bkg); M_{miss, p} [GeV/c^{2}]; Counts / mC ", MM_nbins, MM_xmin, MM_xmax);
  //Missing Mass (after selecting random coincidences, outside the main peak)
  H_MM_K_rand  = new TH1F("H_MM_K_rand", "Kaon Missing Mass (bkg); M_{miss, K^{+}} [GeV/c^{2}]; Counts / mC ", MM_nbins, MM_xmin, MM_xmax);
  H_MM_Pi_rand = new TH1F("H_MM_Pi_rand", "Pion Missing Mass (bkg); M_{miss, #pi^{+}} [GeV/c^{2}]; Counts / mC", MM_nbins, MM_xmin, MM_xmax);
  H_MM_P_rand  = new TH1F("H_MM_P_rand", "Proton Missing Mass (bkg); M_{miss, p} [GeV/c^{2}]; Counts / mC ", MM_nbins, MM_xmin, MM_xmax);
  //Missing Mass (after subtracting the Missing Mass with random coincidences)
  H_MM_K_rand_sub  = new TH1F("H_MM_K_rand_sub", "Kaon Missing Mass (real); M_{miss, K^{+}} [GeV/c^{2}]; Counts / mC", MM_nbins, MM_xmin, MM_xmax);
  H_MM_Pi_rand_sub = new TH1F("H_MM_Pi_rand_sub", "Pion Missing Mass (real); M_{miss, #pi^{+}} [GeV/c^{2}]; Counts / mC", MM_nbins, MM_xmin, MM_xmax);
  H_MM_P_rand_sub  = new TH1F("H_MM_P_rand_sub", "Proton Missing Mass (real); M_{miss, p} [GeV/c^{2}]; Counts / mC ", MM_nbins, MM_xmin, MM_xmax);
  
  // Particle Beta (no pid selection)
  H_Beta_noPID   = new TH1F("H_Beta_noPID",  "#beta (no PID); #beta; Counts / mC", hbeta_nbins, hbeta_xmin, hbeta_xmax);  

  // Particle Beta (pid selection)
  H_Beta_K_PID    = new TH1F("H_Beta_K_PID",  "Kaon #beta (PID); #beta; Counts / mC", hbeta_nbins, hbeta_xmin, hbeta_xmax);  
  H_Beta_Pi_PID   = new TH1F("H_Beta_Pi_PID", "Pion #beta (PID); #beta; Counts / mC", hbeta_nbins, hbeta_xmin, hbeta_xmax);  
  H_Beta_P_PID    = new TH1F("H_Beta_P_PID",  "Proton #beta (PID); #beta; Counts / mC", hbeta_nbins, hbeta_xmin, hbeta_xmax);

  // Particle Beta (real coincidence selection)
  H_Beta_K_real    = new TH1F("H_Beta_K_real",  "Kaon #beta (real+bkg); #beta; Counts / mC", hbeta_nbins, hbeta_xmin, hbeta_xmax);  
  H_Beta_Pi_real   = new TH1F("H_Beta_Pi_real", "Pion #beta (real+bkg); #beta; Counts / mC", hbeta_nbins, hbeta_xmin, hbeta_xmax);  
  H_Beta_P_real    = new TH1F("H_Beta_P_real",  "Proton #beta (real+bkg); #beta; Counts / mC", hbeta_nbins, hbeta_xmin, hbeta_xmax);

  // Particle Beta (random coincidence selection)
  H_Beta_K_rand    = new TH1F("H_Beta_K_rand",  "Kaon #beta (bkg); #beta; Counts / mC", hbeta_nbins, hbeta_xmin, hbeta_xmax);  
  H_Beta_Pi_rand   = new TH1F("H_Beta_Pi_rand", "Pion #beta (bkg); #beta; Counts / mC", hbeta_nbins, hbeta_xmin, hbeta_xmax);  
  H_Beta_P_rand    = new TH1F("H_Beta_P_rand",  "Proton #beta (bkg); #beta; Counts / mC", hbeta_nbins, hbeta_xmin, hbeta_xmax);

  // Particle Beta (random coincidence subtracted)
  H_Beta_K_rand_sub    = new TH1F("H_Beta_K_rand_sub",  "Kaon #beta (bkg); #beta; Counts / mC", hbeta_nbins, hbeta_xmin, hbeta_xmax);  
  H_Beta_Pi_rand_sub   = new TH1F("H_Beta_Pi_rand_sub", "Pion #beta (bkg); #beta; Counts / mC", hbeta_nbins, hbeta_xmin, hbeta_xmax);  
  H_Beta_P_rand_sub    = new TH1F("H_Beta_P_rand_sub",  "Proton #beta (bkg); #beta; Counts / mC", hbeta_nbins, hbeta_xmin, hbeta_xmax);
  
  //2D Histograms

  // Beta vs. Coincidence Time (no PID selection cuts)
  H_Beta_vs_ctime_K_noPID   = new TH2F("H_Beta_vs_ctime_K_noPID",   "#beta vs Kaon Coincidence Time (no PID); Kaon Coincidence Time [ns]; #beta", coin_nbins, coin_xmin, coin_xmax, hbeta_nbins, hbeta_xmin, hbeta_xmax); 
  H_Beta_vs_ctime_Pi_noPID  = new TH2F("H_Beta_vs_ctime_Pi_noPID", "#beta vs Pion Coincidence Time (no PID); Pion Coincidence Time [ns]; #beta", coin_nbins, coin_xmin, coin_xmax, hbeta_nbins, hbeta_xmin, hbeta_xmax); 
  H_Beta_vs_ctime_P_noPID   = new TH2F("H_Beta_vs_ctime_P_noPID",   "#beta vs Proton Coincidence Time (no PID); Proton Coincidence Time [ns]; #beta", coin_nbins, coin_xmin, coin_xmax, hbeta_nbins, hbeta_xmin, hbeta_xmax); 


  // Beta vs. Coincidence Time (PID selection cuts)
  H_Beta_vs_ctime_K_PID   = new TH2F("H_Beta_vs_ctime_K_PID",   "Kaon #beta vs Coincidence Time (PID); Coincidence Time [ns]; #beta_{K}", coin_nbins, coin_xmin, coin_xmax, hbeta_nbins, hbeta_xmin, hbeta_xmax); 
  H_Beta_vs_ctime_Pi_PID  = new TH2F("H_Beta_vs_ctime_Pi_PID", "Pion #beta vs Coincidence Time (PID); Coincidence Time [ns]; #beta_{#pi}", coin_nbins, coin_xmin, coin_xmax, hbeta_nbins, hbeta_xmin, hbeta_xmax); 
  H_Beta_vs_ctime_P_PID   = new TH2F("H_Beta_vs_ctime_P_PID",   "Proton #beta vs Coincidence Time (PID); Coincidence Time [ns]; #beta_{p}", coin_nbins, coin_xmin, coin_xmax, hbeta_nbins, hbeta_xmin, hbeta_xmax); 

  
  //Kaon LT Kinmatic Histograms
  H_Q2_vs_W   = new TH2F("H_Q2_vs_W", "Q^{2} vs. W", W_nbins, W_xmin, W_xmax, Q2_nbins, Q2_xmin, Q2_xmax);
  H_t_vs_phxq = new TH2F("Ht_vs_phxq", "-t Mandelstam vs. #phi_{xq}", phxq_nbins, phxq_xmin, phxq_xmax, MandelT_nbins, MandelT_xmin, MandelT_xmax);

  
  //2D: theta_xq_cm "VS" Phi_xq Histograms (These are important to extract beam asymmetry for "+" (pos) and "-"(neg) beam helicity) assuming Kaons (SHMS) are detected in coincidence with electrons (HMS)
  H_thxqCM_vs_phxq_K_pos = new TH2F("H_thxqCM_vs_phxq_K_pos", "#theta^{cm}_{xq} vs. #phi_{xq} (pos helicity K^{+}); #phi_{xq} [deg]; #theta^{cm}_{xq} [deg]", phxq_nbins, phxq_xmin, phxq_xmax, thxq_cm_nbins, thxq_cm_xmin, thxq_cm_xmax);
  H_thxqCM_vs_phxq_K_neg = new TH2F("H_thxqCM_vs_phxq_K_neg", "#theta^{cm}_{xq} vs. #phi_{xq} (neg helicity K^{+}); #phi_{xq} [deg]; #theta^{cm}_{xq} [deg]", phxq_nbins, phxq_xmin, phxq_xmax, thxq_cm_nbins, thxq_cm_xmin, thxq_cm_xmax);

  //same as above but to be filled with random coincidences
  H_thxqCM_vs_phxq_K_pos_rand = new TH2F("H_thxqCM_vs_phxq_K_pos_rand", "#theta^{cm}_{xq} vs. #phi_{xq} (pos helicity K^{+}, bkg); #phi_{xq} [deg]; #theta^{cm}_{xq} [deg]", phxq_nbins, phxq_xmin, phxq_xmax, thxq_cm_nbins, thxq_cm_xmin, thxq_cm_xmax);
  H_thxqCM_vs_phxq_K_neg_rand = new TH2F("H_thxqCM_vs_phxq_K_neg_rand", "#theta^{cm}_{xq} vs. #phi_{xq} (neg helicity K^{+}, bkg); #phi_{xq} [deg]; #theta^{cm}_{xq} [deg]", phxq_nbins, phxq_xmin, phxq_xmax, thxq_cm_nbins, thxq_cm_xmin, thxq_cm_xmax);

  //same as above
  H_thxqCM_vs_phxq_K_pos_rand_sub = new TH2F("H_thxqCM_vs_phxq_K_pos_rand_sub", "#theta^{cm}_{xq} vs. #phi_{xq} (pos helicity K^{+}, real); #phi_{xq} [deg]; #theta^{cm}_{xq} [deg]", phxq_nbins, phxq_xmin, phxq_xmax, thxq_cm_nbins, thxq_cm_xmin, thxq_cm_xmax);
  H_thxqCM_vs_phxq_K_neg_rand_sub = new TH2F("H_thxqCM_vs_phxq_K_neg_rand_sub", "#theta^{cm}_{xq} vs. #phi_{xq} (neg helicity K^{+}, real); #phi_{xq} [deg]; #theta^{cm}_{xq} [deg]", phxq_nbins, phxq_xmin, phxq_xmax, thxq_cm_nbins, thxq_cm_xmin, thxq_cm_xmax);

  //2D: theta_xq_cm "VS" Phi_xq Histograms (These are important to extract beam asymmetry for "+" (pos) and "-"(neg) beam helicity) assuming Pions (SHMS) are detected in coincidence with electrons (HMS)
  H_thxqCM_vs_phxq_Pi_pos = new TH2F("H_thxqCM_vs_phxq_Pi_pos", "#theta^{cm}_{xq} vs. #phi_{xq} (pos helicity #pi^{+}); #phi_{xq} [deg]; #theta^{cm}_{xq} [deg]", phxq_nbins, phxq_xmin, phxq_xmax, thxq_cm_nbins, thxq_cm_xmin, thxq_cm_xmax);
  H_thxqCM_vs_phxq_Pi_neg = new TH2F("H_thxqCM_vs_phxq_Pi_neg", "#theta^{cm}_{xq} vs. #phi_{xq} (neg helicity #pi^{+}); #phi_{xq} [deg]; #theta^{cm}_{xq} [deg]", phxq_nbins, phxq_xmin, phxq_xmax, thxq_cm_nbins, thxq_cm_xmin, thxq_cm_xmax);

  //same as above but to be filled with random coincidences
  H_thxqCM_vs_phxq_Pi_pos_rand = new TH2F("H_thxqCM_vs_phxq_Pi_pos_rand", "#theta^{cm}_{xq} vs. #phi_{xq} (pos helicity #pi^{+}, bkg); #phi_{xq} [deg]; #theta^{cm}_{xq} [deg]", phxq_nbins, phxq_xmin, phxq_xmax, thxq_cm_nbins, thxq_cm_xmin, thxq_cm_xmax);
  H_thxqCM_vs_phxq_Pi_neg_rand = new TH2F("H_thxqCM_vs_phxq_Pi_neg_rand", "#theta^{cm}_{xq} vs. #phi_{xq} (neg helicity #pi^{+}, bkg); #phi_{xq} [deg]; #theta^{cm}_{xq} [deg]", phxq_nbins, phxq_xmin, phxq_xmax, thxq_cm_nbins, thxq_cm_xmin, thxq_cm_xmax);

  //same as above
  H_thxqCM_vs_phxq_Pi_pos_rand_sub = new TH2F("H_thxqCM_vs_phxq_Pi_pos_rand_sub", "#theta^{cm}_{xq} vs. #phi_{xq} (pos helicity #pi^{+}, real); #phi_{xq} [deg]; #theta^{cm}_{xq} [deg]", phxq_nbins, phxq_xmin, phxq_xmax, thxq_cm_nbins, thxq_cm_xmin, thxq_cm_xmax);
  H_thxqCM_vs_phxq_Pi_neg_rand_sub = new TH2F("H_thxqCM_vs_phxq_Pi_neg_rand_sub", "#theta^{cm}_{xq} vs. #phi_{xq} (neg helicity #pi^{+}, real); #phi_{xq} [deg]; #theta^{cm}_{xq} [deg]", phxq_nbins, phxq_xmin, phxq_xmax, thxq_cm_nbins, thxq_cm_xmin, thxq_cm_xmax);


  //2D: theta_xq_cm "VS" Phi_xq Histograms (These are important to extract beam asymmetry for "+" (pos) and "-"(neg) beam helicity) assuming Protons (SHMS) are detected in coincidence with electrons (HMS)
  H_thxqCM_vs_phxq_P_pos = new TH2F("H_thxqCM_vs_phxq_P_pos", "#theta^{cm}_{xq} vs. #phi_{xq} (pos helicity protons); #phi_{xq} [deg]; #theta^{cm}_{xq} [deg]", phxq_nbins, phxq_xmin, phxq_xmax, thxq_cm_nbins, thxq_cm_xmin, thxq_cm_xmax);
  H_thxqCM_vs_phxq_P_neg = new TH2F("H_thxqCM_vs_phxq_P_neg", "#theta^{cm}_{xq} vs. #phi_{xq} (neg helicity protons); #phi_{xq} [deg]; #theta^{cm}_{xq} [deg]", phxq_nbins, phxq_xmin, phxq_xmax, thxq_cm_nbins, thxq_cm_xmin, thxq_cm_xmax);

  //same as above but to be filled with random coincidences
  H_thxqCM_vs_phxq_P_pos_rand = new TH2F("H_thxqCM_vs_phxq_P_pos_rand", "#theta^{cm}_{xq} vs. #phi_{xq} (pos helicity protons, bkg); #phi_{xq} [deg]; #theta^{cm}_{xq} [deg]", phxq_nbins, phxq_xmin, phxq_xmax, thxq_cm_nbins, thxq_cm_xmin, thxq_cm_xmax);
  H_thxqCM_vs_phxq_P_neg_rand = new TH2F("H_thxqCM_vs_phxq_P_neg_rand", "#theta^{cm}_{xq} vs. #phi_{xq} (neg helicity protons, bkg); #phi_{xq} [deg]; #theta^{cm}_{xq} [deg]", phxq_nbins, phxq_xmin, phxq_xmax, thxq_cm_nbins, thxq_cm_xmin, thxq_cm_xmax);

  //same as above
  H_thxqCM_vs_phxq_P_pos_rand_sub = new TH2F("H_thxqCM_vs_phxq_P_pos_rand_sub", "#theta^{cm}_{xq} vs. #phi_{xq} (pos helicity protons, real); #phi_{xq} [deg]; #theta^{cm}_{xq} [deg]", phxq_nbins, phxq_xmin, phxq_xmax, thxq_cm_nbins, thxq_cm_xmin, thxq_cm_xmax);
  H_thxqCM_vs_phxq_P_neg_rand_sub = new TH2F("H_thxqCM_vs_phxq_P_neg_rand_sub", "#theta^{cm}_{xq} vs. #phi_{xq} (neg helicity protons, real); #phi_{xq} [deg]; #theta^{cm}_{xq} [deg]", phxq_nbins, phxq_xmin, phxq_xmax, thxq_cm_nbins, thxq_cm_xmin, thxq_cm_xmax);

  
  //--Add Histograms from helicityAnalyzer here!--
  
  hel_HList->Add(H_eK_ctime_real);  
  hel_HList->Add(H_ePi_ctime_real); 
  hel_HList->Add(H_ep_ctime_real);  
                   
  hel_HList->Add(H_eK_ctime_rand);  
  hel_HList->Add(H_ePi_ctime_rand); 
  hel_HList->Add(H_ep_ctime_rand);  

  hel_HList->Add(H_MM_K_PID);  
  hel_HList->Add(H_MM_Pi_PID); 
  hel_HList->Add(H_MM_P_PID);
  
  hel_HList->Add(H_MM_K_real);
  hel_HList->Add(H_MM_Pi_real);
  hel_HList->Add(H_MM_P_real);

  hel_HList->Add(H_MM_K_rand);
  hel_HList->Add(H_MM_Pi_rand);
  hel_HList->Add(H_MM_P_rand);

  hel_HList->Add(H_MM_K_rand_sub);
  hel_HList->Add(H_MM_Pi_rand_sub);
  hel_HList->Add(H_MM_P_rand_sub);

  hel_HList->Add(H_Beta_noPID);

  hel_HList->Add(H_Beta_K_PID);
  hel_HList->Add(H_Beta_Pi_PID);
  hel_HList->Add(H_Beta_P_PID);

  hel_HList->Add(H_Beta_K_real);
  hel_HList->Add(H_Beta_Pi_real);
  hel_HList->Add(H_Beta_P_real);

  hel_HList->Add(H_Beta_K_rand);
  hel_HList->Add(H_Beta_Pi_rand);
  hel_HList->Add(H_Beta_P_rand);

  hel_HList->Add(H_Beta_K_rand_sub);
  hel_HList->Add(H_Beta_Pi_rand_sub);
  hel_HList->Add(H_Beta_P_rand_sub);
  
  hel_HList->Add(H_Beta_vs_ctime_K_noPID);
  hel_HList->Add(H_Beta_vs_ctime_Pi_noPID);
  hel_HList->Add(H_Beta_vs_ctime_P_noPID);

  hel_HList->Add(H_Beta_vs_ctime_K_PID);
  hel_HList->Add(H_Beta_vs_ctime_Pi_PID);
  hel_HList->Add(H_Beta_vs_ctime_P_PID);
  
  hel_HList->Add(H_Q2_vs_W);
  hel_HList->Add(H_t_vs_phxq);

  hel_HList->Add(H_thxqCM_vs_phxq_K_pos);
  hel_HList->Add(H_thxqCM_vs_phxq_K_neg);

  hel_HList->Add(H_thxqCM_vs_phxq_K_pos_rand);
  hel_HList->Add(H_thxqCM_vs_phxq_K_neg_rand);

  hel_HList->Add(H_thxqCM_vs_phxq_K_pos_rand_sub);
  hel_HList->Add(H_thxqCM_vs_phxq_K_neg_rand_sub);

  hel_HList->Add(H_thxqCM_vs_phxq_Pi_pos);
  hel_HList->Add(H_thxqCM_vs_phxq_Pi_neg);

  hel_HList->Add(H_thxqCM_vs_phxq_Pi_pos_rand);
  hel_HList->Add(H_thxqCM_vs_phxq_Pi_neg_rand);

  hel_HList->Add(H_thxqCM_vs_phxq_Pi_pos_rand_sub);
  hel_HList->Add(H_thxqCM_vs_phxq_Pi_neg_rand_sub);

  hel_HList->Add(H_thxqCM_vs_phxq_P_pos);
  hel_HList->Add(H_thxqCM_vs_phxq_P_neg);

  hel_HList->Add(H_thxqCM_vs_phxq_P_pos_rand);
  hel_HList->Add(H_thxqCM_vs_phxq_P_neg_rand);

  hel_HList->Add(H_thxqCM_vs_phxq_P_pos_rand_sub);
  hel_HList->Add(H_thxqCM_vs_phxq_P_neg_rand_sub);
}


//_______________________________________________________________________________
void helicityAnalyzer::ReadScalerTree()
{
  cout << "Calling Derived ReadScalerTree() " << endl;

  //As of now, there is nothing to add to the existing method from baseAnalyzer.

  //Call the existing method, as it is needed to read the scaler tree
  baseAnalyzer::ReadScalerTree();

  //To add more scaler variables, one needs to add them to the baseAnalyzer.cpp
  
}

//_______________________________________________________________________________
void helicityAnalyzer::ScalerEventLoop()
{
  cout << "Calling Derived ScalerEventLoop() " << endl;

  //As of now, there is nothing to add to the existing method from baseAnalyzer.

  //Call the existing method, as it is needed to loop over the scaler reads 
  baseAnalyzer::ScalerEventLoop();

  //The ScalerEventLoop(), like the data EventLoop() can either be called from the
  //baseAnalyzer, as the line above, or completely re-defined, as it was done in
  //the helicityAnalyzer::EventLoop(), otherwise, if one tries to add functionality
  //in the event loop after calling the generic baseAnalyzer event loop, there are
  //issues reading the tree variables that were inherited. 
  
}

//_______________________________________________________________________________
void helicityAnalyzer::ReadTree()
{
  cout << "Calling Derived ReadTree() " << endl;


  //Call the existing method, as it is needed to read the generic data tree leafs
  //baseAnalyzer::ReadTree(); //Calling this method to get TTree leafs in the derived method
  //does NOT work.  One needs to explicityly define the TTree leafs from the base method here
  //so thay they may be called properly in the EventLoop()

  baseAnalyzer::ReadTree();

  //The additional helicity tree variablles are added in the baseAnalyzer.cpp, and used
  //via a 'helicity_flag' parameter via the main_controls.inp


  
}

//_______________________________________________________________________________
void helicityAnalyzer::EventLoop()
{

  /*
    Brief: This method should be specific to the derived class in question so thay they may fill the
    necessary histograms by said class. For instace, this helicity class requires fewer histograms to 
    be filled than the generic class. Some variables from the generic class, like detector variables,
    are necessary and are used in this methods via inheritance of the baseAnalyzer class.
  */
  
  cout << "Calling Derived EventLoop() . . . " << endl;
  
  //Loop over Events
  
  if(analysis=="data")
    {
      cout << "Loop over Data Events | nentries -->  " << nentries << endl;
      
      for(int ientry=0; ientry<nentries; ientry++)
	{
	  
	  tree->GetEntry(ientry);
	  
	  //--------------DETERMINE THE 'SIGN' OF THE BEAM HELICITY PER EVENT-------------------
	  /*
	    Brief: The electron beam bunches beam helicity (polarization direction) can be in a
	    '+ - - +' or '- + + -' 4-bunch quartet structure. In addition, changes to the IHWP correspond
	    to changes in the helicity sign. These IHWP changes are done for systematic studies.
	    The helicity sign will be used as a cut to require events that came from either a '+' or '-' 
	    helicity beam bunch, and that way, be able to extract the beam-spin asymmetry: Asy = (N+ - N-) / (N+ + N-)
	  */

	  hel_sign = 0;   //set default sign to 0 (undetermined helicity state)
	  
	  //MPS state = 1 represents the T_settle time, during which beam buch polarization is undefined/undeetermiend
	  //We require MPS < 0.5 (which means, basically, MPS = 0), where the beam helicity is known (either +1 or -1)
	  if(hel_mps < 0.5)
	    { 
	      //set to "+" helicity
	      if(hel>0.5)
		{
		  hel_sign = 1;
		}
	      //set to "-" helicity
	      else if(hel<-0.5)
		{
		  hel_sign = -1;
		}
	      //Require hel_sign flip corresponding to IHWP configuration change for certain run range
	      if(((run > 5033) && (run < 5155)) ||
		 ((run > 5335) && (run < 5668)) ||
		 ((run > 5962) && (run < 6147)))
		{
		  hel_sign = -hel_sign;
		}
	      
	    }
	  
	  
	  //------------------------------------------------------------------------------------

	  
	  
	  //--------------CALCULATED KINEMATICS VARIABLES (IF THEY ARE NOT ALREADY DONE IN HCANA)-----------

	  th_x = xangle - th_e;  //hadron arm central angle for each particle
	  MM2 = MM*MM;           //Missing Mass Squared
 	  ztar_diff = htar_z - etar_z;  //reaction vertex z difference
	  
	  
	  //Missing Mass Squared for Kaons, Pions or Protons (obtained definition from S.A.W beam helicity scripts)
	  MM2_K = pow(Em,2) - pow(Pm,2);
	  MM2_Pi = pow( Em + sqrt( pow(MK, 2) + pow(Pf,2) ) - sqrt(pow(MPi,2) + pow(Pf,2)), 2) - pow(Pm,2);
	  MM2_P = pow( Em + sqrt( pow(MK, 2) + pow(Pf,2) ) - sqrt(pow(MP,2) + pow(Pf,2)), 2) - pow(Pm,2);

	  //Set default to -1000
	  MM_K = -1000.;
	  MM_P = -1000.;
	  MM_P = -1000.;
	  
	  //Require MM2 > 0
	  if(MM2_K>0)  { MM_K = sqrt(MM2_K);   }
	  if(MM2_Pi>0) { MM_Pi = sqrt(MM2_Pi); }
	  if(MM2_P>0)  { MM_P = sqrt(MM2_P);   }

	    
	  //--------------DEFINE CUTS--------------------
	  // ** NOTE ** : The cuts defined here may be copied from the generic analyzer, or re-defined as needed
	  
	  
	  //CUTS USED IN EDTM LIVE TIME CALCULATION
	  c_noedtm = EDTM_tdcTimeRaw == 0.;
	  c_edtm   = EDTM_tdcTimeRaw  > 0.;
	  c_trig1  = TRIG1_tdcTimeRaw > 0.;
	  c_trig2  = TRIG2_tdcTimeRaw > 0.;
	  c_trig3  = TRIG3_tdcTimeRaw > 0.;
	  c_trig4  = TRIG4_tdcTimeRaw > 0.;
	  c_trig5  = TRIG5_tdcTimeRaw > 0.;
	  c_trig6  = TRIG6_tdcTimeRaw > 0.;
	  
	  //=====CUTS USED IN TRACKING EFFICIENCY CALCULATION=====
	  //(See main_controls.inp to determine which cuts file is being loaded, and modify the cuts accordingly)
	  //CUTS: HMS TRACKING EFFICIENCY (May be for e- or hadrons, depending on the limits set in the input file)

	  //Require at least a minimum number of track(s) 
	  if(hdc_ntrk_cut_flag){c_hdc_ntrk = hdc_ntrack >= c_hdc_ntrk_min;}
	  else{
	    cout <<
	      "********************************\n"
	      "TRACKING EFFICIENCY ERROR: \n"
	      "Must set hdc_ntrk_cut_flag = 1 \n"
	      "See " <<  Form("%s", input_CutFileName.Data()) << "\n"
	      "********************************"<< endl;
	    
	    //"*********************************************************" << endl;
	    gSystem->Exit(0);}

	  //Require a "Good Scintillator Hit" in the Fiducial Hodoscopoe Region
	  if(hScinGood_cut_flag){c_hScinGood = hhod_GoodScinHit==1;}
	  else{c_hScinGood=1;} //1 means allow all events (do not put hScinGood cut)

	  //Require HMS Cherenkov Cut (for electron or hadron selection)
	  if(hcer_cut_flag){c_hcer_NPE_Sum = hcer_npesum >= c_hnpeSum_min && hcer_npesum <= c_hnpeSum_max;}
	  else{c_hcer_NPE_Sum=1;}

	  //Require HMS Calorimeter Cut (for additional electron or hadron selection)
	  if(hetotnorm_cut_flag){c_hetotnorm = hcal_etotnorm >= c_hetotnorm_min && hcal_etotnorm <= c_hetotnorm_max;}
	  else{c_hetotnorm=1;}

	  //Require HMS Hodoscope Beta Cut (no track-biased) (for electron or hadron selection)
	  if(hBeta_notrk_cut_flag){c_hBeta_notrk = hhod_beta_ntrk >= c_hBetaNtrk_min && hhod_beta_ntrk <= c_hBetaNtrk_max;}
	  else{c_hBeta_notrk=1;}

	  //electrons (or hadrons) that 'SHOULD' have passed the cuts to form a track
	  good_hms_should = c_hScinGood && c_hcer_NPE_Sum && c_hetotnorm && c_hBeta_notrk;

	  //electrons (or hadrons) that 'DID' passed the cuts to form a track
	  good_hms_did = c_hdc_ntrk && good_hms_should;


	  //CUTS: SHMS TRACKING EFFICIENCY (May be for e- or hadrons, depending on the limits set in the input file)

	  //Require at least a minimum number of track(s) 
	  if(pdc_ntrk_cut_flag){c_pdc_ntrk = pdc_ntrack >= c_pdc_ntrk_min;}
	  else{
	    cout <<
	      "********************************\n"
	      "TRACKING EFFICIENCY ERROR: \n"
	      "Must set pdc_ntrk_cut_flag = 1 \n"
	      "See " <<  Form("%s", input_CutFileName.Data()) << "\n"
	      "********************************"<< endl;
	    
	    //"*********************************************************" << endl;
	    gSystem->Exit(0);}

	  //Require a "Good Scintillator Hit" in the Fiducial Hodoscopoe Region
	  if(pScinGood_cut_flag){c_pScinGood = phod_GoodScinHit==1;}
	  else{c_pScinGood=1;} //1 means allow all events (do not put hScinGood cut)

	  //Require SHMS Noble Gas Cherenkov Cut (for electron or hadron selection)
	  if(pngcer_cut_flag){c_pngcer_NPE_Sum = pngcer_npesum >= c_pngcer_npeSum_min &&  pngcer_npesum <= c_pngcer_npeSum_max;}
	  else{c_pngcer_NPE_Sum=1;}

	  //Require SHMS Heavy Gas Cherenkov Cut (for electron or hadron selection)
	  if(phgcer_cut_flag){c_phgcer_NPE_Sum = phgcer_npesum >= c_phgcer_npeSum_min &&  phgcer_npesum <= c_phgcer_npeSum_max;}
	  else{c_phgcer_NPE_Sum=1;}

	  //WHAT ABOUT AN AEROGEL CUT IN THE TRACKING EFFICIENCY?
	  
	  //Require SHMS Calorimeter Cut (for additional electron or hadron selection)
	  if(petotnorm_cut_flag){c_petotnorm = pcal_etotnorm >= c_petotnorm_min && pcal_etotnorm <= c_petotnorm_max;}
	  else{c_petotnorm=1;}

	  //Require SHMS Hodoscope Beta Cut (no track-biased) (for electron or hadron selection)
	  if(pBeta_notrk_cut_flag){c_pBeta_notrk = phod_beta_ntrk >= c_pBetaNtrk_min && phod_beta_ntrk <= c_pBetaNtrk_max;}
	  else{c_pBeta_notrk=1;}

	  //electrons (or hadrons) that 'SHOULD' have passed the cuts to form a track
	  good_shms_should = c_pScinGood && c_pngcer_NPE_Sum && c_phgcer_NPE_Sum && c_petotnorm && c_pBeta_notrk;

	  //electrons (or hadrons) that 'DID' passed the cuts to form a track
	  good_shms_did = c_pdc_ntrk && good_shms_should;

	  //=====END: CUTS USED IN TRACKING EFFICIENCY CALCULATION=====


	  //Count Accepted EDTM events (no bcm current cut : this is just to compare with counts that have bcm cuts)
	  if(c_edtm){total_edtm_accp++;}
	  
	  //Count Accepted TRIG 1-6 events (no bcm current cut : this is just to compare with counts that have bcm cuts)
	  if(c_trig1){total_trig1_accp++;}
	  if(c_trig2){total_trig2_accp++;}
	  if(c_trig3){total_trig3_accp++;}
	  if(c_trig4){total_trig4_accp++;}
	  if(c_trig5){total_trig5_accp++;}
	  if(c_trig6){total_trig6_accp++;}

	  //----------------------Check If BCM Current is within limits---------------------


	  if(evt_flag_bcm[scal_read]==1)
	    {

	      if(c_edtm){ total_edtm_accp_bcm_cut++;}
	      
	      //Count Accepted TRIG1-6 events (without EDTM and with bcm current cut: to be used in the computer live time calculation)
	      if(c_trig1 && c_noedtm) { total_trig1_accp_bcm_cut++; }
	      if(c_trig2 && c_noedtm) { total_trig2_accp_bcm_cut++; }
	      if(c_trig3 && c_noedtm) { total_trig3_accp_bcm_cut++; }
	      if(c_trig4 && c_noedtm) { total_trig4_accp_bcm_cut++; }
	      if(c_trig5 && c_noedtm) { total_trig5_accp_bcm_cut++; }
	      if(c_trig6 && c_noedtm) { total_trig6_accp_bcm_cut++; }
	      
	      //REQUIRE "NO EDTM" CUT TO FILL DATA HISTOGRAMS
	      if(c_noedtm)
		{

		  //---Calculate HMS Tracking Efficiency Components---

		  //Integrated over +/- helicities
		  if(good_hms_did){ h_did++;}
		  if(good_hms_should){ h_should++; }

		  //HMS: Tracking Efficiency (positive helicity)
		  if(good_hms_did && hel_sign>0){ h_did_pos++;}
		  if(good_hms_should && hel_sign>0){ h_should_pos++; }
		  
		  //HMS: Tracking Efficiency (negative helicity)
		  if(good_hms_did && hel_sign<0){ h_did_neg++;}
		  if(good_hms_should && hel_sign<0){ h_should_neg++; }
		  
		  //---Calculate SHMS Tracking Efficiency Components---

		  //Integrated over +/- helicities
		  if(good_shms_did){ p_did++;}
		  if(good_shms_should){ p_should++; }

		   //SHMS: Tracking Efficiency (positive helicity)
		  if(good_shms_did && hel_sign>0){ p_did_pos++;}
		  if(good_shms_should && hel_sign>0){ p_should_pos++; }

		  //SHMS: Tracking Efficiency (negative helicity)
		  if(good_shms_did && hel_sign<0){ p_did_neg++;}
		  if(good_shms_should && hel_sign<0){ p_should_neg++; }
		  
		  //====DATA ANALYSIS CUTS (MUST BE EXACTLY SAME AS SIMC, except PID CUTS on detectors)====
		  
		  // See main_controls.inp to look at which cuts file is currently being read into the code.
		  // In the cuts file itself, there is a variety of cuts with their flags (ON/OFF) and ranges
		  // which are read into the baseAnalyzer class. Since this class inherits from the baseAnalyzer,
		  // the cuts are set via the method call "baseAnalyzer::SetCuts()" in this code.
		  
		  // Below are a range of cuts needed for the beam asymmetry analysis. The naming scheme is exactly the
		  // same as in the input cuts file to avoid confusion. The cuts below are applied directly without the use
		  // of a flag. See 'set_basic_cuts_KaonLT.inp' to set the cut ranges below.


		  // Define the beam asymmetry boolean flags to define cuts to be used later on. The limits read in here can be searched quickly
		  // via a grep command on the corresponding input cuts file (set_basic_cuts_KaonLT.inp)

		  //**NOTE: When a cut is defaulted to 1, it means that when it is applied, the event will always pass the cut (hence, =1 --> is always true, which allows events to pass and therefore, the CUT is not applied)

		  //----Acceptance Cuts----

		  if(hdelta_cut_flag){c_hdelta = h_delta>=hdelta_min && h_delta<=hdelta_max;} 
		  else{c_hdelta=1;}  
		  		  
		  if(edelta_cut_flag){c_edelta = e_delta>=edelta_min && e_delta<=edelta_max;} 
		  else{c_edelta=1;}
		  
		  if(ztarDiff_cut_flag){c_ztarDiff = ztar_diff>=ztarDiff_min && ztar_diff<=ztarDiff_max;} 
		  else{c_ztarDiff=1;} 

		  //---- Particle Identification Cuts ---

		  //--electron--
		  
		  if(hcer_pidCut_flag){cpid_hcer_NPE_Sum = hcer_npesum>=elec_hcer_npe_thrs;}
		  else{cpid_hcer_NPE_Sum = 1;}
		  
		  if(hetot_trkNorm_pidCut_flag){cpid_hetot_trkNorm =  hcal_etottracknorm>=elec_hcal_thrs;}
		  else{cpid_hetot_trkNorm = 1;}
		  
		  //-- Kaon --
		  if(K_paero_npe_flag) {c_K_paero_npe = paero_npesum>=K_paero_npe_thrs;}
		  else{c_K_paero_npe = 1;}

		  if(K_phgcer_npe_flag) {c_K_phgcer_npe = phgcer_npesum<=K_phgcer_npe_thrs;}
		  else{c_K_phgcer_npe = 1;}

		  if(K_beta_flag) { kaon_beta_cut = phod_gtr_beta>=K_beta_min && phod_gtr_beta<=K_beta_max; }
		  else{ kaon_beta_cut=1; }

		  if(eK_ctime_flag) {
		    eK_ctime_cut = abs(eKCoinTime-K_ctime_offset) <= eK_ctime_thrs;
		    eK_ctime_cut_rand = abs(eKCoinTime-K_ctime_offset) > eK_ctime_thrs && abs(eKCoinTime-K_ctime_offset) <= (eK_mult*eK_ctime_thrs);
		  }
		  else{
		    eK_ctime_cut = 1;
		    eK_ctime_cut_rand = 1;
		  }

		  if(MM_K_cut_flag) { kaon_MM_cut = MM_K >=MM_K_min && MM_K <=MM_K_max; }
		  else { kaon_MM_cut = 1;}

		  //-- Pion --
		  if(Pi_phgcer_npe_flag) {c_Pi_phgcer_npe =  phgcer_npesum>=Pi_phgcer_npe_thrs;}
		  else {c_Pi_phgcer_npe = 1;}
		  
		  if(Pi_beta_flag) { pion_beta_cut = phod_gtr_beta>=Pi_beta_min && phod_gtr_beta<=Pi_beta_max; }
		  else{ pion_beta_cut=1; }

		  if(ePi_ctime_flag) {
		    ePi_ctime_cut = abs(ePiCoinTime-Pi_ctime_offset) <= ePi_ctime_thrs;
		    ePi_ctime_cut_rand = abs(ePiCoinTime-Pi_ctime_offset) > ePi_ctime_thrs && abs(ePiCoinTime-Pi_ctime_offset) <= (ePi_mult*ePi_ctime_thrs);
		  }
		  else{
		    ePi_ctime_cut = 1;
		    ePi_ctime_cut_rand = 1;
		  }

		  if(MM_Pi_cut_flag) { pion_MM_cut = MM_Pi >=MM_Pi_min && MM_Pi <=MM_Pi_max; }
		  else { pion_MM_cut = 1;}

		  //-- Proton --
		  if(P_paero_npe_flag) {c_P_paero_npe = paero_npesum<=P_paero_npe_thrs;}
		  else{c_P_paero_npe = 1;}

		  if(P_phgcer_npe_flag) {c_P_phgcer_npe =  phgcer_npesum<=P_phgcer_npe_thrs;}
		  else {c_P_phgcer_npe = 1;}
		  
		  if(P_beta_flag) { proton_beta_cut = phod_gtr_beta>=P_beta_min && phod_gtr_beta<=P_beta_max; }
		  else{ proton_beta_cut=1; }

		  if(eP_ctime_flag) {
		    eP_ctime_cut = abs(epCoinTime-P_ctime_offset) <= eP_ctime_thrs;
		    eP_ctime_cut_rand = abs(epCoinTime-P_ctime_offset) > eP_ctime_thrs && abs(epCoinTime-P_ctime_offset) <= (eP_mult*eP_ctime_thrs);
		  }
		  else{
		    eP_ctime_cut = 1;
		    eP_ctime_cut_rand = 1;
		  }

		  if(MM_P_cut_flag) { proton_MM_cut = MM_P >=MM_P_min && MM_P <=MM_P_max; }
		  else { proton_MM_cut = 1;}

		  
		  //Acceptance Cuts
		  accp_cuts = c_edelta && c_hdelta && c_ztarDiff;

		  //electron PID Cuts
		  elec_pid = cpid_hcer_NPE_Sum && cpid_hetot_trkNorm;

		  // --- Hadron PID Cuts (using SHMS HGC and AERO Cherenkovs) ----

		  //Kaon Analysis PID Cuts: H(e,e'K+)X,  
		  kaon_pid = c_K_paero_npe && c_K_phgcer_npe;		  
		  		  
		  //Pion Analysis PID Cuts: H(e,e'Pi+)X
		  pion_pid = c_Pi_phgcer_npe;
		  
		  //Proton Analysis PID Cuts: H(e,e'p)X,  
		  proton_pid = c_P_paero_npe && c_P_phgcer_npe;
		  
		  //----------------------Fill DATA Histograms-----------------------


		  //--------------------------------------------------------------------
		  //---------HISTOGRAM CATEGORY: Particle Identification (PID)----------
		  //               INHERITED FROM baseAnalyzer CLASS 
		  //--------------------------------------------------------------------

		  // C.Y. NEED TO: Make Spectrometer Acceptance Cuts here (BASED ON THE 1D ACCEPTANCE HISTOGRAMS FILLED BELOW)
		  
		  //Apply acceptance cut (else, continue to next event within loop)
		  if(accp_cuts) { 		  
		    
		    //Fill all coincidence times to determine the offset for each
		    H_eK_ctime->Fill(eKCoinTime-K_ctime_offset); //this ctime (inherited from baseAnalyzer) is used to set the K_ctime_offset and cuts)
		    H_ePi_ctime->Fill(ePiCoinTime-Pi_ctime_offset);
		    H_ep_ctime->Fill(epCoinTime-P_ctime_offset);
		    
		  
		    //Fill HMS Detectors
		    H_hCerNpeSum->Fill(hcer_npesum);
		    H_hCalEtotNorm->Fill(hcal_etotnorm);
		    H_hCalEtotTrkNorm->Fill(hcal_etottracknorm);
		    H_hHodBetaNtrk->Fill(hhod_beta_ntrk);
		    H_hHodBetaTrk->Fill(hhod_gtr_beta);
		    
		    //Fill SHMS Detectors
		    H_pHGCerNpeSum->Fill(phgcer_npesum);
		    H_pAeroNpeSum->Fill(paero_npesum);
		    H_pCalEtotNorm->Fill(pcal_etotnorm);
		    H_pCalEtotTrkNorm->Fill(pcal_etottracknorm);
		    H_pHodBetaNtrk->Fill(phod_beta_ntrk);
		    H_pHodBetaTrk->Fill(phod_gtr_beta);		  

		    //Fill 2D PID Correlations
		    H_hcal_vs_hcer->Fill(hcal_etottracknorm, hcer_npesum);
		    H_pcal_vs_phgcer->Fill(pcal_etottracknorm, phgcer_npesum);  
		    H_pcal_vs_paero->Fill(pcal_etottracknorm, paero_npesum);   
		    H_paero_vs_phgcer->Fill(paero_npesum, phgcer_npesum); 	      		  
		    
		    //----------------------------------------------------------------------
		    //---------HISTOGRAM CATEGORY: Spectrometer Acceptance  (ACCP)----------
		    //----------------------------------------------------------------------
		    
		    //Fill SPECTROMETER  ACCEPTANCE (focal plane)
		    H_exfp       ->Fill(e_xfp);
		    H_eyfp       ->Fill(e_yfp);
		    H_expfp      ->Fill(e_xpfp);
		    H_eypfp      ->Fill(e_ypfp);
		    
		    H_hxfp       ->Fill(h_xfp);
		    H_hyfp       ->Fill(h_yfp);
		    H_hxpfp      ->Fill(h_xpfp);
		    H_hypfp      ->Fill(h_ypfp);
		    
		    //Fill SPECTROMETER  ACCEPTANCE (reconstructed)		  
		    H_eytar      ->Fill(e_ytar);
		    H_exptar     ->Fill(e_xptar);
		    H_eyptar     ->Fill(e_yptar);
		    H_edelta     ->Fill(e_delta);
		    
		    H_hytar       ->Fill(h_ytar);
		    H_hxptar      ->Fill(h_xptar);
		    H_hyptar      ->Fill(h_yptar);
		    H_hdelta      ->Fill(h_delta);
		    
		    //Fill SPECTROMETER AT REACTION VERTEX
		    H_htar_x       ->Fill(htar_x);
		    H_htar_y       ->Fill(htar_y);
		    H_htar_z       ->Fill(htar_z);
		    H_etar_x       ->Fill(etar_x);
		    H_etar_y       ->Fill(etar_y);
		    H_etar_z       ->Fill(etar_z);


		    //--------------------------------------------------------
		    //---------HISTOGRAM CATEGORY: Kinematics  (KIN)----------
		    //--------------------------------------------------------
		    
		    //Fill Primary Kin Histos
		    H_the    ->Fill(th_e/dtr);
		    H_kf     ->Fill(kf);
		    H_W      ->Fill(W);
		    H_W2     ->Fill(W2);
		    H_Q2     ->Fill(Q2);
		    H_xbj    ->Fill(X);
		    H_nu     ->Fill(nu);
		    H_q      ->Fill(q);
		    H_qx     ->Fill(qx);
		    H_qy     ->Fill(qy);
		    H_qz     ->Fill(qz);
		    H_thq    ->Fill(th_q/dtr);
		    H_phq    ->Fill(ph_q/dtr);
		    H_epsilon->Fill(epsilon); 
		    
		    //Fill Secondary Kin Histos
		    H_Em       ->Fill(Em);
		    H_Pm       ->Fill(Pm);
		    H_Pmx_lab  ->Fill(Pmx_lab);
		    H_Pmy_lab  ->Fill(Pmy_lab);
		    H_Pmz_lab  ->Fill(Pmz_lab);
		    H_Pmx_q    ->Fill(Pmx_q);
		    H_Pmy_q    ->Fill(Pmy_q);
		    H_Pmz_q    ->Fill(Pmz_q);
		    H_Tx       ->Fill(Tx);
		    H_Tr       ->Fill(Tr);
		    H_thx      ->Fill(th_x/dtr);
		    H_Pf       ->Fill(Pf);		  
		    H_Tx_cm    ->Fill(Tx_cm);
		    H_Tr_cm    ->Fill(Tr_cm);		  
		    
		    //pid PID cuts Beta (redundant hsito from above, but with different name)
		    H_Beta_noPID->Fill(phod_gtr_beta);
		    
		    
		    //--------Start (Kaon LT) Helicity Analysis----------
		    
		    
		    //=====================================
		    //  KAON BEAM-SPIN ASYMMETRY ANALYSIS
		    //=====================================
		    
		    //Fill Beta vs. Coin Time (Assuming Kaon Mass in Coin Time calculation)
		    H_Beta_vs_ctime_K_noPID->Fill((eKCoinTime-K_ctime_offset), phod_gtr_beta);
		    
		    // C.Y. NEED TO: Make Kaon selection PID CUT HERE ! (BASED ON THE 2D PID CORRELATIONS FILLED ABOVE)
		    if(kaon_pid && elec_pid) {
		      
		      //C.Y. NEED TO: Make a Beta Cut HERE! (Based on the 1D and 2D Beta plots)
		      if(kaon_beta_cut) {		    
			
			//Fill Beta (Kaon PID)
			H_Beta_K_PID->Fill(phod_gtr_beta);
			H_Beta_vs_ctime_K_PID->Fill((eKCoinTime-K_ctime_offset), phod_gtr_beta);
			
			//Missing Mass (Kaon PID)
			H_MM_K_PID->Fill(MM_K);
			
			//------------------------------------------
			//           TRUE KAON SELECTION
			//------------------------------------------		    
			
			//C.Y. NEED TO: Make coin. time cut to select TRUE eK coin. time peak HERE !
			if(eK_ctime_cut) {			
			  
			  H_eK_ctime_real->Fill(eKCoinTime-K_ctime_offset);			  
			  H_MM_K_real->Fill(MM_K); 
			  
			  //C.Y. NEED TO: Make a Kaon  Missing Mass Cut HERE !
			  if(kaon_MM_cut){			  
			    
			    //Fill Beta (with missing mass cut)
			    H_Beta_K_real->Fill(phod_gtr_beta);
			    
			    //Fill Kaon LT Kinematic Coverage Histos (after selecting Kaon Missing Mass)
			    H_Q2_vs_W->Fill(W, Q2);
			    H_t_vs_phxq->Fill(ph_xq/dtr, -MandelT);
			    H_epsilon->Fill(epsilon);   
			    
			    H_ztar_diff->Fill(ztar_diff);
			    
			    //Check if event has "+" helicity
			    if(hel_sign>0){
			      //Fill Fill In-plane th_xq_cm "VS" Out-of-plane angle, ph_xq (pos helicity) for the Kaons 
			      H_thxqCM_vs_phxq_K_pos->Fill(ph_xq/dtr, th_xq_cm/dtr);   //in deg			    
			    }
			    
			    //Check if event has "-" helicity
			    else if(hel_sign<0){
			      //Fill Fill In-plane th_xq_cm "VS" Out-of-plane angle, ph_xq (neg helicity) for the Kaons
			      H_thxqCM_vs_phxq_K_neg->Fill(ph_xq/dtr, th_xq_cm/dtr);   //in deg			      
			    }
			    
			  } // end missing mass cut around Lambda (1115), to select reaction channel, H(e,e'K+)MM, where MM=Lambda
			  
			} // end eK real coin. selection cut		     
			
			//------------------------------------------
			//         RANDOM BKG KAON SELECTION
			//------------------------------------------
			
			//C.Y. NEED TO: Make coin. time cut to select randoms (sample outside main peak) eK coin. time peak HERE !
			if(eK_ctime_cut_rand) {		      
			  
			  H_eK_ctime_rand->Fill(eKCoinTime-K_ctime_offset);
			  H_MM_K_rand->Fill(MM_K);			
			  
			  //C.Y. NEED TO: Make a Kaon  Missing Mass Cut HERE !
			  if(kaon_MM_cut){			  
			    
			    //Fill Beta (with missing mass cut)
			    H_Beta_K_rand->Fill(phod_gtr_beta);
			    
			    //Check if event has "+" helicity
			    if(hel_sign>0){
			      //Fill Fill In-plane th_xq_cm "VS" Out-of-plane angle, ph_qx (pos helicity) for the Kaons
			      H_thxqCM_vs_phxq_K_pos_rand->Fill(ph_xq/dtr, th_xq_cm/dtr);  //in deg			    
			    }
			    
			  //Check if event has "-" helicity
			    else if(hel_sign<0){
			      //Fill Fill In-plane th_xq_cm "VS" Out-of-plane angle, ph_qx (neg helicity) for the Kaons
			      H_thxqCM_vs_phxq_K_neg_rand->Fill(ph_xq/dtr, th_xq_cm/dtr);  //in deg			    
			    }
			    
			  } // end missing mass cut around Lambda (1115), to select reaction channel, H(e,e'K+)MM, where MM=Lambda
			  
			} // end eK random coin. selection
			
		      } // end Kaon Beta Cut
		      
		    } // end PID selection cut on SHMS (Kaons) and HMS (electrons)
		    
		    
		    //==============================================================================================================================
		    
		    //=====================================
		    //  PION BEAM-SPIN ASYMMETRY ANALYSIS
		    //=====================================
		    
		    //Fill Beta vs. Coin Time (Assuming Pion Mass in Coin Time calculation)
		    H_Beta_vs_ctime_Pi_noPID->Fill((ePiCoinTime-Pi_ctime_offset), phod_gtr_beta);
		    
		    // C.Y. NEED TO: Make Pion selection PID CUT HERE ! (BASED ON THE 2D PID CORRELATIONS FILLED ABOVE)
		    if(pion_pid && elec_pid) {		    
		      
		      //C.Y. NEED TO: Make a Beta Cut HERE! (Based on the 1D and 2D Beta plots)
		      if(pion_beta_cut) {
			
			
			//Fill Beta (Pion PID)
			H_Beta_Pi_PID->Fill(phod_gtr_beta);
			H_Beta_vs_ctime_Pi_PID->Fill((ePiCoinTime-Pi_ctime_offset), phod_gtr_beta);
			
			//Missing Mass (Pion PID)
			H_MM_Pi_PID->Fill(MM_Pi);
			
			//------------------------------------------
			//           TRUE PION SELECTION
			//------------------------------------------		    
			
			//C.Y. NEED TO: Make coin. time cut to select TRUE ePi coin. time peak HERE !
			if(ePi_ctime_cut) {		       
			  
			  H_ePi_ctime_real->Fill(ePiCoinTime-Pi_ctime_offset);
			  H_MM_Pi_real->Fill(MM_Pi); 
			  
			  //C.Y. NEED TO: Make a Pion  Missing Mass Cut HERE !
			  if(pion_MM_cut){			    		    
			    
			    //Fill Beta (with missing mass cut)
			    H_Beta_Pi_real->Fill(phod_gtr_beta);
			    
			    //Check if event has "+" helicity
			    if(hel_sign>0){
			      //Fill Fill In-plane th_xq_cm "VS" Out-of-plane angle, ph_qx (pos helicity) for the Pions 
			      H_thxqCM_vs_phxq_Pi_pos->Fill(ph_xq/dtr, th_xq_cm/dtr);   //in deg			    
			    }
			    
			    //Check if event has "-" helicity
			    else if(hel_sign<0){
			      //Fill Fill In-plane th_xq_cm "VS" Out-of-plane angle, ph_qx (neg helicity) for the Pions
			      H_thxqCM_vs_phxq_Pi_neg->Fill(ph_xq/dtr, th_xq_cm/dtr);   //in deg			    
			    }
			    
			  } // end missing mass cut around neutron (939.565), to select reaction channel, H(e,e'Pi+)MM, where MM=neutron
			  
			} // end ePi real coin. selection cut
			
			
			
			//------------------------------------------
			//         RANDOM BKG PION SELECTION
			//------------------------------------------
			
			//C.Y. NEED TO: Make coin. time cut to select randoms (sample outside main peak) ePi coin. time peak HERE !
			if(ePi_ctime_cut_rand) {		       
			  
			  H_ePi_ctime_rand->Fill(ePiCoinTime-Pi_ctime_offset);			  
			  H_MM_Pi_rand->Fill(MM_Pi);			
			  
			  //C.Y. NEED TO: Make a Pion  Missing Mass Cut HERE !
			  if(pion_MM_cut){			 
			    
			    //Fill Beta (with missing mass cut)
			    H_Beta_Pi_rand->Fill(phod_gtr_beta);
			    
			    //Check if event has "+" helicity
			    if(hel_sign>0){
			      //Fill Fill In-plane th_xq_cm "VS" Out-of-plane angle, ph_qx (pos helicity) for the Pions
			      H_thxqCM_vs_phxq_Pi_pos_rand->Fill(ph_xq/dtr, th_xq_cm/dtr);  //in deg			    
			    }
			    
			    //Check if event has "-" helicity
			    else if(hel_sign<0){
			      //Fill Fill In-plane th_xq_cm "VS" Out-of-plane angle, ph_qx (neg helicity) for the Pions
			      H_thxqCM_vs_phxq_Pi_neg_rand->Fill(ph_xq/dtr, th_xq_cm/dtr);  //in deg			    
			    }
			    
			  }// end missing mass cut around neutron (939.565), to select reaction channel, H(e,e'Pi+)MM, where MM=Pion
			  
			}// end ePi random coin. selection
			
		      } // end Pion Beta Cut
		      
		    }// end PID selection cut on SHMS (Pions) and HMS (electrons)
		    
		    //===================================================================================================================
		    
		    //=======================================
		    //  PROTON BEAM-SPIN ASYMMETRY ANALYSIS
		    //=======================================
		    
		    //Fill Beta vs. Coin Time (Assuming Proton Mass in Coin Time calculation)
		    H_Beta_vs_ctime_P_noPID->Fill((epCoinTime-P_ctime_offset), phod_gtr_beta);
		    
		    // C.Y. NEED TO: Make Proton selection PID CUT HERE ! (BASED ON THE 2D PID CORRELATIONS FILLED ABOVE)
		    if(proton_pid && elec_pid) {		    
		      
		      //C.Y. NEED TO: Make a Beta Cut HERE! (Based on the 1D and 2D Beta plots)
		      if(proton_beta_cut) {		      
			
			//Fill Beta (Proton PID)
			H_Beta_P_PID->Fill(phod_gtr_beta);
			H_Beta_vs_ctime_P_PID->Fill((epCoinTime-P_ctime_offset), phod_gtr_beta);
			
			//Missing Mass (Proton PID)
			H_MM_P_PID->Fill(MM_P);
			
			//------------------------------------------
			//           TRUE PROTON SELECTION
			//------------------------------------------		    
			
			//C.Y. NEED TO: Make coin. time cut to select TRUE eP coin. time peak HERE !
			if(eP_ctime_cut) {
			  
			  H_ep_ctime_real->Fill(epCoinTime-P_ctime_offset);			  
			  H_MM_P_real->Fill(MM_P); 		       
			  
			  //C.Y. NEED TO: Make a Proton  Missing Mass Cut HERE !
			  if(proton_MM_cut){			  			    
			    
			    //Fill Beta (with missing mass cut)
			    H_Beta_P_real->Fill(phod_gtr_beta);
			    
			    //Check if event has "+" helicity
			    if(hel_sign>0){
			      //Fill Fill In-plane th_xq_cm "VS" Out-of-plane angle, ph_qx (pos helicity) for the Protons 
			      H_thxqCM_vs_phxq_P_pos->Fill(ph_xq/dtr, th_xq_cm/dtr);   //in deg			    
			    }
			    
			    //Check if event has "-" helicity
			    else if(hel_sign<0){
			      //Fill Fill In-plane th_xq_cm "VS" Out-of-plane angle, ph_qx (neg helicity) for the Protons
			      H_thxqCM_vs_phxq_P_neg->Fill(ph_xq/dtr, th_xq_cm/dtr);   //in deg			      
			    }
			    
			  } // end missing mass cut around "nothing ... ??" (0), to select reaction channel, H(e,e'p)MM, where MM = 0 [elastic?]
			  
			} // end eP real coin. selection cut						
			
			//------------------------------------------
			//         RANDOM BKG PROTON SELECTION
			//------------------------------------------
			
			//C.Y. NEED TO: Make coin. time cut to select randoms (sample outside main peak) eP coin. time peak HERE !
			if(eP_ctime_cut_rand) {		       
			  
			  H_ep_ctime_rand->Fill(epCoinTime-P_ctime_offset);			  
			  H_MM_P_rand->Fill(MM_P);		      
			  
			  //C.Y. NEED TO: Make a Proton  Missing Mass Cut HERE !
			  if(proton_MM_cut){			  
			    
			    //Fill Beta (with missing mass cut)
			    H_Beta_P_rand->Fill(phod_gtr_beta);
			    
			    //Check if event has "+" helicity
			    if(hel_sign>0){
			      //Fill Fill In-plane th_xq_cm "VS" Out-of-plane angle, ph_qx (pos helicity) for the Protons
			      H_thxqCM_vs_phxq_P_pos_rand->Fill(ph_xq/dtr, th_xq_cm/dtr);  //in deg			      
			    }
			    
			    //Check if event has "-" helicity
			    else if(hel_sign<0){
			      //Fill Fill In-plane th_xq_cm "VS" Out-of-plane angle, ph_qx (neg helicity) for the Protons
			      H_thxqCM_vs_phxq_P_neg_rand->Fill(ph_xq/dtr, th_xq_cm/dtr);  //in deg			    
			    }
			    
			  }// end missing mass cut around "nothing . . .?" (0), to select reaction channel, H(e,e'p)MM, where MM = 0 [elastic?]
			  
			}// end eP random coin. selection
			
		      } // end Proton Beta Cut
		      
		    }// end PID selection cut on SHMS (Protons) and HMS (electrons)

		  } // end accpetance cuts  
		  
		  //===================================================================================================================
		  
		  
		}  //------END: REQUIRE "NO EDTM" CUT TO FILL DATA HISTOGRAMS-----
	      
	    }  //-----END: BCM Current Cut------
	  
	  //Increment Scaler Read if event == scaler_evt_perlimit for that scaler read (See baseAnalyzer::EventLoop() for detailed explanation)
	  
	  gevnum = ientry + 1;
	  
	  if(gevnum==scal_evt_num[scal_read]){ scal_read++; }
	  
	  // print every 100000 events
	  if (ientry % 100000 == 0){
	    cout << "Helicity DataEventLoop: " << ientry << "/" << nentries << "( " << std::setprecision(2) << double(ientry) / nentries * 100. << "  % )" << std::flush << "\r";
	  }
	} //END DATA EVENT LOOP      
      
    }//END DATA ANALYSIS
  
  if(analysis=="simc")
    {
      cout << "SIMC ANALYSIS needs to be done for Kaon LT Analysis (2018) . . . " << endl;
    }
  
} //End helicityAnalyzer::EventLoop()

//_______________________________________________________________________________
void helicityAnalyzer::CalcEff()
{
  cout << "Calling Derived CalcEff() . . . " << endl;
  
  //As of now, there is nothing to add to the existing method from baseAnalyzer.
  
  //Call the existing method, as it is needed to calculate live time / tracking eff.
  // as well as convert rates to kHz and charge to mC.
  baseAnalyzer::CalcEff();
  
  //HERE: Add additional calculations if needed.
  
  //--Calculate HMS Tracking Efficiency--                                                                                                                 

  //HMS + helicity
  hTrkEff_pos = h_did_pos / h_should_pos;                                                                                                                  
  hTrkEff_err_pos = sqrt(h_should_pos-h_did_pos) / h_should_pos;

  //HMS - helicity
  hTrkEff_neg = h_did_neg / h_should_neg;                                                                                                                  
  hTrkEff_err_neg = sqrt(h_should_neg-h_did_neg) / h_should_neg;

  //--Calculate SHMS Tracking Efficiency--

  //SHMS + helicity
  pTrkEff_pos = p_did_pos / p_should_pos; 
  pTrkEff_err_pos = sqrt(p_should_pos-p_did_pos) / p_should_pos;

  //SHMS - helicity
  pTrkEff_neg = p_did_neg / p_should_neg; 
  pTrkEff_err_neg = sqrt(p_should_neg-p_did_neg) / p_should_neg;
}

//_______________________________________________________________________________
void helicityAnalyzer::ApplyWeight()
{
  cout << "Calling Derived ApplyWeight() . . . " << endl;


  /*

  //baseAnalyzer::ApplyWeight();
  
  For now, do not call baseAnalyzer::ApplyWeight, as it defined certain weight factors that are different
  for each experiment, and are used to scale the histos. Need to Scale the histos separately
  from determining the weight factors, so that the method of scaling histos may be re-used
  more effectively. Also, since the derived EventLoop() method decides which histograms to 
  fill, it may be better to determine which histograms to fill and scale by the FullWeight in
  the derived method as not all histograms in the baseAnalyzer may be used. Below, we Scale
  only histograms of interest for the helicity analysis, as well as histos of interest from the
  baseAnalyzer class.
  
  NOTE: It is not clear yet, which asymmetry quantities need to be scaled, as most correction factors
  will cancel out when the ratio is taken.  Asy = ( N+/eff+ - N-/eff- ) / ( N+/eff+ + N-/eff-), where
  eff = total_charge * tracking_efficiency * total_live_tine_efficiency * tgt_boiling * hadron_absorption,
  however, when one separates the total counts N into  N+ and N-, which efficiencies CANCEL OUT, and which
  DO NOT when the ratio is taken?  For now, apply a scaling of 1 (no scaling).

  If scaling is ever done for the pos/neg helicity histograms, separate correction factors for pos/neg
  helicity need to be calculated. For example, for tracking efficiency, one would need determine the
  efficiency for pos/neg helicity separately. Then, one can determine how different are these correction factors.
  If they are not different, then they will cancel out when the ratio is taken.
  */
  
  //HERE: Re-define which weight factors to apply, depending on the specific experiment

  tgtBoil_corr = 1.;  //Assume no target boiling for now
  tgtBoil_corr_err = -1.; 
  hadAbs_corr = 1.;   //Assume no hadron absorption for now
  hadAbs_corr_err = -1.;

  //Modify Full Weight accordingly as Kaon LT group gets results for certain corrections
  //FullWeight = 1. / (total_charge_bcm_cut * hTrkEff * pTrkEff * tLT_trig * tgtBoil_corr * hadAbs_corr);

  //For testing purposes, only normalize by total charge.
  //FullWeight = 1. / total_charge_bcm_cut;

  //Do not apply any scaling for now. 
  FullWeight = 1.;
  
  //Scale Data Histograms by Full Weight (Each run for a particular kinematics can then be combined, once they are scaled by the FullWeight)
  
  //----SCALE HISTOGRAMS BY LOOPING OVER LISTS----
  
  //determine what class types are in the list
  TString class_name;
  
  
  //------------------------------------------
  //Lopp over hel_HList of histogram objects
  //------------------------------------------
  for(int i=0; i<hel_HList->GetEntries(); i++) {
    //Get the class name for each element on the list (either "TH1F" or TH2F")
    class_name = hel_HList->At(i)->ClassName();
    //Read ith histograms in the list from current run
    if(class_name=="TH1F") {
      //Get and scale histogram from the list
      h_i = (TH1F *)hel_HList->At(i); h_i->Scale(FullWeight); 
    }
    if(class_name=="TH2F") {
      //Get and scale histogram from the list
      h2_i = (TH2F *)hel_HList->At(i); h2_i->Scale(FullWeight);
    }   
  } //end loop of hel_HList
  
  //-----------------------------------------------------
  //Lopp over pid_HList of histogram objects (inherited)
  //----------------------------------------------------
  for(int i=0; i<pid_HList->GetEntries(); i++) {
    //Get the class name for each element on the list (either "TH1F" or TH2F")
    class_name = pid_HList->At(i)->ClassName();
    //Read ith histograms in the list from current run
    if(class_name=="TH1F") {
      //Get and scale histogram from the list
      h_i = (TH1F *)pid_HList->At(i); h_i->Scale(FullWeight); 
    }
    if(class_name=="TH2F") {
      //Get and scale histogram from the list
      h2_i = (TH2F *)pid_HList->At(i); h2_i->Scale(FullWeight);
    }   
  }//end loop over pid_HList
  
  //-----------------------------------------------------
  //Lopp over kin_HList of histogram objects (inherited)
  //----------------------------------------------------
  for(int i=0; i<kin_HList->GetEntries(); i++) {
    //Get the class name for each element on the list (either "TH1F" or TH2F")
    class_name = kin_HList->At(i)->ClassName();
    //Read ith histograms in the list from current run
    if(class_name=="TH1F") {
      //Get and scale histogram from the list
      h_i = (TH1F *)kin_HList->At(i); h_i->Scale(FullWeight); 
    }
    if(class_name=="TH2F") {
      //Get and scale histogram from the list
      h2_i = (TH2F *)kin_HList->At(i); h2_i->Scale(FullWeight);
    }   
  }//end loop over kin_HList	
  
  //-----------------------------------------------------
  //Lopp over accp_HList of histogram objects (inherited)
  //----------------------------------------------------
  for(int i=0; i<accp_HList->GetEntries(); i++) {
    //Get the class name for each element on the list (either "TH1F" or TH2F")
    class_name = accp_HList->At(i)->ClassName();
    //Read ith histograms in the list from current run
    if(class_name=="TH1F") {
      //Get and scale histogram from the list
      h_i = (TH1F *)accp_HList->At(i); h_i->Scale(FullWeight); 
    }
    if(class_name=="TH2F") {
      //Get and scale histogram from the list
      h2_i = (TH2F *)accp_HList->At(i); h2_i->Scale(FullWeight);
    }   
  }//end loop over accp_HList
  
  //Call the randoms subtraction methods (after scaling all histograms above)
  RandSub();
  
}

//_______________________________________________________________________________
void helicityAnalyzer::RandSub()
{
  cout << "Calling RandSub() " << endl;

  /*
    Brief: This methods carries out the subtraction of random coincidences (outside coin peak selection) 
    from real coincidences (within coin peak selected) for various histograms 
  */

  // Scale Down (If necessary) the randoms before subtracting it from the reals
  // NOTE: As long as cpid_eKctime_max = eK_mult * eK_ctime_thrs, where eK_mult is a integral multiple (i.e., 2, 3, 4,  . . . ) -- same applies for (Pions, Protons)
  // At least a multiple of 2 is required, otherwise, if multiple =1, min/max cuts are the same which gives 0 scale factor --> infinity 
  // the random coincidences can be scaled down by the following factor:
  
  K_scale_factor =  ((eK_mult*eK_ctime_thrs) - eK_ctime_thrs)     / eK_ctime_thrs;
  Pi_scale_factor = ((ePi_mult*ePi_ctime_thrs) - ePi_ctime_thrs)  / ePi_ctime_thrs;
  P_scale_factor =  ((eP_mult*eP_ctime_thrs) - eP_ctime_thrs)     / eP_ctime_thrs;


  cout << "K_scale_factor = "         << K_scale_factor << endl;
  cout << "Pi_scale_factor = "        << Pi_scale_factor << endl;
  cout << "P_scale_factor = "         << P_scale_factor << endl;
  
  //----Scale Down the random coincidences histograms-----
  // ----(other than the coin. histograms themselves)----

  //Missing Mass Randoms
  H_MM_K_rand->Scale(1./K_scale_factor);
  H_MM_Pi_rand->Scale(1./Pi_scale_factor);
  H_MM_P_rand->Scale(1./P_scale_factor);

  H_Beta_K_rand->Scale(1./K_scale_factor);
  H_Beta_Pi_rand->Scale(1./Pi_scale_factor);
  H_Beta_P_rand->Scale(1./P_scale_factor);
  
  //positive helicity randoms
  H_thxqCM_vs_phxq_K_pos_rand->Scale(1./K_scale_factor);
  H_thxqCM_vs_phxq_Pi_pos_rand->Scale(1./Pi_scale_factor);
  H_thxqCM_vs_phxq_P_pos_rand->Scale(1./P_scale_factor);

  //negative helicity randoms
  H_thxqCM_vs_phxq_K_neg_rand->Scale(1./K_scale_factor);
  H_thxqCM_vs_phxq_Pi_neg_rand->Scale(1./Pi_scale_factor);
  H_thxqCM_vs_phxq_P_neg_rand->Scale(1./P_scale_factor);
  
  // -----Carry out the randoms subtraction------

  //Missing Mass Randoms Subtraction
  H_MM_K_rand_sub->Add(H_MM_K_real, H_MM_K_rand, 1, -1);
  H_MM_Pi_rand_sub->Add(H_MM_Pi_real, H_MM_Pi_rand, 1, -1);
  H_MM_P_rand_sub->Add(H_MM_P_real, H_MM_P_rand, 1, -1);

  H_Beta_K_rand_sub->Add(H_Beta_K_real, H_Beta_K_rand, 1, -1);
  H_Beta_Pi_rand_sub->Add(H_Beta_Pi_real, H_Beta_Pi_rand, 1, -1);
  H_Beta_P_rand_sub->Add(H_Beta_P_real, H_Beta_P_rand, 1, -1);
  
  //positive helicity randoms subtraction
  H_thxqCM_vs_phxq_K_pos_rand_sub->Add(H_thxqCM_vs_phxq_K_pos, H_thxqCM_vs_phxq_K_pos_rand, 1, -1);
  H_thxqCM_vs_phxq_Pi_pos_rand_sub->Add(H_thxqCM_vs_phxq_Pi_pos, H_thxqCM_vs_phxq_Pi_pos_rand, 1, -1);
  H_thxqCM_vs_phxq_P_pos_rand_sub->Add(H_thxqCM_vs_phxq_P_pos, H_thxqCM_vs_phxq_P_pos_rand, 1, -1);
  
  //negative helicity randoms subtraction
  H_thxqCM_vs_phxq_K_neg_rand_sub->Add(H_thxqCM_vs_phxq_K_neg, H_thxqCM_vs_phxq_K_neg_rand, 1, -1);
  H_thxqCM_vs_phxq_Pi_neg_rand_sub->Add(H_thxqCM_vs_phxq_Pi_neg, H_thxqCM_vs_phxq_Pi_neg_rand, 1, -1);
  H_thxqCM_vs_phxq_P_neg_rand_sub->Add(H_thxqCM_vs_phxq_P_neg, H_thxqCM_vs_phxq_P_neg_rand, 1, -1);
  
}

//_______________________________________________________________________________
void helicityAnalyzer::WriteHist()
{

  //Call baseAnalyzer::WriteHist() method to write the inherited TLists (with histos) onto a ROOTfile
  //NOTE: All histograms added in the baseAnalyzer TLists will appear, however, if the histos
  //have NOT been filled in the helicittAnalyzer::EventLoop() method, the histos will be empty.
  baseAnalyzer::WriteHist();
  
  cout << "Calling Derived WriteHist() . . ." << endl;

  //HERE: Make additions of TLists containing histograms defined in this class as needed.
  
  //Write Data Histograms
  if(analysis=="data")
    {
      //Update Output ROOTfile (Since baseAnalyzer::WriteHist() was called firts, it had already
      //created a output root file per run, therefore, it should be open in "UPDATE" mode. 
      outROOT = new TFile(data_OutputFileName, "UPDATE");

      //Make directories to store histograms based on category
      outROOT->mkdir("hel_plots");
      
      //Write HELICITY histos to hel_plots directory
      outROOT->cd("hel_plots");
      hel_HList->Write();

      
      //Close File
      outROOT->Close();
    }

  
}

//_______________________________________________________________________________
void helicityAnalyzer::WriteReport()
{
  
  /*
    Method to write charge, efficiencies, live time and other relevant quantities to a data file
    NOTE: This method could have been inherited from the baseAnalyzer class, howeever, in this
    derived class, we do not yet have set up a system of cuts to apply from the cuts input file.
    Although we do inherit the cuts defined in the input cuts file, we are still in the process of
    figuring out how the cuts will be applied in this class, as we may need to select reals/background
    during a single run, for the purpose of background subtraction.
    
    NOTE: We may need to add additional variables specific to pos/neg helicities, for example tracking eff. for +/- hel.
    to check whether the efficiencies are the same (or at least similar, within statistics), that they would cancel out
    in the asymmetry calculation.
  */

  
  cout << "Calling WriteReport() . . ." << endl;
  
  if(analysis=="data"){

    //---------------------------------------------------------

    //Choose which trigger is relevant
    
    //Check if file already exists
    in_file.open(report_OutputFileName.Data());
    
    if(in_file.fail()){
      
      cout << "Report File does NOT exist, will create one . . . " << endl;
      
      out_file.open(report_OutputFileName);
      out_file << "#-------------------------------------" << endl;
      out_file << "#        Data Analysis Summary        " << endl;
      out_file << "#-------------------------------------" << endl;
      out_file << "#                                     " << endl;
      out_file << Form("# %s  | Beam Current Threshold: > %.2f uA ", bcm_type.Data(), current_thrs_bcm) << endl;
      out_file << "#                                     " << endl;
      out_file << Form("# DAQ Mode: %s | Trigger: %s              ", daq_mode.Data(), trig_type.Data()) << endl;
      out_file << Form("# electron arm: %s                        ", e_arm_name.Data() ) << endl;
      out_file << "#                                              " << endl;
      out_file << "#---Acceptance Cuts--- " << endl;
      if(hdelta_cut_flag)   {out_file << Form("# SHMS Momentum Acceptance (P.gtr.dp): [%.3f, %.3f] %%", hdelta_min, hdelta_max ) << endl;}
      if(edelta_cut_flag)   {out_file << Form("# HMS  Momentum Acceptance (H.gtr.dp): [%.3f, %.3f] %%", edelta_min, edelta_max ) << endl;}
      if(ztarDiff_cut_flag) {out_file << Form("# Z-Target Vertex Difference (P.react.z-H.react.z) : [%.3f, %.3f)", ztarDiff_min, ztarDiff_max) << endl;}      
      out_file << "#                                     " << endl;
      out_file << "# ---Particle Identification (PID) Cuts---" << endl;
      out_file << "# NOTE: Coincidence time background selection is done by averaging a sample (left) and (right) of the main coincidence peak, where the sample is assumed to have the same structure on both sides." << endl;
      out_file << "# The sample background is then scaled to the background underneath the main coin peak, and subtracted." << endl;
      out_file << "#                                              " << endl;
      out_file << "# --> electrons (HMS): " << endl;
      if(hcer_pidCut_flag)           {out_file << Form("# HMS Cherenkov NPE Sum (H.cer.npeSum) >= %.3f", elec_hcer_npe_thrs) << endl; }
      if(hetot_trkNorm_pidCut_flag)  {out_file << Form("# HMS Calorimeter (H.cal.etottracknorm) Etot/Etrak >= %.3f", elec_hcal_thrs) << endl; }
      out_file << "#                                              " << endl;
      out_file << "# --> Kaons     (SHMS): " << endl;
      if(K_paero_npe_flag)   {out_file << Form("# SHMS Aerogel NPE Sum   (P.aero.npeSum)  >= %.3f", K_paero_npe_thrs) << endl; }
      if(K_phgcer_npe_flag)  {out_file << Form("# SHMS Heavy Gas NPE Sum (P.hgcer.npeSum) <= %.3f", K_phgcer_npe_thrs) << endl; }
      if(K_beta_flag)        {out_file << Form("# SHMS Hodoscope Beta    (P.gtr.beta):  [%.3f, %.3f]", K_beta_min, K_beta_max) << endl; }
      if(eK_ctime_flag)
	{
	  out_file << Form("# eK Coincidence Time (CTime.eKCoinTime_ROC2): abs( eK_coin_peak - offset(%.3f) ) <= %.3f ns", K_ctime_offset, eK_ctime_thrs) << endl;
	  out_file << Form("# eK Coincidence Time Background: abs( eK_coin_peak - offset(%.3f) ) > %.3f && abs( eK_coin_peak - offset(%.3f) ) <= (%.3f x %.3f) ns", K_ctime_offset, eK_ctime_thrs, K_ctime_offset, eK_mult, eK_ctime_thrs) << endl;
	}
      if(MM_K_cut_flag)  {out_file << Form("# Kaon Missing Mass: [%.3f, %.3f] GeV/c^2", MM_K_min, MM_K_max) << endl; }
      out_file << "#                                              " << endl;
      out_file << "# --> Pions     (SHMS): " << endl;
      if(Pi_phgcer_npe_flag)  {out_file << Form("# SHMS Heavy Gas NPE Sum (P.hgcer.npeSum) >= %.3f", Pi_phgcer_npe_thrs) << endl; }
      if(Pi_beta_flag)        {out_file << Form("# SHMS Hodoscope Beta    (P.gtr.beta):  [%.3f, %.3f]", Pi_beta_min, Pi_beta_max) << endl; }
      if(ePi_ctime_flag)
	{
	  out_file << Form("# ePi Coincidence Time (CTime.ePiCoinTime_ROC2): abs( ePi_coin_peak - offset(%.3f) ) <= %.3f ns", Pi_ctime_offset, ePi_ctime_thrs) << endl;
	  out_file << Form("# ePi Coincidence Time Background: abs( ePi_coin_peak - offset(%.3f) ) > %.3f && abs( ePi_coin_peak - offset(%.3f) ) <= (%.3f x %.3f) ns", Pi_ctime_offset, ePi_ctime_thrs, Pi_ctime_offset, ePi_mult, ePi_ctime_thrs) << endl;
	}
      if(MM_Pi_cut_flag)  {out_file << Form("# Pion Missing Mass: [%.3f, %.3f] GeV/c^2", MM_Pi_min, MM_Pi_max) << endl; }
      out_file << "#                                              " << endl;
      out_file << "# --> Protons   (SHMS): " << endl;
      if(P_paero_npe_flag)    {out_file << Form("# SHMS Aerogel NPE Sum   (P.aero.npeSum)  <= %.3f", P_paero_npe_thrs) << endl; }
      if(P_phgcer_npe_flag)   {out_file << Form("# SHMS Heavy Gas NPE Sum (P.hgcer.npeSum) <= %.3f", P_phgcer_npe_thrs) << endl; }
      if(P_beta_flag)         {out_file << Form("# SHMS Hodoscope Beta    (P.gtr.beta):  [%.3f, %.3f]", P_beta_min, P_beta_max) << endl; }
      if(eP_ctime_flag)
	{
	  out_file << Form("# eP Coincidence Time (CTime.epCoinTime_ROC2): abs( eP_coin_peak - offset(%.3f) ) <= %.3f ns", P_ctime_offset, eP_ctime_thrs) << endl;
	  out_file << Form("# eP Coincidence Time Background: abs( eP_coin_peak - offset(%.3f) ) > %.3f && abs( eP_coin_peak - offset(%.3f) ) <= (%.3f x %.3f) ns", P_ctime_offset, eP_ctime_thrs, P_ctime_offset, eP_mult, eP_ctime_thrs) << endl;
	}
      if(MM_P_cut_flag)  {out_file << Form("# Proton Missing Mass: [%.3f, %.3f] GeV/c^2", MM_P_min, MM_P_max) << endl; }
      
      out_file << "#                       " << endl;
      out_file << "# Units: charge [mC] | currnet [uA] | rates [kHz] |  efficiencies [fractional form]                       " << endl;
      out_file << "#                       " << endl;
      out_file << std::setw(2) << "#! Run[i,0]/" << std::setw(25) << "charge[f,1]/" << std::setw(25) << "avg_current[f,2]/" << std::setw(25)  << "hTrkEff[f,3]/" << std::setw(25) << "hTrkEff_err[f,4]/" << std::setw(25) << "pTrkEff[f,5]/" << std::setw(25) << "pTrkEff_err[f,6]/" << std::setw(25) << "tgt_boil_factor[f,7]/" << std::setw(30) << "tgt_boil_factor_err[f,8]/" << std::setw(25) << "hadAbs_factor[f,9]/" << std::setw(30) << "hadAbs_factor_err[f,10]/" << std::setw(25) <<  "cpuLT[f,11]/" << std::setw(25) << "cpuLT_err_Bi[f,12]/" << std::setw(25) << "cpuLT_err_Bay[f,13]/" << std::setw(25) << "tLT[f,14]/" << std::setw(25) << "tLT_err_Bi[f,15]/" << std::setw(25) << "tLT_err_Bay[f,16]/" << std::setw(25) <<"S1X_rate[f,17]/" << std::setw(30) << "shms_Ps1_3of4_rate[f,18]/" << std::setw(30) << "hms_Ps3_elreal_rate[f,19]/"<< std::setw(25) << "coin_Ps5_rate[f,18]/"  << std::setw(25) << "edtm_rate[f,19]/"  << std::setw(25) << "Ps1_factor[f,20]/" << std::setw(25) << "Ps3_factor[f,21]/" <<std::setw(25) << "Ps5_factor[f,22]/" <<  std::setw(25) << "hTrkEff_pos[f,23]/" <<  std::setw(25) << "hTrkEff_err_pos[f,24]/" << std::setw(25) << "hTrkEff_neg[f,25]/" <<  std::setw(25) << "hTrkEff_err_neg[f,26]/" <<  std::setw(25) << "pTrkEff_pos[f,27]/" <<  std::setw(25) << "pTrkEff_err_pos[f,28]/" << std::setw(25) << "pTrkEff_neg[f,29]/" <<  std::setw(25) << "pTrkEff_err_neg[f,30]/" << std::setw(25) << "edtm_accp[f,31]/" << std::setw(25) << "edtm_scaler[f,32]/" << std::setw(25) << "trig5_accp[f,33]/" << std::setw(25) << "trig5_scaler[f,34]/" << endl;

      out_file.close();
      in_file.close();

    }

    //Open Report FIle in append mode
    out_file.open(report_OutputFileName, ios::out | ios::app);
    out_file << std::setw(7) << run  << std::setw(25) << total_charge_bcm_cut << std::setw(25) << avg_current_bcm_cut << std::setw(25) << hTrkEff << std::setw(25) << hTrkEff_err << std::setw(25) << pTrkEff << std::setw(25) << pTrkEff_err << std::setw(25) << tgtBoil_corr << std::setw(25) << tgtBoil_corr_err << std::setw(25) << hadAbs_corr << std::setw(25) << hadAbs_corr_err << std::setw(25) << cpuLT_trig << std::setw(25) << cpuLT_trig_err_Bi << std::setw(25) << cpuLT_trig_err_Bay << std::setw(25) << tLT_trig << std::setw(25) << tLT_trig_err_Bi << std::setw(25) << tLT_trig_err_Bay << std::setw(25) << S1XscalerRate_bcm_cut << std::setw(30) << TRIG1scalerRate_bcm_cut << std::setw(30) << TRIG3scalerRate_bcm_cut << std::setw(25) << TRIG5scalerRate_bcm_cut << std::setw(25) << EDTMscalerRate_bcm_cut << std::setw(25) << Ps1_factor << std::setw(25) << Ps3_factor << std::setw(25) << Ps5_factor << std::setw(25) << hTrkEff_pos << std::setw(25) << hTrkEff_err_pos << std::setw(25) << hTrkEff_neg << std::setw(25) << hTrkEff_err_neg << std::setw(25) << pTrkEff_pos << std::setw(25) << pTrkEff_err_pos << std::setw(25) << pTrkEff_neg << std::setw(25) << pTrkEff_err_neg << std::setw(25) << total_edtm_accp_bcm_cut << std::setw(25) << (total_edtm_scaler_bcm_cut / Ps_factor) << std::setw(25) << total_trig_accp_bcm_cut << std::setw(25) << (total_trig_scaler_bcm_cut / Ps_factor) << endl;
    out_file.close();
    
  }

} //End WriteReport()

//_______________________________________________________________________________
void helicityAnalyzer::CombineHistos()
{

  cout << "Calling Derived CombineHistos() " << endl;

  if(combine_histos==0) return;
   
  //Call baseAnalyzer::CombineHistos() method to combine multiple runs from the
  //same kinematics for all histograms defined in the baseAnalyzer.
  baseAnalyzer::CombineHistos();
 
  //Below, add code to combine histograms defined in this helicityAnalyzer class.
  
  //Open existing *_final.root file in UPDATE mode 
  outROOT = new TFile(data_OutputFileName_combined, "UPDATE");

  //check if directory to store helicity plots exits, otherwise, create it and add the histo list from the 1st run
  if(outROOT->GetDirectory("hel_plots")==0)
    {

      cout << "hel_plots direcotry does NOT exits ! Will create one" << endl;
      //If histo directory does not exist, make directories to store histograms based on category
      outROOT->mkdir("hel_plots");
      outROOT->cd("hel_plots");
      
      hel_HList->Write();      
      outROOT->Close();
    }
  
 

  //If helicity histo directory exists, keep adding histo counts for each run
  else 
    {      
      //Open the existing *_final.root file created in the baseAnalyzer:CombineHistos()
      outROOT->ReOpen("READ"); 

      //determine what class types are in the list
      TString class_name;
      
      //Set up histogram names/locations to find
      TString hist_name;
      TString hist_dir;
      
      //------------------------------------------
      //Lopp over hel_HList of histogram objects
      //------------------------------------------
      for(int i=0; i<hel_HList->GetEntries(); i++)
	{	  
	  //Determine object data type (as of now, either TH1F or TH2F are possible)
	  class_name = hel_HList->At(i)->ClassName();
	  
	  //Read ith histograms in the list from current run

	  if(class_name=="TH1F") {
	    //Get histogram from current run
	    h_i = (TH1F *)hel_HList->At(i); 
	    hist_name = h_i->GetName();
	    //full path to histogram must be given, otherwise, histogram will NOT be retreived from ROOTfile
	    hist_dir = "hel_plots/" + hist_name;  
	    //Get cumulative histogram object stored in ROOTfile and add the h_i (histogram from ith run) to h_total (cumulative histogram) 	    
	    outROOT->GetObject(hist_dir, h_total); h_total->Add(h_i); outROOT->ReOpen("UPDATE"); outROOT->cd("hel_plots"); h_total->Write("", TObject::kOverwrite); outROOT->ReOpen("READ");	    

	  }
	  
	  if(class_name=="TH2F") {
	    //Get histogram from current run
	    h2_i = (TH2F *)hel_HList->At(i); 
	    hist_name = h2_i->GetName();
	    //full path to histogram must be given, otherwise, histogram will NOT be retreived from ROOTfile
	    hist_dir = "hel_plots/" + hist_name;  
	    //Get cumulative histogram object stored in ROOTfile and add the h_i (histogram from ith run) to h_total (cumulative histogram) 
	    outROOT->GetObject(hist_dir, h2_total); h2_total->Add(h2_i); outROOT->ReOpen("UPDATE"); outROOT->cd("hel_plots"); h2_total->Write("", TObject::kOverwrite); outROOT->ReOpen("READ");	       
	  
	  }
	  
	}//end loop over hel_HList
            
      outROOT->Close();
      
    }//end else statement
  
}

//_______________________________________________________________________________
void helicityAnalyzer::BinExtraction()
{

  cout << "Calling BinExtraction() " << endl;
  /*
    Brief: This method extracts the bin information from 2D Histograms.
    This is often a very useful application towards the end of an analysis,
    where the user needs to extract the binning information of the cross-section,
    asymmetry, etc. onto a .txt file for further analysis and plotting. 
  */

  //Call the histogram utility function to extract information (see header 'utils/hist_utils.h' in baseAnalyzer.h)

  //Extract Bin Information Per Run
  extract_2d_hist(H_thxqCM_vs_phxq_K_pos_rand_sub,  "phi_xq [deg]", "th_xq_CM [deg]", Form("OUTPUT/pos_hel_Kaon_%d.txt", run));
  extract_2d_hist(H_thxqCM_vs_phxq_Pi_pos_rand_sub, "phi_xq [deg]", "th_xq_CM [deg]", Form("OUTPUT/pos_hel_Pion_%d.txt", run));
  //extract_2d_hist(H_thxqCM_vs_phxq_P_pos_rand_sub,  "phi_xq [deg]", "th_xq_CM [deg]", Form("pos_hel_Proton_%d.txt", run));
   
  extract_2d_hist(H_thxqCM_vs_phxq_K_neg_rand_sub,  "phi_xq [deg]", "th_xq_CM [deg]", Form("OUTPUT/neg_hel_Kaon_%d.txt", run));
  extract_2d_hist(H_thxqCM_vs_phxq_Pi_neg_rand_sub, "phi_xq [deg]", "th_xq_CM [deg]", Form("OUTPUT/neg_hel_Pion_%d.txt", run));
  //extract_2d_hist(H_thxqCM_vs_phxq_P_neg_rand_sub,  "phi_xq [deg]", "th_xq_CM [deg]", Form("neg_hel_Proton_%d.txt", run));

  //Extract Bin Information For the combined runs 
  //Open existing *_final.root file in READ mode (to read in the final combined 2D histos) 

  outROOT = new TFile(data_OutputFileName_combined, "READ");

  //Manually Get each of the four (don't count protons for now) relevant random-subtracted histograms and call the extraction function

  //pos hel
  outROOT->GetObject("hel_plots/H_thxqCM_vs_phxq_K_pos_rand_sub", h2_total);
  extract_2d_hist(h2_total, "phi_xq [deg]", "th_xq_CM [deg]", "OUTPUT/pos_hel_Kaon_combined.txt");

  outROOT->GetObject("hel_plots/H_thxqCM_vs_phxq_Pi_pos_rand_sub", h2_total);
  extract_2d_hist(h2_total, "phi_xq [deg]", "th_xq_CM [deg]", "OUTPUT/pos_hel_Pion_combined.txt");

  //outROOT->GetObject("hel_plots/H_thxqCM_vs_phxq_P_pos_rand_sub", h2_total);
  //extract_2d_hist(h2_total, "phi_xq [deg]", "th_xq_CM [deg]", "pos_hel_Proton_combined.txt");

  //neg hel
  outROOT->GetObject("hel_plots/H_thxqCM_vs_phxq_K_neg_rand_sub", h2_total);
  extract_2d_hist(h2_total, "phi_xq [deg]", "th_xq_CM [deg]", "OUTPUT/neg_hel_Kaon_combined.txt");
  
  outROOT->GetObject("hel_plots/H_thxqCM_vs_phxq_Pi_neg_rand_sub", h2_total);
  extract_2d_hist(h2_total, "phi_xq [deg]", "th_xq_CM [deg]", "OUTPUT/neg_hel_Pion_combined.txt");

  //outROOT->GetObject("hel_plots/H_thxqCM_vs_phxq_P_neg_rand_sub", h2_total);
  //extract_2d_hist(h2_total, "phi_xq [deg]", "th_xq_CM [deg]", "neg_hel_Proton_combined.txt");

  
}
//_______________________________________________________________________________
void helicityAnalyzer::run_helicity_analysis()
{
  
  /*
    Brief:  This method call all the necessary methods to carry out the 
    full helicity analysis. The main controls input parameter is located in: 
    ./main_controls.inp. This file amongst other things, loads the necessary
    input files to set filenames, cuts and histogram binning.
    
    This is supposed to be a derived helicityAnalyzer class which analyzes data in
    and makes used of the existing generic functions/methods of the baseAnalyzer. 
    The analyzer assumes that all the calibrations from the data
    have been done. See each of the methods for details of what the analyzer does. 
    
    Some of the methods called are directly from the baseAnalyzer class. If additional
    functionality is needed, then one can re-define or make additions to the method in 
    the derived class by simply calling the baseAnalyzer directly first, and then, making 
    additions. See, for example, the CreateHist() method in this code.
    
    
  */
  
  SetFileNames();
  SetCuts();
  ReadReport();
  SetHistBins();
  CreateHist();
  
  ReadScalerTree();
  ScalerEventLoop(); 

  ReadTree();
  EventLoop();

  CalcEff();
  ApplyWeight();

  WriteHist();
  WriteReport();
  CombineHistos();

  BinExtraction();
    
}



