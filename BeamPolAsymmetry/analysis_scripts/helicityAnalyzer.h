#ifndef HELICITY_ANALYZER_H
#define HELICITY_ANALYZER_H

#include "baseAnalyzer.h"

//Inheritance: the heli
class helicityAnalyzer : public baseAnalyzer
{
public:
  
  //Constructor / Destructor
  helicityAnalyzer(int irun=-1, string mode="", string earm="", string ana_type="", Bool_t hel_flag=0); 
  ~helicityAnalyzer();

  //MAIN HELICITY ANALYSIS FUNCTIONS
  void run_helicity_analysis();
  
  //Prototypes (obtained from generic baseAnalyzer) -- can be overwritten, or modified
  void SetFileNames();
  void ReadInputFile(string ftype="");
  void SetCuts();
  void ReadReport();
  void SetHistBins();
  void CreateHist();
  void ReadScalerTree();  
  void ScalerEventLoop(); 
  void ReadTree(); 
  void EventLoop(); 
  void CalcEff();
  void ApplyWeight();
  void WriteHist();
  void WriteReport();
  void CombineHistos();

  //Prototypes Specific to this class
  void RandSub();
  void BinExtraction();
  
protected:
  
  //-------------------------------------------------
  //    Kaon LT (2018) HELICITY ANALYSIS SECTION
  //-------------------------------------------------
  
  //SHMS: +hadrons (+Kaons was set in the kinematics file, so the standats missing energy assumes Kaon mass)
  //HMS: -electrons
  //Coincidence trigge: ptrig5 : HMS 3/4 && SHMS 3/4
  //General possible reactions:
  // ---> H(e,e'K) (Lambda or Sigma "missing")
  // ---> H(e,e'Pi) (neutron "missing")
  

  //User-defined variables (to be calculated in event loop)
  Double_t MM2_K;      //missing mass assuming detected particle in SHMS is Kaon
  Double_t MM2_Pi;     // "       " for Pions
  Double_t MM2_P;     //  "       " for Protons

  Double_t MM_K;      //missing mass assuming detected particle in SHMS is Kaon
  Double_t MM_Pi;     // "       " for Pions
  Double_t MM_P;     //  "       " for Protons

  
  //Helicity Scaler-related variables
  //(On Fall 2018, useless helicity scalers)
  //After Kaon LT Fall 2018 ended, helicity scalers hardware was updated
  //but the helicity scaler software is a work in progress . . .

  //Helicity scalers need to be added for beamAsymmetry experiments after Fall 2018.

  
  //Define Helicity Cuts
  Int_t hel_sign;     //helicity sign (either +1 or -1 or 0)
  Bool_t c_pos_hel;   //cut to select ONLY positive helicity events
  Bool_t c_neg_hel;   //cut to select ONLY negative helicity events

  //-----------------------------------------------------------------------------
  //------ HELICITY DATA ANALYSIS PARAMTER CUTS ( see set_basic_KaonLT.inp ) ----
  //-----------------------------------------------------------------------------

  //---accpetance cuts variables---
  Bool_t accp_cuts;
  
  Float_t edelta_min;  
  Float_t edelta_max;  

  Float_t hdelta_min; 
  Float_t hdelta_max;

  Float_t ztarDiff_min;
  Float_t ztarDiff_max;
  
  //----PID cuts variables-----
  //-electron PID-
  Bool_t elec_pid;
  Float_t elec_hcer_npe_thrs;      //threhsold on HMS cherenkov to select electrons
  Float_t elec_hcal_thrs;          //threhsold on HMS calorimeter total energy normalized by track momentum to select electrons
  
  //-Kaon PID-
  Bool_t kaon_pid;

  Bool_t K_paero_npe_flag;
  Bool_t c_K_paero_npe;
  Float_t K_paero_npe_thrs;          //threhsold on aerogel number of photo-electrons
  
  Bool_t K_phgcer_npe_flag;
  Bool_t c_K_phgcer_npe;
  Float_t K_phgcer_npe_thrs;         //threshold on heavy gas cherenkv number of photo-electrons

  
  Bool_t K_beta_flag;
  Bool_t kaon_beta_cut;
  Float_t K_beta_min, K_beta_max;    //particle beta cut

  Bool_t eK_ctime_flag;
  Bool_t eK_ctime_cut;
  Bool_t eK_ctime_cut_rand;
  Float_t eK_ctime_thrs;             //threshold minimum cut on absolute value of electron-Kaon coin time: e.g.,  |ctime| < eK_ctime_min

  Bool_t kaon_MM_cut;

  //-Pion PID-
  Bool_t pion_pid;

  Bool_t Pi_phgcer_npe_flag;
  Bool_t c_Pi_phgcer_npe;
  Float_t Pi_phgcer_npe_thrs;        //threshold on heavy gas cherenkv number of photo-electrons

  
  Bool_t Pi_beta_flag;
  Bool_t pion_beta_cut;
  Float_t Pi_beta_min, Pi_beta_max;  //particle beta cut

  Bool_t ePi_ctime_flag;
  Bool_t ePi_ctime_cut;
  Bool_t ePi_ctime_cut_rand;
  Float_t ePi_ctime_thrs;            //threshold minimum cut on absolute value of electron-Pion coin time: e.g.,  |ctime| < eK_ctime_min

  Bool_t pion_MM_cut;


  //-Proton PID-

  Bool_t proton_pid;

  Bool_t proton_beta_flag;
  Bool_t proton_beta_cut;
  Float_t P_beta_min, P_beta_max;    //particle beta cut

  Bool_t eP_ctime_flag;
  Bool_t eP_ctime_cut;
  Bool_t eP_ctime_cut_rand;

  Bool_t proton_MM_cut;

  Bool_t P_aero_npe_flag;
  Float_t P_paero_npe_thrs;          //threhsold on aerogel number of photo-electrons

  Bool_t P_phgcer_npe_flag;
  Float_t P_phgcer_npe_thrs;         //threshold on heavy gas cherenkv number of photo-electrons

  Bool_t eP_ctime_flag;
  Float_t eP_ctime_thrs;             //threshold minimum cut on absolute value of electron-Proton coin time: e.g.,  |ctime| < eK_ctime_min

  //----Kinematics Cuts----

  //Missing Mass
  Float_t MM_K_min;     //Kaon Missing Mass (missing Lambda 1115)
  Float_t MM_K_max;

  Float_t MM_Pi_min;   //Pion Missing Mass (missing neutron 0.939)
  Float_t MM_Pi_max;
  
  Float_t MM_P_min;    //Proton Missing Mass (missing " . . . " nothing ?)
  Float_t MM_P_max;


  //----Coincidence Time Parameters
  Float_t K_ctime_offset;            //coin time offset parameter to center coin time spectrum at 0 ns.
  Float_t Pi_ctime_offset;           //coin time offset parameter to center coin time spectrum at 0 ns.
  Float_t P_ctime_offset;            //coin time offset parameter to center coin time spectrum at 0 ns.

  Float_t eK_mult;                   //coin time integer multiple of eK_ctime_thrs cut used to select randoms (must be: >=2: i,e, 2, 3, 4, . . .) 
  Float_t ePi_mult;                  //coin time integer multiple of ePi_ctime_thrs cut used to select randoms (must be: >=2: i,e, 2, 3, 4, . . .) 
  Float_t eP_mult;                   //coin time integer multiple of eP_ctime_thrs cut used to select randoms (must be: >=2: i,e, 2, 3, 4, . . .) 

  
  //Scale factor variables (for random coincidence scaling / subtraction)
  Float_t K_scale_factor;     //scale factor used to scale down random eK coincidences before random coincidence subtraction
  Float_t Pi_scale_factor;
  Float_t P_scale_factor;

  //Tracking Efficiency Counter for +/- helicities (passed bcm cuts)
  //HMS + helicity
  Double_t h_did_pos = 0;
  Double_t h_should_pos = 0;

  Double_t hTrkEff_pos;
  Double_t hTrkEff_err_pos;

  //HMS - helicity
  Double_t h_did_neg = 0;
  Double_t h_should_neg = 0;
  
  Double_t hTrkEff_neg;
  Double_t hTrkEff_err_neg;
  
  //SHMS + helicity
  Double_t p_did_pos = 0;
  Double_t p_should_pos = 0;

  Double_t pTrkEff_pos;
  Double_t pTrkEff_err_pos;

  //SHMS - helicity
  Double_t p_did_neg = 0;
  Double_t p_should_neg = 0;

  Double_t pTrkEff_neg;
  Double_t pTrkEff_err_neg;
  
  //------VARIABLES USED TO WRITE HISTOGRAMS TO ROOT FILE-------

  //Create Categorical TLists to store histograms based on caterogy
  TList *hel_HList;    //store helicity analysis histos (all histos declared in this class)

  //---------------------------------------------
  
  //-----------------------------------------------------------------
  //---------Declare Kaon LT (HELICITY ANALYSIS) Histograms----------
  //-----------------------------------------------------------------

  //NOTE: Nomenclature clarification 
  // *_real -> real coincidences selection (within main peak -- may still have bakground underneath),
  // *_rand -> random coincidence selection (outside the main peak)  
  // *_rand_sub -> "true" coincidences after having subtracted the estimated randoms beneath the main peak

  // 1D Histograms

  // Coincidence Time Histograms for real / random coincidences selection
  TH1F *H_eK_ctime_real;
  TH1F *H_ePi_ctime_real;
  TH1F *H_ep_ctime_real;

  TH1F *H_eK_ctime_rand;
  TH1F *H_ePi_ctime_rand;
  TH1F *H_ep_ctime_rand;

  // Missing Mass Histograms (only PID cuts, not coin. cuts)
  TH1F *H_MM_K_PID;
  TH1F *H_MM_Pi_PID;
  TH1F *H_MM_P_PID;
  
  // Missing Mass Histograms from real coincidences
  TH1F *H_MM_K_real;
  TH1F *H_MM_Pi_real;
  TH1F *H_MM_P_real;

  // Missing Mass Histograms from random coincidences selection
  TH1F *H_MM_K_rand;
  TH1F *H_MM_Pi_rand;
  TH1F *H_MM_P_rand;

  // Missing Mass Histograms after random coincidence subtraction from real
  TH1F *H_MM_K_rand_sub;
  TH1F *H_MM_Pi_rand_sub;
  TH1F *H_MM_P_rand_sub;

  // Particle Beta (no PID cuts) -should only be one histo
  //(as there is no specific particle being selected)
  TH1F *H_Beta_noPID;
  
  // Particle Beta (pid selection)
  TH1F *H_Beta_K_PID;
  TH1F *H_Beta_Pi_PID;
  TH1F *H_Beta_P_PID;

  // Particle Beta (real coincidence selection)
  TH1F *H_Beta_K_real;
  TH1F *H_Beta_Pi_real;
  TH1F *H_Beta_P_real;

  // Particle Beta (random coincidence selection)
  TH1F *H_Beta_K_rand;
  TH1F *H_Beta_Pi_rand;
  TH1F *H_Beta_P_rand;

  // Particle Beta (random coincidence subtracted)
  TH1F *H_Beta_K_rand_sub;
  TH1F *H_Beta_Pi_rand_sub;
  TH1F *H_Beta_P_rand_sub;
  
  //2D Histograms

  // Beta vs. Coincidence Time (no PID selection cuts)
  TH2F *H_Beta_vs_ctime_K_noPID;
  TH2F *H_Beta_vs_ctime_Pi_noPID;
  TH2F *H_Beta_vs_ctime_P_noPID;

  // Beta vs. Coincidence Time (PID selection cuts)
  TH2F *H_Beta_vs_ctime_K_PID;
  TH2F *H_Beta_vs_ctime_Pi_PID;
  TH2F *H_Beta_vs_ctime_P_PID;
  
  // Kaon LT Kinematic Histograms (fill with Kaons for now . . .)
  TH2F *H_Q2_vs_W;    //supposed to be a diamond shape
  TH2F *H_t_vs_phxq;  //mandelstam -t vs phi_xq

  
  TH2F *H_Q2_vs_W_real;    
  TH2F *H_t_vs_phxq_real;  

  TH2F *H_Q2_vs_W_rand;   
  TH2F *H_t_vs_phxq_rand;  

  TH2F *H_Q2_vs_W_rand_sub;    
  TH2F *H_t_vs_phxq_rand_sub;  


  //NOTE: The advantage of using the 2D histos below is that one can extract asymmetry information from ph_xq for each th_xq_cm bin
  //easily via a loop over each of the 2D bins

  //==== KAONS ====
  //2D: In-plane angles of X about q, th_xq_cm (in CM frame) "VS"  Out-of-plane angles of X about q, ph_xq (Assume detected X are Kaons)
  TH2F *H_thxqCM_vs_phxq_K_pos;    //+ helicity
  TH2F *H_thxqCM_vs_phxq_K_neg;   //- helicity

  //Same as above, but for randoms (events out side the main coin. peak)
  TH2F *H_thxqCM_vs_phxq_K_pos_rand;  //+ helicity
  TH2F *H_thxqCM_vs_phxq_K_neg_rand;  //- helicity

  //Same as above, but for randoms subtracted
  //Will use these histos to extract bin information for +/- helicity into a data file and calculate asymmetries from data file.
  TH2F *H_thxqCM_vs_phxq_K_pos_rand_sub;  //+ helicity
  TH2F *H_thxqCM_vs_phxq_K_neg_rand_sub;  //- helicity
  
  //==== PIONS ====
  //2D:  In-plane angles of X about q, th_xq_cm (in CM frame) "VS"  Out-of-plane angles of X about q, ph_xq (Assume detected X are Pions)
  TH2F *H_thxqCM_vs_phxq_Pi_pos;  //+ helicity
  TH2F *H_thxqCM_vs_phxq_Pi_neg;  //- helicity

  //Same as above, but for randoms (events outside the main coin. peak)
  TH2F *H_thxqCM_vs_phxq_Pi_pos_rand;  //+ helicity
  TH2F *H_thxqCM_vs_phxq_Pi_neg_rand;  //- helicity

  //Same as above, but for randoms subtracted
  //Will use these histos to extract bin information for +/- helicity into a data file and calculate asymmetries from data file.
  TH2F *H_thxqCM_vs_phxq_Pi_pos_rand_sub;  //+ helicity
  TH2F *H_thxqCM_vs_phxq_Pi_neg_rand_sub;  //- helicity
  
  //==== PROTONS ====
  //2D:  In-plane angles of X about q, th_xq_cm (in CM frame) "VS"  Out-of-plane angles of X about q, ph_xq (Assume detected X are Protons)
  TH2F *H_thxqCM_vs_phxq_P_pos;  //+ helicity
  TH2F *H_thxqCM_vs_phxq_P_neg;  //- helicity

  //Same as above, but for randoms (events outside the main coin. peak)
  TH2F *H_thxqCM_vs_phxq_P_pos_rand;  //+ helicity
  TH2F *H_thxqCM_vs_phxq_P_neg_rand;  //- helicity

  //Same as above, but for randoms subtracted
  //Will use these histos to extract bin information for +/- helicity into a data file and calculate asymmetries from data file.
  TH2F *H_thxqCM_vs_phxq_P_pos_rand_sub;  //+ helicity
  TH2F *H_thxqCM_vs_phxq_P_neg_rand_sub;  //- helicity


  
};

#endif 
