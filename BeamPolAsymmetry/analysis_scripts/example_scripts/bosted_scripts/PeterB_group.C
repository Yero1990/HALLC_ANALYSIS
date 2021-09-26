// Filename: PeterB.C
// Calculates cointime for 3 masses. Makes skim files with needed variables for good el. in HMS
// Needs runnumber, and if target is "h" or "d"
#include <stdio.h>  

void PeterB(Int_t runNumber, Int_t targ=1){

  // Skim files
  TString fileNameSkim = 
"/w/hallc-scifs17exp/c-sidis18/bosted/Skimfiles//Skim" ;
  TString fileNameSkimeff = 
"/w/hallc-scifs17exp/c-sidis18/bosted/Skimfiles//Skimeff" ;
  TString fileNameSkimp = 
"/w/hallc-scifs17exp/c-sidis18/bosted/Skimfiles//Skimp" ;
  TString fileNameChist = 
"/w/hallc-scifs17exp/c-sidis18/bosted/Skimfiles//Chist" ;
  fileNameSkim += runNumber; 
  fileNameSkimp += runNumber; 
  fileNameSkimeff += runNumber; 
  fileNameChist += runNumber; 
  fileNameSkim += ".txt"; 
  fileNameSkimp += ".txt"; 
  fileNameSkimeff += ".txt"; 
  fileNameChist += ".top";
  //fileNameSkim += "_20K.txt"; 
  //fileNameSkimp += "_20K.txt"; 
  FILE *f2 = fopen(fileNameSkim,"w");
  FILE *f5 = fopen(fileNameSkimeff,"w");
  FILE *f3 = fopen(fileNameSkimp,"w");
  //  FILE *f4 = fopen("test.top","w");
  FILE *f6 = fopen(fileNameChist,"w");


  //read the input file from data
  //TString fileNameD = "/w/hallc-scifs17exp/c-sidis18/bosted/ROOTfiles/" ;
  TString fileNameD = "/volatile/hallc/spring17/bosted/ROOTfiles/" ;
 fileNameD += "coin_replay_production_"; //read the root file from data
  fileNameD += runNumber; //read the root file from data
  fileNameD += "_-1.root"; //read the root file from data
  //fileNameD += "_20000.root"; //read the root file from data

  TFile *f1 = new TFile(fileNameD);
  if(f1->GetSize()!=-1){ 

  TTree *tt = (TTree*)f1->Get("T");
  //get the relevant branch
  int nentriesD = tt->GetEntries();
  cout<<"Entries:\t"<<nentriesD<<endl;
  TString fileO;
  fileO += "HISTOGRAMS/COIN/ROOT/run_"; //read the root file from data
  fileO += runNumber; //read the root file from data
  fileO += "_hists_coin.root"; //read the root file from data

  //TFile *fout = new TFile(fileO,"RECREATE");

  gROOT->SetBatch(kTRUE);

  
  //make histograms:

  TH1D *h1_epi_PcointimeROC2    = new TH1D("SHMS ROC2 Corrected epi Coin Time","SHMS ROC2 Corrected epi Coin Time; cointime [ns]",       480, -24, 24); 
  TH1D *h1_epi_PcointimeROC2_R  = new TH1D("SHMS ROC2 Corrected epi Coin Time_R","SHMS ROC2 Corrected epi Coin Time_R; cointime [ns]",   480, -24, 24); 
  TH1D *h1_epi_PcointimeROC2_L  = new TH1D("SHMS ROC2 Corrected epi Coin Time_L","SHMS ROC2 Corrected epi Coin Time_L; cointime [ns]",   480, -24, 24); 
  TH1D *h1_epi_PcointimeROC2_C  = new TH1D("SHMS ROC2 Corrected epi Coin Time_C","SHMS ROC2 Corrected epi Coin Time_C; cointime [ns]",   480, -24, 24); 

  TH1D *h1_ep_PcointimeROC2    = new TH1D("SHMS ROC2 Corrected ep Coin Time","SHMS ROC2 Corrected ep Coin Time; cointime [ns]",       480, -24, 24); 
  TH1D *h1_ep_PcointimeROC2_R  = new TH1D("SHMS ROC2 Corrected ep Coin Time_R","SHMS ROC2 Corrected ep Coin Time_R; cointime [ns]",   480, -24, 24); 
  TH1D *h1_ep_PcointimeROC2_L  = new TH1D("SHMS ROC2 Corrected ep Coin Time_L","SHMS ROC2 Corrected ep Coin Time_L; cointime [ns]",   480, -24, 24); 
  TH1D *h1_ep_PcointimeROC2_C  = new TH1D("SHMS ROC2 Corrected ep Coin Time_C","SHMS ROC2 Corrected ep Coin Time_C; cointime [ns]",   480, -24, 24); 

  TH1D *h1_ek_PcointimeROC2    = new TH1D("SHMS ROC2 Corrected ek Coin Time","SHMS ROC2 Corrected ek Coin Time; cointime [ns]",       480, -24, 24); 
  TH1D *h1_ek_PcointimeROC2_R  = new TH1D("SHMS ROC2 Corrected ek Coin Time_R","SHMS ROC2 Corrected ek Coin Time_R; cointime [ns]",   480, -24, 24); 
  TH1D *h1_ek_PcointimeROC2_L  = new TH1D("SHMS ROC2 Corrected ek Coin Time_L","SHMS ROC2 Corrected ek Coin Time_L; cointime [ns]",   480, -24, 24); 
  TH1D *h1_ek_PcointimeROC2_C  = new TH1D("SHMS ROC2 Corrected ek Coin Time_C","SHMS ROC2 Corrected ek Coin Time_C; cointime [ns]",   480, -24, 24); 

  Double_t HgtrX, HgtrTh, HgtrPh, hdelta, PgtrX, PgtrTh, PgtrPh, pdelta;
  Double_t hdcx,hdcxp,hdcy,hdcyp,pdcx,pdcxp,pdcy,pdcyp;
  Double_t htime1, htime2, hhits1, hhits2,Hgtrxtar ;
  Double_t Hgtryptar, Pgtrxtar,Pgtryptar ;
  Double_t ptime1, ptime2, phits1, phits2,poffsettime,hoffsettime ;
  Double_t HgtrBetaCalc, PgtrBetaCalc, evtType, PgtrP, HgtrP, PhodStatus, PhodStartTime, PhodfpHitsTime;
  Double_t cointime, HhodStatus, HhodStartTime, HhodfpHitsTime, paeronpe,paeropt[99],paeront[99],paeroptd[99],paerontd[99];
  Double_t aero_avtime,ntime ; 
  Double_t pkinW, pEm, pPm, modPm, pbeta, hbeta, hcalepr, hcaletot, hcernpe, pcaletot, pcalepr, pcernpe, pcernpeng;
  Double_t TcoinpTRIG1_ROC1_tdcTimeRaw, TcoinpTRIG4_ROC1_tdcTimeRaw, TcoinpTRIG1_ROC2_tdcTimeRaw;
  Double_t TcoinhTRIG1_ROC1_tdcTimeRaw, TcoinhTRIG1_ROC2_tdcTimeRaw, TcoinhTRIG4_ROC1_tdcTimeRaw;
  Double_t t11,t31,t41,t51,t61,t12,t32,t42,t52,t62,diff1,diff3,diff4,diff5,diff6,tdiff1,tdiff2,td3,td4 ;
  Double_t TcoinhTRIG4_ROC2_tdcTimeRaw, TcoinpTRIG4_ROC2_tdcTimeRaw;
  Double_t TcoinpTRIG3_ROC1_tdcTimeRaw, TcoinpTRIG3_ROC2_tdcTimeRaw;
  Double_t TcoinpTRIG5_ROC1_tdcTimeRaw, TcoinpTRIG5_ROC2_tdcTimeRaw;
  Double_t TcoinpTRIG6_ROC1_tdcTimeRaw, TcoinpTRIG6_ROC2_tdcTimeRaw,sum1,sum2,counters[10];
  Double_t  helmpsa,helnega,helposa,hgtry,pgtry,pcaltot,pgdsc,gdtrk,trkeff,trkeffer,prf, hrf, prfraw, hrfraw, hztar, pztar,trkeffg,trkeffger,trkefft,trkeffter,htof1,htof2 ;
  Double_t  hgdsc,htrkeff,htrkeffer,phbig,ppbig,t1big,t2big,sum3,sum4,sum5,sum6,sum7,sum8 ;
  //  Int_t hel1,helpred,helrep,helmps,helnqrt ;                     
  Double_t hel1,helpred,helrep,helmps,helnqrt,hcaltot ;  
  Int_t  helmpsi,helnegi,helposi,icc,chist[100],chistp[100],chistpt[100],chistpg[100],chistt[100],chistg[100] ;
  Int_t evtypecnr[10],ptdcmulth[15] ;
  Int_t  chistb[20][10],ibeta,ctpih[100],ctbig[100],chistww[1200][2] ;
  Int_t  cch1[100], cch3[100], cch2[100],t41h[1000],goodsch,goodscp,stgoodsch,stgoodscp,iccc ;
  Double_t  cthx[100][10],tsfdiff,tsfdiffh,tsfh[100],tsfhh[100] ;
  Int_t ctrawh[200][2], ct2h[200][2] ;
  Double_t  sthh[102],sthp[102],fphh[102],fphp[102],sttoth,sttotp,fpgoodh,fpgoodp,stgoodh,stgoodp ;
  Double_t  hdchit1,hdchit2,pdchit1,pdchit2,pntrk,hntrk,ngcorr ;
  Double_t  hncombo1[10000],pncombo1[10000] ;
  Double_t  hncombo2[10000],pncombo2[10000],ptdcmult ;
  Int_t  hnumcombo1, pnumcombo1,hnumcombo2, pnumcombo2 ;
  Double_t ctpi1,ctpi2,ctk1,ctk2,ctp1,ctp2,hst,pst,hpst,ctpi3,ctpi4;
  Int_t cntsepi=0, it3;
  Int_t cntsep=0;
  Int_t cntsek=0;
  Int_t cntpos=0;
  Int_t nevnt=0;
  Int_t nevntall=0;
  Int_t nevnt32=0;
  Int_t nevnt32h=0;
  Int_t nevntjust32=0;
  Int_t nevntjust32h=0;
  Int_t nevntjust32st=0;
  Int_t nevntjust32hst=0;

  for(icc=0 ; icc<1200 ;  icc++) {
    chistww[icc][0]=0. ;
    chistww[icc][1]=0. ;
  }
  for(icc=0 ; icc<100 ;  icc++) {
    tsfh[icc]=0. ;
    tsfhh[icc]=0. ;
    chist[icc]=0 ;
    chistp[icc]=0 ;
    chistpt[icc]=0 ;
    chistpg[icc]=0 ;
    chistt[icc]=0 ;
    chistg[icc]=0 ;
    ctpih[icc]=0 ;
    ctbig[icc]=0 ;
    cch1[icc]=0 ;
    cch2[icc]=0 ;
    cch3[icc]=0 ;
  }
  for(icc=0 ; icc<1000 ;  icc++) {
    t41h[icc]=0 ;
  }
  for(icc=0 ; icc<10 ;  icc++) {
    counters[icc]=0 ;
    evtypecnr[icc]=0 ;
  }
  for(icc=0 ; icc<15 ;  icc++) {
    ptdcmulth[icc]=0 ;
  }
  for(icc=0 ; icc<100 ;  icc++) {
   for(iccc=0 ; iccc<10 ;  iccc++) {
     cthx[icc][iccc]=0 ;
   }
  }
  for(icc=0 ; icc<20 ;  icc++) {
    for(ibeta=0 ; ibeta<10 ;  ibeta++) {
      chistb[icc][ibeta]=0 ;
    }
  }
  for(icc=0 ; icc<102 ;  icc++) {
    sthh[icc]=0 ;
    sthp[icc]=0 ;
    fphh[icc]=0 ;
    fphp[icc]=0 ;
  }
  sttoth = 0 ;
  sttotp = 0 ;
  fpgoodh = 0 ;
  fpgoodp = 0 ;
  stgoodh = 0 ;
  stgoodp = 0 ;
  tt->SetBranchAddress("CTime.ePiCoinTime_ROC1",  &ctpi1); 
  tt->SetBranchAddress("CTime.ePiCoinTime_ROC2",  &ctpi2); 
  tt->SetBranchAddress("CTime.eKCoinTime_ROC1",  &ctk1); 
  //tt->SetBranchAddress("CTime.eKCoinTime_ROC2",  &ctk2); 
  tt->SetBranchAddress("CTime.epCoinTime_ROC1",  &ctp1); 
  //tt->SetBranchAddress("CTime.epCoinTime_ROC2",  &ctp2); 
 
  tt->SetBranchAddress("P.dc.x_fp",  &pdcx); 
  tt->SetBranchAddress("P.dc.y_fp",  &pdcy); 
  tt->SetBranchAddress("P.dc.xp_fp", &pdcxp); 
  tt->SetBranchAddress("P.dc.yp_fp", &pdcyp); 
  tt->SetBranchAddress("H.dc.x_fp",  &hdcx); 
  tt->SetBranchAddress("H.dc.y_fp",  &hdcy); 
  tt->SetBranchAddress("H.dc.xp_fp", &hdcxp); 
  tt->SetBranchAddress("H.dc.yp_fp", &hdcyp); 

  tt->SetBranchAddress("P.gtr.y", &pgtry); 
  tt->SetBranchAddress("H.gtr.y", &hgtry); 
  tt->SetBranchAddress("P.gtr.x", &PgtrX); 
  tt->SetBranchAddress("H.gtr.x", &HgtrX); 
  tt->SetBranchAddress("P.gtr.p", &PgtrP); 
  tt->SetBranchAddress("H.gtr.p", &HgtrP); 
  tt->SetBranchAddress("P.gtr.beta", &pbeta);   
  tt->SetBranchAddress("H.gtr.beta", &hbeta); 
  tt->SetBranchAddress("H.gtr.dp", &hdelta);   
  tt->SetBranchAddress("P.gtr.dp", &pdelta);  
  tt->SetBranchAddress("H.gtr.th", &HgtrTh);   
  tt->SetBranchAddress("P.gtr.th", &PgtrTh);  
  tt->SetBranchAddress("H.gtr.ph", &HgtrPh);   
  tt->SetBranchAddress("P.gtr.ph", &PgtrPh);  
  tt->SetBranchAddress("H.react.z", &hztar);  
  tt->SetBranchAddress("P.react.z", &pztar);  
  //tt->SetBranchAddress("H.gtr.xtar", &Hgtrxtar);   
  //tt->SetBranchAddress("P.gtr.xtar", &Pgtrxtar);   
  tt->SetBranchAddress("P.cal.eprtracknorm", &pcalepr);
  tt->SetBranchAddress("P.cal.etottracknorm", &pcaletot); 
  tt->SetBranchAddress("P.cal.etotnorm", &pcaltot); 
  tt->SetBranchAddress("P.ngcer.npeSum", &pcernpeng);          
  tt->SetBranchAddress("P.hgcer.npeSum", &pcernpe);          
  tt->SetBranchAddress("P.aero.npeSum", &paeronpe); 
  tt->SetBranchAddress("P.aero.goodPosAdcPulseTime", &paeropt); 
  tt->SetBranchAddress("P.aero.goodPosAdcTdcDiffTime", &paeroptd); 
  tt->SetBranchAddress("P.aero.goodNegAdcPulseTime", &paeront); 
  tt->SetBranchAddress("P.aero.goodNegAdcTdcDiffTime", &paerontd); 
  tt->SetBranchAddress("H.cal.eprtracknorm", &hcalepr);
  tt->SetBranchAddress("H.cal.etottracknorm", &hcaletot);  
  tt->SetBranchAddress("H.cal.etotnorm", &hcaltot); 
  tt->SetBranchAddress("H.cer.npeSum", &hcernpe); 
  tt->SetBranchAddress(
   "T.coin.pT2_tdcMultiplicity", &ptdcmult);
  tt->SetBranchAddress(
   "T.coin.pHEL_MPS_adcMultiplicity", &helmpsa);
  tt->SetBranchAddress(
   "T.coin.pHEL_NEG_adcMultiplicity", &helnega);
  tt->SetBranchAddress(
   "T.coin.pHEL_POS_adcMultiplicity", &helposa);
  tt->SetBranchAddress("H.hod.goodscinhit", &hgdsc);
  tt->SetBranchAddress("P.hod.goodscinhit", &pgdsc);
  tt->SetBranchAddress("T.helicity.hel",&hel1) ;
  tt->SetBranchAddress("T.helicity.helpred",&helpred) ;
  tt->SetBranchAddress("T.helicity.helrep",&helrep) ; 
  tt->SetBranchAddress("T.helicity.mps",&helmps) ;
  tt->SetBranchAddress("T.helicity.nqrt",&helnqrt) ;
                    
  tt->SetBranchAddress("P.hod.starttime", &PhodStartTime);                                               
  tt->SetBranchAddress("P.hod.offsettime", &poffsettime);                                               
  tt->SetBranchAddress("P.hod.fpHitsTime", &PhodfpHitsTime);                                             
  tt->SetBranchAddress("H.hod.starttime", &HhodStartTime);                                               
  tt->SetBranchAddress("H.hod.offsettime", &hoffsettime);                                               
  tt->SetBranchAddress("H.hod.fpHitsTime", &HhodfpHitsTime); 
  tt->SetBranchAddress("T.coin.pTRIG1_ROC1_tdcTimeRaw", &TcoinpTRIG1_ROC1_tdcTimeRaw);                   
  //tt->SetBranchAddress("T.coin.pTRIG3_ROC1_tdcTimeRaw", &TcoinpTRIG3_ROC1_tdcTimeRaw);                   
  tt->SetBranchAddress("T.coin.pTRIG4_ROC1_tdcTimeRaw", &TcoinpTRIG4_ROC1_tdcTimeRaw);
  //tt->SetBranchAddress("T.coin.pTRIG5_ROC1_tdcTimeRaw", &TcoinpTRIG5_ROC1_tdcTimeRaw);                   
  //tt->SetBranchAddress("T.coin.pTRIG6_ROC1_tdcTimeRaw", &TcoinpTRIG6_ROC1_tdcTimeRaw);                   
  tt->SetBranchAddress("T.coin.pTRIG1_ROC2_tdcTimeRaw", &TcoinpTRIG1_ROC2_tdcTimeRaw);                   
  //tt->SetBranchAddress("T.coin.pTRIG3_ROC2_tdcTimeRaw", &TcoinpTRIG3_ROC2_tdcTimeRaw);                   
  tt->SetBranchAddress("T.coin.pTRIG4_ROC2_tdcTimeRaw", &TcoinpTRIG4_ROC2_tdcTimeRaw);                   
  //tt->SetBranchAddress("T.coin.pTRIG5_ROC2_tdcTimeRaw", &TcoinpTRIG5_ROC2_tdcTimeRaw);                   
  //tt->SetBranchAddress("T.coin.pTRIG6_ROC2_tdcTimeRaw", &TcoinpTRIG6_ROC2_tdcTimeRaw);                   
  tt->SetBranchAddress("T.coin.pRF_tdcTime",&prf) ;   
  tt->SetBranchAddress("T.coin.hRF_tdcTime",&hrf) ;   
  tt->SetBranchAddress("T.coin.pRF_tdcTimeRaw",&prfraw) ;   
  tt->SetBranchAddress("T.coin.hRF_tdcTimeRaw",&hrfraw) ;   
  //tt->SetBranchAddress("T.coin.hSTOF_ROC2",&htof1) ;   
  //tt->SetBranchAddress("T.coin.hSTOF_ROC1",&htof2) ;   
  //tt->SetBranchAddress("P.dc.Ch1.nhit", &pdchit1); 
  //tt->SetBranchAddress("P.dc.Ch2.nhit", &pdchit2); 
  tt->SetBranchAddress("P.dc.ntrack", &pntrk); 
  //tt->SetBranchAddress("H.dc.Ch1.nhit", &hdchit1); 
  //tt->SetBranchAddress("H.dc.Ch2.nhit", &hdchit2); 
  tt->SetBranchAddress("H.dc.ntrack", &hntrk); 
  //tt->SetBranchAddress("P.dc.Ch1.nhit", &pdchit1); 
  //tt->SetBranchAddress("Ndata.H.dc.Ch1.ncombos", &hnumcombo1); 
  //tt->SetBranchAddress("Ndata.H.dc.Ch2.ncombos", &hnumcombo2); 
  //tt->SetBranchAddress("Ndata.P.dc.Ch1.ncombos", &pnumcombo1); 
  //tt->SetBranchAddress("Ndata.P.dc.Ch2.ncombos", &pnumcombo2); 
  //tt->SetBranchAddress("H.dc.Ch1.ncombos", &hncombo1); 
  //tt->SetBranchAddress("H.dc.Ch2.ncombos", &hncombo2); 
  //tt->SetBranchAddress("P.dc.Ch1.ncombos", &pncombo1); 
  //tt->SetBranchAddress("P.dc.Ch2.ncombos", &pncombo2); 
  tt->SetBranchAddress("H.hod.fhodo_time_peak1", &htime1);
  tt->SetBranchAddress("H.hod.fhodo_time_peak2", &htime2);
  tt->SetBranchAddress("H.hod.fhodo_time_hits1", &hhits1);
  tt->SetBranchAddress("H.hod.fhodo_time_hits2", &hhits2);
  tt->SetBranchAddress("P.hod.fhodo_time_peak1", &ptime1);
  tt->SetBranchAddress("P.hod.fhodo_time_peak2", &ptime2);
  tt->SetBranchAddress("P.hod.fhodo_time_hits1", &phits1);
  tt->SetBranchAddress("P.hod.fhodo_time_hits2", &phits2);
 
  TCut hpdelta;
  TCut epiCut;                                                                   



  hpdelta = "P.gtr.dp > -10 && P.gtr.dp < 20 && H.gtr.dp > -10 && H.gtr.dp < 10";
  epiCut = "P.aero.npeSum > 1.0 && P.hgcer.npeSum > 1.0 && P.cal.eprtracknorm < 0.2 && H.cer.npeSum > 1.0 && H.cal.etottracknorm > 0.6 && H.cal.etottracknorm < 2.0 && H.cal.eprtracknorm  > 0.2"; 

  

  //  TCanvas *canvas1 = new TCanvas("canvas1","canvas1");                           
  // tt->Draw("P.hod.starttime >> SHMShodoStartTime", epiCut && hpdelta );  
  // TH1D *h1PhodoStartTime = (TH1D*)gDirectory->Get("SHMShodoStartTime");
  // h1PhodoStartTime->GetXaxis()->SetTitle("SHMS hodo start time [ns]");           
  // Double_t PhodoStartTimeMean = h1PhodoStartTime->GetMean();  
  Double_t PhodoStartTimeMean = 0. ;  


                   
  // TCanvas *canvas2 = new TCanvas("canvas2","canvas2");                           
  // tt->Draw("H.hod.starttime >> HMShodoStartTime", epiCut && hpdelta );  
  // TH1D *h1HhodoStartTime = (TH1D*)gDirectory->Get("HMShodoStartTime");           
  // h1HhodoStartTime->GetXaxis()->SetTitle("HMS hodo start time [ns]");            
  // Double_t HhodoStartTimeMean = h1HhodoStartTime->GetMean();                     
  Double_t HhodoStartTimeMean = 0.;


 
  Double_t HMSpartMass = 0.000510998; // electron mass in GeV/c^2
  Double_t SHMSpartMass;                  
  //  Double_t pOffset = 1.5; //9.5 + 10;  // in ns                     
  Double_t pOffset = 15.5   ; //9.5 + 10;  // in ns                                  
  Double_t hOffset = 335;                                                        
  Double_t speedOfLight = 29.9792; // in cm/ns                                   
  Double_t SHMScentralPathLen = 18.1*100;  // SHMS Target to focal plane path length converted to cm  
  Double_t SHMSpathLength = 18.1*100;      // For now assume that it's same as SHMScentralPathLen  
  Double_t HMScentralPathLen = 22*100;     // HMS Target to focal plane path length converted to cm
  Double_t DeltaHMSpathLength;             // For now assume that it's same as HMScentralPathLen 
  Double_t SHMScoinCorr = 0.0;                                                   
  Double_t HMScoinCorr = 0.0;                                                    
  Double_t SHMSrawCoinTimeROC1;                                                  
  Double_t SHMSrawCoinTimeROC2;                                                  
  Double_t HMSrawCoinTimeROC1;                                                   
  Double_t HMSrawCoinTimeROC2;                                                   
  Double_t ctimepi, ctimepinew, ctimeK, ctimep,ctimeraw,tmp ;                          
                                                                                 
  Double_t SHMScorrCoinTimeROC1;                                                 
  Double_t SHMScorrCoinTimeROC2;                                                 
  Double_t HMScorrCoinTimeROC1;                                                  
  Double_t HMScorrCoinTimeROC2;   
  Double_t  ncoin=0,ncoinf=0,ngdsc1=0.,ngdsc2=0.;
  Double_t  hngdsc1=0.,hngdsc2=0.;
  Double_t  ncoinwtrk=0,ncoinwtrk1=0,ncoinwtrk2=0;
  Double_t  hncoin=0,hncoinf=0,mjcorr;
  Double_t  hncoinwtrk=0,hncoinwtrk1=0,hncoinwtrk2=0;
  Double_t  ncn=0, ncnwtrk=0, hncn=0, hncnwtrk=0 ;   
  Bool_t epievent_cut, epevent_cut, ekevent_cut, positron_cut, event_cut, hpdelta_cut;
  Int_t npassed=0;

  //

  HhodoStartTimeMean = 30. ;
  PhodoStartTimeMean = 56. ;
  //fprintf(f2,"%8.2f %8.2f\n",HhodoStartTimeMean, PhodoStartTimeMean) ;
  // printf("%8.2f %8.2f\n",HhodoStartTimeMean, PhodoStartTimeMean) ;


    for(int iii=0 ; iii<200 ; iii++){
      ctrawh[iii][0]=0 ; 
      ct2h[iii][0]=0 ;
      ctrawh[iii][1]=0 ; 
      ct2h[iii][1]=0 ;
    }
  for (int kk=0; kk<nentriesD;  kk++){
      //  for (int kk=0; kk<100;  kk++){

   tt->GetEntry(kk);
   evtType = tt->GetLeaf("fEvtHdr.fEvtType")->GetValue(); 
   if(evtType < 1. )printf("kk,evtype %d %.1f\n",kk,evtType) ;
   if(evtType < 4.) evtypecnr[0]++ ;
   if(evtType == 4.) evtypecnr[1]++ ;
   if(evtType > 4.) evtypecnr[2]++ ;
   icc=(int)ptdcmult ; 
   if(kk<50) printf(" %.1f %d \n",ptdcmult,icc) ;
   if(icc<0) icc=0 ; if(icc>14) icc=14 ;
   ptdcmulth[icc]++ ; 
   // now only analyzing events with SHMS hodoscope
   // refeence time multiplicity up to 5, not more
    if(evtType > 0.9 && icc<6) {
   //if(evtType > 0 && icc<6) {
     npassed++ ;
    for(int iii=0 ; iii<1000 ; iii++){
      hncombo1[iii]=0 ;
      hncombo2[iii]=0 ;
      pncombo1[iii]=0 ;
      pncombo2[iii]=0 ;
    }

    ctpi1 = ctpi1 - 43. ;
    ctpi2 = ctpi2 - 43. ;
    ctp1 = ctp1 - 43. ;
    ctk1 = ctk1 - 43. ;

    Int_t iii = (int) (ctpi2+100.) ;
    //    if(evtType<4) printf("evtype ct %.1f %.1f\n",evtType,ctpi2) ;
    if(ctpi2>-100 && ctpi2<100){
      if(evtType==2) ct2h[iii][0]++ ;
      if(evtType>3) ct2h[iii][1]++ ;
    }
    if((poffsettime > 100 || poffsettime < -100) &&  pdelta<99.) 
      // printf("poff, aero %6.1f %6.1f\n",
      //poffsettime,paeronpe) ;
    if((hoffsettime > 100 || hoffsettime < -100) &&  pdelta<99.) 
      // printf("hoff,cer %6.1f %6.1f\n",hoffsettime,hcernpe) ;

      if(kk<0) printf("%6.1f %6.1f %6.1f\n",pbeta,hbeta,pdelta) ;

    if(kk<0) printf("%6.1f %6.1f %5.0f %5.0f %6.1f %6.1f %5.0f %5.0f \n",
		      htime1,htime2,hhits1,hhits2,
		      ptime1,ptime2,phits1,phits2) ;
    if(kk<0) printf("%d %6.0f %6.0f %6.0f %6.0f \n",pnumcombo1,pncombo1[0],pncombo1[1],pncombo1[2],pncombo1[3]);

    t11 = TcoinpTRIG1_ROC1_tdcTimeRaw ;   
    t31 = TcoinpTRIG3_ROC1_tdcTimeRaw ;
    t41 = TcoinpTRIG4_ROC1_tdcTimeRaw ;
    t51 = TcoinpTRIG5_ROC1_tdcTimeRaw ;  
    t61 = TcoinpTRIG6_ROC1_tdcTimeRaw ;
    t12 = TcoinpTRIG1_ROC2_tdcTimeRaw ;
    t32 = TcoinpTRIG3_ROC2_tdcTimeRaw ; 
    t42 = TcoinpTRIG4_ROC2_tdcTimeRaw ;
    t52 = TcoinpTRIG5_ROC2_tdcTimeRaw ;
    t62 = TcoinpTRIG6_ROC2_tdcTimeRaw ;                                           
    hst = HhodStartTime ; 
    pst = PhodStartTime ; 
    hpst = hst - pst ;
    

    // for when fp time not in root
    //    HhodfpHitsTime = HhodStartTime ; 
    //PhodfpHitsTime = PhodStartTime ; 


// fix ngcer which has negative values sometime
    ngcorr = pcernpeng ;
    if(ngcorr < 0.) ngcorr = 0. ;

// Lokk at aerogel times
    aero_avtime = 0. ; 
    ntime = 0. ;
    if(paeronpe > 0) {
      if(kk<30)printf("%d %6.1f \n",kk,paeronpe) ;
      for(icc=0 ; icc<6 ;  icc++) {
	if(paeroptd[icc] < 9999.) {
	  ntime++ ;
	  aero_avtime += paeroptd[icc] ; 
	  if(kk<30) printf("p %d %6.1f %6.1f \n",icc,
			    paeropt[icc],paeroptd[icc]) ;
	}
      } 
      // neg side seems to have 200 nsec offset
      for(icc=0 ; icc<6 ;  icc++) {
	if(paerontd[icc] < 9999.) {
	  ntime++ ;
	  aero_avtime += paerontd[icc] - 200. ; 
	  if(kk<30) printf("n %d %6.1f %6.1f \n",icc,
			    paeront[icc]+200.,paerontd[icc]-200.) ;
	}
      } 
      if(ntime>0) aero_avtime = aero_avtime / ntime ;
      if(kk<30) printf("avtime= %6.1f %6.1f %6.1f \n",
		       aero_avtime,PhodfpHitsTime,ctpi2) ;
    }

// Get helicity: differene for sping18 and later
    if(runNumber <  5360) {
   	  helmpsi = 0 ;
	  if(helmpsa > 0) helmpsi = 1 ;
	  helnegi = 0 ;
	  if(helnega > 0) helnegi = 1 ;
	  helposi = 0 ;
	  if(helposa > 0) helposi = 1 ;
    }
    if(runNumber >  5360) {
   	  helmpsi = 0 ;
	  if(helmps > 0) helmpsi = 1 ;
	  helnegi = 0 ;
	  if(hel1 < 0) helnegi = 1 ;
	  helposi = 0 ;
	  if(hel1 > 0) helposi = 1 ;
    }

// Don't Fix fttimes that shouldn't be zero
//    if(HhodfpHitsTime == 0.0) HhodfpHitsTime = HhodStartTime ;
//    if(PhodfpHitsTime == 0.0) PhodfpHitsTime = PhodStartTime ;

    tsfdiff = PhodStartTime -PhodfpHitsTime ; 
    icc = (int)(tsfdiff+50.) ;
    if(icc<0) icc=0;
    if(icc>98) icc=98 ;
    tsfh[icc]++ ;
    tsfh[99]++ ;

    tsfdiffh = HhodStartTime -HhodfpHitsTime ; 
    icc = (int)(tsfdiffh+50.) ;
    if(icc<0) icc=0;
    if(icc>98) icc=98 ;
    tsfhh[icc]++ ;
    tsfhh[99]++ ;

// Don't Use Startime instead of FPtime if diff>1 nsec
// unless starttime=32 (the code for "bad")
//    if(tsfdiff< -1. || tsfdiff>1.) {
//      if(PhodStartTime != 32.) {
//	PhodfpHitsTime = PhodStartTime ;
//      }
//    }

      //if(pdelta<9999. && pdelta>-9999 &&
      //   TcoinpTRIG1_ROC2_tdcTimeRaw == 0)
      //      printf("s %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.1f \n",
      // TcoinpTRIG1_ROC2_tdcTimeRaw*0.0997,
      // PhodStartTime,PhodfpHitsTime,tsfdiff,pdelta,
      //	     pbeta,pgtry,pgdsc) ;

      //if(hdelta<9999. && hdelta>-9999 &&
      // TcoinpTRIG4_ROC2_tdcTimeRaw == 0)
      //      printf("h %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.1f \n",
      //       TcoinpTRIG4_ROC2_tdcTimeRaw*0.0997,
      //      HhodStartTime,HhodfpHitsTime,tsfdiffh,hdelta,
      //hbeta,hgtry,hgdsc) ;

// Define raw coincidence time
    nevntall++; 
    if( 
        HhodfpHitsTime == 32.0 || 
        HhodfpHitsTime > 200 ||
        HhodfpHitsTime < 0 || 
        HhodfpHitsTime == 0.0) {
      if(hdelta > 0. && hdelta < 0.01 ) {
	printf("hms 32! %d %d %.2f %.2f %.2f \n",kk,nevntall,HhodStartTime,
	       HhodfpHitsTime,hdelta) ;
	}
      nevnt32h++ ;
    }
    if(HhodfpHitsTime == 32.0) nevntjust32h++ ;
    if(HhodStartTime == 32.0) nevntjust32hst++ ;

    if(
        PhodfpHitsTime == 32.0 || 
        PhodfpHitsTime > 200. || 
        PhodfpHitsTime < 0.  || 
        PhodfpHitsTime == 0.0) {
      if(pdelta > -0.01 && pdelta < 0.) {
	printf("32! %d %.2f %.2f %.2f \n",kk,PhodStartTime,
           PhodfpHitsTime,pdelta) ;
      }
      nevnt32++ ;
    }
    if(PhodfpHitsTime == 32.0) nevntjust32++ ; 
    if(PhodStartTime == 32.0) nevntjust32st++ ; 

    sttoth = sttoth + 1. ;
    sttotp = sttotp + 1. ;
    icc = (int) HhodStartTime + 1.  ;
    if(icc<0) icc=0 ;
    if(icc>99) icc=99 ;
    sthh[icc] = sthh[icc] + 1. ;
    icc = (int) HhodfpHitsTime + 1.  ;
    if(icc<0) icc=0 ;
    if(icc>99) icc=99 ;
    fphh[icc] = fphh[icc] + 1. ;

    icc = (int) PhodStartTime + 1.  ;
    if(icc<0) icc=0 ;
    if(icc>99) icc=99 ;
    sthp[icc] = sthp[icc] + 1. ;
    icc = (int) PhodfpHitsTime + 1.  ;
    if(icc<0) icc=0 ;
    if(icc>99) icc=99 ;
    fphp[icc] = fphp[icc] + 1. ;

// Set flags to define good times. 
// This no longer includes having  a valid track.
    goodsch = 0. ;
    if(HhodfpHitsTime > 0. && 
       HhodfpHitsTime < 199. &&
       //       hdelta>-99. && hdelta<99. &&
       TcoinpTRIG4_ROC2_tdcTimeRaw > 0 &&
       HhodfpHitsTime !=32.0) {
       fpgoodh = fpgoodh + 1.;
       goodsch = 1. ;
    }
    stgoodsch = 0. ;
    if(HhodStartTime > 0. && 
       HhodStartTime < 99. &&
       //       hdelta>-99. && hdelta<99. &&
       TcoinpTRIG4_ROC2_tdcTimeRaw > 0 &&
       HhodStartTime !=32.0) {
       stgoodh = stgoodh + 1.;
       stgoodsch = 1. ;
    }

    goodscp = 0. ; 
    if(PhodfpHitsTime > 0. && 
       PhodfpHitsTime < 99. &&
       //       pdelta>-99. && pdelta<99. &&
       TcoinpTRIG1_ROC2_tdcTimeRaw > 0 &&
       PhodfpHitsTime != 32.00) {
      fpgoodp = fpgoodp + 1.;
      goodscp = 1 ;
    }

    stgoodscp = 0. ; 
    if(PhodStartTime > 0. && 
       PhodStartTime < 99. &&
       //pdelta>-99. && pdelta<99. &&
       TcoinpTRIG1_ROC2_tdcTimeRaw > 0 &&
       PhodStartTime != 32.00) {
      stgoodp = stgoodp + 1.;
      stgoodscp = 1 ;
    }

    SHMScoinCorr = (PhodoStartTimeMean - PhodfpHitsTime); 
    HMScoinCorr = (HhodoStartTimeMean - HhodfpHitsTime);      

    ctimeraw = 72. + 
         (TcoinpTRIG1_ROC2_tdcTimeRaw*0.0997 - 
         SHMScoinCorr) - 
         (TcoinpTRIG4_ROC2_tdcTimeRaw*0.0997 - 
         HMScoinCorr) - pOffset; 
    if(runNumber <  4860) {
	 ctimeraw = ctimeraw + 31. ;
    }

    iii = (int) (ctimeraw + 100.) ;
    if(ctimeraw>-100 && ctimeraw<100){
      if(evtType==2) ctrawh[iii][0]++ ;
      if(evtType>3) ctrawh[iii][1]++ ;
    }


// time histo versus fptime
    icc = (int) ((ctimeraw + 200.) * 0.25 ) ;
    iccc = (int) (PhodfpHitsTime/10.) ;
    if(iccc<0) iccc=0 ;
    if(iccc>9) iccc=9 ;
    if(icc<0) icc=0 ;
    if(icc>98) icc=98 ;
    cthx[icc][iccc]++ ;
    cthx[99][iccc]++ ;

// time histograms
    icc = (int) ((ctpi1  + 200.) * 0.25 ) ;
    if(icc > -1 && icc < 100) cch1[icc]++ ; 
    if(icc<0) cch1[0]++ ;
    if(icc>99) cch1[99]++ ;
    if(hdelta > -100 && hdelta < 100 && 
       pdelta > -100 && pdelta < 100) {
     if(icc > -1 && icc < 100) cch2[icc]++ ; 
     if(icc<0) cch2[0]++ ;
     if(icc>99) cch2[99]++ ;
    }
    icc = (int) ((ctimeraw  + 200. - 60. ) * 0.25 ) ;
    if(icc > -1 && icc < 100) cch3[icc]++ ; 
    if(icc<0) cch3[0]++ ;
    if(icc>99) cch3[99]++ ;
    icc = (int) ((ctimeraw  + 2000. - 60. ) * 0.025 ) ;
    if(icc > -1 && icc < 100) ctbig[icc]++ ; 
    if(icc<0) ctbig[0]++ ;
    if(icc>99) ctbig[99]++ ;

    phbig = hdelta ;
    if(hdelta > 999.) phbig = 999. ;
    if(hdelta < -99.) phbig = -99. ;
    ppbig - pdelta ;
    if(pdelta > 999.) ppbig= 999. ;
    if(pdelta < -99.) ppbig= -99. ;
    t1big = ctpi1 ;
    if(ctpi1<-9999.) t1big = -9999. ;
    if(ctpi1>99999.) t1big = 99999. ;
    t2big = ctimeraw ;
    if(ctimeraw<-9999.) t2big = -9999. ;
    if(ctimeraw>99999.) t2big = 99999. ;
    if(ctimeraw < -2000 || ctimeraw > 2000.) {
      printf(" %8.1f %8.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n", t1big,
	     t2big,phbig,ppbig,paeronpe,pgdsc,hgdsc,hcernpe) ;
    }

    counters[0] = counters[0]+1 ;

    if(hdelta > -20 && hdelta < 20) {
      counters[1] = counters[1]+1 ;
    }

// start long section of wide cuts used to get ctimepi, k, p
    if(hdelta > -20 && hdelta < 20 && 
       pdelta > -100 && pdelta < 100) {
       counters[2] = counters[2]+1 ;

       SHMSpartMass = 0.1395704; // pion mass in GeV/c^2 
       DeltaHMSpathLength = 12.462*HgtrTh + 
         0.1138*HgtrTh*HgtrX - 
         0.0154*HgtrX - 72.292*HgtrTh*HgtrTh - 
         0.0000544*HgtrX*HgtrX - 
         116.52*HgtrPh*HgtrPh;               
       PgtrBetaCalc = PgtrP/
         sqrt(PgtrP*PgtrP + 
         SHMSpartMass*SHMSpartMass);        
       HgtrBetaCalc = HgtrP/
         sqrt(HgtrP*HgtrP + 
         HMSpartMass*HMSpartMass);          
       SHMScoinCorr = SHMScentralPathLen / 
         (speedOfLight*PgtrBetaCalc) + 
         (SHMSpathLength - SHMScentralPathLen) / 
         (speedOfLight*PgtrBetaCalc) + 
          (PhodoStartTimeMean - PhodfpHitsTime); 
       HMScoinCorr = HMScentralPathLen / 
         (speedOfLight*HgtrBetaCalc) + 
         DeltaHMSpathLength / 
         (speedOfLight*HgtrBetaCalc) + 
         (HhodoStartTimeMean - HhodfpHitsTime);      
       SHMScoinCorr = SHMScentralPathLen / 
         (speedOfLight*PgtrBetaCalc) + 
         (SHMSpathLength - SHMScentralPathLen) / 
         (speedOfLight*PgtrBetaCalc) + 
	 (PhodoStartTimeMean - PhodfpHitsTime); 
       HMScoinCorr = HMScentralPathLen / 
         (speedOfLight*HgtrBetaCalc) + 
         DeltaHMSpathLength / 
         (speedOfLight*HgtrBetaCalc) + 
         (HhodoStartTimeMean - HhodfpHitsTime);      
       SHMScorrCoinTimeROC1 = 
         (TcoinpTRIG1_ROC1_tdcTimeRaw*0.0997 - 
         SHMScoinCorr) - 
         (TcoinpTRIG4_ROC1_tdcTimeRaw*0.0997 - 
         HMScoinCorr) - pOffset; // 0.1 to convert to ns 
       SHMScorrCoinTimeROC2 = 
         (TcoinpTRIG1_ROC2_tdcTimeRaw*0.0997 - 
         SHMScoinCorr) - 
         (TcoinpTRIG4_ROC2_tdcTimeRaw*0.0997 - 
         HMScoinCorr) - pOffset; 

// See if this makes a difference
       mjcorr = ((t11 - t12) - (t12 - t42)) * 0.0997 ;
       if(kk < 0) printf("%.2f %.1f %.1f %.1f %.1f \n",
	 mjcorr,t11,t12,t12,t42) ;

// Coin time using fp times
       ctimepi = -999. ;
       if(goodscp == 1 && goodsch==1) 
        ctimepi = SHMScorrCoinTimeROC2 ;

// Coin time using start times
       ctimepinew = -999. ;
       if(goodscp == 1 && goodsch==1) 
        ctimepinew = SHMScorrCoinTimeROC2 -
	  (PhodfpHitsTime - PhodStartTime) +
	  (HhodfpHitsTime - HhodStartTime) ; 

       tdiff1 = PhodfpHitsTime - HhodfpHitsTime ;
       tdiff2 = (t41 - t42)*0.0997 ;
       SHMSpartMass = 0.93827231; // proton mass in GeV/c^2
       DeltaHMSpathLength = 12.462*HgtrTh + 
         0.1138*HgtrTh*HgtrX - 
         0.0154*HgtrX - 72.292*HgtrTh*HgtrTh - 
         0.0000544*HgtrX*HgtrX - 
         116.52*HgtrPh*HgtrPh;               
       PgtrBetaCalc = PgtrP/
         sqrt(PgtrP*PgtrP + 
         SHMSpartMass*SHMSpartMass);        
       HgtrBetaCalc = HgtrP/
         sqrt(HgtrP*HgtrP + 
         HMSpartMass*HMSpartMass);          
       SHMScoinCorr = SHMScentralPathLen / 
         (speedOfLight*PgtrBetaCalc) + 
         (SHMSpathLength - SHMScentralPathLen) / 
         (speedOfLight*PgtrBetaCalc) + 
         (PhodoStartTimeMean - PhodfpHitsTime); 
       HMScoinCorr = HMScentralPathLen / 
         (speedOfLight*HgtrBetaCalc) + 
         DeltaHMSpathLength / 
         (speedOfLight*HgtrBetaCalc) + 
         (HhodoStartTimeMean - HhodfpHitsTime);      
       SHMScoinCorr = SHMScentralPathLen / 
         (speedOfLight*PgtrBetaCalc) + 
         (SHMSpathLength - SHMScentralPathLen) / 
         (speedOfLight*PgtrBetaCalc) + 
         (PhodoStartTimeMean - PhodfpHitsTime); 
       HMScoinCorr = HMScentralPathLen / 
         (speedOfLight*HgtrBetaCalc) + 
         DeltaHMSpathLength / 
         (speedOfLight*HgtrBetaCalc) + 
         (HhodoStartTimeMean - HhodfpHitsTime);      
       SHMScorrCoinTimeROC1 = 
         (TcoinpTRIG1_ROC1_tdcTimeRaw*0.1 - 
         SHMScoinCorr) - 
         (TcoinpTRIG4_ROC1_tdcTimeRaw*0.1 - 
         HMScoinCorr) - pOffset; // 0.1 to convert to ns 
       SHMScorrCoinTimeROC2 = 
         (TcoinpTRIG1_ROC2_tdcTimeRaw*0.1 - 
         SHMScoinCorr) - 
         (TcoinpTRIG4_ROC2_tdcTimeRaw*0.1 - 
         HMScoinCorr) - pOffset; 
       ctimep = SHMScorrCoinTimeROC2 ;

       SHMSpartMass = 0.497648;// kaon mass in GeV/c^2
       DeltaHMSpathLength = 12.462*HgtrTh + 
         0.1138*HgtrTh*HgtrX - 
         0.0154*HgtrX - 72.292*HgtrTh*HgtrTh - 
         0.0000544*HgtrX*HgtrX - 
         116.52*HgtrPh*HgtrPh;               
       PgtrBetaCalc = PgtrP/
         sqrt(PgtrP*PgtrP + 
         SHMSpartMass*SHMSpartMass);        
       HgtrBetaCalc = HgtrP/
         sqrt(HgtrP*HgtrP + 
         HMSpartMass*HMSpartMass);          
       SHMScoinCorr = SHMScentralPathLen / 
         (speedOfLight*PgtrBetaCalc) + 
         (SHMSpathLength - SHMScentralPathLen) / 
         (speedOfLight*PgtrBetaCalc) + 
         (PhodoStartTimeMean - PhodfpHitsTime); 
       HMScoinCorr = HMScentralPathLen / 
         (speedOfLight*HgtrBetaCalc) + 
         DeltaHMSpathLength / 
         (speedOfLight*HgtrBetaCalc) + 
         (HhodoStartTimeMean - HhodfpHitsTime);      
       SHMScoinCorr = SHMScentralPathLen / 
         (speedOfLight*PgtrBetaCalc) + 
         (SHMSpathLength - SHMScentralPathLen) / 
         (speedOfLight*PgtrBetaCalc) + 
         (PhodoStartTimeMean - PhodfpHitsTime); 
       HMScoinCorr = HMScentralPathLen / 
         (speedOfLight*HgtrBetaCalc) + 
         DeltaHMSpathLength / 
         (speedOfLight*HgtrBetaCalc) + 
         (HhodoStartTimeMean - HhodfpHitsTime);      
       SHMScorrCoinTimeROC1 = 
         (TcoinpTRIG1_ROC1_tdcTimeRaw*0.1 - 
         SHMScoinCorr) - 
         (TcoinpTRIG4_ROC1_tdcTimeRaw*0.1 - 
         HMScoinCorr) - pOffset; // 0.1 to convert to ns 
       SHMScorrCoinTimeROC2 = 
         (TcoinpTRIG1_ROC2_tdcTimeRaw*0.1 - 
         SHMScoinCorr) - 
         (TcoinpTRIG4_ROC2_tdcTimeRaw*0.1 - 
         HMScoinCorr) - pOffset; 
       ctimeK = SHMScorrCoinTimeROC2 ;

       if(runNumber <  4860) {
	 ctimepi = ctimepi + 31. ;
	 ctimeK  = ctimeK + 31. ;
	 ctimep  = ctimep + 31. ;
	 ctpi1 = ctpi1 + 32.2 ;
	 ctk1 = ctk1 + 32.2 ;
	 ctp1 = ctp1 + 32.2 ;
	 ctpi2 = ctpi2 + 32.2 ;
       }
       if(t41>280 && t41<330){
	 tmp = (t41 - 280) * 10. ;
         icc = (int) tmp;
         t41h[icc]++ ;
       }
       ctpi3 = hpst - t11 + t12 + 111.; 
       ctpi4 = hpst - t41 + t42 + 111. - 0.7 ; 
       if(t41 <300.) {
	 ctpi3 = ctpi3 - 135. ;
	 ctpi4 = ctpi4 - 135. ;
       }
       if(ctpi2 > -30. && ctpi2 < 30 && kk<400 && kk>3900) 
	 printf(
        "ct %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f \n",
        hpst,t41,ctpi1,ctpi2,ctimepi,ctimeK,ctk1-ctpi1,ctp1-ctpi1,hdelta,pdelta) ; 
       //        (t11 - t41 - (t12 - t42))) ;
       //       printf("%8.2f %8.2f\n",ctimepi,pbeta) ;
       //       if(kk>30 && kk<20000000) 
       if(hdelta > -10 && hdelta < 10. && 
          pdelta > -10 && pdelta < 20. && 
	 paeronpe>.4 && pcaltot>0.10 &&
         hcaletot>0.6 && hcernpe>2 &&
         ctimepi > -10. && ctimepi < 10) { 
	 tmp = (ctimepi + 10.)*5. ;
	 icc = (int) tmp;
	 ctpih[icc]++ ;
       }
      //      ctpi1 = ctpi1 - (t12-391.) ;
      //ctpi2 = ctpi2 - (t12-391.) + 1.0 ;
// Electron in HMS. Do not require hcer spring18
      if(hdelta > -15 && hdelta < 15. && 
         pdelta > -30 && pdelta < 50. && 
         hcaletot>0.6 && 
         (hcernpe>1 || runNumber<4400)) { 

       if(ctimepi > -149.9 && ctimepi < 149.9) { 
	 tmp = (ctimepi + 150.) / 0.25 ;
         icc = (int) tmp;
	 chistww[icc][0]++ ;
       }
       if(ctimepi < -149.9) chistww[0][0]++ ;
       if(ctimepi > 149.9) chistww[1199][0]++ ; 

       if(ctpi2 > -149.9 && ctpi2 < 149.9) { 
	 tmp = (ctpi2 + 150.) / 0.25 ;
         icc = (int) tmp;
	 chistww[icc][1]++ ;
       }
       if(ctpi2 < -149.9) chistww[0][1]++ ;
       if(ctpi2 > 149.9) chistww[1199][1]++ ; 

       if(ctimepi<-150. && kk>400 && kk<100) printf("%8.1f %8.1f %8.1f %8.1f %8.1f\n",ctimepi,
	 ctpi1,ctpi2, hdelta,pdelta) ;

       counters[3] = counters[3]+1 ;

// using ctpi2 for pi time, but ctpi1 for k, p time diff.
// because k2, p2 not in root file
//       if(ctimepi > -30. && ctimepi < 30) { 
       if(ctpi2 > -30. && ctpi2 < 30) { 
	counters[4] = counters[4]+1 ;
	gdtrk = hgdsc * 10. + pgdsc ;
	fprintf(f2,
"%6.2f %5.2f %5.2f %1d %1d %1d %6.2f %7.4f %7.4f %6.2f %7.4f %7.4f %4.1f %4.2f %4.2f %4.1f %4.1f %4.1f %4.2f %4.2f %5.1f %5.1f %6.3f %6.2f %7.4f %6.2f %7.4f %6.2f %7.4f %6.2f %7.4f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %1d %1d %.3f %.2f %.1f %.0f %.0f\n",
//	   ctimepi,ctimeK-ctimepi,ctimep-ctimepi,
	   ctpi2,ctk1-ctpi1,ctp1-ctpi1,
//	    (int)helmpsa,(int)helnega,(int)helposa,hdelta,HgtrTh,
	   (int)hel1+1,(int)helrep+1,(int)helmps+1,hdelta,HgtrTh,
           HgtrPh,pdelta,PgtrTh,PgtrPh,hcernpe,
	   hcalepr,hcaletot,ngcorr,
	   paeronpe,pcernpe,pcalepr,pcaletot,hgtry,pgtry,pbeta,
	   hdcx,hdcxp,hdcy,hdcyp,pdcx,pdcxp,pdcy,pdcyp,prf,hrf,
		PhodfpHitsTime,HhodfpHitsTime,pztar,hztar,gdtrk,
		goodsch,goodscp,hbeta,ctimepinew,aero_avtime,
		ntime,evtType) ;
	 //	 if (kk % 100 == 0) cout << kk*100/nentriesD << "   % of data done " << pdelta << "  " << paeronpe << "  " << pbeta << " " << pcaltot << endl;
         //if(kk<10) printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f f%8.4f \n",
         //  Hgtrxtar,Pgtrxtar) ;
	 if(hcernpe>1 && hcaletot>0.7) {
	   counters[5] = counters[5]+1 ;
	   if(paeronpe > 2. && pcaletot > 0.02 && 
             ctpi2 > -30. && ctpi2 < 30.) {
	     counters[6] = counters[6]+1 ;
	     td3  = PhodfpHitsTime - prf + 1000. ;
	     it3 = (int) (td3/4.) ;
             td4 = td3 - 4.*it3 ;
	     if(kk>400 && evtType==2) {
	       printf("%7d %4.0f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.3f %6.1f %8.1f\n",
		      kk,pgdsc,paeronpe,pcaltot,ctimeraw,
		      ctimepi,ctpi2,pdelta,pbeta,pgtry,PhodStartTime);
	     }
	   }
 	 }
	 //	 printf("%.2f %.2f %.2f %.2f %.2f \n",hel1,helpred,helrep,helmps,helnqrt) ; 
       } // cointime check
      } // delta, PID cut
// Electron in SHMS
      if(hdelta > -12 && hdelta < 12. && 
          pdelta > -15 && pdelta < 24. && 
         pcaltot >0.6 && paeronpe>1.0 && 
         (ngcorr>1.0 || runNumber>5300) && 
         hcaletot>0.05) { 
       if(ctimepi > -2000. && ctimepi < 2600) { 
	counters[7] = counters[7]+1 ;
       }
       if(ctpi2 > -30. && ctpi2 < 30) { 
	counters[8] = counters[8]+1 ;
	fprintf(f3,
"%6.2f %5.2f %5.2f %1d %1d %1d %6.2f %7.4f %7.4f %6.2f %7.4f %7.4f %4.1f %4.2f %4.2f %4.1f %4.1f %4.1f %4.2f %4.2f %5.1f %5.1f %6.3f %6.2f %7.4f %6.2f %7.4f %6.2f %7.4f %6.2f %7.4f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %1d %1d %.3f %.2f %.1f %.0f\n",
//	   ctimepi,ctimeK-ctimepi,ctimep-ctimepi,
	   ctpi2,ctk1-ctpi1,ctp1-ctpi1,
//	    (int)helmpsa,(int)helnega,(int)helposa,hdelta,HgtrTh,
	   (int)hel1+1,(int)helrep+1,(int)helmps+1,hdelta,HgtrTh,
           HgtrPh,pdelta,PgtrTh,PgtrPh,hcernpe,
	   hcalepr,hcaletot,ngcorr,
	   paeronpe,pcernpe,pcalepr,pcaletot,hgtry,pgtry,pbeta,
	   hdcx,hdcxp,hdcy,hdcyp,pdcx,pdcxp,pdcy,pdcyp,prf,hrf,
		PhodfpHitsTime,HhodfpHitsTime,pztar,hztar,gdtrk,
		goodsch,goodscp,hbeta,ctimepinew,aero_avtime,ntime) ;
       } // cointime check
      } // delta, PID cut
    } //wide pdelta check
    // get SHMS tracking efficiency cut on coin with good e
    // to do add cut on coin time
    // timing with no track info
       //       HMScoinCorr = HMScentralPathLen / 
       //  (speedOfLight*HgtrBetaCalc) + 
       //  DeltaHMSpathLength / 
       //  (speedOfLight*HgtrBetaCalc) + 
       //  (HhodoStartTimeMean - HhodfpHitsTime);      
      //printf("%7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f   \n",ctimeraw,hdelta,pgdsc,paeronpe,pcaltot,hcaletot,hcernpe) ;
    
// get HMS tracking efficiency for HMS
    if(pdelta > -10 && pdelta < 20 && 
       paeronpe > 2.5 && pcaltot > 0.10
       && hcernpe > 0.1 && hcaltot>0.4) {
      hncoin = hncoin + 1. ;
      if(hgdsc>0. && goodsch) { 
       hncoinf = hncoinf + 1. ;
       if(hdelta>-20. && hdelta<20.) {
         hncoinwtrk = hncoinwtrk + 1 ;
 	if(hbeta>0.8 && hbeta<1.2) {
          hncoinwtrk1 = hncoinwtrk1 + 1 ;
 	 if(hgtry>-50 && hgtry< 50.) {
 	   hncoinwtrk2 = hncoinwtrk2 + 1 ;
 	 }
 	}
       }
      }
// get efficiency of hgdsc
      if(hdelta >-5 && hdelta < 5 &&
           hgtry > -2 && hgtry < 2.) {
	  hngdsc1 = hngdsc1 +1 ;
	  if(hgdsc > 0.) hngdsc2 = hngdsc2 + 1 ;
      }
    }
// get SHMS tracking efficiency and other stuff
    if(hdelta > -10 && hdelta < 10 && 
       paeronpe > 0.4 && pcaltot>0.01 &&
       hcaletot > 0.6 && hcernpe > 2) {

      if(ctimeraw > 0 && ctimeraw < 200.) {
	ncoin = ncoin + 1. ;
	if(pgdsc > 0. && goodscp) { 
	  ncoinf = ncoinf + 1. ;
	  if(kk<400) {
	    printf("%5.1f %5.1f %5.1f %5.1f %5.3f %5.1f\n",
		   paeronpe,pcaltot,ctimeraw,
		   pdelta,pbeta,pgtry);
	  }
	  if(pdelta>-30. && pdelta<60.) {
	    ncoinwtrk = ncoinwtrk + 1 ; 
	    if(pbeta>0.8 && pbeta<1.2) {
	      ncoinwtrk1 = ncoinwtrk1 + 1 ;
	      if(pgtry>-50. && pgtry<50.) {
		ncoinwtrk2 = ncoinwtrk2 + 1 ;
	      }
	    }
	  }
	}
      }
// get efficiency of pgdsc
      if(ctimeraw > 0 && ctimeraw < 200.) {
	if(pdelta >-5 && pdelta < 5 &&
           pgtry > -5. && pgtry < 5.) {
	  ngdsc1 = ngdsc1 +1 ;
	  if(pgdsc > 0.) ngdsc2 = ngdsc2 + 1 ;
	}
      } // test on cointimeraw

      icc = (int) ctimeraw -10 ;
      if(icc<0) chist[0]++ ;
      if(icc>99) chist[99]++ ;
      if(icc > -1 && icc < 100) {
       chist[icc]++ ; 
       if(pdelta>-100. && pdelta<100.) chistp[icc]++ ;
       if(pgdsc>0.) {
        chistg[icc]++ ;
   	if(pdelta>-100. && pdelta<100.) chistpg[icc]++ ;
       }
      }
      diff1 = t11 - t12 ;
      diff3 = t31 - t32 ; 
      diff4 = t41 - t42 ;
      diff5 = t51 - t52 ;
      diff6 = t61 - t62 ;

      if(ctimepi > -4.0 && ctimepi < 0.0 && pbeta>0.95 && pbeta<1.05) {
	tmp = 100.*(pbeta-0.949) ;
	ibeta = (int) tmp ;
        tmp = (ctimepi + 4.) * 5. ;
	icc = (int) tmp ;
	if(icc > -1 && icc<20 && ibeta > -1 && ibeta<10) 
          chistb[icc][ibeta]++ ;
        nevnt++ ; 
	/*        fprintf(f4,"%6d %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n",
           nevnt, ctimepi, 
		 SHMScoinCorr,HMScoinCorr, -0.1*diff1, hbeta, pbeta,  
		 diff5 - diff1, diff6-diff1 ) ; */
      } 
      tmp = (ctimepi + 5.) * 10. ;
      icc = (int) tmp ;
      if(icc > -1 && icc < 100) chistt[icc]++ ;
      if(icc > -1 && icc < 100 && pdelta>-100. && pdelta<100.) chistpt[icc]++ ;
    }
   } //test if evttpe==4 AND tdcmult <5
  } // loop over events

  fprintf(f5,"%.0f %.0f %.3f \n",fpgoodh,sttoth,
	  fpgoodh/sttoth) ;
  fprintf(f5,"%.0f %.0f %.3f \n",fpgoodp,sttotp,
	  fpgoodp/sttotp) ;
  fprintf(f5,"%d %d %.3f %d %d %.3f %d %d\n",
          nevnt32, nevntall,(double)nevnt32/(double)nevntall,
          nevnt32h,nevntall,(double)nevnt32h/(double)nevntall,
          npassed,nentriesD) ;
  fprintf(f5,"%8.0f %8.0f %8.0f %8.0f %8.0f\n",
	  hncoin,hncoinf,hncoinwtrk,hncoinwtrk1,hncoinwtrk2) ;
  fprintf(f5,"%8.0f %8.0f %8.0f %8.0f %8.0f\n",
	  ncoin,ncoinf,ncoinwtrk,ncoinwtrk1,ncoinwtrk2) ;
  fprintf(f5,"%7.0f %7.0f %7.3f \n",hngdsc1,hngdsc2,hngdsc2/hngdsc1) ;
  fprintf(f5,"%7.0f %7.0f %7.3f \n",ngdsc1,ngdsc2,ngdsc2/ngdsc1) ;
  fprintf(f5,"%7.3f %7.3f %7.3f %7.3f \n",sum2/sum1,sum4/sum3,sum6/sum5,sum8/sum7) ;
  printf("evtypes %d %d %d \n",evtypecnr[0],evtypecnr[1],evtypecnr[2]) ;

  printf("%8.0f %8.3f %8.3f %8.3f %8.3f\n",
	 ncoin,ncoinf/ncoin,ncoinwtrk/ncoinf,ncoinwtrk1/ncoinf,ncoinwtrk2/ncoinf) ;
  printf("%8.0f %8.3f %8.3f %8.3f %8.3f\n",
	 hncoin,hncoinf/hncoin,hncoinwtrk/hncoinf,
         hncoinwtrk1/hncoinf,hncoinwtrk2/hncoinf) ;

  fprintf(f2,"  0.00  0.00 0.00\n") ;
  fprintf(f2,"%8.0f %8.0f %8.0f %8.0f %8.0f\n",
	  ncoin,ncoinf,ncoinwtrk,ncoinwtrk1,ncoinwtrk2) ;
  fprintf(f2,"%8.0f %8.0f %8.0f %8.0f %8.0f\n",
	  hncoin,hncoinf,hncoinwtrk,hncoinwtrk1,hncoinwtrk2) ;

  for(icc=0 ; icc<100 ;  icc++) {
    //    printf("%5d %5d %6d %6d %6d %6d \n",icc,chist[icc],chistp[icc],
    //	   chistpt[icc],ctpih[icc],ctbig[icc]) ;
    trkeff = 0. ;
    trkeffer = 0. ; 
    if(chistp[icc]>2 && chist[icc]>2) {
      trkeff = (float)chistp[icc] / (float)chist[icc] ;
     trkeffer = sqrt(1./chistp[icc] - 1./chist[icc]) ;
    }
    trkeffg = 0. ;
    trkeffger = 0. ; 
    if(chistpg[icc]>2 && chistg[icc]>2) {
      trkeffg = (float)chistpg[icc] / (float)chistg[icc] ;
     trkeffger = sqrt(1./chistpg[icc] - 1./chistg[icc]) ;
    }
    trkefft = 0. ;
    trkeffter = 0. ; 
    if(chistpt[icc]>2 && chistt[icc]>2) {
      trkefft = (float)chistpt[icc] / (float)chistt[icc] ;
     trkeffter = sqrt(1./chistpt[icc] - 1./chistt[icc]) ;
    }
    fprintf(f2," %5d %5d %6d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
      icc,chist[icc],chistp[icc],trkeff,trkeffer,trkeffg,
      trkeffger,trkefft,trkeffter) ;
    fprintf(f5," %5d %5d %6d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
      icc,chist[icc],chistp[icc],trkeff,trkeffer,trkeffg,
      trkeffger,trkefft,trkeffter) ;
    //printf(" %5d %5d %6d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
    //  icc,chist[icc],chistp[icc],trkeff,trkeffer,trkeffg,
    //  trkeffger,trkefft,trkeffter) ;
  }
  for(ibeta=0 ; ibeta<10 ;  ibeta++) {
    sum1 = 0. ;
    sum2 = 0. ; 
    for(icc=0 ; icc<20 ;  icc++) {
      //fprintf(f4,"%5d %5d\n",icc,chistb[icc][ibeta]) ;
      sum1 = sum1 + chistb[icc][ibeta] ;
      sum2 = sum2 + chistb[icc][ibeta] * icc ;
    }
    //fprintf(f4," hist") ;
    //printf(" ibeta, petak= ,%5d %7.2f \n",ibeta,sum2/sum1) ;
  }
  sum1 = 0. ;
  sum2 = 0. ;
  sum3 = 0. ;
  sum4 = 0. ; 
  sum5 = 0. ;
  sum6 = 0. ; 
  sum7 = 0. ;
  sum8 = 0. ; 

  for(icc=0 ; icc<100 ;  icc++) {
    sum1 = sum1 + cch1[icc] ;
    if(icc==0 || icc==99) sum2 = sum2 + cch1[icc] ;
    sum3 = sum3 + cch2[icc] ;
    if(icc==0 || icc==99) sum4 = sum4 + cch2[icc] ;
    sum5 = sum5 + cch3[icc] ;
    if(icc==0 || icc==99) sum6 = sum6 + cch3[icc] ;
    sum7 = sum7 + ctbig[icc] ;
    if(icc>44 && icc<56) sum8 = sum8 + ctbig[icc] ;
    fprintf(f2," %3d %6d %6d %6d %6d %6d\n",
       icc,cch1[icc],cch2[icc],cch3[icc],chistt[icc],ctbig[icc]) ; 
    fprintf(f5," %3d %6d %6d %6d %6d %6d\n",
       icc,cch1[icc],cch2[icc],cch3[icc],chistt[icc],ctbig[icc]) ; 
    //    printf(" %3d %6d %6d %6d %6d %6d\n",
    //   icc,cch1[icc],cch2[icc],cch3[icc],chistt[icc],ctbig[icc]) ; 
  }
  fprintf(f2,"%7.3f %7.3f %7.3f %7.3f \n",sum2/sum1,sum4/sum3,sum6/sum5,sum8/sum7) ;
  printf(   "%7.3f %7.3f %7.3f %7.3f \n",sum2/sum1,sum4/sum3,sum6/sum5,sum8/sum7) ;
  fprintf(f2,"%7.0f %7.0f %7.3f \n",ngdsc1,ngdsc2,ngdsc2/ngdsc1) ;
  fprintf(f5,"%7.0f %7.0f %7.3f \n",ngdsc1,ngdsc2,ngdsc2/ngdsc1) ;
  printf("%7.0f %7.0f %7.3f \n",ngdsc1,ngdsc2,ngdsc2/ngdsc1) ;
  fprintf(f2,"%7.0f %7.0f %7.3f \n",hngdsc1,hngdsc2,hngdsc2/hngdsc1) ;
  printf("%7.0f %7.0f %7.3f \n",hngdsc1,hngdsc2,hngdsc2/hngdsc1) ;

  for(icc=0 ; icc<0 ;  icc++) {
    printf(" %8.2f %5d\n",
	   280.+0.1*icc,t41h[icc]) ; 
    fprintf(f5," %8.2f %5d\n",
	   280.+0.1*icc,t41h[icc]) ; 
  }
  printf("shms32 %d %d %.3f %.3f %.3f\n",nevnt32,nevntall,
	 (double)nevnt32/(double)nevntall,
	 (double)nevntjust32/(double)nevntall,
	 (double)nevntjust32st/(double)nevntall) ;
  printf("hms32 %d %d %.3f %.3f %.3f\n",nevnt32h,nevntall,
	 (double)nevnt32h/(double)nevntall,
	 (double)nevntjust32h/(double)nevntall,
	 (double)nevntjust32hst/(double)nevntall) ;

  for(icc=0 ; icc<102 ;  icc++) {
    fprintf(f5," %5.1f %7.1f %7.1f %7.1f %7.1f \n",
	   -1. + (float)icc, 
	   100.*sthh[icc]/sttoth,100.*fphh[icc]/sttoth,
	   100.*sthp[icc]/sttotp,100.*fphp[icc]/sttotp) ; 
  }

  for(iccc=0 ; iccc<10 ;  iccc++) {
    if(cthx[99][iccc]<1.) cthx[99][iccc]=1. ;
  }
  for(icc=0 ; icc<100 ;  icc++) {
    fprintf(f5,"%5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f  \n",
	    100.*cthx[icc][0]/cthx[99][0],
	    100.*cthx[icc][1]/cthx[99][1],
	    100.*cthx[icc][2]/cthx[99][2],
	    100.*cthx[icc][3]/cthx[99][3],
	    100.*cthx[icc][4]/cthx[99][4],
	    100.*cthx[icc][5]/cthx[99][5],
	    100.*cthx[icc][6]/cthx[99][6],
	    100.*cthx[icc][7]/cthx[99][7],
	    100.*cthx[icc][8]/cthx[99][8],
	    100.*cthx[icc][9]/cthx[99][9]) ;
  }
  for(icc=0 ; icc<100 ;  icc++) {
    fprintf(f5,"%d %5.1f %5.1f \n",icc,
	    100.*tsfhh[icc]/tsfhh[99],
	    100.*tsfh[icc]/tsfh[99]) ;
  }
  for(icc=0 ; icc<200 ;  icc++) {
    fprintf(f5,"%3d %7d %7d %7d %7d  \n",icc,
	    ct2h[icc][0],ct2h[icc][1],
	    ctrawh[icc][0],ctrawh[icc][1]) ;
  }
  for(icc=0 ; icc<200 ;  icc++) {
    fprintf(f5,"%4d %7d %7d %7d %7d  \n",icc-100,
	    ct2h[icc][0],ct2h[icc][1],
	    ctrawh[icc][0],ctrawh[icc][1]) ;
  }

  for(icc=0 ; icc<9 ; icc++){
    printf(" %7.0f %7.3f %d %.3f\n",counters[icc],
	   counters[icc]/counters[0],ptdcmulth[icc],
           (float)ptdcmulth[icc]/(float)nentriesD) ;
  }
  printf("%d %d %.3f %d %d %.3f %d %d\n",
          nevnt32, nevntall,(double)nevnt32/(double)nevntall,
          nevnt32h,nevntall,(double)nevnt32h/(double)nevntall,
          npassed,nentriesD) ;

  fprintf(f6,"set device postscript ; set scale y log \n") ;
  fprintf(f6,"set window x 1 of 2\n") ;
  for(icc=0 ; icc<1199 ; icc++){
    if(chistww[icc][0]>0) fprintf(f6," %7.2f %d\n",
	-150. + 0.25*icc, chistww[icc][0]) ;
  }
  fprintf(f6,"hist\n") ;
  fprintf(f6,"set window x 2 of 2\n") ;
  for(icc=0 ; icc<1199 ; icc++){
    if(chistww[icc][1]>0) fprintf(f6," %7.2f %d\n",
	-150. + 0.25*icc, chistww[icc][1]) ;
  }
  fprintf(f6,"hist\n") ;

  fclose(f2) ;
  //fclose(f4) ;
  fclose(f5) ;
  fclose(f6) ;
 
  } else {
      printf("--------\n") ;
      printf("run %d not FOUND\n",runNumber) ;
      printf("--------\n") ;

    gSystem->Exit(1) ;
  } 

  gROOT
->SetBatch(kFALSE);

  gSystem->Exit(1) ; 
}

