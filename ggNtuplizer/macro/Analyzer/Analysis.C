#define Analysis_cxx
#include "Analysis.h"
#include <TH2.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
 
void Analysis::Loop(int whichSample, int maxEvents,double weight)
{

  //create a new output file
  std::cout<<"-------> "<<whichSample<<std::endl;
  TFile* fout;
  if(whichSample==0) fout=new TFile("fout_DM_Hgg.root", "RECREATE");
  if(whichSample==1) fout=new TFile("fout_ggF_Hgg.root", "RECREATE");
  if (whichSample==2) fout=new TFile("fout_VH_Hgg.root", "RECREATE");
  if (whichSample==3) fout=new TFile("fout_VBF_Hgg.root", "RECREATE");
  if (whichSample==4) fout=new TFile("fout_ttH_Hgg.root", "RECREATE");
  if(whichSample==5) fout=new TFile("fout_Background_QCD.root", "RECREATE");
  
  //book histos
  TH1F* h_phoE_1 = new TH1F("h_phoE_1", "", 120, 0, 300);
  h_phoE_1->Sumw2();
  TH1F* h_phoEt_1 = new TH1F("h_phoEt_1", "", 120, 0, 300);
  h_phoEt_1->Sumw2();
  TH1F* h_phoEta_1 = new TH1F("h_phoEta_1", "", 60, -3, 3);
  h_phoEta_1->Sumw2();
  TH1F* h_phoSieieEB_Matched_1 = new TH1F("h_phoSieieEB_Matched_1", "", 60, 0, 0.1);
  TH1F* h_phoSieieEB_NotMatched_1 = new TH1F("h_phoSieieEB_NotMatched_1", "", 60, 0, 0.1);
  h_phoSieieEB_Matched_1->Sumw2();
  h_phoSieieEB_NotMatched_1->Sumw2();
  
  //diphoton mass
  TH1F* h_diphotonMass= new TH1F("h_diphotonMass", "", 60, 105, 165);
  h_diphotonMass->Sumw2();
  
  //energy corrections
  TH1F* h_RefinedSCEnergyCorrection_1 = new TH1F("h_RefinedSCEnergyCorrection_1", "", 120,0,5);
  h_RefinedSCEnergyCorrection_1->Sumw2();
  TH2F* h_RefinedSCEnergyCorrectionVsPt_1 = new TH2F("h_RefinedSCEnergyCorrectionVsPt_1", "", 120,0,5, 120, 0, 200);
  h_RefinedSCEnergyCorrectionVsPt_1->Sumw2();
  TH1F* h_ResidualEnergyCorrection_1 = new TH1F("h_ResidualEnergyCorrection_1", "", 120,0,5);
  h_ResidualEnergyCorrection_1->Sumw2();
  TH1F* h_diphotonMassCorr= new TH1F("h_diphotonMassCorr", "", 60, 105, 165);
  h_diphotonMassCorr->Sumw2();

  //isolations - Effective areas study
  TH1F* h_phoPFPhoIso_1 = new TH1F("h_phoPFPhoIso_1", "", 60, -1.,11);
  h_phoPFPhoIso_1->Sumw2();
  TH2F* h_phoPFPhoIsoVsNvtx_1 = new TH2F("h_phoPFPhoIsoVsNvtx_1", "", 60, -1.,11, 50, 0,50);
  h_phoPFPhoIsoVsNvtx_1->Sumw2();
  TH2F* h_rhoVsNvtx = new TH2F("h_rhoVsNvtx_1", "", 120, 0.,500, 50, 0,50);
  h_rhoVsNvtx->Sumw2(); 
  TH1F* h_phoCorrPFPhoIso_1 = new TH1F("h_phoCorrPFPhoIso_1", "", 60, -1.,11);
  h_phoCorrPFPhoIso_1->Sumw2();
  TH2F* h_phoCorrPFPhoIsoVsNvtx_1 = new TH2F("h_phoCorrPFPhoIsoVsNvtx_1", "", 60, -1.,11, 50, 0,50);
  h_phoCorrPFPhoIsoVsNvtx_1->Sumw2();

  TH1F h_passCorrPFPhoIso_nvtx;
  h_passCorrPFPhoIso_nvtx.SetBins(50, 0,50);
  h_passCorrPFPhoIso_nvtx.SetTitle("h_passCorrPFPhoIso_nvtx");
  h_passCorrPFPhoIso_nvtx.SetName("h_passCorrPFPhoIso_nvtx");
  h_passCorrPFPhoIso_nvtx.Sumw2();
  TH1F h_passPFPhoIso_nvtx;
  h_passPFPhoIso_nvtx.SetBins(50, 0,50);
  h_passPFPhoIso_nvtx.SetName("h_passPFPhoIso_nvtx");
  h_passPFPhoIso_nvtx.SetTitle("h_passPFPhoIso_nvtx");
  h_passPFPhoIso_nvtx.Sumw2();
  TH1F h_nvtx;
  h_nvtx.SetBins(50, 0,50);
  h_nvtx.SetName("h_nvtx");
  h_nvtx.SetTitle("h_nvtx");
  h_nvtx.Sumw2();

  //diphoton mass with different IDs
  TH1F* h_diphotonMass_LooseID= new TH1F("h_diphotonMass_LooseID", "", 60, 105, 165);
  h_diphotonMass_LooseID->Sumw2();
  TH1F* h_diphotonMass_MediumID= new TH1F("h_diphotonMass_MediumID", "", 60, 105, 165);
  h_diphotonMass_MediumID->Sumw2();
  TH1F* h_diphotonMass_TightID= new TH1F("h_diphotonMass_TightID", "", 60, 105, 165);
  h_diphotonMass_TightID->Sumw2();

  //DM analysis cut flow
  TH1F* h_DiPhotonMassCorrWeight = new  TH1F("h_DiPhotonMassCorrWeight","",60,105,165);
  h_DiPhotonMassCorrWeight->Sumw2();
  TH1F* h_DiPhotonMassCorrWeight_KinCuts = new TH1F("h_DiPhotonMassCorrWeight_KinCuts","",60,105,165);
  h_DiPhotonMassCorrWeight_KinCuts->Sumw2();
  TH1F* h_DiPhotonMassCorrWeight_KinCuts_MedID = new TH1F("h_DiPhotonMassCorrWeight_KinCuts_MedID","",60,105,165);
  h_DiPhotonMassCorrWeight_KinCuts_MedID->Sumw2();
  TH1F* h_DiPhotonMassCorrWeight_KinCuts_MedID_METCut = new TH1F("h_DiPhotonMassCorrWeight_KinCuts_MedID_METCut","",60,105,165);
  h_DiPhotonMassCorrWeight_KinCuts_MedID_METCut->Sumw2();

  TH1F* h_METWeight =new  TH1F("h_METWeight","",120,0,600);
  h_METWeight->Sumw2();
  TH1F* h_METWeight_KinCuts = new TH1F("h_METWeight_KinCuts","",120,0,600);
  h_METWeight_KinCuts->Sumw2();
  TH1F* h_METWeight_KinCuts_MedID = new TH1F("h_METWeight_KinCuts_MedID","",120,0,600);
  h_METWeight_KinCuts_MedID->Sumw2();
  TH1F* h_METWeight_KinCuts_MedID_METCut = new TH1F("h_METWeight_KinCuts_MedID_METCut","",120,0,600);
  h_METWeight_KinCuts_MedID_METCut->Sumw2();

  
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   //-------------------------- LOOP ------------------------------------------------//
   //   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   if(maxEvents>0) nentries=maxEvents;
 
  for (Long64_t jentry=0; jentry<nentries;jentry++) { //loop on only 1k events
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(ientry%100==0) std::cout<<ientry<<std::endl;
      // if (Cut(ientry) < 0) continue;
    
      //look only events with at least two photon
      double size=phoE->size();
      if(size<2)continue;

      
      //-------------------------- Ex: 1.1 ADD OTHER KINEMATIC DISTRIBUTIONS -----------------------------------------------//

      //-------------------------- Ex: 1.2 IMPLEMENT MC MATCHING AND PLOT VARIABLES FOR (UN)MATCHED PHOTONS ----------------//
      
      //-------------------------- Ex: 1.3 COMPUTE DIPHOTON MASS -----------------------------------------------------------//   


      //-------------------------- Ex: 2.1 PLOT ENERGY CORRECTIONS AS A FUNCTION OF ENERGY ---------------------------------//  

      //-------------------------- Ex: 2.2 EFFECT OF ENERGY CORRECTIONS ON HIGGS PEAK --------------------------------------//  



      //-------------------------- Ex: 3.1 CORRECT PHOTON ISOLATION WITH EFFECTIVE AREAS -----------------------------------//

      //-------------------------- Ex: 3.3 APPLY PHOTON ID USING BOTH HAND IMPLEMENTATION AND VID --------------------------//

      //-------------------------- Ex: 3.4 LOOK AT DIPHOTON MASS  WITH DIFFERENT IDs ---------------------------------------//                                               


      //-------------------------- Ex: 5.2 - 5.4 Perform DM-like Analysis --------------------------------------------------//


      //-------------------------- FILL OTHER HISTOS ---------------------------------------------------------//
      
   }


      //------------------------- COMPUTE ISOLATION EFFICIENCY VS PU W/ AND W/O CORRECTIONS-------------------//



   
   //save output
   fout->cd();
   h_phoE_1->Write();
   h_phoEta_1->Write();
   h_phoSieieEB_Matched_1->Write();
   h_phoSieieEB_NotMatched_1->Write();
   h_diphotonMass->Write();
   h_RefinedSCEnergyCorrection_1->Write();
   h_RefinedSCEnergyCorrectionVsPt_1->Write();
   h_ResidualEnergyCorrection_1->Write();
   h_diphotonMassCorr->Write();
   h_phoPFPhoIsoVsNvtx_1->Write();
   h_rhoVsNvtx->Write();
   h_phoPFPhoIso_1->Write();
   h_phoCorrPFPhoIsoVsNvtx_1->Write();
   h_phoCorrPFPhoIso_1->Write();
   h_nvtx.Write();
   h_passPFPhoIso_nvtx.Write();
   h_passCorrPFPhoIso_nvtx.Write();
   //   PFPhoIso_eff->Write();
   //   CorrPFPhoIso_eff->Write();
   h_diphotonMass_LooseID->Write();
   h_diphotonMass_MediumID->Write();
   h_diphotonMass_TightID->Write();
   h_DiPhotonMassCorrWeight->Write();
   h_DiPhotonMassCorrWeight_KinCuts->Write();
   h_DiPhotonMassCorrWeight_KinCuts_MedID->Write();
   h_DiPhotonMassCorrWeight_KinCuts_MedID_METCut->Write();
   h_METWeight->Write();
   h_METWeight_KinCuts->Write();
   h_METWeight_KinCuts_MedID->Write();
   h_METWeight_KinCuts_MedID_METCut->Write();
   fout->Write();
   fout->Close();
}
