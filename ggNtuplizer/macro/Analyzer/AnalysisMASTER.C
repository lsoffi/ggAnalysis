#define AnalysisMASTER_cxx
#include "AnalysisMASTER.h"
#include <TH2.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
 
void AnalysisMASTER::Loop(int whichSample, int maxEvents,double weight)
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

  TH1F* h_diphotonMass= new TH1F("h_diphotonMass", "", 60, 105, 165);
  h_diphotonMass->Sumw2();
  TH1F* h_RefinedSCEnergyCorrection_1 = new TH1F("h_RefinedSCEnergyCorrection_1", "", 120,0,5);
  h_RefinedSCEnergyCorrection_1->Sumw2();
  TH2F* h_RefinedSCEnergyCorrectionVsPt_1 = new TH2F("h_RefinedSCEnergyCorrectionVsPt_1", "", 120,0,5, 120, 0, 200);
  h_RefinedSCEnergyCorrectionVsPt_1->Sumw2();
  TH1F* h_ResidualEnergyCorrection_1 = new TH1F("h_ResidualEnergyCorrection_1", "", 120,0,5);
  h_ResidualEnergyCorrection_1->Sumw2();
  TH1F* h_diphotonMassCorr= new TH1F("h_diphotonMassCorr", "", 60, 105, 165);
  h_diphotonMassCorr->Sumw2();

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

  TH1F* h_diphotonMass_LooseID= new TH1F("h_diphotonMass_LooseID", "", 60, 105, 165);
  h_diphotonMass_LooseID->Sumw2();
  TH1F* h_diphotonMass_MediumID= new TH1F("h_diphotonMass_MediumID", "", 60, 105, 165);
  h_diphotonMass_MediumID->Sumw2();
  TH1F* h_diphotonMass_TightID= new TH1F("h_diphotonMass_TightID", "", 60, 105, 165);
  h_diphotonMass_TightID->Sumw2();

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

      
      //-------------------------- LOOK AT TWO HIGH PT PHOTONS ------------------------------------------------//
      //save the two highest pt photons
      int pho1_index;
      int pho2_index;
      double maxPt_1st=0;
      double maxPt_2nd=0;
      
      //look at first photon
      for(int ipho=0; ipho<phoE->size();ipho++){
	if((*phoEt)[ipho]>maxPt_1st){
	  maxPt_1st= (*phoEt)[ipho];
	  pho1_index=ipho;
	}
      } 
      
      //look at second photon                                                                                                                           
      for(int ipho=0; ipho<phoE->size();ipho++){
	if(ipho==pho1_index)continue;
	if((*phoEt)[ipho]>maxPt_2nd){
          maxPt_2nd= (*phoEt)[ipho];
	  pho2_index=ipho;
        }
      }
  
      int pho1_genIndex=999;
      int pho2_genIndex=999;
      bool pho1_isMatch=false;
      bool pho2_isMatch=false;

      //Loop on gen particles and check if the photon is matched
      //1st photon
      for(int igen=0;igen<nMC; igen++){
	//if is a photon
	if((*mcPID)[igen]!=22)continue;
	//if is prompt
	if((*mcStatusFlag)[igen]>2)continue;

	//compute DeltaR between recon and gen photon
	double deltaR = sqrt(pow((*phoEta)[pho1_index]-(*mcEta)[igen],2)+ pow((*phoPhi)[pho1_index]-(*mcPhi)[igen],2));
	if(deltaR<0.3){
	  pho1_isMatch = true;
	  pho1_genIndex=igen;
	}
      }

      //2nd photon
      for(int igen=0;igen<nMC; igen++){
        if((*mcPID)[igen]!=22)continue;
	if((*mcStatusFlag)[igen]>2)continue;
        double deltaR = sqrt(pow((*phoEta)[pho2_index]-(*mcEta)[igen],2)+ pow((*phoPhi)[pho2_index]-(*mcPhi)[igen],2));
        if(deltaR<0.3){
	  pho2_genIndex=igen;
	  pho2_isMatch = true;
      }

   }

      //-------------------------- COMPUTE DIPHOTON MASS ------------------------------------------------//   

      TLorentzVector pho1;
      TLorentzVector pho2;
      pho1.SetPtEtaPhiM((*phoEt)[pho1_index], (*phoEta)[pho1_index],(*phoPhi)[pho1_index],0);
      pho2.SetPtEtaPhiM((*phoEt)[pho2_index], (*phoEta)[pho2_index],(*phoPhi)[pho2_index],0);
      double diphotonMass = (pho1+pho2).M();
      h_diphotonMass->Fill(diphotonMass);


      //-------------------------- CORRECTED ENERGIES ------------------------------------------------//  

      //plot energy corrections
      double RefinedSCEnergyCorrection = (*phoSCE)[pho1_index]/(*phoSCRawE)[pho1_index];
      double ResidualEnergyCorrection = (*phoCalibEt)[pho1_index]/(*phoSCE)[pho1_index];

      //plot diphoton mass w/ corrected energies
      TLorentzVector phoCorr1;
      TLorentzVector phoCorr2;
      phoCorr1.SetPtEtaPhiM((*phoCalibEt)[pho1_index], (*phoEta)[pho1_index],(*phoPhi)[pho1_index],0);
      phoCorr2.SetPtEtaPhiM((*phoCalibEt)[pho2_index], (*phoEta)[pho2_index],(*phoPhi)[pho2_index],0);
      double diphotonMassCorr = (phoCorr1+phoCorr2).M();


      h_RefinedSCEnergyCorrection_1->Fill(RefinedSCEnergyCorrection);
      h_RefinedSCEnergyCorrectionVsPt_1->Fill(RefinedSCEnergyCorrection, (*phoEt)[pho1_index]);
      h_ResidualEnergyCorrection_1->Fill(ResidualEnergyCorrection);
      h_diphotonMassCorr->Fill(diphotonMassCorr);


      //-------------------------- ISOLATION CORRECTION WITH EFFECTIVE AREA--------------------//
      //consider only photon isolation for central barrel photons

      h_phoPFPhoIso_1->Fill((*phoPFPhoIso)[pho1_index]);
      h_phoPFPhoIsoVsNvtx_1->Fill((*phoPFPhoIso)[pho1_index],nVtx);
      h_rhoVsNvtx->Fill(rho,nVtx);

      double effArea = 0.121;
      double corrPFPhoIso= (*phoPFPhoIso)[pho1_index]- rho*effArea;
      //correct photon isolation
      h_phoCorrPFPhoIso_1->Fill(corrPFPhoIso);
      h_phoCorrPFPhoIsoVsNvtx_1->Fill(corrPFPhoIso,nVtx);
      
      //consider only photon isolation and fill histos needed to compute the efficiency of the selection on this variables vs nvtx (only EB)
      if(abs((*phoEta)[pho1_index])<1.)h_nvtx.Fill(nVtx);
      bool passPFPhoIso = (*phoPFPhoIso)[pho1_index]<(2.362+0.0047*(*phoEt)[pho1_index])&& abs((*phoEta)[pho1_index])<1.;
      bool passCorrPFPhoIso = ((*phoPFPhoIso)[pho1_index]-rho*effArea)<(2.362+0.0047*(*phoEt)[pho1_index])&& abs((*phoEta)[pho1_index])<1.;
   
      if(passPFPhoIso) h_passPFPhoIso_nvtx.Fill(nVtx);
      if(passCorrPFPhoIso)h_passCorrPFPhoIso_nvtx.Fill(nVtx);



      //-------------------------- APPLY PHOTON ID USING BOTH HAND IMPLEMENTATION AND VID -----//
      //apply medium photon ID in the barrel
      double effAreaPFPhoIsoEGM=0.1210;
      double effAreaPFChIsoEGM=0.0360;
      double effAreaPFNeuIsoEGM=0.0597;
      bool passEtaRange_1=abs((*phoEta)[pho1_index])<1.;
      bool passPFPhoIso_1=((*phoPFPhoIso)[pho1_index]-rho*effAreaPFPhoIsoEGM)<(2.571+0.0047*(*phoEt)[pho1_index]);
      bool passPFChIso_1=((*phoPFChIso)[pho1_index]-rho*effAreaPFChIsoEGM)<0.441;
      bool passPFNeuIso_1=((*phoPFNeuIso)[pho1_index]-rho*effAreaPFNeuIsoEGM)<(2.725+0.0148*(*phoEt)[pho1_index]+0.000017*(*phoEt)[pho1_index]*(*phoEt)[pho1_index]);
      bool passHoE_1=(*phoHoverE)[pho1_index]<0.0396;
      bool passSieie_1= (*phoSigmaIEtaIEtaFull5x5)[pho1_index]<0.01022;
      
      bool passPhID_1 = passEtaRange_1 && passPFPhoIso_1 && passPFChIso_1 && passPFNeuIso_1 && passHoE_1 &&passSieie_1;
      
      bool passPhID_VID_1 = passEtaRange_1 && (((*phoIDbit)[pho1_index] & 2 ) ==2); //medium ID: 2nd bit equal to 1

      if(passPhID_1!=passPhID_VID_1)      std::cout<<passPhID_1<<" "<<passPhID_VID_1<<std::endl;



      //-------------------------- LOOK AT DIPHOTON MASS  WITH DIFFERENT IDs ------------------------------------------------//                                               
      if((((*phoIDbit)[pho1_index] & 1 ) ==1) && (((*phoIDbit)[pho2_index] & 1 ) ==1))      h_diphotonMass_LooseID->Fill(diphotonMass);
      if((((*phoIDbit)[pho1_index] & 2 ) ==2) && (((*phoIDbit)[pho2_index] & 2 ) ==2))      h_diphotonMass_MediumID->Fill(diphotonMass);
      if((((*phoIDbit)[pho1_index] & 4 ) ==4) && (((*phoIDbit)[pho2_index] & 4 ) ==4))      h_diphotonMass_TightID->Fill(diphotonMass);


      //-------------------------- Perform DM-like AnalysisMASTER -----------------------------------------//
      //no cuts
      h_DiPhotonMassCorrWeight->Fill(diphotonMassCorr,weight);
      h_METWeight->Fill(pfMET, weight);

      //applying kinematic cuts
      double pt1OverM = (*phoEt)[pho1_index]/diphotonMassCorr;
      double pt2OverM = (*phoEt)[pho2_index]/diphotonMassCorr;
      if(pt1OverM>0.3 && pt2OverM>0.2) {
	h_DiPhotonMassCorrWeight_KinCuts->Fill(diphotonMassCorr,weight);
	h_METWeight_KinCuts->Fill(pfMET, weight);

	//applying medium ID
	if((((*phoIDbit)[pho1_index] & 2 ) ==2) && (((*phoIDbit)[pho2_index] & 2 ) ==2)) {
	  h_DiPhotonMassCorrWeight_KinCuts_MedID->Fill(diphotonMassCorr,weight);
	  h_METWeight_KinCuts_MedID->Fill(pfMET, weight);

	  //applying MET cut at 130 GeV
	  if(pfMET>130){
	    
	    h_DiPhotonMassCorrWeight_KinCuts_MedID_METCut->Fill(diphotonMassCorr,weight);
	    h_METWeight_KinCuts_MedID_METCut->Fill(pfMET, weight);

	    
	  }
	}
      }
      //-------------------------- FILL OTHER HISTOS ------------------------------------------------//
      h_phoE_1->Fill((*phoE)[pho1_index]);
      h_phoEt_1->Fill((*phoEt)[pho1_index]);
      h_phoEta_1->Fill((*phoEta)[pho1_index]);

      if(pho1_isMatch && abs((*phoEta)[pho1_index])<1.4442) h_phoSieieEB_Matched_1->Fill((*phoSigmaIEtaIEtaFull5x5)[pho1_index]);
      else if(!pho1_isMatch && abs((*phoEta)[pho1_index])<1.4442) h_phoSieieEB_NotMatched_1->Fill((*phoSigmaIEtaIEtaFull5x5)[pho1_index]);
      
   }


   //------------------------- COMPUTE ISOLATION EFFICIENCY VS PU W/ AND W/O CORRECTIONS----//
   TEfficiency* PFPhoIso_eff = new TEfficiency(h_passPFPhoIso_nvtx,h_nvtx);
   PFPhoIso_eff->SetTitle("PFPhoIso_eff");
   PFPhoIso_eff->SetName("PFPhoIso_eff");
   TEfficiency* CorrPFPhoIso_eff = new TEfficiency(h_passCorrPFPhoIso_nvtx,h_nvtx);
   CorrPFPhoIso_eff->SetTitle("CorrPFPhoIso_eff");
   CorrPFPhoIso_eff->SetName("CorrPFPhoIso_eff");
   
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
   PFPhoIso_eff->Write();
   CorrPFPhoIso_eff->Write();
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
