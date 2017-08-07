
// last commit by $Id: analysis.cc,v 1.1 2013/01/31 15:32:00 soffi Exp $
//
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>

#include <TTree.h>
#include <TChain.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TAxis.h>
#include <TGraph.h>

#include "Analysis.C"

using namespace std;



int runAnalysis(int whichSample=1, int maxEvents=1000) {

  double weight;
  //weight = L*xsec*BR/Nevents
  //XSEC 8 BR [fb]: ggF = 110 , VH = 5, VBF = 8, ttH =1, Signal (Baryonic DM , mDM=1 GeV, mMed = 1 TeV) = 0.4 

  if(whichSample==0)weight=1; //signal
  if(whichSample==1)weight=1; //ggF
  if(whichSample==2)weight=1; //VH
  if(whichSample==3)weight=1; //VBF
  if(whichSample==4)weight=1; //ttH
  if(whichSample==5)weight=1.; //QCD Bkg - weight not relevant

      //================ Creating chain 
  


  TFile* fin;
  if(whichSample==0) fin= new TFile("root://eoscms//eos/cms/store/group/phys_egamma/soffi/CMSPOS2017/ggntuples/DM_MonoHgg_Mmed1000GeV_Mdm_1GeV_BaryonicModel.root");
  else if(whichSample==1) fin= new TFile("root://eoscms//eos/cms/store/group/phys_egamma/soffi/CMSPOS2017/ggntuples/GluGluHToGG_M-125_13TeV_powheg_pythia8.root");
  else if(whichSample==2) fin= new TFile("root://eoscms//eos/cms/store/group/phys_egamma/soffi/CMSPOS2017/ggntuples/VHToGG_M-125_13TeV_powheg_pythia8_30K.root");
  else if (whichSample==3) fin= new TFile("root://eoscms//eos/cms/store/group/phys_egamma/soffi/CMSPOS2017/ggntuples/VBFHToGG_M-125_13TeV_powheg_pythia8.root");
  else if (whichSample==4) fin= new TFile("root://eoscms//eos/cms/store/group/phys_egamma/soffi/CMSPOS2017/ggntuples/ttHToGG_M-125_13TeV_powheg_pythia8.root");
  else if (whichSample==5)  fin= new TFile("root://eoscms//eos/cms/store/group/phys_egamma/soffi/CMSPOS2017/ggntuples/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8.root");
  TChain* chain =new TChain("EventTree");

  if(whichSample==0) chain->Add("root://eoscms//eos/cms/store/group/phys_egamma/soffi/CMSPOS2017/ggntuples/DM_MonoHgg_Mmed1000GeV_Mdm_1GeV_BaryonicModel.root/ggNtuplizer/EventTree");
  else if(whichSample==1)  chain->Add("root://eoscms//eos/cms/store/group/phys_egamma/soffi/CMSPOS2017/ggntuples/GluGluHToGG_M-125_13TeV_powheg_pythia8.root/ggNtuplizer/EventTree");
  else if(whichSample==2)  chain->Add("root://eoscms//eos/cms/store/group/phys_egamma/soffi/CMSPOS2017/ggntuples/VHToGG_M-125_13TeV_powheg_pythia8_30K.root/ggNtuplizer/EventTree");
  else if(whichSample==3)  chain->Add("root://eoscms//eos/cms/store/group/phys_egamma/soffi/CMSPOS2017/ggntuples/VBFHToGG_M-125_13TeV_powheg_pythia8.root/ggNtuplizer/EventTree");
  else if(whichSample==4)  chain->Add("root://eoscms//eos/cms/store/group/phys_egamma/soffi/CMSPOS2017/ggntuples/ttHToGG_M-125_13TeV_powheg_pythia8.root/ggNtuplizer/EventTree");
  else if(whichSample==5)  chain->Add("root://eoscms//eos/cms/store/group/phys_egamma/soffi/CMSPOS2017/ggntuples/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8.root/ggNtuplizer/EventTree");

  if(whichSample==0)std::cout<<" Running on DM Signal" <<std::endl;
  else   if(whichSample==1)std::cout<<" Running on  ggF Hgg" <<std::endl;
  else   if(whichSample==2)std::cout<<" Running on  VH Hgg" <<std::endl;
  else   if(whichSample==3)std::cout<<" Running on  VBF Hgg" <<std::endl;
  else   if(whichSample==4)std::cout<<" Running on  ttH Hgg" <<std::endl;
  else   if(whichSample==5)std::cout<<" Running on Background QCD" <<std::endl;

  //================ Run analysis
  Analysis tree( chain );
  tree.Loop(whichSample,maxEvents, weight);
  
  delete fin;
  delete chain;
  return 0;
}





void makeROCcurve(){
  TFile* f_bkg= TFile::Open("fout_Background_QCD.root");
  TFile* f_sig= TFile::Open("fout_Signal_Hgg.root");

  TH1F* bkg = (TH1F*) f_bkg->Get("h_phoSieieEB_NotMatched_1");
  TH1F* sig = (TH1F*) f_sig->Get("h_phoSieieEB_Matched_1");

  bkg->Sumw2();
  sig->Sumw2();

  TGraph *g_e = new TGraph();

  g_e->SetTitle( "ROCs" );

  g_e->GetXaxis()->SetTitle( "signal efficiency" );
  g_e->GetXaxis()->SetRangeUser( 0, 1 );
  g_e->GetYaxis()->SetRangeUser( 0, 1 );
  g_e->GetYaxis()->SetTitle( "background rejection" );

  TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();
  c->SetLogy();
  c->SetLogx();
  g_e->Draw( "apl same" );
  c->SaveAs("~/www/CMSPOS2017/ROCcurve.png");


}
