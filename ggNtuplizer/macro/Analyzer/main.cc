
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

#include "Analysis.C"

using namespace std;



int runAnalysis() {


      //================ Creating chain 
  

  TFile* fin= new TFile("root://eoscms//eos/cms/store/group/phys_egamma/soffi/CMSPOS2017/ggntuples/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
  TChain* chain =new TChain("EventTree");

  chain->Add("root://eoscms//eos/cms/store/group/phys_egamma/soffi/CMSPOS2017/ggntuples/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root/ggNtuplizer/EventTree");  


  // std::cout<<chain->GetEntries()<<std::endl;      
  //================ Run analysis
  Analysis tree( chain );
  tree.Loop(false);
  
  delete fin;
  delete chain;
  return 0;
}
