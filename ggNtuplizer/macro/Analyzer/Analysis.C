#define Analysis_cxx
#include "Analysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void Analysis::Loop(bool isSig)
{
  //create a new output file
  TFile* fout;
  if(isSig) fout=new TFile("fout_HggSignal.root", "RECREATE");
  else fout=new TFile("fout_Background.root", "RECREATE");

  //book histos
  TH1F* h_phoE_1 = new TH1F("h_phoE_1", "", 120, 0, 300);
  h_phoE_1->Sumw2();


   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   //   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   for (Long64_t jentry=0; jentry<1000;jentry++) { //loop on only 1k events
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
    
      //look only events with at least otwo photon
      double size=phoE->size();
      if(size<2)continue;
      h_phoE_1->Fill((*phoE)[0]);

   }


   //save output
   fout->cd();
   h_phoE_1->Write();
   fout->Write();
   fout->Close();
}
