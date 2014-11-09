//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct  8 16:46:53 2014 by ROOT version 5.34/18
// from TTree Tree_2Dtruebkgbkg_template/Tree_2Dtruebkgbkg_template
// found on file: allmc.root
//////////////////////////////////////////////////////////

#ifndef light_h
#define light_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include "TH1F.h"
#include "TRandom3.h"
// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class light {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          event_fileuuid;
   Int_t           event_run;
   Int_t           event_lumi;
   UInt_t          event_number;
   Int_t           dataset_id;
   Float_t         event_luminormfactor;
   Float_t         event_Kfactor;
   Float_t         event_weight;
   Float_t         event_rho;
   Float_t         event_sigma;
   Int_t           event_nPU;
   Int_t           event_nPUtrue;
   Int_t           event_nRecVtx;
   Int_t           event_pass12whoissiglike;
   Float_t         dipho_mgg_photon;
   Float_t         pholead_eta;
   Float_t         photrail_eta;
   Float_t         pholead_phi;
   Float_t         photrail_phi;
   Float_t         pholead_pt;
   Float_t         photrail_pt;
   Float_t         pholead_energy;
   Float_t         photrail_energy;
   Float_t         pholead_SCeta;
   Float_t         photrail_SCeta;
   Float_t         pholead_SCphi;
   Float_t         photrail_SCphi;
   Int_t           pholead_PhoHasPixSeed;
   Int_t           photrail_PhoHasPixSeed;
   Float_t         pholead_r9;
   Float_t         photrail_r9;
   Float_t         pholead_sieie;
   Float_t         photrail_sieie;
   Float_t         pholead_hoe;
   Float_t         photrail_hoe;
   Float_t         pholead_PhoSCRemovalPFIsoCharged;
   Float_t         photrail_PhoSCRemovalPFIsoCharged;
   Float_t         pholead_PhoSCRemovalPFIsoNeutral;
   Float_t         photrail_PhoSCRemovalPFIsoNeutral;
   Float_t         pholead_PhoSCRemovalPFIsoPhoton;
   Float_t         photrail_PhoSCRemovalPFIsoPhoton;
   Float_t         pholead_PhoSCRemovalPFIsoCombined;
   Float_t         photrail_PhoSCRemovalPFIsoCombined;
   Float_t         pholead_PhoIso03Ecal;
   Float_t         pholead_PhoIso03Hcal;
   Float_t         pholead_PhoIso03TrkHollow;
   Float_t         photrail_PhoIso03Ecal;
   Float_t         photrail_PhoIso03Hcal;
   Float_t         photrail_PhoIso03TrkHollow;
   Int_t           pholead_PhoPassConversionVeto;
   Int_t           photrail_PhoPassConversionVeto;
   Float_t         pholead_GenPhotonIsoDR04;
   Float_t         photrail_GenPhotonIsoDR04;
   Int_t           pholead_PhoMCmatchexitcode;
   Int_t           photrail_PhoMCmatchexitcode;
   Float_t         mulead_pt;
   Float_t         mutrail_pt;
   Float_t         mulead_eta;
   Float_t         mutrail_eta;
   Float_t         mulead_phi;
   Float_t         mutrail_phi;
   Float_t         mulead_energy;
   Float_t         mutrail_energy;
   Float_t         mass_mumu;
   Float_t         mass_mumug;
   Float_t         dR1;
   Float_t         dR2;
   Float_t         pholead_m_jet_ptcorr;
   Float_t         pholead_m_jet_dR;
   Float_t         pholead_phopt_footprint_total;
   Float_t         pholead_phopt_footprint_m_frac;
   Float_t         pholead_jetpt_pf;
   Float_t         pholead_jetpt_m_frac;
   Float_t         pholead_jetpt_m_frac_PhoComp;
   Float_t         photrail_m_jet_ptcorr;
   Float_t         photrail_m_jet_dR;
   Float_t         photrail_phopt_footprint_total;
   Float_t         photrail_phopt_footprint_m_frac;
   Float_t         photrail_jetpt_pf;
   Float_t         photrail_jetpt_m_frac;
   Float_t         photrail_jetpt_m_frac_PhoComp;
   Float_t         pholead_pt_closestjet;
   Float_t         pholead_dR_closestjet;
   Float_t         photrail_pt_closestjet;
   Float_t         photrail_dR_closestjet;
   Float_t         phoiso_template_1event_sigsig_1[5];
   Float_t         phoiso_template_1event_sigsig_2[5];
   Float_t         phoiso_template_1event_sigbkg_1[5];
   Float_t         phoiso_template_1event_sigbkg_2[5];
   Float_t         phoiso_template_1event_bkgsig_1[5];
   Float_t         phoiso_template_1event_bkgsig_2[5];
   Float_t         phoiso_template_1event_bkgbkg_1[5];
   Float_t         phoiso_template_1event_bkgbkg_2[5];
   Float_t         phoiso_template_2events_sigsig_1[5];
   Float_t         phoiso_template_2events_sigsig_2[5];
   Float_t         phoiso_template_2events_sigbkg_1[5];
   Float_t         phoiso_template_2events_sigbkg_2[5];
   Float_t         phoiso_template_2events_bkgsig_1[5];
   Float_t         phoiso_template_2events_bkgsig_2[5];
   Float_t         phoiso_template_2events_bkgbkg_1[5];
   Float_t         phoiso_template_2events_bkgbkg_2[5];
   Float_t         rewinfo_template_1event_sigsig_1[30];
   Float_t         rewinfo_template_1event_sigsig_2[30];
   Float_t         rewinfo_template_1event_sigbkg_1[30];
   Float_t         rewinfo_template_1event_sigbkg_2[30];
   Float_t         rewinfo_template_1event_bkgsig_1[30];
   Float_t         rewinfo_template_1event_bkgsig_2[30];
   Float_t         rewinfo_template_1event_bkgbkg_1[30];
   Float_t         rewinfo_template_1event_bkgbkg_2[30];
   Float_t         rewinfo_template_2events_sigsig_1[30];
   Float_t         rewinfo_template_2events_sigsig_2[30];
   Float_t         rewinfo_template_2events_sigbkg_1[30];
   Float_t         rewinfo_template_2events_sigbkg_2[30];
   Float_t         rewinfo_template_2events_bkgsig_1[30];
   Float_t         rewinfo_template_2events_bkgsig_2[30];
   Float_t         rewinfo_template_2events_bkgbkg_1[30];
   Float_t         rewinfo_template_2events_bkgbkg_2[30];
   Int_t           n_jets;
   Float_t         jet_pt[8];   //[n_jets]
   Float_t         jet_eta[8];   //[n_jets]
   Float_t         jet_phi[8];   //[n_jets]
   Float_t         jet_energy[8];   //[n_jets]

   // List of branches
   TBranch        *b_event_fileuuid;   //!
   TBranch        *b_event_run;   //!
   TBranch        *b_event_lumi;   //!
   TBranch        *b_event_number;   //!
   TBranch        *b_dataset_id;   //!
   TBranch        *b_event_luminormfactor;   //!
   TBranch        *b_event_Kfactor;   //!
   TBranch        *b_event_weight;   //!
   TBranch        *b_event_rho;   //!
   TBranch        *b_event_sigma;   //!
   TBranch        *b_event_nPU;   //!
   TBranch        *b_event_nPUtrue;   //!
   TBranch        *b_event_nRecVtx;   //!
   TBranch        *b_event_pass12whoissiglike;   //!
   TBranch        *b_dipho_mgg_photon;   //!
   TBranch        *b_pholead_eta;   //!
   TBranch        *b_photrail_eta;   //!
   TBranch        *b_pholead_phi;   //!
   TBranch        *b_photrail_phi;   //!
   TBranch        *b_pholead_pt;   //!
   TBranch        *b_photrail_pt;   //!
   TBranch        *b_pholead_energy;   //!
   TBranch        *b_photrail_energy;   //!
   TBranch        *b_pholead_SCeta;   //!
   TBranch        *b_photrail_SCeta;   //!
   TBranch        *b_pholead_SCphi;   //!
   TBranch        *b_photrail_SCphi;   //!
   TBranch        *b_pholead_PhoHasPixSeed;   //!
   TBranch        *b_photrail_PhoHasPixSeed;   //!
   TBranch        *b_pholead_r9;   //!
   TBranch        *b_photrail_r9;   //!
   TBranch        *b_pholead_sieie;   //!
   TBranch        *b_photrail_sieie;   //!
   TBranch        *b_pholead_hoe;   //!
   TBranch        *b_photrail_hoe;   //!
   TBranch        *b_pholead_PhoSCRemovalPFIsoCharged;   //!
   TBranch        *b_photrail_PhoSCRemovalPFIsoCharged;   //!
   TBranch        *b_pholead_PhoSCRemovalPFIsoNeutral;   //!
   TBranch        *b_photrail_PhoSCRemovalPFIsoNeutral;   //!
   TBranch        *b_pholead_PhoSCRemovalPFIsoPhoton;   //!
   TBranch        *b_photrail_PhoSCRemovalPFIsoPhoton;   //!
   TBranch        *b_pholead_PhoSCRemovalPFIsoCombined;   //!
   TBranch        *b_photrail_PhoSCRemovalPFIsoCombined;   //!
   TBranch        *b_pholead_PhoIso03Ecal;   //!
   TBranch        *b_pholead_PhoIso03Hcal;   //!
   TBranch        *b_pholead_PhoIso03TrkHollow;   //!
   TBranch        *b_photrail_PhoIso03Ecal;   //!
   TBranch        *b_photrail_PhoIso03Hcal;   //!
   TBranch        *b_photrail_PhoIso03TrkHollow;   //!
   TBranch        *b_pholead_PhoPassConversionVeto;   //!
   TBranch        *b_photrail_PhoPassConversionVeto;   //!
   TBranch        *b_pholead_GenPhotonIsoDR04;   //!
   TBranch        *b_photrail_GenPhotonIsoDR04;   //!
   TBranch        *b_pholead_PhoMCmatchexitcode;   //!
   TBranch        *b_photrail_PhoMCmatchexitcode;   //!
   TBranch        *b_mulead_pt;   //!
   TBranch        *b_mutrail_pt;   //!
   TBranch        *b_mulead_eta;   //!
   TBranch        *b_mutrail_eta;   //!
   TBranch        *b_mulead_phi;   //!
   TBranch        *b_mutrail_phi;   //!
   TBranch        *b_mulead_energy;   //!
   TBranch        *b_mutrail_energy;   //!
   TBranch        *b_mass_mumu;   //!
   TBranch        *b_mass_mumug;   //!
   TBranch        *b_dR1;   //!
   TBranch        *b_dR2;   //!
   TBranch        *b_pholead_m_jet_ptcorr;   //!
   TBranch        *b_pholead_m_jet_dR;   //!
   TBranch        *b_pholead_phopt_footprint_total;   //!
   TBranch        *b_pholead_phopt_footprint_m_frac;   //!
   TBranch        *b_pholead_jetpt_pf;   //!
   TBranch        *b_pholead_jetpt_m_frac;   //!
   TBranch        *b_pholead_jetpt_m_frac_PhoComp;   //!
   TBranch        *b_photrail_m_jet_ptcorr;   //!
   TBranch        *b_photrail_m_jet_dR;   //!
   TBranch        *b_photrail_phopt_footprint_total;   //!
   TBranch        *b_photrail_phopt_footprint_m_frac;   //!
   TBranch        *b_photrail_jetpt_pf;   //!
   TBranch        *b_photrail_jetpt_m_frac;   //!
   TBranch        *b_photrail_jetpt_m_frac_PhoComp;   //!
   TBranch        *b_pholead_pt_closestjet;   //!
   TBranch        *b_pholead_dR_closestjet;   //!
   TBranch        *b_photrail_pt_closestjet;   //!
   TBranch        *b_photrail_dR_closestjet;   //!
   TBranch        *b_phoiso_template_1event_sigsig_1;   //!
   TBranch        *b_phoiso_template_1event_sigsig_2;   //!
   TBranch        *b_phoiso_template_1event_sigbkg_1;   //!
   TBranch        *b_phoiso_template_1event_sigbkg_2;   //!
   TBranch        *b_phoiso_template_1event_bkgsig_1;   //!
   TBranch        *b_phoiso_template_1event_bkgsig_2;   //!
   TBranch        *b_phoiso_template_1event_bkgbkg_1;   //!
   TBranch        *b_phoiso_template_1event_bkgbkg_2;   //!
   TBranch        *b_phoiso_template_2events_sigsig_1;   //!
   TBranch        *b_phoiso_template_2events_sigsig_2;   //!
   TBranch        *b_phoiso_template_2events_sigbkg_1;   //!
   TBranch        *b_phoiso_template_2events_sigbkg_2;   //!
   TBranch        *b_phoiso_template_2events_bkgsig_1;   //!
   TBranch        *b_phoiso_template_2events_bkgsig_2;   //!
   TBranch        *b_phoiso_template_2events_bkgbkg_1;   //!
   TBranch        *b_phoiso_template_2events_bkgbkg_2;   //!
   TBranch        *b_rewinfo_template_1event_sigsig_1;   //!
   TBranch        *b_rewinfo_template_1event_sigsig_2;   //!
   TBranch        *b_rewinfo_template_1event_sigbkg_1;   //!
   TBranch        *b_rewinfo_template_1event_sigbkg_2;   //!
   TBranch        *b_rewinfo_template_1event_bkgsig_1;   //!
   TBranch        *b_rewinfo_template_1event_bkgsig_2;   //!
   TBranch        *b_rewinfo_template_1event_bkgbkg_1;   //!
   TBranch        *b_rewinfo_template_1event_bkgbkg_2;   //!
   TBranch        *b_rewinfo_template_2events_sigsig_1;   //!
   TBranch        *b_rewinfo_template_2events_sigsig_2;   //!
   TBranch        *b_rewinfo_template_2events_sigbkg_1;   //!
   TBranch        *b_rewinfo_template_2events_sigbkg_2;   //!
   TBranch        *b_rewinfo_template_2events_bkgsig_1;   //!
   TBranch        *b_rewinfo_template_2events_bkgsig_2;   //!
   TBranch        *b_rewinfo_template_2events_bkgbkg_1;   //!
   TBranch        *b_rewinfo_template_2events_bkgbkg_2;   //!
   TBranch        *b_n_jets;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_energy;   //!

   //own varialbles
  TH1F *EB_temp;
  TH1F *EB_si;
 TH1F *m1EE_temp ;
 TH1F *EE_temp;
 TH1F *m1EE_si;
  TH1F *EE_si;
 TH1F *m1EE_pt;
 TH1F *m1EE_diphopt;
  TH1F *m1EE_diphomass;
 TH1F *EE_diphomass;
  TH1F *EB_diphomass;
  TH1F *h_weight;
   Bool_t do2dstd; 
   Bool_t do2dff;
   Bool_t do2dpp;
   Bool_t do2d1side;
   Bool_t do1p1f;
   Bool_t do2dside;
   Bool_t dorcone;
   Bool_t isdata;
   Float_t roovar1;
   Float_t roovar2;
   Float_t roopt1;
   Float_t roopt2;
   Float_t roosieie1;
   Float_t roosieie2;
   Float_t rooeta1;
   Float_t rooeta2;
   Float_t roorho;
   Float_t roosigma;
   Float_t roonvtx;
   Float_t rooweight;
  Int_t nbins; 
 TRandom3 *randomgen;


 
   light(TString, TString, Bool_t);
   virtual ~light();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   
  
};

#endif


#ifdef light_cxx
light::light(TString filename, TString treename, Bool_t realdata) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  // if (tree == 0) {
 //     TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("allmc.root");
   //   if (!f || !f->IsOpen()) {
       TFile*  f = new TFile(filename.Data());
    //  }
TTree *tree;
      f->GetObject(treename.Data(),tree);
      std::cout << treename << std::endl;
  // }
   Init(tree);
  TH1F::SetDefaultSumw2(kTRUE);


  Float_t xstart=400.;
        nbins=46;
        Float_t xbins[46];
        xbins[0]=0.;
        for(int i=0 ; i <nbins ; i++){
        if (xbins[i] < xstart)  xbins[i+1]= xbins[i]+10;
        else if ( xstart <= xbins[i] < 1000)  xbins[i+1]=xbins[i]+100;

}

 EB_temp= new TH1F("EB_temp","EB_temp",90, 0.,9.);
 EB_si= new TH1F("EB_si","EB_si", 400, 0, 0.04);
  m1EE_temp= new TH1F("m1EE_temp","m1EE_temp", 90, 0., 9.);
  EE_temp= new TH1F("EE_temp","EE_temp", 90, 0., 9.);
 m1EE_si= new TH1F("m1EE_si","m1EE_si", 400, 0, 0.04);
  EE_si= new TH1F("EE_si","EE_si", 400, 0, 0.04);
 m1EE_pt= new TH1F("m1EE_pt","m1EE_pt", 200, 0., 200.);
m1EE_diphopt= new TH1F("m1EE_diphopt","m1EE_diphopt", 600, 0., 600.);
  m1EE_diphomass= new TH1F("m1EE_diphomass","m1EE_diphomass", nbins, xbins);
  EE_diphomass= new TH1F("EE_diphomass","EE_diphomass", nbins, xbins);
  EB_diphomass= new TH1F("EB_diphomass","EB_diphomass", nbins, xbins);
  h_weight =new TH1F("h_weight","h_weight", 500, 0.,50.); 
  do2dstd=false;
   do2dff=false;
   do2dpp=false;
   do2d1side=false;
   do2dside=false;
   do1p1f=false;
   dorcone=false;
   isdata= realdata;
   if(treename== "Tree_2Dtruebkgbkg_template"){do2dff=true;}
   else if(treename== "Tree_2Dstandard_selection"){do2dstd=true;}
   else if(treename== "Tree_2Dsideband_template"){do2dside=true;}
   else if(treename== "Tree_2Dtruesigsig_template"){do2dpp=true; }
   else if(treename== "Tree_2Drandomcone_template"){dorcone=true;}
   else if(treename== "Tree_2Drandomconesideband_template"){do2d1side=true;}
   //141024 keep in mind no gen prompt anymore!
   else if(treename== "Tree_2Dgenpromptplussideband_template"){do2d1side=true;}
   else if(treename== "Tree_2Dtruesigbkg_template"){do1p1f=true;isdata=false;}
   //TFile *f ;
//TTree *tree ;
//tree->Branch("roodset_signal",&roovar1,"roovar1/F:roovar1:roovar2:rooeta1:rooeta2:roopt1:roopt2:roosieie1:roosieie2:rooweight:roorho:roosigma:rooweight");
//TBranch* b1;
 randomgen = new TRandom3(0);

  roovar1 = -999;
  roovar2 = -999;
  rooeta1 = -999;
  rooeta2 = -999;
  roopt1 = -999;
  roopt2 = -999;
  roosieie1 = -999;
  roosieie2 = -999;
  roorho = -999;
  roosigma = -999;
  roonvtx = -999;
  rooweight = -999;


}

light::~light()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   delete randomgen;
}

Int_t light::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t light::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void light::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event_fileuuid", &event_fileuuid, &b_event_fileuuid);
   fChain->SetBranchAddress("event_run", &event_run, &b_event_run);
   fChain->SetBranchAddress("event_lumi", &event_lumi, &b_event_lumi);
   fChain->SetBranchAddress("event_number", &event_number, &b_event_number);
   fChain->SetBranchAddress("dataset_id", &dataset_id, &b_dataset_id);
   fChain->SetBranchAddress("event_luminormfactor", &event_luminormfactor, &b_event_luminormfactor);
   fChain->SetBranchAddress("event_Kfactor", &event_Kfactor, &b_event_Kfactor);
   fChain->SetBranchAddress("event_weight", &event_weight, &b_event_weight);
   fChain->SetBranchAddress("event_rho", &event_rho, &b_event_rho);
   fChain->SetBranchAddress("event_sigma", &event_sigma, &b_event_sigma);
   fChain->SetBranchAddress("event_nPU", &event_nPU, &b_event_nPU);
   fChain->SetBranchAddress("event_nPUtrue", &event_nPUtrue, &b_event_nPUtrue);
   fChain->SetBranchAddress("event_nRecVtx", &event_nRecVtx, &b_event_nRecVtx);
   fChain->SetBranchAddress("event_pass12whoissiglike", &event_pass12whoissiglike, &b_event_pass12whoissiglike);
   fChain->SetBranchAddress("dipho_mgg_photon", &dipho_mgg_photon, &b_dipho_mgg_photon);
   fChain->SetBranchAddress("pholead_eta", &pholead_eta, &b_pholead_eta);
   fChain->SetBranchAddress("photrail_eta", &photrail_eta, &b_photrail_eta);
   fChain->SetBranchAddress("pholead_phi", &pholead_phi, &b_pholead_phi);
   fChain->SetBranchAddress("photrail_phi", &photrail_phi, &b_photrail_phi);
   fChain->SetBranchAddress("pholead_pt", &pholead_pt, &b_pholead_pt);
   fChain->SetBranchAddress("photrail_pt", &photrail_pt, &b_photrail_pt);
   fChain->SetBranchAddress("pholead_energy", &pholead_energy, &b_pholead_energy);
   fChain->SetBranchAddress("photrail_energy", &photrail_energy, &b_photrail_energy);
   fChain->SetBranchAddress("pholead_SCeta", &pholead_SCeta, &b_pholead_SCeta);
   fChain->SetBranchAddress("photrail_SCeta", &photrail_SCeta, &b_photrail_SCeta);
   fChain->SetBranchAddress("pholead_SCphi", &pholead_SCphi, &b_pholead_SCphi);
   fChain->SetBranchAddress("photrail_SCphi", &photrail_SCphi, &b_photrail_SCphi);
   fChain->SetBranchAddress("pholead_PhoHasPixSeed", &pholead_PhoHasPixSeed, &b_pholead_PhoHasPixSeed);
   fChain->SetBranchAddress("photrail_PhoHasPixSeed", &photrail_PhoHasPixSeed, &b_photrail_PhoHasPixSeed);
   fChain->SetBranchAddress("pholead_r9", &pholead_r9, &b_pholead_r9);
   fChain->SetBranchAddress("photrail_r9", &photrail_r9, &b_photrail_r9);
   fChain->SetBranchAddress("pholead_sieie", &pholead_sieie, &b_pholead_sieie);
   fChain->SetBranchAddress("photrail_sieie", &photrail_sieie, &b_photrail_sieie);
   fChain->SetBranchAddress("pholead_hoe", &pholead_hoe, &b_pholead_hoe);
   fChain->SetBranchAddress("photrail_hoe", &photrail_hoe, &b_photrail_hoe);
   fChain->SetBranchAddress("pholead_PhoSCRemovalPFIsoCharged", &pholead_PhoSCRemovalPFIsoCharged, &b_pholead_PhoSCRemovalPFIsoCharged);
   fChain->SetBranchAddress("photrail_PhoSCRemovalPFIsoCharged", &photrail_PhoSCRemovalPFIsoCharged, &b_photrail_PhoSCRemovalPFIsoCharged);
   fChain->SetBranchAddress("pholead_PhoSCRemovalPFIsoNeutral", &pholead_PhoSCRemovalPFIsoNeutral, &b_pholead_PhoSCRemovalPFIsoNeutral);
   fChain->SetBranchAddress("photrail_PhoSCRemovalPFIsoNeutral", &photrail_PhoSCRemovalPFIsoNeutral, &b_photrail_PhoSCRemovalPFIsoNeutral);
   fChain->SetBranchAddress("pholead_PhoSCRemovalPFIsoPhoton", &pholead_PhoSCRemovalPFIsoPhoton, &b_pholead_PhoSCRemovalPFIsoPhoton);
   fChain->SetBranchAddress("photrail_PhoSCRemovalPFIsoPhoton", &photrail_PhoSCRemovalPFIsoPhoton, &b_photrail_PhoSCRemovalPFIsoPhoton);
   fChain->SetBranchAddress("pholead_PhoSCRemovalPFIsoCombined", &pholead_PhoSCRemovalPFIsoCombined, &b_pholead_PhoSCRemovalPFIsoCombined);
   fChain->SetBranchAddress("photrail_PhoSCRemovalPFIsoCombined", &photrail_PhoSCRemovalPFIsoCombined, &b_photrail_PhoSCRemovalPFIsoCombined);
   fChain->SetBranchAddress("pholead_PhoIso03Ecal", &pholead_PhoIso03Ecal, &b_pholead_PhoIso03Ecal);
   fChain->SetBranchAddress("pholead_PhoIso03Hcal", &pholead_PhoIso03Hcal, &b_pholead_PhoIso03Hcal);
   fChain->SetBranchAddress("pholead_PhoIso03TrkHollow", &pholead_PhoIso03TrkHollow, &b_pholead_PhoIso03TrkHollow);
   fChain->SetBranchAddress("photrail_PhoIso03Ecal", &photrail_PhoIso03Ecal, &b_photrail_PhoIso03Ecal);
   fChain->SetBranchAddress("photrail_PhoIso03Hcal", &photrail_PhoIso03Hcal, &b_photrail_PhoIso03Hcal);
   fChain->SetBranchAddress("photrail_PhoIso03TrkHollow", &photrail_PhoIso03TrkHollow, &b_photrail_PhoIso03TrkHollow);
   fChain->SetBranchAddress("pholead_PhoPassConversionVeto", &pholead_PhoPassConversionVeto, &b_pholead_PhoPassConversionVeto);
   fChain->SetBranchAddress("photrail_PhoPassConversionVeto", &photrail_PhoPassConversionVeto, &b_photrail_PhoPassConversionVeto);
   fChain->SetBranchAddress("pholead_GenPhotonIsoDR04", &pholead_GenPhotonIsoDR04, &b_pholead_GenPhotonIsoDR04);
   fChain->SetBranchAddress("photrail_GenPhotonIsoDR04", &photrail_GenPhotonIsoDR04, &b_photrail_GenPhotonIsoDR04);
   fChain->SetBranchAddress("pholead_PhoMCmatchexitcode", &pholead_PhoMCmatchexitcode, &b_pholead_PhoMCmatchexitcode);
   fChain->SetBranchAddress("photrail_PhoMCmatchexitcode", &photrail_PhoMCmatchexitcode, &b_photrail_PhoMCmatchexitcode);
/*   fChain->SetBranchAddress("mulead_pt", &mulead_pt, &b_mulead_pt);
   fChain->SetBranchAddress("mutrail_pt", &mutrail_pt, &b_mutrail_pt);
   fChain->SetBranchAddress("mulead_eta", &mulead_eta, &b_mulead_eta);
   fChain->SetBranchAddress("mutrail_eta", &mutrail_eta, &b_mutrail_eta);
   fChain->SetBranchAddress("mulead_phi", &mulead_phi, &b_mulead_phi);
   fChain->SetBranchAddress("mutrail_phi", &mutrail_phi, &b_mutrail_phi);
   fChain->SetBranchAddress("mulead_energy", &mulead_energy, &b_mulead_energy);
   fChain->SetBranchAddress("mutrail_energy", &mutrail_energy, &b_mutrail_energy);
   fChain->SetBranchAddress("mass_mumu", &mass_mumu, &b_mass_mumu);
   fChain->SetBranchAddress("mass_mumug", &mass_mumug, &b_mass_mumug);
   fChain->SetBranchAddress("dR1", &dR1, &b_dR1);
   fChain->SetBranchAddress("dR2", &dR2, &b_dR2);*/
   fChain->SetBranchAddress("pholead_m_jet_ptcorr", &pholead_m_jet_ptcorr, &b_pholead_m_jet_ptcorr);
   fChain->SetBranchAddress("pholead_m_jet_dR", &pholead_m_jet_dR, &b_pholead_m_jet_dR);
   fChain->SetBranchAddress("pholead_phopt_footprint_total", &pholead_phopt_footprint_total, &b_pholead_phopt_footprint_total);
   fChain->SetBranchAddress("pholead_phopt_footprint_m_frac", &pholead_phopt_footprint_m_frac, &b_pholead_phopt_footprint_m_frac);
   fChain->SetBranchAddress("pholead_jetpt_pf", &pholead_jetpt_pf, &b_pholead_jetpt_pf);
   fChain->SetBranchAddress("pholead_jetpt_m_frac", &pholead_jetpt_m_frac, &b_pholead_jetpt_m_frac);
   fChain->SetBranchAddress("pholead_jetpt_m_frac_PhoComp", &pholead_jetpt_m_frac_PhoComp, &b_pholead_jetpt_m_frac_PhoComp);
   fChain->SetBranchAddress("photrail_m_jet_ptcorr", &photrail_m_jet_ptcorr, &b_photrail_m_jet_ptcorr);
   fChain->SetBranchAddress("photrail_m_jet_dR", &photrail_m_jet_dR, &b_photrail_m_jet_dR);
   fChain->SetBranchAddress("photrail_phopt_footprint_total", &photrail_phopt_footprint_total, &b_photrail_phopt_footprint_total);
   fChain->SetBranchAddress("photrail_phopt_footprint_m_frac", &photrail_phopt_footprint_m_frac, &b_photrail_phopt_footprint_m_frac);
   fChain->SetBranchAddress("photrail_jetpt_pf", &photrail_jetpt_pf, &b_photrail_jetpt_pf);
   fChain->SetBranchAddress("photrail_jetpt_m_frac", &photrail_jetpt_m_frac, &b_photrail_jetpt_m_frac);
   fChain->SetBranchAddress("photrail_jetpt_m_frac_PhoComp", &photrail_jetpt_m_frac_PhoComp, &b_photrail_jetpt_m_frac_PhoComp);
   fChain->SetBranchAddress("pholead_pt_closestjet", &pholead_pt_closestjet, &b_pholead_pt_closestjet);
   fChain->SetBranchAddress("pholead_dR_closestjet", &pholead_dR_closestjet, &b_pholead_dR_closestjet);
   fChain->SetBranchAddress("photrail_pt_closestjet", &photrail_pt_closestjet, &b_photrail_pt_closestjet);
   fChain->SetBranchAddress("photrail_dR_closestjet", &photrail_dR_closestjet, &b_photrail_dR_closestjet);
   fChain->SetBranchAddress("phoiso_template_1event_sigsig_1", phoiso_template_1event_sigsig_1, &b_phoiso_template_1event_sigsig_1);
   fChain->SetBranchAddress("phoiso_template_1event_sigsig_2", phoiso_template_1event_sigsig_2, &b_phoiso_template_1event_sigsig_2);
   fChain->SetBranchAddress("phoiso_template_1event_sigbkg_1", phoiso_template_1event_sigbkg_1, &b_phoiso_template_1event_sigbkg_1);
   fChain->SetBranchAddress("phoiso_template_1event_sigbkg_2", phoiso_template_1event_sigbkg_2, &b_phoiso_template_1event_sigbkg_2);
   fChain->SetBranchAddress("phoiso_template_1event_bkgsig_1", phoiso_template_1event_bkgsig_1, &b_phoiso_template_1event_bkgsig_1);
   fChain->SetBranchAddress("phoiso_template_1event_bkgsig_2", phoiso_template_1event_bkgsig_2, &b_phoiso_template_1event_bkgsig_2);
   fChain->SetBranchAddress("phoiso_template_1event_bkgbkg_1", phoiso_template_1event_bkgbkg_1, &b_phoiso_template_1event_bkgbkg_1);
   fChain->SetBranchAddress("phoiso_template_1event_bkgbkg_2", phoiso_template_1event_bkgbkg_2, &b_phoiso_template_1event_bkgbkg_2);
   fChain->SetBranchAddress("phoiso_template_2events_sigsig_1", phoiso_template_2events_sigsig_1, &b_phoiso_template_2events_sigsig_1);
   fChain->SetBranchAddress("phoiso_template_2events_sigsig_2", phoiso_template_2events_sigsig_2, &b_phoiso_template_2events_sigsig_2);
   fChain->SetBranchAddress("phoiso_template_2events_sigbkg_1", phoiso_template_2events_sigbkg_1, &b_phoiso_template_2events_sigbkg_1);
   fChain->SetBranchAddress("phoiso_template_2events_sigbkg_2", phoiso_template_2events_sigbkg_2, &b_phoiso_template_2events_sigbkg_2);
   fChain->SetBranchAddress("phoiso_template_2events_bkgsig_1", phoiso_template_2events_bkgsig_1, &b_phoiso_template_2events_bkgsig_1);
   fChain->SetBranchAddress("phoiso_template_2events_bkgsig_2", phoiso_template_2events_bkgsig_2, &b_phoiso_template_2events_bkgsig_2);
   fChain->SetBranchAddress("phoiso_template_2events_bkgbkg_1", phoiso_template_2events_bkgbkg_1, &b_phoiso_template_2events_bkgbkg_1);
   fChain->SetBranchAddress("phoiso_template_2events_bkgbkg_2", phoiso_template_2events_bkgbkg_2, &b_phoiso_template_2events_bkgbkg_2);
   fChain->SetBranchAddress("rewinfo_template_1event_sigsig_1", rewinfo_template_1event_sigsig_1, &b_rewinfo_template_1event_sigsig_1);
   fChain->SetBranchAddress("rewinfo_template_1event_sigsig_2", rewinfo_template_1event_sigsig_2, &b_rewinfo_template_1event_sigsig_2);
   fChain->SetBranchAddress("rewinfo_template_1event_sigbkg_1", rewinfo_template_1event_sigbkg_1, &b_rewinfo_template_1event_sigbkg_1);
   fChain->SetBranchAddress("rewinfo_template_1event_sigbkg_2", rewinfo_template_1event_sigbkg_2, &b_rewinfo_template_1event_sigbkg_2);
   fChain->SetBranchAddress("rewinfo_template_1event_bkgsig_1", rewinfo_template_1event_bkgsig_1, &b_rewinfo_template_1event_bkgsig_1);
   fChain->SetBranchAddress("rewinfo_template_1event_bkgsig_2", rewinfo_template_1event_bkgsig_2, &b_rewinfo_template_1event_bkgsig_2);
   fChain->SetBranchAddress("rewinfo_template_1event_bkgbkg_1", rewinfo_template_1event_bkgbkg_1, &b_rewinfo_template_1event_bkgbkg_1);
   fChain->SetBranchAddress("rewinfo_template_1event_bkgbkg_2", rewinfo_template_1event_bkgbkg_2, &b_rewinfo_template_1event_bkgbkg_2);
   fChain->SetBranchAddress("rewinfo_template_2events_sigsig_1", rewinfo_template_2events_sigsig_1, &b_rewinfo_template_2events_sigsig_1);
   fChain->SetBranchAddress("rewinfo_template_2events_sigsig_2", rewinfo_template_2events_sigsig_2, &b_rewinfo_template_2events_sigsig_2);
   fChain->SetBranchAddress("rewinfo_template_2events_sigbkg_1", rewinfo_template_2events_sigbkg_1, &b_rewinfo_template_2events_sigbkg_1);
   fChain->SetBranchAddress("rewinfo_template_2events_sigbkg_2", rewinfo_template_2events_sigbkg_2, &b_rewinfo_template_2events_sigbkg_2);
   fChain->SetBranchAddress("rewinfo_template_2events_bkgsig_1", rewinfo_template_2events_bkgsig_1, &b_rewinfo_template_2events_bkgsig_1);
   fChain->SetBranchAddress("rewinfo_template_2events_bkgsig_2", rewinfo_template_2events_bkgsig_2, &b_rewinfo_template_2events_bkgsig_2);
   fChain->SetBranchAddress("rewinfo_template_2events_bkgbkg_1", rewinfo_template_2events_bkgbkg_1, &b_rewinfo_template_2events_bkgbkg_1);
   fChain->SetBranchAddress("rewinfo_template_2events_bkgbkg_2", rewinfo_template_2events_bkgbkg_2, &b_rewinfo_template_2events_bkgbkg_2);
   fChain->SetBranchAddress("n_jets", &n_jets, &b_n_jets);
   fChain->SetBranchAddress("jet_pt", jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_energy", jet_energy, &b_jet_energy);
   Notify();
}

Bool_t light::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void light::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t light::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef light_cxx
