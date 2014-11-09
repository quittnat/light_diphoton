
#define light_cxx
#include "light.h"
#include <iostream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TFormula.h"
#include <TH1F.h>
#include "TString.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
using namespace std;
const char*  mc_name ="allmc";
Float_t pixel_lumi= 19.7;
//Float_t pixel_lumi= 1.;

void light::Loop()
{
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
     	if (fChain == 0) return;
TFile *f = new TFile("forroofit/1legsideband_swaped_allmc_EBEB_forroofit.root","RECREATE");
TTree *tree = new TTree("for_roofit","for_roofit");
tree->Branch("roovar1",&roovar1, "roovar1/F");
tree->Branch("roovar2",&roovar2, "roovar2/F");
tree->Branch("rooeta1",&rooeta1, "rooeta1/F");
tree->Branch("rooeta2",&rooeta2, "rooeta2/F");
tree->Branch("roopt1",&roopt1, "roopt1/F");
tree->Branch("roopt2",&roopt2, "roopta2/F");
tree->Branch("roosieie1",&roosieie1, "roosieie1/F");
tree->Branch("roosieie2",&roosieie2, "roosieie2/F");
tree->Branch("roorho",&roorho, "roorho/F");
tree->Branch("roosigma",&roosigma, "roosigma/F");
tree->Branch("roonvtx",&roonvtx, "roonvtx/F");
tree->Branch("rooweight",&rooweight, "rooweight/F");
   Long64_t nentries = fChain->GetEntriesFast();


 Float_t weight=0.; 
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (pholead_pt<40. || photrail_pt<25.) continue;     
      if (rooeta1>1.4442 && rooeta1<1.566) continue;
      if (rooeta2>1.4442 && rooeta2<1.566) continue;
//define variables   
     if(!isdata) {
         weight=pixel_lumi*event_luminormfactor*event_Kfactor*event_weight; 
      }

      else if(isdata){
       weight=1.;
      }

       float in1=pholead_PhoSCRemovalPFIsoCharged;
        float in2=photrail_PhoSCRemovalPFIsoCharged;
        float ptin1=pholead_pt;
        float ptin2=photrail_pt;
        float sieiein1=pholead_sieie;
        float sieiein2=photrail_sieie;
        float etain1=fabs(pholead_SCeta);
        float etain2=fabs(photrail_SCeta);
 bool doswap=false;
//swap leading and subleading photon randomly
       if(do2d1side || do1p1f){
//event_pass12whoissiglike==0 leading photon signal region, subleading sideban
//for swaping -> all photons in signal band are in roovar1, all photons in sideband are in roovar2 
//                 if(event_pass12whoissiglike==0) { doswap=false;}
//                  if(event_pass12whoissiglike==1) { doswap=true;}
        }     
          
        if ( (randomgen->Uniform(0.,1.)>0.5)) doswap=true;

      if(!doswap){
      roovar1= pholead_PhoSCRemovalPFIsoCharged;
      roovar2= photrail_PhoSCRemovalPFIsoCharged;
//for eta division      
      rooeta1=fabs(pholead_SCeta);
      rooeta2=fabs(photrail_SCeta);
      roopt1=pholead_pt;
      roopt2=photrail_pt;
      roosieie1 = pholead_sieie;
      roosieie2 = photrail_sieie;
 }
     rooweight= weight;
//     h_weight->Fill(rooweight);
//furher variables for rootfile
    roorho=event_rho;
    roosigma=event_sigma;
    roonvtx=event_nRecVtx;
//swap
   if (doswap){
        float temp;
         temp=in1; in1=in2; in2=temp;
         temp=ptin1; ptin1=ptin2; ptin2=temp;
         temp=sieiein1; sieiein1=sieiein2; sieiein2=temp;
         temp=etain1; etain1=etain2; etain2=temp;
         event_pass12whoissiglike=!event_pass12whoissiglike;         
   }    

     roovar1=in1;
     roovar2=in2;
     roopt1=ptin1;
     roopt2=ptin2;
     roosieie1=sieiein1;
     roosieie2=sieiein2;
     rooeta1=etain1;
     rooeta2=etain2;


     TLorentzVector pho1;
     TLorentzVector pho2;
     TLorentzVector dipho;   
     Float_t dipho_pt=0.;
     Float_t dipho_mass=0;
//not in dead area
//compute diphoton momentum      
     pho1.SetPtEtaPhiE(pholead_pt,pholead_eta,pholead_phi, pholead_energy);
     pho2.SetPtEtaPhiE(photrail_pt,photrail_eta,photrail_phi, photrail_energy);
	     
     dipho=pho1+pho2;
     dipho_pt=dipho.Pt();
     dipho_mass=(pho1+pho2).M();

//EBEB
     if (rooeta1 < 1.4442 && rooeta2 < 1.4442){ 
         EB_diphomass->Fill(dipho_mass,weight);	     
         if(do2dff || do2dpp ||do2dside || do2dstd || dorcone){
	               if(roovar1 >= 0.) {    EB_temp->Fill(roovar1,weight);}
		       if(roovar2 >= 0.) {    EB_temp->Fill(roovar2,weight);}
}
	 tree->Fill();
	
/*	 if(do2d1side || do1p1f){
          	 if(event_pass12whoissiglike==0) {   EB_temp->Fill(roovar2,weight); EB_si->Fill(roosieie2,weight);tree->Fill();}
                  if(event_pass12whoissiglike==1) {    EB_temp->Fill(roovar1,weight);EB_si->Fill(roosieie1,weight);tree->Fill();}
        }*/
     }	
//one leg in EE     
      if (rooeta1 > 1.566 || rooeta2 > 1.566){
	      m1EE_diphopt->Fill(dipho_pt,weight);
	        m1EE_diphomass->Fill(dipho_mass,weight);
       	if(do2dff || do2dpp || do2dside || do2dstd ||dorcone){
	 	if(roovar1 >= 0. && rooeta1 > 1.566 ) {
		  	m1EE_temp->Fill(roovar1,weight);    
			m1EE_pt->Fill(roopt1,weight);
		}
	        if(roovar2 >= 0. && rooeta2 > 1.566 ) {
			m1EE_temp->Fill(roovar2,weight);  
		      	m1EE_pt->Fill(roopt2,weight);
		}

	}    
     	if(do2d1side || do1p1f){
//if pt fake > than pt true event_pass12whoissiglike==1 
     		if(event_pass12whoissiglike==1  && rooeta1 > 1.566) {    m1EE_temp->Fill(roovar1,weight); m1EE_si->Fill(pholead_sieie,weight);  m1EE_pt->Fill(pholead_pt,weight);
		 m1EE_diphomass->Fill(dipho_mass,weight);
		}
		if(event_pass12whoissiglike==0  && rooeta2 > 1.566) {    m1EE_temp->Fill(roovar2,weight); m1EE_si->Fill(roosieie2,weight);  m1EE_pt->Fill(roopt2,weight);
		 m1EE_diphomass->Fill(dipho_mass,weight);
		}
     	}
        	
     }
//EE EE
     if (rooeta1 > 1.566 && rooeta2 > 1.566){
	 EE_diphomass->Fill(dipho_mass,weight);
     	 if(do2dff || do2dpp || do2dside || do2dstd || dorcone){
       		if(roovar1 >= 0.) {    EE_temp->Fill(roovar1,weight);}
      		if(roovar2 >= 0.) {    EE_temp->Fill(roovar2,weight);}
	 }
        if(do2d1side || do1p1f){
        	if(event_pass12whoissiglike==0) {    EE_temp->Fill(roovar2,weight); EE_si->Fill(roosieie2,weight);}
  	        if(event_pass12whoissiglike==1) {    EE_temp->Fill(roovar1,weight);EE_si->Fill(roosieie1,weight);}
     	} 
     }
     pho1.Clear();
     pho2.Clear();
   }  
f->Write();
//h_weight->Draw();
 /*  
for(int bin=0; bin<nbins; bin++){ EB_diphomass->SetBinContent(bin+1,EB_diphomass->GetBinContent(bin+1)/(EB_diphomass->GetBinWidth(bin+1)));
//EB_diphomass->SetBinError(bin+1,EB_diphomass->GetBinError(bin+1)/(EB_diphomass->GetBinWidth(bin+1)));}
for(int bin=0; bin<nbins; bin++){ m1EE_diphomass->SetBinContent(bin+1,m1EE_diphomass->GetBinContent(bin+1)/(m1EE_diphomass->GetBinWidth(bin+1)));
	m1EE_diphomass->SetBinError(bin+1,m1EE_diphomass->GetBinError(bin+1)/(m1EE_diphomass->GetBinWidth(bin+1)));}
for(int bin=0; bin<nbins; bin++){ EE_diphomass->SetBinContent(bin+1,EE_diphomass->GetBinContent(bin+1)/(EE_diphomass->GetBinWidth(bin+1)));
	EE_diphomass->SetBinError(bin+1,EE_diphomass->GetBinError(bin+1)/(EE_diphomass->GetBinWidth(bin+1)));}
*/
//	EB_diphomass->Draw();
//plot
//  EB_temp->SaveAs("1dtruesigsig_signalmc_roovar2.root");
//  EB_temp->Scale(1.0/EB_temp->Integral()); EB_si->Scale(1.0/EB_si->Integral()); EB_diphomass->Scale(1.0/EB_diphomass->Integral());
   //Form() is wrapper from printf
//  const char* outfileEB=Form("rootfiles/%s_%s_EB.root",((do2d1side) ? "2d1side" : "1p1f"), ((isdata)? "data" : mc_name));  EB_temp->SaveAs(outfileEB);
/*  const char* diphomass_outfileEB=Form("rootfiles/%s_%s_EB_diphomass.root",((do2dstd) ? "2dstd" : "2dpp"), ((isdata)? "data" : mc_name));  EB_diphomass->SaveAs(diphomass_outfileEB);
  const char* outfilem1EE=Form("rootfiles/%s_%s_m1EE.root",((do2dstd) ? "2dstd" : "2dpp"), ((isdata)? "data" : mc_name));   m1EE_temp->SaveAs(outfilem1EE); 
   const char* diphomass_outfilem1EE=Form("rootfiles/%s_%s_m1EE_diphomass.root",((do2dstd) ? "2dstd" : "2dpp"), ((isdata)? "data" : mc_name));   m1EE_diphomass->SaveAs(diphomass_outfilem1EE);
//  EE_temp->Scale(1.0/EE_temp->Integral()); EE_si->Scale(1.0/EE_si->Integral());EE_diphomass->Scale(1.0/EE_diphomass->Integral());
   const char* outfileEE=Form("rootfiles/%s_%s_EE.root",((do2dstd) ? "2dstd" : "2dpp"), ((isdata)? "data" : mc_name));   EE_temp->SaveAs(outfileEE);
   const char* diphomass_outfileEE=Form("rootfiles/%s_%s_EE_diphomass.root",((do2dstd) ? "2dstd" : "2dpp"), ((isdata)? "data" : mc_name));  EE_diphomass->SaveAs(diphomass_outfileEE);
*/
}

