#define light_cxx

#include "light.h"
#include <assert.h>
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
//Float_t pixel_lumi= 19.7;
Float_t pixel_lumi= 1.;
bool fillroodata=false;
bool fill1dtemplate=true; bool areascaled=false;
bool doswap=false;
Bool_t save= true;
Bool_t b_diphomass = true;
TFile *f;
TTree *tree;
Float_t xbins[66];
void light::Loop()
{
if (fChain == 0) return;
randomgen = new TRandom3(0);
  TH1F::SetDefaultSumw2(kTRUE);

if(fill1dtemplate){
  Float_t xstart=400.;
        nbins=66;
        xbins[0]=0.;
        for(int i=0 ; i <nbins ; i++){
        if (xbins[i] < xstart)  xbins[i+1]= xbins[i]+10;
        else if ( xstart <= xbins[i] && xbins[i]< 3000)  xbins[i+1]=xbins[i]+100;
}
}
EBEB=false;m1EE=false; EEEE=false;fulletarange=false;
if(etarange=="EBEB") {seta="EBEB";  EBEB=true;}
else if(etarange=="m1EE") {seta="m1EE";  m1EE=true;}
else if(etarange=="EEEE") {seta="EEEE";  EEEE=true;}
else if(etarange=="fulletarange") {seta="fulletarange"; fulletarange=true;}
else {cout << "etarange not correct "<< endl;}
if(fill1dtemplate){
 h_iso= new TH1F("h_iso","h_iso",90, 0.,9.);
 h_sieie= new TH1F("h_sieie","h_sieie", 400, 0, 0.04);
 h_pt= new TH1F("h_pt","h_pt", 200, 0., 400.);
 h_diphopt= new TH1F("h_diphopt","h_diphopt", 1000, 0., 1000.);
 h_diphomass= new TH1F("h_diphomass","h_diphomass", nbins, xbins);
 h_weight =new TH1F("h_weight","h_weight", 500, 0.,50.); 
}
 do2dstd=false;do2dff=false; do2dpp=false; do2d1side=false;do2dside=false;do1p1f=false; do2drcone=false;
 kSignal_l=-999;kElectron_l=-999; kSignal_t=-999;kElectron_t=-999;
 if(realdata=="data")isdata=true;
 else if(realdata=="mc")isdata=false;
 else {cout << "no data or mc" << endl;}
 if(sel=="2dstd") do2dstd=true;
 else if(sel=="2dpp") do2dpp=true;
 else if(sel=="2dff") do2dff=true;
 else if(sel=="1p1f") do1p1f=true;
 else if(sel=="2drcone") do2drcone=true;
else if(sel=="2dside") do2dside=true;
else if(sel=="2d1side") do2d1side=true;
cout << sel << endl;
/*   if(treename== "Tree_2Dtruebkgbkg_template"){do2dff=true;isdata=false;name="2dff";}
//   else if(treename== "Tree_2Dstandard_selection"){do2dstd=true;isdata=false;name="2dstd";}
   else if(treename== "Tree_2Dstandard_selection"){do2dpp=true;name="2dpp";}
   else if(treename== "Tree_2Dsideband_template"){do2dside=true;isdata=false;name="2dside";}
   else if(treename== "Tree_2Drandomcone_template"){dorcone=true;isdata=false;name="2drcone";}
   else if(treename== "Tree_2Drandomconesideband_template"){do2d1side=true;isdata=false;name="2drconesideband";}
   //141024 keep in mind no gen prompt anymore!
   else if(treename== "Tree_2Dgenpromptplussideband_template"){do2d1side=true;isdata=false;name="2dsideband";}
   else if(treename== "Tree_2Dtruesigbkg_template"){do1p1f=true;isdata=false;name="1p1f";}
*/
if(fillroodata){
	 const char* fitfilename=Form("forroofit/%s%s_%s_%s_template.root", ((isdata)? "data" : "mc"),(sel.Data()), (seta.Data()), (doswap)? "subleadpho":"leadpho"); 
	f = new TFile(fitfilename,"RECREATE");
cout << fitfilename << endl;
	tree = new TTree("for_roofit","for_roofit");
	tree->Branch("roovar1",&roovar1, "roovar1/F");
	tree->Branch("roovar2",&roovar2, "roovar2/F");
	if(do2dff||do2dpp||do1p1f){
		tree->Branch("rootruth1",&rootruth1, "rootruth1/I");
	        tree->Branch("rootruth2",&rootruth2, "rootruth2/I");
	}

	tree->Branch("rooeta1",&rooeta1, "rooeta1/F");
	tree->Branch("rooeta2",&rooeta2, "rooeta2/F");
	tree->Branch("roopt1",&roopt1, "roopt1/F");
	tree->Branch("roopt2",&roopt2, "roopta2/F");
	tree->Branch("roosieie1",&roosieie1, "roosieie1/F");
	tree->Branch("roosieie2",&roosieie2, "roosieie2/F");
	tree->Branch("roodiphopt",&roodiphopt, "roodiphopt/F");
	tree->Branch("roorho",&roorho, "roorho/F");
	tree->Branch("roosigma",&roosigma, "roosigma/F");
	tree->Branch("roonvtx",&roonvtx, "roonvtx/F");
	tree->Branch("rooweight",&rooweight, "rooweight/F");

	roovar1 = -999;
    	roovar2 = -999;
    	rootruth1 = -999;
  	rootruth2 = -999;
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
 Long64_t nentries = fChain->GetEntriesFast();


 Float_t weight=0.; 
 Long64_t nbytes = 0, nb = 0;
 for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (pholead_pt<40. || photrail_pt<25.) continue;     
      if (fabs(pholead_SCeta)>1.4442 && fabs(pholead_SCeta)<1.566) continue;
      if (fabs(photrail_SCeta)>1.4442 && fabs(photrail_SCeta)<1.566) continue;
//define variables  
      
     if(!isdata) {
         weight=pixel_lumi*event_luminormfactor*event_Kfactor*event_weight; 
      }
     else if(isdata){
       weight=1.;
      }
     if(fillroodata){
     if(do2d1side || do1p1f){
//event_pass12whoissiglike==0 leading photon signal region, subleading sideban
//only for truth match for swaping -> all photons in signal band are in roovar1, all photons in sideband are in roovar2 
             if(event_pass12whoissiglike==0) { doswap=true;}
             if(event_pass12whoissiglike==1) { doswap=false;}
      }     
//   if ( (randomgen->Uniform(0.,1.)>0.5)) doswap=true;
      }
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
   if (doswap){
      roovar1= photrail_PhoSCRemovalPFIsoCharged;
      roovar2= pholead_PhoSCRemovalPFIsoCharged;
//for eta division      
      rooeta1=fabs(photrail_SCeta);
      rooeta2=fabs(pholead_SCeta);
      roopt1=photrail_pt;
      roopt2=pholead_pt;
      roosieie1 = photrail_sieie;
      roosieie2 = pholead_sieie;
 }
      rooweight= weight; roorho=event_rho;roosigma=event_sigma;  roonvtx=event_nRecVtx;
     rootruth1=-999;rootruth2=-999;
     TLorentzVector pho1; TLorentzVector pho2;TLorentzVector dipho;   
     Float_t diphopt=0.;  Float_t diphomass=0;
//compute diphoton momentum      
     pho1.SetPtEtaPhiE(pholead_pt,pholead_eta,pholead_phi, pholead_energy);
     pho2.SetPtEtaPhiE(photrail_pt,photrail_eta,photrail_phi, photrail_energy);
     diphopt=(pho1+pho2).Pt();
     diphomass=(pho1+pho2).M();
     roodiphopt=diphopt;
//EBEB
     if ((EBEB && (rooeta1 < 1.4442 && rooeta2 < 1.4442)) ||(m1EE && (rooeta1 > 1.566 || rooeta2 > 1.566))||(EEEE && (rooeta1 > 1.566 && rooeta2 > 1.566))|| fulletarange){
//if (rooeta1 < 1.4442 && rooeta2 < 1.4442){  
     if(fillroodata) {
		 int gen_code_l=-999;int gen_code_t=-999;float geniso_l=-999.;float geniso_t=-999.;
	     if(do2dff||do2dpp||do1p1f ){
	     assert (!isdata);
	     gen_code_l = pholead_PhoMCmatchexitcode;
             geniso_l = pholead_GenPhotonIsoDR04;
             if ((gen_code_l==1 || gen_code_l==2) && (geniso_l<5)) kSignal_l=1;
             else if (gen_code_l==4) kElectron_l=2;
             else kSignal_l=0;
    	     gen_code_t = photrail_PhoMCmatchexitcode;
             geniso_t = photrail_GenPhotonIsoDR04;
             if ((gen_code_t==1 || gen_code_t==2) && (geniso_t<5)) kSignal_t=1;
             else if (gen_code_t==4) kElectron_t=2;
             else kSignal_t=0;
	     rootruth1=kSignal_l;
	     rootruth2=kSignal_t;
 	     if(doswap){rootruth1=kSignal_t;rootruth2=kSignal_l;}
	     if(do2dff && (rootruth1==0 && rootruth2==0)) {tree->Fill();}
	     else if((do2dpp && rootruth1==1 && rootruth2==1)) {tree->Fill();}
	     else if(do1p1f && ((rootruth1==1 && rootruth2==0)|| (rootruth1==0 && rootruth2==1))){ tree->Fill();}
	     }


	     else tree->Fill();
      	}
      	if(fill1dtemplate){
         	h_diphomass->Fill(diphomass,weight); h_diphopt->Fill(diphopt,weight); 
      		if(do2dff || do2dpp ||do2dside || do2dstd || do2drcone){
        		if(roovar1 >= 0.) {    h_iso->Fill(roovar1,weight);}
                	if(roovar2 >= 0.) {   h_iso->Fill(roovar2,weight);}
     	 	}	
      	
        	 if(do2d1side || do1p1f){
         		if(event_pass12whoissiglike==0) {   h_iso->Fill(roovar2,weight); h_sieie->Fill(roosieie2,weight);}
           		if(event_pass12whoissiglike==1) {    h_iso->Fill(roovar1,weight);h_sieie->Fill(roosieie1,weight);}
    		    }	
        }
	}
     pho1.Clear();
     pho2.Clear();
 }
if(fillroodata) {f->Write();f->Close();}

//save output
if(fill1dtemplate) {
  if(b_diphomass){
       for(int bin=0; bin<nbins; bin++){ h_diphomass->SetBinContent(bin+1,h_diphomass->GetBinContent(bin+1)/(h_diphomass->GetBinWidth(bin+1)));
       h_diphomass->SetBinError(bin+1,h_diphomass->GetBinError(bin+1)/(h_diphomass->GetBinWidth(bin+1)));}
  }
  if(areascaled){
        h_iso->Scale(1.0/h_iso->Integral()); h_sieie->Scale(1.0/h_sieie->Integral()); h_diphomass->Scale(1.0/h_diphomass->Integral());h_diphopt->Scale(1.0/h_diphopt->Integral());
  } 
   //Form() is wrapper from printf
if(save){
  const char* outfile=Form("rootfiles/%s%s_%s_isotemp.root", ((isdata)? "data" : "mc"),(name.Data()), (seta.Data()));  h_iso->SaveAs(outfile);
  const char* diphomass_outfile=Form("rootfiles/%s%s_%s_diphomass.root",((isdata)? "data" : "mc"),(name.Data()), (seta.Data()));h_diphomass->SaveAs(diphomass_outfile);
  }
}
}

