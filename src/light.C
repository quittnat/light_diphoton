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
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TVector3.h"
using namespace std;
//Float_t pixel_lumi= 19.7;
Float_t pixel_lumi= 1.;
bool areascaled=false;bool eff_calc=false;
//bool doswap=false;
TFile *f;TFile *f_h;
TTree *tree;
Float_t xbins[66];
Long64_t nentries;
int u=0;int v=0;int y=0;int w=0;
TH1F* h_pvz; TH1F* h_h2ggvz; TH1F * h_allz;
TH1F* h_pv;TH1F* h_h2ggv;TH1F* h_truth; TH1F* h_all; TH2F* h_eratio;TH1F* h_geniso;TH1F* h_genisotruth; TH1F* h_isowv;


void PrintProgress(Long64_t entry){
    int step = 10; 
    // Adapt step in powers of 10 (every 10 below 100, every 100 below 1000, etc.)
// Method that prints the progress at reasonable frequency
    Long64_t power = 1;
    for ( size_t i=1; i<10; ++i ) { // up to 10^10...
        power *= 10; 
        if ( !(entry/power) ) break;
        step = power;
    }   
    if( !(entry%step) ) cout << ">>> Processing event # " << entry << endl;
    else if ( entry==nentries-1 )
        cout << ">>> Processing last event # " << entry << endl;
}



void light::Loop()
{
	if (fChain == 0) return;
	randomgen = new TRandom3(0);
  	TH1F::SetDefaultSumw2(kTRUE);
	//variable binning
	Float_t xstart=400.;
    nbins=66;
    xbins[0]=0.;
    for(int i=0 ; i <nbins ; i++){
        if (xbins[i] < xstart)  xbins[i+1]= xbins[i]+10;
        else if ( xstart <= xbins[i] && xbins[i]< 3000)  xbins[i+1]=xbins[i]+100;
	}
	EBEB=false;m1EE=false; EEEE=false;fulletarange=false;
	if(etarange=="EBEB") {seta="EBEB";  EBEB=true;}
	else if(etarange=="m1EE") {seta="m1EE";  m1EE=true;}
	else if(etarange=="EEEE") {seta="EEEE";  EEEE=true;}
	else if(etarange=="fulletarange") {seta="fulletarange"; fulletarange=true;}
	else {cout << "etarange not correct "<< endl;}
    cout << "doswap " << doswap << endl;
	const char* fitfilename=Form("../forroofit/March03/mc%s_%s%s.root",(sel.Data()), (seta.Data()),(doswap) ? "_rev" : "");
//	const char* fitfilename=Form("test.root");
	f = new TFile(fitfilename,"RECREATE");
	cout << fitfilename << endl;
	h_iso= new TH1F("h_iso","h_iso",90, 0.,9.);
	h_geniso= new TH1F("h_geniso","h_geniso",2000, -1000.,1000.);
	h_genisotruth= new TH1F("h_genisotruth","h_geniso",90, 0.,9.);
	h_truth=new TH1F("h_truth","h_truth",2,0.,2.);
	h_all=new TH1F("h_all","h_all",2,0.,2.);
	h_pv=new TH1F("h_pv","h_pv",2,0.,2.);
	h_h2ggv=new TH1F("h_h2ggv","h_h2ggv",2,0.,2.);
	h_allz=new TH1F("h_allz","h_allz",2,0.,2.);
	h_pvz=new TH1F("h_pvz","h_pvz",2,0.,2.);
	h_h2ggvz=new TH1F("h_h2ggvz","h_h2ggvz",2,0.,2.);
	h_isopv= new TH1F("h_isopv","h_isopv",90, 0.,9.);
	h_sieie= new TH1F("h_sieie","h_sieie", 400, 0, 0.04);
	h_isowv= new TH1F("h_isowv","h_isowv",90, 0.,9.);
	h_isogen= new TH1F("h_isogen","h_isogen",90, 0.,9.);
	h_pt= new TH1F("h_pt","h_pt", 200, 0., 400.);
	h_diphopt= new TH1F("h_diphopt","h_diphopt",1000,0.,1000.);
	h_h2ggv_diphopt= new TH1F("h_h2ggv_diphopt","h_h2ggv_diphopt",1000,0.,1000.);
	eff_h2ggv_diphopt= new TGraphAsymmErrors();
	eff_pv_diphopt= new TGraphAsymmErrors();
	h_pv_diphopt= new TH1F("h_pv_diphopt","h_pv_diphopt",1000,0.,1000.);
	h_diphomass= new TH1F("h_diphomass","h_diphomass", nbins, xbins);
	h_weight =new TH1F("h_weight","h_weight", 500, 0.,50.); 
    h_eratio=new TH2F("h_eratio","h_eratio",100, 0.,1.,2000,-1000.,1000.);
	do2dstd=false;do2dff=false; do2dpp=false; do2d1side=false;do2dside=false;do1p1f=false; do2drcone=false;
	kSignal_l=-999;kElectron_l=-999; kSignal_t=-999;kElectron_t=-999;
	if(realdata=="data")isdata=true;
	else if(realdata=="mc"){isdata=false;cout << "mc " << endl;}
	else {cout << "no data or mc" << endl;}
	
	tree = new TTree("for_roofit","for_roofit");
    tree->Branch("roovar1",&roovar1, "roovar1/F");
	tree->Branch("roovar2",&roovar2, "roovar2/F");
	tree->Branch("rootruth1",&rootruth1, "rootruth1/I");
    tree->Branch("rootruth2",&rootruth2, "rootruth2/I");
	tree->Branch("rooisopv1",&rooisopv1, "rooisopv1/F");
	tree->Branch("rooisopv2",&rooisopv2, "rooisopv2/F");

	tree->Branch("rooisowv1",&rooisowv1, "rooisowv1/F");
	tree->Branch("rooisowv2",&rooisowv2, "rooisowv2/F");
	tree->Branch("rooeta1",&rooeta1, "rooeta1/F");
	tree->Branch("rooeta2",&rooeta2, "rooeta2/F");
	tree->Branch("roopt1",&roopt1, "roopt1/F");
	tree->Branch("roopt2",&roopt2, "roopta2/F");
	tree->Branch("roosieie1",&roosieie1, "roosieie1/F");
	tree->Branch("roosieie2",&roosieie2, "roosieie2/F");
	tree->Branch("roodiphopt",&roodiphopt, "roodiphopt/F");
	tree->Branch("roodiphomass",&roodiphomass, "roodiphomass/F");
	tree->Branch("roorho",&roorho, "roorho/F");
	tree->Branch("roosigma",&roosigma, "roosigma/F");
	tree->Branch("roonvtx",&roonvtx, "roonvtx/F");
	tree->Branch("rooweight",&rooweight, "rooweight/F");

	roovar1 = -999; 	roovar2 = -999;
	rooisopv1= -999;	rooisopv2= -999;
	rootruth1 = -999;	rootruth2 = -999;
	rooeta1 = -999;		rooeta2 = -999;
	roopt1 = -999;		roopt2 = -999;
	roosieie1 = -999;	roosieie2 = -999;
	roorho = -999;		roosigma = -999;
	roonvtx = -999;		rooweight = -999;
    rooisowv1= -999;	rooisowv2= -999;
    roodiphomass= -999;	roodiphopt= -999;


	int gen_code_l=-999.;int gen_code_t=-999.;float geniso_l=-999.;float geniso_t=-999.;
	bool pv_match=false; bool h2ggv_match=false;
	Float_t weight=-999.;
    
	nentries = fChain->GetEntriesFast();
    cout << nentries << endl;

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
		  Long64_t ientry = LoadTree(jentry);
		  if (ientry < 0) break;
	
		  nb = fChain->GetEntry(jentry);   nbytes += nb;
     //     PrintProgress(jentry);
		  if(!isdata) {
			weight=pixel_lumi*event_luminormfactor*event_Kfactor*event_weight; 
		  }
		 else if(isdata){
		   	weight=1.;
		  }
	//get gen level information
		 if(!isdata){
			 gen_code_l = pholead_PhoMCmatchexitcode;
			 geniso_l = pholead_GenPhotonIsoDR04;
	//also valide for sherpa as all set to status 2
			 if ((gen_code_l==1 || gen_code_l==2) && ((geniso_l<5.)&& (geniso_l >=0.))) {kSignal_l=1;}
			 else if (gen_code_l==4) {kElectron_l=2; }
			 else{kSignal_l=0;}
			 gen_code_t = photrail_PhoMCmatchexitcode;
			 geniso_t = photrail_GenPhotonIsoDR04;
			 if ((photrail_PhoMCmatchexitcode==1 ||photrail_PhoMCmatchexitcode==2)&&(geniso_t<5.)&& (geniso_l >=0.)){kSignal_t=1;}
			 else if (gen_code_t==4) {kElectron_t=2; }
			 else {kSignal_t=0;}
			 rootruth1=kSignal_l; 
			 rootruth2=kSignal_t;
			// cout << "rootruth1 " << rootruth1 << "rootruth2 " << rootruth2 << endl;
		 }

//selection setting
	     if(sel=="2dstd") do2dstd=true;
		 else if((sel=="2dpp") && rootruth1==1 && rootruth2==1) {do2dpp=true;} 
         else if((sel=="1p1f") && ((rootruth1==1 && rootruth2==0)|| (rootruth1==0 && rootruth2==1))){ do1p1f=true;}
         else if((sel=="1f1p") && ((rootruth1==1 && rootruth2==0)|| (rootruth1==0 && rootruth2==1))){ do1f1p=true;}
		 else if(sel=="2drcone") do2drcone=true;
		 else if(sel=="2drconeside") do2drconeside=true;
	     else if(sel=="2dside") do2dside=true;
	     else if(sel=="2d1side") do2d1side=true;
		 else if((sel=="2dff")&&(rootruth1==0 && rootruth2==0)){do2dff=true;}
		 else {cout << "no selection " << endl;
			 cout << "rootruth1 " << rootruth1 << "rootruth2 " << rootruth2 << endl;}
		 if((do2d1side || do1f1p)){
	//event_pass12whoissiglike==0 leading photon signal region, subleading sideban
	//only for truth match for swaping -> all photons in signal band are in roovar1, all photons in sideband are in roovar2 
				 if(event_pass12whoissiglike==0) { doswap=true;}
				 if(event_pass12whoissiglike==1) { doswap=false;}
		 }     
		 else if(do2drconeside||do1p1f){//get rcone from case where one passes photon selection, one in sideband. put passing photons in first branch
				 if(event_pass12whoissiglike==0) { doswap=false;}
				 if(event_pass12whoissiglike==1) { doswap=true;}
		 }     
        else if(do2dpp ||do2dstd ||do2dff||do2dside ||do2drcone ) { if((randomgen->Uniform(0.,1.)>0.5)) doswap=true;}
	
		 if(!doswap){
			  if(do2dpp || do2dstd || do1p1f || do2dff || do2d1side || do2dside || do1f1p){
				  rooisopv1= pholead_PhoSCRemovalPFIsoCharged;
				  rooisopv2= photrail_PhoSCRemovalPFIsoCharged;
				  rooisowv1=pholead_worstiso;
				  rooisowv2=photrail_worstiso;
				  roovar1= pholead_isoh2ggvtx;
				  roovar2= photrail_isoh2ggvtx;
			  }
			  else if(do2drcone || do2drconeside){
				  //! be sure that Tree random cone selected
				  rooisopv1= pholead_PhoSCRemovalPFIsoCharged;
				  rooisopv2= photrail_PhoSCRemovalPFIsoCharged;
				  rooisowv1=pholead_worstrcone;
				  rooisowv2=photrail_worstrcone;
				  roovar1= pholead_rconeh2ggvtx;
				  roovar2= photrail_rconeh2ggvtx;
			  }	 
			  else cout << "selection wrong " << endl; 
			  rooeta1=fabs(pholead_SCeta);
			  rooeta2=fabs(photrail_SCeta);
			  roopt1=pholead_pt;
			  roopt2=photrail_pt;
			  roosieie1 = pholead_sieie;
			  roosieie2 = photrail_sieie;
		 }
		 if (doswap){
			 if(do1p1f || do2dpp || do2dff || do2dstd || do2d1side || do2dside || do1f1p ){
				 rooisopv1= photrail_PhoSCRemovalPFIsoCharged;
				  rooisopv2= pholead_PhoSCRemovalPFIsoCharged;
				  rooisowv1=photrail_worstiso;
				  rooisowv2=pholead_worstiso;
				  roovar1= photrail_isoh2ggvtx;
				  roovar2= pholead_isoh2ggvtx;
			 }
			 else if(do2drconeside || do2drcone){
				  rooisopv1= photrail_PhoSCRemovalPFIsoCharged;
				  rooisopv2= pholead_PhoSCRemovalPFIsoCharged;
				  rooisowv1=photrail_worstrcone;
				  rooisowv2=pholead_worstrcone;
				  roovar1= photrail_rconeh2ggvtx;
				  roovar2= pholead_rconeh2ggvtx;
			  }
			  else cout << "wrong selection " << endl; 
			  rooeta1=fabs(photrail_SCeta);
			  rooeta2=fabs(pholead_SCeta);
			  roopt1=photrail_pt;
			  roopt2=pholead_pt;
			  roosieie1 = photrail_sieie;
			  roosieie2 = pholead_sieie;
			  rootruth1=kSignal_t;rootruth2=kSignal_l;
		 }
		 rooweight= weight; roorho=event_rho;roosigma=event_sigma;  roonvtx=event_nRecVtx;
		 TLorentzVector pho1; TLorentzVector pho2;TLorentzVector dipho;   
		 Float_t diphopt=0.;  Float_t diphomass=0;
	//compute diphoton momentum      
		 pho1.SetPtEtaPhiE(pholead_pt,pholead_eta,pholead_phi, pholead_energy);
		 pho2.SetPtEtaPhiE(photrail_pt,photrail_eta,photrail_phi, photrail_energy);
		 diphopt=(pho1+pho2).Pt();
		 diphomass=(pho1+pho2).M();
		 roodiphopt=diphopt;
		 roodiphomass=diphomass;
		
		 //see if primary vertex and/or h2ggvtx matches gen vertex
		 if(!isdata){
				if(fabs((pholead_GenVz-primVtxz)) < 1.){ pv_match=true;}
				else {pv_match=false;w++;}
				if(fabs((pholead_GenVz-diphoton_h2ggvtx_Vz)) < 1.) {h2ggv_match=true;}
				else {h2ggv_match=false;y++;}
		 }
		 if ((EBEB && (rooeta1 < 1.4442 && rooeta2 < 1.4442)) ||(m1EE && (rooeta1 > 1.566 || rooeta2 > 1.566))||(EEEE && (rooeta1 > 1.566 && rooeta2 > 1.566))|| fulletarange){  
				 f->cd();
		      	 h_weight->Fill(weight);
				 h_diphomass->Fill(diphomass,weight);
				 h_diphopt->Fill(diphopt,weight); 
				 if(diphopt < 250.){
				 h_allz->Fill(1.,weight); 
				 }
				 if(pv_match && diphopt<250.){h_pvz->Fill(1.,weight);}
				 if(h2ggv_match && diphopt < 250.){h_h2ggvz->Fill(1.,weight);}
				 if(pv_match){h_pv_diphopt->Fill(diphopt,weight);u++;h_pv->Fill(1.,weight);}
				 if(h2ggv_match){h_h2ggv_diphopt->Fill(diphopt,weight);v++;h_h2ggv->Fill(1.,weight);}
				//vertex studies
				if((do2dff || do2dpp ||do2dside || do2dstd || do2drcone)){
						if((do2dpp || do2drcone) &&(randomgen->Uniform(0.,1.)>0.9)){
							tree->Fill();
						}
						h_all->Fill(1.,weight);
						 h_geniso->Fill(photrail_GenPhotonIsoDR04,weight);
						 h_iso->Fill(roovar1,weight);
						 h_isopv->Fill(rooisopv1,weight);
						 h_isowv->Fill(rooisowv1,weight);
				}	
				else if(do2d1side || do1f1p ||do2drconeside ||do1p1f){
					if(randomgen->Uniform(0.,1.)>0.5)	{
						tree->Fill();
					}
				    h_all->Fill(1.,weight);
					//tree->Fill();	
   					h_iso->Fill(roovar1,weight);h_isopv->Fill(rooisopv1,weight);h_isowv->Fill(rooisowv1,weight);
				}	 
 		 }

// eta range finished
// clear for next entry
	pho1.Clear();
	pho2.Clear();
	}
//LOOP over events finished
	
	cout << "  h_all->Integral()  " << h_all->Integral()  <<  endl;
	
//	  cout << "h_pv integral " << h_pv->Integral() << endl;
 // cout << "______________________" << endl;
//  	  cout << " h_hwggv  " << h_truth->GetBinContent(2)  <<  endl;
	   cout << "h_h2ggv integral "  << h_h2ggv->Integral() << endl;
	   cout << "h2ggv eff " << h_h2ggv->Integral()/h_all->Integral() << endl;
//	  cout << "h_genisotruth integral " << h_genisotruth->Integral() << endl;
	  cout << "h_pv integral " << h_pv->Integral() << endl;
	  cout << "pv  eff " << h_pv->Integral()/h_all->Integral() << endl;

  cout << "______________________" << endl;
	   cout << "h2ggv eff phopt < 250 GeV " << h_h2ggvz->Integral()/h_allz->Integral() << endl;
	  cout << "pv  eff phopt < 250 GeV" << h_pvz->Integral()/h_allz->Integral() << endl;
//	  cout <<" # 2dstd selection " << u <<" # truth matching " << v <<endl; //"x " << x << "y " << y << endl;
	  cout <<" eff: #pvmatch true  " << u <<" pv match false " << w  <<" h2ggmatch true  " << v <<" h2gg match false " << y << endl; //"x " << x << "y " << y << endl;
	if(eff_calc){
		 eff_pv_diphopt->Draw("AP");
		 eff_h2ggv_diphopt->Draw("SAME");
	}	 

		 	for(int bin=0; bin<nbins; bin++){ 
				h_diphopt->SetBinContent(bin+1,h_diphopt->GetBinContent(bin+1)/(h_diphopt->GetBinWidth(bin+1)));
				h_diphopt->SetBinError(bin+1,h_diphopt->GetBinError(bin+1)/(h_diphopt->GetBinWidth(bin+1)));
				h_h2ggv_diphopt->SetBinContent(bin+1,h_h2ggv_diphopt->GetBinContent(bin+1)/(h_h2ggv_diphopt->GetBinWidth(bin+1)));
				h_h2ggv_diphopt->SetBinError(bin+1,h_h2ggv_diphopt->GetBinError(bin+1)/(h_h2ggv_diphopt->GetBinWidth(bin+1)));
				h_pv_diphopt->SetBinContent(bin+1,h_pv_diphopt->GetBinContent(bin+1)/(h_pv_diphopt->GetBinWidth(bin+1)));
				h_pv_diphopt->SetBinError(bin+1,h_pv_diphopt->GetBinError(bin+1)/(h_pv_diphopt->GetBinWidth(bin+1)));
			}  
			if(eff_calc){
				eff_pv_diphopt->Divide(h_pv_diphopt,h_diphopt,"cl=0.683 b(1,1) mode");
				eff_h2ggv_diphopt->Divide(h_h2ggv_diphopt,h_diphopt,"cl=0.683 b(1,1) mode");
				eff_pv_diphopt->SetName("eff_pv_diphopt"); eff_pv_diphopt->Write();eff_h2ggv_diphopt->SetName("eff_h2ggv_diphopt");eff_h2ggv_diphopt->Write();                
				TCanvas * c1= new TCanvas("c1","c1");
				c1->cd();	
				eff_pv_diphopt->Draw("AP");
				TCanvas * c2= new TCanvas("c2","c2");
				c2->cd();
				eff_h2ggv_diphopt->Draw("AP");	
			}
//correct errors for diphomass
				for(int bin=0; bin<nbins; bin++){ 
					  h_diphomass->SetBinContent(bin+1,h_diphomass->GetBinContent(bin+1)/(h_diphomass->GetBinWidth(bin+1)));
					  h_diphomass->SetBinError(bin+1,h_diphomass->GetBinError(bin+1)/(h_diphomass->GetBinWidth(bin+1)));}
			
			if(areascaled){
				h_iso->Scale(1.0/h_iso->Integral()); h_sieie->Scale(1.0/h_sieie->Integral()); h_diphomass->Scale(1.0/h_diphomass->Integral());h_diphopt->Scale(1.0/h_diphopt->Integral());
				h_pv_diphopt->Scale(1.0/h_pv_diphopt->Integral()); h_h2ggv_diphopt->Scale(1.0/h_h2ggv_diphopt->Integral());
				h_isopv->Scale(1.0/h_isopv->Integral()); h_isogen->Scale(1.0/h_isogen->Integral());h_isowv->Scale(1.0/h_isowv->Integral());
			}	 
 f->Write();f->Close();

}



