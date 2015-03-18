#define counttruth_cxx

#include "counttruth.h"
#include "TGraphErrors.h"
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
#include "TLegend.h"
#include "TPad.h"
#include "TVector3.h"
using namespace std;
//Float_t pixel_lumi= 19.7;
Float_t pixel_lumi= 1.;
//bool doswap=false;
TFile *f;TFile *f_h;
TTree *tree;
Float_t xbins[66];
Long64_t nentries;

TH1F* h_diphomass; TH1F* h_all;
TH1F* h_diphomassEB; TH1F* h_allEB;
TH1F* h_diphomassnEB; TH1F* h_allnEB;
TH1F* h_pEB; TH1F* h_pEBc;
TH1F* h_pnEB; TH1F* h_pnEBc;
TH1F* h_npEB; TH1F* h_npEBc;
TH1F* h_npnEB; TH1F* h_npnEBc;
int pEB=0; int pnEB=0; int npEB=0; int npnEB=0;

Float_t dpmass[]={0., 73.334,79.8851, 82.5001,85.0449, 87.5896,90.1619, 93.227,96.2921, 99.3572,103.387, 107.673,112.509, 118.,124.897, 133.411,144.482, 159.143,182.111, 222.319,2900.};

void counttruth::Loop()
{

	if (fChain == 0) return;
  	TH1F::SetDefaultSumw2(kTRUE);
	//variable binning
	const char* fitfilename="../forroofit/March16/150316counttruth.root";
	f = new TFile(fitfilename,"RECREATE");
	cout << fitfilename << endl;
	const int mbins=20;
	h_diphomass= new TH1F("h_diphomass","h_diphomass",mbins,dpmass);
	h_all=new TH1F("h_all","h_all",2,0.,2.);
	h_diphomassEB= new TH1F("h_diphomassEB","h_diphomassEB",mbins,dpmass);
//	h_allEB= new TH1F("h_allEB","h_allEB",mbins,dpmass);
	h_allEB=new TH1F("h_allEB","h_allEB",2,0.,2.);
	h_diphomassnEB= new TH1F("h_diphomassnEB","h_diphomassnEB",mbins,dpmass);
	h_allnEB=new TH1F("h_allnEB","h_allnEB",2,0.,2.);
	h_pEB= new TH1F("h_pEB","h_pEB",mbins,dpmass);
	h_pEBc= new TH1F("h_pEBc","h_pEBc",mbins,dpmass);
	h_pnEB= new TH1F("h_pnEB","h_pnEB",mbins,dpmass);
	h_pnEBc= new TH1F("h_pnEBc","h_pnEBc",mbins,dpmass);
	h_npEB= new TH1F("h_npEB","h_npEB",mbins,dpmass);
	h_npEBc= new TH1F("h_npEBc","h_npEBc",mbins,dpmass);
	h_npnEB= new TH1F("h_npnEB","h_npnEB",mbins,dpmass);
	h_npnEBc= new TH1F("h_npnEBc","h_npnEBc",mbins,dpmass);

	kSignal_l=-999; kSignal_t=-999;

	int gen_code_l=-999.;int gen_code_t=-999.;float geniso_l=-999.;float geniso_t=-999.;
	Float_t weight=-999.;
    
	nentries = fChain->GetEntriesFast();
    cout << nentries << endl;

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++)
    {

//    for (Long64_t jentry=0; jentry<2000;jentry++) {
		  Long64_t ientry = LoadTree(jentry);
		  if (ientry < 0) break;
	
		  nb = fChain->GetEntry(jentry);   nbytes += nb;
     //     PrintProgress(jentry);
			weight=pixel_lumi*event_luminormfactor*event_Kfactor*event_weight; 
	//get gen level informationi, kSignal_l ->true photon or kSignal_t true photon=setAlias in line
		 if(!isdata)
		 {
			 gen_code_l = pholead_PhoMCmatchexitcode;
			 geniso_l = pholead_GenPhotonIsoDR04;
	//also valide for sherpa as all set to status 2
			 if ((gen_code_l==1 || gen_code_l==2) && ((geniso_l<5.)&& (geniso_l >=0.))) {kSignal_l=1;}
			 else{kSignal_l=0;}
			 gen_code_t = photrail_PhoMCmatchexitcode;
			 geniso_t = photrail_GenPhotonIsoDR04;
			 if ((photrail_PhoMCmatchexitcode==1 ||photrail_PhoMCmatchexitcode==2)&&(geniso_t<5.)&& (geniso_l >=0.)){kSignal_t=1;}
			 else {kSignal_t=0;}
			// cout << "kSignal_l " << kSignal_l << "kSignal_t " << kSignal_t << endl;
		 }
         //TODO load tgrapherror files and get N
		 
		 TLorentzVector pho1; TLorentzVector pho2;TLorentzVector dipho;   
	     Float_t diphomass=0;
	//compute diphoton momentum      
		 pho1.SetPtEtaPhiE(pholead_pt,pholead_eta,pholead_phi, pholead_energy);
		 pho2.SetPtEtaPhiE(photrail_pt,photrail_eta,photrail_phi, photrail_energy);
		 diphomass=(pho1+pho2).M();
		 roodiphomass=diphomass;
		
		 rooeta1=fabs(pholead_SCeta);
		 rooeta2=fabs(photrail_SCeta);
		 //cut overflow away
		 if((pholead_isoh2ggvtx < 8.9999999) && (photrail_isoh2ggvtx < 8.9999999))
		 {
			 h_diphomass->Fill(diphomass,weight);
			 h_all->Fill(1.,weight);

			 if (rooeta1 < 1.4442 && rooeta2 < 1.4442)
			 {
			 			h_diphomassEB->Fill(diphomass,weight);
			 			h_allEB->Fill(1.,weight);
						//if GJet, diphotons or double fakes
						 if((kSignal_l==1 && kSignal_t==0) ||(kSignal_l==0 && kSignal_t==1))
						 {
							pEB++;
							//weight * 0.5 as only once looped over events and only one photon in it
							//otherwise double counted
							h_pEBc->Fill(1.,weight*0.5);
							h_pEB->Fill(diphomass,weight*0.5);
							npEB++;	
							h_npEBc->Fill(1.,weight*0.5);
							h_npEB->Fill(diphomass,weight*0.5);
						 }
						 else if(kSignal_l==1 && kSignal_t==1)
						 {
							pEB++;
							h_pEBc->Fill(1.,weight);
							h_pEB->Fill(diphomass,weight);
						 }
						 else if(kSignal_l==0 && kSignal_t==0)
						 {
							npEB++;	
							h_npEBc->Fill(1.,weight);
							h_npEB->Fill(diphomass,weight);
						 }
						 
					}
			 
				else if(rooeta1 > 1.566 || rooeta2 > 1.566)
				{
			 			h_diphomassnEB->Fill(diphomass,weight);
			 			h_allnEB->Fill(1.,weight);

						if((kSignal_l==1 && kSignal_t==0) ||(kSignal_l==0 && kSignal_t==1))
						{
							pnEB++;
							h_pnEBc->Fill(1.,weight*0.5);
							h_pnEB->Fill(diphomass,weight*0.5);
							npnEB++;	
							h_npnEBc->Fill(1.,weight*0.5);
							h_npnEB->Fill(diphomass,weight*0.5);
						 }
						 else if(kSignal_l==1 && kSignal_t==1)
						 {
							pnEB++;
							h_pnEBc->Fill(1.,weight);
							h_pnEB->Fill(diphomass,weight);
						 }
						 //else{
						 else if(kSignal_l==0 && kSignal_t==0)
						 {
							npnEB++;	
							h_npnEBc->Fill(1.,weight);
							h_npnEB->Fill(diphomass,weight);
						 }
			}
		 }	
    pho1.Clear();
	pho2.Clear();
	}
//LOOP over events finished
	
	cout << " all weighted events: h_all->Integral()  " << h_all->Integral()  <<  endl;
	cout << " all weighted events in diphomass   " << h_diphomass->Integral()  <<  endl;
    cout << "_________EBEB_____________" << endl;
	cout << "h_pEB integral "  << h_pEB->Integral() << endl;
	cout << "h_pEBc integral "  << h_pEBc->Integral() << endl;
	cout << "h_allEB integral "  << h_allEB->Integral() << endl;
	cout << "overall purity EB " << h_pEB->Integral()/h_allEB->Integral() << endl;
    cout << "__________sainty check____________" << endl;
	cout << "h_npEB integral "  << h_npEB->Integral() << endl;
	cout << "overall fakes EB " << h_npEB->Integral()/h_allEB->Integral() << endl;
    cout << "_________nEBEB_____________" << endl;
	cout << "h_pnEB integral "  << h_pnEB->Integral() << endl;
	cout << "h_allnEB integral "  << h_allnEB->Integral() << endl;
	cout << "overall purity nEB " << h_pnEB->Integral()/h_allnEB->Integral() << endl;
    cout << "__________sainty check____________" << endl;
	cout << "h_npnEB integral "  << h_npnEB->Integral() << endl;
	cout << "overall fakes nEB " << h_npnEB->Integral()/h_allnEB->Integral() << endl;
    cout << "______________________" << endl;
	cout <<" # pEB  " << pEB <<" # npEB " << npEB  <<" pnEB  " << pnEB <<" npnEB " << npnEB << endl;
	/*
	h_pEB->Scale(1.0/h_pEB->Integral());
	h_pnEB->Scale(1.0/h_pnEB->Integral());
	h_npnEB->Scale(1.0/h_npnEB->Integral());
	h_npEB->Scale(1.0/h_npEB->Integral());
	h_diphomass->Scale(1.0/h_diphomass->Integral());
	*/
	/*
	h_pEB->Scale(1.0/h_allEB->Integral());
	h_pnEB->Scale(1.0/h_allnEB->Integral());
	h_npnEB->Scale(1.0/h_allnEB->Integral());
	h_npEB->Scale(1.0/h_allEB->Integral());
	h_diphomass->Scale(1.0/h_diphomass->Integral());
	*/
	TCanvas * cmass= new TCanvas("cmass","cmass");
	 TLegend *leg2 = new TLegend(0.15,0.6,0.35,0.7);
	 cmass->cd();
	 cmass->SetLogx();
	 h_pEB->SetMarkerStyle(20);
	 h_pnEB->SetMarkerStyle(20);
	 h_diphomass->SetMarkerStyle(20);
	 h_diphomass->SetMarkerColor(kBlack);
	 h_pEB->SetMarkerColor(kRed);
	 h_pnEB->SetMarkerColor(kBlue);
	 h_diphomass->SetMarkerColor(kBlack);
	 h_pEB->SetLineColor(kRed);
	 h_pnEB->SetLineColor(kBlue);
	 h_diphomass->SetLineColor(kBlack);
	 h_pEB->GetYaxis()->SetRange(0.,1.);
 	 h_pEB->GetYaxis()->SetTitle("# Events");
     h_pEB->GetXaxis()->SetTitle("diphoton mass [GeV]");
	 leg2->SetFillColor(10);
	 leg2->SetLineColor(10);
	 leg2->SetTextSize(0.03);
	 leg2->SetTextFont(42);
	 leg2->AddEntry(h_pEB,"purity EBEB" ,"pl");
   	 leg2->AddEntry(h_pnEB,"purity nEBEB" ,"pl");
	 leg2->AddEntry(h_diphomass,"overall diphomass" ,"pl");
	 
	 h_diphomass->Draw();
     h_pnEB->Draw("SAME");
   	 h_pEB->Draw("SAME");
     leg2->Draw();
     cmass->SaveAs("../plots/diphomassall.png") ;
     cmass->SaveAs("../plots/diphomassall.root") ;
//TODO calculate purities per bin,store them in array an plot them with TGraphErrors
/*	Float_t purtruthEB[mbins+1];
    Float_t purtruthnEB[mbins+1];
	purtruthEB[0]=0;
	purtruthnEB[0]=0;
   for(int k=1; k< (mbins+1); k++){
	  purtruthEB[k]=h_pEB->GetBinContent(k)/h_diphomassEB->GetBinContent(k);
      purtruthnEB[k]=h_pnEB->GetBinContent(k)/h_diphomassnEB->GetBinContent(k);

	  cout << k << " purtruthEB[k] "<< purtruthEB[k] << endl;
	 cout <<  "dpmass[k]" << dpmass[k] << endl;
	}
    for(int k=1; k< (mbins+1); k++){
	  cout << k << " purtruthnEB[k] "<< purtruthnEB[k] << endl;
	}
     cout << "purtruthEB[20]" << purtruthEB[20] << endl; 
     cout << "purtruthnEB[20]" << purtruthnEB[20] << endl; 
	 */
	 cout << __LINE__ << endl;
    TGraphAsymmErrors* grEB=new TGraphAsymmErrors();
    TGraphAsymmErrors* grnEB=new TGraphAsymmErrors();
	 TCanvas *c1 = new TCanvas("c1","c1",0,0,1024,768);
	TLegend* leg = new TLegend(0.6,0.7,0.7,0.8);
	 leg->SetFillColor(10);
	 leg->SetLineColor(10);
	 leg->SetTextSize(0.03);
	 leg->SetTextFont(42);
	 c1->Clear();
	 c1->cd(1);
	 cout << __LINE__ << endl;
	 //TGraph *grEB = new TGraphErrors(mbins+1, dpmass, purtruthEB);
	 grEB->SetName("grEB");
	 //TGraph *grnEB = new TGraphErrors(mbins+1, dpmass, purtruthnEB);
	 grnEB->SetName("grnEB");
	 cout << __LINE__ << endl;
     grEB->Divide(h_pEB,h_diphomassEB,"cl=0.683 b(1,1) mode");
     grnEB->Divide(h_pnEB,h_diphomassnEB,"cl=0.683 b(1,1) mode");
	 TPad *pad1 = new TPad("pad1","",0,0,1,1);
     pad1->Draw();
     pad1->cd(); 
	 pad1->SetGridx();
	 pad1->SetGridy();
	 pad1->SetLogx();
 	 grEB->SetTitle("purity from counting");
	 grEB->SetMarkerStyle(20);
	 grEB->SetMarkerSize(1.3);
	 grEB->GetYaxis()->SetTitle("purity");
	 grEB->GetXaxis()->SetTitle("diphoton mass [GeV]");
	 grEB->GetYaxis()->SetRangeUser(0.2,1.);
	 grEB->GetXaxis()->SetLimits(0.,3000.);
	 grEB->GetYaxis()->SetTitleOffset(1.2);
	 grnEB->SetMarkerColor(kRed);
	 grnEB->SetLineColor(kRed);
	 grnEB->SetMarkerStyle(20);
	 grEB->SetMarkerColor(kBlue);
	 grEB->SetLineColor(kBlue);
	 grnEB->SetMarkerSize(1.3); 
	 pad1->SetTicks(0,2);
	 leg->AddEntry(grEB,"purity EB ","p");    
	 leg->AddEntry(grnEB,"purity nEB","p");    
	 grEB->Draw("AP");
	 grnEB->Draw("SAME P");
	 leg->Draw();
	 grEB->SaveAs("../plots/March16/grEB.root");
	 grnEB->SaveAs("../plots/March16/grnEB.root");
 	 c1->Print("../plots/March16/150316_purity_counting.root","root");
	 c1->Print("../plots/March16/150316_purity_counting.png","png");
//TODO Tgrapherrors give name
	f->Write();
	f->Close();
}



