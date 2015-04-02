#define counttruth_cxx

#include "counttruth.h"
#include "TGraphErrors.h"
#include <assert.h>
#include <iostream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TFormula.h"
#include <TH1D.h>
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
Double_t xbins[66];
Long64_t nentries;

TH1D* h_diphomass=NULL; TH1F* h_all=NULL;
TH1D* h_diphomassEB; TH1F* h_allEB;
TH1D* h_diphomassnEB; TH1F* h_allnEB;
TH1D* h_pEB; TH1F* h_pEBc;
TH1D* h_pnEB; TH1F* h_pnEBc;
TH1D* h_npEB; TH1F* h_npEBc;
TH1D* h_npnEB; TH1F* h_npnEBc;
//int pEB=0; int pnEB=0; int npEB=0; int npnEB=0;
float pEB=0; int pnEB=0; float npEB=0; int npnEB=0;

//Float_t dpmass[]={0., 73.334,79.8851, 82.5001,85.0449, 87.5896,90.1619, 93.227,96.2921, 99.3572,103.387, 107.673,112.509, 118.,124.897, 133.411,144.482, 159.143,182.111, 222.319,8000.};
Double_t dpmass[]={0., 73.3073,79.8459, 82.5331,85.1072, 87.6912,90.3288, 93.4158,96.5028, 99.5897,103.699, 107.964,112.849, 118.299,125.262, 133.815,144.895, 159.528,182.622, 222.844,2900.};

void counttruth::Loop()
{

	if (fChain == 0) return;
  	TH1D::SetDefaultSumw2(kTRUE);
	//variable binning
	const char* fitfilename="../plots/March30/150330counttruth_nooverflow.root";
	f = new TFile(fitfilename,"RECREATE");
	cout << fitfilename << endl;
	const int mbins=20;
//	h_diphomass= new TH1D("h_diphomass","h_diphomass",mbins,dpmass);
	h_diphomass= new TH1D("h_diphomass","h_diphomass",mbins,dpmass);
	h_all=new TH1F("h_all","h_all",3,0.,3.);
	h_diphomassEB= new TH1D("h_diphomassEB","h_diphomassEB",mbins,dpmass);
//	h_allEB= new TH1D("h_allEB","h_allEB",mbins,dpmass);
	h_allEB=new TH1F("h_allEB","h_allEB",2,0.,2.);
	h_diphomassnEB= new TH1D("h_diphomassnEB","h_diphomassnEB",mbins,dpmass);
	h_allnEB=new TH1F("h_allnEB","h_allnEB",2,0.,2.);
	h_pEB= new TH1D("h_pEB","h_pEB",mbins,dpmass);
	h_pEBc= new TH1F("h_pEBc","h_pEBc",2,0.,2.);
	h_pnEB= new TH1D("h_pnEB","h_pnEB",mbins,dpmass);
	h_pnEBc= new TH1F("h_pnEBc","h_pnEBc",mbins,dpmass);
	h_npEB= new TH1D("h_npEB","h_npEB",2,0.,2.);
	h_npEBc= new TH1F("h_npEBc","h_npEBc",2,0.,2.);
	h_npnEB= new TH1D("h_npnEB","h_npnEB",mbins,dpmass);
	h_npnEBc= new TH1F("h_npnEBc","h_npnEBc",2,0.,2.);

	kSignal_l=-999; kSignal_t=-999;

	int gen_code_l=-999.;int gen_code_t=-999.;float geniso_l=-999.;float geniso_t=-999.;
	Double_t weight=-999.;
    Double_t sumweights=0;
    Float_t sumweightsf=0;
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
			 gen_code_l=0; geniso_l=0;
			 gen_code_l = pholead_PhoMCmatchexitcode;
			 geniso_l = pholead_GenPhotonIsoDR04;
	//also valide for sherpa as all set to status 2
			 if ((gen_code_l==1 || gen_code_l==2) && ((geniso_l<5.)&& (geniso_l >=0.))) {kSignal_l=1;}
			 else{kSignal_l=0;}
			 gen_code_t=0; geniso_t=0;
			 gen_code_t = photrail_PhoMCmatchexitcode;
			 geniso_t = photrail_GenPhotonIsoDR04;
			 if ((photrail_PhoMCmatchexitcode==1 ||photrail_PhoMCmatchexitcode==2)&&(geniso_t<5.)&& (geniso_t >=0.)){kSignal_t=1;}
			 else {kSignal_t=0;}
		 }
         //TODO load tgrapherror files and get N
		 
		 TLorentzVector pho1; TLorentzVector pho2;TLorentzVector dipho;   
	     Double_t diphomass=0;
	//compute diphoton momentum      
		 pho1.SetPtEtaPhiE(pholead_pt,pholead_eta,pholead_phi, pholead_energy);
		 pho2.SetPtEtaPhiE(photrail_pt,photrail_eta,photrail_phi, photrail_energy);
		 diphomass=(pho1+pho2).M();
		 roodiphomass=diphomass;
		
		 rooeta1=fabs(pholead_SCeta);
		 rooeta2=fabs(photrail_SCeta);
		 //cut overflow away
		 if((pholead_isoh2ggvtx < (9.-1.e-5)) && (photrail_isoh2ggvtx < 9.-1.e-5))
		 {
			//TODO why different from h_all? make histo without binning
			 h_diphomass->Fill(diphomass,weight);
			 h_all->Fill(1,weight); //same as if directly in diphoton minitree

			 if (rooeta1 < 1.4442 && rooeta2 < 1.4442)
			 {
			 			h_diphomassEB->Fill(diphomass,weight);
			 			h_allEB->Fill(1.,weight);
//						sumweights+=weight;
//						sumweightsf+=weight;
						//if GJet, diphotons or double fakes
						 if((kSignal_l==1 && kSignal_t==0 ) ||(kSignal_l==0 && kSignal_t==1 ))
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
						 else if(kSignal_l==0 && kSignal_t==0  )
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
							h_pnEBc->Fill(1.,weight*.5);
							h_pnEB->Fill(diphomass,weight*.5);
							npnEB++;	
							h_npnEBc->Fill(1.,weight*.5);
						    h_npnEB->Fill(diphomass,weight*.5);
						 }
						 else if(kSignal_l==1 && kSignal_t==1)
						 {
							pnEB++;
							h_pnEBc->Fill(1.,weight);
							h_pnEB->Fill(diphomass,weight);
						 }
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
	
	cout << " all events  " << h_all->GetEntries()  <<  endl;
	cout << " all weighted events: h_all->Integral()  " << h_all->Integral()  <<  endl;
	cout << " all weights double " << sumweights  <<  endl;
	cout << " all weights  float" << sumweightsf  <<  endl;
	cout << " all events in diphomass   " << h_diphomass->GetEntries()  <<  endl;
	cout << " all weighted events in diphomass   " << h_diphomass->Integral()  <<  endl;
    
	cout << "_________EBEB_____________" << endl;
	cout << "h_pEB integral double "  << h_pEB->Integral() << endl;
	cout << "h_pEBc integral float "  << h_pEBc->Integral() << endl;
	cout << "h_allEB integral float "  << h_allEB->Integral() << endl;
	cout << "h_diphomassEB integral double "  << h_diphomassEB->Integral() << endl;
	cout << "overall purity EB double " << h_pEB->Integral()/h_diphomassEB->Integral() << endl;
	cout << "overall purity EB float " << h_pEBc->Integral()/h_allEB->Integral() << endl;
    cout << "_________nEBEB_____________" << endl;
	cout << "h_pnEB integral double "  << h_pnEB->Integral() << endl;
	cout << "h_pnEBc integral float "  << h_pnEBc->Integral() << endl;
	cout << "h_allnEB integral float "  << h_allnEB->Integral() << endl;
	cout << "h_allnEB integral double "  << h_diphomassnEB->Integral() << endl;
	cout << "overall purity nEB float " << h_pnEBc->Integral()/h_allnEB->Integral() << endl;
	cout << "h_diphomass!EB integral double "  << h_diphomassnEB->Integral() << endl;
	cout << "overall purity !EB double " << h_pnEB->Integral()/h_diphomassnEB->Integral() << endl;
    cout << "______not for comparison with other routines:_________" << endl;
	cout <<" # pEB  " << pEB <<" # npEB " << npEB  <<" pnEB  " << pnEB <<" npnEB " << npnEB << endl;
	
	
	
	//TCanvas * integral= new TCanvas("integral","integral");
//	integral->cd();   
//	h_all->Draw();
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
     
	 cmass->SaveAs("../plots/March30/diphomassall.png") ;
     cmass->SaveAs("../plots/March30/diphomassall.root") ;
    
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
	 grEB->SetName("grEB");
	 grnEB->SetName("grnEB");
     grEB->Divide(h_pEB,h_diphomassEB,"cl=0.683 b(1,1) mode");
     grnEB->Divide(h_pnEB,h_diphomassnEB,"cl=0.683 b(1,1) mode");
	 Double_t * yval=grEB->GetY();
     for(int k=1; k< (mbins+1); k++){
     cout << " k "<<  k <<" dpmass[k] " << dpmass[k] <<  " # of all EB *2" << (h_diphomassEB->GetBinContent(k)*2) << " # pure in EB *2" << (h_pEB->GetBinContent(k)*2)  << endl;
		 cout  <<  " ratio " <<  yval[k-1] << endl;
	 }
	 Double_t * nyval=grnEB->GetY();
     for(int k=1; k< (mbins+1); k++){
     	cout << " k "<<  k <<" dpmass[k] " << dpmass[k] <<  " # of all nEB *2" << (h_diphomassnEB->GetBinContent(k)*2) << " # pure in nEB *2" << (h_pnEB->GetBinContent(k)*2)  << endl;
		 cout  <<  " ratio " <<  nyval[k-1] << endl;
	 } 
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
	 grEB->SaveAs("../plots/March30/grEB.root");
	 grnEB->SaveAs("../plots/March30/grnEB.root");
 	 c1->Print("../plots/March30/150330_purity_counting.root","root");
	 c1->Print("../plots/March30/150330_purity_counting.png","png");
   	f->Write();
	f->Close();
	
}



