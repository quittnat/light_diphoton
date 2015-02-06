#include <iostream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TH1F.h>
#include <TFile.h>
#include <TLegend.h>
 #include <TAttLine.h>
#include <TAttMarker.h>
#include <TAttFill.h>
#include "TLine.h"
#include <TLatex.h>
#include "TGraphAsymmErrors.h"

TCanvas* c4;

Bool_t logplot=kTRUE;
Bool_t  iso=kTRUE;
Bool_t eff=kFALSE;
void draw_light(){
	gStyle->SetOptStat(0);
    TCanvas* c1= new TCanvas("c1","c1");
    TLatex b;
  	TFile* f1=NULL;TFile* f2=NULL; 
	TH1F *h1=NULL;TH1F *h2=NULL;TGraphAsymmErrors *h5=NULL;TGraphAsymmErrors *h4=NULL;
	 f1=new TFile("../rootfiles/mc2dstd_EBEB__DiPhoJets.root","read");
	f2=new TFile("../rootfiles/mc1p1f_EBEB__.root","read");
	assert(f1);
    if(iso){
	   h1=(TH1F*)f1->Get("h_iso");
	   h2=(TH1F*)f7->Get("h_isogen");
	   TLegend* leg = new TLegend(0.55, 0.65, .9, .9);
	   b.SetNDC();b.SetTextSize(0.06);b.SetTextColor(kRed);
	   b.DrawLatex(0.55,0.5,"PRELIMINARY");
	   c1->cd();
	   if(logplot){
			c1->SetLogy();
		}   
	   c1->SetGridx();
	   h1->SetMarkerColor(kBlue+1);h1->SetLineColor(kBlue+1);h1->SetMarkerStyle(20);
	//   h1->SetTitle("both photons in EB, true pp photons in DiPhotonBox+ DiPhotonJets madgraph");  
	     h1->SetTitle("EBEB, fake photon from Gamma Jets, true photon from DiPhotonBox+ DiPhotonJets madgraph");  
	     h1->GetXaxis()->SetTitle("Charged Iso (GeV)");
	   h2->SetMarkerColor(kRed+1);h2->SetLineColor(kRed+1); h2->SetMarkerStyle(20);
	   h1->Draw();
	   h2->Draw("SAME");
	   leg->AddEntry(h1,"true p with h2gg vertex" ,"p");
	   leg->AddEntry(h2,"true f with primary vertex", "p");
	   leg->Draw();
	}
	if(eff){
	  c4= new TCanvas("c4","c4");
      c4->cd();
//	   TCanvas* c2= (TCanvas*)f1->Get("eff_h2ggv_diphopt");
	   h4=(TGraphAsymmErrors*)f1->Get("eff_h2ggv_diphopt");
	 //  TCanvas* c3= (TCanvas*)f1->Get("eff_pv_diphopt");
	  // h5=(TGraphAsymmErrors*)f1->GetPrimitive("eff_pv_diphopt");
	   TLegend* leg2 = new TLegend(0.55, 0.65, .9, .9);
     
      
	   c4->SetGridx();
	   h4->SetMarkerColor(kBlue+1);h4->SetLineColor(kBlue+1);h4->SetMarkerStyle(20);
	   h4->SetTitle("both photons in EB, true pp photons in DiPhotonBox+ DiPhotonJets madgraph");h4->GetXaxis()->SetTitle("Diphoton pt (GeV)"); 
	   c4->Draw();
        h4->Draw("AP");
	  // h5->Draw("P SAME");
	   leg2->AddEntry(h4,"vertex by legacy algorithm" ,"p");
	   leg2->AddEntry(h5,"primary vertex", "p");
	   leg2->Draw();
       b.DrawLatex(0.55,0.5,"PRELIMINARY");
       c4->Update();
	}
	  /////////////////////////////
		   

//Save
    const char* outfileiso=Form("../plots/1p1f2ppisolation_pvvertex_comp_EB_%s.root", ((logplot)? "log" : "lin"));  c1->SaveAs(outfileiso);  
   const char* outfileiso2=Form("../plots/1p1f2ppisolation_pvvertex_comp_EB_%s.pdf", ((logplot)? "log" : "lin"));  c1->SaveAs(outfileiso2);
   if(eff){	const char* outfileeff="../plots/efficiency_vertex_comp_EB.root";  c4->SaveAs(outfileeff);}

}
