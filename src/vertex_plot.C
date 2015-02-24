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
Bool_t  iso=kFALSE;
Bool_t eff=kTRUE;
void draw_light(){
	gStyle->SetOptStat(0);
    TCanvas* c1= new TCanvas("c1","c1");
    TLatex b;
    TLatex res;
  	TFile* f1=NULL;TFile* f2=NULL; 
	TH1F *h1=NULL;TH1F *h2=NULL;TGraphAsymmErrors *h5=NULL;TGraphAsymmErrors *h4=NULL;
//	 f1=new TFile("../rootfiles/mc2dstd_EBEB__DiPhoJets.root","read");
 	// f1=new TFile("../forroofit/mc2dpp_fulletarange.root","read");
 	 f1=new TFile("../forroofit/mc2dpp_fulletarange_madgraph.root","read");
//	 f1=new TFile("../forroofit/mc1f1p_fulletarange.root","read");
	 f2=new TFile("../forroofit/mc2dpp_fulletarange.root","read");
//	f2=new TFile("../rootfiles/mc1p1f_EBEB__.root","read");
	assert(f1);assert(f2);
    if(iso){
	   h1=(TH1F*)f1->Get("h_iso");
	   h2=(TH1F*)f2->Get("h_isogen");
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
	   h5=(TGraphAsymmErrors*)f1->Get("eff_pv_diphopt");
	 //  TCanvas* c3= (TCanvas*)f1->Get("eff_pv_diphopt");
	  // h5=(TGraphAsymmErrors*)f1->GetPrimitive("eff_pv_diphopt");
	   TLegend* leg2 = new TLegend(0.12, 0.12, .4, .2);
      leg2->SetFillColor(kWhite);
      
	   c4->SetGridx();
	   h4->SetMarkerColor(kBlue+1);h4->SetLineColor(kBlue+1);h4->SetMarkerStyle(20);
	   //h4->GetXaxis()->SetRangeUser(0.,500.);
	   h5->SetMarkerColor(kRed+1);h5->SetLineColor(kRed+1);h5->SetMarkerStyle(20);
	   //h4->SetTitle("both photons in EB, true pp photons in DiPhotonBox+ DiPhotonJets madgraph");h4->GetXaxis()->SetTitle("Diphoton pt (GeV)"); 
	  // h4->SetTitle("true pp photons in DiPhotonBox+ Jets Sherpa, matching within 1 cm in z");h4->GetXaxis()->SetTitle("Diphoton pt (GeV)"); 
	  // h4->SetTitle("1 true 1 fake photon in GJets MC, matching within 1 cm in z");h4->GetXaxis()->SetTitle("Diphoton pt (GeV)"); 
	   h4->SetTitle("true pp photons in DiPhotonBox+ Jets madgraph, matching within 1 cm in z");h4->GetXaxis()->SetTitle("Diphoton pt (GeV)"); 
	   c4->Draw();
        h4->Draw("AP");
	   h5->Draw("P SAME");
	   leg2->AddEntry(h4,"vertex by legacy algorithm" ,"p");
	   leg2->AddEntry(h5,"primary vertex", "p");
	   leg2->Draw();
	   b.SetNDC();b.SetTextSize(0.06);b.SetTextColor(kRed);
       b.DrawLatex(0.12,0.4,"PRELIMINARY");
 	   res.SetNDC();res.SetTextSize(0.04);res.SetTextColor(kBlack);
     //  res.DrawLatex(0.12,0.3,"hggv eff 0.75, pv eff 0.53");
   //     res.DrawLatex(0.12,0.3,"hggv eff 0.79, pv eff 0.65");
       res.DrawLatex(0.12,0.3,"hggv eff 0.77, pv eff 0.62");
       c4->Update();
	}
	  /////////////////////////////
		   

//Save
    const char* outfileiso=Form("../plots/1p1f2ppisolation_pvvertex_comp_EB_%s.root", ((logplot)? "log" : "lin"));  c1->SaveAs(outfileiso);  
   const char* outfileiso2=Form("../plots/1p1f2ppisolation_pvvertex_comp_EB_%s.pdf", ((logplot)? "log" : "lin"));  c1->SaveAs(outfileiso2);
   if(eff){	const char* outfileeff="../plots/efficiency_vertex_comp_madgraph_fulletarange.root";  c4->SaveAs(outfileeff);
   	   const char* outfileeff2="../plots/efficiency_vertex_comp_madgraph_fulletarange.png";  c4->SaveAs(outfileeff2);}

}
