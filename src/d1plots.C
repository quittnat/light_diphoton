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
Bool_t logplot=kTRUE;
void draw_light(){
	gStyle->SetOptStat(1111111);
        TCanvas* c1= new TCanvas("c1","c1");
      	TCanvas* c2= new TCanvas("c2","c2");
      	TCanvas* c3= new TCanvas("c3","c3");
        TLatex a;
  	TFile* f1=NULL; TFile* f2=NULL;TFile* f3=NULL;TFile* f4=NULL; //TFile* f5=NULL;TFile* f6=NULL; TFile* f7=NULL; TFile* f8=NULL;TFile* f9=NULL; 
	TH1F *h1=NULL;TH1F *h2=NULL;TH1F *h11=NULL;TH1F *h4=NULL;TH1F *h5=NULL;TH1F *h22=NULL;TH1F *h7=NULL;TH1F *h8=NULL;TH1F *h44=NULL;TH1F *h55=NULL;TH1F *h88=NULL;TH1F *h77=NULL;
	 f1=new TFile("../forroofit/mc2dpp_EBEB.root","read");
	 f2=new TFile("../forroofit/mc1p1f_EBEB.root","read");
	 f3=new TFile("../forroofit/mc2drcone_EBEB.root","read");
	 f4=new TFile("../forroofit/mc2d1side_EBEB.root","read");
	 h1=(TH1F*)f1->Get("h_iso");
	 h2=(TH1F*)f2->Get("h_iso");
	 h4=(TH1F*)f1->Get("h_isopv");
	 h5=(TH1F*)f2->Get("h_isopv");
	 h7=(TH1F*)f1->Get("h_isowv");
	 h8=(TH1F*)f2->Get("h_isowv");
	
	 h11=(TH1F*)f3->Get("h_iso");
	 h22=(TH1F*)f4->Get("h_iso");
	 h44=(TH1F*)f3->Get("h_isopv");
	 h55=(TH1F*)f4->Get("h_isopv");
	 h77=(TH1F*)f3->Get("h_isowv");
	 h88=(TH1F*)f4->Get("h_isowv");
	 TLegend* leg = new TLegend(0.4, 0.7, .75, .9);
     leg->SetFillColor(kWhite);
   c1->cd();
   if(logplot){
       	c1->SetLogy();
        }   
   c1->SetGridx();
   h2->SetMarkerColor(kBlue+1);h1->SetMarkerStyle(20);
   h1->SetTitle("Both photons in EB, hgg vertex");h1->GetXaxis()->SetTitle("Charged Iso (GeV)");
   h1->SetMarkerColor(kRed+1); h2->SetMarkerStyle(20);
   h22->SetLineColor(kRed+1);h11->SetLineWidth(2);
   h11->SetLineColor(kBlue+1);h22->SetLineWidth(2);
   h1->Draw();
   h2->Draw("SAME");
   h11->Draw("SAMEHE2");
   h22->Draw("SAME HE2");
   leg->AddEntry(h1,"true pp photons in Sherpa MC" ,"p");
   leg->AddEntry(h2,"true f photons in GJet MC", "p");
   leg->AddEntry(h11,"random cone in Sherpa MC" ,"l");
   leg->AddEntry(h22,"sideband in GJet MC (npmatch)", "l");
   leg->Draw();

   c2->cd();
   if(logplot){
       c2->SetLogy();
	}
   c2->SetGridx();
   h4->SetMarkerColor(kBlue+1);h5->SetMarkerStyle(20);
   h4->SetTitle("Both photons in EB, primary vertex");    
   h4->GetXaxis()->SetTitle("Charged Iso (GeV)");
   h5->SetMarkerColor(kRed+1);h5->SetLineColor(kRed+1); h4->SetMarkerStyle(20);
   
   h55->SetLineColor(kRed+1);h44->SetLineWidth(2);
   h44->SetLineColor(kBlue+1);h55->SetLineWidth(2);
   h4->Draw();
   h5->Draw("SAME");
   h44->Draw("SAME H E2 ");
   h55->Draw("SAME H  E2 ");
   leg->Draw();
/////////////////////////////////////////   
   c3->cd();
   if(logplot){
       c3->SetLogy();
   }
   c3->SetGridx();
   h7->SetMarkerColor(kBlue+1); h8->SetMarkerStyle(20);h7->SetLineColor(kBlue+1);
   h7->SetTitle("Both photons in EB, worst vertex");  h7->GetXaxis()->SetTitle("Charged Iso (GeV)");
   h8->SetMarkerColor(kRed+1); h7->SetMarkerStyle(20);h8->SetLineColor(kRed+1);
   h88->SetLineColor(kRed+1);h77->SetLineWidth(2);
   h77->SetLineColor(kBlue+1);h88->SetLineWidth(2);
   h7->Draw();
   h8->Draw("SAME");
   h77->Draw("SAME HE2");
   h88->Draw("SAME HE2");
   leg->Draw();

  /////////////////////////////
       

//Save
  const char* outfilehgg=Form("../plots/rsb_comp_EB_hgg_%s.root", ((logplot)? "log" : "lin"));  c1->SaveAs(outfilehgg);  
   const char* outfilepv=Form("../plots/rsb_comp_EB_pv_%s.root", ((logplot)? "log" : "lin"));  c2->SaveAs(outfilepv);
  const char* outfilewv=Form("../plots/rsb_comp_EB_wv_%s.root", ((logplot)? "log" : "lin"));  c3->SaveAs(outfilewv);

}
