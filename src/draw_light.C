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
#include "tdrstyle.C"
#include <TAttFill.h>
#include "TLine.h"
#include <TLatex.h>
Bool_t logplot=kTRUE;
Bool_t diphomass=kFALSE;
void draw_light(){
	gStyle->SetOptStat(0);
//if(!diphomass){
	TFile* f1=new TFile("1d2dstd_allmc_roovar2.root","read");
	TFile* f2=new TFile("1dtruesigsig_signalmc_roovar2.root","read");
//	TFile* f3=new TFile("rootfiles/141031/2dstd_data_EB.root","read");
	 TH1F *h1=(TH1F*)f1->Get("EB_temp");
	 TH1F *h2=(TH1F*)f2->Get("EB_temp");
//	 TH1F *h3=(TH1F*)f3->Get("EB_temp");
/*
        TFile* f4=new TFile("rootfiles/rcone_allmc_m1EE.root","read");
        TFile* f5=new TFile("rootfiles/2dpp_signalmc_m1EE.root","read");
        TFile* f6=new TFile("rootfiles/rcone_data_m1EE.root","read");

         TH1F *h4=(TH1F*)f4->Get("m1EE_temp");
          TH1F *h5=(TH1F*)f5->Get("m1EE_temp");
          TH1F *h6=(TH1F*)f6->Get("m1EE_temp");

        TFile* f7=new TFile("rootfiles/rcone_allmc_EE.root","read");
        TFile* f8=new TFile("rootfiles/2dpp_signalmc_EE.root","read");
        TFile* f9=new TFile("rootfiles/rcone_data_EE.root","read");

         TH1F *h7=(TH1F*)f7->Get("EE_temp");
          TH1F *h8=(TH1F*)f8->Get("EE_temp");
          TH1F *h9=(TH1F*)f9->Get("EE_temp");	  
//}*/
/*
else {
	 TFile* f1=new TFile("rootfiles/2dstd_allmc_EB_diphomass.root","read");
	 TFile* f2=new TFile("rootfiles/2dpp_signalmc_EB_diphomass.root","read");
         TFile* f3=new TFile("rootfiles/2dstd_data_EB_diphomass.root","read");
	
	TH1F *h1=(TH1F*)f1->Get("EB_diphomass");
	TH1F *h2=(TH1F*)f2->Get("EB_diphomass");
        TH1F *h3=(TH1F*)f3->Get("EB_diphomass");
           
         TFile* f4=new TFile("rootfiles/2dstd_allmc_m1EE_diphomass.root","read");
         TFile* f5=new TFile("rootfiles/2dpp_signalmc_m1EE_diphomass.root","read");
         TFile* f6=new TFile("rootfiles/2dstd_data_m1EE_diphomass.root","read");
        TH1F *h4=(TH1F*)f4->Get("m1EE_diphomass");
        TH1F *h5=(TH1F*)f5->Get("m1EE_diphomass");
        TH1F *h6=(TH1F*)f6->Get("m1EE_diphomass");
         TFile* f7=new TFile("rootfiles/2dstd_allmc_EE_diphomass.root","read");
         TFile* f8=new TFile("rootfiles/2dpp_signalmc_EE_diphomass.root","read");
         TFile* f9=new TFile("rootfiles/2dstd_data_EE_diphomass.root","read");
        TH1F *h7=(TH1F*)f7->Get("EE_diphomass");
        TH1F *h8=(TH1F*)f8->Get("EE_diphomass");
        TH1F *h9=(TH1F*)f9->Get("EE_diphomass");
	 
	}
*/
 // ratio plots
        TH1F *r1 = (TH1F*)h2->Clone("r1");
	r1->Divide(h1);
        r1->SetMarkerStyle(20);
        r1->SetTitle("");
	TCanvas* c1= new TCanvas("c1","c1");
//	TCanvas* c2= new TCanvas("c2","c2");
//	TCanvas* c3= new TCanvas("c3","c3");
        TLegend* leg = new TLegend(0.55, 0.65, .9, .9);

	c1->cd();
	if(logplot){
	c1->SetLogy();
	}
 TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
   pad1->SetBottomMargin(0); // Upper and lower plot are joined
   pad1->SetGridx();         // Vertical grid
   pad1->Draw();             // Draw the upper pad: pad1
   pad1->cd();  
if(logplot){pad1->SetLogy();}            // pad1 becomes the current pad
     if(diphomass){  r1->GetXaxis()->SetTitle("Diphoton mass (GeV)"); h1->SetTitle("both photons in EB and signal band, diphoton mass");}
     else {  r1->GetXaxis()->SetTitle("Charged Iso (GeV)"); h1->SetTitle("isolation for both photons in EB and in signal band");}
//       h1->GetXaxis()->SetRangeUser(-100.,1000.);
      h1->GetYaxis()->SetRangeUser(10,100000);
      h1->SetMarkerColor(kBlack); h1->SetMarkerStyle(20);h1->SetLineColor(kBlack);
      h2->SetMarkerColor(kRed+1); h2->SetMarkerStyle(20);h2->SetLineColor(kRed+1);
 //     h3->SetMarkerStyle(20);  h3->SetMarkerColor(kGreen-3);h3->SetLineColor(kGreen-3);  
// get mean
      Double_t meanh1=h1->Integral(0,h1->GetNbinsX());
      Double_t meanh2=h2->Integral(0,h2->GetNbinsX());      
      Double_t ratio=meanh2/meanh1;
      cout << "h1 " << h1->GetNbinsX() <<"h2 " << h2->GetNbinsX() << endl;
      cout << "meanh1 " << meanh1 << "meanh2 " << meanh2 << "ratio " << ratio << endl;
      h1->Draw();
      h2->Draw("SAME");
    //  h2->Draw("SAME");
 
       leg->SetFillColor(0);
/*       leg->AddEntry(h1,"2 Sieie sideband in full  MC" ,"p");
       leg->AddEntry(h2,"2ff fake photons in QCD/DYJets MC", "p");
       leg->AddEntry(h3,"2 Sieie sideband in  DATA", "p");
*/
       leg->AddEntry(h1,"leading photon in full MC" ,"p");
       leg->AddEntry(h2,"leading true photon in GJets MC", "p");
  //     leg->AddEntry(h3,"2 photons  in  DATA", "p");
    
       leg->Draw();
        TLatex a;
        a.SetNDC();
        a.SetTextSize(0.03);
	a.DrawLatex(0.55,0.6,"#splitline{CMS Internal}{#sqrt{s} = 8 TeV L = 19.7 fb^{-1}}");
     // Do not draw the Y axis label on the upper plot and redraw a small
   // axis instead, in order to avoid the first label (0) to be clipped.
//   h1->GetYaxis()->SetLabelSize(0.);
//   TGaxis *axis = h1->GetXaxis();
 //  axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
 //  axis->SetLabelSize(15);
 //  axis->Draw();
// lower plot will be in pad
   c1->cd();          // Go back to the main canvas before defining pad2
   TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
   pad2->SetTopMargin(0);
   pad2->SetBottomMargin(0.4);
   pad2->SetTicky();
   pad2->SetGridx(); // vertical grid
   pad2->Draw();
   pad2->cd();
   TLine *r = new TLine(0.,ratio,9.,ratio);
   r->SetLineColor(kRed);
   r1->Draw("ep");
   r->Draw("SAME");
  // r1->SetMinimum(0.0);  // Define Y ..
 //  r1->SetMaximum(10.); // .. range
r1->GetYaxis()->SetTitle("MCtruth/fullMC");
   r1->GetYaxis()->SetNdivisions(505);
   r1->GetYaxis()->SetTitleSize(18);
   r1->GetYaxis()->SetTitleFont(43);
   r1->GetYaxis()->SetTitleOffset(1.05);
   r1->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   r1->GetYaxis()->SetLabelSize(15);
 // X axis ratio plot settings
   r1->GetXaxis()->SetTitleSize(20);
   r1->GetXaxis()->SetTitleFont(43);
   r1->GetXaxis()->SetTitleOffset(4.);
   r1->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   r1->GetXaxis()->SetLabelSize(15);


//////////////////////


/*

   c2->cd();
	if(logplot){
       c2->SetLogy();
	}
c2->SetGridx();
       h4->SetMarkerColor(kBlack);h4->SetLineColor(kBlack);h4->SetMarkerStyle(20);
      if(diphomass){h4->SetTitle("min. one photon in EE, diphoton mass if both in signal band");h4->GetXaxis()->SetTitle("Diphoton mass (GeV)"); }
      else{  h4->SetTitle("min. one photon in EE, isolation for 2 candidates in signal band");    h4->GetXaxis()->SetTitle("Charged Iso (GeV)");}
       h5->SetMarkerColor(kRed+1);h5->SetLineColor(kRed+1); h5->SetMarkerStyle(20);
       h6->SetMarkerColor(kGreen-3); h6->SetMarkerStyle(20);h6->SetLineColor(kGreen-3);
       h6->Draw();
       h5->Draw("SAME");
       h4->Draw("SAME");
leg->Draw();
  a.DrawLatex(0.55,0.6,"#splitline{CMS Internal}{#sqrt{s} = 8 TeV L = 19.7 fb^{-1}}");

       c3->cd();
       if(logplot){
       c3->SetLogy();
       }
c3->SetGridx();
       h7->SetMarkerColor(kBlack); h7->SetMarkerStyle(20);h7->SetLineColor(kBlack);
       if(diphomass){h7->SetTitle("both photons in EE and signal band, diphoton mass");h7->GetXaxis()->SetTitle("Diphoton mass (GeV)");}
      else {h7->SetTitle("isolation for both photons in EE and in signal band");  h7->GetXaxis()->SetTitle("Charged Iso (GeV)");}
       h8->SetMarkerColor(kRed+1); h8->SetMarkerStyle(20);h8->SetLineColor(kRed+1);
       h9->SetMarkerColor(kGreen-3);h9->SetMarkerStyle(20);h9->SetLineColor(kGreen-3);
       h9->Draw();
       h8->Draw("SAME");
       h7->Draw("SAME");
leg->Draw();
  a.DrawLatex(0.55,0.6,"#splitline{CMS Internal}{#sqrt{s} = 8 TeV L = 19.7 fb^{-1}}");
*/

  /////////////////////////////
       

//Save
//  const char* outfileEB=Form("plots/sig_2dstd_EB_%s_%s.root",((diphomass) ? "Lumi_diphomass" : "Lumi_isotemp"), ((logplot)? "log" : "lin"));  c1->SaveAs(outfileEB);  
 //  const char* outfilem1EE=Form("plots/sig_2dstd_m1EE_%s_%s.root",((diphomass) ? "Lumi_diphomass" : "Lumi_isotemp"), ((logplot)? "log" : "lin"));  c2->SaveAs(outfilem1EE);
 //  const char* outfileEE=Form("plots/sig_2dstd_EE_%s_%s.root",((diphomass) ? "Lumi_diphomass" : "Lumi_isotemp"), ((logplot)? "log" : "lin"));  c3->SaveAs(outfileEE);


	c1->SaveAs("truesig_subleadphotonsignalmc_EB.png");
  //      c2->SaveAs("bkg_2dside_m1EE_isotemp_lin.png");
    //    c3->SaveAs("bkg_2dside_EE_isotemp_lin.png");
}
