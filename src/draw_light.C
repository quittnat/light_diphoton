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
Bool_t logplot=kFALSE;
Boo_t compplots=kFALSE;
Bool_t diphomass=kFALSE;
Bool_t correlation=kTrue;
Bool_t isdata=kFALSE;
Bool_t ratioplot=kFALSE;
void draw_light(){
	gStyle->SetOptStat(0);
        TCanvas* c1= new TCanvas("c1","c1");
      	TCanvas* c2= new TCanvas("c2","c2");
      	TCanvas* c3= new TCanvas("c3","c3");
        TLatex a;
  	TFile* f1=NULL; TFile* f2=NULL;TFile* f3=NULL;TFile* f4=NULL; TFile* f5=NULL;TFile* f6=NULL; TFile* f7=NULL; TFile* f8=NULL;TFile* f9=NULL; 
	TH1F *h1=NULL;TH1F *h2=NULL;TH1F *h3=NULL;TH1F *h4=NULL;TH1F *h5=NULL;TH1F *h6=NULL;TH1F *h7=NULL;TH1F *h8=NULL;TH1F *h9=NULL;
if(compplots){
	 f1=new TFile("rootfiles/2dff_allmc_EBswapped.root","read");
	 f2=new TFile("rootfiles/2dpp_allmc_EBswapped.root","read");
	 h1=(TH1F*)f1->Get("EB_temp");
	 h2=(TH1F*)f2->Get("EB_temp");
         f4=new TFile("rootfiles/2dff_allmc_m1EEswapped.root","read");
         f5=new TFile("rootfiles/2dpp_allmc_m1EEswapped.root","read");
	 h4=(TH1F*)f4->Get("m1EE_temp");
         h5=(TH1F*)f5->Get("m1EE_temp");
	 f7=new TFile("rootfiles/2dff_allmc_EEswapped.root","read");
         f8=new TFile("rootfiles/2dpp_allmc_EEswapped.root","read");
	 h7=(TH1F*)f7->Get("EE_temp");
         h8=(TH1F*)f8->Get("EE_temp");
else if(diphomass){
	 f1=new TFile("rootfiles/2dstd_allmc_EB_diphomass.root","read");
	 f2=new TFile("rootfiles/2dpp_signalmc_EB_diphomass.root","read");
         f3=new TFile("rootfiles/2dstd_data_EB_diphomass.root","read");
	
	 h1=(TH1F*)f1->Get("EB_diphomass");
	 h2=(TH1F*)f2->Get("EB_diphomass");
         h3=(TH1F*)f3->Get("EB_diphomass");
           
         f4=new TFile("rootfiles/2dstd_allmc_m1EE_diphomass.root","read");
         f5=new TFile("rootfiles/2dpp_signalmc_m1EE_diphomass.root","read");
         f6=new TFile("rootfiles/2dstd_data_m1EE_diphomass.root","read");
         h4=(TH1F*)f4->Get("m1EE_diphomass");
         h5=(TH1F*)f5->Get("m1EE_diphomass");
         h6=(TH1F*)f6->Get("m1EE_diphomass");
         f7=new TFile("rootfiles/2dstd_allmc_EE_diphomass.root","read");
         f8=new TFile("rootfiles/2dpp_signalmc_EE_diphomass.root","read");
         f9=new TFile("rootfiles/2dstd_data_EE_diphomass.root","read");
         h7=(TH1F*)f7->Get("EE_diphomass");
         h8=(TH1F*)f8->Get("EE_diphomass");
         h9=(TH1F*)f9->Get("EE_diphomass"); 
	}
if(corrplots){}
if(ratiolot){
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
  //    h1->GetYaxis()->SetRangeUser(10,100000);
      h1->SetMarkerColor(kBlue+1); h1->SetMarkerStyle(20);h1->SetLineColor(kBlue+1);
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
       leg->AddEntry(h1,"2 Sieie sideband in full  MC" ,"p");
       leg->AddEntry(h2,"2ff fake photons in QCD/DYJets MC", "p");
       leg->AddEntry(h3,"2 Sieie sideband in  DATA", "p");

       leg->AddEntry(h1,"true pp photons in full MC" ,"p");
       leg->AddEntry(h2,"true ff photons in full  MC", "p");
  //     leg->AddEntry(h3,"2 photons  in  DATA", "p");
    
       leg->Draw();
        TLatex a;
        a.SetNDC();
        a.SetTextSize(0.03);
	TLatex b;b.SetNDC();b.SetTextSize(0.06);b.SetTextColor(kRed);
	b.DrawLatex(0.55,0.5,"PRELIMINARY  EB");
//	a.DrawLatex(0.55,0.6,"#splitline{CMS Internal}{#sqrt{s} = 8 TeV L = 19.7 fb^{-1}}");
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
//   TLine *r = new TLine(0.,ratio,9.,ratio);
  // r->SetLineColor(kRed);
   r1->Draw("ep");
//   r->Draw("SAME");
  // r1->SetMinimum(0.0);  // Define Y ..
 //  r1->SetMaximum(10.); // .. range
r1->GetYaxis()->SetTitle("MCtruth/MCfake");
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

}
else {
   TLegend* leg = new TLegend(0.55, 0.65, .9, .9);
   TLatex b;b.SetNDC();b.SetTextSize(0.06);b.SetTextColor(kRed);
   b.DrawLatex(0.55,0.5,"PRELIMINARY");

   c1->cd();
   if(logplot){
       	c1->SetLogy();
        }   
   c1->SetGridx();
   h1->SetMarkerColor(kBlue+1);h1->SetLineColor(kBlue+1);h1->SetMarkerStyle(20);
   if(diphomass){h1->SetTitle("both photons in EB");h1->GetXaxis()->SetTitle("Diphoton mass (GeV)"); }
   else{  h1->SetTitle("both photons in EB");    h1->GetXaxis()->SetTitle("Charged Iso (GeV)");}
   h2->SetMarkerColor(kRed+1);h2->SetLineColor(kRed+1); h2->SetMarkerStyle(20);
  //     h3->SetMarkerColor(kGreen-3); h3->SetMarkerStyle(20);h3->SetLineColor(kGreen-3);
   h1->Draw();
   h2->Draw("SAME");
    //   h3->Draw("SAME");
   leg->AddEntry(h1,"true pp photons in full MC" ,"p");
   leg->AddEntry(h2,"true ff photons in full  MC", "p");
   leg->Draw();
   if(isdata){   
 	 a.DrawLatex(0.55,0.6,"#splitline{CMS Internal}{#sqrt{s} = 8 TeV L = 19.7 fb^{-1}}");
   }

   c2->cd();
   if(logplot){
       c2->SetLogy();
	}
   c2->SetGridx();
   h4->SetMarkerColor(kBlue+1);h4->SetLineColor(kBlue+1);h4->SetMarkerStyle(20);
   if(diphomass){h4->SetTitle("min. one photon in EE, diphoton mass if both in signal band");h4->GetXaxis()->SetTitle("Diphoton mass (GeV)"); }
   else{  h4->SetTitle("min. one photon in EE, isolation for 2 candidates in signa");    h4->GetXaxis()->SetTitle("Charged Iso (GeV)");}
   h5->SetMarkerColor(kRed+1);h5->SetLineColor(kRed+1); h5->SetMarkerStyle(20);
   //h6->SetMarkerColor(kGreen-3); h6->SetMarkerStyle(20);h6->SetLineColor(kGreen-3);
   //  h6->Draw();
   h5->Draw();
   h4->Draw("SAME");
   leg->Draw();
   if(isdata){
  	a.DrawLatex(0.55,0.6,"#splitline{CMS Internal}{#sqrt{s} = 8 TeV L = 19.7 fb^{-1}}");
   }
   
   c3->cd();
   if(logplot){
       c3->SetLogy();
   }
   c3->SetGridx();
   h7->SetMarkerColor(kBlue+1); h7->SetMarkerStyle(20);h7->SetLineColor(kBlue+1);
   if(diphomass){h7->SetTitle("both photons in EE and signal band, diphoton mass");h7->GetXaxis()->SetTitle("Diphoton mass (GeV)");}
   else {h7->SetTitle("isolation for both photons in EE");  h7->GetXaxis()->SetTitle("Charged Iso (GeV)");}
   h8->SetMarkerColor(kRed+1); h8->SetMarkerStyle(20);h8->SetLineColor(kRed+1);
    //   h9->SetMarkerColor(kGreen-3);h9->SetMarkerStyle(20);h9->SetLineColor(kGreen-3);
 //      h9->Draw();
   h8->Draw();
   h7->Draw("SAME");
   leg->Draw();
   if(isdata){
  	a.DrawLatex(0.55,0.6,"#splitline{CMS Internal}{#sqrt{s} = 8 TeV L = 19.7 fb^{-1}}");
   }
}
  /////////////////////////////
       

//Save
  const char* outfileEB=Form("plots/sigbkg_comp_EB_%s_%s.root",((diphomass) ? "Lumi_diphomass" : "Lumi_isotemp"), ((logplot)? "log" : "lin"));  c1->SaveAs(outfileEB);  
   const char* outfilem1EE=Form("plots/sigbkg_comp_m1EE_%s_%s.root",((diphomass) ? "Lumi_diphomass" : "Lumi_isotemp"), ((logplot)? "log" : "lin"));  c2->SaveAs(outfilem1EE);
  const char* outfileEE=Form("plots/sigbkg_comp_EE_%s_%s.root",((diphomass) ? "Lumi_diphomass" : "Lumi_isotemp"), ((logplot)? "log" : "lin"));  c3->SaveAs(outfileEE);

}
