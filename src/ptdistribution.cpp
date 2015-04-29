//
//M Quittnat 
//For plotting of purity vs mass bins
//
#include "TROOT.h"
#include <iostream>
#include "TFile.h"
#include <TStyle.h>
#include "RooDataHist.h"
#include "RooDataHist.h"
#include "TCanvas.h"
#include "TPad.h"
#include  "TLegend.h"
#include "TAxis.h"
#include "TString.h"
#include "TH2F.h"
#include "TH1F.h"
#include <Riostream.h>
#include "TTree.h"
#include <stdlib.h>
#include "TF1.h"
#include <string.h>
#include "RooDataHist.h"
#define DTTMFMT "%Y-%m-%d"
#define DTTMSZ 11
TString dir;
const char * eta=NULL;
using namespace std; 

static char *getDtTm (char *buff) {
	   time_t t = time (0);
	   strftime (buff, DTTMSZ, DTTMFMT, localtime (&t));
	   return buff;
}
TString date;
int nbins=0;
int tnbins=0;
const int split=2;




//main script
int plot (TString eta_q="EBEB");
int plot(TString eta_q )
{	
gROOT->Reset();
gROOT->SetStyle("Plain"); // possibilities: Default, Plain, Bold, Video, Pub
gStyle->SetOptStat(111111);
gStyle->SetOptTitle(0);

TCanvas *cres = new TCanvas("cres","Graph",0,0,1024,768);
TLegend* leg = new TLegend(0.5,0.5,0.7,0.9);
leg->SetFillColor(10);
leg->SetLineColor(10);
leg->SetTextSize(0.03);
leg->SetTextFont(42);

/////////////count purity directly from diphotonminitree file
char buff[DTTMSZ];
date=getDtTm(buff);
eta=eta_q.Data();
TH1D* grcp=NULL;

dir=Form("../plots/April29/%s/massbinned20_4bins/",eta);
if(eta_q=="EBEB")
{
	TFile * fcp=new TFile(Form("%sptdata_wholemassrange.root",dir.Data()),"READ");
	assert(fcp);
	grcp=(TH1D*)fcp->Get("hallpt__roopt1");
}
assert(grcp);
    for(int bin=0; bin< grcp->GetNbinsX();bin++)
{
	grcp->SetBinContent(bin+1,grcp->GetBinContent(bin+1)/(grcp->GetBinWidth(bin+1)));
}
grcp->Scale(1.0/grcp->Integral());
// ************************read in mass bins and purity
TH1D* hpt0[5];
TH1D* hpt5[5];
TH1D* hpt10[5];
TH1D* hpt15[5];
//*********************************getpt***********************************//

//TODO make more flexible for 100 bins 
//f1=new TFile(Form("%s%s_0_10.root",dir,path),"READ");
for(int i=0;i<5;i++)
{
		TFile* f0=new TFile(Form("%shpt0%u.root",dir.Data(),i+1),"READ");
		assert(f0);
		TH1D* hpttemp;
//		htemp->SetName(histname0);
        hpttemp=(TH1D*)f0->Get("hpt__roopt1");
    //    assert(hpttemp);
	   hpt0[i]=hpttemp;
		TFile* f5=new TFile(Form("%shpt5%u.root",dir.Data(),i+1),"READ");
		assert(f5);
        TH1D* hpttemp5=(TH1D*)f5->Get("hpt__roopt1");
        assert(hpttemp5);
		hpt5[i]=hpttemp5;
		TFile* f10=new TFile(Form("%shpt10%u.root",dir.Data(),i+1),"READ");
		assert(f10);
        TH1D* hpttemp10=(TH1D*)f10->Get("hpt__roopt1");
        assert(hpttemp10);
		hpt10[i]=hpttemp10;
		TFile* f15=new TFile(Form("%shpt15%u.root",dir.Data(),i+1),"READ");
		assert(f15);
        TH1D* hpttemp15=(TH1D*)f15->Get("hpt__roopt1");
        assert(hpttemp15);
		hpt15[i]=hpttemp15;
       
	   	
        hpt0[i]->Scale(1.0/hpt0[i]->Integral());
        hpt5[i]->Scale(1.0/hpt5[i]->Integral());
        hpt10[i]->Scale(1.0/hpt10[i]->Integral());
        hpt15[i]->Scale(1.0/hpt15[i]->Integral());
    for(int bin=0; bin< hpt0[i]->GetNbinsX();bin++)
{
	hpt0[i]->SetBinContent(bin+1,hpt0[i]->GetBinContent(bin+1)/(hpt0[i]->GetBinWidth(bin+1)));
	hpt5[i]->SetBinContent(bin+1,hpt5[i]->GetBinContent(bin+1)/(hpt5[i]->GetBinWidth(bin+1)));
	hpt10[i]->SetBinContent(bin+1,hpt10[i]->GetBinContent(bin+1)/(hpt10[i]->GetBinWidth(bin+1)));
	hpt15[i]->SetBinContent(bin+1,hpt15[i]->GetBinContent(bin+1)/(hpt15[i]->GetBinWidth(bin+1)));
}
	  //TODO multiply by bin width and colorwheel different colors 
		}

//*************************plot *********************************//
cres->cd();
cres->SetLogx();
TPad *pad1 = new TPad("pad1", "pad1", 0, 0., 1, 1.0);
pad1->SetGridx();         // Vertical grid
pad1->SetGridy();
pad1->Draw();
pad1->cd(); 
pad1->SetLogy();

grcp->GetXaxis()->SetTitle("single photon pt [GeV]");
//grcp->GetYaxis()->SetRangeUser(10,1e5);
grcp->GetYaxis()->SetTitleOffset(1.2);
 
grcp->SetMarkerStyle(20);grcp->SetMarkerColor(kBlack); grcp->SetLineColor(kBlack); grcp->SetMarkerSize(1.3);

//leg->AddEntry(grcp,Form("truth %s",eta),"p");    
//leg->AddEntry(gr1,Form("truth fit rew %s",eta),"p");    
//leg->AddEntry(gr11,Form("template fit not rew%s",eta),"p");    
grcp->SetTitle(Form("single pt %s purity",eta));

grcp->Draw();

leg->AddEntry(grcp,Form("pt whole masssp. %s",eta),"p");    
for(int i=0;i<2;i++)
{
	hpt0[i]->SetMarkerStyle(20);hpt0[i]->SetMarkerColor(kRed-2*i); hpt0[i]->SetLineColor(kRed-2*i); hpt0[i]->SetMarkerSize(0.9);
	hpt5[i]->SetMarkerStyle(20);hpt5[i]->SetMarkerColor(kBlue-2*i); hpt5[i]->SetLineColor(kBlue-2*i); hpt5[i]->SetMarkerSize(0.9);
	hpt10[i]->SetMarkerStyle(20);hpt10[i]->SetMarkerColor(kCyan-2*i); hpt10[i]->SetLineColor(kCyan-2*i); hpt10[i]->SetMarkerSize(0.9);
	hpt15[i]->SetMarkerStyle(20);hpt15[i]->SetMarkerColor(kGreen-2*i); hpt15[i]->SetLineColor(kGreen-2*i); hpt15[i]->SetMarkerSize(0.9);
   leg->AddEntry(hpt0[i],Form("mbin 0+%u",i),"pi");    
   leg->AddEntry(hpt5[i],Form("mbin 5+%u",i),"pi");    
   leg->AddEntry(hpt10[i],Form("mbin 10+%u",i),"pi");    
   leg->AddEntry(hpt15[i],Form("mbin 15+%u",i),"pi");    
	hpt0[i]->Draw("SAME");
    hpt5[i]->Draw("SAME");
    hpt10[i]->Draw("SAME");
    hpt15[i]->Draw("SAME");
}

leg->Draw();
cres->Print(Form("%sptdistribution_comparison_%s_truth_comp.png",dir.Data(),eta),"png");
cres->Print(Form("%sptdistribution_comparison_%s_truth_comp.root",dir.Data(),eta),"root");

return 0;
}

