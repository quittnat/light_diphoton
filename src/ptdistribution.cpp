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
TLegend* leg = new TLegend(0.6,0.4,0.9,0.7);
leg->SetFillColor(10);
leg->SetLineColor(10);
leg->SetTextSize(0.03);
leg->SetTextFont(42);

/////////////count purity directly from diphotonminitree file
char buff[DTTMSZ];
date=getDtTm(buff);
eta=eta_q.Data();
TH1D* grcp=NULL;
TH1D* hside=NULL;

dir=Form("../plots/April29/%s/massbinned20_4bins_onlyrhosigmaetarew/",eta);
if(eta_q=="EBEB")
{
	TFile * fcp=new TFile(Form("%sptdata_wholemassrange.root",dir.Data()),"READ");
	assert(fcp);
	grcp=(TH1D*)fcp->Get("hallpt__roopt1");
	TFile * fside=new TFile(Form("%sptsideb_wholemassrange.root",dir.Data()),"READ");
	assert(fside);
	hside=(TH1D*)fside->Get("hallpt__roopt1");
}
assert(grcp);
assert(hside);
    for(int bin=0; bin< grcp->GetNbinsX();bin++)
{
	grcp->SetBinContent(bin+1,grcp->GetBinContent(bin+1)/(grcp->GetBinWidth(bin+1)));
	hside->SetBinContent(bin+1,hside->GetBinContent(bin+1)/(hside->GetBinWidth(bin+1)));
}
grcp->Scale(1.0/grcp->Integral());
hside->Scale(1.0/hside->Integral());
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
		TFile* f0=new TFile(Form("%shptdata0%u.root",dir.Data(),i+1),"READ");
		assert(f0);
		TH1D* hpttemp;
//		htemp->SetName(histname0);
        hpttemp=(TH1D*)f0->Get("hpt__roopt1");
    //    assert(hpttemp);
	   hpt0[i]=hpttemp;
		TFile* f5=new TFile(Form("%shptdata5%u.root",dir.Data(),i+1),"READ");
		assert(f5);
        TH1D* hpttemp5=(TH1D*)f5->Get("hpt__roopt1");
        assert(hpttemp5);
		hpt5[i]=hpttemp5;
		TFile* f10=new TFile(Form("%shptdata10%u.root",dir.Data(),i+1),"READ");
		assert(f10);
        TH1D* hpttemp10=(TH1D*)f10->Get("hpt__roopt1");
        assert(hpttemp10);
		hpt10[i]=hpttemp10;
		TFile* f15=new TFile(Form("%shptdata15%u.root",dir.Data(),i+1),"READ");
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
TPad *pad1 = new TPad("pad1", "pad1", 0, 0., 1, 1.0);
pad1->SetGridx();         // Vertical grid
pad1->SetGridy();
pad1->Draw();
pad1->cd(); 
pad1->SetLogy();
pad1->SetLogx();

grcp->GetXaxis()->SetTitle("single photon pt [GeV]");
//grcp->GetYaxis()->SetRangeUser(10,1e5);
grcp->GetYaxis()->SetTitleOffset(1.2);
 
grcp->SetMarkerStyle(20);grcp->SetMarkerColor(kBlack); grcp->SetLineColor(kBlack); grcp->SetMarkerSize(1.3);
hside->SetMarkerStyle(20);hside->SetMarkerColor(kYellow+2); hside->SetLineColor(kYellow+2); hside->SetMarkerSize(1.3);

//leg->AddEntry(grcp,Form("truth %s",eta),"p");    
//leg->AddEntry(gr1,Form("truth fit rew %s",eta),"p");    
//leg->AddEntry(gr11,Form("template fit not rew%s",eta),"p");    
grcp->SetTitle(Form("single pt %s purity",eta));

grcp->Draw("HIST");
hside->Draw("SAME HIST");

leg->AddEntry(grcp,Form("pt whole masssp. %s data",eta),"p");    
leg->AddEntry(hside,Form("pt whole masssp. %s sideband",eta),"p");    
for(int i=0;i<1;i++)
{
	hpt0[i]->SetMarkerStyle(20);hpt0[i]->SetMarkerColor(kRed-2*i); hpt0[i]->SetLineColor(kRed-2*i); hpt0[i]->SetMarkerSize(0.9);
	hpt5[i]->SetMarkerStyle(20);hpt5[i]->SetMarkerColor(kBlue-2*i); hpt5[i]->SetLineColor(kBlue-2*i); hpt5[i]->SetMarkerSize(0.9);
	hpt10[i]->SetMarkerStyle(20);hpt10[i]->SetMarkerColor(kCyan-2*i); hpt10[i]->SetLineColor(kCyan-2*i); hpt10[i]->SetMarkerSize(0.9);
	hpt15[i]->SetMarkerStyle(20);hpt15[i]->SetMarkerColor(kGreen-2*i); hpt15[i]->SetLineColor(kGreen-2*i); hpt15[i]->SetMarkerSize(0.9);
   leg->AddEntry(hpt0[i],Form("mbin 0+%u",i),"pi");    
   leg->AddEntry(hpt5[i],Form("mbin 5+%u",i),"pi");    
   leg->AddEntry(hpt10[i],Form("mbin 10+%u",i),"pi");    
   leg->AddEntry(hpt15[i],Form("mbin 15+%u",i),"pi");    
	hpt0[i]->Draw("SAME HIST");
    hpt5[i]->Draw("SAME HIST");
    hpt10[i]->Draw("SAME HIST");
    hpt15[i]->Draw("SAME HIST");
}

leg->Draw();
cres->Print(Form("%sptdistributiondata5stepsHIST_comparison_%s.png",dir.Data(),eta),"png");
cres->Print(Form("%sptdistributiondataa5stepsHIST_comparison_%s.root",dir.Data(),eta),"root");

return 0;
}

