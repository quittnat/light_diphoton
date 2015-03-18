//
//M Quittnat 
//For plotting of purity vs mass bins
//
#include "TROOT.h"
#include <iostream>
#include "TFile.h"
#include <TStyle.h>
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
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
using namespace std;
//massbins
int nbins=0;
const int split=2;

int plot_pull(double dpmass[], double purity[], double w2err[],Bool_t truthcheck)  
{
//plot projections-should be gaussian
//plot as a function of diphoton mass    
     TFile* fpull;
	fpull=new TFile("../plots/March16/grEB.root","READ");
   	assert(fpull);
    TGraphAsymmErrors* grpull=(TGraphAsymmErrors*)fpull->Get("grEB");
	assert(grpull);
	//const int nbinspull=grpull->GetN();
//	Double_t pull[nbins];
//	Double_t errxl[nbins];
//	Double_t errxh[nbins];
	double *ppull=new double[nbins];
	double  *perrxl=new double[nbins];
	double *perrxh=new double[nbins];

    cout <<"# entries of counting true photons " <<grpull->GetN() << endl;

    if((grpull->GetN())==(nbins))
	{
    
		for(int k=0; k< nbins;k++)
		{
			ppull[k]=((grpull->GetY()[k])-purity[k])/w2err[k];
		    perrxl[k]=grpull->GetErrorXlow(k);
		    perrxh[k]=grpull->GetErrorXhigh(k);
			cout <<" grpull->GetY()[k] "<< grpull->GetY()[k] << " purity[k] " << purity[k] << " w2err[k] " << w2err[k] << " pull[k] " << ppull[k] << " errxl " <<perrxl[k] << " errxh" << perrxh[k] <<endl;
		}
    }
	TString canname;
	canname=Form("cpull_%s",(truthcheck)? "truth": "rcone_side");
	TCanvas *cpull=new TCanvas(canname,canname);
	cpull->Divide(2,1);
	cpull->cd(1);

	TPad *padpull = new TPad("padpull","",0,0,1,1);
	padpull->Draw();
	padpull->cd(); 
	padpull->SetGridx();
	padpull->SetGridy();
	padpull->SetLogx();
	 //TODO correct errors with error propagation for tgraph and th1 
	TGraphAsymmErrors *grpm = new TGraphAsymmErrors(nbins,dpmass, ppull,perrxl, perrxh, w2err);
//	grpm->SetTitle(Form("mgg %s",eta_q.Data()));
	grpm->SetTitle("mgg EBEB");
	grpm->SetMarkerStyle(20);
	grpm->SetMarkerSize(1.3);
	grpm->GetXaxis()->SetLimits(50.,3000.);
	grpm->GetYaxis()->SetTitle("Pull");
	grpm->GetXaxis()->SetTitle("diphoton mass [GeV]");
//		grpm->GetYaxis()->SetRangeUser(0.,0.1);
	grpm->GetYaxis()->SetTitleOffset(1.2);
	padpull->SetTicks(0,2);
	grpm->Draw("AP");

	cpull->cd(2);
	gStyle->SetOptFit(1);
    TPad *padpull2 = new TPad("padpull2","",0,0,1,1);
	padpull2->Draw();
	padpull2->cd();
    TH1F* proj=new TH1F("proj","proj",nbins,-5.,5.);
    TF1 * g 	 =new TF1("g","gaus",-3.,3.);
    for(int i=0; i<=nbins;i++)
    {
	   proj->Fill(ppull[i]);
    }

   // proj->GetXaxis()->SetRangeUser(-5.,5.);
	proj->Fit("g","L R");
	proj->Draw("");
	g->Draw("SAME");
//	cpull->Print(Form("../plots/March18/massbinned20/150318_truth_pull_EBEB_%u.root",nbins),"root");
//	cpull->Print(Form("../plots/March18/massbinned20/150318_truth_pull_EBEB_%u.png",nbins),"png");
	cpull->Print(Form("../plots/March18/massbinned20_4bins/150318_%s_pull_EBEB_%u.root",(truthcheck)? "truth": "rcone_side",nbins),"root");
	cpull->Print(Form("../plots/March18/massbinned20_4bins/150318_%s_pull_EBEB_%u.png",(truthcheck)? "truth": "rcone_side",nbins),"png");
	return 0;
}
//TODO etGGa_q definieren



//main script
int plot (Bool_t truth=kFALSE);
int plot(Bool_t truth )
{	
gROOT->Reset();
gROOT->SetStyle("Plain"); // possibilities: Default, Plain, Bold, Video, Pub
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);
gStyle->SetStatFont(63); 
gStyle->SetStatFontSize(30); 
gStyle->SetPadLeftMargin(0.15);
gStyle->SetPadBottomMargin(0.17);
// Define Canvas

TCanvas *cres = new TCanvas("cres","Graph",0,0,1024,768);
TLegend* leg = new TLegend(0.6,0.8,0.7,0.88);
leg->SetFillColor(10);
leg->SetLineColor(10);
leg->SetTextSize(0.03);
leg->SetTextFont(42);

/////////////count purity directly from diphotonminitree file
TFile* fcp;TFile* f1;TFile* f2;TFile* f3;TFile* f4;

fcp=new TFile("../plots/March16/grEB.root","READ");
assert(fcp);
TGraphAsymmErrors* grcp=(TGraphAsymmErrors*)fcp->Get("grEB");
assert(grcp);
// ************************read in mass bins and purity
const char *path=Form("150318_truth_purity_sumw2erron_EBEB_massbin_20_range");
const char *dir=Form("../plots/March18/massbinned20_4bins/");
TGraphErrors* gr_n1=new TGraphErrors();
TGraphErrors* gr_n2=new TGraphErrors();
//const char *dir2=Form("../plots/March16/massbinned20/");
int	nbins1=0;
int	nbins2=0;
//*********************************truthfit***********************************//

//TODO make more flexible for 100 bins 
f1=new TFile(Form("%s%s_0_11.root",dir,path),"READ");
assert(f1);
TGraphErrors* gr1=(TGraphErrors*)f1->Get("Graph");
assert(gr1); 

f2=new TFile(Form("%s%s_10_21.root",dir,path),"READ");
assert(f2);
TGraphErrors* gr2=(TGraphErrors*)f2->Get("Graph");
assert(gr2);
if(truth)
{
	gr_n1=gr1;
	gr_n2=gr2;
	nbins1=(gr1->GetN())-1;
	nbins2=(gr2->GetN())-1;
	cout << "nbins1 truth " << nbins1<< "nbins2 truth" << nbins2<<endl;
}
//*************************rcone and sideband *********************************//
const char *path2=Form("150318_rcone_sideb_purity_sumw2erron_EBEB_massbin_20_range");
f3=new TFile(Form("%s%s_0_11.root",dir,path2),"READ");
assert(f3);
TGraphErrors* gr3=(TGraphErrors*)f3->Get("Graph");
assert(gr3);

f4=new TFile(Form("%s%s_10_21.root",dir,path2),"READ");
assert(f4);
TGraphErrors* gr4=(TGraphErrors*)f4->Get("Graph");
assert(gr4);
if(!truth)
{
	gr_n1=gr3;
	gr_n2=gr4;
	nbins1=(gr2->GetN())-1;
	nbins2=(gr3->GetN())-1;
	cout << "nbins1 rcone sideband " << nbins1<< "nbins2" << nbins2<<endl;
}

nbins=nbins1+nbins2;
//Double_t dpmass_all[nbins];
//Double_t purity_all[nbins];
//Double_t w2err_all[nbins];
double *pdpmass_all=new double[nbins];
double *ppurity_all=new double[nbins];
double *pw2err_all=new double[nbins];
for(int k=0;k<nbins1; k++ )
{
	gr_n1->GetPoint(k,pdpmass_all[k],ppurity_all[k]);
	pw2err_all[k]=gr_n1->GetErrorY(k);
	cout << "k " << k << "dpmass " << pdpmass_all[k] << "purity " << ppurity_all[k] << "w2err_all[k] " << pw2err_all[k] <<endl;
}

cout << "second filling " << endl;
for(int k=0;k<nbins2; k++ )
{
	gr_n2->GetPoint(k,pdpmass_all[k+nbins1],ppurity_all[k+nbins1]);
	pw2err_all[k+nbins1]=gr4->GetErrorY(k);
	cout << "k " << k+nbins1 << "dpmass " << pdpmass_all[k+nbins1] <<  "purity " << ppurity_all[k+nbins1] << "w2err_all[k] " << pw2err_all[k+nbins1] <<endl;
} 

cres->cd();
TPad *pad1 = new TPad("pad1","",0,0,1,1);
pad1->Draw();
pad1->cd(); 
pad1->SetGridy();
pad1->SetGridy();
pad1->SetLogx();
gr1->GetYaxis()->SetTitle("purity");
gr1->GetXaxis()->SetTitle("diphoton mass [GeV]");
gr1->GetXaxis()->SetLimits(10.,3000.);
//gr1->GetXaxis()->SetRange(50,3000);
gr1->GetYaxis()->SetRangeUser(.1,1.);
gr1->GetYaxis()->SetTitleOffset(1.2);
 
grcp->SetMarkerStyle(20);grcp->SetMarkerColor(kBlack); grcp->SetLineColor(kBlack); grcp->SetMarkerSize(1.3);
gr1->SetMarkerStyle(20);gr1->SetMarkerColor(kRed); gr1->SetLineColor(kRed); gr1->SetMarkerSize(1.3);
gr2->SetMarkerStyle(20);gr2->SetMarkerColor(kRed); gr2->SetLineColor(kRed); gr2->SetMarkerSize(1.3);
gr3->SetMarkerStyle(20);gr3->SetMarkerColor(kBlue); gr3->SetLineColor(kBlue); gr3->SetMarkerSize(1.3);
gr4->SetMarkerStyle(20);gr4->SetMarkerColor(kBlue); gr4->SetLineColor(kBlue); gr4->SetMarkerSize(1.3);
leg->AddEntry(grcp,"truth  EBEB","p");    
leg->AddEntry(gr2,"truth fit EBEB","p");    
leg->AddEntry(gr3,"template fit EBEB","p");    
gr1->SetTitle("mgg EBEB purity");
gr1->Draw("AP");
grcp->Draw("P");
//gr1->Draw("AP");
gr2->Draw(" P");
gr3->Draw("P");
gr4->Draw("P");
leg->Draw();
pad1->SetTicks(0,2);
cres->Print("../plots/March18/massbinned20_4bins/150318_purity_EBEB_truth_rconecomp.root","root");
cres->Print("../plots/March18/massbinned20_4bins/150318_purity_EBEB_truth_rconecomp.png","png");

plot_pull(pdpmass_all,ppurity_all,pw2err_all,truth);

return 0;
}

