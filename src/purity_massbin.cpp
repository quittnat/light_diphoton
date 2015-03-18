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
const int nbins=20;
const int split=2;
Double_t dpmass_all[nbins];
Double_t purity_all[nbins];
Double_t w2err_all[nbins];
int plot_pull()  
{
//plot projections-should be gaussian
//plot as a function of diphoton mass    
     TFile* fpull;
	fpull=new TFile("../plots/March16/grEB.root","READ");
   	assert(fpull);
    TGraphAsymmErrors* grpull=(TGraphAsymmErrors*)fpull->Get("grEB");
	assert(grpull);
	Double_t pull[nbins];
    cout <<"grpull->GetN() " <<grpull->GetN() << endl;
//	if((grpull->GetN())==(nbins-1)){
    
	for(int k=0; k< nbins;k++)
	{
			cout <<" grpull->GetY()[k] "<< grpull->GetY()[k] << "purity_all[k] " << purity_all[k] << "w2err_all[k] " << w2err_all[k] << endl;
			pull[k]=((grpull->GetY()[k])-purity_all[k])/w2err_all[k];
			cout << pull[k] << endl;
	}
//	}
	TCanvas* cpull =new TCanvas("cpull","cpull");
	cpull->Divide(2,1);
	cpull->cd(1);
	TPad *padpull = new TPad("padpull","",0,0,1,1);
	padpull->Draw();
	padpull->cd(); 
	padpull->SetGridx();
	padpull->SetGridy();
	padpull->SetLogx();
	  
	TGraphAsymmErrors *grpm = new TGraphAsymmErrors(nbins,dpmass_all, pull);
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
	TPad *padpull2 = new TPad("padpull2","",0,0,1,1);
	padpull2->Draw();
    TH1F* proj=new TH1F("proj","proj",nbins,-5.,5.);
    TF1 * g 	 =new TF1("g","gaus",-3.,3.);
    for(int i=0; i<=nbins;i++)
    {
	   proj->Fill(pull[i]);
    }

   // proj->GetXaxis()->SetRangeUser(-5.,5.);
	proj->Fit("g","R L");
	proj->Draw("");
	g->Draw("SAME");
//	cpull->Print(Form("../plots/March17/massbinned20/150317_truth_pull_EBEB_%u.root",nbins),"root");
//	cpull->Print(Form("../plots/March17/massbinned20/150317_truth_pull_EBEB_%u.png",nbins),"png");
	cpull->Print(Form("../plots/March17/massbinned20/150317_rcone_pull_EBEB_%u.root",nbins),"root");
	cpull->Print(Form("../plots/March17/massbinned20/150317_rcone_pull_EBEB_%u.png",nbins),"png");
	return 0;
}
//TODO etGGa_q definieren



//main script

int plot()
{
gROOT->Reset();
gROOT->SetStyle("Plain"); // possibilities: Default, Plain, Bold, Video, Pub
gStyle->SetOptStat(0);
gStyle->SetOptFit(1);
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
const char *path=Form("150317_truth_purity_sumw2erron_EBEB_massbin_20_range");
const char *dir=Form("../plots/March17/massbinned20/");
const char *dir2=Form("../plots/March16/massbinned20/");
//*********************************truthfit***********************************//
f1=new TFile(Form("%s%s_0_12.root",dir,path),"READ");
assert(f1);
TGraphErrors* gr1=(TGraphErrors*)f1->Get("Graph");
assert(gr1); 
const int	nbins1=gr1->GetN();
cout << "nbins1 " << nbins1<<endl;

f2=new TFile(Form("%s%s_10_21.root",dir,path),"READ");
assert(f2);
TGraphErrors* gr2=(TGraphErrors*)f2->Get("Graph");
assert(gr2);
//TODO for Pull get tgraph with sumw2 error on, but only if ready and not over canvas anymore
const int	nbins2=gr2->GetN();
cout << "nbins2 " << nbins2<<endl;

//*************************rcone and sideband *********************************//
const char *path2=Form("150316_rcone_sideb_purity_sumw2erron_EBEB_massbin_20_range");
f3=new TFile(Form("%s%s_0_11.root",dir2,path2),"READ");
assert(f3);
TGraphErrors* gr3=(TGraphErrors*)f3->Get("Graph");
assert(gr3);

f4=new TFile(Form("%s%s_10_21.root",dir2,path2),"READ");
assert(f4);
TGraphErrors* gr4=(TGraphErrors*)f4->Get("Graph");
assert(gr4);

//TODO build new array with all dpmass and purity from variable arrays from tgrapherrors after each other
//int end1=11;
//
//for(int k=0;k<nbins1; k++ )
/*
for(int k=0;k<nbins1-1; k++ )
{
	gr1->GetPoint(k,dpmass_all[k],purity_all[k]);
	w2err_all[k]=gr1->GetErrorY(k);
	cout << "k " << k << "dpmass " << dpmass_all[k] << "purity " << purity_all[k] << "w2err_all[k] " << w2err_all[k] <<endl;
}

cout << "second filling " << endl;
for(int k=1;k<nbins2-1; k++ )
{
	gr2->GetPoint(k,dpmass_all[k+nbins1-2],purity_all[k+nbins1-2]);
	w2err_all[k+nbins1-2]=gr2->GetErrorY(k);
	cout << "k " << k+nbins1-2 << "dpmass " << dpmass_all[k+nbins1-2] <<  "purity " << purity_all[k+nbins1-2] << "w2err_all[k] " << w2err_all[k+nbins1-2] <<endl;
} 
*/

for(int k=0;k<nbins1-2; k++ )
{
	gr3->GetPoint(k+1,dpmass_all[k],purity_all[k]);
	w2err_all[k]=gr3->GetErrorY(k+1);
	cout << "k " << k << "dpmass " << dpmass_all[k] << "purity " << purity_all[k] << "w2err_all[k] " << w2err_all[k] <<endl;
}

cout << "second filling " << endl;
for(int k=1;k<nbins2; k++ )
{
	gr4->GetPoint(k,dpmass_all[k+nbins1-3],purity_all[k+nbins1-3]);
	w2err_all[k+nbins1-3]=gr4->GetErrorY(k);
	cout << "k " << k+nbins1-3 << "dpmass " << dpmass_all[k+nbins1-3] <<  "purity " << purity_all[k+nbins1-3] << "w2err_all[k] " << w2err_all[k+nbins1-3] <<endl;
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
cres->Print("../plots/March17/150317_purity_EBEB_truth_rconecomp.root","root");
cres->Print("../plots/March17/150317_purity_EBEB_truth_rconecomp.png","png");

plot_pull();
return 0;
}

