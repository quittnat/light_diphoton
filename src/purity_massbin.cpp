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
#define DTTMFMT "%Y-%m-%d"
#define DTTMSZ 11
TString dir;
TString dir2;
const char * eta=NULL;
using namespace std; 

static char *getDtTm (char *buff) {
	   time_t t = time (0);
	   strftime (buff, DTTMSZ, DTTMFMT, localtime (&t));
	   return buff;
}
TString date;
int nbins=0;
const int split=2;

int plot_pull(TGraphAsymmErrors* grpull,double dpmass[], double purity[], double w2errL[],double w2errH[],Bool_t truthcheck)  
{
	  TH1F::SetDefaultSumw2(kTRUE);
//plot projections-should be gaussian
//plot as a function of diphoton mass    
	assert(grpull);
	//const int nbinspull=grpull->GetN();
	double *ppull=new double[nbins];
	double  *perrxl=new double[nbins];
	double *perrxh=new double[nbins];

    cout <<"# entries of counting true photons " <<grpull->GetN() << endl;

    if((grpull->GetN())==(nbins))
	{
    
		for(int k=0; k< nbins;k++)
		{   
			if(purity[k] <= (grpull->GetY()[k]))
			{
				ppull[k]=((grpull->GetY()[k])-purity[k])/w2errH[k];
			}
			else if(purity[k] > (grpull->GetY()[k]))
			{
				ppull[k]=((grpull->GetY()[k])-purity[k])/w2errL[k];
			}
		    perrxl[k]=grpull->GetErrorXlow(k);
		    perrxh[k]=grpull->GetErrorXhigh(k);
		//	cout <<" grpull->GetY()[k] "<< grpull->GetY()[k] << " purity[k] " << purity[k] << " w2errL[k] " << w2errL[k] << " w2errH[k] " << w2errH[k] << " pull[k] " << ppull[k] << " errxl " <<perrxl[k] << " errxh" << perrxh[k] <<endl;
			cout << "k " << k << " dpmass " << dpmass[k] << " grpull->GetY()[k] "<< grpull->GetY()[k] << " pull[k] " << ppull[k]  <<endl;
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
	TGraphErrors *grpm = new TGraphErrors(nbins,dpmass, ppull);
	grpm->SetTitle(Form("mgg %s",eta));
	grpm->SetMarkerStyle(20);
	grpm->SetMarkerSize(1.3);
	grpm->GetXaxis()->SetLimits(50.,3000.);
	grpm->GetYaxis()->SetTitle("Pull");
	grpm->GetXaxis()->SetTitle("diphoton mass [GeV]");
	grpm->GetYaxis()->SetTitleOffset(1.2);
	padpull->SetTicks(0,2);
	grpm->Draw("AP");

	cpull->cd(2);
	gStyle->SetOptFit(1);
    TPad *padpull2 = new TPad("padpull2","",0,0,1,1);
	padpull2->Draw();
	padpull2->cd();
    TH1F* proj=new TH1F("proj","proj",20,-20.,20.);
    TF1 * g 	 =new TF1("g","gaus",-20.,20.);
    for(int i=0; i<=nbins;i++)
    {
	   proj->Fill(ppull[i]);
    }
 proj->Sumw2();
   // proj->GetXaxis()->SetRangeUser(-5.,5.);
	proj->Fit("g","L ");
	proj->Draw("HIST");
	g->Draw("SAME");
	cpull->Print(Form("%s%s_%s_pull_%s_template_%u_masscut.root",dir2.Data(),(truthcheck)? "truth": "rcone_side",date.Data(),eta,nbins),"root");
	cpull->Print(Form("%s%s_%s_pull_%s_template_%u_masscut.png",dir2.Data(),(truthcheck)? "truth": "rcone_side",date.Data(),eta,nbins),"png");
	return 0;
}
//TODO etGGa_q definieren



//main script
int plot (Bool_t truth=kFALSE,TString eta_q="EBEB");
int plot(Bool_t truth, TString eta_q )
{	
gROOT->Reset();
gROOT->SetStyle("Plain"); // possibilities: Default, Plain, Bold, Video, Pub
gStyle->SetOptStat(111111);
gStyle->SetOptTitle(0);

TCanvas *cres = new TCanvas("cres","Graph",0,0,1024,768);
TLegend* leg = new TLegend(0.6,0.8,0.7,0.88);
leg->SetFillColor(10);
leg->SetLineColor(10);
leg->SetTextSize(0.03);
leg->SetTextFont(42);

/////////////count purity directly from diphotonminitree file
char buff[DTTMSZ];
date=getDtTm(buff);
TFile* fcp;TFile* f1;TFile* f2;TFile* f3;TFile* f4;
TFile* f11;TFile* f22;TFile* f33;TFile* f44;
eta=eta_q.Data();
TGraphAsymmErrors* grcp=NULL;

if(eta_q=="EBEB")
{
	fcp=new TFile("../plots/March30/grEB.root","READ");
	assert(fcp);
	grcp=(TGraphAsymmErrors*)fcp->Get("grEB");
}
else 	
{
	fcp=new TFile("../plots/March30/grnEB.root","READ");
	assert(fcp);
	grcp=(TGraphAsymmErrors*)fcp->Get("grnEB");
}
assert(grcp);
// ************************read in mass bins and purity

//const char *path=Form("purity_%s_truth_sumw2erron_%s_massbin_20_range",date.Data(),eta);
//const char *path=Form("purity_2015-03-31_truth_sumw2erron_%s_massbin_20_range",eta);
const char *path=Form("purity_2015-04-01_truth_sumw2erron_%s_massbin_20_range",eta);
const char *path2=Form("purity_2015-04-01_truth_sumw2erron_%s_massbin_20_range",eta);
//const char *path2=Form("purity_%s_truth_sumw2erron_%s_massbin_20_range",date.Data(),eta);
//const char *path2=Form("purity_%s_rcone_sideb_sumw2erron_%s_massbin_20_range",date.Data(),eta);
//const char *path2=Form("purity_2015-03-31_truth_sumw2erron_%s_massbin_20_range",eta);
//const char *path=Form("purity_2015-03-31_rcone_sideb_sumw2erron_%s_massbin_20_range",eta);
dir=Form("../plots/April1/%s/massbinned20_4bins_truthreweight/",eta);
//dir=Form("../plots/March31/%s/massbinned20_4bins_reweight/",eta);
//dir=Form("../plots/March31/%s/massbinned20_4bins_5binspacing_reweight/",eta);
dir2=Form("../plots/April1/%s/massbinned20_4bins_truthreweight/",eta);
TGraphAsymmErrors* gr_n1=new TGraphAsymmErrors(); TGraphAsymmErrors* gr_n2=new TGraphAsymmErrors();
TGraphAsymmErrors* gr_n3=new TGraphAsymmErrors(); TGraphAsymmErrors* gr_n4=new TGraphAsymmErrors();
int	nbins1=0; int	nbins2=0;
int	nbins3=0; int	nbins4=0;
//*********************************truthfit***********************************//

//TODO make more flexible for 100 bins 
//f1=new TFile(Form("%s%s_0_10.root",dir,path),"READ");
f1=new TFile(Form("%s%s_0_5.root",dir2.Data(),path2),"READ");
assert(f1);
TGraphAsymmErrors* gr1=(TGraphAsymmErrors*)f1->Get("Graph");
assert(gr1); 

f2=new TFile(Form("%s%s_5_10.root",dir2.Data(),path2),"READ");
//f2=new TFile(Form("%s%s_10_20.root",dir,path),"READ");
assert(f2);
TGraphAsymmErrors* gr2=(TGraphAsymmErrors*)f2->Get("Graph");
assert(gr2);

f3=new TFile(Form("%s%s_10_15.root",dir2.Data(),path2),"READ");
assert(f3);
TGraphAsymmErrors* gr3=(TGraphAsymmErrors*)f3->Get("Graph");
assert(gr3);
f4=new TFile(Form("%s%s_15_20.root",dir2.Data(),path2),"READ");
assert(f4);
TGraphAsymmErrors* gr4=(TGraphAsymmErrors*)f4->Get("Graph");
assert(gr4);

if(truth)
{
	gr_n1=gr1;
	gr_n2=gr2;
	gr_n3=gr3;
	gr_n4=gr4;
	nbins1=(gr1->GetN());
	nbins2=(gr2->GetN());
	nbins3=(gr3->GetN());
	nbins4=(gr4->GetN());
	cout << "nbins1 truth " << nbins1<< "nbins2 truth" << nbins2<<endl;
}
//*************************rcone and sideband *********************************//
f11=new TFile(Form("%s%s_0_5.root",dir2.Data(),path2),"READ");
//f11=new TFile(Form("%s%s_0_10.root",dir2,path2),"READ");
assert(f11);
TGraphAsymmErrors* gr11=(TGraphAsymmErrors*)f11->Get("Graph");
assert(gr11);
f22=new TFile(Form("%s%s_5_10.root",dir2.Data(),path2),"READ");
//f22=new TFile(Form("%s%s_0_10.root",dir2,path2),"READ");
assert(f22);
TGraphAsymmErrors* gr22=(TGraphAsymmErrors*)f22->Get("Graph");
assert(gr22);

f33=new TFile(Form("%s%s_10_15.root",dir2.Data(),path2),"READ");
assert(f33);
TGraphAsymmErrors* gr33=(TGraphAsymmErrors*)f33->Get("Graph");
assert(gr33);
f44=new TFile(Form("%s%s_15_20.root",dir2.Data(),path2),"READ");
assert(f44);
TGraphAsymmErrors* gr44=(TGraphAsymmErrors*)f44->Get("Graph");
assert(gr44);

if(!truth)
{
	gr_n1=gr11;
	gr_n2=gr22;
	gr_n3=gr33;
	gr_n4=gr44;
	nbins1=(gr11->GetN());
	nbins2=(gr22->GetN());
	nbins3=(gr33->GetN());
	nbins4=(gr44->GetN());
	cout << "nbins1 rcone sideband " << nbins1<< "nbins2" << nbins2<<endl;
}


nbins=nbins1+nbins2+nbins3+nbins4;
cout << "nbins " << nbins << endl;
//Double_t dpmass_all[nbins];
//Double_t purity_all[nbins];
//Double_t w2err_all[nbins];
double *pdpmass_all=new double[nbins];
double *ppurity_all=new double[nbins];
//double *pw2err_all=new double[nbins];
double *pw2errL_all=new double[nbins];
double *pw2errH_all=new double[nbins];
for(int k=0;k<nbins1; k++ )
{
	gr_n1->GetPoint(k,pdpmass_all[k],ppurity_all[k]);
	pw2errH_all[k]=gr_n1->GetErrorYhigh(k);
    pw2errL_all[k]=gr_n1->GetErrorYlow(k);
	cout << "k " << k << "dpmass " << pdpmass_all[k] << "purity " << ppurity_all[k] << "w2errH_all[k] " << pw2errH_all[k] << "w2errL_all[k] " << pw2errL_all[k] <<endl;
}
//TODO get low and high error
cout << "second filling " << endl;
for(int k=0;k<nbins2; k++ )
{
	gr_n2->GetPoint(k,pdpmass_all[k+nbins1],ppurity_all[k+nbins1]);
	pw2errH_all[k+nbins1]=gr_n2->GetErrorYhigh(k);
	pw2errL_all[k+nbins1]=gr_n2->GetErrorYlow(k);
	cout << "k " << k+nbins << "dpmass " << pdpmass_all[k+nbins1] << "purity " << ppurity_all[k+nbins1] << "w2errH_all[k] " << pw2errH_all[k+nbins1] << "w2errL_all[k] " << pw2errL_all[k+nbins1] <<endl;
} 

cout << "third filling " << endl;
for(int k=0;k<nbins3; k++ )
{
	gr_n3->GetPoint(k,pdpmass_all[k+nbins1+nbins2],ppurity_all[k+nbins1+nbins2]);
	pw2errH_all[k+nbins1+nbins2]=gr_n2->GetErrorYhigh(k);
	pw2errL_all[k+nbins1+nbins2]=gr_n2->GetErrorYlow(k);
	cout << "k " << k+nbins1+nbins2 << "dpmass " << pdpmass_all[k+nbins1+nbins2] << "purity " << ppurity_all[k+nbins1+nbins2] << "w2errH_all[k] " << pw2errH_all[k+nbins1+nbins2] << "w2errL_all[k] " << pw2errL_all[k+nbins1+nbins2] <<endl;
} 
cout << "fourth filling " << endl;
for(int k=0;k<nbins4; k++ )
{
	gr_n4->GetPoint(k,pdpmass_all[k+nbins1+nbins2+nbins3],ppurity_all[k+nbins1+nbins2+nbins3]);
	pw2errH_all[k+nbins1+nbins2+nbins3]=gr_n4->GetErrorYhigh(k);
	pw2errL_all[k+nbins1+nbins2+nbins3]=gr_n4->GetErrorYlow(k);
	cout << "k " << k+nbins1+nbins2+nbins3 << "dpmass " << pdpmass_all[k+nbins1+nbins2+nbins3] << "purity " << ppurity_all[k+nbins1+nbins2+nbins3] << "w2errH_all[k] " << pw2errH_all[k+nbins1+nbins2+nbins3] << "w2errL_all[k] " << pw2errL_all[k+nbins1+nbins2+nbins3] <<endl;
} 

cres->cd();

TPad *pad1 = new TPad("pad1","",0,0,1,1);
pad1->Draw();
pad1->cd(); 
pad1->SetGridx();
pad1->SetGridy();
pad1->SetLogx();

grcp->GetYaxis()->SetTitle("purity");
grcp->GetXaxis()->SetTitle("diphoton mass [GeV]");
grcp->GetXaxis()->SetLimits(10.,3000.);
//gr1->GetXaxis()->SetRange(50,3000);
grcp->GetYaxis()->SetRangeUser(.1,1.);
grcp->GetYaxis()->SetTitleOffset(1.2);
 
grcp->SetMarkerStyle(20);grcp->SetMarkerColor(kBlack); grcp->SetLineColor(kBlack); grcp->SetMarkerSize(1.3);
gr1->SetMarkerStyle(20);gr1->SetMarkerColor(kRed); gr1->SetLineColor(kRed); gr1->SetMarkerSize(1.3);
gr2->SetMarkerStyle(20);gr2->SetMarkerColor(kRed); gr2->SetLineColor(kRed); gr2->SetMarkerSize(1.3);
gr3->SetMarkerStyle(20);gr3->SetMarkerColor(kRed); gr3->SetLineColor(kRed); gr3->SetMarkerSize(1.3);
gr4->SetMarkerStyle(20);gr4->SetMarkerColor(kRed); gr4->SetLineColor(kRed); gr4->SetMarkerSize(1.3);
gr11->SetMarkerStyle(20);gr11->SetMarkerColor(kBlue); gr11->SetLineColor(kBlue); gr11->SetMarkerSize(1.3);
gr22->SetMarkerStyle(20);gr22->SetMarkerColor(kBlue); gr22->SetLineColor(kBlue); gr22->SetMarkerSize(1.3);
gr33->SetMarkerStyle(20);gr33->SetMarkerColor(kBlue); gr33->SetLineColor(kBlue); gr33->SetMarkerSize(1.3);
gr44->SetMarkerStyle(20);gr44->SetMarkerColor(kBlue); gr44->SetLineColor(kBlue); gr44->SetMarkerSize(1.3);
leg->AddEntry(grcp,Form("truth %s",eta),"p");    
//leg->AddEntry(gr1,"truth fit EBEB","p");    
//leg->AddEntry(gr11,"template fit EBEB","p");    
//leg->AddEntry(gr1,Form("truth fit no reweight %s",eta),"p");    
leg->AddEntry(gr11,"truth fit EBEB reweighted","p");    
//leg->AddEntry(gr1,Form("template fit EBEB not reweighted %s",eta),"p");    
grcp->SetTitle(Form("mgg %s purity",eta));

grcp->Draw("AP");
/*
gr1->Draw("P");
gr2->Draw(" P");
gr3->Draw("P");
gr4->Draw("P");
*/
gr11->Draw("P");
gr22->Draw(" P");
gr33->Draw("P");
gr44->Draw("P");

leg->Draw();
pad1->SetTicks(0,2);
pad1->Update();

cres->Update();
cres->Print(Form("%spurity_%s_%s_template_truth_weightingcomp.root",dir2.Data(),date.Data(),eta),"root");
cres->Print(Form("%spurity_%s_%s_template_truth_weightingcomp.png",dir2.Data(),date.Data(),eta),"png");
//cres->Print(Form("%spurity_%s_%s_truth_reweight_comp.root",dir,date.Data(),eta),"root");
//cres->Print(Form("%spurity_%s_%s_truth_rweweight_comp.png",dir,date.Data(),eta),"png");
plot_pull(grcp,pdpmass_all,ppurity_all,pw2errL_all,pw2errH_all,truth);
delete[] pdpmass_all;
delete[] ppurity_all;
//delete[] pw2err_all;
delete[] pw2errL_all;
delete[] pw2errH_all;

return 0;
}

