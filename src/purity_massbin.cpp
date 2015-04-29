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
int tnbins=0;
const int split=2;

//int plot_pull(TGraphAsymmErrors* grpull,double dpmass[], double purity[], double w2errL[],double w2errH[],Bool_t truthcheck)  
int plot_pull(TGraphAsymmErrors* grpull,double dpmass[], double purity[], double w2err[],Bool_t truthcheck)  
{
	cout << "truthcheck "<<  truthcheck << endl;
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
		    ppull[k]=((grpull->GetY()[k])-purity[k])/w2err[k];
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
    TH1F* proj=new TH1F("proj","proj",10,-5.,5.);
    TF1 * g 	 =new TF1("g","gaus",-5.,5.);
    for(int i=0; i<=nbins;i++)
    {
	   proj->Fill(ppull[i]);
    }
 proj->Sumw2();
   // proj->GetXaxis()->SetRangeUser(-5.,5.);
	proj->Fit("g","L ");
	proj->Draw("HIST");
	g->Draw("SAME");
//	cpull->Print(Form("%s%s_%s_pull_%s_template_%u_masscut.root",dir2.Data(),(truthcheck)? "truth": "rcone_side",date.Data(),eta,nbins),"root");
//	cpull->Print(Form("%s%s_%s_pull_%s_template_%u_masscut.png",dir2.Data(),(truthcheck)? "truth": "rcone_side",date.Data(),eta,nbins),"png");
cpull->Print(Form("%s%s_%s_pull_%s_template_%u.root",dir2.Data(),(truthcheck)? "truth": "rcone_side",date.Data(),eta,nbins),"root");
cpull->Print(Form("%s%s_%s_pull_%s_template_%u.png",dir2.Data(),(truthcheck)? "truth": "rcone_side",date.Data(),eta,nbins),"png");
//	cpull->Print(Form("%srcone_sidebtemp_%s_%s_pull_%s_%u.root",dir2.Data(),(truthcheck)? "nonrew_withtruthinfo": "rew_andnonrew",date.Data(),eta,nbins),"root");
//	cpull->Print(Form("%srcone_sidebtemp_%s_%s_pull_%s_%u.png",dir2.Data(),(truthcheck)? "nonrew_withtruthinfo": "rew_andnonrew",date.Data(),eta,nbins),"png");
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
	fcp=new TFile("../plots/April8/grEB.root","READ");
	assert(fcp);
	grcp=(TGraphAsymmErrors*)fcp->Get("grEB");
}
else 	
{
	fcp=new TFile("../plots/April8/grnEB.root","READ");
	assert(fcp);
	grcp=(TGraphAsymmErrors*)fcp->Get("grnEB");
}
assert(grcp);
// ************************read in mass bins and purity
//truth templates:
const char *path=Form("purity_2015-04-08_truth_sumw2erron_%s_massbin_20_range",eta);
dir=Form("../plots/April8/%s/massbinned20_4bins_truthreweight/",eta);
/// if template unreweighted
//const char *path=Form("purity_2015-04-29_rcone_sideb_sumw2erron_%s_massbin_20_range",eta);
//dir=Form("../plots/April29/%s/massbinned20_4bins_nonweighted/",eta);
//reweighted templates
const char *path2=Form("purity_2015-04-29_rcone_sideb_sumw2erron_%s_massbin_20_range",eta);
dir2=Form("../plots/April29/%s/massbinned20_4bins/",eta);

TGraphErrors* gr_n1=new TGraphErrors(); TGraphErrors* gr_n2=new TGraphErrors();
TGraphErrors* gr_n3=new TGraphErrors(); TGraphErrors* gr_n4=new TGraphErrors();
int	nbins1=0; int	nbins2=0;
int	nbins3=0; int	nbins4=0;
TGraphErrors* gr_tn1=new TGraphErrors(); TGraphErrors* gr_tn2=new TGraphErrors();
TGraphErrors* gr_tn3=new TGraphErrors(); TGraphErrors* gr_tn4=new TGraphErrors();
int	nbinst1=0; int	nbinst2=0;
int	nbinst3=0; int	nbinst4=0;
//*********************************truthfit***********************************//

//TODO make more flexible for 100 bins 
//f1=new TFile(Form("%s%s_0_10.root",dir,path),"READ");
f1=new TFile(Form("%s%s_0_5.root",dir.Data(),path),"READ");
assert(f1);
//TGraphAsymmErrors* gratio=(TGraphAsymmErrors*)f1->Get("Graph");
TGraphErrors* gr1=(TGraphErrors*)f1->Get("Graph");
assert(gr1); 

f2=new TFile(Form("%s%s_5_10.root",dir.Data(),path),"READ");
//f2=new TFile(Form("%s%s_10_20.root",dir,path),"READ");
assert(f2);
//TGraphAsymmErrors* gr2=(TGraphAsymmErrors*)f2->Get("Graph");
TGraphErrors* gr2=(TGraphErrors*)f2->Get("Graph");
assert(gr2);

f3=new TFile(Form("%s%s_10_15.root",dir.Data(),path),"READ");
assert(f3);
//TGraphAsymmErrors* gr3=(TGraphAsymmErrors*)f3->Get("Graph");
TGraphErrors* gr3=(TGraphErrors*)f3->Get("Graph");
assert(gr3);
f4=new TFile(Form("%s%s_15_20.root",dir.Data(),path),"READ");
assert(f4);
//TGraphAsymmErrors* gr4=(TGraphAsymmErrors*)f4->Get("Graph");
TGraphErrors* gr4=(TGraphErrors*)f4->Get("Graph");
assert(gr4);

gr_tn1=gr1;
gr_tn2=gr2;
gr_tn3=gr3;
gr_tn4=gr4;
nbinst1=(gr1->GetN());
nbinst2=(gr2->GetN());
nbinst3=(gr3->GetN());
nbinst4=(gr4->GetN());
cout << "nbins1 truth " << nbins1<< "nbins2 truth" << nbins2<<endl;

//*************************rcone and sideband *********************************//
f11=new TFile(Form("%s%s_0_5.root",dir2.Data(),path2),"READ");
//f11=new TFile(Form("%s%s_0_10.root",dir2,path2),"READ");
assert(f11);
TGraphErrors* gr11=(TGraphErrors*)f11->Get("Graph");
assert(gr11);
f22=new TFile(Form("%s%s_5_10.root",dir2.Data(),path2),"READ");
//f22=new TFile(Form("%s%s_0_10.root",dir2,path2),"READ");
assert(f22);
TGraphErrors* gr22=(TGraphErrors*)f22->Get("Graph");
assert(gr22);

f33=new TFile(Form("%s%s_10_15.root",dir2.Data(),path2),"READ");
assert(f33);
TGraphErrors* gr33=(TGraphErrors*)f33->Get("Graph");
assert(gr33);
f44=new TFile(Form("%s%s_15_20.root",dir2.Data(),path2),"READ");
assert(f44);
TGraphErrors* gr44=(TGraphErrors*)f44->Get("Graph");
assert(gr44);

gr_n1=gr11;
gr_n2=gr22;
gr_n3=gr33;
gr_n4=gr44;
nbins1=(gr11->GetN());
nbins2=(gr22->GetN());
nbins3=(gr33->GetN());
nbins4=(gr44->GetN());
cout << "nbins1 rcone sideband " << nbins1<< "nbins2" << nbins2<<endl;

nbins=nbins1+nbins2+nbins3+nbins4;
tnbins=nbinst1+nbinst2+nbinst3+nbinst4;
if(nbins!=nbins) {return -1;}
cout << "nbins " << nbins << endl;
double *pdpmasst_all=new double[nbins];
double *masstrue=new double[nbins];
double *ptr=new double[nbins];
double *ratio=new double[nbins];
double *ppurityt_all=new double[nbins];
double *pw2errt_all=new double[nbins];
double *pdpmass_all=new double[nbins];
double *ppurity_all=new double[nbins];
double *pw2err_all=new double[nbins];
for(int k=0;k<nbins1; k++ )
{
	gr_n1->GetPoint(k,pdpmass_all[k],ppurity_all[k]);
    pw2err_all[k]=gr_n1->GetErrorY(k);
	cout << "k " << k << "dpmass " << pdpmass_all[k] << "purity " << ppurity_all[k] << "w2err_all[k] " << pw2err_all[k] <<endl;
}
//TODO get low and high error
cout << "second filling " << endl;
for(int k=0;k<nbins2; k++ )
{
	gr_n2->GetPoint(k,pdpmass_all[k+nbins1],ppurity_all[k+nbins1]);
	pw2err_all[k+nbins1]=gr_n2->GetErrorY(k);
	cout << "k " << k+nbins << "dpmass " << pdpmass_all[k+nbins1] << "purity " << ppurity_all[k+nbins1] << " w2err_all[k] " << pw2err_all[k+nbins1] <<endl;
} 

cout << "third filling " << endl;
for(int k=0;k<nbins3; k++ )
{
	gr_n3->GetPoint(k,pdpmass_all[k+nbins1+nbins2],ppurity_all[k+nbins1+nbins2]);
	pw2err_all[k+nbins1+nbins2]=gr_n3->GetErrorY(k);
	cout << "k " << k+nbins1+nbins2 << "dpmass " << pdpmass_all[k+nbins1+nbins2] << "purity " << ppurity_all[k+nbins1+nbins2] << "w2err_all[k] " << pw2err_all[k+nbins1+nbins2] <<endl;
} 
cout << "fourth filling " << endl;
for(int k=0;k<nbins4; k++ )
{
	gr_n4->GetPoint(k,pdpmass_all[k+nbins1+nbins2+nbins3],ppurity_all[k+nbins1+nbins2+nbins3]);
//	pw2errH_all[k+nbins1+nbins2+nbins3]=gr_n4->GetErrorYhigh(k);
	pw2err_all[k+nbins1+nbins2+nbins3]=gr_n4->GetErrorY(k);
	cout << "k " << k+nbins1+nbins2+nbins3 << "dpmass " << pdpmass_all[k+nbins1+nbins2+nbins3] << "purity " << ppurity_all[k+nbins1+nbins2+nbins3]  << "w2err_all[k] " << pw2err_all[k+nbins1+nbins2+nbins3] <<endl;
}
//truth templates
//
for(int k=0;k<nbinst1; k++ )
{
	gr_tn1->GetPoint(k,pdpmasst_all[k],ppurityt_all[k]);
    pw2errt_all[k]=gr_tn1->GetErrorY(k);
}
for(int k=0;k<nbinst2; k++ )
{
	gr_tn2->GetPoint(k,pdpmasst_all[k+nbinst1],ppurityt_all[k+nbinst1]);
	pw2errt_all[k+nbinst1]=gr_tn2->GetErrorY(k);
} 
for(int k=0;k<nbinst3; k++ )
{
	gr_tn3->GetPoint(k,pdpmasst_all[k+nbinst1+nbinst2],ppurityt_all[k+nbinst1+nbinst2]);
	pw2errt_all[k+nbinst1+nbinst2]=gr_n3->GetErrorY(k);
	cout << "k " << k+nbinst1+nbinst2 << "dpmass " << pdpmasst_all[k+nbinst1+nbinst2] << "purity " << ppurityt_all[k+nbinst1+nbinst2] << "w2err_all[k] " << pw2errt_all[k+nbinst1+nbinst2] <<endl;
} 
cout << "fourth filling " << endl;
for(int k=0;k<nbinst4; k++ )
{
	gr_tn4->GetPoint(k,pdpmasst_all[k+nbinst1+nbinst2+nbinst3],ppurityt_all[k+nbinst1+nbinst2+nbinst3]);
	pw2errt_all[k+nbinst1+nbinst2+nbinst3]=gr_tn4->GetErrorY(k);
	cout << "k " << k+nbinst1+nbinst2+nbinst3 << "dpmass " << pdpmasst_all[k+nbinst1+nbinst2+nbinst3] << "purity " << ppurityt_all[k+nbinst1+nbinst2+nbinst3]  << "w2err_all[k] " << pw2errt_all[k+nbinst1+nbinst2+nbinst3] <<endl;
}

for(int k=0; k< nbins ;k++)
{
	grcp->GetPoint(k,masstrue[k],ptr[k]);
	cout << "ptr[k]" << ptr[k] << endl;
	if(truth)
	{
	ratio[k]=(ptr[k]-ppurityt_all[k])/ptr[k];
	}
	else if(!truth)
	{
		ratio[k]=(ptr[k]-ppurity_all[k])/ptr[k];
//
//     for template comp
//		ratio[k]=(ppurityt_all[k]-ppurity_all[k])/ppurityt_all[k];
	}
}

TGraphErrors *gratio= new TGraphErrors(nbins,masstrue,ratio);
        gratio->SetMarkerStyle(20);
        gratio->SetTitle("");
cres->cd();

TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.15);
TPad *pad1 = new TPad("pad1", "pad1", 0., 0.15, 1., 1.0);
pad1->SetBottomMargin(0); // Upper and lower plot are joined
pad1->SetGridx();         // Vertical grid
pad1->SetGridy();
pad1->Draw();
pad1->cd(); 
pad1->SetLogx();



gratio->GetXaxis()->SetTitle("Diphoton mass [GeV]");

grcp->GetYaxis()->SetTitle("purity");
grcp->GetXaxis()->SetTitle("diphoton mass [GeV]");
//grcp->GetXaxis()->SetLimits(10.,7000.);
gratio->GetXaxis()->SetLimits(50,3000);
gratio->GetYaxis()->SetRangeUser(-0.5,0.5);
grcp->GetXaxis()->SetLimits(50,3000);
grcp->GetYaxis()->SetRangeUser(0.,0.9);
grcp->GetYaxis()->SetTitleOffset(1.2);
 
grcp->SetMarkerStyle(20);grcp->SetMarkerColor(kBlack); grcp->SetLineColor(kBlack); grcp->SetMarkerSize(1.3);
gr1->SetMarkerStyle(20);gr1->SetMarkerColor(kBlue); gr1->SetLineColor(kBlue); gr1->SetMarkerSize(1.3);
gr2->SetMarkerStyle(20);gr2->SetMarkerColor(kBlue); gr2->SetLineColor(kBlue); gr2->SetMarkerSize(1.3);
gr3->SetMarkerStyle(20);gr3->SetMarkerColor(kBlue); gr3->SetLineColor(kBlue); gr3->SetMarkerSize(1.3);
gr4->SetMarkerStyle(20);gr4->SetMarkerColor(kBlue); gr4->SetLineColor(kBlue); gr4->SetMarkerSize(1.3);
gr11->SetMarkerStyle(20);gr11->SetMarkerColor(kRed); gr11->SetLineColor(kRed); gr11->SetMarkerSize(1.3);
gr22->SetMarkerStyle(20);gr22->SetMarkerColor(kRed); gr22->SetLineColor(kRed); gr22->SetMarkerSize(1.3);
gr33->SetMarkerStyle(20);gr33->SetMarkerColor(kRed); gr33->SetLineColor(kRed); gr33->SetMarkerSize(1.3);gr44->SetMarkerStyle(20);gr44->SetMarkerColor(kRed); gr44->SetLineColor(kRed); gr44->SetMarkerSize(1.3);
leg->AddEntry(grcp,Form("truth %s",eta),"p");    
leg->AddEntry(gr1,Form("truth fit rew %s",eta),"p");    
//leg->AddEntry(gr1,Form("template fit non rew %s",eta),"p");    
leg->AddEntry(gr11,Form("template fit rew %s",eta),"p");    
grcp->SetTitle(Form("mgg %s purity",eta));

grcp->Draw("AP");

gr1->Draw("P");
gr2->Draw(" P");
gr3->Draw("P");
gr4->Draw("P");

gr11->Draw("P");
gr22->Draw(" P");
gr33->Draw("P");
gr44->Draw("P");

leg->Draw();
pad1->SetTicks(0,2);
pad1->Update();

		 
 //  pad2->SetTopMargin(0);
   pad2->SetBottomMargin(0.3);
   pad2->SetTicky();
   pad2->SetGridx(); // vertical grid
   pad2->SetGridy(); // vertical grid
   pad2->SetLogx(); // vertical grid
   pad2->Draw();
   pad2->cd();
   pad2->Update();
gratio->Draw("AP");
pad2->Update();
if(truth)
{
gratio->GetYaxis()->SetTitle("(ptr-ptrrew)/ptr");
//gratio->GetYaxis()->SetTitle("(ptr-ptpnrew)/ptr");
}
else {gratio->GetYaxis()->SetTitle("(ptnrew-ptprew)/ptnrew");}
//else {gratio->GetYaxis()->SetTitle("(ptr-ptprew)/ptr");}
//gratio->GetYaxis()->SetTitle("(ptpr-ptpnr)/ptr");
   gratio->GetYaxis()->SetNdivisions(505);
   gratio->GetYaxis()->SetTitleSize(14);
   gratio->GetYaxis()->SetTitleFont(43);
   gratio->GetYaxis()->SetTitleOffset(1.05);
   gratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   gratio->GetYaxis()->SetLabelSize(15);
 // X axis ratio plot settings
   gratio->GetXaxis()->SetTitleSize(18);
   gratio->GetXaxis()->SetTitleFont(43);
   gratio->GetXaxis()->SetTitleOffset(4.);
   gratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   gratio->GetXaxis()->SetLabelSize(15);
   //   TLine *r = new TLine(0.,ratio,9.,ratio);
   pad2->Update();
  //cres->Update();
pad2->Update(); 


//cres->Print(Form("%spurity_%s_%s_template_comp.root",dir2.Data(),date.Data(),eta),"root");
//cres->Print(Form("%spurity_%s_%s_template_comp.png",dir2.Data(),date.Data(),eta),"png");
cres->Print(Form("%spurity_%s_%s_template_truthcomp_ratio%s.root",dir2.Data(),date.Data(),eta,(truth)? "truthtemp": "rconesidebtemp"),"root");
cres->Print(Form("%spurity_%s_%s_template_truthcomp_ratio%s.png",dir2.Data(),date.Data(),eta,(truth)? "truthtemp": "rconesidebtemp"),"png");
//cres->Print(Form("%spurity_%s_%s_template_truthcomp_ratio%s.root",dir2.Data(),date.Data(),eta,(truth)? "templatenonrew": "temptotempnonrew"),"root");
//cres->Print(Form("%spurity_%s_%s_template_truthcomp_ratio%s.png",dir2.Data(),date.Data(),eta,(truth)? "templatenonrew": "temptotempnonrew"),"png");
if(truth)
{
plot_pull(grcp,pdpmasst_all,ppurityt_all,pw2errt_all,truth);

}
else if(!truth)
{
plot_pull(grcp,pdpmass_all,ppurity_all,pw2err_all,truth);
}
delete[] pdpmass_all;
delete[] ppurity_all;
delete[] pw2err_all;

return 0;
}

