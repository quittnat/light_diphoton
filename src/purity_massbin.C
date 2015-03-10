//
//M Quittnat 
//For plotting of purity vs mass bins
//
#include "TROOT.h"
#include <iostream>
#include "TFile.h"
#include <TStyle.h>
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include  "TLegend.h"
#include "TAxis.h"
#include "TString.h"
#include "TH2F.h"
#include <Riostream.h>
#include "TTree.h"
#include <stdlib.h>
#include <string.h>
using namespace std;
//massbins
const int nbins=20;
//main script
int plot()
{
gROOT->Reset();
gROOT->SetStyle("Plain"); // possibilities: Default, Plain, Bold, Video, Pub
gStyle->SetOptStat(0);
gStyle->SetOptFit(0);
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
//

TFile* fcp;TFile* f1;TFile* f2;TFile* f3;TFile* f4;
fcp=new TFile("../plots/March09/150309_purity_counting.root","READ");
assert(fcp);
TCanvas *ccp = (TCanvas*)fcp->Get("c1");
assert(ccp);
TPad* padcp=(TPad*)ccp->GetPrimitive("pad1");
//TGraphErrors* grcp=(TGraphErrors*)padcp->GetPrimitive("grEB");
TGraphErrors* grcp=(TGraphErrors*)padcp->GetPrimitive("grnEB");
assert(grcp);
// TGraphErrors* grEB=(TGraphErrors*)ccp->GetPrimitive("grEB");
// TGraphErrors* grnEB=(TGraphErrors*)ccp->GetPrimitive("grnEB");
///////////// 
const char *path=Form("150306_truth_purity_nEBEB_massbin_20_range");
//read in mass bins and purity
//truthfit
f1=new TFile(Form("../plots/March05/truthfit/nEBEB/%s_0_11.root",path),"READ");
assert(f1);
TCanvas *c1 = (TCanvas*)f1->Get("cgr");
assert(c1);
TGraphErrors* gr1=(TGraphErrors*)c1->GetPrimitive("Graph");
assert(gr1); 
const int	nbins1=gr1->GetN();
cout << "nbins1 " << nbins1<<endl;
//for (i=0;i< nbins1;i++){
//mvec[i]=gr1->GetXaxis()->GetBinUpEdge(i);
//cout << "mvec[i]" << mvec[i] << endl;
// }
f2=new TFile(Form("../plots/March05/truthfit/nEBEB/%s_10_21.root",path),"READ");
assert(f2);
TCanvas *c2 = (TCanvas*)f2->Get("cgr");
assert(c2);
TGraphErrors* gr2=(TGraphErrors*)c2->GetPrimitive("Graph");
assert(gr2);

const int	nbins2=gr2->GetN();
cout << "nbins2 " << nbins2<<endl;
//rcone
//const char *path2=Form("150306_truth_purity_nEBEB_massbin_20_range");
const char *path2=Form("150306_rcone_sideb_purity_nEBEB_massbin_20_range");
f3=new TFile(Form("../plots/March05/templates/nEBEB/%s_0_11.root",path2),"READ");
//f3=new TFile("../plots/March05/templates/nEBEB/150306_rcone_sideb_purity_nEBEB_massbin_20_range_0_11.root","READ");
//f3=new TFile(Form("../plots/March05/truthfit/nEBEB/%s_0_11.root",path2),"READ");
assert(f3);
TCanvas *c3 = (TCanvas*)f3->Get("cgr");
assert(c3);
TGraphErrors* gr3=(TGraphErrors*)c3->GetPrimitive("Graph");
assert(gr3);

//f4=new TFile("../plots/March05/templates/nEBEB/150306_rcone_sideb_purity_nEBEB_massbin_20_range_10_21.root","READ");
f4=new TFile(Form("../plots/March05/templates/nEBEB/%s_10_21.root",path2),"READ");
//f4=new TFile(Form("../plots/March05/truthfit/EBEB/%s_10_21.root",path2),"READ");
assert(f4);
TCanvas *c4 = (TCanvas*)f4->Get("cgr");
assert(c4);
TGraphErrors* gr4=(TGraphErrors*)c4->GetPrimitive("Graph");
assert(gr4);
//plot
cres->cd();
TPad *pad1 = new TPad("pad1","",0,0,1,1);
pad1->Draw();
pad1->cd(); 
pad1->SetGridy();
pad1->SetGridy();
pad1->SetLogx();
//gr1->SetTitle("mgg nEBEB");
//gr1->SetTitle("mgg nEBEB purity, truth & template fit");
gr1->GetYaxis()->SetTitle("purity");
gr1->GetXaxis()->SetTitle("diphoton mass [GeV]");
//gr1->GetXaxis()->SetRangeUser(50.,3000.);
gr1->GetXaxis()->SetLimits(50.,3000.);
//gr1->GetXaxis()->SetRange(50,3000);
gr1->GetYaxis()->SetRangeUser(.1,1.);
gr1->GetYaxis()->SetTitleOffset(1.2);
 //gr2->GetXaxis()->SetNdivisions(504);
 // gr2->GetYaxis()->SetNdivisions(503);
 
cout << __LINE__ << endl;
 grcp->SetMarkerStyle(20);grcp->SetMarkerColor(kBlack); grcp->SetLineColor(kBlack); grcp->SetMarkerSize(1.3);
cout << __LINE__ << endl;
 gr1->SetMarkerStyle(20);gr1->SetMarkerColor(kRed); gr1->SetLineColor(kRed); gr1->SetMarkerSize(1.3);
cout << __LINE__ << endl;
 gr2->SetMarkerStyle(20);gr2->SetMarkerColor(kRed); gr2->SetLineColor(kRed); gr2->SetMarkerSize(1.3);
 cout << __LINE__ << endl;
 gr3->SetMarkerStyle(20);gr3->SetMarkerColor(kBlue); gr3->SetLineColor(kBlue); gr3->SetMarkerSize(1.3);
 cout << __LINE__ << endl;
 gr4->SetMarkerStyle(20);gr4->SetMarkerColor(kBlue); gr4->SetLineColor(kBlue); gr4->SetMarkerSize(1.3);
cout << __LINE__ << endl;
leg->AddEntry(grcp,"truth  nEBEB","p");    
leg->AddEntry(gr2,"truth fit nEBEB","p");    
leg->AddEntry(gr3,"template fit nEBEB","p");    
cout << __LINE__ << endl;
gr1->SetTitle("mgg nEBEB purity");
gr1->Draw("AP");
grcp->Draw("P");
gr2->Draw(" P");
gr3->Draw("P");
gr4->Draw("P");
leg->Draw();
pad1->SetTicks(0,2);
//cres->Print("../plots/March10/150310_purity_truth_etacomp.root","root");
//cres->Print("../plots/March10/150310_purity_truth_etacomp.png","png");
cres->Print("../plots/March10/150310_purity_nEBEB_truth_rconecomp.root","root");
cres->Print("../plots/March10/150310_purity_nEBEB_truth_rconecomp.png","png");
// cres->SaveAs(Form("../plots/March05/%s.root",path));
// cres->SaveAs(Form("../plots/March05/%s.png",path));

return 0;
}
