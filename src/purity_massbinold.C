//
//M Quittnat 
//For plotting of purity vs mass bins
//
#include "TROOT.h"
#include <iostream>
#include <TStyle.h>
#include "TGraphErrors.h"
#include "TCanvas.h"
#include  "TLegend.h"
#include "TAxis.h"
#include "TH2F.h"
#include <Riostream.h>
//#include <fstream.h>
#include <stdlib.h>
//#include <iomanip.h>
#include <string.h>
using namespace std;
int nbins=10;
// Define User fit function: 2 exponentials + 1 constant

//Double_t fitf(Double_t *x, Double_t *par)
//{  
//Double_t fitval=fabs(par[0])*exp(-x[0]/par[3])+fabs(par[1])*exp(-x[0]/par[4])+fabs(par[2]);
//return fitval;
//}

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

 Double_t mvec[10] ;
 Double_t pvec[10] ;
 Double_t errvec[10] ;
 Double_t xvec[10];
 Double_t mvect[10] ;
 Double_t pvect[10] ;
 Double_t errvect[10] ;
 Double_t xvect[10];
 Double_t mvecr[10] ;
 Double_t pvecr[10] ;
 Double_t errvecr[10] ;
 Double_t xvecr[10];

 char str[100];
 double value1;
 double value2;
 double value3;
 double value1t;
 double value2t;
 double value3t;
 double value1r;
 double value2r;
 double value3r;
 int i=0;

// Define Canvas

 TCanvas *c1 = new TCanvas("c1","Graph",0,0,1024,768);
TLegend* leg = new TLegend(0.6,0.7,0.7,0.8);
 leg->SetFillColor(10);
 leg->SetLineColor(10);
 leg->SetTextSize(0.03);
 leg->SetTextFont(42);
 c1->Clear();
   c1->cd(1);

//read in mass bins and purity
 ifstream infileu("../plots/March02/EBEB/rcone_sideband/massbinnned/purity_rcone_sideband_EBEB.txt",ios::in);
 if(!infileu.is_open()){
   cout<<"No inputfile1 to open!"<<endl;
   return -1;
 } 
 for (i=0;i<nbins;i++){
   infileu.getline(str,256);
   sscanf(str,"%lf %lf %lf", &value1, &value2, &value3);
  mvec[i] = value1; pvec[i] = value2;errvec[i] = value3;
  xvec[i]=0.;
 cout << i <<" mass bins : " << mvec[i] << " purity+ error "<< pvec[i] << " " << errvec[i] << endl;
 }
 infileu.close();
/*
 ifstream infiler("../plots/March03/EBEB/rebinned/purity_rcone_sideband_EBEB_rebinned.txt",ios::in);
 if(!infiler.is_open()){
   cout<<"No inputfile2 to open!"<<endl;
   return -1;
 } 
 for (i=0;i<nbins;i++){
   infiler.getline(str,256);
   sscanf(str,"%lf %lf %lf", &value1r, &value2r, &value3r);
  mvecr[i] = value1r; pvecr[i] = value2r;errvecr[i] = value3r;
  xvecr[i]=0.;
 cout << i <<" mass bins : " << mvecr[i] << " purity+ error "<< pvecr[i] << " " << errvecr[i] << endl;
 }
*/

 ifstream infilet("../plots/March02/EBEB/truth_fit/massbinned/purity_truth_EBEB.txt",ios::in);
 if(!infilet.is_open()){
   cout<<"No inputfile3 to open!"<<endl;
   return -1;
 }
//TH2F* gr2=new TH2F("gr2","gr2",nbins,0.,3000.,nbins,0.,1.);
 for (i=0;i<nbins;i++){
   infilet.getline(str,256);
   sscanf(str,"%lf %lf %lf", &value1t, &value2t, &value3t);
  mvect[i] = value1t; pvect[i] = value2t;errvect[i] = value3t;
  xvect[i]=0.;
/*   if(i!=0){
	   mvec[i]=mvec[i]-(mvec[i]-mvec[i-1])/2;
	   xvec[i]=(mvec[i]-mvec[i-1])/2;
   }
   else if(i==0) {
	   mvec[i]=79/2.;
	   xvec[0]=79./2.;}
*/	   
//   gr2->Fill(mvec[i],pvec[i]);
 //  gr2->SetBinError(i,errvec[i]);
 cout << i <<" mass bins : " << mvect[i] << " purity+ error "<< pvect[i] << " " << errvect[i] << endl;
 }
 infilet.close();

// draw graph with crystal  data

 TGraph *gr2 = new TGraphErrors(nbins, mvec, pvec,xvec,errvec);
 TGraph *gr = new TGraphErrors(nbins, mvect, pvect,xvect,errvect);
 //TGraph *grr = new TGraphErrors(nbins, mvecr, pvecr,xvecr,errvecr);
c1->SetLogx();
c1->SetGridy();
 gr2->SetTitle("mgg EBEB");
gr2->SetMarkerStyle(20);
 gr2->SetMarkerSize(1.3);
 gr2->GetYaxis()->SetTitle("purity");
 gr2->GetXaxis()->SetTitle("diphoton mass [GeV]");
//  gr2->GetXaxis()->SetLimits(0.,3000.);
 //gr2->GetXaxis()->SetRangeUser(0.01,1000.);
gr2->GetYaxis()->SetRangeUser(.4,1.);
 gr2->GetYaxis()->SetTitleOffset(1.2);
 //gr2->GetXaxis()->SetNdivisions(504);
 // gr2->GetYaxis()->SetNdivisions(503);
 gr2->SetMarkerStyle(20);
 gr2->SetMarkerColor(kRed);
 gr2->SetLineColor(kRed);
 gr2->SetMarkerSize(1.); 
 gr->SetMarkerStyle(20);
 gr->SetMarkerColor(kBlue);
 gr->SetLineColor(kBlue);
 gr->SetMarkerSize(1.); 
/* grr->SetMarkerStyle(20);
 grr->SetMarkerColor(kBlack);
 grr->SetLineColor(kBlack);
 grr->SetMarkerSize(1.); 
 */
  leg->AddEntry(gr2,"template fit purity ","p");    
//  leg->AddEntry(grr,"template fit purity rebinned","p");    
  leg->AddEntry(gr," truth fit purity ","p");    
 gr2->Draw("APL");
 gr->Draw("SAME PL");
// grr->Draw("SAME PL");
 leg->Draw();
// gPad->Modified();

 //c1->Print("150203_truth_purity_massbin.root","root");
// c1->Print("150203_truth_purity_massbin.png","png");
 c1->Print("../plots/March03/nEBEB/150303_purity_massbin_nEBEB.root","root");
 c1->Print("../plots/March03/nEBEB/150303_purity_massbin_nEBEB.png","png");
return 0;
}
