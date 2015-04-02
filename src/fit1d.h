
#include <assert.h>
#include "TLine.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TStyle.h"
#include "RooFitResult.h"
#include "TString.h"
#include "TFile.h"
#include "TLatex.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom3.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooHistPdf.h"
#include "RooFormulaVar.h"
#include "RooRealVar.h"
#include "RooRealConstant.h"
#include "RooPlot.h"
#include "RooMinuit.h"
#include "RooMCStudy.h"
#include "RooBinning.h"
#include "RooGaussian.h"
#include "RooExtendPdf.h"
#include "RooGenericPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include <stdio.h>
#include "RooNLLVar.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooMinimizer.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooKeysPdf.h"
#include "RooNDKeysPdf.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooClassFactory.h"
#include "TCanvas.h"
#include "RooConstraintSum.h"
#include "RooAddition.h"
#include "RooAbsDataStore.h"
#include "RooCachedPdf.h"
#include "RooThresholdCategory.h"
#include "TF1.h"
#include "TF2.h"
#include "TLegend.h"
#include "TSystem.h"

Bool_t dovzero=kFALSE;
Bool_t dohggv=kTRUE;
Bool_t dowv=kFALSE;

//declaration of functions and variables
RooRealVar *roovar1=NULL;
RooRealVar *roovar2=NULL;
RooRealVar *rooisopv1=NULL;
RooRealVar *rooisopv2=NULL;
RooRealVar *rooisowv1=NULL;
RooRealVar *rooisowv2=NULL;
RooRealVar *roopt1=NULL;
RooRealVar *roopt2=NULL;
RooRealVar *roodiphopt=NULL;
RooRealVar *roodiphomass=NULL;
RooRealVar *roosieie1=NULL;
RooRealVar *roosieie2=NULL;
RooRealVar *rooeta1=NULL;      
RooRealVar *rooeta2=NULL;
RooRealVar *roorho=NULL;
RooRealVar *roosigma=NULL;
RooRealVar *roonvtx=NULL;
RooRealVar *rooweight=NULL;
RooDataSet *dataset_data = NULL;
RooDataSet *dataset_rcone = NULL;
RooDataSet *dataset_sideband = NULL;
RooDataSet *dataset_testsig = NULL;
RooDataSet *dataset_test1 = NULL;
RooDataSet *dataset_testbkg = NULL;
RooDataSet *dataset_test2 = NULL;
RooDataSet *dset_data = NULL;
RooDataSet *dset_rcone = NULL;
RooDataSet *dset_testsig1 = NULL;
RooDataSet *dset_sideb = NULL;
RooDataSet *dset_testbkg1 = NULL;
RooRealVar *binning_roovar1=NULL;
TString eta_q;
TString dir;
Int_t startbin_q=-100;
Int_t  endbin_q=-100;
Int_t ntot_q=-100;
Int_t nq_q=-100;
bool truthfit;
bool massbinned=kTRUE;
Bool_t isdata_q=kTRUE;
Bool_t EBEB;
bool debug=kFALSE;
bool check=kFALSE;

//bins for reweihgting
RooBinning tbins (0,9); 
const int n_ptbins_forreweighting = 5;
Float_t ptbins_forreweighting[n_ptbins_forreweighting+1]={20,40,60,100,150,999};
const int n_diphoptbins_forreweighting = 4;
Float_t diphoptbins_forreweighting[n_diphoptbins_forreweighting+1]={0,40,50,200,2000};
//const int n_ptbins_forreweighting =6;
//Float_t ptbins_forreweighting[n_ptbins_forreweighting+1]={20,35,50,80,100,200,999};
//const int n_diphoptbins_forreweighting = 6;
//Float_t diphoptbins_forreweighting[n_diphoptbins_forreweighting+1]={0,20,40,50,60,100,2000};
const int n_etabins_forreweighting = 13;
Float_t etabins_forreweighting[n_etabins_forreweighting+1]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.3,1.9,2.5};
const int n_rhobins_forreweighting = 4;
Float_t rhobins_forreweighting[n_rhobins_forreweighting+1]={0,10,15,20,60};
const int n_sigmabins_forreweighting = 5;
Float_t sigmabins_forreweighting[n_sigmabins_forreweighting+1]={0,2,3,4,6,20};


//histograms for reweightings
double wxmin=0.;
double wxmax=10.;
int wbins=100;
int wbins2=200;
double wxmin2=0.;
double wxmax2=20.;

TH2F*  w_rhornr=new TH2F("w_rhornr","oldw/neww",n_rhobins_forreweighting ,0.,n_rhobins_forreweighting+1,wbins2,wxmin2,20.);
TH2F*  w_rhor2o=new TH2F("w_rhor2o","oldw",n_rhobins_forreweighting,0.,n_rhobins_forreweighting+1,wbins2,wxmin2,25.);
TH2F*  w_rhorn=new TH2F("w_rhorn","neww",n_rhobins_forreweighting ,0.,n_rhobins_forreweighting+1,wbins,wxmin,10.);
TH2F*  w_rhor2=new TH2F("w_rhor2","oldw/neww",50,0.,50.,wbins2,wxmin2,15);
TH1F *w_rhor= new TH1F("w_rhor","neww",wbins ,wxmin,20.);
TH1F *w_rhoro= new TH1F("w_rhoro","oldw",wbins ,wxmin,25.);
TH2F*  w_sigmarnr=new TH2F("w_sigmarnr","w_sigmarnr",n_sigmabins_forreweighting,0.,n_sigmabins_forreweighting+1,wbins2,wxmin2,wxmax2);
TH2F*  w_sigmar2o=new TH2F("w_sigmar2o","w_sigmar2o",n_sigmabins_forreweighting ,0.,n_sigmabins_forreweighting+1,wbins2,wxmin2,25.);
TH2F*  w_sigmarn=new TH2F("w_sigmarn","w_sigmarn",n_sigmabins_forreweighting ,0.,n_sigmabins_forreweighting+1,wbins,wxmin,wxmax2);
TH2F*  w_sigmar2=new TH2F("w_sigmar2","w_sigmar2",50,0.,50.,wbins2,wxmin2,15);


TH2F*  w_etanr=new TH2F("w_etanr","old/new weights",n_etabins_forreweighting ,0.,n_etabins_forreweighting+1,wbins2,wxmin2,3.);
TH2F*  w_eta2o=new TH2F("w_eta2o","old weights",n_etabins_forreweighting ,0.,n_etabins_forreweighting+1,wbins2,wxmin2,10.);
TH2F*  w_etan=new TH2F("w_etan","new weights",n_etabins_forreweighting ,0.,n_etabins_forreweighting+1,wbins,wxmin,15.);
TH2F*  w_eta2=new TH2F("w_eta2","old/new weights",n_etabins_forreweighting,0.,2.5,wbins2,wxmin2,2.);
TH1F *w_eta= new TH1F("w_eta","new weights",wbins ,wxmin,15.);
TH1F *w_etao= new TH1F("w_etao","old weights",wbins ,wxmin,10.);
TH2F*  w_etarnr=new TH2F("w_etarnr","old/neww",n_etabins_forreweighting ,0.,n_etabins_forreweighting+1,wbins2,wxmin2,3.);
TH2F*  w_etar2o=new TH2F("w_etar2o","oldw",n_etabins_forreweighting ,0.,n_etabins_forreweighting+1,wbins2,wxmin2,10.);
TH2F*  w_etarn=new TH2F("w_etarn","neww",n_etabins_forreweighting ,0.,n_etabins_forreweighting+1,wbins,wxmin,15.);
TH2F*  w_etar2=new TH2F("w_etar2","oldw/neww",n_etabins_forreweighting,0.,2.5,wbins2,wxmin2,2.);
TH1F *w_etar= new TH1F("w_etar","neww",wbins ,wxmin,15.);
TH1F *w_etaro= new TH1F("w_etaro","oldw",wbins ,wxmin,10.);

TH2F*  w_rhonr=new TH2F("w_rhonr","old/new w",n_rhobins_forreweighting ,0.,n_rhobins_forreweighting+1,20,wxmin2,2.);
TH2F*  w_rho2o=new TH2F("w_rho2o","old w",n_rhobins_forreweighting ,0.,n_rhobins_forreweighting+1,10,wxmin2,10.);
TH2F*  w_rhon=new TH2F("w_rhon","new w",n_rhobins_forreweighting ,0.,n_rhobins_forreweighting+1,10,wxmin,10.);
TH2F*  w_rho2=new TH2F("w_rho2","oldw/neww",n_rhobins_forreweighting,0.,n_rhobins_forreweighting+1,20,wxmin2,2.);
TH1F *w_rho= new TH1F("w_rho","new weights",wbins ,wxmin,10.);
TH1F *w_rhoo= new TH1F("w_rhoo","oldw",wbins ,wxmin,10.);
TH2F*  w_sigmanr=new TH2F("w_sigmanr","oldw/neww",n_sigmabins_forreweighting ,0.,n_sigmabins_forreweighting+1,20,wxmin2,2.);
TH2F*  w_sigma2o=new TH2F("w_sigma2o","oldw",n_sigmabins_forreweighting ,0.,n_sigmabins_forreweighting+1,10,wxmin2,10.);
TH2F*  w_sigman=new TH2F("w_sigman","new weights",n_sigmabins_forreweighting ,0.,n_sigmabins_forreweighting+1,10,wxmin,10.);
TH2F*  w_sigma2=new TH2F("w_sigma2","oldw/neww",n_sigmabins_forreweighting,0.,n_sigmabins_forreweighting+1,wbins2,wxmin2,2.);
TH2F*  w_diphoptnr=new TH2F("w_diphoptnr","oldw/neww",n_diphoptbins_forreweighting ,0.,n_diphoptbins_forreweighting+1,wbins2,wxmin2,20.);
TH2F*  w_diphopt2o=new TH2F("w_diphopt2o","oldw",n_diphoptbins_forreweighting ,0.,n_diphoptbins_forreweighting+1,wbins2,wxmin2,25.);
TH2F*  w_diphoptn=new TH2F("w_diphoptn","new weights vs pt bins",n_diphoptbins_forreweighting ,0.,n_diphoptbins_forreweighting+1,wbins,wxmin,100.);
TH2F*  w_diphopt2=new TH2F("w_diphopt2","oldw/neww",100,0.,1000,wbins2,wxmin2,15);
TH2F*  w_ptn=new TH2F("w_ptn","new weights vs pt bins",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1,wbins,wxmin,25.);
TH2F*  w_ptnr=new TH2F("w_ptnr","oldw/neww",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1, wbins2,wxmin2,10.);
TH2F*  w_pt2o=new TH2F("w_pt2o","oldw",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1, wbins2,wxmin2,15.);
TH2F*  w_pt2=new TH2F("w_pt2","oldw/neww",500 ,0.,500,wbins2,wxmin2,5.);
TH1F *w_pt= new TH1F("w_pt","neww",wbins ,wxmin,wxmax);
TH1F *w_diphopt= new TH1F("w_diphopt","new weights",wbins ,wxmin,wxmax);
TH1F *w_diphopto= new TH1F("w_diphopto","oldw",wbins ,wxmin,wxmax);
TH1F *w_pto= new TH1F("w_pto","oldw",wbins ,wxmin,wxmax);
TH2F*  w_diphoptrnr=new TH2F("w_diphoptrnr","w_diphoptrnr",n_diphoptbins_forreweighting ,0.,n_diphoptbins_forreweighting+1,wbins2,wxmin2,20.);
TH2F*  w_diphoptr2o=new TH2F("w_diphoptr2o","w_diphoptr2o",n_diphoptbins_forreweighting ,0.,n_diphoptbins_forreweighting+1,wbins2,wxmin2,25.);
TH2F*  w_diphoptrn=new TH2F("w_diphoptrn","w_diphoptrn",n_diphoptbins_forreweighting ,0.,n_diphoptbins_forreweighting+1,wbins,wxmin,100.);
TH2F*  w_diphoptr2=new TH2F("w_diphoptr2","w_diphoptr2",100,0.,1000,wbins2,wxmin2,15);
TH2F*  w_ptrn=new TH2F("w_ptrn","w_ptrn",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1,wbins,wxmin,25.);
TH2F*  w_ptrnr=new TH2F("w_ptrnr","w_ptrnr",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1, wbins2,wxmin2,10.);
TH2F*  w_ptr2o=new TH2F("w_ptr2o","w_ptr2o",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1, wbins2,wxmin2,15.);
TH2F*  w_ptr2=new TH2F("w_ptr2","w_ptr2",500 ,0.,500,wbins2,wxmin2,wxmax);
TH1F *w_ptr= new TH1F("w_ptr","w_ptr",wbins ,wxmin,wxmax);
TH1F *w_diphoptr= new TH1F("w_diphoptr","w_diphoptr",wbins ,wxmin,wxmax);
TH1F *w_diphoptro= new TH1F("w_diphoptro","w_diphoptro",wbins ,wxmin,wxmax);
TH1F *w_ptro= new TH1F("w_ptro","w_ptro",wbins ,wxmin,12.);


TH2F*  tw_rhornr=new TH2F("tw_rhornr","oldw/neww",n_rhobins_forreweighting ,0.,n_rhobins_forreweighting+1,wbins2,wxmin2,20.);
TH2F*  tw_rhor2o=new TH2F("tw_rhor2o","oldw",n_rhobins_forreweighting,0.,n_rhobins_forreweighting+1,wbins2,wxmin2,25.);
TH2F*  tw_rhorn=new TH2F("tw_rhorn","neww",n_rhobins_forreweighting ,0.,n_rhobins_forreweighting+1,wbins,wxmin,10.);
TH2F*  tw_rhor2=new TH2F("tw_rhor2","oldw/neww",50,0.,50.,wbins2,wxmin2,15);
TH1F *tw_rhor= new TH1F("tw_rhor","neww",wbins ,wxmin,20.);
TH1F *tw_rhoro= new TH1F("tw_rhoro","oldw",wbins ,wxmin,25.);
TH2F*  tw_sigmarnr=new TH2F("tw_sigmarnr","tw_sigmarnr",n_sigmabins_forreweighting,0.,n_sigmabins_forreweighting+1,wbins2,wxmin2,wxmax2);
TH2F*  tw_sigmar2o=new TH2F("tw_sigmar2o","tw_sigmar2o",n_sigmabins_forreweighting ,0.,n_sigmabins_forreweighting+1,wbins2,wxmin2,25.);
TH2F*  tw_sigmarn=new TH2F("tw_sigmarn","tw_sigmarn",n_sigmabins_forreweighting ,0.,n_sigmabins_forreweighting+1,wbins,wxmin,wxmax2);
TH2F*  tw_sigmar2=new TH2F("tw_sigmar2","tw_sigmar2",50,0.,50.,wbins2,wxmin2,15);


TH2F*  tw_etanr=new TH2F("tw_etanr","old/new weights",n_etabins_forreweighting ,0.,n_etabins_forreweighting+1,wbins2,wxmin2,3.);
TH2F*  tw_eta2o=new TH2F("tw_eta2o","old weights",n_etabins_forreweighting ,0.,n_etabins_forreweighting+1,wbins2,wxmin2,10.);
TH2F*  tw_etan=new TH2F("tw_etan","new weights",n_etabins_forreweighting ,0.,n_etabins_forreweighting+1,wbins,wxmin,15.);
TH2F*  tw_eta2=new TH2F("tw_eta2","old/new weights",n_etabins_forreweighting,0.,2.5,wbins2,wxmin2,2.);
TH1F *tw_eta= new TH1F("tw_eta","new weights",wbins ,wxmin,15.);
TH1F *tw_etao= new TH1F("tw_etao","old weights",wbins ,wxmin,10.);
TH2F*  tw_etarnr=new TH2F("tw_etarnr","old/neww",n_etabins_forreweighting ,0.,n_etabins_forreweighting+1,wbins2,wxmin2,3.);
TH2F*  tw_etar2o=new TH2F("tw_etar2o","oldw",n_etabins_forreweighting ,0.,n_etabins_forreweighting+1,wbins2,wxmin2,10.);
TH2F*  tw_etarn=new TH2F("tw_etarn","neww",n_etabins_forreweighting ,0.,n_etabins_forreweighting+1,wbins,wxmin,15.);
TH2F*  tw_etar2=new TH2F("tw_etar2","oldw/neww",n_etabins_forreweighting,0.,2.5,wbins2,wxmin2,2.);
TH1F *tw_etar= new TH1F("tw_etar","neww",wbins ,wxmin,15.);
TH1F *tw_etaro= new TH1F("tw_etaro","oldw",wbins ,wxmin,10.);

TH2F*  tw_rhonr=new TH2F("tw_rhonr","old/new w",n_rhobins_forreweighting ,0.,n_rhobins_forreweighting+1,20,wxmin2,2.);
TH2F*  tw_rho2o=new TH2F("tw_rho2o","old w",n_rhobins_forreweighting ,0.,n_rhobins_forreweighting+1,10,wxmin2,10.);
TH2F*  tw_rhon=new TH2F("tw_rhon","new w",n_rhobins_forreweighting ,0.,n_rhobins_forreweighting+1,10,wxmin,10.);
TH2F*  tw_rho2=new TH2F("tw_rho2","oldw/neww",n_rhobins_forreweighting,0.,n_rhobins_forreweighting+1,20,wxmin2,2.);
TH1F *tw_rho= new TH1F("tw_rho","new weights",wbins ,wxmin,10.);
TH1F *tw_rhoo= new TH1F("tw_rhoo","oldw",wbins ,wxmin,10.);
TH2F*  tw_sigmanr=new TH2F("tw_sigmanr","oldw/neww",n_sigmabins_forreweighting ,0.,n_sigmabins_forreweighting+1,20,wxmin2,2.);
TH2F*  tw_sigma2o=new TH2F("tw_sigma2o","oldw",n_sigmabins_forreweighting ,0.,n_sigmabins_forreweighting+1,10,wxmin2,10.);
TH2F*  tw_sigman=new TH2F("tw_sigman","new weights",n_sigmabins_forreweighting ,0.,n_sigmabins_forreweighting+1,10,wxmin,10.);
TH2F*  tw_sigma2=new TH2F("tw_sigma2","oldw/neww",n_sigmabins_forreweighting,0.,n_sigmabins_forreweighting+1,wbins2,wxmin2,2.);
TH2F*  tw_diphoptnr=new TH2F("tw_diphoptnr","oldw/neww",n_diphoptbins_forreweighting ,0.,n_diphoptbins_forreweighting+1,wbins2,wxmin2,20.);
TH2F*  tw_diphopt2o=new TH2F("tw_diphopt2o","oldw",n_diphoptbins_forreweighting ,0.,n_diphoptbins_forreweighting+1,wbins2,wxmin2,25.);
TH2F*  tw_diphoptn=new TH2F("tw_diphoptn","new weights vs pt bins",n_diphoptbins_forreweighting ,0.,n_diphoptbins_forreweighting+1,wbins,wxmin,100.);
TH2F*  tw_diphopt2=new TH2F("tw_diphopt2","oldw/neww",100,0.,1000,wbins2,wxmin2,15);
TH2F*  tw_ptn=new TH2F("tw_ptn","new weights vs pt bins",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1,wbins,wxmin,25.);
TH2F*  tw_ptnr=new TH2F("tw_ptnr","oldw/neww",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1, wbins2,wxmin2,10.);
TH2F*  tw_pt2o=new TH2F("tw_pt2o","oldw",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1, wbins2,wxmin2,15.);
TH2F*  tw_pt2=new TH2F("tw_pt2","oldw/neww",500 ,0.,500,wbins2,wxmin2,5.);
TH1F *tw_pt= new TH1F("tw_pt","neww",wbins ,wxmin,wxmax);
TH1F *tw_diphopt= new TH1F("tw_diphopt","new weights",wbins ,wxmin,wxmax);
TH1F *tw_diphopto= new TH1F("tw_diphopto","oldw",wbins ,wxmin,wxmax);
TH1F *tw_pto= new TH1F("tw_pto","oldw",wbins ,wxmin,wxmax);
TH2F*  tw_diphoptrnr=new TH2F("tw_diphoptrnr","tw_diphoptrnr",n_diphoptbins_forreweighting ,0.,n_diphoptbins_forreweighting+1,wbins2,wxmin2,20.);
TH2F*  tw_diphoptr2o=new TH2F("tw_diphoptr2o","tw_diphoptr2o",n_diphoptbins_forreweighting ,0.,n_diphoptbins_forreweighting+1,wbins2,wxmin2,25.);
TH2F*  tw_diphoptrn=new TH2F("tw_diphoptrn","tw_diphoptrn",n_diphoptbins_forreweighting ,0.,n_diphoptbins_forreweighting+1,wbins,wxmin,100.);
TH2F*  tw_diphoptr2=new TH2F("tw_diphoptr2","tw_diphoptr2",100,0.,1000,wbins2,wxmin2,15);
TH2F*  tw_ptrn=new TH2F("tw_ptrn","tw_ptrn",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1,wbins,wxmin,25.);
TH2F*  tw_ptrnr=new TH2F("tw_ptrnr","tw_ptrnr",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1, wbins2,wxmin2,10.);
TH2F*  tw_ptr2o=new TH2F("tw_ptr2o","tw_ptr2o",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1, wbins2,wxmin2,15.);
TH2F*  tw_ptr2=new TH2F("tw_ptr2","tw_ptr2",500 ,0.,500,wbins2,wxmin2,wxmax);
TH1F *tw_ptr= new TH1F("tw_ptr","tw_ptr",wbins ,wxmin,wxmax);
TH1F *tw_diphoptr= new TH1F("tw_diphoptr","tw_diphoptr",wbins ,wxmin,wxmax);
TH1F *tw_diphoptro= new TH1F("tw_diphoptro","tw_diphoptro",wbins ,wxmin,wxmax);
TH1F *tw_ptro= new TH1F("tw_ptro","tw_ptro",wbins ,wxmin,12.);
  
//define functions for fit
void get_roodset_from_ttree(TDirectoryFile *f, TString treename, RooDataSet* &roodset);
//int do1dfit(RooDataSet *dset_data_axis1,RooDataSet *dset_sig_axis1,RooDataSet *dset_bkg_axis1,int im);
int do1dfit(RooDataSet *dset_data_axis1,RooDataSet* dset_sigcut, RooDataSet* dset_bkgcut,RooDataSet *dset_sig_axis1,RooDataSet *dset_bkg_axis1,RooDataSet* dset_testsig,RooDataSet* dset_testbkg, int im, double pu[], double puerr[],double puerrLo[],double puerrHi[],double pullerr[],double pullerrLo[],double pullerrHi[]);
/*
void reweight_rhosigma(TH1F* weight_rho, TH1F* weight_rhoo,TH2F*weight_rhon,TH2F*weight_rho2o,TH2F* weight_rhonr, TH2F* weight_rho2,TH2F*weight_sigman,TH2F* weight_sigma2o,TH2F* weight_sigmanr, TH2F* weight_sigma2,RooDataSet **dset, RooDataSet *dsetdestination, bool deleteold = kTRUE);
void reweight_pt_1d(TH1F* weight_pt,TH1F* weight_pto,TH2F*weight_ptn,TH2F* weight_pt2o,TH2F* weight_ptnr,TH2F*weight_pt2,RooDataSet **dset, RooDataSet *dsetdestination, int numvar);
void reweight_diphotonpt_1d(TH1F* weight_diphopt,TH1F* weight_diphopto,TH2F*weight_diphoptn,TH2F* weight_diphopt2o,TH2F*weight_diphoptnr,TH2F*weight_diphopt2,RooDataSet **dset, RooDataSet *dsetdestination);
void reweight_eta_1d(TH1F* weight_eta,TH1F* weight_etao,TH2F*weight_etan,TH2F* weight_eta2o,TH2F* weight_etanr,TH2F*weight_eta2,RooDataSet **dset, RooDataSet *dsetdestination, int numvar);
*/

void reweight_rhosigma(TH1F* weight_rho=NULL, TH1F* weight_rhoo=NULL,TH2F*weight_rhon=NULL,TH2F*weight_rho2o=NULL,TH2F* weight_rhonr=NULL, TH2F* weight_rho2=NULL,TH2F*weight_sigman=NULL,TH2F* weight_sigma2o=NULL,TH2F* weight_sigmanr=NULL, TH2F* weight_sigma2=NULL,RooDataSet **dset=NULL, RooDataSet *dsetdestination=NULL, bool deleteold = kTRUE);
void reweight_pt_1d(TH1F* weight_pt=NULL,TH1F* weight_pto=NULL,TH2F*weight_ptn=NULL,TH2F* weight_pt2o=NULL,TH2F* weight_ptnr=NULL,TH2F*weight_pt2=NULL,RooDataSet **dset=NULL, RooDataSet *dsetdestination=NULL, int numvar=0);
void reweight_diphotonpt_1d(TH1F* weight_diphopt=NULL,TH1F* weight_diphopto=NULL,TH2F*weight_diphoptn=NULL,TH2F* weight_diphopt2o=NULL,TH2F*weight_diphoptnr=NULL,TH2F*weight_diphopt2=NULL,RooDataSet **dset=NULL, RooDataSet *dsetdestination=NULL);


void ratioplot(RooDataSet *data_axis1, RooDataSet *truth_axis1, Bool_t logplot);
//double quantiles(TH1F* mass, int ntot, int nqtemp,int startn, int endn, TString eta_q,double _proptemp[], double dpmqtemp[]);
int quantiles(TH1F* mass,double _proptemp[], double dpmqtemp[]);
int plot_purity_massbin(double mass[],double masserr[],double pur[],double purerr[],double puerrLo[],double puerrHi[], bool errb) ; 
int prep_dataset(Bool_t isdata=kTRUE,TString eta="EBEB", Bool_t truthf=kFALSE,int ntotbins=10,int startbin=0, int endbin=2);
 void print_mem();
const int numcpu=1;
TString cut_s;
TString date;
ProcInfo_t procinfo;
