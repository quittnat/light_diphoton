
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
//RooDataSet *dataset_testsig = NULL;
RooDataSet *dataset_test1 = NULL;
//RooDataSet *dataset_testbkg = NULL;
RooDataSet *dataset_test2 = NULL;
//RooDataSet *dset_data = NULL;
//RooDataSet *dset_rcone = NULL;
RooDataSet *dset_testsig1 = NULL;
//RooDataSet *dset_sideb = NULL;
RooDataSet *dset_testbkg1 = NULL;
RooRealVar *binning_roovar1=NULL;
TString eta_q;
TString dir;
Int_t startbin_q=-100;
Int_t  endbin_q=-100;
Int_t ntot_q=-100;
Int_t nq_q=-100;
bool truthfit=kFALSE;
bool massbinned=kTRUE;
Bool_t isdata_q=kTRUE;
Bool_t EBEB=kTRUE;
bool debug=kFALSE;
bool check=kTRUE;

//bins for reweihgting
RooBinning tbins (0,9); 
const int n_ptbins_forreweighting = 5;
Float_t ptbins_forreweighting[n_ptbins_forreweighting+1]={20,40,60,100,150,999};
const int n_diphoptbins_forreweighting = 4;
Float_t diphoptbins_forreweighting[n_diphoptbins_forreweighting+1]={0,40,50,200,7000};
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


TH2F*  wrhoratiobin=new TH2F("wrhoratiobin","old/new weights",n_rhobins_forreweighting ,0.,n_rhobins_forreweighting+1,wbins2,wxmin2,6.);
TH2F*  wrhoratio=new TH2F("wrhoratio","old/new weights",50,0.,  50.,wbins2,wxmin2,6.);
TH2F*  wrho_eta=new TH2F("wrho_eta","eta",100,0.,2.5,wbins2,wxmin2,6.);
TH2F*  wrho_pt=new TH2F("wrho_pt","pt",100,0.,2000.,wbins2,wxmin2,6.);
TH2F* wrho_diphopt= new TH2F("wrho_diphopt","diphopt",100,0.,2000.,wbins2,wxmin2,6.);
TH2F*  wrhorratiobin=new TH2F("wrrhoratiobin","old/new weights",n_rhobins_forreweighting ,0.,n_rhobins_forreweighting+1,wbins2,wxmin2,6.);
TH2F*  wrhorratio=new TH2F("wrrhoratio","old/new weights",50,0., 50.,wbins2,wxmin2,6.);
TH2F*  wrhor_eta=new TH2F("wrhor_eta","eta",100,0.,2.5,wbins2,wxmin2,6.);
TH2F*  wrhor_pt=new TH2F("wrhor_pt","pt",100,0.,2000.,wbins2,wxmin2,6.);
TH2F* wrhor_diphopt= new TH2F("wrhor_diphopt","diphopt",100,0.,2000.,wbins2,wxmin2,6.);

TH2F*  wsigmaratiobin=new TH2F("wsigmaratiobin","old/new weights",n_sigmabins_forreweighting ,0.,n_sigmabins_forreweighting+1,wbins2,wxmin2,6.);
TH2F*  wsigmaratio=new TH2F("wsigmaratio","old/new weights",20,0.,  20.,wbins2,wxmin2,6.);
TH2F*  wsigmarratiobin=new TH2F("wrsigmaratiobin","old/new weights",n_sigmabins_forreweighting ,0.,n_sigmabins_forreweighting+1,wbins2,wxmin2,6.);
TH2F*  wsigmarratio=new TH2F("wrsigmaratio","old/new weights",20,0., 20.,wbins2,wxmin2,6.);


TH2F*  wetaratiobin=new TH2F("wetaratiobin","old/new weights",n_etabins_forreweighting ,0.,n_etabins_forreweighting+1,wbins2,wxmin2,6.);
TH2F*  wetaratio=	new TH2F("wetaratio","old/new weights",25,0.,2.5,wbins2,wxmin2,6.);
TH2F*  weta_rho=	new TH2F("weta_rho","rho",100,0.,50.,wbins2,wxmin2,6.);
TH2F*  weta_sigma=	new TH2F("weta_sigma","sigma",20,0.,20.,wbins2,wxmin2,6.);
TH2F*  weta_pt=		new TH2F("weta_pt","pt",100,0.,2000.,wbins2,wxmin2,6.);
TH2F* weta_diphopt= new TH2F("weta_diphopt","diphopt",100,0.,2000.,wbins2,wxmin2,6.);
TH2F*  wetarratiobin=new TH2F("wetarratiobin","old/new weights",n_etabins_forreweighting ,0.,n_etabins_forreweighting+1,wbins2,wxmin2,6.);
TH2F*  wetarratio=	new TH2F("wetarratio","old/new weights",25,0.,2.5,wbins2,wxmin2,6.);
TH2F*  wetar_rho=	new TH2F("wetar_rho","rho",100,0.,50.,wbins2,wxmin2,6.);
TH2F*  wetar_sigma=	new TH2F("wetar_sigma","sigma",20,0.,20.,wbins2,wxmin2,6.);
TH2F*  wetar_pt=	new TH2F("wetar_pt","pt",100,0.,2000.,wbins2,wxmin2,6.);
TH2F* wetar_diphopt= new TH2F("wetar_diphopt","diphopt",100,0.,2000.,wbins2,wxmin2,6.);

TH2F*  wpt_rho=new TH2F("wpt_rho","rho bkg", 100,0.,50.,wbins2,wxmin2,6.);
TH2F*  wpt_sigma=new TH2F("wpt_sigma","sigma bkg",50.,0.,20.,wbins2,wxmin2,6.);
TH2F * wpt_eta= new TH2F("wpt_eta","eta bkg",50,0.,2.5,wbins2,wxmin2,6.);
TH2F * wpt_diphopt= new TH2F("wpt_diphopt","diphopt bkg",100 ,wxmin,2000.,wbins2,wxmin2,6.);
TH2F*  wptratiobin=new TH2F("wptratiobin","oldw/neww",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1, wbins2,wxmin2,6.);
TH2F*  wptratio=new TH2F("wptratio","oldw/neww",100 ,0.,1000., wbins2,wxmin2,6.);


TH2F*  wdiphopt_rho=new TH2F("wdiphopt_rho","rho bkg", 100,0.,50.,wbins2,wxmin2,6.);
TH2F*  wdiphopt_sigma=new TH2F("wdiphopt_sigma","sigma bkg",50.,0.,20.,wbins2,wxmin2,6.);
TH2F *wdiphopt_eta= new TH2F("wdiphopt_eta","eta bkg",50,0.,2.5,wbins2,wxmin2,6.);
TH2F *wdiphopt_pt= new TH2F("wdiphopt_pt","pt bkg",wbins ,wxmin,wxmax,wbins2,wxmin2,6.);
TH2F*  wdiphoptratiobin=new TH2F("wdiphoptratiobin","oldw/neww",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1, wbins2,wxmin2,6.);
TH2F*  wdiphoptratio=new TH2F("wdiphoptratio","oldw/neww",100 ,0.,1000., wbins2,wxmin2,6.);

//random cone
TH2F*  wptr_rho=new TH2F("wptr_rho","rho sig",100,0.,50.,wbins2,wxmin2,6.);
TH2F*  wptr_sigma=new TH2F("wptr_sigma","sigma sig",50.,0.,20.,wbins2,wxmin2,6.);
TH2F *wptr_eta= new TH2F("wptr_eta","eta sig",50.,0.,2.5,wbins2,wxmin2,6.);
TH2F *wptr_diphopt= new TH2F("wptr_diphopt","diphopt sig",100,wxmin,2000.,wbins2,wxmin2,6.);
TH2F*  wptrratiobin=new TH2F("wptrratiobin","oldw/neww",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1, wbins2,wxmin2,6.);
TH2F*  wptrratio=new TH2F("wptrratio","oldw/neww",100,0.,1000., wbins2,wxmin2,6.);


TH2F*  wdiphoptr_rho=new TH2F("wdiphoptr_rho","rho sig",100,0.,50.,wbins2,wxmin2,6.);
TH2F*  wdiphoptr_sigma=new TH2F("wdiphoptr_sigma","sigma sig",50.,0.,20.,wbins2,wxmin2,6.);
TH2F *wdiphoptr_eta= new TH2F("wdiphoptr_eta","eta sig",50.,0.,2.5,wbins2,wxmin2,6.);
TH2F *wdiphoptr_pt= new TH2F("wdiphoptr_pt","pt sig",wbins ,wxmin,wxmax,wbins2,wxmin2,6.);
TH2F*  wdiphoptrratiobin=new TH2F("wdiphoptrratiobin","oldw/neww",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1, wbins2,wxmin2,6.);
TH2F*  wdiphoptrratio=new TH2F("wdiphoptrratio","oldw/neww",100,0.,1000., wbins2,wxmin2,6.);



TH2F*  twrhoratiobin=new TH2F("twrhoratiobin","old/new weights",n_rhobins_forreweighting ,0.,n_rhobins_forreweighting+1,wbins2,wxmin2,6.);
TH2F*  twrhoratio=new TH2F("twrhoratio","old/new weights",50,0.,  50.,wbins2,wxmin2,6.);
TH2F*  twrho_eta=new TH2F("twrho_eta","eta",100,0.,2.5,wbins2,wxmin2,6.);
TH2F*  twrho_pt=new TH2F("twrho_pt","pt",100,0.,2000.,wbins2,wxmin2,6.);
TH2F* twrho_diphopt= new TH2F("twrho_diphopt","diphopt",100,0.,2000.,wbins2,wxmin2,6.);
TH2F*  twrhorratiobin=new TH2F("twrhorratiobin","old/new weights",n_rhobins_forreweighting ,0.,n_rhobins_forreweighting+1,wbins2,wxmin2,20.);
TH2F*  twrhorratio=new TH2F("twrhorratio","old/new weights",50,0., 50.,wbins2,wxmin2,20.);
TH2F*  twrhor_eta=new TH2F("twrhor_eta","eta",100,0.,2.5,wbins2,wxmin2,6.);
TH2F*  twrhor_pt=new TH2F("twrhor_pt","pt",100,0.,2000.,wbins2,wxmin2,6.);
TH2F* twrhor_diphopt= new TH2F("twrhor_diphopt","diphopt",100,0.,2000.,wbins2,wxmin2,6.);

TH2F*  twsigmaratiobin=new TH2F("twsigmaratiobin","old/new weights",n_sigmabins_forreweighting ,0.,n_sigmabins_forreweighting+1,wbins2,wxmin2,6.);
TH2F*  twsigmaratio=new TH2F("twsigmaratio","old/new weights",20,0.,  20.,wbins2,wxmin2,6.);
TH2F*  twsigmarratiobin=new TH2F("twsigmarratiobin","old/new weights",n_sigmabins_forreweighting ,0.,n_sigmabins_forreweighting+1,wbins2,wxmin2,6.);
TH2F*  twsigmarratio=new TH2F("twsigmarratio","old/new weights",20,0., 20.,wbins2,wxmin2,6.);








TH2F*  twetaratiobin=new TH2F("twetaratiobin","old/new weights",n_etabins_forreweighting ,0.,n_etabins_forreweighting+1,wbins2,wxmin2,6.);
TH2F*  twetaratio=	new TH2F("twetaratio","old/new weights",25,0.,2.5,wbins2,wxmin2,6.);
TH2F*  tweta_rho=	new TH2F("tweta_rho","rho",100,0.,50.,wbins2,wxmin2,6.);
TH2F*  tweta_sigma=	new TH2F("tweta_sigma","sigma",20,0.,20.,wbins2,wxmin2,6.);
TH2F*  tweta_pt=	new TH2F("tweta_pt","pt",100,0.,2000.,wbins2,wxmin2,6.);
TH2F*  tweta_diphopt= new TH2F("tweta_diphopt","diphopt",100,0.,2000.,wbins2,wxmin2,6.);
TH2F*  twetarratiobin=new TH2F("twetarratiobin","old/new weights",n_etabins_forreweighting ,0.,n_etabins_forreweighting+1,wbins2,wxmin2,6.);
TH2F*  twetarratio=	new TH2F("twetarratio","old/new weights",25,0.,2.5,wbins2,wxmin2,2.);
TH2F*  twetar_rho=	new TH2F("twetar_rho","rho",100,0.,50.,wbins2,wxmin2,6.);
TH2F*  twetar_sigma=new TH2F("twetar_sigma","sigma",20.,0.,20.,wbins2,wxmin2,6.);
TH2F*  twetar_pt=	new TH2F("twetar_pt","pt",100,0.,2000.,wbins2,wxmin2,6.);
TH2F*  twetar_diphopt= new TH2F("twetar_diphopt","diphopt",100,0.,2000.,wbins2,wxmin2,6.);




TH2F*  twpt_rho=new TH2F("twpt_rho","rho bkg",100.,0.,50.,wbins2,wxmin2,6.);
TH2F*  twpt_sigma=new TH2F("twpt_sigma","sigma bkg",20,0.,20.,wbins2,wxmin2,6.);
TH2F *twpt_eta= new TH2F("twpt_eta","eta bkg",50,0.,2.5,wbins2,wxmin2,6.);
TH2F *twpt_diphopt= new TH2F("twpt_diphopt","diphopt bkg",100 ,wxmin,2000.,wbins2,wxmin2,6.);
TH2F*  twptratiobin=new TH2F("twptratiobin","oldw/neww",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1, wbins2,wxmin2,6.);
TH2F*  twptratio=new TH2F("twptratio","oldw/neww",100 ,0.,1000., wbins2,wxmin2,6.);


TH2F*  twptr_rho=new TH2F("twptr_rho","rho sig",100,0.,50.,wbins2,wxmin2,6.);
TH2F*  twptr_sigma=new TH2F("twptr_sigma","sigma sig",20.,0.,20.,wbins2,wxmin2,6.);
TH2F *twptr_eta= new TH2F("twptr_eta","eta sig",50.,0.,2.5,wbins2,wxmin2,6.);
TH2F *twptr_diphopt= new TH2F("twptr_diphopt","diphopt sig",100,wxmin,2000.,wbins2,wxmin2,6.);
TH2F*  twptrratiobin=new TH2F("twptrratiobin","oldw/neww",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1, wbins2,wxmin2,6.);
TH2F*  twptrratio=new TH2F("twptrratio","oldw/neww",100.,0.,1000., wbins2,wxmin2,6.);


TH2F*  twdiphopt_rho=new TH2F("twdiphopt_rho","rho bkg",100.,0.,50.,wbins2,wxmin2,6.);
TH2F*  twdiphopt_sigma=new TH2F("twdiphopt_sigma","sigma bkg",20,0.,20.,wbins2,wxmin2,6.);
TH2F *twdiphopt_eta= new TH2F("twdiphopt_eta","eta bkg",50,0.,2.5,wbins2,wxmin2,6.);
TH2F *twdiphopt_pt= new TH2F("twdiphopt_diphopt","pt bkg",wbins ,wxmin,wxmax,wbins2,wxmin2,6.);
TH2F*  twdiphoptratiobin=new TH2F("twdiphoptratiobin","oldw/neww",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1, wbins2,wxmin2,6.);
TH2F*  twdiphoptratio=new TH2F("twdiphoptratio","oldw/neww",100 ,0.,1000., wbins2,wxmin2,6.);


TH2F*  twdiphoptr_rho=new TH2F("twdiphoptr_rho","rho ",100,0.,50.,wbins2,wxmin2,6.);
TH2F*  twdiphoptr_sigma=new TH2F("twdiphoptr_sigma","sigma sig",20.,0.,20.,wbins2,wxmin2,6.);
TH2F *twdiphoptr_eta= new TH2F("twdiphoptr_eta","eta sig",50.,0.,2.5,wbins2,wxmin2,6.);
TH2F *twdiphoptr_pt= new TH2F("twdiphoptr_pt","pt sig",wbins ,wxmin,wxmax,wbins2,wxmin2,6.);
TH2F*  twdiphoptrratiobin=new TH2F("twdiphoptrratiobin","oldw/neww",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1, wbins2,wxmin2,6.);
TH2F*  twdiphoptrratio=new TH2F("twdiphoptrratio","oldw/neww",100.,0.,1000., wbins2,wxmin2,6.);

//define functions for fit
void get_roodset_from_ttree(TDirectoryFile *f, TString treename, RooDataSet* &roodset);
//int do1dfit(RooDataSet *dset_data_axis1,RooDataSet *dset_sig_axis1,RooDataSet *dset_bkg_axis1,int im);
int do1dfit(RooDataSet *dset_data_axis1,RooDataSet* dset_sigcut, RooDataSet* dset_bkgcut,RooDataSet *dset_sig_axis1,RooDataSet *dset_bkg_axis1,RooDataSet* dset_testsig,RooDataSet* dset_testbkg, int im, double pu[], double puerr[],double pullerr[],double data_en[]);
/*
void reweight_rhosigma(RooDataSet **dset, RooDataSet *dsetdestination, bool deleteold = kTRUE);
void reweight_pt_1d(RooDataSet **dset, RooDataSet *dsetdestination, int numvar);
void reweight_diphotonpt_1d(RooDataSet **dset, RooDataSet *dsetdestination);
void reweight_eta_1d(RooDataSet **dset, RooDataSet *dsetdestination, int numvar);
*/

void reweight_rhosigma(TH2F* weight_eta=NULL, TH2F* weight_pt=NULL,TH2F*weight_diphopt=NULL,TH2F*weightratiobin_rho=NULL,TH2F* weightratio_rho=NULL, TH2F* weightratiobin_sigma=NULL,TH2F*weightratio_sigma=NULL,RooDataSet **dset=NULL, RooDataSet *dsetdestination=NULL, bool deleteold = kTRUE);
void reweight_eta_1d(TH2F* weight_rho=NULL,TH2F* weight_sigma=NULL,TH2F*weight_pt=NULL,TH2F* weight_diphopt=NULL,TH2F* weightratiobin=NULL,TH2F*weightratio=NULL,RooDataSet **dset=NULL, RooDataSet *dsetdestination=NULL, int numvar=0);
void reweight_pt_1d(TH2F* weight_rho=NULL,TH2F* weight_sigma=NULL,TH2F*weight_eta=NULL,TH2F* weight_diphopt=NULL,TH2F* weightratiobin_pt=NULL,TH2F*weightratio_pt=NULL,RooDataSet **dset=NULL, RooDataSet *dsetdestination=NULL, int numvar=0);
void reweight_diphotonpt_1d(TH2F* weight_rho=NULL,TH2F* weight_sigma=NULL,TH2F*weight_eta=NULL,TH2F* weight_pt=NULL,TH2F* weightratiobin=NULL,TH2F*weightratio=NULL,RooDataSet **dset=NULL, RooDataSet *dsetdestination=NULL);


void ratioplot(RooDataSet *data_axis1, RooDataSet *truth_axis1, Bool_t logplot);
//double quantiles(TH1F* mass, int ntot, int nqtemp,int startn, int endn, TString eta_q,double _proptemp[], double dpmqtemp[]);
int quantiles(TH1D* mass,double _proptemp[], double dpmqtemp[]);
int plot_purity_massbin(double mass[],double masserr[],double pur[],double purerr[], bool errb) ; 
int prep_dataset(Bool_t isdata=kTRUE,TString eta="EBEB", Bool_t truthf=kFALSE,int ntotbins=10,int startbin=0, int endbin=2);
 void print_mem();
const int numcpu=1;
TString cut_s;
TString date;
ProcInfo_t procinfo;
