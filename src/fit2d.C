#include <assert.h>
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"
#include "RooFitResult.h"
#include "TLatex.h"
#include "TString.h"
#include "TFile.h"
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
#include "RooFitResult.h"
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
using namespace std; 
using namespace RooFit;
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
RooDataSet *dataset_test = NULL;
RooDataSet *dset_data = NULL;
RooDataSet *dset_rcone = NULL;
RooDataSet *dset_sideb = NULL;

RooRealVar *binning_roovar1=NULL;

bool isdata=kTRUE;
//bins for reweihgting
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
//const int n_rhobins_forreweighting = 4;
//Float_t rhobins_forreweighting[n_rhobins_forreweighting+1]={0,10,15,20,60};
const int n_rhobins_forreweighting = 3;
Float_t rhobins_forreweighting[n_rhobins_forreweighting+1]={0,10,15,60};
const int n_sigmabins_forreweighting = 3;
Float_t sigmabins_forreweighting[n_sigmabins_forreweighting+1]={0,2,3,20};


//histograms for reweightings
double wxmin=0.;
double wxmax=10.;
int wbins=20;
int wbins2=20;
double wxmin2=0.;
double wxmax2=50.;
/*TH2F*  w_rhonr=new TH2F("w_rhonr","w_rhonr",n_rhobins_forreweighting ,0.,n_rhobins_forreweighting+1,wbins2,wxmin2,20.);
TH2F*  w_rho2o=new TH2F("w_rho2o","w_rho2o",n_rhobins_forreweighting ,0.,n_rhobins_forreweighting+1,wbins2,wxmin2,25.);
TH2F*  w_rhon=new TH2F("w_rhon","w_rhon",n_rhobins_forreweighting ,0.,n_rhobins_forreweighting+1,wbins,wxmin,100.);
TH2F*  w_rho2=new TH2F("w_rho2","w_rho2",50,0.,50.,wbins2,wxmin2,15);
TH1F *w_rho= new TH1F("w_rho","w_rho",wbins ,wxmin,100.);
TH1F *w_rhoo= new TH1F("w_rhoo","w_rhoo",wbins ,wxmin,25.);
TH2F*  w_sigmanr=new TH2F("w_sigmanr","w_sigmanr",n_sigmabins_forreweighting ,0.,n_sigmabins_forreweighting+1,wbins2,wxmin2,20.);
TH2F*  w_sigma2o=new TH2F("w_sigma2o","w_sigma2o",n_sigmabins_forreweighting ,0.,n_sigmabins_forreweighting+1,wbins2,wxmin2,25.);
TH2F*  w_sigman=new TH2F("w_sigman","w_sigman",n_sigmabins_forreweighting ,0.,n_sigmabins_forreweighting+1,wbins,wxmin,100.);
TH2F*  w_sigma2=new TH2F("w_sigma2","w_sigma2",50,0.,50.,wbins2,wxmin2,15);
TH1F *w_sigma= new TH1F("w_sigma","w_sigma",wbins ,wxmin,100.);
TH1F *w_sigmao= new TH1F("w_sigmao","w_sigmao",wbins ,wxmin,25.);
*/
TH2F*  w_rhornr=new TH2F("w_rhornr","w_rhornr",n_rhobins_forreweighting ,0.,n_rhobins_forreweighting+1,wbins2,wxmin2,20.);
TH2F*  w_rhor2o=new TH2F("w_rhor2o","w_rhor2o",n_rhobins_forreweighting ,0.,n_rhobins_forreweighting+1,wbins2,wxmin2,25.);
TH2F*  w_rhorn=new TH2F("w_rhorn","w_rhorn",n_rhobins_forreweighting ,0.,n_rhobins_forreweighting+1,wbins,wxmin,100.);
TH2F*  w_rhor2=new TH2F("w_rhor2","w_rhor2",50,0.,50.,wbins2,wxmin2,15);
TH1F *w_rhor= new TH1F("w_rhor","w_rhor",wbins ,wxmin,100.);
TH1F *w_rhoro= new TH1F("w_rhoro","w_rhoro",wbins ,wxmin,25.);
TH2F*  w_sigmarnr=new TH2F("w_sigmarnr","w_sigmarnr",n_sigmabins_forreweighting ,0.,n_sigmabins_forreweighting+1,wbins2,wxmin2,20.);
TH2F*  w_sigmar2o=new TH2F("w_sigmar2o","w_sigmar2o",n_sigmabins_forreweighting ,0.,n_sigmabins_forreweighting+1,wbins2,wxmin2,25.);
TH2F*  w_sigmarn=new TH2F("w_sigmarn","w_sigmarn",n_sigmabins_forreweighting ,0.,n_sigmabins_forreweighting+1,wbins,wxmin,100.);
TH2F*  w_sigmar2=new TH2F("w_sigmar2","w_sigmar2",50,0.,50.,wbins2,wxmin2,15);


TH2F*  w_etanr=new TH2F("w_etanr","w_etanr",n_etabins_forreweighting ,0.,n_etabins_forreweighting+1,wbins2,wxmin2,3.);
TH2F*  w_eta2o=new TH2F("w_eta2o","w_eta2o",n_etabins_forreweighting ,0.,n_etabins_forreweighting+1,wbins2,wxmin2,10.);
TH2F*  w_etan=new TH2F("w_etan","w_etan",n_etabins_forreweighting ,0.,n_etabins_forreweighting+1,wbins,wxmin,15.);
TH2F*  w_eta2=new TH2F("w_eta2","w_eta2",n_etabins_forreweighting,0.,2.5,wbins2,wxmin2,2.);
TH1F *w_eta= new TH1F("w_eta","w_eta",wbins ,wxmin,15.);
TH1F *w_etao= new TH1F("w_etao","w_etao",wbins ,wxmin,10.);
TH2F*  w_etarnr=new TH2F("w_etarnr","w_etarnr",n_etabins_forreweighting ,0.,25,wbins2,wxmin2,3.);
TH2F*  w_etar2o=new TH2F("w_etar2o","w_etar2o",n_etabins_forreweighting ,0.,25,wbins2,wxmin2,10.);
TH2F*  w_etarn=new TH2F("w_etarn","w_etarn",n_etabins_forreweighting ,0.,25,wbins,wxmin,15.);
TH2F*  w_etar2=new TH2F("w_etar2","w_etar2",n_etabins_forreweighting,0.,2.5,wbins2,wxmin2,2.);
TH1F *w_etar= new TH1F("w_etar","w_etar",wbins ,wxmin,15.);
TH1F *w_etaro= new TH1F("w_etaro","w_etaro",wbins ,wxmin,10.);



TH2F*  w_rhonr=new TH2F("w_rhonr","w_rhonr",100 ,0.,100+1,20,wxmin2,2.);
TH2F*  w_rho2o=new TH2F("w_rho2o","w_rho2o",100 ,0.,100+1,10,wxmin2,10.);
TH2F*  w_rhon=new TH2F("w_rhon","w_rhon",100 ,0.,100+1,10,wxmin,10.);
TH2F*  w_rho2=new TH2F("w_rho2","w_rho2",100,0.,100.+1,20,wxmin2,2.);
TH1F *w_rho= new TH1F("w_rho","w_rho",wbins ,wxmin,10.);
TH1F *w_rhoo= new TH1F("w_rhoo","w_rhoo",wbins ,wxmin,10.);
TH2F*  w_sigmanr=new TH2F("w_sigmanr","w_sigmanr",20 ,0.,20+1,20,wxmin2,2.);
TH2F*  w_sigma2o=new TH2F("w_sigma2o","w_sigma2o",20 ,0.,20+1,10,wxmin2,10.);
TH2F*  w_sigman=new TH2F("w_sigman","w_sigman",20 ,0.,20+1,10,wxmin,10.);
TH2F*  w_sigma2=new TH2F("w_sigma2","w_sigma2",20,0.,20.,wbins2,wxmin2,2.);
//TH1F*  w_sigmasigma=new TH1F("w_sigmasigma","w_sigmasigma",wbins,wxmin,wxmax);
//TH2F*  w_rhorsigmar2=new TH2F("w_rhorsigmar2","w_rhorsigmar2",50,0.,50.,wbins2,wxmin2,wxmax2);
//TH1F*  w_sigmarsigmar=new TH1F("w_sigmarsigmar","w_sigmarsigmar",wbins,wxmin,wxmax);
//TH2F*  w_rhosigmar2=new TH2F("w_rhosigmar2","w_rhosigmar2",50,0.,50.,wbins2,wxmin2,wxmax2);
//TH2F*  w_rhosigmar3=new TH2F("w_rhosigmar3","w_rhosigmar3",20,0.,20.,wbins2,wxmin2,wxmax2);
//TH2F*  w_rhosigmar3=new TH2F("w_rhosigmar3","w_rhosigmar3",50,0.,50.,20,0.,20.);
//TH1F* w_eta=new TH1F("w_eta","w_eta",wbins,wxmin,wxmax);
//TH1F* w_eta=new TH1F("w_eta","w_eta",25,0.,2.5);
//TH2F*  w_eta2=new TH2F("w_eta2","w_eta2",25,0.,2.5,wbins2,wxmin2,wxmax2);
TH2F*  w_diphoptnr=new TH2F("w_diphoptnr","w_diphoptnr",n_diphoptbins_forreweighting ,0.,n_diphoptbins_forreweighting+1,wbins2,wxmin2,20.);
TH2F*  w_diphopt2o=new TH2F("w_diphopt2o","w_diphopt2o",n_diphoptbins_forreweighting ,0.,n_diphoptbins_forreweighting+1,wbins2,wxmin2,25.);
TH2F*  w_diphoptn=new TH2F("w_diphoptn","w_diphoptn",n_diphoptbins_forreweighting ,0.,n_diphoptbins_forreweighting+1,wbins,wxmin,100.);
TH2F*  w_diphopt2=new TH2F("w_diphopt2","w_diphopt2",100,0.,1000,wbins2,wxmin2,15);
TH2F*  w_ptn=new TH2F("w_ptn","w_ptn",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1,wbins,wxmin,25.);
TH2F*  w_ptnr=new TH2F("w_ptnr","w_ptnr",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1, wbins2,wxmin2,10.);
TH2F*  w_pt2o=new TH2F("w_pt2o","w_pt2o",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1, wbins2,wxmin2,15.);
TH2F*  w_pt2=new TH2F("w_pt2","w_pt2",500 ,0.,500,wbins2,wxmin2,5.);
TH1F *w_pt= new TH1F("w_pt","w_pt",wbins ,wxmin,25.);
TH1F *w_diphopt= new TH1F("w_diphopt","w_diphopt",wbins ,wxmin,100.);
TH1F *w_diphopto= new TH1F("w_diphopto","w_diphopto",wbins ,wxmin,25.);
TH1F *w_pto= new TH1F("w_pto","w_pto",wbins ,wxmin,12.);
TH2F*  w_diphoptrnr=new TH2F("w_diphoptrnr","w_diphoptrnr",n_diphoptbins_forreweighting ,0.,n_diphoptbins_forreweighting+1,wbins2,wxmin2,20.);
TH2F*  w_diphoptr2o=new TH2F("w_diphoptr2o","w_diphoptr2o",n_diphoptbins_forreweighting ,0.,n_diphoptbins_forreweighting+1,wbins2,wxmin2,25.);
TH2F*  w_diphoptrn=new TH2F("w_diphoptrn","w_diphoptrn",n_diphoptbins_forreweighting ,0.,n_diphoptbins_forreweighting+1,wbins,wxmin,100.);
TH2F*  w_diphoptr2=new TH2F("w_diphoptr2","w_diphoptr2",100,0.,1000,wbins2,wxmin2,15);
TH2F*  w_ptrn=new TH2F("w_ptrn","w_ptrn",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1,wbins,wxmin,25.);
TH2F*  w_ptrnr=new TH2F("w_ptrnr","w_ptrnr",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1, wbins2,wxmin2,10.);
TH2F*  w_ptr2o=new TH2F("w_ptr2o","w_ptr2o",n_ptbins_forreweighting ,0.,n_ptbins_forreweighting+1, wbins2,wxmin2,15.);
TH2F*  w_ptr2=new TH2F("w_ptr2","w_ptr2",500 ,0.,500,wbins2,wxmin2,5.);
TH1F *w_ptr= new TH1F("w_ptr","w_ptr",wbins ,wxmin,25.);
TH1F *w_diphoptr= new TH1F("w_diphoptr","w_diphoptr",wbins ,wxmin,100.);
TH1F *w_diphoptro= new TH1F("w_diphoptro","w_diphoptro",wbins ,wxmin,25.);
TH1F *w_ptro= new TH1F("w_ptro","w_ptro",wbins ,wxmin,12.);
  
//define functions for fit
void get_roodset_from_ttree(TDirectoryFile *f, TString treename, RooDataSet* &roodset);
int do1dfit();
int do2dfit();
//void reweight_rhosigma(TH1F* w_rhosigma,TH2F*w_rhosigma2,TH2F*w_rhosigma3,RooDataSet **dset, RooDataSet *dsetdestination, bool deleteold = kTRUE);
void reweight_rhosigma(TH1F* weight_rho, TH1F* weight_rhoo,TH2F*weight_rhon,TH2F*weight_rho2o,TH2F* weight_rhonr, TH2F* weight_rho2,TH2F*weight_sigman,TH2F* weight_sigma2o,TH2F* weight_sigmanr, TH2F* weight_sigma2,RooDataSet **dset, RooDataSet *dsetdestination, bool deleteold = kTRUE);
void reweight_pt_1d(TH1F* weight_pt,TH1F* weight_pto,TH2F*weight_ptn,TH2F* weight_pt2o,TH2F* weight_ptnr,TH2F*weight_pt2,RooDataSet **dset, RooDataSet *dsetdestination, int numvar);
void reweight_diphotonpt_1d(TH1F* weight_diphopt,TH1F* weight_diphopto,TH2F*weight_diphoptn,TH2F* weight_diphopt2o,TH2F*weight_diphoptnr,TH2F*weight_diphopt2,RooDataSet **dset, RooDataSet *dsetdestination);
void reweight_eta_1d(TH1F* weight_eta,TH1F* weight_etao,TH2F*weight_etan,TH2F* weight_eta2o,TH2F* weight_etanr,TH2F*weight_eta2,RooDataSet **dset, RooDataSet *dsetdestination, int numvar);
void ratioplot(RooDataSet *data_axis1, RooDataSet *truth_axis1, Bool_t logplot);
void quantiles(TH1F* mass);

void quantiles(TH1F* mass) {
        // demo for quantiles
        const Int_t nq = 10;
        mass->Scale(1.0/mass->Integral());
        Double_t prob[nq];  // position where to compute the quantiles in [0,1]
        Double_t yq[nq];  // array to contain the quantiles
        for (Int_t i=0;i<nq;i++) 
        //xq[i] = Float_t(i+1)/nq;
    //    prob[i] =  (0.1+ i*0.1);
        prob[i] =  Float_t(i+1)/nq;
        mass->GetQuantiles(nq,yq,prob);

        //show the original histogram in the top pad
        TCanvas *cq = new TCanvas("cq","demo quantiles",10,10,700,900);
        cq->Divide(1,2);
        cq->cd(1);
        mass->Draw();

        // show the quantiles in the bottom pad
        cq->cd(2);
        gPad->SetGrid();
        TGraph *gr = new TGraph(nq,prob,yq);
        gr->SetMarkerStyle(21);
        gr->Draw("alp");
     }



//get the RooDataSet
void get_roodset_from_ttree(TDirectoryFile *f, TString treename, RooDataSet* &roodset){
  cout << "Creating roodset from file: " << f->GetName() << " with tree: " << treename.Data() << endl;

  TTree *t = NULL;
  assert(roodset==NULL);
  f->GetObject(treename.Data(),t);
  if (!t) {cout << "Impossible to find TTree " << treename.Data() << endl; return;}
  TObjArray *objs = t->GetListOfBranches();
  //disables all branches
  t->SetBranchStatus("*",0);


  float v_rooisopv1;
  float v_rooisopv2;

  float v_rooisowv1;
  float v_rooisowv2;
  float v_roovar1;
  float v_roovar2;
  float v_roopt1;
  float v_roosieie1;
  float v_rooeta1;
  float v_roopt2;
  float v_roosieie2;
  float v_rooeta2;
  float v_roodiphopt;
  float v_roorho;
  float v_roosigma;
  float v_roonvtx;
  float v_rooweight;

  TBranch *b_roovar1;
  TBranch *b_roovar2;

  TBranch *b_rooisopv1;
  TBranch *b_rooisopv2;
  TBranch *b_rooisowv1;
  TBranch *b_rooisowv2;
  TBranch *b_roopt1;
  TBranch *b_roosieie1;
  TBranch *b_rooeta1;
  TBranch *b_roopt2;
  TBranch *b_roosieie2;
  TBranch *b_rooeta2;
  TBranch *b_roodiphopt;
  TBranch *b_roorho;
  TBranch *b_roosigma;
  TBranch *b_roonvtx;
  TBranch *b_rooweight;

  const int nvars = 17;
  float* ptrs[nvars]={&v_roovar1,&v_roovar2,&v_rooisopv1,&v_rooisopv2,&v_rooisowv1,&v_rooisowv2,&v_roopt1,&v_roosieie1,&v_rooeta1,&v_roopt2,&v_roosieie2,&v_rooeta2,&v_roodiphopt,&v_roorho,&v_roosigma,&v_roonvtx,&v_rooweight};
  TBranch** branches[nvars]={&b_roovar1,&b_roovar2,&b_rooisopv1,&b_rooisopv2,&b_rooisowv1,&b_rooisowv2,&b_roopt1,&b_roosieie1,&b_rooeta1,&b_roopt2,&b_roosieie2,&b_rooeta2,&b_roodiphopt,&b_roorho,&b_roosigma,&b_roonvtx,&b_rooweight};
  RooRealVar* rooptrs[nvars]={roovar1,roovar2,rooisopv1,rooisopv2,rooisowv1,rooisowv2,roopt1,roosieie1,rooeta1,roopt2,roosieie2,rooeta2,roodiphopt,roorho,roosigma,roonvtx,rooweight};
  bool status[nvars];
  RooArgSet args;
  for (int i=0; i<nvars; i++){ 
    status[i]=0;
    TString name = rooptrs[i]->GetName();
    TObject *obj = objs->FindObject(name.Data());
    if (!obj) continue;
    t->SetBranchStatus(name.Data(),1);
    status[i]=1;
    t->SetBranchAddress(name.Data(),ptrs[i],branches[i]);
 
   args.add(*(rooptrs[i]));
  }

  TString newname = Form("roo_%s",t->GetName());
  roodset = new RooDataSet(newname.Data(),newname.Data(),args,WeightVar(*rooweight) );
  for (int j=0; j<t->GetEntries(); j++){
    t->GetEntry(j);
    for (int i=0; i<nvars; i++){
      if (!status[i]) continue;
      rooptrs[i]->setVal(*(ptrs[i]));
    }
    roodset->add(args,v_rooweight);
  }

  cout << "Imported roodset " << newname.Data() << " from TTree " << t->GetName() << endl;
  roodset->Print();

}


void reweight_rhosigma(TH1F* weight_rho, TH1F* weight_rhoo,TH2F*weight_rhon,TH2F*weight_rho2o,TH2F* weight_rhonr, TH2F* weight_rho2,TH2F*weight_sigman,TH2F*weight_sigma2o,TH2F* weight_sigmanr, TH2F* weight_sigma2,RooDataSet **dset, RooDataSet *dsetdestination, bool deleteold){
  if (!(*dset)) return;
  (*dset)->Print();

//  TH2F *hnum = new TH2F("hnum","hnum",n_rhobins_forreweighting,rhobins_forreweighting,n_sigmabins_forreweighting,sigmabins_forreweighting);
//  TH2F *hden = new TH2F("hden","hden",n_rhobins_forreweighting,rhobins_forreweighting,n_sigmabins_forreweighting,sigmabins_forreweighting);
 // TH2F *hnum = new TH2F("hnum","hnum",60,0,60,15,0,15);
 // TH2F *hden = new TH2F("hden","hden",60,0,60,15,0,15);
  TH2F *hnum = new TH2F("hnum","hnum",100,0,100,20,0,20);
  TH2F *hden = new TH2F("hden","hden",100,0,100,20,0,20);
  hnum->Sumw2();
  hden->Sumw2();
  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill(fabs((*dset)->get(i)->getRealValue("roorho")),fabs((*dset)->get(i)->getRealValue("roosigma")),(*dset)->store()->weight(i));
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(fabs(dsetdestination->get(i)->getRealValue("roorho")),fabs(dsetdestination->get(i)->getRealValue("roosigma")),dsetdestination->store()->weight(i));
  }


  hnum->Scale(1.0/hnum->Integral());
  //TH2F *h_data = hnum;
  hden->Scale(1.0/hden->Integral());
//data/MC
  hnum->Divide(hden);
  TH2F *h = hnum;

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_rhosigmarew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float rho = args.getRealValue("roorho");
    float sigma = args.getRealValue("roosigma");
    float neww = oldw*h->GetBinContent(h->FindBin(rho,sigma));
   //     std::cout << rho << " " << sigma << std::endl;
    //    std::cout << oldw << " " << neww << std::endl;
//    weight_rhosigma->Fill(neww);	
 //   weight_rhosigma2->Fill(rho,neww/oldw);
   // weight_rhosigma3->Fill(sigma,neww/oldw);
  //  weight_rhosigma3->Fill(rho,sigma);
//weight_rhosigma2->Fill(neww, args.getRealValue("roovar1"));
    weight_rho->Fill(neww);
    weight_rhoo->Fill(oldw);
   weight_rho2o->Fill(h->GetXaxis()->FindBin(rho),oldw);	
   weight_rhon->Fill(h->GetXaxis()->FindBin(rho),neww);	
if(oldw!=0)weight_rhonr->Fill(h->GetXaxis()->FindBin(rho),oldw/neww);
 else {weight_rhonr->Fill(-10,1);}//cout << "dipho weight old 0" << endl;}
if(oldw!=0)weight_rho2->Fill(rho,oldw/neww);
   weight_sigma2o->Fill(h->GetYaxis()->FindBin(sigma),oldw);	
   weight_sigman->Fill(h->GetYaxis()->FindBin(sigma),neww);	
if(oldw!=0)weight_sigmanr->Fill(h->GetYaxis()->FindBin(sigma),oldw/neww);
 else {weight_sigmanr->Fill(-10,1);}//cout << "dipho weight old 0" << endl;}
if(oldw!=0)weight_sigma2->Fill(sigma,oldw/neww);
	newdset->add(args,neww);
  }
  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

//  weight_rhosigma3->Scale(1.0/weight_rhosigma3->Integral());
  //weight_rhosigma3->Divide(h_data);
  delete hnum; delete hden;
  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "RhoSigma2D rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;
  (*dset)->Print();
  if (deleteold) delete old_dset;
};


void reweight_pt_1d(TH1F* weight_pt, TH1F* weight_pto,TH2F*weight_ptn,TH2F*weight_pt2o,TH2F* weight_ptnr, TH2F* weight_pt2,RooDataSet **dset, RooDataSet *dsetdestination, int numvar){

  if (!(*dset)) return;
///numerator and denominator
  TH1F *hnum = new TH1F("hnum","hnum",n_ptbins_forreweighting,ptbins_forreweighting);
  TH1F *hden = new TH1F("hden","hden",n_ptbins_forreweighting,ptbins_forreweighting);
  hnum->Sumw2();
  hden->Sumw2();

  TH1F *h_data = new TH1F("h_data","h_data",n_ptbins_forreweighting,ptbins_forreweighting);
  h_data->Sumw2();
  const char* ptname=Form("roopt%d",numvar);
// RooAbsData->get*() Load a given row of data
 //getRealValue Get value of a RooAbsReal stored in set with given name. If none is found, value of defVal is returned.
  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill(fabs((*dset)->get(i)->getRealValue(ptname)),(*dset)->store()->weight(i));
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(fabs(dsetdestination->get(i)->getRealValue(ptname)),dsetdestination->store()->weight(i));
  }
  
 // for (int i=0; i<dsetdestination->numEntries(); i++){
  //  h_data->Fill(fabs(dsetdestination->get(i)->getRealValue(ptname)));
 // }
//MQ
 //normalize to one 
  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());
  hnum->Divide(hden);
  TH1F *h = hnum;
  for (int i=0; i<dsetdestination->numEntries(); i++){
   h_data->Fill(fabs(dsetdestination->get(i)->getRealValue(ptname)));
  }


  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_ptrew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float pt = args.getRealValue(ptname);
    float neww = oldw*h->GetBinContent(h->FindBin(fabs(pt)));
    weight_pt->Fill(neww);	
	weight_pto->Fill(oldw);
 //   weight_ptn->Fill((h->GetBinContent(h->FindBin(fabs(pt)))),neww);	
    weight_pt2o->Fill(h->FindBin(fabs(pt)),oldw);
    weight_ptn->Fill(h->FindBin(fabs(pt)),neww);	
//    weight_ptnr->Fill(h->FindBin(fabs(pt)),oldw/neww);	
   if(oldw!=0 && neww!=0)weight_ptnr->Fill(h->FindBin(fabs(pt)),oldw/neww);
  else {weight_ptnr->Fill(-10,1);}
   // weight_pt2->Fill(pt,neww/oldw);
   if(oldw!=0 && neww!=0)weight_pt2->Fill(pt,oldw/neww);
  else {weight_pt2->Fill(-10,1);}
//weight_pt2->Fill(neww, args.getRealValue("roovar1"));
//    std::cout << oldw << " " << neww << std::endl;
  // cout << __LINE__ << endl; 
    newdset->add(args,neww);
  }

// for (int i=0; i<n_ptbins_forreweighting; i++){
// h_data->Fill(h->GetBinContent(i),);
// }


  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());
  //weight_pt->Scale(1.0/weight_pt->Integral());
// weight_pt->Divide(h_data);
  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Pt 1d rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;

  delete old_dset;

};

void reweight_diphotonpt_1d(TH1F*weight_diphopt,TH1F*weight_diphopto,TH2F*weight_diphoptn,TH2F*weight_diphopt2o,TH2F*weight_diphoptnr, TH2F* weight_diphopt2,RooDataSet **dset, RooDataSet *dsetdestination){

  if (!(*dset)) return;
if (!(dsetdestination)) return;
///numerator and denominator
  TH1F *hnum = new TH1F("hnum","hnum",n_diphoptbins_forreweighting,diphoptbins_forreweighting);
  TH1F *hden = new TH1F("hden","hden",n_diphoptbins_forreweighting,diphoptbins_forreweighting);
  hnum->Sumw2();
  hden->Sumw2();

     
// RooAbsData->get*() Load a given row of data

 //getRealValue Get value of a RooAbsReal stored in set with given name. If none is found, value of defVal is returned.
 for (int i=0; i<(*dset)->numEntries(); i++){
      hden->Fill(fabs((*dset)->get(i)->getRealValue("roodiphopt")),(*dset)->store()->weight(i));
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(fabs(dsetdestination->get(i)->getRealValue("roodiphopt")),dsetdestination->store()->weight(i));
  }
//MQ
 //normalize to one 
  hnum->Scale(1.0/hnum->Integral());
  //TH1F *h_data = hnum;
  hden->Scale(1.0/hden->Integral());
  hnum->Divide(hden);
  TH1F *h = hnum;
  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_diphoptrew",(*dset)->GetName()));
  newdset->reset();
  assert(newdset);
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float diphopt = fabs((*dset)->get(i)->getRealValue("roodiphopt"));
    float neww = oldw*h->GetBinContent(h->FindBin(fabs(diphopt)));
    newdset->add(args,neww);
// cout << __LINE__ << endl; 
    weight_diphopt->Fill(neww);
    weight_diphopto->Fill(oldw);
    //weight_diphopt->Fill(neww);	
   // weight_diphopt->Fill(diphopt);	
  //  weight_diphoptnr->Fill((h->GetBinContent(h->FindBin(fabs(diphopt)))),neww/oldw);	
 //   weight_diphoptn->Fill((h->GetBinContent(h->FindBin(fabs(diphopt)))),neww);	
   weight_diphopt2o->Fill((h->FindBin(fabs(diphopt))),oldw);	
   weight_diphoptn->Fill((h->FindBin(fabs(diphopt))),neww);	
 ///   weight_diphopt2->Fill(diphopt, neww/oldw);
if(oldw!=0)weight_diphoptnr->Fill(h->FindBin(fabs(diphopt)),oldw/neww);
 else {weight_diphoptnr->Fill(-10,1);}//cout << "dipho weight old 0" << endl;}
if(oldw!=0)weight_diphopt2->Fill(diphopt,oldw/neww);
 else {weight_diphopt2->Fill(-10,1);}//cout << "dipho weight old 0" << endl;}
  // weight_diphopt2->Fill(neww, args.getRealValue("roovar1"));
  }

  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());
 // weight_diphopt->Scale(1.0/weight_diphopt->Integral());
 // weight_diphopt->Divide(h_data);
  delete hnum; delete hden;
 
  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Diphoton Pt 1d rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;

  delete old_dset;

};








void reweight_eta_1d(TH1F* weight_eta,TH1F* weight_etao,TH2F*weight_etan,TH2F* weight_eta2o,TH2F* weight_etanr,TH2F*weight_eta2,RooDataSet **dset, RooDataSet *dsetdestination, int numvar){

  if (!(*dset)) return;

  TH1F *hnum = new TH1F("hnum","hnum",n_etabins_forreweighting,etabins_forreweighting);
 TH1F *hden = new TH1F("hden","hden",n_etabins_forreweighting,etabins_forreweighting);
//  TH1F *hnum = new TH1F("hnum","hnum",25,0.,2.5);
//  TH1F *hden = new TH1F("hden","hden",25,0.,2.5);
  hnum->Sumw2();
  hden->Sumw2();

  const char* etaname=Form("rooeta%d",numvar);

  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill(fabs((*dset)->get(i)->getRealValue(etaname)),(*dset)->store()->weight(i));
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(fabs(dsetdestination->get(i)->getRealValue(etaname)),dsetdestination->store()->weight(i));
  }


  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

 // TH1F *h_data = hnum;
  hnum->Divide(hden);
  TH1F *h = hnum;

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_etarew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float eta = args.getRealValue(etaname);
    float neww = oldw*h->GetBinContent(h->FindBin(fabs(eta)));
  
  	//    std::cout << oldw << " " << neww << std::endl;
  	 //   std::cout << h->FindBin(fabs(eta)) << " " << fabs(eta) << std::endl;
  
    weight_eta->Fill(neww);	
	weight_etao->Fill(oldw);
	weight_etan->Fill(h->FindBin(fabs(eta)),neww);	
    weight_eta2o->Fill(h->FindBin(fabs(eta)),oldw);	
   if(oldw!=0 && neww!=0)weight_etanr->Fill(h->FindBin(fabs(eta)),oldw/neww);
  else {weight_etanr->Fill(-10,1);}
   // weight_pt2->Fill(pt,neww/oldw);
   if(oldw!=0 && neww!=0)weight_eta2->Fill(fabs(eta),oldw/neww);
  else {weight_eta2->Fill(-10,1);}
  	newdset->add(args,neww);
  
 //   weight_eta->Fill(neww);	
  //  weight_eta->Fill(eta);	
   // weight_eta2->Fill(eta,neww/oldw);
    //weight_eta2->Fill(neww, args.getRealValue("roovar1"));	
  }


  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

//  weight_eta->Scale(1.0/weight_eta->Integral());
//  weight_eta->Divide(h_data);
  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Eta 1d rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;

  delete old_dset;

};
void ratioplot(RooDataSet *data_axis1, RooDataSet* truth_axis1,Bool_t logplot){
   assert(data_axis1);
   assert(truth_axis1);
 TH1F::SetDefaultSumw2(kTRUE);
  TH1F* EB_data= new TH1F("EB_data","EB_data",100, 0.,9.);
  for(int i=0;i<(*data_axis1).numEntries();i++) {
     EB_data->Fill(fabs((*data_axis1).get(i)->getRealValue("roovar1")),(*data_axis1).store()->weight(i));
        }   
  TH1F* EB_truth= new TH1F("EB_truth","EB_truth",100, 0.,9.);
  for(int i=0;i<(*truth_axis1).numEntries();i++) {
EB_truth->Fill(fabs((*truth_axis1).get(i)->getRealValue("roovar1")),(*truth_axis1).store()->weight(i));
        } 
TCanvas* c5= new TCanvas("c5","c5"); 
c5->cd();
EB_data->SetLineColor(kRed);EB_data->SetMarkerSize(20);
EB_truth->SetLineColor(kBlue);EB_truth->SetMarkerSize(20);
EB_data->Draw();
// ratio plots
        TH1F *r1 = (TH1F*)EB_truth->Clone("r1");
        r1->Divide(EB_data);
     //   r1->SetMarkerStyle(20);
        r1->SetTitle("");
        TCanvas* c6= new TCanvas("c6","c6");
        if(logplot){
        c6->SetLogy();
        }
 TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
   pad1->SetBottomMargin(0); // Upper and lower plot are joined
   pad1->SetGridx();         // Vertical grid
   pad1->Draw();             // Draw the upper pad: pad1
   pad1->cd();
if(logplot){pad1->SetLogy();}            // pad1 becomes the current pad
    r1->GetXaxis()->SetTitle("Charged Iso (GeV)");
 EB_data->SetTitle("MC prompt photon comparison  reweighted");
//       h1->GetXaxis()->SetRangeUser(-100.,1000.);
    //  EB_data->GetYaxis()->SetRangeUser(10,100000);
   //   h1->SetMarkerColor(kBlack); h1->SetMarkerStyle(20);h1->SetLineColor(kBlack);
   //   h2->SetMarkerColor(kRed+1); h2->SetMarkerStyle(20);h2->SetLineColor(kRed+1);
// get mean
//Double_t err1;
//Double_t err2;
//      Double_t meanh1=EB_data->IntegralAndError(0,EB_data->GetNbinsX(),err1);
//      Double_t meanh2=EB_truth->IntegralAndError(0,EB_truth->GetNbinsX(),err2);
Double_t meanh1=EB_data->Integral(0,EB_data->GetNbinsX());
Double_t meanh2=EB_truth->Integral(0,EB_truth->GetNbinsX());  
    Double_t ratio=meanh2/meanh1;
//      Double_t err_ratio= ratio*ratio*pow(((err1/meanh1)*(err1/meanh1)+(err2/meanh2)*(err2/meanh2)),0.5);
     // Double_t ratio_error=
      cout << "data " << EB_data->GetNbinsX() <<" truth " << EB_truth->GetNbinsX() << endl;
 //     cout << "meanh1 " << meanh1 << " +/- " << err1 << "meanh2 " << meanh2 << " +/- " <<err2<<"ratio " << ratio << "+/-"<< err_ratio << endl;
  cout << " meanh1 " << meanh1  << " meanh2 " << meanh2<<" ratio " << ratio  << endl;  
    EB_data->Draw();
      EB_truth->Draw("SAME");
   TLegend* leg = new TLegend(0.55, 0.65, .9, .9);
       leg->SetFillColor(0);
       leg->AddEntry(EB_data," data " ,"pl");
       leg->AddEntry(EB_truth,"p  photon", "pl");

       leg->Draw();
        TLatex a;
        a.SetNDC();
//const char * ratio_out= Form("ratio of truth MC to data MC: %f +/-%f",ratio,err_ratio);
const char * ratio_out= Form("ratio of truth MC to data MC: %f",ratio);  
      a.DrawLatex(0.55,0.6,ratio_out);
   c6->cd();          // Go back to the main canvas before defining pad2
   TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  // pad2->SetLogy();
   gStyle->SetOptStat(0);
   pad2->SetTopMargin(0);
   pad2->SetBottomMargin(0.4);
   pad2->SetTicky();
   pad2->SetGridx(); // vertical grid
   pad2->Draw();
   pad2->cd();
   TLine *r = new TLine(0.,ratio,9.,ratio);
   r->SetLineColor(kRed);
   r1->Draw("ep");
   r->Draw("SAME");
  // r1->SetMinimum(0.0);  // Define Y ..
 //  r1->SetMaximum(10.); // .. range
   r1->GetYaxis()->SetTitle("MCrew/MC");
   r1->GetYaxis()->SetNdivisions(505);
   r1->GetYaxis()->SetTitleSize(18);
   r1->GetYaxis()->SetTitleFont(43);
   r1->GetYaxis()->SetTitleOffset(1.05);
   r1->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   r1->GetYaxis()->SetLabelSize(15);
 // X axis ratio plot settings
   r1->GetXaxis()->SetTitleSize(20);
   r1->GetXaxis()->SetTitleFont(43);
   r1->GetXaxis()->SetTitleOffset(4.);
   r1->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   r1->GetXaxis()->SetLabelSize(15);
  // c2->SaveAs("../plots/d1fit/mc_truth_EBEB.png");
 //  c2->SaveAs("../plots/d1fit/mc_truth_EBEB.root");
 //  c2->SaveAs("../plots/d1fit/mc_truth_EBEB.pdf");

}
/////////////////////////////////////////////////////////////////////
int do1dfit(){
//1d fit following dolightcomp from Marcos code
 //taking isolation variable
	  gStyle->SetOptStat(111111);
  Float_t leftrange=0. ;
  Float_t rightrange=9.;
  TH1F::SetDefaultSumw2(kTRUE);
  roovar1 = new RooRealVar("roovar1","roovar1",leftrange,rightrange);
  roovar2 = new RooRealVar("roovar2","roovar2",leftrange,rightrange);
  roovar1->setRange(leftrange,rightrange);
  roovar2->setRange(leftrange,rightrange);
   RooBinning tbins (0,9); 
   tbins.addUniform(90,0,9) ;
//   tbins.addUniform(5,0,1) ;
//   tbins.addUniform(30,1,3);
//     tbins.addUniform(15,3,9) ;
//     tbins.addUniform(10,0.,1.) ;
  //   tbins.addUniform(20,1.,3.);
   //    tbins.addUniform(10,3.,6.) ;
   //    tbins.addUniform(3,6.,9.) ;
         roovar1->setBinning(tbins);
  roovar1->SetTitle("Iso_{1}");
  roovar2->SetTitle("Iso_{2}"); 
  rooisopv1 = new RooRealVar("rooisopv1","rooisopv1",leftrange,rightrange);
  rooisopv2 = new RooRealVar("rooisopv2","rooisopv2",leftrange,rightrange);
  rooisopv1->setRange(leftrange,rightrange);
  rooisopv2->setRange(leftrange,rightrange);
  rooisopv1->setBinning(tbins);
  rooisopv1->SetTitle("IsoVzero_{1}");
  rooisopv2->SetTitle("IsoVzero_{2}"); 
  rooisowv1 = new RooRealVar("rooisowv1","rooisowv1",leftrange,rightrange);
  rooisowv2 = new RooRealVar("rooisowv2","rooisowv2",leftrange,rightrange);
  rooisowv1->setRange(leftrange,rightrange);
  rooisowv2->setRange(leftrange,rightrange);
  rooisowv1->setBinning(tbins);
  rooisowv1->SetTitle("IsoWV_{1}");
  rooisowv2->SetTitle("IsoWV_{2}"); 
  rooeta1 = new RooRealVar("rooeta1","rooeta1",0,2.5);
  rooeta2 = new RooRealVar("rooeta2","rooeta2",0,2.5);
  //roopt1 = new RooRealVar("roopt1","roopt1",100.);
 // roopt2 = new RooRealVar("roopt2","roopt2",100.);


  Float_t ptleftrange=0. ;
  Float_t ptrightrange=1000.;
//  Float_t ptrightrange=400.;
  RooBinning ptbins (ptleftrange,ptrightrange); 
  ptbins.addUniform(40,20,100) ;
  ptbins.addUniform(40,100,200) ;
  ptbins.addUniform(10,200,400) ;
  ptbins.addUniform(1,400,1000) ;
  Float_t dptleftrange=0. ;
  Float_t dptrightrange=1000.;
  RooBinning dptbins (dptleftrange,dptrightrange); 
  dptbins.addUniform(80,0,80) ;
  dptbins.addUniform(40,80,200) ;
  dptbins.addUniform(20,200,400) ;
  dptbins.addUniform(1,400,2000) ;


  roopt1 = new RooRealVar("roopt1","roopt1",ptleftrange,ptrightrange);
  roopt2 = new RooRealVar("roopt2","roopt2",ptleftrange,ptrightrange);
  roopt1->setRange(ptleftrange,ptrightrange);
  roopt2->setRange(ptleftrange,ptrightrange);
  roopt1->setBinning(ptbins);
  roopt1->SetTitle("pt_{1}");
  roodiphopt = new RooRealVar("roodiphopt","roodiphopt",dptleftrange,dptrightrange);
  roodiphopt->setRange(dptleftrange,dptrightrange);
  roodiphopt->setBinning(dptbins);
  roodiphopt->SetTitle("dipho_{1}");
  
  roosieie1 = new RooRealVar("roosieie1","roosieie1",0,0.045);
  roosieie2 = new RooRealVar("roosieie2","roosieie2",0,0.045);
  //roodiphopt = new RooRealVar("roodiphopt","roodiphopt",100.);
  //roorho = new RooRealVar("roorho","roorho",20.);
  //roosigma = new RooRealVar("roosigma","roosigma",2.);
  roorho = new RooRealVar("roorho","roorho",0.,100.);
  roosigma = new RooRealVar("roosigma","roosigma",0.,20.);
  roonvtx = new RooRealVar("roonvtx","roonvtx",20.);
  rooweight = new RooRealVar("rooweight","rooweight",1.);
  assert (roovar1);
  assert (roovar1);
  assert (rooisopv1);
  assert (rooisopv2);
  assert (rooisowv1);
  assert (rooisowv2);
  assert (roopt1);
  assert (roosieie1);
  assert (rooeta1);
  assert (roopt2);
  assert (roosieie2);
  assert (rooeta2);
  assert (roodiphopt);
  assert (roorho);
  assert (roosigma);
  assert (roonvtx);
  assert (rooweight);
  TH1::SetDefaultSumw2(kTRUE);
  //Axis 1 and 2 randomly swapped
TFile* fdata=NULL; TFile* frcone=NULL; TFile * fsideband=NULL; //TFile* testfile=NULL;
 	 //fdata_std = new TFile("../forroofit/fullmc_fulletarange/mc2dstd_fulletarange.root","read"); 
 	 fdata = new TFile("../forroofit/Feb13/data.root","read"); 
 	 get_roodset_from_ttree(fdata,"for_roofit",dataset_data);
//both leading and subleading photon 
 	 //rcone = new TFile("../forroofit/fullmc_fulletarange/sig_fulletarange.root","read"); 
 //	frcone = new TFile("../forroofit/Feb13/rcone.root","read"); 
   	frcone = new TFile("../forroofit/Feb13/sig.root","read"); 
	 get_roodset_from_ttree(frcone,"for_roofit",dataset_rcone);
	 //sideband= new TFile("../forroofit/fullmc_fulletarange/bkg_fulletarange.root","read");
//	 fsideband= new TFile("../forroofit/Feb13/sideband.root","read");
	 fsideband= new TFile("../forroofit/Feb13/bkg.root","read");
	 get_roodset_from_ttree(fsideband,"for_roofit",dataset_sideband);


  assert(fdata); assert(frcone);  assert(fsideband);
//define mass bins
  TH1F* diphomass_data=(TH1F*)fdata->Get("h_diphomass") ;
  assert(diphomass_data);
  quantiles(diphomass_data);
//mass->Scale(1.0/h_diphomass->Integral())
//get_roodset_from_ttree(testfile,"for_roofit",dataset_test);
//cout << "before cuts " << endl;
//  TString cut_string= Form("roovar1<%f",rightrange-1e-5);
    // TString cut_string= Form("roovar1<%f && rooeta1<1.5 && rooeta2 < 1.5",rightrange-1e-5);
  // TString cut_string= Form("roovar1<%f && roowisov1<%f && rooisopv1<%f",rightrange-1e-5,rightrange-1e-5,rightrange-1e-5);
  // TString cut_string= Form("roovar1<%f && roowisov1<%f && rooisopv1<%f",rightrange-1e-5,rightrange-1e-5,rightrange-1e-5);
   TString cut_string= Form("roovar1<%f && rooeta1<1.5 && rooeta2 < 1.5",rightrange-1e-5);
  if (dataset_data) dset_data = (RooDataSet*)(dataset_data->reduce(Name("dset_data"),Cut(cut_string)));
  if (dataset_rcone) dset_rcone = (RooDataSet*)(dataset_rcone->reduce(Name("dset_rcone"),Cut(cut_string)));
  if (dataset_sideband) dset_sideb = (RooDataSet*)(dataset_sideband->reduce(Name("dset_sideb"),Cut(cut_string)));
////create subset of dataset e.g. for 1d signal template, here leading photon 
 
 
 
  //RooDataSet *dset_data_axis1 = (RooDataSet*)(dset_data->reduce(Name("dset_data_axis1"),SelectVars(RooArgList(*roovar1,*rooisopv1,*rooisowv1,*roopt1,*roosieie1,*roodiphopt,*rooeta1,*roorho,*roosigma,*roonvtx))));
 // RooDataSet *dset_rcone_axis1 = (RooDataSet*)(dset_rcone->reduce(Name("dset_rcone_axis1"),SelectVars(RooArgList(*roovar1,*rooisopv1, *rooisowv1,*roopt1,*roosieie1,*rooeta1,*roodiphopt,*roorho,*roosigma,*roonvtx))));
 // RooDataSet *dset_sideb_axis1 = (RooDataSet*)(dset_sideb->reduce(Name("dset_sideb_axis1"),SelectVars(RooArgList(*roovar1,*rooisopv1,*rooisowv1,*roopt1,*roosieie1,*rooeta1,*roorho,*roodiphopt,*roosigma,*roonvtx))))
  RooDataSet *dset_data_axis1 = (RooDataSet*)(dset_data->reduce(Name("dset_data_axis1"),SelectVars(RooArgList(*roovar1,*rooisopv1,*rooisowv1,*roopt1,*roosieie1,*roodiphopt,*rooeta1,*roorho,*roosigma))));
  RooDataSet *dset_rcone_axis1 = (RooDataSet*)(dset_rcone->reduce(Name("dset_rcone_axis1"),SelectVars(RooArgList(*roovar1,*rooisopv1, *rooisowv1,*roopt1,*roosieie1,*rooeta1,*roodiphopt,*roorho,*roosigma))));
  RooDataSet *dset_sideb_axis1 = (RooDataSet*)(dset_sideb->reduce(Name("dset_sideb_axis1"),SelectVars(RooArgList(*roovar1,*rooisopv1,*rooisowv1,*roopt1,*roosieie1,*rooeta1,*roorho,*roodiphopt,*roosigma))));
  assert(dset_data_axis1);assert(dset_rcone_axis1);assert(dset_sideb_axis1); 

  cout << "ratio Entries truth/data weighted " << dset_rcone_axis1->sumEntries()/ dset_data_axis1->sumEntries()  << endl;
  //reweight always reweight to data in signal region
  reweight_rhosigma(w_rhor,w_rhoro,w_rhorn,w_rhor2o,w_rhornr,w_rhor2,w_sigmarn,w_sigmar2o,w_sigmarnr,w_sigmar2,&dset_sideb_axis1,dset_data_axis1); 
  reweight_rhosigma(w_rho,w_rhoo,w_rhon,w_rho2o,w_rhonr,w_rho2,w_sigman,w_sigma2o,w_sigmanr,w_sigma2,&dset_sideb_axis1,dset_data_axis1); 
  reweight_eta_1d(w_etar,w_etaro,w_etarn,w_etar2o,w_etarnr,w_etar2,&dset_rcone_axis1,dset_data_axis1,1);
  reweight_eta_1d(w_eta,w_etao,w_etan,w_eta2o,w_etanr,w_eta2,&dset_sideb_axis1,dset_data_axis1,1);
  reweight_pt_1d(w_pt,w_pto,w_ptn,w_pt2o,w_ptnr,w_pt2,&dset_sideb_axis1,dset_data_axis1,1);
  reweight_pt_1d(w_ptr,w_ptro,w_ptrn,w_ptr2o,w_ptrnr,w_ptr2,&dset_rcone_axis1,dset_data_axis1,1);
  reweight_diphotonpt_1d(w_diphopt,w_diphopto,w_diphoptn,w_diphopt2o,w_diphoptnr,w_diphopt2,&dset_sideb_axis1,dset_data_axis1); 
  reweight_diphotonpt_1d(w_diphopt,w_diphopto,w_diphoptn,w_diphopt2o,w_diphoptnr,w_diphopt2,&dset_sideb_axis1,dset_data_axis1); 
  

w_pt->SetTitle("new weights ");
w_ptn->SetTitle("new weights vs pt bins");
w_pto->SetTitle("old weights");
w_ptnr->SetTitle("old weights/new weights vs pt bins");
w_pt2o->SetTitle("old weights vs pt bins");
w_pt2->SetTitle("new weights vs pt");
w_diphoptn->SetTitle("new weights vs diphopt bins");
w_diphopt2o->SetTitle("old weights vs diphopt bins");
w_diphoptnr->SetTitle("old weights/new weights vs diphopt bins");
w_diphopt->SetTitle("new weights");
w_diphopt2->SetTitle("new weights vs diphopt");
w_diphopto->SetTitle("old weights");
w_etan->SetTitle("new weights vs eta bins");
w_eta2o->SetTitle("old weights vs eta bins");
w_etanr->SetTitle("old weights/new weights vs eta bins");
w_eta->SetTitle("new weights");
w_eta2->SetTitle("new weights vs eta");
w_etao->SetTitle("old weights");
w_rhon->SetTitle("new weights vs rho bins");
w_rho2o->SetTitle("old weights vs rho bins");
w_rhonr->SetTitle("old weights/new weights vs rho bins");
w_rho->SetTitle("new weights");
w_rho2->SetTitle("new weights vs rho");
w_rhoo->SetTitle("old weights");
w_sigman->SetTitle("new weights vs sigma bins");
w_sigma2o->SetTitle("old weights vs sigma bins");
w_sigmanr->SetTitle("old weights/new weights vs sigma bins");
	TCanvas * c_reweights=new TCanvas("c_reweights","c_reweights");
  c_reweights->Divide(3,2);  c_reweights->cd(1); w_ptn->Draw("colz");
  c_reweights->cd(2);  w_pt2o->Draw("colz");  c_reweights->cd(3);   w_ptnr->Draw("colz");
  c_reweights->cd(4);  w_pt2->Draw("colz"); c_reweights->cd(5);w_pt->Draw();
  c_reweights->cd(6);   w_pto->Draw(); c_reweights->SaveAs("../plots/c_reweights.root");c_reweights->SaveAs("../plots/c_reweights.png");
  
 TCanvas * c_reweights2=new TCanvas("c_reweights2","c_reweights2");
  c_reweights2->Divide(3,2); c_reweights2->cd(1);  w_diphoptn->Draw("colz");  c_reweights2->cd(2);    w_diphopt2o->Draw("colz");
  c_reweights2->cd(3);    w_diphoptnr->Draw("colz");c_reweights2->cd(4);  w_diphopt2->Draw("colz");
  c_reweights2->cd(5);w_diphopt->Draw();  c_reweights2->cd(6);    w_diphopto->Draw();
 c_reweights2->SaveAs("../plots/c_reweights2.root");
 c_reweights2->SaveAs("../plots/c_reweights2.png");

  TCanvas * c_rhoreweights=new TCanvas("c_rhoreweights","c_rhoreweights");
  c_rhoreweights->Divide(3,2);  c_rhoreweights->cd(1);  w_rhon->Draw("colz");
  c_rhoreweights->cd(2);    w_rho2o->Draw("colz");  c_rhoreweights->cd(3);    w_rhonr->Draw("colz");
  c_rhoreweights->cd(4);  w_rho2->Draw("colz");  c_rhoreweights->cd(5);
  w_rho->Draw();c_rhoreweights->cd(6);    w_rhoo->Draw(); c_rhoreweights->SaveAs("../plots/c_rhoreweights.root");
  TCanvas * c_sigmareweights=new TCanvas("c_sigmareweights","c_sigmareweights");
  c_sigmareweights->Divide(2,2);  c_sigmareweights->cd(1);  w_sigman->Draw("colz");  c_sigmareweights->cd(2); 
   w_sigma2o->Draw("colz");  c_sigmareweights->cd(3);    w_sigmanr->Draw("colz");  c_sigmareweights->cd(4);
  w_sigma2->Draw("colz");   c_sigmareweights->SaveAs("../plots/c_sigmareweights.root");
  
  TCanvas * c_etareweights=new TCanvas("c_etareweights","c_etareweights");
  c_etareweights->Divide(3,2);  c_etareweights->cd(1);  w_etan->Draw("colz");
  c_etareweights->cd(2);    w_eta2o->Draw("colz");  c_etareweights->cd(3);    w_etanr->Draw("colz");
  c_etareweights->cd(4);  w_eta2->Draw("colz");  c_etareweights->cd(5);
  w_eta->Draw();c_etareweights->cd(6);   w_etao->Draw();  c_etareweights->SaveAs("../plots/c_etareweights.root");
 
 /*
  cout << "after reweighting" << endl;
	  cout << "sig " << dset_rcone_axis1->sumEntries() << endl;
  cout << "bkg " << dset_sideb_axis1->sumEntries() << endl;
  cout << "data " << dset_data_axis1->sumEntries() << endl;
  cout << "signal+ fakes " << dset_rcone_axis1->sumEntries()+ dset_sideb_axis1->sumEntries()<< endl;
  cout << "one? " << (dset_rcone_axis1->sumEntries()+ dset_sideb_axis1->sumEntries())/dset_data_axis1->sumEntries()<< endl;
  cout << "ratio Entries truth/data weighted " << dset_rcone_axis1->sumEntries()/dset_data_axis1->sumEntries()  << endl;
*/

// ratioplot(dset_data_axis1, dset_rcone_axis1,kTRUE);
  //fitting
 cout << "fitting starts" << endl;
//first step determine single-photon purities on Iso1 and Iso2 axes
//for integral 0-90 GeV, fit should not depend on initial values!  
//give initial values for purity fraction -j1 overall purity for leg 1 
  RooRealVar *j1=NULL;RooFormulaVar *fsig1=NULL;
  j1 = new RooRealVar("j1","j1",0.3,0,1);
  fsig1 = new RooFormulaVar("fsig1","fsig1","j1",RooArgList(*j1));
  
  //produces binned RooDataHist also when RooRealVar input and weights are unbinned
  RooArgSet var;
  if(dohggv){var.add(*roovar1);}
  else if(dovzero) {var.add(*rooisopv1);}
  else if(dowv) {var.add(*rooisowv1);}
  RooDataHist *drconehist_axis1 = new RooDataHist("drconehist_axis1","drconehist_axis1",var);
  drconehist_axis1->add(*dset_rcone_axis1);
  drconehist_axis1->Print("V");
  RooDataHist *dsidebhist_axis1 = new RooDataHist("dsidebhist_axis1","dsidebhist_axis1",var);
  dsidebhist_axis1->add(*dset_sideb_axis1);
  dsidebhist_axis1->Print("V"); 
  RooHistPdf *drconepdf_axis1 = new RooHistPdf("drconepdf_axis1","drconepdf_axis1",var,*drconehist_axis1);
  RooHistPdf *dsidebpdf_axis1 = new RooHistPdf("dsidebpdf_axis1","dsidebpdf_axis1",var,*dsidebhist_axis1);
  //normalizes histograms over all observables
  //RooHistPdf *drconepdf_axis1 = new RooHistPdf("drconepdf_axis1","drconepdf_axis1",RooArgList(*roovar1),*drconehist_axis1);
  //RooHistPdf *dsidebpdf_axis1 = new RooHistPdf("dsidebpdf_axis1","dsidebpdf_axis1",RooArgList(*roovar1),*dsidebhist_axis1);
  
  /*
 //debugging reweighting 
  RooDataHist *test_dpt = new RooDataHist("test_dpt","test_dpt",RooArgSet(*roodiphopt));
  test_dpt->add(*dset_data_axis1);
  RooDataHist *test_spt = new RooDataHist("test_spt","test_spt",RooArgSet(*roodiphopt));
  test_spt->add(*dset_sideb_axis1);
  RooHistPdf *test_datapt = new RooHistPdf("test_datapt","test_datapt",RooArgList(*roodiphopt),*test_dpt);
  RooHistPdf *test_sidept = new RooHistPdf("test_sidept","test_sidept",RooArgList(*roodiphopt),*test_spt);
  TCanvas *c3=new TCanvas("c3","diphopt after reweighting rho, eta , single pt");
  c3->cd();
  RooPlot *testpt = roodiphopt->frame(Title("diphopt after reweighting rho, eta , single pt"));
  test_datapt->plotOn(testpt,LineStyle(kDashed),LineColor(kRed),Name("data"));
  test_sidept->plotOn(testpt,LineStyle(kDashed),LineColor(kBlack),Name("sideband"));
  testpt->Draw(); 
 */
/*  TCanvas *c8=new TCanvas("c8","c8");
  c8->cd();
  RooPlot *frame1bla2 = roovar1->frame(Title("roohistpdfs"));
  drconehist_axis1->plotOn(frame1bla2,LineStyle(kDashed),MarkerColor(kRed),Name("sig"));
  dsidebhist_axis1->plotOn(frame1bla2,LineStyle(kDashed),MarkerColor(kBlack),Name("bkg"));
  frame1bla2->Draw();
*/ 


 
 // RooAddPdf is an efficient implementation of a sum of PDFs of the form  c_1*PDF_1 + c_2*PDF_2 + ... (1-sum(c_1...c_n-1))*PDF_n
  //the sum of the coefficients is enforced to be one,and the coefficient of the last PDF is calculated from that condition.
  //ArgList of coefficients -1 of roohistpdf
  RooAddPdf *model_axis1=NULL;
  model_axis1 = new RooAddPdf("model_axis1","model_axis1",RooArgList(*drconepdf_axis1,*dsidebpdf_axis1),RooArgList(*fsig1),kFALSE);
  //if does not converge take NumCPU strategy 2
  RooFitResult *firstpass1;
 //hesse run by default
 // firstpass1 = model_axis1->fitTo(*test, NumCPU(8), Extended(false),SumW2Error(kTRUE),Save(kTRUE));
  firstpass1 = model_axis1->fitTo(*dset_data_axis1, NumCPU(8), Extended(false),SumW2Error(kTRUE),Save(kTRUE));
  TFile resFile("mctruth_nooverflow_reweight.root","RECREATE"); 
  firstpass1->Write();
  resFile.Write();  
  resFile.Close();
  firstpass1->Print();

  TLatex b;b.SetNDC();b.SetTextSize(0.06);b.SetTextColor(kRed);
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  TLegend *leg = new TLegend(0.15,0.8,0.35,0.9);
  c1->cd(1);
  TString title; 
  title= "1d fit for hgg vertex, EBEB, reweighted";

  TString titlevzero; 
  titlevzero= "1d fit for vertex[0], full eta range, reweighted";
  TString titlewv; 
  titlewv= "1d fit for worst isolation, full eta range, reweighted";
  //lin plot
  RooPlot *frame1bla;
  if(dohggv) {frame1bla = roovar1->frame(Title(title.Data()));}
  else  if(dovzero) {frame1bla = rooisopv1->frame(Title(titlevzero.Data()));}
  else  if(dowv) {frame1bla = rooisowv1->frame(Title(titlewv.Data()));}
 
 
  dset_data_axis1->plotOn(frame1bla,Binning(tbins),Name("data"));
  // dset_data_axis1->plotOn(frame1bla,Binning(tbins),Name("data"));
  model_axis1->plotOn(frame1bla,Name("fit"));
  model_axis1->plotOn(frame1bla,Components("drconepdf_axis1"),LineStyle(kDashed),LineColor(kRed),Name("signal"));
  model_axis1->plotOn(frame1bla,Components("dsidebpdf_axis1"),LineStyle(kDashed),LineColor(kBlack),Name("background"));
  frame1bla->Draw();
  leg->AddEntry("fit","fit","l");
//  leg->AddEntry("signal","rcone MC","l");
//  leg->AddEntry("background","sideband  MC","l");
  leg->AddEntry("signal","signal MC","l");
  leg->AddEntry("background","background  MC","l");
  leg->SetFillColor(kWhite); 
  leg->Draw();
  b.DrawLatex(0.55,0.7,"PRELIMINARY");

//log plot
  TCanvas *c2 = new TCanvas("c2","c2",1200,800);
  c2->cd(1);
  gPad->SetLogy();
  RooPlot *frame1logbla;
  if(dohggv) {frame1logbla = roovar1->frame(Title(title.Data()));}
  else if(dovzero) {frame1logbla = rooisopv1->frame(Title(titlevzero.Data()));}
  else if(dowv) {frame1logbla = rooisowv1->frame(Title(titlewv.Data()));}
  dset_data_axis1->plotOn(frame1logbla,Binning(tbins),Name("data"));
  model_axis1->plotOn(frame1logbla,Name("fit"));
  model_axis1->plotOn(frame1logbla,Components("drconepdf_axis1"),LineStyle(kDashed),LineColor(kRed),Name("signal"));
  model_axis1->plotOn(frame1logbla,Components("dsidebpdf_axis1"),LineStyle(kDashed),LineColor(kBlack),Name("background"));
  frame1logbla->Draw();
  leg->Draw();
  b.DrawLatex(0.55,0.7,"PRELIMINARY");
  model_axis1->Print();
  dset_data_axis1->Print();

  //save

  TString titleout;

//  if(dohggv) titleout= "EBEB_reweight_rcone_sideband_hggv";
  if(dohggv) titleout= "EBEB_reweight_truth_hggv";
  else if(dovzero) titleout= "nooverflow_reweight_fulletarange_rcone_sideband_vzero";
  else if(dowv) titleout= "nooverflow_reweight_fulletarange_rcone_sideband_wv";
  //titleout="rcone_sideb";
  const char* outfile=Form("../plots/mc_%s", titleout.Data());
  c1->SaveAs(Form("%s_lin.png",outfile));c1->SaveAs(Form("%s_lin.root",outfile));c1->SaveAs(Form("%s_lin.pdf",outfile));
  c2->SaveAs(Form("%s_log.png",outfile));c2->SaveAs(Form("%s_log.root",outfile));c2->SaveAs(Form("%s_log.pdf",outfile));
return 0;
  }

//////////////////2d fit

int do2dfit(){
/*
	gStyle->SetOptStat(111111);
  Float_t leftrange=0. ;
  Float_t rightrange=9.;
  Int_t roovarbins=90;
  TH1F::SetDefaultSumw2(kTRUE);
  roovar1 = new RooRealVar("roovar1","roovar1",leftrange,rightrange);
  roovar2 = new RooRealVar("roovar2","roovar2",leftrange,rightrange);
  roovar1->setRange(leftrange,rightrange);
  roovar2->setRange(leftrange,rightrange);
  roovar1->setBins(roovarbins);
  roovar2->setBins(roovarbins);
  roovar1->SetTitle("Iso_{1}");
  roovar2->SetTitle("Iso_{2}"); 
  rooeta1 = new RooRealVar("rooeta1","rooeta1",0,2.5);
  rooeta2 = new RooRealVar("rooeta2","rooeta2",0,2.5);
  roopt1 = new RooRealVar("roopt1","roopt1",25,1000);
  roopt2 = new RooRealVar("roopt2","roopt2",25,1000);
  roosieie1 = new RooRealVar("roosieie1","roosieie1",0,0.045);
  roosieie2 = new RooRealVar("roosieie2","roosieie2",0,0.045);
 roodiphopt = new RooRealVar("roodiphopt","roodiphopt",200.);
  roorho = new RooRealVar("roorho","roorho",20.);
  roosigma = new RooRealVar("roosigma","roosigma",1.);
  roonvtx = new RooRealVar("roonvtx","roonvtx",20.);
  rooweight = new RooRealVar("rooweight","rooweight",1.);
  assert (roovar1);
  assert (roovar1);
  assert (roopt1);
  assert (roosieie1);
  assert (rooeta1);
  assert (roopt2);
assert (roodiphopt);
  assert (roosieie2);
  assert (rooeta2);
  assert (roorho);
  assert (roosigma);
  assert (roonvtx);
  assert (rooweight);
  TH1::SetDefaultSumw2(kTRUE);
  //Axis 1 and 2 randomly swapped
  //TFile *fdata_std = new TFile("forroofit/2dstd_data_EBEB_forroofit.root","read");
  TFile *fdata_std = new TFile("forroofit/2dstd_swaped_allmc_EBEB_forroofit.root","read");
  get_roodset_from_ttree(fdata_std,"for_roofit",dataset_data);
//both leading and subleading photon 
// TFile *rcone = new TFile("forroofit/1drcone_template.root","read");
 TFile *rcone = new TFile("forroofit/2drcone_swaped_allmc_EBEB_forroofit.root","read");
  get_roodset_from_ttree(rcone,"for_roofit",dataset_rcone);
//first fake then true photon
//  TFile *fdatasideb_b= new TFile("forroofit/side_signalband_data_EBEB_forroofit.root","read");
TFile *fdatasideb_b= new TFile("forroofit/2drconesideband_swaped_allmc_EBEB_forroofit.root","read");
  //signal+ sideband region first variable signal, second sideband
// TFile *fdatasideb_b= new TFile("forroofit/1legsideband_swaped_allmc_EBEB_forroofit.root","read");
//  TFile *fdatasideb_b= new TFile("forroofit/1dtemp_bkgtruth.root","read");
 get_roodset_from_ttree(fdatasideb_b,"for_roofit",dataset_datasideb_b);
//both in sideband, randomly swapped
  TFile *sideband= new TFile("forroofit/2dsideband_swaped_allmc_EBEB_forroofit.root","read");
//TFile *sideband= new TFile("forroofit/1dsideb_template.root","read");
  get_roodset_from_ttree(sideband,"for_roofit",dataset_sideband);

  assert(fdata_std); assert(rcone); assert(fdatasideb_b);  assert(sideband);
//cut away overflow bin
  if (dataset_data) dset_data1 = (RooDataSet*)(dataset_data->reduce(Name("dset_data1"),Cut(Form("roovar1<%f",rightrange-1e-5))));
  dset_data = (RooDataSet*)(dset_data1->reduce(Name("dset_data"),Cut(Form("roovar2<%f",rightrange-1e-5))));
  if (dataset_rcone) dset_rcone1 = (RooDataSet*)(dataset_rcone->reduce(Name("dset_rcone1"),Cut(Form("roovar2<%f",rightrange-1e-5))));
  dset_rcone = (RooDataSet*)(dset_rcone1->reduce(Name("dset_rcone"),Cut(Form("roovar1<%f",rightrange-1e-5))));
if (dataset_datasideb_b) dset_sideb1 = (RooDataSet*)(dataset_datasideb_b->reduce(Name("dset_sideb1"),Cut(Form("roovar1<%f",rightrange-1e-5))));
  dset_sideb = (RooDataSet*)(dset_sideb1->reduce(Name("dset_sideb"),Cut(Form("roovar2<%f",rightrange-1e-5))));
  if (dataset_sideband) dset_sideb2s1 = (RooDataSet*)(dataset_sideband->reduce(Name("dset_sideb2s1"),Cut(Form("roovar1<%f",rightrange-1e-5))));
  dset_sideb2s = (RooDataSet*)(dset_sideb2s1->reduce(Name("dset_sideb2s"),Cut(Form("roovar2<%f",rightrange-1e-5))));

    const float pp_init = 0.30;
        const float pf_init = 0.23;//0.33
//give initial values for purity fraction -j1 overall purity for leg 1 
  RooRealVar *pp = new RooRealVar("pp","pp",pp_init,0,1);

   RooRealVar *j1 = new RooRealVar("j1","j1",pp_init+pf_init,0,1);

   //  RooRealVar *j1_pf = new RooRealVar("j1_pf","j1_pf",pf_init,0,1);
//  RooRealVar *j2 = new RooRealVar("j2","j2",pp_init+pf_init,0,1);
//  RooRealVar *j2_pf = new RooRealVar("j2_pf","j2_pf",pf_init,0,1);
  RooFormulaVar *fsig1 = new RooFormulaVar("fsig1","fsig1","j1",RooArgList(*j1));
  //RooFormulaVar *fsig2 = new RooFormulaVar("fsig2","fsig2","j2",RooArgList(*j2));
  RooFormulaVar *fsigsig = new RooFormulaVar("fsigsig","fsigsig","pp",RooArgList(*pp));
  RooFormulaVar *fsigbkg = new RooFormulaVar("fsigbkg","fsigbkg","fsig1-pp",RooArgList(*fsig1,*pp));
RooFormulaVar *fbkgbkg = new RooFormulaVar("fbkgbkg","fbkgbkg","1-fsigsig-fsigbkg",RooArgList(*fsigsig,*fsigbkg));
 RooDataHist *sigsigdhist = new RooDataHist("sigsigdhist","sigsigdhist",RooArgList(*roovar1,*roovar2),*dset_rcone);
    RooDataHist *sigbkgdhist = new RooDataHist("sigbkgdhist","sigbkgdhist",RooArgList(*roovar1,*roovar2),*dset_sideb);
    RooDataHist *bkgbkgdhist = new RooDataHist("bkgbkgdhist","bkgbkgdhist",RooArgList(*roovar1,*roovar2),*dset_sideb2s);
  RooAbsPdf*    sigsigpdf = new RooHistPdf("sigsigpdf","sigsigpdf",RooArgList(*roovar1,*roovar2),*sigsigdhist);
  RooAbsPdf*    sigbkgpdf = new RooHistPdf("sigbkgpdf","sigbkgpdf",RooArgList(*roovar1,*roovar2),*sigbkgdhist);
   RooAbsPdf*   bkgbkgpdf = new RooHistPdf("bkgbkgpdf","bkgbkgpdf",RooArgList(*roovar1,*roovar2),*bkgbkgdhist);
// TFile *res_firstpass = new TFile("result1d.root","read");
      //    float lowerbounds[4]={0,f1p-1,f2p-1,f1p+f2p-1};
//    float upperbounds[4]={1,f1l,f2l,f1l+f2l};
    float lowerbounds[4]={0,5.2698e-01-1,5.2698e-01-1,5.2698e-01+5.2698e-01-1};
    float upperbounds[4]={1,5.2698e-01,5.2698e-01,5.2698e-01+5.2698e-01};

    float minpp = TMath::MaxElement(4,lowerbounds);
    float maxpp = TMath::MinElement(4,upperbounds);
    pp->setVal((minpp+maxpp)/2);
RooFitResult *secondpass;
 //   RooGaussian *constraint_gaussian_j1 = new RooGaussian("constraint_gaussian_j1","constraint_gaussian_j1",*j1,RooRealConstant::value(j1->getVal()),RooRealConstant::value(j1->getPropagatedError(*firstpass)));
RooGaussian *constraint_gaussian_j1 = new RooGaussian("constraint_gaussian_j1","constraint_gaussian_j1",*j1,RooRealConstant::value(j1->getVal()),RooRealConstant::value(j1->getError()));
RooArgSet *constraint_pdf_set = new RooArgSet(*constraint_gaussian_j1);
    RooArgSet *constraint_parameters_set = new RooArgSet(*j1);
RooConstraintSum *constraint_gaussian_nll = new RooConstraintSum("constraint_gaussian_nll","constraint_gaussian_nll",*constraint_pdf_set,*constraint_parameters_set);


    RooAddPdf *model_2D_uncorrelated = new RooAddPdf("model_2D_uncorrelated","model_2D_uncorrelated",RooArgList(*sigsigpdf,*sigbkgpdf,*bkgbkgpdf),RooArgList(*fsigsig,*fsigbkg,*fbkgbkg),kFALSE);
    RooNLLVar *model_2D_uncorrelated_noextended_nll = new RooNLLVar("model_2D_uncorrelated_noextended_nll","model_2D_uncorrelated_noextended_nll",*model_2D_uncorrelated,*dset_data,NumCPU(8));
    RooAddition *model_2D_uncorrelated_noextended_nll_constraint = new RooAddition("model_2D_uncorrelated_noextended_nll_constraint","model_2D_uncorrelated_noextended_nll_constraint",RooArgSet(*model_2D_uncorrelated_noextended_nll,*constraint_gaussian_nll));


    RooMinimizer *minuit_secondpass_constraint = new RooMinimizer(*model_2D_uncorrelated_noextended_nll_constraint);
    minuit_secondpass_constraint->migrad();
    minuit_secondpass_constraint->hesse();
    RooFitResult *secondpass_constraint;
    secondpass_constraint = minuit_secondpass_constraint->save("secondpass_constraint","secondpass_constraint");
    secondpass_constraint->Print();  
RooMinimizer *minuit_secondpass = new RooMinimizer(*model_2D_uncorrelated_noextended_nll);
    minuit_secondpass->migrad();
    minuit_secondpass->hesse();
secondpass = minuit_secondpass->save("secondpass","secondpass");
     secondpass->Print();
   cout << "pp " << fsigsig->getVal() << " " << fsigsig->getPropagatedError(*secondpass) << std::endl;
 cout << "pf " << fsigbkg->getVal() << " " << fsigbkg->getPropagatedError(*secondpass) << std::endl; 
     cout << "ff " << fbkgbkg->getVal() << " " << fbkgbkg->getPropagatedError(*secondpass) << std::endl;



      TCanvas *c4 = new TCanvas("c4","c4",1200,800);
      c4->Divide(2,2);
        c4->cd(1);
        RooPlot *frame1final = roovar1->frame(Title("Fit axis 1 - binned"));
        dset_data->plotOn(frame1final,Name("data"));
        model_2D_uncorrelated->plotOn(frame1final,Name("fit"));
        sigsigpdf->plotOn(frame1final,Normalization(fsigsig->getVal(),RooAbsPdf::Relative),Name("plot_sigsig_axis1"),LineStyle(kDashed),LineColor(kRed));
        sigbkgpdf->plotOn(frame1final,Normalization(fsigbkg->getVal(),RooAbsPdf::Relative),Name("plot_sigbkg_axis1"),LineStyle(kDashed),LineColor(kGreen));
        bkgbkgpdf->plotOn(frame1final,Normalization(fbkgbkg->getVal(),RooAbsPdf::Relative),Name("plot_bkgbkg_axis1"),LineStyle(kDashed),LineColor(kBlack));
        frame1final->Draw();
        //    c2->GetPad(1)->SetLogy(1);
        TLegend *leg = new TLegend(0.18,0.74,0.38,0.94);
        leg->AddEntry("data","data","lp");
        leg->AddEntry("fit","fit","l");
        leg->AddEntry("plot_sigsig_axis1","prompt-prompt","l");
        leg->AddEntry("plot_sigbkg_axis1","prompt-fake","l");
        leg->AddEntry("plot_bkgbkg_axis1","fake-fake","l");
        leg->SetFillColor(kWhite);
        leg->Draw();

        c4->cd(2);
        RooPlot *frame2final = roovar2->frame(Title("Fit axis 2 - binned"));
        dset_data->plotOn(frame2final);
        model_2D_uncorrelated->plotOn(frame2final);
        sigsigpdf->plotOn(frame2final,Normalization(fsigsig->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kRed));
        sigbkgpdf->plotOn(frame2final,Normalization(fsigbkg->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kGreen));
        bkgbkgpdf->plotOn(frame2final,Normalization(fbkgbkg->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kBlack));
        frame2final->Draw();
        //    c2->GetPad(2)->SetLogy(1);
        leg->Draw();






  c4->SaveAs("../plots/2dfitrcone_log.png");
  c4->SaveAs("../plots/2dfitrcone_log.pdf");
  c4->SaveAs("../plots/2dfitrcone_log.root");
*/
	return 0;
}









