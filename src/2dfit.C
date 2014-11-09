#include <assert.h>
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

using namespace std; 
using namespace RooFit;
//declaration of functions and variables
RooRealVar *roovar1=NULL;
RooRealVar *roovar2=NULL;
RooRealVar *roovar1_s=NULL;
RooRealVar *roovar2_s=NULL;
RooRealVar *roovar1_b=NULL;
RooRealVar *roovar2_b=NULL;
RooRealVar *roopt1=NULL;
RooRealVar *roopt2=NULL;
RooRealVar *roosieie1=NULL;
RooRealVar *roosieie2=NULL;
RooRealVar *rooeta1=NULL;      
RooRealVar *rooeta2=NULL;
RooRealVar *roorho=NULL;
RooRealVar *roosigma=NULL;
RooRealVar *roonvtx=NULL;
RooRealVar *rooweight=NULL;
RooDataSet *dset_mctrue_s = NULL;
RooDataSet *dset_mctrue_std = NULL;
RooDataSet *dataset_data_std = NULL;
RooDataSet *dataset_datarcone_s = NULL;
RooDataSet *dset_data_std1 = NULL;
RooDataSet *dset_datarcone_s1 = NULL;
RooDataSet *dataset_datasideb_b = NULL;
 RooDataSet *dset_datasideb_b1 = NULL;
RooDataSet *dataset_datasideb_b2s = NULL;
RooDataSet *dset_datasideb_b2s1 = NULL;
RooDataSet *dset_data_std = NULL;
RooDataSet *dset_datarcone_s = NULL;
RooDataSet *dset_datasideb_b = NULL;
RooDataSet *dset_datasideb_b2s = NULL;

RooRealVar *binning_roovar1=NULL;

bool isdata=kTRUE;
const int n_ptbins_forreweighting = 4;
Float_t ptbins_forreweighting[n_ptbins_forreweighting+1]={20,35,50,80,999};



void get_roodset_from_ttree(TDirectoryFile *f, TString treename, RooDataSet* &roodset);
void do1dfit();
void reweight_rhosigma(RooDataSet **dset, RooDataSet *dsetdestination, bool deleteold = kTRUE);
void reweight_pt_1d(RooDataSet **dset, RooDataSet *dsetdestination, int numvar);
 void reweight_eta_1d(RooDataSet **dset, RooDataSet *dsetdestination, int numvar);








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

  float v_roovar1;
  float v_roovar2;
  float v_roopt1;
  float v_roosieie1;
  float v_rooeta1;
  float v_roopt2;
  float v_roosieie2;
  float v_rooeta2;
  float v_roorho;
  float v_roosigma;
  float v_roonvtx;
  float v_rooweight;

  TBranch *b_roovar1;
  TBranch *b_roovar2;
  TBranch *b_roopt1;
  TBranch *b_roosieie1;
  TBranch *b_rooeta1;
  TBranch *b_roopt2;
  TBranch *b_roosieie2;
  TBranch *b_rooeta2;
  TBranch *b_roorho;
  TBranch *b_roosigma;
  TBranch *b_roonvtx;
  TBranch *b_rooweight;

  const int nvars = 12;
  float* ptrs[nvars]={&v_roovar1,&v_roovar2,&v_roopt1,&v_roosieie1,&v_rooeta1,&v_roopt2,&v_roosieie2,&v_rooeta2,&v_roorho,&v_roosigma,&v_roonvtx,&v_rooweight};
  TBranch** branches[nvars]={&b_roovar1,&b_roovar2,&b_roopt1,&b_roosieie1,&b_rooeta1,&b_roopt2,&b_roosieie2,&b_rooeta2,&b_roorho,&b_roosigma,&b_roonvtx,&b_rooweight};
  RooRealVar* rooptrs[nvars]={roovar1,roovar2,roopt1,roosieie1,rooeta1,roopt2,roosieie2,rooeta2,roorho,roosigma,roonvtx,rooweight};
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


void reweight_rhosigma(RooDataSet **dset, RooDataSet *dsetdestination, bool deleteold){
  if (!(*dset)) return;
  (*dset)->Print();
  TH2F *hnum = new TH2F("hnum","hnum",50,0,50,15,0,15);
  TH2F *hden = new TH2F("hden","hden",50,0,50,15,0,15);
  hnum->Sumw2();
  hden->Sumw2();
  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill(fabs((*dset)->get(i)->getRealValue("roorho")),fabs((*dset)->get(i)->getRealValue("roosigma")),(*dset)->store()->weight(i));
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(fabs(dsetdestination->get(i)->getRealValue("roorho")),fabs(dsetdestination->get(i)->getRealValue("roosigma")),dsetdestination->store()->weight(i));
  }


  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

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
    //    std::cout << oldw << " " << neww << std::endl;
    newdset->add(args,neww);
  }
  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());
  delete hnum; delete hden;
  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "RhoSigma2D rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;
  (*dset)->Print();
  if (deleteold) delete old_dset;
};


void reweight_pt_1d(RooDataSet **dset, RooDataSet *dsetdestination, int numvar){

  if (!(*dset)) return;
///numerator and denominator
  TH1F *hnum = new TH1F("hnum","hnum",n_ptbins_forreweighting,ptbins_forreweighting);
  TH1F *hden = new TH1F("hden","hden",n_ptbins_forreweighting,ptbins_forreweighting);
  hnum->Sumw2();
  hden->Sumw2();

  const char* ptname=Form("roopt%d",numvar);
// RooAbsData->get*() Load a given row of data
 //getRealValue Get value of a RooAbsReal stored in set with given name. If none is found, value of defVal is returned.
  for (int i=0; i<(*dset)->numEntries(); i++){
    hden->Fill(fabs((*dset)->get(i)->getRealValue(ptname)),(*dset)->store()->weight(i));
  }
  for (int i=0; i<dsetdestination->numEntries(); i++){
    hnum->Fill(fabs(dsetdestination->get(i)->getRealValue(ptname)),dsetdestination->store()->weight(i));
  }
//MQ
 //normalize to one 
  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());

  hnum->Divide(hden);
  TH1F *h = hnum;

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_ptrew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float pt = args.getRealValue(ptname);
    float neww = oldw*h->GetBinContent(h->FindBin(fabs(pt)));
  //    std::cout << oldw << " " << neww << std::endl;
    newdset->add(args,neww);
  }


  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Pt 1d rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;

  delete old_dset;

};

void reweight_eta_1d(RooDataSet **dset, RooDataSet *dsetdestination, int numvar){

  if (!(*dset)) return;

  TH1F *hnum = new TH1F("hnum","hnum",25,0,2.5);
  TH1F *hden = new TH1F("hden","hden",25,0,2.5);
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
    newdset->add(args,neww);
  }


  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Eta 1d rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;

  delete old_dset;

};
/////////////////////////////////////////////////////////////////////
void do1dfit(){
//TODO implement switch background or signal template
//start with signal template
//1d fit following dolightcomp from Marcos code
 //taking isolation variable
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
  roorho = new RooRealVar("roorho","roorho",0,50);
  roosigma = new RooRealVar("roosigma","roosigma",0,50);
  roonvtx = new RooRealVar("roonvtx","roonvtx",0,50);
  rooweight = new RooRealVar("rooweight","rooweight",0.,100.);
  assert (roovar1);
  assert (roovar1);
  assert (roopt1);
  assert (roosieie1);
  assert (rooeta1);
  assert (roopt2);
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
  get_roodset_from_ttree(fdata_std,"for_roofit",dataset_data_std);
//both leading and subleading photon 
// TFile *fdatarcone_s = new TFile("forroofit/2drcone_data_EBEB_forroofit.root","read");
 TFile *fdatarcone_s = new TFile("forroofit/2drancomcone_swaped_allmc_EBEB_forroofit.root","read");
  get_roodset_from_ttree(fdatarcone_s,"for_roofit",dataset_datarcone_s);
//first fake then true photon
//  TFile *fdatasideb_b= new TFile("forroofit/side_signalband_data_EBEB_forroofit.root","read");
//TFile *fdatasideb_b= new TFile("forroofit/side_signalband_allmc_EBEB_forroofit.root","read"); 
  //signal+ sideband region first variable signal, second sideband
 TFile *fdatasideb_b= new TFile("forroofit/1legsideband_swaped_allmc_EBEB_forroofit.root","read");
//   TFile *fdatasideb_b= new TFile("forroofit/fullbkg_roovar1_EBEB_forroofit.root","read");
  get_roodset_from_ttree(fdatasideb_b,"for_roofit",dataset_datasideb_b);
//both in sideband, randomly swapped
  TFile *fdatasideb_b2s= new TFile("forroofit/2dsideband_swaped_allmc_EBEB_forroofit.root","read");
  get_roodset_from_ttree(fdatasideb_b2s,"for_roofit",dataset_datasideb_b2s);

  assert(fdata_std); assert(fdatarcone_s); assert(fdatasideb_b); assert(fdatasideb_b2s);
  //cut away overflow bin
  if (dataset_data_std) dset_data_std1 = (RooDataSet*)(dataset_data_std->reduce(Name("dset_data_std1"),Cut(Form("roovar1<%f",rightrange-1e-5))));
  dset_data_std = (RooDataSet*)(dset_data_std1->reduce(Name("dset_data_std"),Cut(Form("roovar2<%f",rightrange-1e-5))));
  if (dataset_datarcone_s) dset_datarcone_s1 = (RooDataSet*)(dataset_datarcone_s->reduce(Name("dset_datarcone_s1"),Cut(Form("roovar2<%f",rightrange-1e-5))));
  dset_datarcone_s = (RooDataSet*)(dset_datarcone_s1->reduce(Name("dset_datarcone_s"),Cut(Form("roovar1<%f",rightrange-1e-5))));
  if (dataset_datasideb_b) dset_datasideb_b1 = (RooDataSet*)(dataset_datasideb_b->reduce(Name("dset_datasideb_b1"),Cut(Form("roovar1<%f",rightrange-1e-5))));
  dset_datasideb_b = (RooDataSet*)(dset_datasideb_b1->reduce(Name("dset_datasideb_b"),Cut(Form("roovar2<%f",rightrange-1e-5))));
  if (dataset_datasideb_b2s) dset_datasideb_b2s1 = (RooDataSet*)(dataset_datasideb_b2s->reduce(Name("dset_datasideb_b2s1"),Cut(Form("roovar1<%f",rightrange-1e-5))));
  dset_datasideb_b2s = (RooDataSet*)(dset_datasideb_b2s1->reduce(Name("dset_datasideb_b2s"),Cut(Form("roovar2<%f",rightrange-1e-5))));
  cout << "datasets" << endl;
//  std::cout << "MC datasets" << std::endl;
  dset_data_std->Print();
  dset_datarcone_s->Print();
  dset_datasideb_b->Print();
  dset_datasideb_b2s->Print();
////create subset of dataset e.g. for 1d signal template, here leading photon 
  RooDataSet *dset_datastd_axis1 = (RooDataSet*)(dset_data_std->reduce(Name("dset_datastd_axis1"),SelectVars(RooArgList(*roovar1,*roopt1,*roosieie1,*rooeta1,*roorho,*roosigma,*roonvtx))));
  RooDataSet *dset_datastd_axis2 = (RooDataSet*)(dset_data_std->reduce(Name("dset_datastd_axis2"),SelectVars(RooArgList(*roovar2,*roopt2,*roosieie2,*rooeta2,*roorho,*roosigma,*roonvtx))));
  RooDataSet *dset_datarcone_axis1 = (RooDataSet*)(dset_datarcone_s->reduce(Name("dset_datarcone_axis1"),SelectVars(RooArgList(*roovar1,*roopt1,*roosieie1,*rooeta1,*roorho,*roosigma,*roonvtx))));
  RooDataSet *dset_datarcone_axis2 = (RooDataSet*)(dset_datarcone_s->reduce(Name("dset_datarcone_axis2"),SelectVars(RooArgList(*roovar2,*roopt2,*roosieie2,*rooeta2,*roorho,*roosigma,*roonvtx))));
  RooDataSet *dset_sigside_axis1 = (RooDataSet*)(dset_datasideb_b->reduce(Name("dset_sigside_axis1"),SelectVars(RooArgList(*roovar1,*roopt1,*roosieie1,*rooeta1,*roorho,*roosigma,*roonvtx))));
  RooDataSet *dset_sigside_axis2 = (RooDataSet*)(dset_datasideb_b->reduce(Name("dset_sigside_axis2"),SelectVars(RooArgList(*roovar2,*roopt2,*roosieie2,*rooeta2,*roorho,*roosigma,*roonvtx))));
  RooDataSet *dset_sideb_axis1 = (RooDataSet*)(dset_datasideb_b2s->reduce(Name("dset_sideb_axis1"),SelectVars(RooArgList(*roovar1,*roopt1,*roosieie1,*rooeta1,*roorho,*roosigma,*roonvtx))));
  RooDataSet *dset_sideb_axis2 = (RooDataSet*)(dset_datasideb_b2s->reduce(Name("dset_sideb_axis2"),SelectVars(RooArgList(*roovar2,*roopt2,*roosieie2,*rooeta2,*roorho,*roosigma,*roonvtx))));
  assert(dset_datastd_axis1);
  assert(dset_datastd_axis2);
  assert(dset_datarcone_axis1);
  assert(dset_datarcone_axis2);
  assert(dset_sigside_axis1);    
  assert(dset_sigside_axis2);
  assert(dset_sideb_axis1); 
  assert(dset_sideb_axis2);
 // RooDataSet *dsettry_axis1 = (RooDataSet*)(dset_datarcone_axis2->reduce(Name("dsettry_axis1"),SelectVars(RooArgList(*roovar2))));
 // assert(dsettry_axis1);
/*TH1F* EB_temp= new TH1F("EB_temp","EB_temp",100, 0.,9.);
			  for(int i=0;i<(*dsettry_axis1).numEntries();i++) {
				  RooArgSet args = *((*dsettry_axis1).get(i));
				      float iso = args.getRealValue("roovar2");
   EB_temp->Fill(iso);
}
TCanvas *c3 = new TCanvas("c3","c3",1200,800);
c3->cd(1);
c3->SetLogy();
cout << EB_temp->GetBinContent(101) << endl;
EB_temp->Draw();*/
  //reweight always reweight to data in signal region
  reweight_rhosigma(&dset_datarcone_axis1,dset_datastd_axis1);
  reweight_rhosigma(&dset_datarcone_axis2,dset_datastd_axis2);
  reweight_rhosigma(&dset_sideb_axis1,dset_datastd_axis1);
  reweight_rhosigma(&dset_sideb_axis2,dset_datastd_axis2);
 reweight_rhosigma(&dset_sigside_axis2,dset_datastd_axis2);
 reweight_rhosigma(&dset_sigside_axis1,dset_datastd_axis1); 
 //MC
 // if(!isdata){
//   reweight_pt_1d(&dset_datarcone_axis1,dset_datastd_axis1,1);
//  reweight_pt_1d(&dset_datarcone_axis2,dset_datastd_axis2,2);
//  }
  //
 reweight_pt_1d(&dset_sigside_axis1,dset_datastd_axis1,1);
 reweight_pt_1d(&dset_sigside_axis2,dset_datastd_axis2,2);
  reweight_pt_1d(&dset_sideb_axis1,dset_datastd_axis1,1);
  reweight_pt_1d(&dset_sideb_axis2,dset_datastd_axis2,2);
  reweight_eta_1d(&dset_datarcone_axis1,dset_datastd_axis1,1);
  reweight_eta_1d(&dset_datarcone_axis2,dset_datastd_axis2,2);
  reweight_eta_1d(&dset_sideb_axis1,dset_datastd_axis1,1);
  reweight_eta_1d(&dset_sideb_axis2,dset_datastd_axis2,2);
 reweight_eta_1d(&dset_sigside_axis1,dset_datastd_axis1,1);
 reweight_eta_1d(&dset_sigside_axis2,dset_datastd_axis2,2);


  dset_datarcone_axis1->Print();
  dset_datarcone_axis2->Print();
  dset_sideb_axis1->Print();
  dset_sideb_axis2->Print();
  dset_sigside_axis1->Print();
  dset_sigside_axis2->Print();
/*  RooDataSet *dsettry_axis2 = (RooDataSet*)(dset_datarcone_axis2->reduce(Name("dsettry_axis2"),SelectVars(RooArgList(*roovar2))));
   assert(dsettry_axis2);
  RooDataSet *dsettry_axisb = (RooDataSet*)(dset_datasideb_axis2->reduce(Name("dsettry_axisb"),SelectVars(RooArgList(*roovar2))));
   assert(dsettry_axisb);
  TH1F* EB_temp2= new TH1F("EB_temp2","EB_temp2",100, 0.,9.);
  for(int i=0;i<(*dsettry_axis2).numEntries();i++) {
     RooArgSet args = *((*dsettry_axis2).get(i));
     float iso = args.getRealValue("roovar2");
//   float weight = args.getRealValue("rooweight");
     EB_temp2->Fill(iso);
	}
  TH1F* EB_temp= new TH1F("EB_temp","EB_temp",100, 0.,9.);
  for(int i=0;i<(*dsettry_axisb).numEntries();i++) {
     RooArgSet args = *((*dsettry_axisb).get(i));
     float iso = args.getRealValue("roovar2");
//   float weight = args.getRealValue("rooweight");
   EB_temp->Fill(iso);
        } 
 TCanvas *c3 = new TCanvas("c3","c3",1200,800);
 c3->cd(1);
 c3->SetLogy();
 cout << EB_temp2->GetBinContent(101) << endl;
EB_temp2->SetLineColor(kRed);
 EB_temp2->Draw();
EB_temp->SetMarkerColor(kBlue);
 EB_temp->Draw("SAME");
 */
  //fitting
// RooAddPdf is an efficient implementation of a sum of PDFs of the form  c_1*PDF_1 + c_2*PDF_2 + ... (1-sum(c_1...c_n-1))*PDF_n
//the sum of the coefficients is enforced to be one,
// and the coefficient of the last PDF is calculated from that condition.
//first step determine single-photon purities on Iso1 and Iso2 axes
//j1=pp1 +fp fraction of photons populating axis 1
//j2= pp2+fp fraction of photons populating axis 2
//for integral 0-90 GeV, fit should not depend on initial values!  
     const float pp_init = 0.15;
        const float pf_init = 0.20;//0.33
//give initial values for purity fraction -j1 overall purity for leg 1 
  RooRealVar *pp = new RooRealVar("pp","pp",pp_init,0,1);
  RooRealVar *j1 = new RooRealVar("j1","j1",pp_init+pf_init,0,1);
//  RooRealVar *j1_pf = new RooRealVar("j1_pf","j1_pf",pf_init,0,1);
  RooRealVar *j2 = new RooRealVar("j2","j2",pp_init+pf_init,0,1);
//  RooRealVar *j2_pf = new RooRealVar("j2_pf","j2_pf",pf_init,0,1);
  RooFormulaVar *fsig1 = new RooFormulaVar("fsig1","fsig1","j1",RooArgList(*j1));
  RooFormulaVar *fsig2 = new RooFormulaVar("fsig2","fsig2","j2",RooArgList(*j2));
//  RooFormulaVar *fbkgbkg1 = new RooFormulaVar("fbkgbkg1","fbkgbkg1","1-fsigsig1-fsigbkg1",RooArgList(*fsigsig1,*fsigbkg1));
   // RooFormulaVar *fsigsig2 = new RooFormulaVar("fsigsig2","fsigsig2","j2_pp",RooArgList(*j2_pp));
    //  RooFormulaVar *fsigbkg2 = new RooFormulaVar("fsigbkg2","fsigbkg2","fsigsig2-j2_pp",RooArgList(*fsigsig2,*j2_pp));
  //      RooFormulaVar *fbkgbkg2 = new RooFormulaVar("fbkgbkg2","fbkgbkg2","1-fsigsig2-fsigbkg2",RooArgList(*fsigsig2,*fsigbkg2));
  RooDataHist *drconehist_axis1 = new RooDataHist("drconehist_axis1","drconehist_axis1",RooArgList(*roovar1),*dset_datarcone_axis1);
  RooDataHist *drconehist_axis2 = new RooDataHist("drconehist_axis2","drconehist_axis2",RooArgList(*roovar2),*dset_datarcone_axis2);
  RooDataHist *dsidebhist_axis1 = new RooDataHist("dsidebhist_axis1","dsidebhist_axis1",RooArgList(*roovar1),*dset_sideb_axis1);
  RooDataHist *dsidebhist_axis2 = new RooDataHist("dsidebhist_axis2","dsidebhist_axis2",RooArgList(*roovar2),*dset_sideb_axis2);
  RooDataHist *dsigsidehist_axis1 = new RooDataHist("dsigsidehist_axis1","dsigsidehist_axis1",RooArgList(*roovar1),*dset_sigside_axis1);
  RooDataHist *dsigsidehist_axis2 = new RooDataHist("dsigsidehist_axis2","dsigsidehist_axis2",RooArgList(*roovar2),*dset_sigside_axis2);
  //normalizes histograms
  RooHistPdf *drconepdf_axis1 = new RooHistPdf("drconepdf_axis1","drconepdf_axis1",RooArgList(*roovar1),*drconehist_axis1);
  RooHistPdf *drconepdf_axis2 = new RooHistPdf("drconepdf_axis2","drconepdf_axis2",RooArgList(*roovar2),*drconehist_axis2);
  RooHistPdf *dsidebpdf_axis1 = new RooHistPdf("dsidebpdf_axis1","dsidebpdf_axis1",RooArgList(*roovar1),*dsidebhist_axis1);
  RooHistPdf *dsidebpdf_axis2 = new RooHistPdf("dsidebpdf_axis2","dsidebpdf_axis2",RooArgList(*roovar2),*dsidebhist_axis2);
 RooHistPdf *dsigsidepdf_axis1 = new RooHistPdf("dsigsidepdf_axis1","dsigsidepdf_axis1",RooArgList(*roovar1),*dsigsidehist_axis1);
 RooHistPdf *dsigsidepdf_axis2 = new RooHistPdf("dsigsidepdf_axis2","dsigsidepdf_axis2",RooArgList(*roovar2),*dsigsidehist_axis2);
/* TCanvas *c2 = new TCanvas("c2","c2",1200,800);
  c2->cd();
    gPad->SetLogy();
   RooPlot *test = roovar2->frame(Title("test axis 2"));
   drconehist_axis2->plotOn(test,MarkerColor(kRed),Name("signal region test"));
   dsidebhist_axis2->plotOn(test,MarkerColor(kBlue),Name("side test"));
   dset_datastd_axis2->plotOn(test,Name("data"));
  // drconepdf_axis2->plotOn(test,LineStyle(kDashed),LineColor(kRed),Name("signal region test"));
  // dsidebpdf_axis2->plotOn(test,LineStyle(kDashed),LineColor(kBlue),Name("side region test"));
    test->Draw(); 
EB_temp2->Draw("SAME");
EB_temp->Draw("SAME");*/
 // RooAddPdf is an efficient implementation of a sum of PDFs of the form  c_1*PDF_1 + c_2*PDF_2 + ... (1-sum(c_1...c_n-1))*PDF_n
  //the sum of the coefficients is enforced to be one,and the coefficient of the last PDF is calculated from that condition.
  //ArgList of coefficients -1 of roohistpdf
  RooAddPdf *model_axis1 = new RooAddPdf("model_axis1","model_axis1",RooArgList(*drconepdf_axis1,*dsidebpdf_axis1),RooArgList(*fsig1),kFALSE);
  RooAddPdf *model_axis2 = new RooAddPdf("model_axis2","model_axis2",RooArgList(*drconepdf_axis2,*dsidebpdf_axis2),RooArgList(*fsig2),kFALSE);
  RooAbsReal * nll1;RooAbsReal * nll2;
// Construct representation of -log(L) of PDFwith given dataset. If dataset is binned, a binned likelihood is constructed.
//if does not converge take NumCPU strategy 2
  nll1 = model_axis1->createNLL(*dset_datastd_axis1, NumCPU(8), Extended(false));
  nll2 = model_axis2->createNLL(*dset_datastd_axis2, NumCPU(8), Extended(false));
  RooFitResult *firstpass1;
  RooFitResult *firstpass2;
  //RooMinimizer can minimize any RooAbsReal function with respect to its parameters. Usual choices for minimization are RooNLLVar
  RooMinimizer *minuit_firstpass1 = new RooMinimizer(*nll1);
  RooMinimizer *minuit_firstpass2 = new RooMinimizer(*nll2);
//find function minimum, calcuates function gradient, follow to local minimum, recalculate gradient, iterate until minimum found
  minuit_firstpass1->migrad();
//calculate errors by explicit finding points wherre delta -log L=0.5
  minuit_firstpass1->minos();
  minuit_firstpass2->migrad();
  minuit_firstpass2->minos();
    //minuit_firstpass->hesse();
  firstpass1 = minuit_firstpass1->save("firstpass1","firstpass1");
  firstpass2 = minuit_firstpass2->save("firstpass2","firstpass2");
    firstpass1->Print();
    firstpass2->Print();
    //out->fsig1_firstpass=fsig1->getVal();
   // out->fsig2_firstpass=fsig2->getVal();
   // out->fsig1_firstpass_err=fsig1->getPropagatedError(*firstpass);
  //  out->fsig2_firstpass_err=fsig2->getPropagatedError(*firstpass2);
// RooNLLVar *model_axis1_noextended_nll = new RooNLLVar("model_axis1_noextended_nll","model_axis1_noextended_nll",*model_axis1,*dataset_axis1,NumCPU(numcpu>1 ? numcpu/2 : 1));

  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  c1->Divide(2);
  c1->cd(1);
  gPad->SetLogy();
  RooPlot *frame1bla = roovar1->frame(Title("Fit axis 1"));
  dset_datastd_axis1->plotOn(frame1bla,Name("data"));
  model_axis1->plotOn(frame1bla,Name("fit"));
  model_axis1->plotOn(frame1bla,Components("drconepdf_axis1"),LineStyle(kDashed),LineColor(kRed),Name("signal"));
  model_axis1->plotOn(frame1bla,Components("dsigsidepdf_axis1"),LineStyle(kDashed),LineColor(kGreen-3),Name("sigside"));
  model_axis1->plotOn(frame1bla,Components("dsidebpdf_axis1"),LineStyle(kDashed),LineColor(kBlack),Name("background"));
  frame1bla->Draw();

  TLegend *leg = new TLegend(0.2,0.8,0.45,0.9);
  leg->AddEntry("data","data","lp");
  leg->AddEntry("fit","fit","l");
  leg->AddEntry("signal","2 rcone","l");
   leg->AddEntry("sigside","sig + side","l");
  leg->AddEntry("background","2 sideb","l");
  leg->SetFillColor(kWhite);
  leg->Draw();
  c1->cd(2);
  gPad->SetLogy();
  RooPlot *frame2bla = roovar2->frame(Title("Fit axis 2 "));
  dset_datastd_axis2->plotOn(frame2bla,Name("data"));
  model_axis2->plotOn(frame2bla,Name("fit"));
  model_axis2->plotOn(frame2bla,Components("drconepdf_axis2"),LineStyle(kDashed),LineColor(kRed),Name("signal"));
  model_axis2->plotOn(frame2bla,Components("dsigsidepdf_axis2"),LineStyle(kDashed),LineColor(kGreen-3),Name("sigside"));
  model_axis2->plotOn(frame2bla,Components("dsidebpdf_axis2"),LineStyle(kDashed),LineColor(kBlack),Name("background"));
  frame2bla->Draw();
  leg->Draw();

  model_axis1->Print();
  model_axis2->Print();
  dset_datastd_axis1->Print();
  dset_datastd_axis2->Print();

  c1->SaveAs("plots/fittingplot_mctest_log.png");
  c1->SaveAs("plots/fittingplot_mctest_log.pdf");
  c1->SaveAs("plots/fittingplot_mctest_log.root");

}









