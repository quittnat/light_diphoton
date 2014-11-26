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


bool fit2comp=true;
bool fit3comp=false;
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
RooRealVar *roodiphopt=NULL;
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
const int n_diphoptbins_forreweighting = 5;
Float_t diphoptbins_forreweighting[n_diphoptbins_forreweighting+1]={0,20,35,50,80,999};


void get_roodset_from_ttree(TDirectoryFile *f, TString treename, RooDataSet* &roodset);
void do1dfit();
void do2dfit();
void reweight_rhosigma(RooDataSet **dset, RooDataSet *dsetdestination, bool deleteold = kTRUE);
void reweight_pt_1d(RooDataSet **dset, RooDataSet *dsetdestination, int numvar);
void reweight_diphotonpt_1d(RooDataSet **dset, RooDataSet *dsetdestination);
void reweight_eta_1d(RooDataSet **dset, RooDataSet *dsetdestination, int numvar);
void ratioplot(RooDataSet *data_axis1, RooDataSet *truth_axis1, Bool_t logplot);







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
float v_roodiphopt;
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
TBranch *b_roodiphopt;
  TBranch *b_roorho;
  TBranch *b_roosigma;
  TBranch *b_roonvtx;
  TBranch *b_rooweight;

  const int nvars = 13;
  float* ptrs[nvars]={&v_roovar1,&v_roovar2,&v_roopt1,&v_roosieie1,&v_rooeta1,&v_roopt2,&v_roosieie2,&v_rooeta2,&v_roodiphopt,&v_roorho,&v_roosigma,&v_roonvtx,&v_rooweight};
  TBranch** branches[nvars]={&b_roovar1,&b_roovar2,&b_roopt1,&b_roosieie1,&b_rooeta1,&b_roopt2,&b_roosieie2,&b_rooeta2,&b_roodiphopt,&b_roorho,&b_roosigma,&b_roonvtx,&b_rooweight};
  RooRealVar* rooptrs[nvars]={roovar1,roovar2,roopt1,roosieie1,rooeta1,roopt2,roosieie2,rooeta2,roodiphopt,roorho,roosigma,roonvtx,rooweight};
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

void reweight_diphotonpt_1d(RooDataSet **dset, RooDataSet *dsetdestination){

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
  }

  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());

  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Diphoton Pt 1d rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;

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
        TCanvas* c2= new TCanvas("c2","c2");
        if(logplot){
        c2->SetLogy();
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
      cout << "data " << EB_data->GetNbinsX() <<"truth " << EB_truth->GetNbinsX() << endl;
 //     cout << "meanh1 " << meanh1 << " +/- " << err1 << "meanh2 " << meanh2 << " +/- " <<err2<<"ratio " << ratio << "+/-"<< err_ratio << endl;
  cout << "meanh1 " << meanh1  << "meanh2 " << meanh2<<"ratio " << ratio  << endl;  
    EB_data->Draw();
      EB_truth->Draw("SAME");
   TLegend* leg = new TLegend(0.55, 0.65, .9, .9);
       leg->SetFillColor(0);
       leg->AddEntry(EB_data,"2p photons " ,"pl");
       leg->AddEntry(EB_truth,"1p (1f) photon", "pl");

       leg->Draw();
        TLatex a;
        a.SetNDC();
//const char * ratio_out= Form("ratio of truth MC to data MC: %f +/-%f",ratio,err_ratio);
const char * ratio_out= Form("ratio of truth MC to data MC: %f",ratio);  
      a.DrawLatex(0.55,0.6,ratio_out);
   c2->cd();          // Go back to the main canvas before defining pad2
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
   c2->SaveAs("plots/template_comp/MCprompt_iso_fulleta_rew.png");
   c2->SaveAs("plots/template_comp/MCprompt_iso_fulleta_rew.root");
   c2->SaveAs("plots/template_comp/MCprompt_iso_fulleta_rew.pdf");

}
/////////////////////////////////////////////////////////////////////
void do1dfit(){
//TODO implement switch background or signal template
//start with signal template
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
     tbins.addUniform(10,0,1.) ;
     tbins.addUniform(5,1.,3);
       tbins.addUniform(2,3,9) ;
         roovar1->setBinning(tbins);
  roovar1->SetTitle("Iso_{1}");
  roovar2->SetTitle("Iso_{2}"); 
  rooeta1 = new RooRealVar("rooeta1","rooeta1",0,2.5);
  rooeta2 = new RooRealVar("rooeta2","rooeta2",0,2.5);
  roopt1 = new RooRealVar("roopt1","roopt1",100.);
  roopt2 = new RooRealVar("roopt2","roopt2",100.);
  roosieie1 = new RooRealVar("roosieie1","roosieie1",0,0.045);
  roosieie2 = new RooRealVar("roosieie2","roosieie2",0,0.045);
  roodiphopt = new RooRealVar("roodiphopt","roodiphopt",100.);
  roorho = new RooRealVar("roorho","roorho",20.);
  roosigma = new RooRealVar("roosigma","roosigma",2.);
  roonvtx = new RooRealVar("roonvtx","roonvtx",20.);
  rooweight = new RooRealVar("rooweight","rooweight",1.);
  assert (roovar1);
  assert (roovar1);
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
TFile* fdata_std=NULL; TFile* fdatarcone_s=NULL; TFile * fdatasideb_b2s=NULL; TFile *fdatasideb_b=NULL;  
  if(fit2comp){
    //   	fdata_std = new TFile("forroofit/sig2dpp_fulleta_truth.root","read");
 fdata_std = new TFile("forroofit/sig2dpp_fulleta_truth.root","read"); 
 	  get_roodset_from_ttree(fdata_std,"for_roofit",dataset_data_std);
//both leading and subleading photon 
 	 fdatarcone_s = new TFile("forroofit/mc1p1f_fulletarange_leadpho_template.root","read"); 
 //	 fdatarcone_s = new TFile("forroofit/comp2rcone.root ","read");
//	  fdatarcone_s = new TFile("forroofit/mc1p1f_EBEB_leadpho_template.root ","read");
	 get_roodset_from_ttree(fdatarcone_s,"for_roofit",dataset_datarcone_s);
//	 fdatasideb_b2s= new TFile("forroofit/bkgtruth_template.root","read");
	 fdatasideb_b2s= new TFile("forroofit/mc1p1f_fulletarange_subleadpho_template.root","read");
	 get_roodset_from_ttree(fdatasideb_b2s,"for_roofit",dataset_datasideb_b2s);
  }

  else if(fit3comp){
        fdata_std = new TFile("forroofit/mc2dstd_EBEB.root","read");
        get_roodset_from_ttree(fdata_std,"for_roofit",dataset_data_std);
//both leading and subleading photon 
        fdatarcone_s = new TFile("forroofit/pp_template.root ","read");
        get_roodset_from_ttree(fdatarcone_s,"for_roofit",dataset_datarcone_s);
//first fake then true photon
	fdatasideb_b= new TFile("forroofit/pf_template.root","read");
	get_roodset_from_ttree(fdatasideb_b,"for_roofit",dataset_datasideb_b);
        fdatasideb_b2s= new TFile("forroofit/ff_template.root","read");
        get_roodset_from_ttree(fdatasideb_b2s,"for_roofit",dataset_datasideb_b2s);
  }
else {cout << "no input files" << endl;}


//   TString cut_string= Form("roovar1<%f",1000000000.);
   TString cut_string= Form("roovar1<%f",rightrange-1e-5);
//  TString cut_string= Form("roovar1<%f && rooeta1<1.5 && rooeta2 < 1.5",rightrange-1e-5);
  assert(fdata_std); assert(fdatarcone_s);  assert(fdatasideb_b2s);
 if(fit3comp)assert(fdatasideb_b); 
 //cut away overflow bin
  if (dataset_data_std) dset_data_std = (RooDataSet*)(dataset_data_std->reduce(Name("dset_data_std1"),Cut(cut_string)));
  if (dataset_datarcone_s) dset_datarcone_s = (RooDataSet*)(dataset_datarcone_s->reduce(Name("dset_datarcone_s1"),Cut(cut_string)));
  if (dataset_datasideb_b2s) dset_datasideb_b2s = (RooDataSet*)(dataset_datasideb_b2s->reduce(Name("dset_datasideb_b2s1"),Cut(cut_string)));
  if(fit3comp){
  	if (dataset_datasideb_b) dset_datasideb_b = (RooDataSet*)(dataset_datasideb_b->reduce(Name("dset_datasideb_b1"),Cut(cut_string)));
  }
////create subset of dataset e.g. for 1d signal template, here leading photon 
  RooDataSet *dset_datastd_axis1 = (RooDataSet*)(dset_data_std->reduce(Name("dset_datastd_axis1"),SelectVars(RooArgList(*roovar1,*roopt1,*roosieie1,*roodiphopt,*rooeta1,*roorho,*roosigma,*roonvtx))));
  RooDataSet *dset_datarcone_axis1 = (RooDataSet*)(dset_datarcone_s->reduce(Name("dset_datarcone_axis1"),SelectVars(RooArgList(*roovar1,*roopt1,*roosieie1,*rooeta1,*roodiphopt,*roorho,*roosigma,*roonvtx))));
  RooDataSet *dset_sideb_axis1 = (RooDataSet*)(dset_datasideb_b2s->reduce(Name("dset_sideb_axis1"),SelectVars(RooArgList(*roovar1,*roopt1,*roosieie1,*rooeta1,*roorho,*roodiphopt,*roosigma,*roonvtx))));
  RooDataSet *dset_rcsideb_axis1=NULL;
  if(fit3comp){
	dset_rcsideb_axis1 = (RooDataSet*)(dset_datasideb_b->reduce(Name("dset_rcsideb_axis1"),SelectVars(RooArgList(*roovar1,*roopt1,*roosieie1,*rooeta1,*roorho,*roodiphopt,*roosigma,*roonvtx))));
}


  assert(dset_datastd_axis1);assert(dset_datarcone_axis1);assert(dset_sideb_axis1); 
  cout << "pp " << dset_datarcone_axis1->sumEntries() << endl;
  if(fit3comp)  cout << "pf " << dset_rcsideb_axis1->sumEntries() << endl;
  cout << "ff " << dset_sideb_axis1->sumEntries() << endl;
  cout << "data " << dset_datastd_axis1->sumEntries() << endl;
  cout << "ratio " << dset_datarcone_axis1->sumEntries()/ dset_datastd_axis1->sumEntries()  << endl;
  if(fit3comp) cout << "signal+ fakes " << dset_datarcone_axis1->sumEntries()+ dset_sideb_axis1->sumEntries() + dset_rcsideb_axis1->sumEntries()<< endl;

//reweight always reweight to data in signal region
  //MC
  
  reweight_rhosigma(&dset_datarcone_axis1,dset_datastd_axis1);
  reweight_rhosigma(&dset_sideb_axis1,dset_datastd_axis1);

  reweight_eta_1d(&dset_datarcone_axis1,dset_datastd_axis1,1);
  reweight_eta_1d(&dset_sideb_axis1,dset_datastd_axis1,1);
  reweight_pt_1d(&dset_sideb_axis1,dset_datastd_axis1,1);
  if(fit2comp){  
  	reweight_pt_1d(&dset_datarcone_axis1,dset_datastd_axis1,1); 
	}
  reweight_diphotonpt_1d(&dset_datarcone_axis1,dset_datastd_axis1);
  reweight_diphotonpt_1d(&dset_sideb_axis1,dset_datastd_axis1); 
 
  /* if(fit3comp){
	 
 reweight_rhosigma(&dset_rcsideb_axis1,dset_datastd_axis1);
	reweight_eta_1d(&dset_rcsideb_axis1,dset_datastd_axis1,1);
	reweight_pt_1d(&dset_rcsideb_axis1,dset_datastd_axis1,1);
 	reweight_diphotonpt_1d(&dset_rcsideb_axis1,dset_datastd_axis1);
	}	
*/
ratioplot(dset_datastd_axis1,dset_datarcone_axis1,true);
  dset_datarcone_axis1->Print();dset_sideb_axis1->Print();
/*
  if(fit3comp){dset_rcsideb_axis1->Print();}

// ratioplot(dset_datastd_axis1, dset_datarcone_axis1,kTRUE);
  //fitting
// RooAddPdf is an efficient implementation of a sum of PDFs of the form  c_1*PDF_1 + c_2*PDF_2 + ... (1-sum(c_1...c_n-1))*PDF_n
//the sum of the coefficients is enforced to be one,
// and the coefficient of the last PDF is calculated from that condition.
//first step determine single-photon purities on Iso1 and Iso2 axes
//j1=pp1 +fp fraction of photons populating axis 1
//j2= pp2+fp fraction of photons populating axis 2
//for integral 0-90 GeV, fit should not depend on initial values!  
  const float pp_init = 0.23;
  const float pf_init = 0.50;//0.33
//give initial values for purity fraction -j1 overall purity for leg 1 
RooRealVar *pp=NULL ;RooRealVar *pf=NULL; RooRealVar *j1=NULL;RooFormulaVar *fsig1=NULL;RooFormulaVar *fsig2=NULL;
  if(fit3comp){
 	 pp = new RooRealVar("pp","pp",pp_init,0,1);
         pf = new RooRealVar("pf","pf",pf_init,0,1);
//	pf->setVal(0.);
//	 pf->setConstant();
        fsig1 = new RooFormulaVar("fsig1","fsig1","pp",RooArgList(*pp));
        fsig2 = new RooFormulaVar("fsig2","fsig2","pf",RooArgList(*pf));
	}
  else if(fit2comp){
 	 j1 = new RooRealVar("j1","j1",pp_init+pf_init,0,1);
	 fsig1 = new RooFormulaVar("fsig1","fsig1","j1",RooArgList(*j1));
	}
  //produces binned RooDataHist also when RooRealVar input and weights are unbinned
  RooDataHist *drconehist_axis1 = new RooDataHist("drconehist_axis1","drconehist_axis1",RooArgSet(*roovar1));
  drconehist_axis1->add(*dset_datarcone_axis1);
  drconehist_axis1->Print("V");
  RooDataHist *dsidebhist_axis1 = new RooDataHist("dsidebhist_axis1","dsidebhist_axis1",RooArgSet(*roovar1));
dsidebhist_axis1->add(*dset_sideb_axis1);
  
  //normalizes histograms over all observables
  RooHistPdf *drconepdf_axis1 = new RooHistPdf("drconepdf_axis1","drconepdf_axis1",RooArgList(*roovar1),*drconehist_axis1);
  RooHistPdf *dsidebpdf_axis1 = new RooHistPdf("dsidebpdf_axis1","dsidebpdf_axis1",RooArgList(*roovar1),*dsidebhist_axis1);
  RooDataHist *drcsidebhist_axis1=NULL;RooHistPdf *drcsidebpdf_axis1=NULL;
  if(fit3comp){
  	drcsidebhist_axis1 = new RooDataHist("drcsidebhist_axis1","drcsidebhist_axis1",RooArgSet(*roovar1));
	drcsidebhist_axis1->add(*dset_rcsideb_axis1);
  	drcsidebpdf_axis1 = new RooHistPdf("drcsidebpdf_axis1","drcsidebpdf_axis1",RooArgList(*roovar1),*drcsidebhist_axis1);
  }
TCanvas *c8=new TCanvas("c8","c8");
c8->cd();
RooPlot *frame1bla2 = roovar1->frame(Title("test"));
drconehist_axis1->plotOn(frame1bla2,LineStyle(kDashed),LineColor(kRed),Name("sig"));
//drcsidebhist_axis1->plotOn(frame1bla2,LineStyle(kDashed),LineColor(kViolet+7),Name("sigbkg"));
dsidebhist_axis1->plotOn(frame1bla2,LineStyle(kDashed),LineColor(kBlack),Name("bkg"));
frame1bla2->Draw();
 // RooAddPdf is an efficient implementation of a sum of PDFs of the form  c_1*PDF_1 + c_2*PDF_2 + ... (1-sum(c_1...c_n-1))*PDF_n
  //the sum of the coefficients is enforced to be one,and the coefficient of the last PDF is calculated from that condition.
  //ArgList of coefficients -1 of roohistpdf
RooAddPdf *model_axis1=NULL;
  if(fit2comp){
  	model_axis1 = new RooAddPdf("model_axis1","model_axis1",RooArgList(*drconepdf_axis1,*dsidebpdf_axis1),RooArgList(*fsig1),kFALSE);
  }
  else if(fit3comp){
	model_axis1 = new RooAddPdf("model_axis1","model_axis1",RooArgList(*drconepdf_axis1,*drcsidebpdf_axis1,*dsidebpdf_axis1),RooArgList(*fsig1, *fsig2),kFALSE);
	}
//if does not converge take NumCPU strategy 2
  RooFitResult *firstpass1;
 //hesse run by default
  firstpass1 = model_axis1->fitTo(*dset_datastd_axis1, NumCPU(8), Extended(false),SumW2Error(kTRUE),Save(kTRUE));
  TFile resFile("mctruth_nooverflow_reweight.root","RECREATE"); 
  firstpass1->Write();
  resFile.Write();  
  resFile.Close();
//  firstpass2 = minuit_firstpass2->save("firstpass2","firstpass2");
  firstpass1->Print();
 //   firstpass2->Print();

  TLatex b;b.SetNDC();b.SetTextSize(0.06);b.SetTextColor(kRed);
  
TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  TLegend *leg = new TLegend(0.15,0.8,0.35,0.9);
  //leg->AddEntry("data","data","lp");
//  leg->AddEntry("fit","fit","l");
//  leg->AddEntry("background","sideband from gj + jj MC","l");
// c1->Divide(2);
  c1->cd(1);
//  gPad->SetLogy();
 TString title; 
if(fit2comp){
 title= "1d fit, mc truth no overflow, reweight";
}
else if (fit3comp){
title="1d fit, mc random cone and sideband";
}

  RooPlot *frame1bla = roovar1->frame(Title(title.Data()));
  dset_datastd_axis1->plotOn(frame1bla,Binning(tbins),Name("data"));
  model_axis1->plotOn(frame1bla,Name("fit"));
  model_axis1->plotOn(frame1bla,Components("drconepdf_axis1"),LineStyle(kDashed),LineColor(kRed),Name("signal"));
  if(fit3comp){
  	model_axis1->plotOn(frame1bla,Components("drcsidebpdf_axis1"),LineStyle(kDashed),LineColor(kViolet+7),Name("sigbkg"));
  }
  frame1bla->Draw();
  leg->AddEntry("fit","fit","l");
  if(fit3comp){
	 leg->AddEntry("signal","rcone from gg  MC","l");
	 leg->AddEntry("sigbkg","rcone + sideband from gj MC","l");
	 leg->AddEntry("background","sideband from jj MC","l");
	} 
  else if(fit2comp){leg->AddEntry("signal","signal MC","l");
 	 leg->AddEntry("background","fake  MC","l");
 	}
   leg->SetFillColor(kWhite); 
   leg->Draw();
   b.DrawLatex(0.55,0.7,"PRELIMINARY");


  TCanvas *c2 = new TCanvas("c2","c2",1200,800);
  c2->cd(1);
  gPad->SetLogy();
  RooPlot *frame1logbla = roovar1->frame(Title(title.Data()));
  dset_datastd_axis1->plotOn(frame1logbla,Binning(tbins),Name("data"));
  model_axis1->plotOn(frame1logbla,Name("fit"));
  model_axis1->plotOn(frame1logbla,Components("drconepdf_axis1"),LineStyle(kDashed),LineColor(kRed),Name("signal"));
  model_axis1->plotOn(frame1logbla,Components("dsidebpdf_axis1"),LineStyle(kDashed),LineColor(kBlack),Name("background"));
 if(fit3comp){
        model_axis1->plotOn(frame1logbla,Components("drcsidebpdf_axis1"),LineStyle(kDashed),LineColor(kViolet+7),Name("sigbkg"));
  }
  frame1logbla->Draw();
  leg->Draw();
  b.DrawLatex(0.55,0.7,"PRELIMINARY");
  model_axis1->Print();
//  model_axis2->Print();
  dset_datastd_axis1->Print();
//  dset_datastd_axis2->Print();
 TString titleout;
if(fit2comp){
 titleout= "nooverflow_reweight";
}
else if (fit3comp){
titleout="rcone_sideb";
}
 const char* outfile=Form("plots/mc_EBEB_%s", titleout.Data());
  c1->SaveAs(Form("%s_lin.png",outfile));c1->SaveAs(Form("%s_lin.root",outfile));c1->SaveAs(Form("%s_lin.pdf",outfile));
  c2->SaveAs(Form("%s_log.png",outfile));c2->SaveAs(Form("%s_log.root",outfile));c2->SaveAs(Form("%s_log.pdf",outfile));
*/
  }

//////////////////2d fit
void do2dfit(){
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
  get_roodset_from_ttree(fdata_std,"for_roofit",dataset_data_std);
//both leading and subleading photon 
// TFile *fdatarcone_s = new TFile("forroofit/1drcone_template.root","read");
 TFile *fdatarcone_s = new TFile("forroofit/2drcone_swaped_allmc_EBEB_forroofit.root","read");
  get_roodset_from_ttree(fdatarcone_s,"for_roofit",dataset_datarcone_s);
//first fake then true photon
//  TFile *fdatasideb_b= new TFile("forroofit/side_signalband_data_EBEB_forroofit.root","read");
TFile *fdatasideb_b= new TFile("forroofit/2drconesideband_swaped_allmc_EBEB_forroofit.root","read");
  //signal+ sideband region first variable signal, second sideband
// TFile *fdatasideb_b= new TFile("forroofit/1legsideband_swaped_allmc_EBEB_forroofit.root","read");
//  TFile *fdatasideb_b= new TFile("forroofit/1dtemp_bkgtruth.root","read");
 get_roodset_from_ttree(fdatasideb_b,"for_roofit",dataset_datasideb_b);
//both in sideband, randomly swapped
  TFile *fdatasideb_b2s= new TFile("forroofit/2dsideband_swaped_allmc_EBEB_forroofit.root","read");
//TFile *fdatasideb_b2s= new TFile("forroofit/1dsideb_template.root","read");
  get_roodset_from_ttree(fdatasideb_b2s,"for_roofit",dataset_datasideb_b2s);

  assert(fdata_std); assert(fdatarcone_s); assert(fdatasideb_b);  assert(fdatasideb_b2s);
//cut away overflow bin
  if (dataset_data_std) dset_data_std1 = (RooDataSet*)(dataset_data_std->reduce(Name("dset_data_std1"),Cut(Form("roovar1<%f",rightrange-1e-5))));
  dset_data_std = (RooDataSet*)(dset_data_std1->reduce(Name("dset_data_std"),Cut(Form("roovar2<%f",rightrange-1e-5))));
  if (dataset_datarcone_s) dset_datarcone_s1 = (RooDataSet*)(dataset_datarcone_s->reduce(Name("dset_datarcone_s1"),Cut(Form("roovar2<%f",rightrange-1e-5))));
  dset_datarcone_s = (RooDataSet*)(dset_datarcone_s1->reduce(Name("dset_datarcone_s"),Cut(Form("roovar1<%f",rightrange-1e-5))));
if (dataset_datasideb_b) dset_datasideb_b1 = (RooDataSet*)(dataset_datasideb_b->reduce(Name("dset_datasideb_b1"),Cut(Form("roovar1<%f",rightrange-1e-5))));
  dset_datasideb_b = (RooDataSet*)(dset_datasideb_b1->reduce(Name("dset_datasideb_b"),Cut(Form("roovar2<%f",rightrange-1e-5))));
  if (dataset_datasideb_b2s) dset_datasideb_b2s1 = (RooDataSet*)(dataset_datasideb_b2s->reduce(Name("dset_datasideb_b2s1"),Cut(Form("roovar1<%f",rightrange-1e-5))));
  dset_datasideb_b2s = (RooDataSet*)(dset_datasideb_b2s1->reduce(Name("dset_datasideb_b2s"),Cut(Form("roovar2<%f",rightrange-1e-5))));

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
 RooDataHist *sigsigdhist = new RooDataHist("sigsigdhist","sigsigdhist",RooArgList(*roovar1,*roovar2),*dset_datarcone_s);
    RooDataHist *sigbkgdhist = new RooDataHist("sigbkgdhist","sigbkgdhist",RooArgList(*roovar1,*roovar2),*dset_datasideb_b);
    RooDataHist *bkgbkgdhist = new RooDataHist("bkgbkgdhist","bkgbkgdhist",RooArgList(*roovar1,*roovar2),*dset_datasideb_b2s);
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
    RooNLLVar *model_2D_uncorrelated_noextended_nll = new RooNLLVar("model_2D_uncorrelated_noextended_nll","model_2D_uncorrelated_noextended_nll",*model_2D_uncorrelated,*dset_data_std,NumCPU(8));
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
        dset_data_std->plotOn(frame1final,Name("data"));
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
        dset_data_std->plotOn(frame2final);
        model_2D_uncorrelated->plotOn(frame2final);
        sigsigpdf->plotOn(frame2final,Normalization(fsigsig->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kRed));
        sigbkgpdf->plotOn(frame2final,Normalization(fsigbkg->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kGreen));
        bkgbkgpdf->plotOn(frame2final,Normalization(fbkgbkg->getVal(),RooAbsPdf::Relative),LineStyle(kDashed),LineColor(kBlack));
        frame2final->Draw();
        //    c2->GetPad(2)->SetLogy(1);
        leg->Draw();





  c4->SaveAs("plots/2dfitrcone_log.png");
  c4->SaveAs("plots/2dfitrcone_log.pdf");
  c4->SaveAs("plots/2dfitrcone_log.root");

}









