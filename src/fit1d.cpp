//TODO on name for current date changes everywere
//Change counting a bit
#include <assert.h>
#include "TLine.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPad.h"
#include "TStyle.h"
#include "RooFitResult.h"
#include "TString.h"
#include "TFile.h"
#include "TLatex.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
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
#include "RooKeysPdf.h"
#include "RooNDKeysPdf.h"
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
#include "fit1d.h"
#include <ctime> 
#define DTTMFMT "%Y-%m-%d"
#define DTTMSZ 11

using namespace std; 
using namespace RooFit;

static char *getDtTm (char *buff) {
	   time_t t = time (0);
	   strftime (buff, DTTMSZ, DTTMFMT, localtime (&t));
	   return buff;
}
//check how much memory the jobs use, 4 GB normal
 void print_mem(){
	 gSystem->GetProcInfo(&procinfo);
	 std::cout << "Resident mem (kB): " << procinfo.fMemResident << std::endl;
	 std::cout << "Virtual mem (kB):  " << procinfo.fMemVirtual << std::endl;
	 gSystem->Sleep((UInt_t)1e3);
	  };

//computes mass quantiles for diphoton mass binning and plots the result
int quantiles(TH1D* mass,double probtemp[], double dpmqtemp[] ) {
		cout << "define mass bins " << endl;
    //    mass->Scale(1.0/mass->Integral());     
        for (Int_t i =0;i<=nq_q;i++){ 
		//first entry e.g. 0.05 - diphopt 73 to end up with 10 bins in the end	
	    probtemp[i] =  Float_t(i+startbin_q)/ntot_q;
		}
        mass->GetQuantiles(nq_q+1,dpmqtemp,probtemp);
        //show the original histogram in the top pad
			TCanvas *cq = new TCanvas("cq","mass quantiles",10,10,700,900);
			cq->Divide(1,2);
			cq->cd(1);
			cq->SetLogy();
			cq->Update();
			mass->Draw();

			// show the quantiles in the bottom pad
			cq->cd(2);
			gPad->SetGrid();
			TGraph *gr = new TGraph(nq_q+1,probtemp,dpmqtemp);
			gr->SetMarkerStyle(21);
			gr->Draw("alp");
			cq->SaveAs(Form("%smassquantiles_%s_range_%u_%u_%s.root",dir.Data(),(truthfit)? "truth": "rcone_sideb",startbin_q,endbin_q,eta_q.Data()));
			cq->SaveAs(Form("%smassquantiles_%s_range_%u_%u_%s.png",dir.Data(),(truthfit)? "truth": "rcone_sideb",startbin_q,endbin_q,eta_q.Data()));
//   	
	delete cq;
	delete gr;
	return 0;
	}

//final purity plot vs diphoton mass bins
//TODO plot as bin on the x-axis, not as marker point
int plot_purity_massbin(double mass[], double masserr[],double pur[],double purerr[],bool berr)  {
	TCanvas* cgr =new TCanvas("cgr","cgr");
 	TLegend* leg = new TLegend(0.6,0.7,0.7,0.8);
   	leg->SetFillColor(10);
    leg->SetLineColor(10);
    leg->SetTextSize(0.03);
    leg->SetTextFont(42);
   // TGraphAsymmErrors *gr = new TGraphAsymmErrors(nq_q,mass, pur,masserr,masserr,pullerrLo,purerrHi);
    TGraphErrors *gr = new TGraphErrors(nq_q,mass, pur,masserr,purerr);
	cgr->SetGridy();
	gr->SetTitle(Form("mgg %s",eta_q.Data()));
	gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.3);
    gr->GetYaxis()->SetTitle("purity");
    gr->GetXaxis()->SetTitle("diphoton mass [GeV]");
    gr->GetYaxis()->SetRangeUser(0.,1.);
    gr->GetYaxis()->SetTitleOffset(1.2);
    leg->AddEntry(gr,"template fit purity ","p");    
    gr->Draw("AP");
    leg->Draw();
    gr->SaveAs(Form("%spurity_%s_%s_%s_%s_massbin_%u_range_%u_%u.root",dir.Data(),date.Data(),(truthfit) ? "truth": "rcone_sideb",(berr)? "sumw2erron": "sumw2erroff",eta_q.Data(),ntot_q,startbin_q,endbin_q));
    cgr->Print(Form("%spurity_%s_%s_%s_%s_massbin_%u_range_%u_%u.png",dir.Data(),date.Data(),(truthfit) ? "truth": "rcone_sideb",(berr)? "sumw2erron": "sumw2erroff",eta_q.Data(),ntot_q,startbin_q,endbin_q),"png");
    
	delete cgr;
	delete leg;
    delete gr;
	return 0;
}
//get the RooDataSet, multiply each roorealvar by the event_weight
void get_roodset_from_ttree(TDirectoryFile *f, TString treename, RooDataSet* &roodset){
  cout << "Creating roodset from file: " << f->GetName() << " with tree: " << treename.Data() << endl;

  TTree *t = NULL;
  assert(roodset==NULL);
  f->GetObject(treename.Data(),t);
  if (!t) {cout << "Impossible to find TTree " << treename.Data() << endl; return;}
  TObjArray *objs = t->GetListOfBranches();
  //disables all branches
  t->SetBranchStatus("*",0);


  double v_rooisopv1;
  double v_rooisopv2;
  double v_rooisowv1;
  double v_rooisowv2;
  double v_roovar1;
  double v_roovar2;
  double v_roopt1;
  double v_roosieie1;
  double v_rooeta1;
  double v_roopt2;
  double v_roosieie2;
  double v_rooeta2;
  double v_roodiphopt;
  double v_roodiphomass;
  double v_roorho;
  double v_roosigma;
  double v_roonvtx;
  double v_rooweight;

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
  TBranch *b_roodiphomass;
  TBranch *b_roorho;
  TBranch *b_roosigma;
  TBranch *b_roonvtx;
  TBranch *b_rooweight;

  const int nvars = 18;
  double* ptrs[nvars]={&v_roovar1,&v_roovar2,&v_rooisopv1,&v_rooisopv2,&v_rooisowv1,&v_rooisowv2,&v_roopt1,&v_roosieie1,&v_rooeta1,&v_roopt2,&v_roosieie2,&v_rooeta2,&v_roodiphopt,&v_roodiphomass,&v_roorho,&v_roosigma,&v_roonvtx,&v_rooweight};
  TBranch** branches[nvars]={&b_roovar1,&b_roovar2,&b_rooisopv1,&b_rooisopv2,&b_rooisowv1,&b_rooisowv2,&b_roopt1,&b_roosieie1,&b_rooeta1,&b_roopt2,&b_roosieie2,&b_rooeta2,&b_roodiphopt,&b_roodiphomass,&b_roorho,&b_roosigma,&b_roonvtx,&b_rooweight};
  RooRealVar* rooptrs[nvars]={roovar1,roovar2,rooisopv1,rooisopv2,rooisowv1,rooisowv2,roopt1,roosieie1,rooeta1,roopt2,roosieie2,rooeta2,roodiphopt,roodiphomass,roorho,roosigma,roonvtx,rooweight};
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

//2d reweighting of rho and its sigma
void reweight_rhosigma(TH2F* weight_eta, TH2F* weight_pt,TH2F*weight_diphopt,TH2F*weightratiobin_rho,TH2F* weightratio_rho, TH2F* weightratiobin_sigma,TH2F*weightratio_sigma,RooDataSet **dset, RooDataSet *dsetdestination, bool deleteold){
  if (!(*dset)) return;
  TH2F *hnum = new TH2F("hnum","hnum",n_rhobins_forreweighting,rhobins_forreweighting,n_sigmabins_forreweighting,sigmabins_forreweighting);
  TH2F *hden = new TH2F("hden","hden",n_rhobins_forreweighting,rhobins_forreweighting,n_sigmabins_forreweighting,sigmabins_forreweighting);
//  TH2F *hnum = new TH2F("hnum","hnum",100,0,100,20,0,20);
 // TH2F *hden = new TH2F("hden","hden",100,0,100,20,0,20);
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
//data/MC
  hnum->Divide(hden);
  TH2F *h = hnum;

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_rhosigmarew",(*dset)->GetName()));
  newdset->reset();
  int zero=0;
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    double oldw = (*dset)->store()->weight(i);
    double rho = args.getRealValue("roorho");
    double sigma = args.getRealValue("roosigma");
    double eta = fabs(args.getRealValue("rooeta1"));
    double pt = fabs(args.getRealValue("roopt1"));
    double diphopt = fabs(args.getRealValue("roodiphopt"));
    double neww = oldw*h->GetBinContent(h->FindBin(rho,sigma));
	//if(debug){
		if(oldw!=0)
		{
		weight_diphopt->Fill(diphopt, oldw/neww);
		weight_eta->Fill(eta, oldw/neww);
		weight_pt->Fill(pt, oldw/neww);

		}
	    if(oldw!=0)weightratiobin_rho->Fill(h->GetXaxis()->FindBin(rho),oldw/neww);
		if(oldw!=0)weightratio_rho->Fill(rho,oldw/neww);
		if(oldw!=0)weightratiobin_sigma->Fill(h->GetYaxis()->FindBin(sigma),oldw/neww);
		if(oldw!=0)weightratio_sigma->Fill(sigma,oldw/neww);
		else{zero++;}
//	}
		newdset->add(args,neww);
	  }
	  newdset->SetName((*dset)->GetName());
	  newdset->SetTitle((*dset)->GetTitle());
	  delete hnum; delete hden;
	  RooDataSet *old_dset = *dset;
	  *dset=newdset;
	  std::cout << "RhoSigma2D rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << " # oldw = 0: " << zero << std::endl;
	  if (deleteold) delete old_dset;
	};


void reweight_pt_1d(TH2F* weight_rho, TH2F* weight_sigma,TH2F*weight_eta,TH2F*weight_diphopt,TH2F* weightratiobin_pt, TH2F* weightratio_pt,RooDataSet **dset, RooDataSet *dsetdestination, int numvar){

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
  
 //normalize to one 
  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());
  hnum->Divide(hden);
  TH1F *h = hnum;

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_ptrew",(*dset)->GetName()));
  newdset->reset();
  int zero=0;
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    double oldw = (*dset)->store()->weight(i);
    double pt = args.getRealValue(ptname);
    double neww = oldw*h->GetBinContent(h->FindBin(fabs(pt)));
    if(oldw!=0){
	    double diphopt=fabs(args.getRealValue("roodiphopt"));
		double rho=args.getRealValue("roorho");
		double sigma=fabs(args.getRealValue("roosigma"));
		double eta=fabs(args.getRealValue("rooeta1"));
		weight_diphopt->Fill(diphopt, oldw/neww);
		weight_eta->Fill(eta, oldw/neww);
		weight_rho->Fill(rho, oldw/neww);
		weight_sigma->Fill(sigma, oldw/neww);
	}
	
		if(oldw!=0)weightratiobin_pt->Fill(h->FindBin(fabs(pt)),oldw/neww);
		if(oldw!=0)weightratio_pt->Fill(pt,oldw/neww);
		else{zero++;}
	 newdset->add(args,neww);
  }

  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());
  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Pt 1d rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << " # oldw=0: " << zero << std::endl;
  delete old_dset;
};

void reweight_diphotonpt_1d(TH2F*weight_rho,TH2F*weight_sigma,TH2F*weight_eta,TH2F* weight_pt,TH2F*weight_ratiobin, TH2F* weight_ratio,RooDataSet **dset, RooDataSet *dsetdestination){

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
 //normalize to one 
  hnum->Scale(1.0/hnum->Integral());
  hden->Scale(1.0/hden->Integral());
  hnum->Divide(hden);
  TH1F *h = hnum;
  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_diphoptrew",(*dset)->GetName()));
  newdset->reset();
  assert(newdset);
  int zero=0;
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    double oldw = (*dset)->store()->weight(i);
    double diphopt = fabs((*dset)->get(i)->getRealValue("roodiphopt"));
    double neww = oldw*h->GetBinContent(h->FindBin(fabs(diphopt)));
    newdset->add(args,neww);
   // if(debug){
	    double eta=fabs(args.getRealValue("rooeta1"));
		double pt=fabs(args.getRealValue("roopt1"));
		double rho=args.getRealValue("roorho");
		double sigma=fabs(args.getRealValue("roosigma"));
		if(oldw!=0){
			weight_pt->Fill(pt,oldw/neww);
			weight_eta->Fill(eta, oldw/neww);
			weight_sigma->Fill(sigma, oldw/neww);
			weight_rho->Fill(rho, oldw/neww);
		}
		if(oldw!=0)weight_ratiobin->Fill(h->FindBin(fabs(diphopt)),oldw/neww);
//		else {weight_ratiobin->Fill(-10,1);}//cout << "dipho weight old 0" << endl;}
		if(oldw!=0)weight_ratio->Fill(diphopt,oldw/neww);
 		else{zero++;}
  }

  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());
  delete hnum; delete hden;
 
  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Diphoton Pt 1d rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << "# oldw=0 : " << zero << std::endl;

  delete old_dset;

};



void reweight_eta_1d(TH2F* weight_rho,TH2F* weight_sigma,TH2F*weight_pt,TH2F* weight_diphopt,TH2F* weight_ratiobin,TH2F*weight_ratio,RooDataSet **dset, RooDataSet *dsetdestination, int numvar){

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

  hnum->Divide(hden);
  TH1F *h = hnum;
  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_etarew",(*dset)->GetName()));
  newdset->reset();
  int zero=0;
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    double oldw = (*dset)->store()->weight(i);
    double eta = args.getRealValue(etaname);
    double neww = oldw*h->GetBinContent(h->FindBin(fabs(eta)));
  //  if(debug){ 
		double pt=fabs(args.getRealValue("roopt1"));
		double diphopt=fabs(args.getRealValue("roodiphopt"));
		double rho=args.getRealValue("roorho");
		double sigma=fabs(args.getRealValue("roosigma"));
		if(oldw!=0){
		weight_diphopt->Fill(diphopt, oldw/neww);	
		weight_pt->Fill(pt,oldw/neww);
		weight_rho->Fill(rho, oldw/neww);
		weight_sigma->Fill(sigma, oldw/neww);
		}
		if(oldw!=0)weight_ratiobin->Fill(h->FindBin(fabs(eta)),oldw/neww);
	   if(oldw!=0)weight_ratio->Fill(fabs(eta),oldw/neww);
	   else{zero++;}
	   newdset->add(args,neww);
  }

  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());
  delete hnum; delete hden;

  RooDataSet *old_dset = *dset;
  *dset=newdset;
  std::cout << "Eta 1d rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << " # oldw=0 "  << zero << std::endl;
  delete old_dset;

};
//plot to compare two distributions and plot ther ratio
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
   gStyle->SetOptStat(111111);
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
int prep_dataset(Bool_t isdata,TString eta,Bool_t truthf, int ntotbins, int startbin, int endbin){
 //template from MC truth or with random cone/sideband?
//  dir=Form("../plots/May5/%s/massbinned20_4bins_nonweighted/",eta.Data());
//  dir=Form("../plots/April29/%s/massbinned20_4bins_onlyrhosigmaetarew/",eta.Data());
  dir=Form("../plots/May5/test/");
  truthfit=truthf;
  // EBEB or nEBEB
  eta_q=eta;
  isdata_q=isdata;
  //startbin for massquantiles, important if massbins splitted in different jobs
  startbin_q=startbin;
  //endbin for massquantiles! loop starts with i=1 such that the [0] entry of array is 0 and not the diphoptbin of first quantile. Therefore endbin=#total bins (ntotbins) +1
  endbin_q=endbin;
  ntot_q=ntotbins;

  if(eta=="EBEB"){EBEB=true;}
  else if (eta=="nEBEB"){EBEB=false;}
  else cout << "wrong eta range " << endl;
  const Int_t nq=endbin-startbin;
  nq_q=nq;
  char buff[DTTMSZ];
  date=getDtTm(buff);
  cout << getDtTm(buff)<< endl;
  cout << "ntot bins" << ntotbins << endl;
  cout << "massbins " << nq << " startbin " << startbin << " endbin " << endbin << endl;
  
//100 mass bins 10-15 % statistics on the purity fit are sufficient
  Double_t prob[nq+1];  // position where to compute the quantiles in [0,1]
  Double_t dpmq[nq+1];  // array to contain the quantiles
  Double_t data_entries[nq];  
  //arrays to hold results from 1d fit in different massbins with and without sumw2error enabled
  Double_t purity[nq];  
  Double_t purityerr[nq]; 
  Double_t pullerror[nq];  
  //define roovariables
  gStyle->SetOptStat(111111);
  Float_t leftrange=0. ;
  Float_t rightrange=9.;
  TH1F::SetDefaultSumw2(kTRUE);
  roovar1 = new RooRealVar("roovar1","roovar1",leftrange,rightrange);
  roovar2 = new RooRealVar("roovar2","roovar2",leftrange,rightrange);
  roovar1->setRange(leftrange,rightrange);
  roovar2->setRange(leftrange,rightrange);
   //tbins.addUniform(90,0,9) ;
   //6 bins for purity fit in different massbins
  if(massbinned){
  	tbins.addUniform(1, 0.,0.1);
  	tbins.addUniform(1, 0.1,1);
  	tbins.addUniform(1, 1,5);
  	tbins.addUniform(1, 5,9);
  }

  //only one diphoton mass-binning a bit finer
  if (!massbinned)
  {
	  
  	tbins.addUniform(1, 0.,0.1);
  	tbins.addUniform(1, 0.1,1);
  	tbins.addUniform(1, 1,5);
  	tbins.addUniform(1, 5,9);
/*
  	  
    tbins.addUniform(3,0.,.5) ;
  	tbins.addUniform(2,0.5,1.) ;
  	tbins.addUniform(6,1.,3.) ;
	tbins.addUniform(3,3.,6.) ;
	tbins.addUniform(3,6.,9.);
 */
   	} 
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
  
 //pt binning mainly for debugging and to compare distributions 
  Float_t ptleftrange=0. ;
  Float_t ptrightrange=3000.;
  RooBinning ptbins (ptleftrange,ptrightrange); 
  ptbins.addUniform(20,20,40) ;
  ptbins.addUniform(10,40,100) ;
  ptbins.addUniform(5,100,200) ;
  ptbins.addUniform(2,200,400) ;
  ptbins.addUniform(1,400,3000) ;
  Float_t dptleftrange=0. ;
  Float_t dptrightrange=1000.;
  RooBinning dptbins (dptleftrange,dptrightrange); 
  dptbins.addUniform(80,0,80) ;
  dptbins.addUniform(40,80,200) ;
  dptbins.addUniform(20,200,400) ;
  dptbins.addUniform(1,400,3000) ;

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
  Float_t dpmleftrange=0. ;
  Float_t dpmrightrange=7000.;
  RooBinning dpmbins (dpmleftrange,dpmrightrange); 
  dpmbins.addUniform(80,0.,400) ;
  dpmbins.addUniform(26,400,3000) ;
  dpmbins.addUniform(1,3000,7000) ;
  roodiphomass = new RooRealVar("roodiphomass","roodiphomass",dpmleftrange,dpmrightrange);
  roodiphomass->setRange(dpmleftrange,dpmrightrange);
  roodiphomass->setBinning(dpmbins);
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
  assert (roodiphomass);
  assert (roorho);
  assert (roosigma);
  assert (roonvtx);
  assert (rooweight);
  TFile* fdata=NULL; TFile* frcone=NULL; TFile * fsideband=NULL; TFile* testfilesig=NULL; TFile* testfilebkg=NULL;
  
  if(isdata_q)
  {
	   fdata = new TFile("../forroofit/April9/data.root","read"); 
//	    fdata = new TFile("../forroofit/April1/datamc.root","read"); 
	   get_roodset_from_ttree(fdata,"for_roofit",dataset_data);
	//both leading and subleading photon 
  	   frcone = new TFile("../forroofit/April9/rconedata.root","read"); 	
  	   get_roodset_from_ttree(frcone,"for_roofit",dataset_rcone);
  	   fsideband= new TFile("../forroofit/April9/sidebanddata.root","read");
  	   get_roodset_from_ttree(fsideband,"for_roofit",dataset_sideband);
  	if(check) 
		{	
  			testfilesig= new TFile("../forroofit/April1/rconemc.root","read");
			get_roodset_from_ttree(testfilesig,"for_roofit",dataset_test1);
  			testfilebkg= new TFile("../forroofit/April1/sidebandmc.root","read");
			get_roodset_from_ttree(testfilebkg,"for_roofit",dataset_test2);
  		}
  }		
  if(!isdata_q)
  {
	  if(!truthfit)
	  { 
	    fdata = new TFile("../forroofit/April1/datamc.root","read"); 
	    get_roodset_from_ttree(fdata,"for_roofit",dataset_data);
  		frcone = new TFile("../forroofit/April1/rconemc.root","read"); 	
  		get_roodset_from_ttree(frcone,"for_roofit",dataset_rcone);
  		fsideband= new TFile("../forroofit/April1/sidebandmc.root","read");
  		get_roodset_from_ttree(fsideband,"for_roofit",dataset_sideband);
  		if(check) 
		{	
  			testfilesig= new TFile("../forroofit/April1/sig.root","read");
			get_roodset_from_ttree(testfilesig,"for_roofit",dataset_test1);
  			testfilebkg= new TFile("../forroofit/April1/bkg.root","read");
			get_roodset_from_ttree(testfilebkg,"for_roofit",dataset_test2);
  		}
 	 }
	 if(truthfit)
	 { 
 	   fdata = new TFile("../forroofit/April1/datamc.root","read"); 
	   get_roodset_from_ttree(fdata,"for_roofit",dataset_data);
  		frcone = new TFile("../forroofit/April1/sig.root","read"); 	
  		get_roodset_from_ttree(frcone,"for_roofit",dataset_rcone);
    	fsideband= new TFile("../forroofit/April1/bkg.root","read");
  		get_roodset_from_ttree(fsideband,"for_roofit",dataset_sideband);
	  	if(check) {	
  			testfilesig= new TFile("../forroofit/April1/sig.root","read");
			get_roodset_from_ttree(testfilesig,"for_roofit",dataset_test1);
  			testfilebkg= new TFile("../forroofit/April1/bkg.root","read");
			get_roodset_from_ttree(testfilebkg,"for_roofit",dataset_test2);
			assert(testfilesig);assert(testfilebkg);
	  	}
    }
  }
   assert(fdata); assert(frcone);  assert(fsideband);
   	cout << "Entries truth/data weighted " << dataset_rcone->sumEntries()/ dataset_data->sumEntries()  << endl;

  // cut datasets to etarange and cut away overflow bin
  TString cut_string; 
  if(!EBEB){
   	cut_string= Form("(roovar1<%f && roovar2<%f)  && (rooeta1 > 1.566 || rooeta2 > 1.566)",rightrange-1e-5, rightrange-1e-5);
  }
  else if(EBEB){
	cut_string= Form(" (roovar1<%f && roovar2<%f) && (rooeta1<1.4442 && rooeta2 < 1.4442)",rightrange-1e-5,rightrange-1e-5);
  }

  RooDataSet *dset_data = NULL;
  RooDataSet *dset_rcone = NULL;
  RooDataSet *dset_sideb = NULL;
  RooDataSet *dataset_testsig = NULL;
  RooDataSet *dataset_testbkg = NULL;
  RooDataSet *dset_rcone_allm_axis1 = NULL;
  RooDataSet *dset_sideb_allm_axis1 = NULL;
  if (dataset_data)dset_data = (RooDataSet*)(dataset_data->reduce(Name("dset_data"),Cut(cut_string)));
  if (dataset_rcone) dset_rcone = (RooDataSet*)(dataset_rcone->reduce(Name("dset_rcone"),Cut(cut_string)));
  if (dataset_sideband)dset_sideb = (RooDataSet*)(dataset_sideband->reduce(Name("dset_sideb"),Cut(cut_string)));
  if (dataset_test1) dataset_testsig = (RooDataSet*)(dataset_test1->reduce(Name("dataset_testsig"),Cut(cut_string)));
  if (dataset_test2) dataset_testbkg = (RooDataSet*)(dataset_test2->reduce(Name("dataset_testbkg"),Cut(cut_string)));
////create subset of dataset e.g. for 1d signal template, here leading photon 
  RooDataSet *dset_data_allm_axis1 = (RooDataSet*)(dset_data->reduce(Name("dset_data_allm_axis1"),SelectVars(RooArgList(*roovar1,*rooisopv1,*roopt1,*rooeta1,*roosieie1,*roodiphopt,*roodiphomass,*roorho,*roosigma))));
 if(dataset_rcone)  dset_rcone_allm_axis1 = (RooDataSet*)(dset_rcone->reduce(Name("dset_rcone_allm_axis1"),SelectVars(RooArgList(*roovar1,*rooisopv1,*roopt1,*roosieie1,*rooeta1,*roodiphopt,*roodiphomass,*roorho,*roosigma))));
 if(dataset_sideband) dset_sideb_allm_axis1 = (RooDataSet*)(dset_sideb->reduce(Name("dset_sideb_allm_axis1"),SelectVars(RooArgList(*roovar1,*rooisopv1,*roopt1,*roosieie1,*rooeta1,*roorho,*roodiphopt,*roodiphomass,*roosigma))));
  dset_data_allm_axis1->Print();
 if(dataset_rcone) dset_rcone_allm_axis1->Print();
 if(dataset_sideband) dset_sideb_allm_axis1->Print();
if(dataset_rcone)  cout << "ratio after eta and overflow cut Entries truth/data weighted " << dset_rcone_allm_axis1->sumEntries()/ dset_data_allm_axis1->sumEntries()  << endl;
if(check)
{
	dset_testbkg1 = (RooDataSet*)(dataset_testbkg->reduce(Name("dset_testbkg1"),SelectVars(RooArgList(*roovar1,*rooisopv1,*roopt1,*roosieie1,*rooeta1,*roorho,*roodiphopt,*roodiphomass,*roosigma))));
	  dset_testsig1 = (RooDataSet*)(dataset_testsig->reduce(Name("dset_testsig1"),SelectVars(RooArgList(*roovar1,*rooisopv1,*roopt1,*roosieie1,*rooeta1,*roorho,*roodiphopt,*roodiphomass,*roosigma))));
}
	  //define mass bins
 
  RooArgSet massarg;
  massarg.add(*roodiphomass); 
  RooDataHist* mass_data = new RooDataHist("mass_data","mass_data",massarg);
  mass_data->add(*dset_data_allm_axis1);
  RooHistPdf* mass_datapdf = new RooHistPdf("mass_datapdf","mass_datapdf",massarg,*mass_data);
  TH1D* diphomass_data=(TH1D*)(mass_datapdf->createHistogram("diphomass_data",*roodiphomass,Binning(dpmbins)));
  TCanvas* cbla=new TCanvas("cbla","cbla");
  cbla->cd();
  diphomass_data->Draw();
  assert(diphomass_data);
 /* 
  RooArgSet isoarg;
  isoarg.add(*roovar1); 
  RooDataHist* iso_data = new RooDataHist("iso_data","iso_data",isoarg);
  iso_data->add(*dset_rcone_allm_axis1);
  RooHistPdf* iso_datapdf = new RooHistPdf("iso_datapdf","iso_datapdf",isoarg,*iso_data);
   TH1D* tiso=(TH1D*)(iso_datapdf->createHistogram("tiso",*roovar1,Binning(tbins)));
  RooDataHist* iso_test = new RooDataHist("iso_test","iso_test",isoarg);
  iso_test->add(*dset_testsig1);
  RooHistPdf* iso_testpdf = new RooHistPdf("iso_testpdf","iso_testpdf",isoarg,*iso_test);
   TH1D* tisotest=(TH1D*)(iso_testpdf->createHistogram("tiso",*roovar1,Binning(tbins)));
  TCanvas* ciso=new TCanvas("ciso","ciso");
  ciso->cd();
  tiso->SetLineColor(kRed);
  tiso->Draw();
  tisotest->Draw("SAME");
  */
  // compute mass quantiles to have the same number of weighted entries per massbin
  if(massbinned)
  {
  	quantiles(diphomass_data,prob,dpmq); 
    for(int k=0;k<=nq;k++)
	cout << "prob " << prob[k] << " diphomass " << dpmq[k]  << endl; 
  }
 //compare templates bkg and sideband to rcone and sideband to see if they are described correctly
if(check) { 
//	dset_testbkg1 = (RooDataSet*)(dataset_testbkg->reduce(Name("dset_testbkg1"),SelectVars(RooArgList(*roovar1,*rooisopv1,*roopt1,*roosieie1,*rooeta1,*roorho,*roodiphopt,*roodiphomass,*roosigma))));
//	  dset_testsig1 = (RooDataSet*)(dataset_testsig->reduce(Name("dset_testsig1"),SelectVars(RooArgList(*roovar1,*rooisopv1,*roopt1,*roosieie1,*rooeta1,*roorho,*roodiphopt,*roodiphomass,*roosigma))));
	  TCanvas *c5=new TCanvas("c5","isolation comparison");
	  c5->cd();
	  c5->SetLogy();
	  TLegend *leg5 = new TLegend(0.6,0.8,0.9,0.9);
	  leg5->SetFillColor(kWhite); 
	  RooPlot *comp = roovar1->frame(Title("comparison before reweighting"));
	 cout << "dset_rcone_allm_axis1->sumEntries()  " << dset_rcone_allm_axis1->sumEntries() << " dset_rcone_allm_axis1->sumEntries(0.,0.1) " << dset_rcone_allm_axis1->sumEntries("roovar1","") 
	<< " dset_rcone_allm_axis1->sumEntries(0.,9) " << dset_rcone_allm_axis1->sumEntries("roovar1","< 9.") << endl;

	 //	  iso_datapdf->plotOn(comp,LineColor(kBlack),Name("rcone"));
     //  iso_testpdf->plotOn(comp,LineColor(kBlue),Name("sigtruth"));
	  dset_rcone_allm_axis1->plotOn(comp,Rescale(1./(dset_rcone_allm_axis1->sumEntries())),Binning(tbins),MarkerStyle(20),MarkerColor(kBlack),LineColor(kBlack),Name("rcone"));
	  dset_testsig1->plotOn(comp,Rescale(1./(dset_testsig1->sumEntries())),Binning(tbins),MarkerStyle(20),MarkerColor(kBlue),LineColor(kBlue),Name("sigtruth"));
	  dset_sideb_allm_axis1->plotOn(comp,Rescale(1./(dset_sideb_allm_axis1->sumEntries())),Binning(tbins),MarkerStyle(20),MarkerColor(kRed),LineColor(kRed),Name("sideband"));
	  dset_testbkg1->plotOn(comp,Rescale(1./(dset_testbkg1->sumEntries())),Binning(tbins),MarkerStyle(20),MarkerColor(kMagenta),LineColor(kMagenta),Name("bkgtruth"));
	  comp->Draw();
	  comp->SetAxisRange(0.,1e-3,"Y");
	  if(!isdata_q)
      {
	  	leg5->AddEntry("rcone","random cone","l");
	   	leg5->AddEntry("sigtruth","sig truth","l");
	    leg5->AddEntry("bkgtruth","bkg truth","l");
	    leg5->AddEntry("sideband","sideband","l");
	  }
	  else if(isdata_q)
	  {
	  	leg5->AddEntry("rcone","random cone data","l");
	    leg5->AddEntry("sideband","sideband data","l");
	   	leg5->AddEntry("sigtruth","random cone mc","l");
	    leg5->AddEntry("bkgtruth","sideband mc","l");

	  }
	  leg5->Draw();
	  c5->SaveAs(Form("%sfullmassrange_template_comparison_noreweight_%slog.root",dir.Data(),eta_q.Data())) ;
	  c5->SaveAs(Form("%sfullmassrange_template_comparison_noreweight_%slog.png",dir.Data(),eta_q.Data())) ;
}

	RooDataHist* h_allpt = new RooDataHist("h_allpt","h_allpt",*roopt1);
//	h_allpt->add(*dset_data_allm_axis1);
	h_allpt->add(*dset_sideb_allm_axis1);
	TH1D* hallpt=(TH1D*)(h_allpt->createHistogram("hallpt",*roopt1));
 //   hallpt->SaveAs(Form("%sptdata_wholemassrange.root",dir.Data())) ;
    hallpt->SaveAs(Form("%sptsideb_wholemassrange.root",dir.Data())) ;
	  Double_t dpmq2[nq];  
  Double_t masserror[nq];
	  // ratioplot(dset_data_axis1, dset_rconeallm_axis1,kTRUE);
 if(massbinned){
	   
	 for(int i=0; i<nq ;i++)
	 {
		RooDataSet *dset_rcone_temp=new RooDataSet(*dset_rcone_allm_axis1,Form("dset_rcone%u_%u_axis1",i,startbin));
		RooDataSet *dset_sideb_temp=new RooDataSet(*dset_sideb_allm_axis1,Form("dset_sideb%u_%u_axis1",i,startbin));
    	RooDataSet *dset_massc = NULL;
    	RooDataSet *dset_sigmassc = NULL;
    	RooDataSet *dset_bkgmassc = NULL;
    	RooDataSet *dset_testsigmassc = NULL;
    	RooDataSet *dset_testbkgmassc = NULL;
		TString cut_s2= Form("roodiphomass>%f && roodiphomass<%f",dpmq[i],dpmq[i+1]);
		cut_s= Form("%1.0f_%2.0f",dpmq[i],dpmq[i+1]);
	 	cout << cut_s << endl;
  
		dset_massc = (RooDataSet*)(dset_data_allm_axis1->reduce(Name(Form("dset_massc%u_%u",i,startbin)),Cut(cut_s2)));
cout  << "entries overall data set " <<  dset_data_allm_axis1->sumEntries() << " "  << diphomass_data->GetEntries() <<  " entries of masscutted dataset " << dset_massc->sumEntries() << endl;
		dset_sigmassc = (RooDataSet*)(dset_rcone_allm_axis1->reduce(Name(Form("dset_sigmassc%u_%u",i,startbin)),Cut(cut_s2)));
		dset_bkgmassc = (RooDataSet*)(dset_sideb_allm_axis1->reduce(Name(Form("dset_bkgmassc%u_%u",i,startbin)),Cut(cut_s2)));
        if(check)
		{
			dset_testsigmassc = (RooDataSet*)(dset_testsig1->reduce(Name(Form("dset_testsigmassc%u_%u",i,startbin)),Cut(cut_s2)));
			dset_testbkgmassc = (RooDataSet*)(dset_testbkg1->reduce(Name(Form("dset_testbkgmassc%u_%u",i,startbin)),Cut(cut_s2)));
		}
	do1dfit(dset_massc,dset_sigmassc,dset_bkgmassc,dset_rcone_temp,dset_sideb_temp,dset_testsigmassc,dset_testbkgmassc,i, purity, purityerr, pullerror,data_entries);
    //  do1dfit(dset_massc,dset_sigmassc,dset_bkgmassc,dset_sigmassc,dset_bkgmassc,dset_testsigmassc,dset_testbkgmassc,i, purity, purityerr, pullerror,data_entries);
	
		cut_s.Clear();
		cut_s2.Clear();
		//dont want first 0 element
		dpmq2[i]     =(dpmq[i]+dpmq[i+1])/2.;
     	masserror[i]=(dpmq[i]-dpmq[i+1])/2.;

	 }
        
	   for(int k=0;k<nq;k++) 
	  {
	  cout << k << " entries  " << data_entries[k] << " purity " <<	purity[k] << " propagated error sumw2erroff " << purityerr[k] << " propagated error sumw2erron " <<	pullerror[k]  << endl;
	  }
	 plot_purity_massbin(dpmq2,masserror,purity,purityerr,false) ; 
	 if(!isdata_q)
	 {
	 plot_purity_massbin(dpmq2,masserror,purity,pullerror, true) ; 
	 }
 }
 else{
//no mass bining, one fit for all diphoton masses, just but high number a the end of title (10000)
   int i= 10000;
   do1dfit(dset_data_allm_axis1,dset_rcone_allm_axis1,dset_sideb_allm_axis1,dset_rcone_allm_axis1,dset_sideb_allm_axis1,dset_testsig1,dset_testbkg1,i, purity, purityerr,pullerror, data_entries);
 }
  delete mass_datapdf; 
  delete mass_data;
  delete diphomass_data;
   delete  fdata;delete frcone;delete fsideband; delete testfilesig; delete testfilebkg;
   return 0;
    
}


int do1dfit(RooDataSet *dset_data_axis1,RooDataSet * dset_sigcut,RooDataSet* dset_bkgcut,RooDataSet *dset_rconeallm_axis1,RooDataSet *dset_sideballm_axis1,RooDataSet* dset_testsig,RooDataSet*dset_testbkg,int im, double pu[], double puerr[],double pullerr[],double data_en[])
{
//reweighting
	cout << "massbin " << im << endl;
	dset_data_axis1->Print();
    data_en[im]=dset_data_axis1->sumEntries();
	dset_rconeallm_axis1->Print();
	dset_sideballm_axis1->Print();
   	cout << "ratio in massbin  truth/data weighted " << dset_rconeallm_axis1->sumEntries()/ dset_data_axis1->sumEntries()  << endl;
	RooDataSet* dset_sigorg=(RooDataSet*)dset_rconeallm_axis1->Clone("dset_sigorg");
	RooDataSet* dset_bkgorg=(RooDataSet*)dset_sideballm_axis1->Clone("dset_bkgorg");
	if(check) dset_testsig->Print();
	if(check) dset_testbkg->Print(); 
    if(check) { 
	  TCanvas *c3=new TCanvas("c3","isolation comparison");
	  c3->cd();
	  c3->SetLogy();
	  TLegend *leg3 = new TLegend(0.6,0.8,0.9,0.9);
	  leg3->SetFillColor(kWhite); 
	  RooPlot *testp = roovar1->frame(Title("comparison before reweighting"));
	  dset_sigcut->plotOn(testp,Rescale(1./dset_sigcut->sumEntries()),Binning(tbins),MarkerStyle(20),MarkerColor(kBlack),LineColor(kBlue),Name("rcone"));
	  dset_testsig->plotOn(testp,Rescale(1./dset_testsig->sumEntries()),Binning(tbins),MarkerStyle(20),MarkerColor(kBlue),LineColor(kBlue),Name("sigtruth"));
	  dset_bkgcut->plotOn(testp,Rescale(1./dset_bkgcut->sumEntries()),Binning(tbins),MarkerStyle(20),MarkerColor(kRed),LineColor(kRed),Name("sideband"));
	  dset_testbkg->plotOn(testp,Rescale(1./dset_testbkg->sumEntries()),Binning(tbins),MarkerStyle(20),MarkerColor(kMagenta),LineColor(kMagenta),Name("bkgtruth"));
	  testp->Draw();
	  testp->SetAxisRange(1e-3,20.,"Y");
	  if(!isdata_q)
      {
	  	leg3->AddEntry("rcone","random cone of massbin","l");
	   	leg3->AddEntry("sigtruth","sig truth of massbin","l");
	    leg3->AddEntry("bkgtruth","bkg truth of massbin","l");
	    leg3->AddEntry("sideband","sideband of massbin","l");
	  }
	  else if(isdata_q)
	  {
	  	leg3->AddEntry("rcone","random cone data","l");
	    leg3->AddEntry("sideband","sideband data","l");
	   	leg3->AddEntry("sigtruth","random cone mc","l");
	    leg3->AddEntry("bkgtruth","sideband mc","l");

	  }
	  leg3->Draw();
      c3->SaveAs(Form("%stemplate_comparison_noreweight_%s_%s_%u_%ulog.root",dir.Data(),(truthfit)? "truth": "rcone_sideb",cut_s.Data(),startbin_q,im)) ;
	  c3->SaveAs(Form("%stemplate_comparison_noreweight_%s_%s_%u_%ulog.png",dir.Data(),(truthfit)? "truth": "rcone_sideb",cut_s.Data(),startbin_q,im)) ;
	// 	cout << "ratio Entries truth/data weighted " << dset_rconeallm_axis1->sumEntries()/ dset_data_axis1->sumEntries()  << endl;
	}
	RooDataHist* h_pt = new RooDataHist("h_pt","h_pt",*roopt1);
	h_pt->add(*dset_data_axis1);
	TH1D* hpt=(TH1D*)(h_pt->createHistogram("hpt",*roopt1));
    hpt->SaveAs(Form("%shptdata%u%u.root",dir.Data(),startbin_q,im)) ;
	
/*
 	cout << __LINE__ << endl;	
	RooDataHist* hsig_eta = new RooDataHist("hsig_eta","hsig_eta",*rooeta1);
	hsig_eta->add(*dset_rconeallm_axis1);
	TH1D* heta=(TH1D*)(hsig_eta->createHistogram("heta",*rooeta1));
	RooDataHist* hsig_rho = new RooDataHist("hsig_rho","hsig_rho",*roorho);
	hsig_rho->add(*dset_rconeallm_axis1);
	TH1D* hrho=(TH1D*)(hsig_rho->createHistogram("hrho",*roorho));
	RooDataHist* hsig_sigma = new RooDataHist("hsig_sigma","hsig_sigma",*roosigma);
	hsig_sigma->add(*dset_rconeallm_axis1);
	TH1D* hsigma=(TH1D*)(hsig_sigma->createHistogram("hsigma",*roosigma));

	RooDataHist* hsig_pt = new RooDataHist("hsig_pt","hsig_pt",*roopt1);
	hsig_pt->add(*dset_rconeallm_axis1);
	TH1D* hpt=(TH1D*)(hsig_pt->createHistogram("hpt",*roopt1));
	RooDataHist* hsig_diphopt = new RooDataHist("hsig_diphopt","hsig_diphopt",*roodiphopt);
	hsig_diphopt->add(*dset_rconeallm_axis1);
	TH1D* hdiphopt=(TH1D*)(hsig_diphopt->createHistogram("hdiphopt",*roodiphopt));
    
   cout << __LINE__ << endl;	
	TCanvas *ccsig=new TCanvas("ccsig","random cone distributions before reweighting");
	ccsig->Divide(3,2);
	ccsig->cd(1);
    heta->Draw();
	ccsig->cd(2);
    hrho->Draw();
	ccsig->cd(3);
    hsigma->Draw();
	ccsig->cd(4);
    hpt->Draw();
	ccsig->cd(5);
    hdiphopt->Draw();
	
	ccsig->Draw();
    
   cout << __LINE__ << endl;	
	TString titlesig2= Form("sigtemplatecheck_%s_%s_sigtemplatecomp_%s_startb%u_%u",(EBEB)? "EBEB": "notEBEB",(truthfit)? "truth": "rcone_sideb",cut_s.Data(),startbin_q,im);
	ccsig->SaveAs(Form("%smassbinned20_4bins/%s.root",dir.Data(),titlesig2.Data()));
	ccsig->SaveAs(Form("%smassbinned20_4bins/%s.png",dir.Data(),titlesig2.Data()));  
	
   cout << __LINE__ << endl;	
	RooDataHist* hbkg_eta = new RooDataHist("hbkg_eta","hbkg_eta",*rooeta1);
	hbkg_eta->add(*dset_sideballm_axis1);
	TH1D* hetas=(TH1D*)(hbkg_eta->createHistogram("hetas",*rooeta1));
	RooDataHist* hbkg_rho = new RooDataHist("hbkg_rho","hbkg_rho",*roorho);
	hbkg_rho->add(*dset_sideballm_axis1);
	TH1D* hrhos=(TH1D*)(hbkg_rho->createHistogram("hrhos",*roorho));
	RooDataHist* hbkg_sigma = new RooDataHist("hbkg_sigma","hbkg_sigma",*roosigma);
	hbkg_sigma->add(*dset_sideballm_axis1);
	TH1D* hsigmas=(TH1D*)(hbkg_sigma->createHistogram("hsigmas",*roosigma));
	RooDataHist* hbkg_pt = new RooDataHist("hbkg_pt","hbkg_pt",*roopt1);
	hbkg_pt->add(*dset_sideballm_axis1);
	TH1D* hpts=(TH1D*)(hbkg_pt->createHistogram("hpts",*roopt1));
	RooDataHist* hbkg_diphopt = new RooDataHist("hbkg_diphopt","hbkg_diphopt",*roodiphopt);
	hbkg_diphopt->add(*dset_sideballm_axis1);
	TH1D* hdiphopts=(TH1D*)(hbkg_diphopt->createHistogram("hdiphopts",*roodiphopt));
    
	TCanvas *ccbkg=new TCanvas("ccbkg","sideband distributions before reweighting");
	ccbkg->Divide(3,2);
	ccbkg->cd(1);
    hetas->Draw();
	ccbkg->cd(2);
    hrhos->Draw();
	ccbkg->cd(3);
    hsigmas->Draw();
	ccbkg->cd(4);
    hpts->Draw();
	ccbkg->cd(5);
    hdiphopts->Draw();
    ccbkg->Draw();
    TString titlebkg2= Form("bkgtemplatecheck_%s_%s_sigtemplatecomp_%s_startb%u_%u",(EBEB)? "EBEB": "notEBEB",(truthfit)? "truth": "rcone_sideb",cut_s.Data(),startbin_q,im);
    ccbkg->SaveAs(Form("%smassbinned20_4bins/%s.root",dir.Data(),titlebkg2.Data()));
	ccbkg->SaveAs(Form("%smassbinned20_4bins/%s.png",dir.Data(),titlebkg2.Data()));  
*/

  	reweight_rhosigma(wrhor_eta,wrhor_pt,wrhor_diphopt,wrhorratiobin,wrhorratio,wsigmarratiobin,wsigmarratio,&dset_rconeallm_axis1,dset_data_axis1); 
  	reweight_rhosigma(wrho_eta,wrho_pt,wrho_diphopt,wrhoratiobin,wrhoratio,wsigmaratiobin,wsigmaratio,&dset_sideballm_axis1,dset_data_axis1); 
    reweight_eta_1d(weta_rho,weta_sigma,weta_pt,weta_diphopt,wetarratiobin,wetarratio,&dset_rconeallm_axis1,dset_data_axis1,1);
  	reweight_eta_1d(wetar_rho,wetar_sigma,wetar_pt,wetar_diphopt,wetaratiobin,wetaratio,&dset_sideballm_axis1,dset_data_axis1,1);

	reweight_pt_1d(wpt_rho,wpt_sigma,wpt_eta,wpt_diphopt,wptratiobin,wptratio,&dset_sideballm_axis1,dset_data_axis1,1); 
	reweight_diphotonpt_1d(wdiphopt_rho,wdiphopt_sigma,wdiphopt_eta,wdiphopt_pt,wdiphoptratiobin,wdiphoptratio,&dset_sideballm_axis1,dset_data_axis1); 
//no reweighting for pt if random cone applied

	
if(truthfit){
		  reweight_pt_1d(wptr_rho,wptr_sigma,wptr_eta,wptr_diphopt,wptrratiobin,wptrratio,&dset_rconeallm_axis1,dset_data_axis1,1);
		  reweight_diphotonpt_1d(wdiphoptr_rho,wdiphoptr_sigma,wdiphoptr_eta,wdiphoptr_pt,wdiphoptrratiobin,wdiphoptrratio,&dset_rconeallm_axis1,dset_data_axis1);
	}
								
//	if(check) { 
       if(isdata_q)
	   {
  	      reweight_rhosigma(twrhor_eta,twrhor_pt,twrhor_diphopt,twrhorratiobin,twrhorratio,twsigmarratiobin,twsigmarratio,&dset_testsig,dset_data_axis1); 
  	    reweight_rhosigma(twrho_eta,twrho_pt,twrho_diphopt,wrhoratiobin,wrhoratio,wsigmaratiobin,wsigmaratio,&dset_testbkg,dset_data_axis1); 
	 	  reweight_eta_1d(tweta_rho,tweta_sigma,tweta_pt,tweta_diphopt,twetaratiobin,twetaratio,&dset_testbkg,dset_data_axis1,1);
		  reweight_pt_1d(twptr_rho,twptr_sigma,twptr_eta,twptr_diphopt,twptrratiobin,twptrratio,&dset_testbkg,dset_data_axis1,1);
		  reweight_diphotonpt_1d(twdiphoptr_rho,twdiphoptr_sigma,twdiphoptr_eta,twdiphoptr_pt,twdiphoptrratiobin,twdiphoptrratio,&dset_testbkg,dset_data_axis1);
	 	  reweight_eta_1d(twetar_rho,twetar_sigma,twetar_pt,twetar_diphopt,twetarratiobin,twetarratio,&dset_testsig,dset_data_axis1,1);
	   }

		  TCanvas *c4=new TCanvas("c4","isolation comparison after reweighting");
		  c4->cd();
		  c4->SetLogy();
		  TLegend *leg4 = new TLegend(0.6,0.8,0.9,0.9);
		  leg4->SetFillColor(kWhite); 
		  RooPlot *testp2 = roovar1->frame(Title("comparison after reweighting"));
	      dset_rconeallm_axis1->plotOn(testp2,Rescale(1./dset_rconeallm_axis1->sumEntries()),Binning(tbins),MarkerStyle(20),MarkerColor(kBlack),LineColor(kBlack),Name("rcone1"));
	      dset_testsig->plotOn(testp2,Rescale(1./dset_testsig->sumEntries()),Binning(tbins),MarkerStyle(20),MarkerColor(kBlue),LineColor(kBlue),Name("sigtruth1"));
	      dset_sideballm_axis1->plotOn(testp2,Rescale(1./dset_sideballm_axis1->sumEntries()),Binning(tbins),MarkerStyle(20),MarkerColor(kRed),LineColor(kRed),Name("sideband1"));
	      dset_testbkg->plotOn(testp2,Rescale(1./dset_testbkg->sumEntries()),Binning(tbins),MarkerStyle(20),MarkerColor(kMagenta+2),LineColor(kMagenta+2),Name("bkgtruth1"));
		  testp2->Draw(); 
	      testp2->SetAxisRange(1e-3,20.,"Y");
		  if(!isdata_q)
		  {
		  	leg4->AddEntry("rcone1","random cone rew","l");
		  	leg4->AddEntry("sigtruth1","sig truth of massbin not rew","l");
			leg4->AddEntry("bkgtruth1","bkg truth of massbin not rew","l");
		  	leg4->AddEntry("sideband1","sideband rew","l");
		  }
		  else if(isdata_q)
		  {

		  	leg4->AddEntry("rcone1"," data random cone of massbin","l");
		  	leg4->AddEntry("sideband1","data sideband of massbin","l");
		  	leg4->AddEntry("sigtruth1","mc random cone of massbin","l");
			leg4->AddEntry("bkgtruth1","mc sideband of massbin","l");
		  }
		  leg4->Draw();
		  c4->SaveAs(Form("%stemplate_comparison_reweight_%s_%s_%u_%ulog.root",dir.Data(),(truthfit)? "truth": "rcone_sideb",cut_s.Data(),startbin_q,im)) ;
		  c4->SaveAs(Form("%stemplate_comparison_reweight_%s_%s_%u_%ulog.png",dir.Data(),(truthfit)? "truth": "rcone_sideb",cut_s.Data(),startbin_q,im)) ;
//	}	 

		  //check distributions, the one over whole massspectrum, in massbin and after reweighting for signal template
		 
		  Double_t sumentries=0; Double_t sumentries_cut=0; Double_t sumentries_rew=0; 
		  sumentries=dset_sigorg->sumEntries(); 
		  sumentries_cut=dset_sigcut->sumEntries(); 
		  sumentries_rew=dset_rconeallm_axis1->sumEntries(); 
		  cout << "sumentries " << sumentries << " sumentries_cut "  << sumentries_cut << " sumentries_rew " << sumentries_rew << endl;
		  TCanvas *csig=new TCanvas("csig","signal template comparison");
		  TLegend *leg1 = new TLegend(0.6,0.8,0.9,0.9);
		  leg1->SetFillColor(kWhite); 
		  csig->cd();
		   csig->SetLogy();
		  RooPlot *sig = roovar1->frame(Title("signal template comparison"));
		  dset_sigorg->plotOn(sig,Rescale(1./dset_sigorg->sumEntries()),Binning(tbins),MarkerStyle(20),MarkerColor(kBlack),LineColor(kBlack),Name("sigall"));
		  dset_sigcut->plotOn(sig,Rescale(1./dset_sigcut->sumEntries()),Binning(tbins),MarkerStyle(20),MarkerColor(kBlue),LineColor(kBlue),Name("sigmassbin"));
		  dset_rconeallm_axis1->plotOn(sig,Rescale(1./dset_rconeallm_axis1->sumEntries()),Binning(tbins),MarkerStyle(20),MarkerColor(kRed),LineColor(kRed),Name("sigrew"));
		  sig->Draw();
	      sig->SetAxisRange(1e-3,20.,"Y");
		  leg1->AddEntry("sigall","whole mass spectrum","l");
		  leg1->AddEntry("sigmassbin","in massbin","l");
		  leg1->AddEntry("sigrew","after reweight","l");
		  leg1->Draw();
		  TString titlesig;
		  titlesig= Form("sigtemplate_%s_%s_sigtemplatecomp_%s_startb%u_%u",(EBEB)? "EBEB": "nEBEB",(truthfit)? "truth": "rcone_sideb",cut_s.Data(),startbin_q,im);
		  //titlesig= Form("sigtemplate_%s_%s_sigtemplatecomp_startb%u_%u",(EBEB)? "EBEB": "nEBEB",(truthfit)? "truth": "rcone_sideb",startbin_q,im);
		  csig->SaveAs(Form("%s%s.root",dir.Data(),titlesig.Data()));
		  csig->SaveAs(Form("%s%s.png",dir.Data(),titlesig.Data())); 
		  //compare the background template to check if binning correct
		  Double_t sumentriesb=0; Double_t sumentries_cutb=0; Double_t sumentries_rewb=0; 
		  sumentriesb=dset_bkgorg->sumEntries(); 
		  sumentries_cutb=dset_bkgcut->sumEntries(); 
		  sumentries_rewb=dset_sideballm_axis1->sumEntries(); 
		  cout << "sumentries bkg " << sumentriesb << " sumentries_cut bkg "  << sumentries_cutb << " sumentries_rew bkg" << sumentries_rewb << endl;
		  TCanvas *cbkg=new TCanvas("cbkg","bkgnal template comparison");
		  cbkg->cd();
		  cbkg->SetLogy();
		  RooPlot *bkg = roovar1->frame(Title("bkgnal template comparison"));
		  dset_bkgorg->plotOn(bkg,Rescale(1./dset_bkgorg->sumEntries()),Binning(tbins),MarkerStyle(20),MarkerColor(kBlack),LineColor(kBlack),Name("bkgall"));
		  dset_bkgcut->plotOn(bkg,Rescale(1./dset_bkgcut->sumEntries()),Binning(tbins),MarkerStyle(20),MarkerColor(kBlue),LineColor(kBlue),Name("bkgmassbin"));
		  dset_sideballm_axis1->plotOn(bkg,Rescale(1./dset_sideballm_axis1->sumEntries()),Binning(tbins),MarkerStyle(20),MarkerColor(kRed),LineColor(kRed),Name("bkgrew"));
		  bkg->Draw();
	      bkg->SetAxisRange(1e-3,20.,"Y");
		  leg1->Draw();
		  TString titlebkg;
		  titlebkg= Form("bkg_template_%s_%s_bkgtemplatecomp_%s_startb%u_%u",(EBEB)? "EBEB": "nEBEB",(truthfit)? "truth": "rcone_sideb",cut_s.Data(),startbin_q,im);
		  //titlebkg= Form("bkg_template_%s_%s_bkgtemplatecomp_startb%u_%u",(EBEB)? "EBEB": "nEBEB",(truthfit)? "truth": "rcone_sideb",startbin_q,im);
		  cbkg->SaveAs(Form("%s%s.root",dir.Data(),titlebkg.Data()));
		  cbkg->SaveAs(Form("%s%s.png",dir.Data(),titlebkg.Data()));
//plots of reweighting proceduyre sidebands
		
        gStyle->SetOptStat(11111111);
	
		TCanvas * c_reweights=new TCanvas("c_reweights","c_reweights",2000,2000);
		c_reweights->Divide(3,2);  c_reweights->cd(1); 
		wptratiobin->Draw("colz");
    	c_reweights->cd(2);   wptratio->Draw("colz");
		c_reweights->cd(3);  wpt_rho->Draw("colz");
		c_reweights->cd(4);  wpt_sigma->Draw("colz");
		c_reweights->cd(5); wpt_eta->Draw("colz");
		c_reweights->cd(6); wpt_diphopt->Draw("colz");
		TString titlerew; 
		titlerew= Form("%sreweight_%s_%s_%s_startbin%u_%u",dir.Data(),(EBEB)? "EBEB": "nEBEB",(truthfit)? "truth": "rcone_sideb",cut_s.Data(),startbin_q,im);
		 c_reweights->SaveAs(Form("%s_ptbkg.root",titlerew.Data()));
		 c_reweights->SaveAs(Form("%s_ptbkg.png",titlerew.Data()));
		 
		TCanvas * c_reweights2=new TCanvas("c_reweights2","c_reweights2");
		c_reweights2->Divide(3,2);  c_reweights2->cd(1); 
		wdiphoptratiobin->Draw("colz");
    	c_reweights2->cd(2);   wdiphoptratio->Draw("colz");
		c_reweights2->cd(3);  wdiphopt_rho->Draw("colz");
		c_reweights2->cd(4);  wdiphopt_sigma->Draw("colz");
		c_reweights2->cd(5); wdiphopt_eta->Draw("colz");
		c_reweights2->cd(6); wdiphopt_pt->Draw("colz");
		c_reweights2->SaveAs(Form("%s_diphoptbkg.root",titlerew.Data()));
		c_reweights2->SaveAs(Form("%s_dphoptbkg.png",titlerew.Data()));

		TCanvas * c_rhoreweights=new TCanvas("c_rhoreweights","c_rhoreweights");
		c_rhoreweights->Divide(4,2);  c_rhoreweights->cd(1); 
		wrhoratiobin->Draw("colz");
    	c_rhoreweights->cd(2);   wrhoratio->Draw("colz");
		  c_rhoreweights->cd(3);  wsigmaratiobin->Draw("colz");  
		  c_rhoreweights->cd(4);  wsigmaratio->Draw("colz");
		c_rhoreweights->cd(5);  wrho_eta->Draw("colz");
		c_rhoreweights->cd(6);  wrho_pt->Draw("colz");
		c_rhoreweights->cd(7); wrho_diphopt->Draw("colz");
	  	  c_rhoreweights->SaveAs(Form("%s_rhosigmabkg.root",titlerew.Data()));
	  	  c_rhoreweights->SaveAs(Form("%s_rhosigmabkg.png",titlerew.Data()));
		  
		  TCanvas * c_etareweights=new TCanvas("c_etareweights","c_etareweights");
		  c_etareweights->Divide(3,2);  c_etareweights->cd(1);  wetaratiobin->Draw("colz");
	   	c_etareweights->cd(2);  wetaratio->Draw("colz");
		c_etareweights->cd(3); weta_rho->SetTitle("rhon/rhoo"); weta_rho->Draw("colz");
		c_etareweights->cd(4); weta_sigma->SetTitle("sigman/sigmao"); weta_sigma->Draw("colz");
		   c_etareweights->cd(5);    weta_pt->Draw("colz");
		  c_etareweights->cd(6);  weta_diphopt->Draw("colz"); 
	  	  c_etareweights->SaveAs(Form("%s_etabkg.root",titlerew.Data()));
	  	  c_etareweights->SaveAs(Form("%s_etabkg.png",titlerew.Data()));
//signal template
       if(truthfit)
	   {
  		  TCanvas * c_reweightsr=new TCanvas("c_reweightsr","c_reweightsr");
		c_reweightsr->Divide(3,2); 
    	c_reweightsr->cd(1);   wptrratiobin->Draw("colz");
		c_reweightsr->cd(2);  wptrratio->Draw("colz");
		c_reweightsr->cd(3); wptr_rho->Draw("colz");
		c_reweightsr->cd(4); wptr_sigma->Draw("colz");
		c_reweightsr->cd(5); wptr_eta->Draw("colz");
		c_reweightsr->cd(6); wptr_diphopt->Draw("colz");
	  	c_reweightsr->SaveAs(Form("%s_ptsig.root",titlerew.Data()));
	  	c_reweightsr->SaveAs(Form("%s_ptsig.png",titlerew.Data()));
		  
		TCanvas * c_reweightsr2=new TCanvas("c_reweightsr2","c_reweightsr2");
		c_reweightsr2->Divide(3,2); 
    	c_reweightsr2->cd(1);   wdiphoptrratiobin->Draw("colz");
		c_reweightsr2->cd(2);  wdiphoptrratio->Draw("colz");
		c_reweightsr2->cd(3); wdiphoptr_rho->Draw("colz");
		c_reweightsr2->cd(4); wdiphoptr_sigma->Draw("colz");
		c_reweightsr2->cd(5); wdiphoptr_eta->Draw("colz");
		c_reweightsr2->cd(6); wdiphoptr_pt->Draw("colz");
	  	c_reweightsr2->SaveAs(Form("%s_diphoptsig.root",titlerew.Data()));
	  	c_reweightsr2->SaveAs(Form("%s_diphoptsig.png",titlerew.Data()));
       }
		  
		 TCanvas * c_rhoreweightsr=new TCanvas("c_rhoreweightsr","c_rhoreweightsr");
		c_rhoreweightsr->Divide(4,2);  c_rhoreweightsr->cd(1); 
		wrhorratiobin->Draw("colz");
    	c_rhoreweightsr->cd(2);   wrhoratio->Draw("colz");
		  c_rhoreweightsr->cd(3);  wsigmarratiobin->Draw("colz");
		  c_rhoreweightsr->cd(4);  wsigmarratio->Draw("colz");
		c_rhoreweightsr->cd(5);  wrhor_eta->Draw("colz");
		c_rhoreweightsr->cd(6);  wrhor_pt->Draw("colz");
		c_rhoreweightsr->cd(7); wrhor_diphopt->Draw("colz");
	  	  c_rhoreweightsr->SaveAs(Form("%s_rhosigmasig.root",titlerew.Data()));
	  	  c_rhoreweightsr->SaveAs(Form("%s_rhosigmasig.png",titlerew.Data()));
		  TCanvas * c_etareweightsr=new TCanvas("c_etareweightsr","c_etareweightsr");
		  c_etareweightsr->Divide(3,2);  c_etareweightsr->cd(1);  wetarratiobin->Draw("colz");
	   	c_etareweightsr->cd(2);  wetarratio->Draw("colz");
		c_etareweightsr->cd(3); wetar_rho->SetTitle("rhon/rhoo"); wetar_rho->Draw("colz");
		c_etareweightsr->cd(4); wetar_sigma->SetTitle("sigman/sigmao"); wetar_sigma->Draw("colz");
		   c_etareweightsr->cd(5);    wetar_pt->Draw("colz");
		  c_etareweightsr->cd(6);  wetar_diphopt->Draw("colz"); 
	  	  c_etareweightsr->SaveAs(Form("%s_etasig.root",titlerew.Data()));
	    	c_etareweightsr->SaveAs(Form("%s_etasig.png",titlerew.Data()));
 
  if(check) {
	  cout << "after reweighting" << endl; 
	  cout << "sig " << dset_rconeallm_axis1->sumEntries() << endl;
	  cout << "bkg " << dset_sideballm_axis1->sumEntries() << endl;
	  cout << "data " << dset_data_axis1->sumEntries() << endl;
	  cout << "signal+ fakes " << dset_rconeallm_axis1->sumEntries()+ dset_sideballm_axis1->sumEntries()<< endl;
	  cout << "one? " << (dset_rconeallm_axis1->sumEntries()+ dset_sideballm_axis1->sumEntries())/dset_data_axis1->sumEntries()<< endl;
	  cout << "ratio Entries truth/data weighted " << dset_rconeallm_axis1->sumEntries()/dset_data_axis1->sumEntries()  << endl;
	  }
//fitting
	  cout << "fitting starts" << endl;
	//first step determine single-photon purities on Iso1 and Iso2 axes
	//for integral 0-90 GeV, fit should not depend on initial values!  
	//give initial values for purity fraction -j1 overall purity for leg 1 
	  RooRealVar *j1=NULL;RooFormulaVar *fsig1=NULL;
	  RooRealVar *j2=NULL;RooFormulaVar *fsig2=NULL;
	  j1 = new RooRealVar("j1","j1",0.5,0,1);
	  j2 = new RooRealVar("j2","j2",0.3,0,1);
	  fsig1 = new RooFormulaVar("fsig1","fsig1","j1",RooArgList(*j1));
	  fsig2 = new RooFormulaVar("fsig2","fsig2","j2",RooArgList(*j2));
	  
	  //produces binned RooDataHist also when RooRealVar input and weights are unbinned
	  RooArgSet var;
	  if(dohggv){var.add(*roovar1);}
	  else if(dovzero) {var.add(*rooisopv1);}
	  else if(dowv) {var.add(*rooisowv1);}
	  RooDataHist *drconehist_axis1=NULL;
	  RooDataHist *dsidebhist_axis1=NULL;
	  RooHistPdf *drconepdf_axis1 =NULL;
	  RooHistPdf *dsidebpdf_axis1 =NULL;
	  drconehist_axis1 = new RooDataHist("drconehist_axis1","drconehist_axis1",var);
	  drconehist_axis1->add(*dset_rconeallm_axis1);
	  if(debug){
	  	drconehist_axis1->Print("V");
	  }
	  dsidebhist_axis1 = new RooDataHist("dsidebhist_axis1","dsidebhist_axis1",var);
	  dsidebhist_axis1->add(*dset_sideballm_axis1);
	  if(debug){
		  dsidebhist_axis1->Print("V"); 
      }
      drconepdf_axis1 = new RooHistPdf("drconepdf_axis1","drconepdf_axis1",var,*drconehist_axis1);
	  dsidebpdf_axis1 = new RooHistPdf("dsidebpdf_axis1","dsidebpdf_axis1",var,*dsidebhist_axis1);
	  //normalizes histograms over all observables
	 // RooAddPdf is an efficient implementation of a sum of PDFs of the form  c_1*PDF_1 + c_2*PDF_2 + ... (1-sum(c_1...c_n-1))*PDF_n
	  //the sum of the coefficients is enforced to be one,and the coefficient of the last PDF is calculated from that condition.
	  //ArgList of coefficients -1 of roohistpdf
	  RooAddPdf *model_axis1=NULL;
	  // do ML fit of pdf to unbinned data
	  // weights are supported in unbinned dataset
	  model_axis1 = new RooAddPdf("model_axis1","model_axis1",RooArgList(*drconepdf_axis1,*dsidebpdf_axis1),RooArgList(*fsig1),kFALSE);
      TCanvas* c_deb=new TCanvas("c_deb","c_deb");
      c_deb->cd();
      RooPlot *cdeb;
      cdeb = roovar1->frame(Title("compare pdfs before fit"));
	  drconepdf_axis1->plotOn(cdeb,Components("drconepdf_axis1"),LineStyle(kDashed),LineColor(kRed),Name("signalpdf"));
      dsidebpdf_axis1->plotOn(cdeb,Components("dsidebpdf_axis1"),LineStyle(kDashed),LineColor(kBlack),Name("bkgpdf"));
	  model_axis1->plotOn(cdeb,Components("drconepdf_axis1"),LineColor(kRed-2),Name("signal"));
	  model_axis1->plotOn(cdeb,Components("dsidebpdf_axis1"),LineColor(kBlue+2),Name("background"));
      cdeb->Draw();
	  c_deb->SaveAs(Form("%spdfcomp_%s.png",dir.Data(),cut_s.Data()));c_deb->SaveAs(Form("%spdfcomp_%s.root",dir.Data(),cut_s.Data()));
	  //if does not converge take NumCPU strategy 2
	  RooFitResult *firstpass1;
	  TFile resFile(Form("%sfitresult_%s_%s_sumw2erroff_%u_%s_range_%u_%u.root",dir.Data(),(truthfit)? "truth": "rcone_sideb",(EBEB) ? "EBEB": "nEBEB",im,cut_s.Data(),startbin_q,endbin_q),"RECREATE"); 
	  //second fit with sumw2error on for PULL fct= (ptrue-pmeas)/sigma_meas
	  //if does not converge take NumCPU strategy 2
	 //ML fit to weighted dataset: SumW2Error takes statistics of dataset into account, scales with number of events in dataset
	 // By default the fit is executed through the MINUIT commands MINOS, HESSE and MINOS in succession 
	 //fit to unbinned! data (roodataset)
	  cout << "firstpass1 sumw2error false" << endl;
	  resFile.cd();
	  firstpass1 = model_axis1->fitTo(*dset_data_axis1, NumCPU(8), Extended(false),SumW2Error(kFALSE),Verbose(kFALSE),Save(kTRUE));
	  firstpass1->Print();
	  firstpass1->Write();
	  resFile.Write();  
	  resFile.Close();
      RooAddPdf *model2_axis1=NULL;
      RooFitResult *firstpass2=NULL;
	  if(!isdata_q)
	  {
		  //second fit with sumw2error on for PULL fct= (ptrue-pmeas)/sigma_meas
		  model2_axis1 = new RooAddPdf("model2_axis1","model2_axis1",RooArgList(*drconepdf_axis1,*dsidebpdf_axis1),RooArgList(*fsig2),kFALSE);
		  TFile resFile2(Form("%sfitresult_%s_%s_sumw2erron_%u_%s_range_%u_%u.root",dir.Data(),(truthfit)? "truth": "rcone_sideb",(EBEB) ? "EBEB": "nEBEB",im,cut_s.Data(),startbin_q,endbin_q),"RECREATE"); 
		  cout << "firstpass2 sumw2error true " << endl;
		  resFile2.cd();
		  firstpass2 = model2_axis1->fitTo(*dset_data_axis1, NumCPU(8), Extended(false),SumW2Error(kTRUE),Verbose(kFALSE),Save(kTRUE));
		  firstpass2->Print();
		  firstpass2->Write();
		  resFile2.Write();  
		  resFile2.Close();

		  TLatex b2;b2.SetNDC();b2.SetTextSize(0.06);b2.SetTextColor(kRed);
		  TCanvas *c12 = new TCanvas("c12","c12",1200,800);
		  TLegend *leg2 = new TLegend(0.15,0.8,0.35,0.9);
		  c12->cd(1);
		  TString title2; 
		  if(dohggv){
			title2= Form("1d fit for hgg vertex SumW2Error, %s, reweighted,%u mass: %s",(EBEB)? "EBEB": "nEBEB" ,im,cut_s.Data());
		  }
	   //lin plot
		  RooPlot *frame1bla2;
		  frame1bla2 = roovar1->frame(Title(title2.Data()));
		  dset_data_axis1->plotOn(frame1bla2,Binning(tbins),Name("data"));
		  model2_axis1->plotOn(frame1bla2,Name("fit"));
		  model2_axis1->plotOn(frame1bla2,Components("drconepdf_axis1"),LineStyle(kDashed),LineColor(kRed),Name("signal"));
		  model2_axis1->plotOn(frame1bla2,Components("dsidebpdf_axis1"),LineStyle(kDashed),LineColor(kBlack),Name("background"));
		  frame1bla2->Draw();
		  leg2->AddEntry("fit","fit","l");
		  if(!truthfit){
			leg2->AddEntry("signal","rcone MC","l");
			leg2->AddEntry("background","sideband  MC","l");
		  }
		  else if(truthfit){
			  leg2->AddEntry("signal","signal MC","l");
			  leg2->AddEntry("background","background  MC","l");
		  }
			  leg2->SetFillColor(kWhite); 
		  leg2->Draw();
		  b2.DrawLatex(0.55,0.7,"PRELIMINARY");

//log plot
		  TCanvas *c22 = new TCanvas("c22","c22",1200,800);
		  c22->cd(1);
		  gPad->SetLogy();
		  RooPlot *frame1logbla2;
		  frame1logbla2 = roovar1->frame(Title(title2.Data()));
		  dset_data_axis1->plotOn(frame1logbla2,Binning(tbins),Name("data"));
		  model2_axis1->plotOn(frame1logbla2,Name("fit"));
		  model2_axis1->plotOn(frame1logbla2,Components("drconepdf_axis1"),LineStyle(kDashed),LineColor(kRed),Name("signal"));
		  model2_axis1->plotOn(frame1logbla2,Components("dsidebpdf_axis1"),LineStyle(kDashed),LineColor(kBlack),Name("background"));
		  frame1logbla2->Draw();
		  leg2->Draw();
		  b2.DrawLatex(0.55,0.7,"PRELIMINARY");
		  TString titleout2;
		  if(dohggv) titleout2= Form("%s_reweight_%s_hggvtx_SUMW2_%s_startb_%u_%u_sigcomp",(EBEB)? "EBEB": "nEBEB",(truthfit)? "truth": "rcone_sideb",cut_s.Data(),startbin_q,im);
		  const char* outfile2=Form("%sroofit_%s_%s",dir.Data(), (isdata_q)? "data": "mc",titleout2.Data());
		  c12->SaveAs(Form("%s_lin.png",outfile2));c12->SaveAs(Form("%s_lin.root",outfile2));c12->SaveAs(Form("%s_lin.pdf",outfile2));
		  c22->SaveAs(Form("%s_log.png",outfile2));c22->SaveAs(Form("%s_log.root",outfile2));c22->SaveAs(Form("%s_log.pdf",outfile2));
      }
	  if(massbinned)
	  {
	  	pu[im]=fsig1->getVal();
		puerr[im]=fsig1->getPropagatedError(*firstpass1);
	 	if(!isdata_q)
		{	
		pullerr[im]=fsig2->getPropagatedError(*firstpass2);
		}
	  }
	  TLatex b;b.SetNDC();b.SetTextSize(0.06);b.SetTextColor(kRed);
	  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
	  TLegend *leg = new TLegend(0.15,0.8,0.35,0.9);
	  c1->cd(1);
	  TString title; 
	  if(dohggv){
	  //	title= Form("1d fit for hgg vertex no SumW2Error, %s, reweighted,%u rangebins: %u -  %u",(EBEB)? "EBEB": "nEBEB" ,im, startbin_q, endbin_q);
	  	title= Form("1d fit for hgg vertex no SumW2Error, %s, reweighted,%u mass: %s",(EBEB)? "EBEB": "nEBEB" ,im,cut_s.Data());
      }
	  else if(dovzero){
	  	title= Form("1d fit for vertex[0], %s, reweighted,%u",(EBEB)? "EBEB": "nEBEB" ,im);
	  }
	  else if(dowv){
	  	title= Form("1d fit for worst isolation, %s, reweighted,%u",(EBEB)? "EBEB": "nEBEB" ,im);
	  }
   //lin plot
  	  RooPlot *frame1bla;
 	  frame1bla = roovar1->frame(Title(title.Data()));
	  dset_data_axis1->plotOn(frame1bla,Binning(tbins),Name("data"));
	  model_axis1->plotOn(frame1bla,Name("fit"));
	  model_axis1->plotOn(frame1bla,Components("drconepdf_axis1"),LineStyle(kDashed),LineColor(kRed),Name("signal"));
	  model_axis1->plotOn(frame1bla,Components("dsidebpdf_axis1"),LineStyle(kDashed),LineColor(kBlack),Name("background"));
	  frame1bla->Draw();
	  leg->AddEntry("fit","fit","l");
	  if(!truthfit){
	  	leg->AddEntry("signal","rcone MC","l");
	  	leg->AddEntry("background","sideband  MC","l");
	  }
	  else if(truthfit){
		  leg->AddEntry("signal","signal MC","l");
	      leg->AddEntry("background","background  MC","l");
	  }
		  leg->SetFillColor(kWhite); 
	  leg->Draw();
	  b.DrawLatex(0.55,0.7,"PRELIMINARY");

//log plot
	  TCanvas *c2 = new TCanvas("c2","c2",1200,800);
	  c2->cd(1);
	  gPad->SetLogy();
	  RooPlot *frame1logbla;
	  frame1logbla = roovar1->frame(Title(title.Data()));
	  dset_data_axis1->plotOn(frame1logbla,Binning(tbins),Name("data"));
	  model_axis1->plotOn(frame1logbla,Name("fit"));
	  model_axis1->plotOn(frame1logbla,Components("drconepdf_axis1"),LineStyle(kDashed),LineColor(kRed),Name("signal"));
	  model_axis1->plotOn(frame1logbla,Components("dsidebpdf_axis1"),LineStyle(kDashed),LineColor(kBlack),Name("background"));
	  frame1logbla->Draw();
	  leg->Draw();
	  b.DrawLatex(0.55,0.7,"PRELIMINARY");
//	  model_axis1->Print();
	 // dset_data_axis1->Print();
	  //save fit plots
	  TString titleout;
	  if(dohggv) titleout= Form("%s_reweight_%s_hggvtx_noSUMW2_%s_startb_%u_%u_sigcomp",(EBEB)? "EBEB": "nEBEB",(truthfit)? "truth": "rcone_sideb",cut_s.Data(),startbin_q,im);
	  else if(dovzero) titleout= "nooverflow_reweight_fulletarange_rcone_sideband_vzero";
	  else if(dowv) titleout= "nooverflow_reweight_fulletarange_rcone_sideband_wv";
	  const char* outfile=Form("%sroofit_%s_%s",dir.Data(), (isdata_q)? "data": "mc",titleout.Data());
	  c1->SaveAs(Form("%s_lin.png",outfile));c1->SaveAs(Form("%s_lin.root",outfile));c1->SaveAs(Form("%s_lin.pdf",outfile));
	  c2->SaveAs(Form("%s_log.png",outfile));c2->SaveAs(Form("%s_log.root",outfile));c2->SaveAs(Form("%s_log.pdf",outfile));
	delete dset_data_axis1;
	delete dset_rconeallm_axis1;
	delete dset_sideballm_axis1;
	delete dset_sigorg;
	delete dset_bkgorg;
	delete dset_testbkg;
	delete drconehist_axis1;
	delete dsidebhist_axis1;    
	delete drconepdf_axis1;  
	delete dsidebpdf_axis1;
	dset_data_axis1=NULL;
	 dset_rconeallm_axis1=NULL;
	 dset_sideballm_axis1=NULL;
	 dset_sigorg=NULL;
	 dset_bkgorg=NULL;
	 dset_testbkg=NULL;
	 drconehist_axis1=NULL;
	 dsidebhist_axis1=NULL;    
	 drconepdf_axis1=NULL;  

 	 delete c1;
	delete c2;
	delete csig;
	delete cbkg;
	
	delete c_reweights;
	delete c_reweights2;
	delete c_rhoreweights;
	delete c_etareweights;

//	delete c_reweightsr;
//	delete c_reweightsr2;
	delete c_rhoreweightsr;
	delete c_etareweightsr;

    delete model_axis1;	
    delete model2_axis1;
/*
	if(check){
		delete w_rhornr; delete w_rhor2o; delete w_rhorn; delete w_rhor2; delete w_rhor; delete w_rhoro;
		delete w_rhonr; delete w_rho2o; delete w_rhon; delete w_rho2; delete w_rho; delete w_rhoo;
		delete w_sigmarnr; delete w_sigmar2o; delete w_sigmarn; delete w_sigmar2;
		delete w_sigmanr; delete w_sigma2o; delete w_sigman; delete w_sigma2; 
		delete w_ptrnr; delete w_ptr2o; delete w_ptrn; delete w_ptr2; delete w_ptr; delete w_ptro;
		delete w_ptnr; delete w_pt2o; delete w_ptn; delete w_pt2; delete w_pt; delete w_pto;
		delete w_etarnr; delete w_etar2o; delete w_etarn; delete w_etar2; delete w_etar; delete w_etaro;
		delete w_etanr; delete w_eta2o; delete w_etan; delete w_eta2; delete w_eta; delete w_etao;
		delete w_diphoptrnr; delete w_diphoptr2o; delete w_diphoptrn; delete w_diphoptr2; delete w_diphoptr; delete w_diphoptro;
		delete w_diphoptnr; delete w_diphopt2o; delete w_diphoptn; delete w_diphopt2; delete w_diphopt; delete w_diphopto;
		
		delete tw_rhornr; delete tw_rhor2o; delete tw_rhorn; delete tw_rhor2; delete tw_rhor; delete tw_rhoro;
		delete tw_rhonr; delete tw_rho2o; delete tw_rhon; delete tw_rho2; delete tw_rho; delete tw_rhoo;
		delete tw_sigmarnr; delete tw_sigmar2o; delete tw_sigmarn; delete tw_sigmar2;
		delete tw_sigmanr; delete tw_sigma2o; delete tw_sigman; delete tw_sigma2; 
		delete tw_ptrnr; delete tw_ptr2o; delete tw_ptrn; delete tw_ptr2; delete tw_ptr; delete tw_ptro;
		delete tw_ptnr; delete tw_pt2o; delete tw_ptn; delete tw_pt2; delete tw_pt; delete tw_pto;
		delete tw_etarnr; delete tw_etar2o; delete tw_etarn; delete tw_etar2; delete tw_etar; delete tw_etaro;
		delete tw_etanr; delete tw_eta2o; delete tw_etan; delete tw_eta2; delete tw_eta; delete tw_etao;
		delete tw_diphoptrnr; delete tw_diphoptr2o; delete tw_diphoptrn; delete tw_diphoptr2; delete tw_diphoptr; delete tw_diphoptro;
	}
*/
	print_mem();
	return 0;

}


