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

using namespace std; 
using namespace RooFit;


//check how much memory the jobs use, 4 GB normal
 void print_mem(){
	 gSystem->GetProcInfo(&procinfo);
	 std::cout << "Resident mem (kB): " << procinfo.fMemResident << std::endl;
	 std::cout << "Virtual mem (kB):  " << procinfo.fMemVirtual << std::endl;
	 gSystem->Sleep((UInt_t)1e3);
	  };

//computes mass quantiles for diphoton mass binning and plots the result
int quantiles(TH1F* mass,double probtemp[], double dpmqtemp[] ) {
		cout << "define mass bins " << endl;
        mass->Scale(1.0/mass->Integral());     
        for (Int_t i =0;i<=nq_q;i++){ 
        //probtemp[i] =  Float_t(i+startbin_q)/ntot_q;
        
		//first entry e.g. 0.05 - diphopt 73 to end up with 10 bins in the end	
	    probtemp[i] =  Float_t(i+startbin_q)/ntot_q;
        mass->GetQuantiles(nq_q,dpmqtemp,probtemp);
		cout << "quantile  "<<probtemp[i] << " diphopt  " << dpmqtemp[i]  << endl;
		}
        //show the original histogram in the top pad
//        if(debug){ 
			TCanvas *cq = new TCanvas("cq","mass quantiles",10,10,700,900);
			cq->Divide(1,2);
			cq->cd(1);
			cq->SetLogy();
			cq->Update();
			mass->Draw();

			// show the quantiles in the bottom pad
			cq->cd(2);
			gPad->SetGrid();
			TGraph *gr = new TGraph(nq_q,probtemp,dpmqtemp);
			gr->SetMarkerStyle(21);
			gr->Draw("alp");
			cq->SaveAs(Form("../plots/March18/massbinned20_4bins/massquantiles_%s_range_%u_%u_%s.root",(truthfit)? "truth": "rcone_sideb",startbin_q,endbin_q,eta_q.Data()));
			cq->SaveAs(Form("../plots/March18/massbinned20_4bins/massquantiles_%s_range_%u_%u_%s.png",(truthfit)? "truth": "rcone_sideb",startbin_q,endbin_q,eta_q.Data()));
//		}
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
	
    //for(int i=0;i< nq;i++) {
//		mass
//	}

	TGraphErrors *gr = new TGraphErrors(nq_q,mass, pur,masserr,purerr);
 //	TGraphAsymmErrors *gr = new TGraphAsymmErrors(nq_q,mass, pur,masserrl,masserrr,purerr,purerr);
//	cgr->SetLogx();
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
    gr->SaveAs(Form("../plots/March18/massbinned20_4bins/150318_%s_purity_%s_%s_massbin_%u_range_%u_%u.root",(truthfit) ? "truth": "rcone_sideb",(berr)? "sumw2erron": "sumw2erroff",eta_q.Data(),ntot_q,startbin_q,endbin_q));
    cgr->Print(Form("../plots/March18/massbinned20_4bins/150318_%s_purity_%s_%s_massbin_%u_range_%u_%u.png",(truthfit) ? "truth": "rcone_sideb",(berr)? "sumw2erron": "sumw2erroff",eta_q.Data(),ntot_q,startbin_q,endbin_q),"png");
    
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
  float v_roodiphomass;
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
  TBranch *b_roodiphomass;
  TBranch *b_roorho;
  TBranch *b_roosigma;
  TBranch *b_roonvtx;
  TBranch *b_rooweight;

  const int nvars = 18;
  float* ptrs[nvars]={&v_roovar1,&v_roovar2,&v_rooisopv1,&v_rooisopv2,&v_rooisowv1,&v_rooisowv2,&v_roopt1,&v_roosieie1,&v_rooeta1,&v_roopt2,&v_roosieie2,&v_rooeta2,&v_roodiphopt,&v_roodiphomass,&v_roorho,&v_roosigma,&v_roonvtx,&v_rooweight};
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
void reweight_rhosigma(TH1F* weight_rho, TH1F* weight_rhoo,TH2F*weight_rhon,TH2F*weight_rho2o,TH2F* weight_rhonr, TH2F* weight_rho2,TH2F*weight_sigman,TH2F*weight_sigma2o,TH2F* weight_sigmanr, TH2F* weight_sigma2,RooDataSet **dset, RooDataSet *dsetdestination, bool deleteold){
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
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float rho = args.getRealValue("roorho");
    float sigma = args.getRealValue("roosigma");
    float neww = oldw*h->GetBinContent(h->FindBin(rho,sigma));
	//if(debug){
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
//	}
		newdset->add(args,neww);
	  }
	  newdset->SetName((*dset)->GetName());
	  newdset->SetTitle((*dset)->GetTitle());
	  delete hnum; delete hden;
	  RooDataSet *old_dset = *dset;
	  *dset=newdset;
	  std::cout << "RhoSigma2D rew: norm from " << old_dset->sumEntries() << " to " << newdset->sumEntries() << std::endl;
	  if (deleteold) delete old_dset;
	};


void reweight_pt_1d(TH1F* weight_pt, TH1F* weight_pto,TH2F*weight_ptn,TH2F*weight_pt2o,TH2F* weight_ptnr, TH2F* weight_pt2,RooDataSet **dset, RooDataSet *dsetdestination, int numvar){

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
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float pt = args.getRealValue(ptname);
    float neww = oldw*h->GetBinContent(h->FindBin(fabs(pt)));
  //  if(debug){
		weight_pt->Fill(neww);	
		weight_pto->Fill(oldw);
		weight_pt2o->Fill(h->FindBin(fabs(pt)),oldw);
		weight_ptn->Fill(h->FindBin(fabs(pt)),neww);	
		if(oldw!=0 && neww!=0)weight_ptnr->Fill(h->FindBin(fabs(pt)),oldw/neww);
		else {weight_ptnr->Fill(-10,1);}
		if(oldw!=0 && neww!=0)weight_pt2->Fill(pt,oldw/neww);
		else {weight_pt2->Fill(-10,1);}
    // }
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
   // if(debug){
		weight_diphopt->Fill(neww);
		weight_diphopto->Fill(oldw);
		weight_diphopt2o->Fill((h->FindBin(fabs(diphopt))),oldw);	
		weight_diphoptn->Fill((h->FindBin(fabs(diphopt))),neww);	
		if(oldw!=0)weight_diphoptnr->Fill(h->FindBin(fabs(diphopt)),oldw/neww);
		else {weight_diphoptnr->Fill(-10,1);}//cout << "dipho weight old 0" << endl;}
		if(oldw!=0)weight_diphopt2->Fill(diphopt,oldw/neww);
		else {weight_diphopt2->Fill(-10,1);}//cout << "dipho weight old 0" << endl;}
	 //}
		}

  newdset->SetName((*dset)->GetName());
  newdset->SetTitle((*dset)->GetTitle());
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

  hnum->Divide(hden);
  TH1F *h = hnum;

  RooDataSet *newdset = new RooDataSet(**dset,Form("%s_etarew",(*dset)->GetName()));
  newdset->reset();
  for (int i=0; i<(*dset)->numEntries(); i++){
    RooArgSet args = *((*dset)->get(i));
    float oldw = (*dset)->store()->weight(i);
    float eta = args.getRealValue(etaname);
    float neww = oldw*h->GetBinContent(h->FindBin(fabs(eta)));
  //  if(debug){ 
		weight_eta->Fill(neww);	
		weight_etao->Fill(oldw);
		weight_etan->Fill(h->FindBin(fabs(eta)),neww);	
		weight_eta2o->Fill(h->FindBin(fabs(eta)),oldw);
		if(oldw!=0 && neww!=0)weight_etanr->Fill(h->FindBin(fabs(eta)),oldw/neww);
		else {weight_etanr->Fill(-10,1);}
	   // weight_pt2->Fill(pt,neww/oldw);
	   if(oldw!=0 && neww!=0)weight_eta2->Fill(fabs(eta),oldw/neww);
	   else {weight_eta2->Fill(-10,1);}
//	}  
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
int prep_dataset(TString eta,Bool_t truthf, int ntotbins, int startbin, int endbin){
 //template from MC truth or with random cone/sideband?
  truthfit=truthf;
  // EBEB or nEBEB
  eta_q=eta;
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
  cout << "ntot bins" << ntotbins << endl;
  cout << "massbins " << nq << " startbin " << startbin << " endbin " << endbin << endl;
  
//100 mass bins 10-15 % statistics on the purity fit are sufficient
Double_t prob[nq];  // position where to compute the quantiles in [0,1]
Double_t dpmq[nq];  // array to contain the quantiles
Double_t purity[nq];  // array to contain the quantiles
Double_t purityerr[nq];  // array to contain the quantiles
Double_t pullerror[nq];  // array to contain the quantiles
//Double_t masserror[nq];  // array to contain the quantiles

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
 // 	tbins.addUniform(1, 0.,0.1);
//  	tbins.addUniform(2, 0.1,2);
//  	tbins.addUniform(3, 2,9);
  	tbins.addUniform(1,0.,0.1);
 	tbins.addUniform(1,0.1,1.);
  	tbins.addUniform(2,1.,9.);
  }
  //only one diphoton mass-binning a bit finer
  if (!massbinned){
    tbins.addUniform(5,0.,.5) ;
  	tbins.addUniform(3,0.5,1.) ;
  	tbins.addUniform(9,1.,3.) ;
	tbins.addUniform(3,3.,6.) ;
	tbins.addUniform(3,6.,9.);
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
  Float_t ptrightrange=1000.;
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
  roodiphomass = new RooRealVar("roodiphomass","roodiphomass",0.,4000.);
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
  TH1::SetDefaultSumw2(kTRUE);
  //Axis 1 and 2 randomly swapped
  TFile* fdata=NULL; TFile* frcone=NULL; TFile * fsideband=NULL; TFile* testfile=NULL;
  RooDataSet* dset_test=NULL;
//  fdata = new TFile("../forroofit/Feb27/data.root","read"); 
  fdata = new TFile("../forroofit/March09/data.root","read"); 
  get_roodset_from_ttree(fdata,"for_roofit",dataset_data);
//both leading and subleading photon 
  if(!truthfit){ 
//  	frcone = new TFile("../forroofit/Feb27/rcone.root","read"); 	
  	frcone = new TFile("../forroofit/March09/rcone.root","read"); 	
  	get_roodset_from_ttree(frcone,"for_roofit",dataset_rcone);
  	//fsideband= new TFile("../forroofit/Feb27/sideband.root","read");
  	fsideband= new TFile("../forroofit/March09/sideband.root","read");
  	get_roodset_from_ttree(fsideband,"for_roofit",dataset_sideband);
  	if(check) {	
	testfile= new TFile("../forroofit/Feb27/sideband.root","read");
		get_roodset_from_ttree(testfile,"for_roofit",dataset_test1);
  	}
  }
  if(truthfit){ 
 // 	frcone = new TFile("../forroofit/Feb27/sig.root","read"); 	
  	frcone = new TFile("../forroofit/March09/sig.root","read"); 	
  	get_roodset_from_ttree(frcone,"for_roofit",dataset_rcone);
  //	fsideband= new TFile("../forroofit/Feb27/bkg.root","read");
    fsideband= new TFile("../forroofit/March09/bkg.root","read");
  	get_roodset_from_ttree(fsideband,"for_roofit",dataset_sideband);
  	if(check) {	
  		testfile= new TFile("../forroofit/Feb27/bkg.root","read");
		get_roodset_from_ttree(testfile,"for_roofit",dataset_test1);
  	}
   }
   assert(fdata); assert(frcone);  assert(fsideband);
  



//define mass bins
  TH1F* diphomass_data=(TH1F*)fdata->Get("h_diphomass") ;
  assert(diphomass_data);
  quantiles(diphomass_data,prob,dpmq); 
  TString cut_string; 
  if(!EBEB){
   	cut_string= Form("roovar1<%f && (rooeta1 > 1.5 || rooeta2 > 1.5)",rightrange-1e-5);
  }
  else if(EBEB){
	cut_string= Form("roovar1<%f && (rooeta1<1.5 && rooeta2 < 1.5)",rightrange-1e-5);
  }
  if (dataset_data) dset_data = (RooDataSet*)(dataset_data->reduce(Name("dset_data"),Cut(cut_string)));
  if (dataset_rcone) dset_rcone = (RooDataSet*)(dataset_rcone->reduce(Name("dset_rcone"),Cut(cut_string)));
  if (dataset_sideband) dset_sideb = (RooDataSet*)(dataset_sideband->reduce(Name("dset_sideb"),Cut(cut_string)));
  if (dataset_test1) dataset_test = (RooDataSet*)(dataset_test1->reduce(Name("dataset_test"),Cut(cut_string)));
////create subset of dataset e.g. for 1d signal template, here leading photon 
  RooDataSet *dset_data_allm_axis1 = (RooDataSet*)(dset_data->reduce(Name("dset_data_allm_axis1"),SelectVars(RooArgList(*roovar1,*rooisopv1,*roopt1,*rooeta1,*roosieie1,*roodiphopt,*roodiphomass,*roorho,*roosigma))));
  RooDataSet *  dset_rcone_allm_axis1 = (RooDataSet*)(dset_rcone->reduce(Name("dset_rcone_allm_axis1"),SelectVars(RooArgList(*roovar1,*rooisopv1,*roopt1,*roosieie1,*rooeta1,*roodiphopt,*roodiphomass,*roorho,*roosigma))));
  RooDataSet *  dset_sideb_allm_axis1 = (RooDataSet*)(dset_sideb->reduce(Name("dset_sideb_allm_axis1"),SelectVars(RooArgList(*roovar1,*rooisopv1,*roopt1,*roosieie1,*rooeta1,*roorho,*roodiphopt,*roodiphomass,*roosigma))));
  assert(dset_data_allm_axis1);assert(dset_rcone_allm_axis1);assert(dset_sideb_allm_axis1); 

  

 //compare two distributions, e.g. bkg and sideband
if(check) { 
	  dset_test = (RooDataSet*)(dataset_test->reduce(Name("dset_test"),SelectVars(RooArgList(*roovar1,*rooisopv1,*roopt1,*roosieie1,*rooeta1,*roorho,*roodiphopt,*roodiphomass,*roosigma))));

	  RooDataHist *test_side = new RooDataHist("test_side","test_side",RooArgSet(*roovar1));
	  test_side->add(*dset_sideb_allm_axis1);
//	  test_side->add(*dset_rcone_allm_axis1);
	  RooDataHist *test_bkg = new RooDataHist("test_bkg","test_bkg",RooArgSet(*roovar1));
	  test_bkg->add(*dset_test);
	  RooHistPdf *test_sidepdf = new RooHistPdf("test_sidepdf","test_sidepdf",RooArgList(*roovar1),*test_side);
	  RooHistPdf *test_bkgpdf = new RooHistPdf("test_bkgpdf","test_bkgpdf",RooArgList(*roovar1),*test_bkg);
	  
//	  TCanvas *c3=new TCanvas("c3","isolation comparison");
	  TCanvas *c3=new TCanvas("c3","isolation comparison");
	  c3->cd();
	  c3->cd();
	  c3->SetLogy();
	//  RooPlot *testp = roovar1->frame(Title("sideband bkg comparison before reweighting"));
	  RooPlot *testp = roovar1->frame(Title("sig comparison before reweighting"));
//	  test_sidepdf->plotOn(testp,LineStyle(kDashed),LineColor(kRed),Name("sideband"));
//	  test_bkgpdf->plotOn(testp,LineStyle(kDashed),LineColor(kBlack),Name("bkg"));
	  test_sidepdf->plotOn(testp,LineStyle(kDashed),LineColor(kRed),Name("sig new"));
	  test_bkgpdf->plotOn(testp,LineStyle(kDashed),LineColor(kBlack),Name("sig old"));
	  testp->Draw(); 
	  //c3->SaveAs("../plots/sideband_bkg_comparison_noreweight.root") ;
//	  c3->SaveAs("../plots/sideband_bkg_comparison_noreweight.png") ;
	  c3->SaveAs("../plots/side_comparison_noreweight.root") ;
	  c3->SaveAs("../plots/side_comparison_noreweight.png") ;
  }
      Double_t dpmq2[nq];  
      //Double_t masserrorl[nq];
      Double_t masserror[nq];
	  
      
	  // ratioplot(dset_data_axis1, dset_rcone_axis1,kTRUE);
 if(massbinned){
	   
	 for(int i=1; i<=nq ;i++)
	 {
		RooDataSet *dset_rcone_temp=NULL;
		dset_rcone_temp=new RooDataSet(*dset_rcone_allm_axis1,Form("dset_rcone%u_%u_axis1",i,startbin));
		RooDataSet *dset_sideb_temp=NULL;
		dset_sideb_temp=new RooDataSet(*dset_sideb_allm_axis1,Form("dset_sideb%u_%u_axis1",i,startbin));
    	RooDataSet *dset_massc = NULL;
    	RooDataSet *dset_sigmassc = NULL;
    	RooDataSet *dset_bkgmassc = NULL;
		
		cut_s= Form("dpmass>%f_dpmass<%f",dpmq[i-1],dpmq[i]);
	 	cout << cut_s << endl;
  
		dset_massc = (RooDataSet*)(dset_data_allm_axis1->reduce(Name(Form("dset_massc%u_%u",i,startbin)),Cut(cut_s)));
		dset_sigmassc = (RooDataSet*)(dset_rcone_allm_axis1->reduce(Name(Form("dset_sigmassc%u_%u",i,startbin)),Cut(cut_s)));
		dset_bkgmassc = (RooDataSet*)(dset_sideb_allm_axis1->reduce(Name(Form("dset_bkgmassc%u_%u",i,startbin)),Cut(cut_s)));
		do1dfit(dset_massc,dset_sigmassc,dset_bkgmassc,dset_rcone_temp,dset_sideb_temp,dset_test,i, purity, purityerr, pullerror);
		//masserror[i]=0.;
	
		cut_s.Clear();
		//dont want first 0 element
		//dpmq2[i-1]     =dpmq2[i]-(dpmq[i]-dpmq[i-1])/2;
		dpmq2[i-1]     =(dpmq[i]+dpmq[i-1])/2;
	//	masserrorl[i]=dpmq2[i]-dpmq[i-1];
	//	masserrorr[i]=dpmq[i]-dpmq2[i];
     	masserror[i]=(dpmq[i]-dpmq[i-1])/2;

		
	 }
        
//	 plot_purity_massbin(dpmq2,masserrorl,masserrorr,purity,purityerr, true) ; 
//	 plot_purity_massbin(dpmq2,masserrorl,masserrorr,purity,pullerror,false) ; 
	 plot_purity_massbin(dpmq2,masserror,purity,purityerr, true) ; 
	 plot_purity_massbin(dpmq2,masserror,purity,pullerror,false) ; 
 }
 else{
//no mass bining, one fit for all diphoton masses, just but high number a the end of title (10000)
   int i= 10000;
   do1dfit(dset_data_allm_axis1,dset_rcone_allm_axis1,dset_sideb_allm_axis1,dset_rcone_allm_axis1,dset_sideb_allm_axis1,dset_test,i, purity, purityerr,pullerror);
 }
   return 0;
    
}


int do1dfit(RooDataSet *dset_data_axis1,RooDataSet * dset_sigcut,RooDataSet* dset_bkgcut,RooDataSet *dset_rcone_axis1,RooDataSet *dset_sideb_axis1,RooDataSet*dset_testbkg,int im, double pu[], double puerr[],double pullerr[]){
//reweighting
	cout << "massbin " << im << endl;
	dset_data_axis1->Print();
	dset_rcone_axis1->Print();
	dset_sideb_axis1->Print();
//	RooDataSet* dset_sigorg=NULL;
//	RooDataSet* dset_bkgorg=NULL;
	RooDataSet* dset_sigorg=(RooDataSet*)dset_rcone_axis1->Clone("dset_sigorg");
	RooDataSet* dset_bkgorg=(RooDataSet*)dset_sideb_axis1->Clone("dset_bkgorg");
	//RooDataSet* dset_bkgorg=new RooDataSet(*dset_sideb_axis1,"dset_bkgorg");
//	dset_bkgorg=dset_sideb_axis1;
	if(debug) dset_testbkg->Print();
  // 	cout << "ratio Entries truth/data weighted " << dset_rcone_axis1->sumEntries()/ dset_data_axis1->sumEntries()  << endl;
	   reweight_rhosigma(w_rhor,w_rhoro,w_rhorn,w_rhor2o,w_rhornr,w_rhor2,w_sigmarn,w_sigmar2o,w_sigmarnr,w_sigmar2,&dset_rcone_axis1,dset_data_axis1); 
		reweight_rhosigma(w_rho,w_rhoo,w_rhon,w_rho2o,w_rhonr,w_rho2,w_sigman,w_sigma2o,w_sigmanr,w_sigma2,&dset_sideb_axis1,dset_data_axis1); 
		reweight_eta_1d(w_etar,w_etaro,w_etarn,w_etar2o,w_etarnr,w_etar2,&dset_rcone_axis1,dset_data_axis1,1);
		reweight_eta_1d(w_eta,w_etao,w_etan,w_eta2o,w_etanr,w_eta2,&dset_sideb_axis1,dset_data_axis1,1);
		reweight_pt_1d(w_pt,w_pto,w_ptn,w_pt2o,w_ptnr,w_pt2,&dset_sideb_axis1,dset_data_axis1,1); 
//no reweighting for pt if random cone applied
		if(truthfit){
		  reweight_pt_1d(w_ptr,w_ptro,w_ptrn,w_ptr2o,w_ptrnr,w_ptr2,&dset_rcone_axis1,dset_data_axis1,1);
		  reweight_diphotonpt_1d(w_diphoptr,w_diphoptro,w_diphoptrn,w_diphoptr2o,w_diphoptrnr,w_diphoptr2,&dset_rcone_axis1,dset_data_axis1); 
		}
		reweight_diphotonpt_1d(w_diphopt,w_diphopto,w_diphoptn,w_diphopt2o,w_diphoptnr,w_diphopt2,&dset_sideb_axis1,dset_data_axis1); 

     if(check){
	 	  reweight_rhosigma(w_rhor,w_rhoro,w_rhorn,w_rhor2o,w_rhornr,w_rhor2,w_sigmarn,w_sigmar2o,w_sigmarnr,w_sigmar2,&dset_testbkg,dset_data_axis1);  
	 	  reweight_eta_1d(w_etar,w_etaro,w_etarn,w_etar2o,w_etarnr,w_etar2,&dset_testbkg,dset_data_axis1,1);
		  reweight_pt_1d(w_ptr,w_ptro,w_ptrn,w_ptr2o,w_ptrnr,w_ptr2,&dset_testbkg,dset_data_axis1,1);
		  reweight_diphotonpt_1d(w_diphoptr,w_diphoptro,w_diphoptrn,w_diphoptr2o,w_diphoptrnr,w_diphoptr2,&dset_testbkg,dset_data_axis1); 
		  RooDataHist *test_side2 = new RooDataHist("test_side2","test_side2",RooArgSet(*roovar1));
		  test_side2->add(*dset_sideb_axis1);
		  RooDataHist *test_bkg2 = new RooDataHist("test_bkg2","test_bkg2",RooArgSet(*roovar1));
		  test_bkg2->add(*dset_testbkg);
		  RooHistPdf *test_sidepdf2 = new RooHistPdf("test_sidepdf2","test_sidepdf2",RooArgList(*roovar1),*test_side2);
		  RooHistPdf *test_bkgpdf2 = new RooHistPdf("test_bkgpdf2","test_bkgpdf2",RooArgList(*roovar1),*test_bkg2);
		  TCanvas *c4=new TCanvas("c4","isolation comparison after reweighting");
		  c4->cd();
		  c4->SetLogy();
		  RooPlot *testp2 = roovar1->frame(Title("sideband bkg comparison after reweighting"));
		  test_sidepdf2->plotOn(testp2,LineStyle(kDashed),LineColor(kRed),Name("sideband"));
		  test_bkgpdf2->plotOn(testp2,LineStyle(kDashed),LineColor(kBlack),Name("bkg"));
		  testp2->Draw(); 
	//	  c4->SaveAs("../plots/sideband_bkg_comparison_reweight.root") ;
//		  c4->SaveAs("../plots/sideband_bkg_comparison_reweight.png") ;
		  c4->SaveAs("../plots/sig_comparison_reweight.root") ;
		  c4->SaveAs("../plots/sig_comparison_reweight.png") ;
	 }
		 //check distributions, the one over whole massspectrum, in massbin and after reweighting for signal template
		 
	      Double_t sumentries=0; Double_t sumentries_cut=0; Double_t sumentries_rew=0; 
	      sumentries=dset_sigorg->sumEntries(); 
		  sumentries_cut=dset_sigcut->sumEntries(); 
		  sumentries_rew=dset_rcone_axis1->sumEntries(); 
		  cout << "sumentries " << sumentries << " sumentries_cut "  << sumentries_cut << " sumentries_rew " << sumentries_rew << endl;
		  TCanvas *csig=new TCanvas("csig","signal template comparison");
	      TLegend *leg1 = new TLegend(0.6,0.8,0.9,0.9);
		  leg1->SetFillColor(kWhite); 
		  csig->cd();
		  RooPlot *sig = roovar1->frame(Title("signal template comparison"));
		  dset_sigorg->plotOn(sig,Binning(tbins),MarkerStyle(20),MarkerColor(kBlack),LineColor(kBlack),Name("sigall"));
		  dset_sigcut->plotOn(sig,Rescale(sumentries/sumentries_cut),Binning(tbins),MarkerStyle(20),MarkerColor(kBlue),LineColor(kBlue),Name("sigmassbin"));
		  dset_rcone_axis1->plotOn(sig,Rescale(sumentries/sumentries_rew),Binning(tbins),MarkerStyle(20),MarkerColor(kRed),LineColor(kRed),Name("sigrew"));
		   
		  sig->Draw();
	  	  leg1->AddEntry("sigall","whole mass spectrum","l");
	  	  leg1->AddEntry("sigmassbin","in massbin","l");
	  	  leg1->AddEntry("sigrew","after reweight","l");
		  leg1->Draw();
		  TString titlesig;
		  titlesig= Form("sigtemplate_%s_%s_sigtemplatecomp_%s_startb%u_%u",(EBEB)? "EBEB": "notEBEB",(truthfit)? "truth": "rcone_sideb",cut_s.Data(),startbin_q,im);
	  	  csig->SaveAs(Form("../plots/March18/massbinned20_4bins/%s.root",titlesig.Data()));
	  	  csig->SaveAs(Form("../plots/March18/massbinned20_4bins/%s.png",titlesig.Data())); 
		  //compare the background template to check if binning correct
	      Double_t sumentriesb=0; Double_t sumentries_cutb=0; Double_t sumentries_rewb=0; 
	      sumentriesb=dset_bkgorg->sumEntries(); 
		  sumentries_cutb=dset_bkgcut->sumEntries(); 
		  sumentries_rewb=dset_sideb_axis1->sumEntries(); 
		  cout << "sumentries bkg " << sumentriesb << " sumentries_cut bkg "  << sumentries_cutb << " sumentries_rew bkg" << sumentries_rewb << endl;
		  TCanvas *cbkg=new TCanvas("cbkg","bkgnal template comparison");
		  cbkg->cd();
		  cbkg->SetLogy();
		  RooPlot *bkg = roovar1->frame(Title("bkgnal template comparison"));
		  dset_bkgorg->plotOn(bkg,Binning(tbins),MarkerStyle(20),MarkerColor(kBlack),LineColor(kBlack),Name("bkgall"));
		  dset_bkgcut->plotOn(bkg,Rescale(sumentriesb/sumentries_cutb),Binning(tbins),MarkerStyle(20),MarkerColor(kBlue),LineColor(kBlue),Name("bkgmassbin"));
		  dset_sideb_axis1->plotOn(bkg,Rescale(sumentriesb/sumentries_rewb),Binning(tbins),MarkerStyle(20),MarkerColor(kRed),LineColor(kRed),Name("bkgrew"));
		  bkg->Draw();
		  leg1->Draw();
	  	  TString titlebkg;
		  titlebkg= Form("bkg_template_%s_%s_bkgtemplatecomp_%s_startb%u_%u",(EBEB)? "EBEB": "notEBEB",(truthfit)? "truth": "rcone_sideb",cut_s.Data(),startbin_q,im);
	  	  cbkg->SaveAs(Form("../plots/March18/massbinned20_4bins/%s.root",titlebkg.Data()));
	  	  cbkg->SaveAs(Form("../plots/March18/massbinned20_4bins/%s.png",titlebkg.Data()));
         
	   
		  //  if(massbinned) {
	    TCanvas * c_reweights=new TCanvas("c_reweights","c_reweights");
		c_reweights->Divide(3,2);  c_reweights->cd(1); w_ptn->Draw("colz");
		c_reweights->cd(2);  w_pt2o->Draw("colz");  c_reweights->cd(3);   w_ptnr->Draw("colz");
	    c_reweights->cd(4);  w_pt2->Draw("colz"); c_reweights->cd(5);w_pt->Draw();
		c_reweights->cd(6);   w_pto->Draw(); 
		TString titlerew; 
		titlerew= Form("../plots/March18/massbinned20_4bins/reweight_%s_%s_%s_startbin%u_%u",(EBEB)? "EBEB": "notEBEB",(truthfit)? "truth": "rcone_sideb",cut_s.Data(),startbin_q,im);
	  	 c_reweights->SaveAs(Form("%s_ptbkg.root",titlerew.Data()));
	  	 c_reweights->SaveAs(Form("%s_ptbkg.png",titlerew.Data()));
		  
		TCanvas * c_reweights2=new TCanvas("c_reweights2","c_reweights2");
		c_reweights2->Divide(3,2); c_reweights2->cd(1);  w_diphoptn->Draw("colz");  c_reweights2->cd(2);    w_diphopt2o->Draw("colz");
		c_reweights2->cd(3);    w_diphoptnr->Draw("colz");c_reweights2->cd(4);  w_diphopt2->Draw("colz");
		c_reweights2->cd(5);w_diphopt->Draw();  c_reweights2->cd(6);    w_diphopto->Draw();
	  	c_reweights2->SaveAs(Form("%s_diphoptbkg.root",titlerew.Data()));
	  	c_reweights2->SaveAs(Form("%s_dphoptbkg.png",titlerew.Data()));

		 TCanvas * c_rhoreweights=new TCanvas("c_rhoreweights","c_rhoreweights");
		  c_rhoreweights->Divide(3,2);  c_rhoreweights->cd(1);  w_rhon->Draw("colz");
		  c_rhoreweights->cd(2);    w_rho2o->Draw("colz");  c_rhoreweights->cd(3);    w_rhonr->Draw("colz");
		  c_rhoreweights->cd(4);  w_rho2->Draw("colz");  c_rhoreweights->cd(5);
		  w_rho->Draw();c_rhoreweights->cd(6);    w_rhoo->Draw();
	  	  c_rhoreweights->SaveAs(Form("%s_rhobkg.root",titlerew.Data()));
	  	  c_rhoreweights->SaveAs(Form("%s_rhobkg.png",titlerew.Data()));
		  TCanvas * c_sigmareweights=new TCanvas("c_sigmareweights","c_sigmareweights");
		  c_sigmareweights->Divide(3,2);  c_sigmareweights->cd(1);  w_sigman->Draw("colz");  c_sigmareweights->cd(2); 
		   w_sigma2o->Draw("colz");  c_sigmareweights->cd(3);    w_sigmanr->Draw("colz");  c_sigmareweights->cd(4);
		  w_sigma2->Draw("colz"); 
	  	  c_sigmareweights->SaveAs(Form("%s_sigmabkg.root",titlerew.Data()));
	  	  c_sigmareweights->SaveAs(Form("%s_sigmabkg.png",titlerew.Data()));
		  
		  TCanvas * c_etareweights=new TCanvas("c_etareweights","c_etareweights");
		  c_etareweights->Divide(3,2);  c_etareweights->cd(1);  w_etan->Draw("colz");
		  c_etareweights->cd(2);    w_eta2o->Draw("colz");  c_etareweights->cd(3);    w_etanr->Draw("colz");
		  c_etareweights->cd(4);  w_eta2->Draw("colz");  c_etareweights->cd(5);
		  w_eta->Draw();c_etareweights->cd(6);   w_etao->Draw(); 
	  	  c_etareweights->SaveAs(Form("%s_etabkg.root",titlerew.Data()));
	  	  c_etareweights->SaveAs(Form("%s_etabkg.png",titlerew.Data()));
//signal template
  		  TCanvas * c_reweightsr=new TCanvas("c_reweightsr","c_reweightsr");
		c_reweightsr->Divide(3,2);  c_reweightsr->cd(1); w_ptn->Draw("colz");
		c_reweightsr->cd(2);  w_pt2o->Draw("colz");  c_reweightsr->cd(3);   w_ptnr->Draw("colz");
	    c_reweightsr->cd(4);  w_pt2->Draw("colz"); c_reweightsr->cd(5);w_pt->Draw();
		c_reweightsr->cd(6);   w_pto->Draw();
	  	c_reweightsr->SaveAs(Form("%s_ptsig.root",titlerew.Data()));
	  	c_reweightsr->SaveAs(Form("%s_ptsig.png",titlerew.Data()));
		  
		TCanvas * c_reweightsr2=new TCanvas("c_reweightsr2","c_reweightsr2");
		c_reweightsr2->Divide(3,2); c_reweightsr2->cd(1);  w_diphoptn->Draw("colz");  c_reweightsr2->cd(2);    w_diphopt2o->Draw("colz");
		c_reweightsr2->cd(3);    w_diphoptnr->Draw("colz");c_reweightsr2->cd(4);  w_diphopt2->Draw("colz");
		c_reweightsr2->cd(5);w_diphopt->Draw();  c_reweightsr2->cd(6);    w_diphopto->Draw();
	  	c_reweightsr2->SaveAs(Form("%s_diphoptsig.root",titlerew.Data()));
	  	c_reweightsr2->SaveAs(Form("%s_diphoptsig.png",titlerew.Data()));

		  TCanvas * c_rhoreweightsr=new TCanvas("c_rhoreweightsr","c_rhoreweightsr");
		  c_rhoreweightsr->Divide(3,2);  c_rhoreweightsr->cd(1);  w_rhon->Draw("colz");
		  c_rhoreweightsr->cd(2);    w_rho2o->Draw("colz");  c_rhoreweightsr->cd(3);    w_rhonr->Draw("colz");
		  c_rhoreweightsr->cd(4);  w_rho2->Draw("colz");  c_rhoreweightsr->cd(5);
		  w_rho->Draw();c_rhoreweightsr->cd(6);    w_rhoo->Draw();
	  	c_rhoreweightsr->SaveAs(Form("%s_rhosig.root",titlerew.Data()));
	  	c_rhoreweightsr->SaveAs(Form("%s_rhosig.png",titlerew.Data()));
		  TCanvas * c_sigmareweightsr=new TCanvas("c_sigmareweightsr","c_sigmareweightsr");
		  c_sigmareweightsr->Divide(3,2);  c_sigmareweightsr->cd(1);  w_sigman->Draw("colz");  c_sigmareweightsr->cd(2); 
		   w_sigma2o->Draw("colz");  c_sigmareweightsr->cd(3);    w_sigmanr->Draw("colz");  c_sigmareweightsr->cd(4);
		  w_sigma2->Draw("colz"); 
	  	c_sigmareweightsr->SaveAs(Form("%s_sigmasig.root",titlerew.Data()));
	  	c_sigmareweightsr->SaveAs(Form("%s_sigmasig.png",titlerew.Data()));
		  
		  TCanvas * c_etareweightsr=new TCanvas("c_etareweightsr","c_etareweightsr");
		  c_etareweightsr->Divide(3,2);  c_etareweightsr->cd(1);  w_etan->Draw("colz");
		  c_etareweightsr->cd(2);    w_eta2o->Draw("colz");  c_etareweightsr->cd(3);    w_etanr->Draw("colz");
		  c_etareweightsr->cd(4);  w_eta2->Draw("colz");  c_etareweightsr->cd(5);
		  w_eta->Draw();c_etareweightsr->cd(6);   w_etao->Draw(); 
	  	  c_etareweightsr->SaveAs(Form("%s_etasig.root",titlerew.Data()));
	    	c_etareweightsr->SaveAs(Form("%s_etasig.png",titlerew.Data()));
  // }
  if(check) {
	  cout << "after reweighting" << endl; 
	  cout << "sig " << dset_rcone_axis1->sumEntries() << endl;
	  cout << "bkg " << dset_sideb_axis1->sumEntries() << endl;
	  cout << "data " << dset_data_axis1->sumEntries() << endl;
	  cout << "signal+ fakes " << dset_rcone_axis1->sumEntries()+ dset_sideb_axis1->sumEntries()<< endl;
	  cout << "one? " << (dset_rcone_axis1->sumEntries()+ dset_sideb_axis1->sumEntries())/dset_data_axis1->sumEntries()<< endl;
	  cout << "ratio Entries truth/data weighted " << dset_rcone_axis1->sumEntries()/dset_data_axis1->sumEntries()  << endl;
	  }
	  if(check){
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
	   }
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
	  RooDataHist *drconehist_axis1=NULL;
	  RooDataHist *dsidebhist_axis1=NULL;
	  RooHistPdf *drconepdf_axis1 =NULL;
	  RooHistPdf *dsidebpdf_axis1 =NULL;
	  drconehist_axis1 = new RooDataHist("drconehist_axis1","drconehist_axis1",var);
	  drconehist_axis1->add(*dset_rcone_axis1);
	  if(debug){
	  	drconehist_axis1->Print("V");
	  }
	  dsidebhist_axis1 = new RooDataHist("dsidebhist_axis1","dsidebhist_axis1",var);
	  dsidebhist_axis1->add(*dset_sideb_axis1);
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
	  model_axis1 = new RooAddPdf("model_axis1","model_axis1",RooArgList(*drconepdf_axis1,*dsidebpdf_axis1),RooArgList(*fsig1),kFALSE);
	  //if does not converge take NumCPU strategy 2
	  RooFitResult *firstpass1;
	  TFile resFile(Form("../plots/March18/massbinned20_4bins/fitresult_%s_%s_sumw2erroff_%u_%s_range_%u_%u.root",(truthfit)? "truth": "rcone_sideb",(EBEB) ? "EBEB": "notEBEB",im,cut_s.Data(),startbin_q,endbin_q),"RECREATE"); 
	  //second fit with sumw2error on for PULL fct= (ptrue-pmeas)/sigma_meas
	  //if does not converge take NumCPU strategy 2
	 //SumW2Error takes statistics of MC (weights) into account
	 // By default the fit is executed through the MINUIT commands MIGRAD, HESSE and MINOS in succession 
	 //firstpass1 = model_axis1->fitTo(*dset_data_axis1, NumCPU(8), Extended(false),SumW2Error(kTRUE),Save(kTRUE));
	 cout << "firstpass1 sumw2error false" << endl;
	  resFile.cd();
	  firstpass1 = model_axis1->fitTo(*dset_data_axis1, NumCPU(8), Extended(false),SumW2Error(kFALSE),Save(kTRUE));
	  firstpass1->Print();
	  firstpass1->Write();
	  resFile.Write();  
	  resFile.Close();
	  //second fit with sumw2error on for PULL fct= (ptrue-pmeas)/sigma_meas
	  RooAddPdf *model2_axis1=NULL;
	  model2_axis1 = new RooAddPdf("model2_axis1","model2_axis1",RooArgList(*drconepdf_axis1,*dsidebpdf_axis1),RooArgList(*fsig1),kFALSE);
	  RooFitResult *firstpass2;
	  TFile resFile2(Form("../plots/March18/massbinned20_4bins/fitresult_%s_%s_sumw2erron_%u_%s_range_%u_%u.root",(truthfit)? "truth": "rcone_sideb",(EBEB) ? "EBEB": "notEBEB",im,cut_s.Data(),startbin_q,endbin_q),"RECREATE"); 
      cout << "firstpass2 sumw2error true " << endl;
  	  resFile.cd();
	  firstpass2 = model2_axis1->fitTo(*dset_data_axis1, NumCPU(8), Extended(false),SumW2Error(kTRUE),Save(kTRUE));
	  firstpass2->Print();
	  firstpass2->Write();
	  resFile2.Write();  
	  resFile2.Close();
	  if(massbinned)
	  {
	  	pu[im-1]=fsig1->getVal();
		puerr[im-1]=fsig1->getPropagatedError(*firstpass1);
	  	pullerr[im-1]=fsig1->getPropagatedError(*firstpass2);
	  }
	  TLatex b;b.SetNDC();b.SetTextSize(0.06);b.SetTextColor(kRed);
	  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
	  TLegend *leg = new TLegend(0.15,0.8,0.35,0.9);
	  c1->cd(1);
	  TString title; 
	  if(dohggv){
	  //	title= Form("1d fit for hgg vertex no SumW2Error, %s, reweighted,%u rangebins: %u -  %u",(EBEB)? "EBEB": "notEBEB" ,im, startbin_q, endbin_q);
	  	title= Form("1d fit for hgg vertex no SumW2Error, %s, reweighted,%u mass: %s",(EBEB)? "EBEB": "notEBEB" ,im,cut_s.Data());
      }
	  else if(dovzero){
	  	title= Form("1d fit for vertex[0], %s, reweighted,%u",(EBEB)? "EBEB": "notEBEB" ,im);
	  }
	  else if(dowv){
	  	title= Form("1d fit for worst isolation, %s, reweighted,%u",(EBEB)? "EBEB": "notEBEB" ,im);
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
	  model_axis1->Print();
	  dset_data_axis1->Print();
	  //save fit plots
	  TString titleout;
//	  if(dohggv) titleout= Form("%s_reweight_%s_hggvtx_noSUMW2_varbinning%u_range%u_%u",(EBEB)? "EBEB": "notEBEB",(truthfit)? "truth": "rcone_sideb",im,startbin_q, endbin_q);
	 // if(dohggv) titleout= Form("%s_reweight_%s_hggvtx_noSUMW2_varbinning_startb_%u_%u_sigcomp",(EBEB)? "EBEB": "notEBEB",(truthfit)? "truth": "rcone_sideb",startbin_q,im);
	  if(dohggv) titleout= Form("%s_reweight_%s_hggvtx_noSUMW2_%s_startb_%u_%u_sigcomp",(EBEB)? "EBEB": "notEBEB",(truthfit)? "truth": "rcone_sideb",cut_s.Data(),startbin_q,im);
	  else if(dovzero) titleout= "nooverflow_reweight_fulletarange_rcone_sideband_vzero";
	  else if(dowv) titleout= "nooverflow_reweight_fulletarange_rcone_sideband_wv";
	  const char* outfile=Form("../plots/March18/massbinned20_4bins/roofit_%s_%s", (isdata)? "data": "mc",titleout.Data());
	  c1->SaveAs(Form("%s_lin.png",outfile));c1->SaveAs(Form("%s_lin.root",outfile));c1->SaveAs(Form("%s_lin.pdf",outfile));
	  c2->SaveAs(Form("%s_log.png",outfile));c2->SaveAs(Form("%s_log.root",outfile));c2->SaveAs(Form("%s_log.pdf",outfile));
	delete dset_data_axis1;
	delete dset_rcone_axis1;
	delete dset_sideb_axis1;
	delete dset_sigorg;
	delete dset_bkgorg;
	delete dset_testbkg;
	delete drconehist_axis1;
	delete dsidebhist_axis1;
    
	delete drconepdf_axis1;  
	delete dsidebpdf_axis1;
	delete c1;
	delete c2;
	delete csig;
	delete cbkg;
	delete c_reweights;
	delete c_reweights2;
	delete c_rhoreweights;
	delete c_sigmareweights;
	delete c_etareweights;
	delete c_reweightsr;
	delete c_reweightsr2;
	delete c_rhoreweightsr;
	delete c_sigmareweightsr;
	delete c_etareweightsr;

	if(!massbinned){
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
		
	}

	print_mem();
	return 0;

}


