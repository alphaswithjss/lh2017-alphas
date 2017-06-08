//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldExample.cxx 348 2014-08-08 22:18:23Z T.J.Adye@rl.ac.uk $
//
// Description:
//      Simple example usage of the RooUnfold package using toy MC.
//
// Authors: Tim Adye <T.J.Adye@rl.ac.uk> and Fergus Wilson <fwilson@slac.stanford.edu>
//
//==============================================================================

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
//#include "RooUnfoldSvd.h"
//#include "RooUnfoldTUnfold.h"

#include "SoftDropHelper.h"
#include "TRandom3.h"

#endif

//==============================================================================
// Example Unfolding
//==============================================================================

double chi2(TH1D* data, TH1D* mytemplate){
  double out = 0.;
  //mytemplate->Scale(data->Integral()/mytemplate->Integral());
  for (int i=1; i<=data->GetNbinsX(); i++){
    out+=pow(data->GetBinContent(i)-mytemplate->GetBinContent(mytemplate->FindBin(data->GetBinCenter(i))),2);
  }
  return out;
}

//---------------------------------------------------------------------------------------------------
//define constants:

//Z mass
double mz=91.18;

//pi
double pi=3.141592654;

//Color charges
double CF=4.0/3.0;
double CA=3.0;

//definitions associated with alpha_s
double alpha_s_mZ_standard=0.118;
double beta0=23.0/(12.0*pi);
double Bg=-11.0/12.0+5.0/18.0;
double Bq=-3.0/4.0;

//running of alpha_s
double alpha_s(const double kappa, const double asmz){

  return asmz/(1.0+2.0*asmz*beta0*log(kappa/mz));

}


//---------------------------------------------------------------------------------------------------
//define several functions that I need

//define shorthand for xlog(x)
double xlog(const double x){

        return x*log(fabs(x));

}

//define step function
double step(const double x){

        double var = 0.0;
        
        if(x > 0){
          var = 1.0;}         

        return var;

}


//---------------------------------------------------------------------------------------------------
//define different components of the mass distribution

//define radiator
double Rad(const double e2, const double Ci, const double Bi, const double beta0, const double as, const double zcut){
        
        double log_arg1=1-2*as*beta0*log(1/e2); 
        double log_arg2=1-as*beta0*log(1/e2); 
        double log_arg3=1-as*beta0*log(1/e2)-as*beta0*log(1/zcut);
        double log_arg4=1-2*as*beta0*log(1/zcut);

        double part1=(Ci/(2*pi*as*beta0*beta0))*(xlog(log_arg1)-2*xlog(log_arg2)-2*as*beta0*Bi*log(fabs(log_arg2)));
        
        double part2=(Ci/(2*pi*as*beta0*beta0))*(-xlog(log_arg4)-2*xlog(log_arg2)-2*as*beta0*Bi*log(fabs(log_arg2))+2*xlog(log_arg3));

        return part1*step(e2-zcut)+part2*step(zcut-e2);
}


//define Sudakov
double sudakov(const double e2, const double Ci, const double Bi, const double zcut, const double as_mz, const double pt){

        return exp(-Rad(e2,Ci,Bi,beta0,alpha_s(pt*0.8,as_mz),zcut)) ;

}

//Derivative of Sudakov (here this is just a simple numerical derivative with step delta)
double diffsudakov(const double e2, const double Ci, const double Bi, const double zcut, const double as_mz, const double pt, const double delta){

  return (1/delta)*(sudakov(e2+delta,Ci,Bi,zcut,as_mz,pt)-sudakov(e2,Ci,Bi,zcut,as_mz,pt));

}


//cutoff
double cutoff(const double e2, const double zcut, const double pt){
        
        return step(e2-1/(zcut*pt*pt));

}


//-------------------------------------------------------------------------------------------------------------
//result to plot (this is e2*dsigma/de2)
//the arguments are:
//e2 (self explanatory)
//Ci=color charge, set to CF for quark, CA for gluon
//Bi= set to Bq for quark, Bg for gluon
//zcut (self explanatory)
//as_mz= alpha_s value
//pt (self explanatory)
//delta= step taken in the numerical derivative
double soft_drop_mass(const double e2, const double Ci, const double Bi, const double zcut, const double as_mz, const double pt, const double delta){

  return e2*diffsudakov(e2,Ci,Bi,zcut,as_mz,pt,delta)*cutoff(e2,zcut,pt)/(diffsudakov(zcut+0.0001,Ci,Bi,zcut,as_mz,pt,delta));

}

void RooUnfoldExample(){

  TRandom3* myrand = new TRandom3(1);

   vector<double> bins;
   bins.push_back(0.05);
   while (bins[bins.size()-1] > 1e-4){
    bins.push_back(bins[bins.size()-1]/1.1);
   }
   std::reverse(bins.begin(), bins.end());
   TH1F* gluons_pdf = new TH1F("","",bins.size()-1,&bins[0]);
   TH1F* quarks_pdf = new TH1F("","",bins.size()-1,&bins[0]);

   vector<double> bins2;
   bins2.push_back(0.05);
   while (bins2[bins2.size()-1] > 1e-4){
    bins2.push_back(bins2[bins2.size()-1]/2);
   }
   std::reverse(bins2.begin(), bins2.end());

   for (int i=1; i<=gluons_pdf->GetNbinsX(); i++){
      double e2 = gluons_pdf->GetXaxis()->GetBinCenter(i);
      gluons_pdf->SetBinContent(i, soft_drop_mass(e2,CA,Bg,0.1,0.10,500.0,0.000001)/e2);
      quarks_pdf->SetBinContent(i, soft_drop_mass(e2,CF,Bq,0.1,0.10,500.0,0.000001)/e2);      
   }

  //double q_frac = 0.0;
  double uJMS = 0.05;
  double uJMR = 0.2;

  double alpha_s_min = 0.01;
  double alpha_s_max = 0.25;
  double dalpha_s = (alpha_s_max-alpha_s_min)/1000;

  TH1D* myscan = new TH1D("","",30,0,1.1);
  for (int q=1; q<myscan->GetNbinsX(); q++){
    double q_frac = myscan->GetXaxis()->GetBinCenter(q);
    double uncert_JMS = 0.;
    double uncert_JMR = 0.;

     vector<double> reco_vals;
     vector<double> true_vals;
     TH1D* Truth = new TH1D("Truth","Truth",bins2.size()-1,&bins2[0]);
     TH1D* Reco = new TH1D("Reco","Reco",bins2.size()-1,&bins2[0]);
     TH2D* Response = new TH2D("Response","Response",bins2.size()-1,&bins2[0],bins2.size()-1,&bins2[0]);     
     int ntoys = 1000000; //this is MC stats
     for (int ntoy = 0; ntoy < ntoys; ntoy++){
        double coin = myrand->Uniform(0,1);
        double e2_toy = 0.;
        if (coin > q_frac) e2_toy = gluons_pdf->GetRandom();
        else e2_toy = e2_toy = quarks_pdf->GetRandom();
        Truth->Fill(e2_toy);  
        double e2_toy_reco = e2_toy*myrand->Gaus(1.0,0.2);
        Reco->Fill(e2_toy_reco);  
        Response->Fill(e2_toy_reco,e2_toy);  
        reco_vals.push_back(e2_toy_reco);
        true_vals.push_back(e2_toy);
     }

    RooUnfoldResponse myresponse (Reco,Truth,Response);  
    RooUnfoldBayes   unfold_hold2 (&myresponse, Reco, 4);
    TH1D* myout2 = (TH1D*) unfold_hold2.Hreco();
    double min_val = 9999.;
    double Min_as = -1;
    myout2->Scale(1./myout2->Integral("width"));
    for (double as = alpha_s_min; as<alpha_s_max; as+=dalpha_s){
        TH1D* mytemplate = new TH1D(TString::Format("0.2%f_%i",as,q),TString::Format("0.2%f_%i",as,q),bins.size()-1,&bins[0]);
        for (int i=1; i<=mytemplate->GetNbinsX(); i++){
          double e2 = mytemplate->GetXaxis()->GetBinCenter(i);
          mytemplate->SetBinContent(i, (1.-q_frac)*soft_drop_mass(e2,CA,Bg,as,0.10,500.0,0.000001)/e2+q_frac*soft_drop_mass(e2,CF,Bq,as,0.10,500.0,0.000001)/e2); 
        }
        mytemplate->Scale(1./mytemplate->Integral("width"));
        double mychi2 = chi2(myout2,mytemplate);
        if (mychi2 < min_val){
          min_val = mychi2;
          Min_as = as;
        }
     }
    if (1==1){
      double JMS = 1.+uJMS;
      TH1D* Reco2 = new TH1D("Reco","Reco",bins2.size()-1,&bins2[0]);
      for (int ntoy = 0; ntoy < ntoys; ntoy++){
          double e2_toy_reco = reco_vals[ntoy]*JMS;
          Reco2->Fill(e2_toy_reco);  
      }
      RooUnfoldBayes   unfold_hold (&myresponse, Reco2, 4);
      TH1D* myout = (TH1D*) unfold_hold.Hreco();
      myout->Scale(1./myout->Integral("width"));
      myout->Draw("");
      double min_val = 9999.;
      for (double as = alpha_s_min; as<alpha_s_max; as+=dalpha_s){
        TH1D* mytemplate = new TH1D(TString::Format("s0.2%f_%i",as,q),TString::Format("s0.2%f_%i",as,q),bins.size()-1,&bins[0]);
        for (int i=1; i<=mytemplate->GetNbinsX(); i++){
          double e2 = mytemplate->GetXaxis()->GetBinCenter(i);
          mytemplate->SetBinContent(i, (1.-q_frac)*soft_drop_mass(e2,CA,Bg,as,0.10,500.0,0.000001)/e2+q_frac*soft_drop_mass(e2,CF,Bq,as,0.10,500.0,0.000001)/e2); 
        }
        mytemplate->Scale(1./mytemplate->Integral("width"));
        double mychi2 = chi2(myout,mytemplate);
        if (mychi2 < min_val){
          min_val = mychi2;
          uncert_JMS = fabs(as-Min_as)/Min_as;
        }
     }
    }
    if (1==1){
      double JMS = 1.+uJMR;
      TH1D* Reco2 = new TH1D("Reco","Reco",bins2.size()-1,&bins2[0]);
      for (int ntoy = 0; ntoy < ntoys; ntoy++){
          double e2_toy_reco = true_vals[ntoy]+(reco_vals[ntoy]-true_vals[ntoy])*JMS;
          Reco2->Fill(e2_toy_reco);  
      }
      RooUnfoldBayes   unfold_hold (&myresponse, Reco2, 4);
      TH1D* myout = (TH1D*) unfold_hold.Hreco();
      myout->Scale(1./myout->Integral("width"));
      myout->Draw("");
      double min_val = 9999.;
      for (double as = alpha_s_min; as<alpha_s_max; as+=dalpha_s){
        TH1D* mytemplate = new TH1D(TString::Format("r0.2%f_%i",as,q),TString::Format("r0.2%f_%i",as,q),bins.size()-1,&bins[0]);
        for (int i=1; i<=mytemplate->GetNbinsX(); i++){
          double e2 = mytemplate->GetXaxis()->GetBinCenter(i);
          mytemplate->SetBinContent(i, (1.-q_frac)*soft_drop_mass(e2,CA,Bg,as,0.10,500.0,0.000001)/e2+q_frac*soft_drop_mass(e2,CF,Bq,as,0.10,500.0,0.000001)/e2); 
        }
        mytemplate->Scale(1./mytemplate->Integral("width"));
        double mychi2 = chi2(myout,mytemplate);
        if (mychi2 < min_val){
          min_val = mychi2;
          uncert_JMR = fabs(as-Min_as)/Min_as;
        }
     }
    }

    double total_uncert = sqrt(uncert_JMS*uncert_JMS+uncert_JMR*uncert_JMR);
    myscan->SetBinContent(q,total_uncert*100);
  }

  TCanvas* c1 = new TCanvas("","",500,500);
  myscan->GetYaxis()->SetRangeUser(0,30);
  gStyle->SetOptStat(0);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.15);
  myscan->GetXaxis()->SetRangeUser(0,1);
  myscan->SetTitle("");
  myscan->GetXaxis()->SetTitle("quark fraction");
  myscan->GetYaxis()->SetTitle("Total Systematic Uncertainty [%]");
  myscan->SetLabelFont(42);
  myscan->SetTitleFont(42);
  myscan->GetYaxis()->SetLabelFont(42);
  myscan->GetYaxis()->SetTitleFont(42);  
  myscan->GetYaxis()->SetNdivisions(505);
  myscan->GetXaxis()->SetTitleOffset(1.4);
  myscan->GetYaxis()->SetTitleOffset(1.6);   
  myscan->Draw();
  TLatex l  = TLatex(); 
  l.SetNDC();
  l.SetTextColor(1);
  l.DrawLatex(0.2,0.85,"#bf{#scale[0.5]{Softdrop @ MLL, A. Larkoski et al. 1402.2657}}"); 
  l.DrawLatex(0.2,0.8,"#bf{#scale[0.5]{p_{T} = 500 GeV, z_{cut} = 0.1}}");   
  l.DrawLatex(0.2,0.75,"#bf{#scale[0.5]{JMS, JMR uncert = "+TString::Format("%0.2f",uJMS)+", "+TString::Format("%0.2f",uJMR)+"}}");     
  c1->Print("plots/total_qgscan.pdf");  

}

#ifndef __CINT__
int main () { RooUnfoldExample(); return 0; }  // Main program when run stand-alone
#endif
