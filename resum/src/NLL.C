#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <map>
#include <iostream>

//NLL
#include "resum-SD.hh"
#include "resum-mMDT.hh"

//Plotting stuff
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TString.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TMarker.h"

using namespace std;
double alphamax = 0.15;

TH1F* mycopy(TH1F* input){
  TH1F* output = new TH1F("","",100,0.05,alphamax);
  for (int i=1; i<=input->GetNbinsX(); i++){
    output->SetBinContent(i,input->GetBinContent(i));
  }
  return output;
}

TH1F* mycopy2(TH1F* input){

  vector<double> bins2;
  bins2.push_back(input->GetXaxis()->GetBinCenter(1)-0.5*input->GetXaxis()->GetBinWidth(1));
  for (int i=1; i<=input->GetNbinsX(); i++){
    bins2.push_back(input->GetXaxis()->GetBinCenter(i)+0.5*input->GetXaxis()->GetBinWidth(i));
  }
  //for (double xx = 3.; xx < 8; xx+=0.25){
  //  bins2.push_back(xx);
  //}

  TH1F* output = new TH1F("","",bins2.size()-1,&bins2[0]);
  for (int i=1; i<=input->GetNbinsX(); i++){
    output->SetBinContent(i,input->GetBinContent(i));
  }
  return output;
}

double chi2(TH1F* data, TH1F* mytemplate){
  double out = 0.;
  mytemplate->Scale(data->Integral()/mytemplate->Integral());
  for (int i=1; i<=data->GetNbinsX(); i++){
    if (data->GetBinContent(i) > 0) out+=pow((data->GetBinContent(i)-mytemplate->GetBinContent(i))/data->GetBinError(i),2);
  }
  //std::cout << "here " << TMath::Prob(out,data->GetNbinsX()-1) << std::endl;
  return TMath::Prob(out,data->GetNbinsX()-1); //out/(data->GetNbinsX()-1);
}

//define constants

//Z mass
double mz=91.18;

//pi
double pi=3.141592654;

//Color charges
double myCF=4.0/3.0;
double myCA=3.0;

//definitions associated with alpha_s
double alpha_s_mZ=0.118;
double beta0=23.0/(12.0*pi);
double myBg=-11.0/12.0+5.0/18.0;
double myBq=-3.0/4.0;

double alpha_s(const double kappa, const double asmz){

	return asmz/(1.0+2.0*asmz*beta0*log(kappa/mz));

}


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

//define radiator
double Rad(const double e2, const double Ci, const double Bi, const double beta0, const double as, const double zcut, const double sdbeta){
        
        double log_arg1=1-2*as*beta0*log(1/e2); 
        double log_arg2=1-as*beta0*log(1/e2); 
        double log_arg3=1-2*(1+sdbeta)/(2+sdbeta)*as*beta0*log(1/e2)-2/(2+sdbeta)*as*beta0*log(1/zcut);
        double log_arg4=1-2*as*beta0*log(1/zcut);

        double part1=(Ci/(2*pi*as*beta0*beta0))*(xlog(log_arg1)-2*xlog(log_arg2)-2*as*beta0*Bi*log(fabs(log_arg2)));
    
        //only part 2 is modified by the sdbeta value
        double part2=(Ci/(2*pi*as*beta0*beta0))*(-1/(1+sdbeta)*xlog(log_arg4)-2*xlog(log_arg2)-2*as*beta0*Bi*log(fabs(log_arg2))+(2+sdbeta)/(1+sdbeta)*xlog(log_arg3));

        return part1*step(e2-zcut)+part2*step(zcut-e2);
        //return (Ci/(2*pi*as*beta0*beta0))*(-xlog(log_arg4)-2*xlog(log_arg2)-2*as*beta0*Bi*log(fabs(log_arg2)));
        //return Bi*log(fabs(log_arg2));
        //return Bi;
}


//define Sudakov
double sudakov(const double e2, const double Ci, const double Bi, const double zcut, const double sdbeta, const double as_mz, const double pt){

        return exp(-Rad(e2,Ci,Bi,beta0,alpha_s(pt*0.8,as_mz),zcut,sdbeta)) ;

}

//Derivative of Sudakov (here this is just a simple numerimyCAl derivative with step delta)
double diffsudakov(const double e2, const double Ci, const double Bi, const double zcut, const double sdbeta, const double as_mz, const double pt, const double delta){

	return (1/delta)*(sudakov(e2+delta,Ci,Bi,zcut,sdbeta,as_mz,pt)-sudakov(e2,Ci,Bi,zcut,sdbeta,as_mz,pt));

}


//cutoff
double cutoff(const double e2, const double zcut, const double pt){
        
        return step(e2-1/(zcut*pt*pt));

}



//result to plot
double soft_drop_mass(const double e2, const double Ci, const double Bi, const double zcut, const double sdbeta, const double as_mz, const double pt, const double delta){

	return e2*diffsudakov(e2,Ci,Bi,zcut,sdbeta,as_mz,pt,delta)*cutoff(e2,zcut,pt)/(diffsudakov(zcut+0.0001,Ci,Bi,zcut,sdbeta,as_mz,pt,delta));

}

//-------------------------------------------------------------------------------------------------------------
//test evaluating functions at several points

int main() {

  double gfrac1 = 0.8;
  double gfrac2 = 0.2;

  int myseed = 0; //2345

  const int nzcuts = 1; //3
  const int nbetas = 1; //2
  const int nalphas = 3;
  double zcuts[nzcuts] = {0.1};//{0.05,0.1,0.2};
  double betas[nbetas] = {1.}; //0.,1.
  double alphas[nalphas] = {0.5,1.,2.};

  std::vector<double> weights_q_resum;
  std::vector<double> weights_g_resum;
  std::vector<double> weights_q_lo;
  std::vector<double> weights_g_lo;
  std::vector<double> weights_q_nlo; 
  std::vector<double> weights_g_nlo;

  double pt = 500;
  double R = 0.8;
  double alphasMZ = 0.1;
  double muR = 1.;
  double muC = 1.;
  double mufr = 1.; //GeV
  double endpoint = 0.279303; //1/4 for g -> qq; this is the maximal mass / pTR.  Small correction from that for not small angle approx.
  double LambdaNP = 3.; // GeV

  std::cout << " Let's calculate the bounds " << std::endl;

  vector<double> mylogmaxs;
  vector<double> mylogmins;
  for (int iz = 0; iz<nzcuts; iz++){
    for (int ib = 0; ib<nbetas; ib++){
      for (int ia = 0; ia<nalphas; ia++){   
        double mylogmax = -3. + log(zcuts[iz]/0.1) + (alphas[ia]-2.)*log(R);
        double mylogmin = ((alphas[ia]+betas[ib]) / (1. + betas[ib]))*log(LambdaNP / (pt*R)) + ((1. - alphas[ia])/(1. + betas[ib])) * log(zcuts[iz]) + alphas[ia]*log(R);
        mylogmaxs.push_back(mylogmax);
        mylogmins.push_back(mylogmin);
        std::cout << zcuts[iz] << " " << betas[ib] << " " << alphas[ia] << " " << mylogmin << " " << mylogmax << std::endl;
      }
    }
  }

  vector<vector<double> > bins;
  for (unsigned int i=0; i<mylogmaxs.size(); i++){
    vector<double> bins_hold;
    for (double xx = -mylogmaxs[i]; xx < -mylogmins[i]; xx+=fabs(mylogmaxs[i]-mylogmins[i])/20.){
      bins_hold.push_back(xx);
    }    
    bins.push_back(bins_hold);
  }
  //std::reverse(bins.begin(), bins.end());  
  //const std::vector<double> logrhos = bins;

  TLatex l  = TLatex(); 
  l.SetNDC();
  l.SetTextColor(1);

  TCanvas *c1 = new TCanvas("","",500,500);
  //gPad->SetLogx();
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.15);
  gPad->SetTopMargin(0.1);
  gStyle->SetOptStat(0);

  vector<TH2F*> all_histosq;
  vector<TH2F*> all_histosg;

  int mycounter = -1;
   for (int iz = 0; iz<nzcuts; iz++){
    for (int ib = 0; ib<nbetas; ib++){
      for (int ia = 0; ia<nalphas; ia++){

        mycounter++;
        const std::vector<double> logrhos = bins[mycounter];

      double zcut = zcuts[iz];
      double beta = betas[ib];
      double alpha = alphas[ia];

      weights_SD(pt, R, zcut, beta, alpha, alphasMZ, muR, muC, mufr,
              logrhos, endpoint, 
              weights_q_resum,
              weights_g_resum,
              weights_q_lo,
              weights_g_lo,
              weights_q_nlo, 
              weights_g_nlo);  

      if (beta==0){
          weights_mMDT(pt, R, zcut, alpha, alphasMZ, muR, muC, mufr,
                  logrhos, endpoint, 
                  weights_q_resum,
                  weights_g_resum,
                  weights_q_lo,
                  weights_g_lo,
                  weights_q_nlo, 
                  weights_g_nlo);  
      }    
   
       TH1F* ungluons = new TH1F("","",bins[mycounter].size()-1,&bins[mycounter][0]);
       TH1F* unquarks = new TH1F("","",bins[mycounter].size()-1,&bins[mycounter][0]);
       for (int i=1; i<=ungluons->GetNbinsX(); i++){
          ungluons->SetBinContent(i, weights_g_resum[i-1]);
          unquarks->SetBinContent(i, weights_q_resum[i-1]);    
       }

       ungluons->Scale(1./ungluons->Integral());
       unquarks->Scale(1./unquarks->Integral());

       TH1F* gluons = mycopy2(ungluons);
       TH1F* quarks = mycopy2(unquarks);

       gluons->GetYaxis()->SetNdivisions(505);
       gluons->GetXaxis()->SetTitleOffset(1.4);
       gluons->GetYaxis()->SetTitleOffset(1.6);
       gluons->GetYaxis()->SetRangeUser(0,gluons->GetMaximum()*1.5);
       gluons->GetXaxis()->SetTitle("-log((m/Rp_{T})^{2})");
       gluons->GetYaxis()->SetTitle("d#sigma / dlog(#rho)");
       gluons->SetLineColor(2);
       gluons->SetLineWidth(3);
       quarks->SetLineWidth(3);       
       gluons->Draw();
       quarks->SetLineColor(4);
       quarks->SetLineStyle(3);
       quarks->Draw("same");

       //l.DrawLatex(0.16,0.95,"#bf{#scale[0.95]{Softdrop @ NLL, 1704.02210 }}"); 
       l.DrawLatex(0.16,0.93,"#bf{#scale[0.7]{p_{T} = 500 GeV, #alpha = "+TString::Format("%0.2f",alpha)+", z_{cut} = "+TString::Format("%0.2f",zcut)+" and #beta = "+TString::Format("%0.1f",beta)+"}}"); 

       //l.DrawLatex(0.2,0.85,"#bf{#scale[0.5]{Softdrop @ NLL, 1704.02210 }}"); 
       //l.DrawLatex(0.2,0.8,"#bf{#scale[0.5]{p_{T} = 500 GeV, z_{cut} = "+TString::Format("%0.2f",zcut)+" , #beta =" +TString::Format("%0.2f",beta)+"}}");    
       TLegend* leg = new TLegend(.7,.7,0.85,.9);
       leg->SetTextFont(42);
       leg->SetHeader("");
       leg->SetNColumns(1);
       leg->AddEntry(gluons,"#scale[1.5]{gluons}","l");
       leg->AddEntry(quarks,"#scale[1.5]{quarks}","l");
       leg->SetFillStyle(0);
       leg->SetFillColor(0);
       leg->SetBorderSize(0);
       leg->Draw();    
       c1->Print("PDFs_alpha_"+TString::Format("%0.2f",alpha)+"zcut"+TString::Format("%0.2f",zcut)+"_beta_"+TString::Format("%0.2f",beta)+TString::Format("%i",myseed)+".pdf");
    }
  }
}

  mycounter = -1;
 for (int iz = 0; iz<nzcuts; iz++){
  for (int ib = 0; ib<nbetas; ib++){
    for (int ia = 0; ia<nalphas; ia++){

        mycounter++;
        const std::vector<double> logrhos = bins[mycounter];

      double zcut = zcuts[iz];
      double beta = betas[ib];
      double alpha = alphas[ia];

      weights_SD(pt, R, zcut, beta, alpha, alphasMZ, muR, muC, mufr,
              logrhos, endpoint, 
              weights_q_resum,
              weights_g_resum,
              weights_q_lo,
              weights_g_lo,
              weights_q_nlo, 
              weights_g_nlo);  

      if (beta==0){
          weights_mMDT(pt, R, zcut, alpha, alphasMZ, muR, muC, mufr,
                  logrhos, endpoint, 
                  weights_q_resum,
                  weights_g_resum,
                  weights_q_lo,
                  weights_g_lo,
                  weights_q_nlo, 
                  weights_g_nlo);  
      }    
   
       TH1F* gluons = new TH1F("","",bins[mycounter].size()-1,&bins[mycounter][0]);
       TH1F* quarks = new TH1F("","",bins[mycounter].size()-1,&bins[mycounter][0]);
       for (int i=1; i<=gluons->GetNbinsX(); i++){
          gluons->SetBinContent(i, weights_g_resum[i-1]);
          quarks->SetBinContent(i, weights_q_resum[i-1]);    
       }

   //Now, let's start with statistical uncertainty.
   int ntoys = 10000;
   //For each pseudo-experiment, fit alpha_s.
   //First step: generate templates in alpha_s.  Let's start by using only gluons.
   std::map<double,TH1F*> alpha_s_maps;
   std::map<double,TH1F*> alpha_s_maps_quark;
   for (double as = 0.05; as<alphamax; as+=0.001){
      TH1F* gluon_template = new TH1F(TString::Format("0.2%f",as),TString::Format("0.2%f",as),bins[mycounter].size()-1,&bins[mycounter][0]);
      TH1F* quark_template = new TH1F(TString::Format("0.2%f",as),TString::Format("0.2%f",as),bins[mycounter].size()-1,&bins[mycounter][0]);
    
      weights_SD(pt, R, zcut, beta, alpha, as, muR, muC, mufr,
              logrhos, endpoint, 
              weights_q_resum,
              weights_g_resum,
              weights_q_lo,
              weights_g_lo,
              weights_q_nlo, 
              weights_g_nlo);  

      if (beta==0){
          weights_mMDT(pt, R, zcut, alpha, as, muR, muC, mufr,
                  logrhos, endpoint, 
                  weights_q_resum,
                  weights_g_resum,
                  weights_q_lo,
                  weights_g_lo,
                  weights_q_nlo, 
                  weights_g_nlo);  
      }    
   
      for (int i=1; i<=gluon_template->GetNbinsX(); i++){
        gluon_template->SetBinContent(i, weights_g_resum[i-1]); 
        quark_template->SetBinContent(i, weights_q_resum[i-1]);  
        gluon_template->SetBinError(i, 0.); 
        quark_template->SetBinError(i, 0.);  
      }
      alpha_s_maps[as]=gluon_template;
      alpha_s_maps_quark[as]=quark_template;
   }

   //Now, the psuedo-experiments.
   TRandom3 *myrand = new TRandom3(0); //try 0?
   TH1F* gluons_toy = new TH1F("","",bins[mycounter].size()-1,&bins[mycounter][0]);
   TH1F* quarks_toy = new TH1F("","",bins[mycounter].size()-1,&bins[mycounter][0]);
   //ntoys*=100;
   for (int ntoy = 0; ntoy < ntoys; ntoy++){
    double e1_toyg = gluons->GetRandom();
    double e1_toyq = quarks->GetRandom();
    double e2_toyg = gluons->GetRandom();
    double e2_toyq = quarks->GetRandom();
    if (myrand->Uniform(0,1) > gfrac1) gluons_toy->Fill(e1_toyq);
    else gluons_toy->Fill(e1_toyg);
    if (myrand->Uniform(0,1) > gfrac2) quarks_toy->Fill(e2_toyq);
    else quarks_toy->Fill(e2_toyg);  
   }
   for (int i=1; i<=gluons_toy->GetNbinsX(); i++){
    //std::cout << i << " " << gluons_toy->GetBinContent(i) << " " << ntoys*(gluons->GetBinContent(i)*gfrac1/gluons->Integral()+(1.-gfrac1)*quarks->GetBinContent(i)/quarks->Integral()) << std::endl;
    double myn1 = ntoys*(gluons->GetBinContent(i)*gfrac1/gluons->Integral()+(1.-gfrac1)*quarks->GetBinContent(i)/quarks->Integral());
    double myn2 = ntoys*(gluons->GetBinContent(i)*gfrac2/gluons->Integral()+(1.-gfrac2)*quarks->GetBinContent(i)/quarks->Integral());
    gluons_toy->SetBinContent(i,myrand->Poisson(myn1));
    quarks_toy->SetBinContent(i,myrand->Poisson(myn2));
    gluons_toy->SetBinError(i,sqrt(gluons_toy->GetBinContent(i)));
    quarks_toy->SetBinError(i,sqrt(quarks_toy->GetBinContent(i)));
   }   
   /*
   for (int i=1; i<=gluons_toy->GetNbinsX(); i++){
    std::cout << i << " " << gluons_toy->GetBinContent(i) << " " << ntoys*(gluons->GetBinContent(i)*gfrac1/gluons->Integral()+(1.-gfrac1)*quarks->GetBinContent(i)/quarks->Integral()) << std::endl;
    gluons_toy->SetBinContent(i,ntoys*(gluons->GetBinContent(i)*gfrac1/gluons->Integral()+(1.-gfrac1)*quarks->GetBinContent(i)/quarks->Integral()));
    quarks_toy->SetBinContent(i,ntoys*(gluons->GetBinContent(i)*gfrac2/gluons->Integral()+(1.-gfrac2)*quarks->GetBinContent(i)/quarks->Integral()));
    gluons_toy->SetBinError(i,10);
    quarks_toy->SetBinError(i,10);
   }
   exit(1);
   */

   TH2F* chi2plot = new TH2F("","",100,0.05,alphamax,99,0,1);
   TH2F* chi2plotq = new TH2F("","",100,0.05,alphamax,99,0,1);
   TH1F* chi2plot_combined = new TH1F("","",100,0.05,alphamax);
   TH1F* chi2plot_combinedq = new TH1F("","",100,0.05,alphamax);
   TH1F* chi2plot_combinedg = new TH1F("","",100,0.05,alphamax);
   for (double gfrac = 0.; gfrac <= 1.0001; gfrac +=0.01){
     for (double as = 0.05; as<alphamax; as+=0.001){
      TH1F* mixture = new TH1F("","",bins[mycounter].size()-1,&bins[mycounter][0]);
      for (int i=1; i<=mixture->GetNbinsX(); i++){
        mixture->SetBinContent(i,gfrac * alpha_s_maps[as]->GetBinContent(i) / alpha_s_maps[as]->Integral() + (1.-gfrac) * alpha_s_maps_quark[as]->GetBinContent(i) / alpha_s_maps_quark[as]->Integral());
      }

      if (1==2 && (fabs(gfrac - gfrac1) < 0.0001) && (fabs(as - 0.1) < 0.0001)){

        std::cout << "HERE " << std::endl;
        std::cout << chi2(gluons_toy,mixture) << std::endl;
        gluons_toy->Scale(1./gluons_toy->Integral());
        mixture->Scale(1./mixture->Integral());
        for (int q=1; q<=gluons_toy->GetNbinsX(); q++){
          std::cout << "   " << q << " " << gluons_toy->GetBinContent(q) << " " << mixture->GetBinContent(q) << " " << alpha_s_maps[as]->GetBinContent(q) << " " << gluons->GetBinContent(q) <<" " << alpha_s_maps_quark[as]->GetBinContent(q) << " " << quarks->GetBinContent(q) << std::endl;
        }
        exit(1);

      }
      chi2plot->SetBinContent(chi2plot->GetXaxis()->FindBin(as),chi2plot->GetYaxis()->FindBin(gfrac),chi2(gluons_toy,mixture));
     }
   }
   for (double gfrac = 0.; gfrac <= 1.0001; gfrac +=0.01){
     for (double as = 0.05; as<alphamax; as+=0.001){
      TH1F* mixture = new TH1F("","",bins[mycounter].size()-1,&bins[mycounter][0]);
      for (int i=1; i<=mixture->GetNbinsX(); i++){
        mixture->SetBinContent(i,gfrac * alpha_s_maps[as]->GetBinContent(i) / alpha_s_maps[as]->Integral() + (1.-gfrac) * alpha_s_maps_quark[as]->GetBinContent(i) / alpha_s_maps_quark[as]->Integral());
      }
      chi2plotq->SetBinContent(chi2plotq->GetXaxis()->FindBin(as),chi2plotq->GetYaxis()->FindBin(gfrac),chi2(quarks_toy,mixture));
     }
   }   

   gPad->SetLogx(0);
   gPad->SetRightMargin(0.15);
   gPad->SetTopMargin(0.12);
   //gPad->SetLogz();
   chi2plot->GetZaxis()->SetRangeUser(0.01,1);
   chi2plot->GetXaxis()->SetTitle("#alpha_{s}");
   chi2plot->GetYaxis()->SetTitleOffset(1.6);
   chi2plot->GetZaxis()->SetTitleOffset(1.5);
   chi2plot->GetXaxis()->SetNdivisions(505);
   chi2plot->GetYaxis()->SetTitle("Gluon Fraction");
   chi2plot->GetZaxis()->SetTitle("#chi^{2} compatibility (probability)");
   chi2plot->Draw("colz");
   //l.DrawLatex(0.2,0.95,"#bf{#scale[0.5]{Softdrop @ NLL, 1704.02210 }}"); 
   //l.DrawLatex(0.2,0.92,"#bf{#scale[0.5]{p_{T} = 500 GeV, z_{cut} = "+TString::Format("%0.2f",zcut)+", #alpha_{s} = "+TString::Format("%0.1f",alphasMZ)+", f_{g,1} = "+TString::Format("%0.0f%%",100*gfrac1)+", f_{g,2} = "+TString::Format("%0.0f%%",100*gfrac2)+", 100k events}}");     
   
   //l.DrawLatex(0.16,0.95,"#bf{#scale[0.95]{Softdrop @ NLL, 1704.02210 }}"); 
   l.DrawLatex(0.16,0.95,"#bf{#scale[0.7]{p_{T} = 500 GeV, f_{g,1} = "+TString::Format("%0.0f%%",100*gfrac1)+", f_{g,2} = "+TString::Format("%0.0f%%",100*gfrac2)+", 100k events}}");     
   l.DrawLatex(0.16,0.9,"#bf{#scale[0.7]{#alpha = "+TString::Format("%0.2f",alpha)+", z_{cut} = "+TString::Format("%0.2f",zcut)+" and #beta = "+TString::Format("%0.1f",beta)+"}}"); 

   chi2plotq->GetZaxis()->SetRangeUser(0.01,1);
   chi2plotq->Draw("colzsame");
   TLine *myline = new TLine(0.1,0,0.1,1);
   myline->SetLineStyle(3);
   myline->Draw();
   TMarker *m = new TMarker(0.1, 0.2, 23); 
   m->SetMarkerColor(2);
   m->SetMarkerSize(1.5);
   m->Draw();

   TMarker *m2 = new TMarker(0.1, 0.8, 22); 
   m2->SetMarkerColor(2);
   m2->SetMarkerSize(1.5);
   m2->Draw();

   gPad->RedrawAxis();
   c1->Print("banana_alpha_"+TString::Format("%0.2f",alpha)+"beta_"+TString::Format("%0.2f",beta)+"_zcut_"+TString::Format("%0.2f",zcut)+TString::Format("%i",myseed)+".pdf");

   all_histosq.push_back(chi2plotq);
   all_histosg.push_back(chi2plot);

   for (int i=1; i<=chi2plot->GetNbinsX(); i++){
    double maxprobg = 0;
    double maxprobq = 0;
    for (int j=1; j<=chi2plot->GetNbinsY(); j++){
      if (j==chi2plot->GetYaxis()->FindBin(gfrac2)) chi2plot_combinedq->SetBinContent(i,pow(chi2plotq->GetBinContent(i,j),2));
      if (j==chi2plot->GetYaxis()->FindBin(gfrac1)) chi2plot_combinedg->SetBinContent(i,pow(chi2plot->GetBinContent(i,j),2));
      if (chi2plot->GetBinContent(i,j) > maxprobg) maxprobg = chi2plot->GetBinContent(i,j);
      if (chi2plotq->GetBinContent(i,j) > maxprobq) maxprobq = chi2plotq->GetBinContent(i,j);
    }
    chi2plot_combined->SetBinContent(i,maxprobg*maxprobq); 
   }

   chi2plot_combined->Scale(1./chi2plot_combined->Integral());
   chi2plot_combinedg->Scale(1./chi2plot_combinedg->Integral());
   chi2plot_combinedq->Scale(1./chi2plot_combinedq->Integral());

   TH1F* chi2plot_combined_copy = mycopy(chi2plot_combined);
   TH1F* chi2plot_combinedg_copy = mycopy(chi2plot_combinedg);
   TH1F* chi2plot_combinedq_copy = mycopy(chi2plot_combinedq);

   chi2plot_combined_copy->SetLineColor(1);
   chi2plot_combined_copy->GetYaxis()->SetRangeUser(0.0,chi2plot_combinedg_copy->GetMaximum()*1.2);
   chi2plot_combined_copy->GetXaxis()->SetTitle("#alpha_{s}");
   chi2plot_combined_copy->GetYaxis()->SetTitleOffset(1.6);
   chi2plot_combined_copy->GetYaxis()->SetTitle("p(alpha_{s})");
   chi2plot_combined_copy->GetXaxis()->SetNdivisions(505);
   chi2plot_combined_copy->Draw();
   chi2plot_combinedg_copy->SetLineColor(2);
   chi2plot_combinedq_copy->SetLineColor(4);
   chi2plot_combinedq_copy->SetLineWidth(3);
   chi2plot_combined_copy->SetLineWidth(3);
   chi2plot_combinedg_copy->SetLineWidth(3);   
   chi2plot_combinedg_copy->Draw("same");
   chi2plot_combinedq_copy->Draw("same");

   //l.DrawLatex(0.16,0.95,"#bf{#scale[0.95]{Softdrop @ NLL, 1704.02210 }}"); 
   l.DrawLatex(0.16,0.95,"#bf{#scale[0.7]{p_{T} = 500 GeV, f_{g,1} = "+TString::Format("%0.0f%%",100*gfrac1)+", f_{g,2} = "+TString::Format("%0.0f%%",100*gfrac2)+", 100k events}}");     
   l.DrawLatex(0.16,0.9,"#bf{#scale[0.7]{#alpha = "+TString::Format("%0.2f",alpha)+", z_{cut} = "+TString::Format("%0.2f",zcut)+" and #beta = "+TString::Format("%0.1f",beta)+"}}");  

   //l.DrawLatex(0.2,0.95,"#bf{#scale[0.5]{Softdrop @ NLL, 1704.02210 }}"); 
   //l.DrawLatex(0.2,0.92,"#bf{#scale[0.5]{p_{T} = 500 GeV, z_{cut} = "+TString::Format("%0.2f",zcut)+", #alpha_{s} = "+TString::Format("%0.1f",alphasMZ)+", f_{g,1} = "+TString::Format("%0.0f%%",100*gfrac1)+", f_{g,2} = "+TString::Format("%0.0f%%",100*gfrac2)+", 100k events}}");     
   
   TLegend* leg2 = new TLegend(.2,.6,0.55,.85);
   leg2->SetTextFont(42);
   leg2->SetHeader("");
   leg2->SetNColumns(1);
   leg2->AddEntry(chi2plot_combined_copy,"Unknown Combined","l");
   leg2->AddEntry(chi2plot_combinedg_copy,"Known Gluon","l");
   leg2->AddEntry(chi2plot_combinedq_copy,"Known Quark","l");   
   leg2->SetFillStyle(0);
   leg2->SetFillColor(0);
   leg2->SetBorderSize(0);
   leg2->Draw();    

   c1->Print("palpha_alpha_"+TString::Format("%0.1f",alpha)+"beta_"+TString::Format("%0.1f",beta)+"_zcut_"+TString::Format("%0.3f",zcut)+TString::Format("%i",myseed)+".pdf");
 }
}
}

TH1F* chi2plot_combinedq = new TH1F("","",100,0.05,alphamax);
TH1F* chi2plot_combinedg = new TH1F("","",100,0.05,alphamax);
TH1F* chi2plot_combined = new TH1F("","",100,0.05,alphamax);

for (int ii=1; ii<=chi2plot_combined->GetNbinsX(); ii++){
  chi2plot_combinedq->SetBinContent(ii,1);
  chi2plot_combinedg->SetBinContent(ii,1);
  chi2plot_combined->SetBinContent(ii,1);
}

for (int i=0; i<all_histosq.size(); i++){
  for (int j=0; j<all_histosg.size(); j++){
    for (int ii=1; ii<=all_histosq[i]->GetNbinsX(); ii++){
        double maxprobg = 0;
        double maxpromybq = 0;
        for (int jj=1; jj<=all_histosq[i]->GetNbinsY(); jj++){
          if (jj==all_histosq[i]->GetYaxis()->FindBin(gfrac2)) chi2plot_combinedq->SetBinContent(ii,chi2plot_combinedq->GetBinContent(ii)*pow(all_histosq[i]->GetBinContent(ii,jj),2));
          if (jj==all_histosg[j]->GetYaxis()->FindBin(gfrac1)) chi2plot_combinedg->SetBinContent(ii,chi2plot_combinedg->GetBinContent(ii)*pow(all_histosg[j]->GetBinContent(ii,jj),2));
          if (all_histosg[j]->GetBinContent(ii,jj) > maxprobg) maxprobg = all_histosg[j]->GetBinContent(ii,jj);
          if (all_histosq[i]->GetBinContent(ii,jj) > maxpromybq) maxpromybq = all_histosq[i]->GetBinContent(ii,jj);
        }
        chi2plot_combined->SetBinContent(ii,maxprobg*maxpromybq*chi2plot_combined->GetBinContent(ii)); 
       }
  }
}

   chi2plot_combined->Scale(1./chi2plot_combined->Integral());
   chi2plot_combinedg->Scale(1./chi2plot_combinedg->Integral());
   chi2plot_combinedq->Scale(1./chi2plot_combinedq->Integral());

   TH1F* chi2plot_combined_copy = mycopy(chi2plot_combined);
   TH1F* chi2plot_combinedg_copy = mycopy(chi2plot_combinedg);
   TH1F* chi2plot_combinedq_copy = mycopy(chi2plot_combinedq);

   gPad->SetTopMargin(0.12);
   chi2plot_combined_copy->SetLineColor(1);
   chi2plot_combined_copy->GetYaxis()->SetRangeUser(0.0,chi2plot_combinedg_copy->GetMaximum()*1.2);
   chi2plot_combined_copy->GetXaxis()->SetTitle("#alpha_{s}");
   chi2plot_combined_copy->GetYaxis()->SetTitleOffset(1.6);
   chi2plot_combined_copy->GetYaxis()->SetTitle("p(alpha_{s})");
   chi2plot_combined_copy->GetXaxis()->SetNdivisions(505);
   chi2plot_combined_copy->Draw();
   chi2plot_combinedg_copy->SetLineColor(2);
   chi2plot_combinedq_copy->SetLineColor(4);
   chi2plot_combinedq_copy->SetLineWidth(3);
   chi2plot_combined_copy->SetLineWidth(3);
   chi2plot_combinedg_copy->SetLineWidth(3);
   chi2plot_combinedg_copy->Draw("same");
   chi2plot_combinedq_copy->Draw("same");
   //l.DrawLatex(0.16,0.95,"#bf{#scale[0.95]{Softdrop @ NLL, 1704.02210 }}"); 
   l.DrawLatex(0.16,0.95,"#bf{#scale[0.7]{p_{T} = 500 GeV, f_{g,1} = "+TString::Format("%0.0f%%",100*gfrac1)+", f_{g,2} = "+TString::Format("%0.0f%%",100*gfrac2)+", 100k events}}");     
   l.DrawLatex(0.16,0.9,"#bf{#scale[0.7]{sum over z_{cut} and #beta, #alpha_{s} = "+TString::Format("%0.1f",alphasMZ)+"}}");        

   TLegend* leg2 = new TLegend(.2,.6,0.55,.85);
   leg2->SetTextFont(42);
   leg2->SetHeader("");
   leg2->SetNColumns(1);
   leg2->AddEntry(chi2plot_combined_copy,"#scale[1.5]{Unknown Combined}","l");
   leg2->AddEntry(chi2plot_combinedg_copy,"#scale[1.5]{Uknown Gluon}","l");
   leg2->AddEntry(chi2plot_combinedq_copy,"#scale[1.5]{Uknown Quark}","l");   
   leg2->SetFillStyle(0);
   leg2->SetFillColor(0);
   leg2->SetBorderSize(0);
   leg2->Draw();  

   c1->Print("combination"+TString::Format("%i",myseed)+".pdf");  

}




