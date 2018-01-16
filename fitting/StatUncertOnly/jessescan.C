//Formulae for the Soft drop mass.
//Distinction between quarks and gluons is encoded in the arguments C_i and B_i for the different functions.
//For quarks, C_i->CF, B_i->Bq,
//For gluons, C_i->CA, B_i->Bg,
//all of which are constants defined below.
//Several test calls are shown at the bottom for some specific values.
//---------------------------------------------------------------------------------------------------

#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <map>
#include <iostream>

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

TH1F* mycopy(TH1F* input){
  TH1F* output = new TH1F("","",100,0.05,0.2);
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

//-------------------------------------------------------------------------------------------------------------
//test evaluating functions at several points

void mymain() {

   //test alpha_s
   printf("alpha_s = %f\n", alpha_s(150,alpha_s_mZ_standard));  

   double gfrac1 = 0.8;
   double gfrac2 = 0.2;
   double myalphas = 0.1;

   const int nzcuts = 3;
   const int nbetas = 1;
   double zcuts[nzcuts] = {0.05,0.1,0.2};
   double betas[nbetas] = {0.};

   vector<double> bins;
   bins.push_back(0.05);
   while (bins[bins.size()-1] > 1e-4){
    bins.push_back(bins[bins.size()-1]/1.1);
   }
   std::reverse(bins.begin(), bins.end());

   TLatex l  = TLatex(); 
   l.SetNDC();
   l.SetTextColor(1);

   TCanvas *c1 = new TCanvas("","",500,500);
   gPad->SetLogx();
   gPad->SetBottomMargin(0.15);
   gPad->SetLeftMargin(0.15);
   gStyle->SetOptStat(0);

   vector<TH2F*> all_histosq;
   vector<TH2F*> all_histosg;

   for (int iz = 0; iz<nzcuts; iz++){
    for (int ib = 0; ib<nbetas; ib++){

      double zcut = zcuts[iz];
      double beta = betas[ib];
   
       TH1F* gluons = new TH1F("","",bins.size()-1,&bins[0]);
       TH1F* quarks = new TH1F("","",bins.size()-1,&bins[0]);
       for (int i=1; i<=gluons->GetNbinsX(); i++){
          double e2 = gluons->GetXaxis()->GetBinCenter(i);
          gluons->SetBinContent(i, soft_drop_mass(e2,CA,Bg,zcut,myalphas,500.0,0.000001));
          quarks->SetBinContent(i, soft_drop_mass(e2,CF,Bq,zcut,myalphas,500.0,0.000001));    
       }

       gluons->GetYaxis()->SetNdivisions(505);
       gluons->GetXaxis()->SetTitleOffset(1.4);
       gluons->GetYaxis()->SetTitleOffset(1.6);
       gluons->GetYaxis()->SetRangeUser(0,gluons->GetMaximum()*1.5);
       gluons->GetXaxis()->SetTitle("log((m/p_{T})^{2})");
       gluons->GetYaxis()->SetTitle("d#sigma / dlog(e_{2})");
       gluons->SetLineColor(2);
       gluons->Draw();
       quarks->SetLineColor(4);
       quarks->SetLineStyle(3);
       quarks->Draw("same");
       l.DrawLatex(0.2,0.85,"#bf{#scale[0.5]{Softdrop @ MLL, A. Larkoski et al. 1402.2657}}"); 
       l.DrawLatex(0.2,0.8,"#bf{#scale[0.5]{p_{T} = 500 GeV, z_{cut} = "+TString::Format("%0.1f",zcut)+"}}");    
       TLegend* leg = new TLegend(.7,.75,0.85,.95);
       leg->SetTextFont(42);
       leg->SetHeader("");
       leg->SetNColumns(1);
       leg->AddEntry(gluons,"gluons","l");
       leg->AddEntry(quarks,"quarks","l");
       leg->SetFillStyle(0);
       leg->SetFillColor(0);
       leg->SetBorderSize(0);
       leg->Draw();    
       c1->Print("PDFs_"+TString::Format("%0.1f",zcut)+"_"+TString::Format("%0.3f",zcut)+".pdf");
    }
  }

 for (int iz = 0; iz<nzcuts; iz++){
  for (int ib = 0; ib<nbetas; ib++){

    double zcut = zcuts[iz];
    double beta = betas[ib];

   TH1F* gluons = new TH1F("","",bins.size()-1,&bins[0]);
   TH1F* quarks = new TH1F("","",bins.size()-1,&bins[0]);
   for (int i=1; i<=gluons->GetNbinsX(); i++){
      double e2 = gluons->GetXaxis()->GetBinCenter(i);
      gluons->SetBinContent(i, soft_drop_mass(e2,CA,Bg,zcut,myalphas,500.0,0.000001));
      quarks->SetBinContent(i, soft_drop_mass(e2,CF,Bq,zcut,myalphas,500.0,0.000001));    
   }

   //Now, let's start with statistical uncertainty.
   int ntoys = 100000;
   //For each pseudo-experiment, fit alpha_s.
   //First step: generate templates in alpha_s.  Let's start by using only gluons.
   std::map<double,TH1F*> alpha_s_maps;
   std::map<double,TH1F*> alpha_s_maps_quark;
   for (double as = 0.05; as<0.2; as+=0.001){
      TH1F* gluon_template = new TH1F(TString::Format("0.2%f",as),TString::Format("0.2%f",as),bins.size()-1,&bins[0]);
      TH1F* quark_template = new TH1F(TString::Format("0.2%f",as),TString::Format("0.2%f",as),bins.size()-1,&bins[0]);
      for (int i=1; i<=gluon_template->GetNbinsX(); i++){
        double e2 = gluon_template->GetXaxis()->GetBinCenter(i);
        gluon_template->SetBinContent(i, soft_drop_mass(e2,CA,Bg,zcut,as,500.0,0.000001)); 
        quark_template->SetBinContent(i, soft_drop_mass(e2,CF,Bq,zcut,as,500.0,0.000001));  
        gluon_template->SetBinError(i, 0.); 
        quark_template->SetBinError(i, 0.);  
      }
      alpha_s_maps[as]=gluon_template;
      alpha_s_maps_quark[as]=quark_template;
   }
   //Now, the psuedo-experiments.
   TRandom3 *myrand = new TRandom3(2345);
   TH1F* gluons_toy = new TH1F("","",bins.size()-1,&bins[0]);
   TH1F* quarks_toy = new TH1F("","",bins.size()-1,&bins[0]);
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
   TH2F* chi2plot = new TH2F("","",100,0.05,0.2,99,0,1);
   TH2F* chi2plotq = new TH2F("","",100,0.05,0.2,99,0,1);
   TH1F* chi2plot_combined = new TH1F("","",100,0.05,0.2);
   TH1F* chi2plot_combinedq = new TH1F("","",100,0.05,0.2);
   TH1F* chi2plot_combinedg = new TH1F("","",100,0.05,0.2);
   for (double gfrac = 0.; gfrac <= 1; gfrac +=0.01){
     for (double as = 0.05; as<0.2; as+=0.001){
      TH1F* mixture = new TH1F("","",bins.size()-1,&bins[0]);
      for (int i=1; i<=mixture->GetNbinsX(); i++){
        mixture->SetBinContent(i,gfrac * alpha_s_maps[as]->GetBinContent(i) / alpha_s_maps[as]->Integral() + (1.-gfrac) * alpha_s_maps_quark[as]->GetBinContent(i) / alpha_s_maps_quark[as]->Integral());
      }
      chi2plot->SetBinContent(chi2plot->GetXaxis()->FindBin(as),chi2plot->GetYaxis()->FindBin(gfrac),chi2(gluons_toy,mixture));
     }
   }
   for (double gfrac = 0.; gfrac <= 1; gfrac +=0.01){
     for (double as = 0.05; as<0.2; as+=0.001){
      TH1F* mixture = new TH1F("","",bins.size()-1,&bins[0]);
      for (int i=1; i<=mixture->GetNbinsX(); i++){
        mixture->SetBinContent(i,gfrac * alpha_s_maps[as]->GetBinContent(i) / alpha_s_maps[as]->Integral() + (1.-gfrac) * alpha_s_maps_quark[as]->GetBinContent(i) / alpha_s_maps_quark[as]->Integral());
      }
      chi2plotq->SetBinContent(chi2plotq->GetXaxis()->FindBin(as),chi2plotq->GetYaxis()->FindBin(gfrac),chi2(quarks_toy,mixture));
     }
   }   
   gPad->SetLogx(0);
   gPad->SetRightMargin(0.15);
   //gPad->SetLogz();
   chi2plot->GetZaxis()->SetRangeUser(0.01,1);
   chi2plot->GetXaxis()->SetTitle("#alpha_{s}");
   chi2plot->GetYaxis()->SetTitleOffset(1.6);
   chi2plot->GetZaxis()->SetTitleOffset(1.5);
   chi2plot->GetYaxis()->SetTitle("Gluon Fraction");
   chi2plot->GetZaxis()->SetTitle("#chi^{2} compatibility (probability)");
   chi2plot->Draw("colz");
   l.DrawLatex(0.2,0.95,"#bf{#scale[0.5]{Softdrop @ MLL, 1402.2657}}"); 
   l.DrawLatex(0.2,0.92,"#bf{#scale[0.5]{p_{T} = 500 GeV, z_{cut} = "+TString::Format("%0.1f",zcut)+", #alpha_{s} = "+TString::Format("%0.1f",myalphas)+", f_{g,1} = "+TString::Format("%0.0f%%",100*gfrac1)+", f_{g,2} = "+TString::Format("%0.0f%%",100*gfrac2)+", 100k events}}");     
   chi2plotq->GetZaxis()->SetRangeUser(0.01,1);
   chi2plotq->Draw("colzsame");
   TLine *myline = new TLine(0.1,0,0.1,1);
   myline->SetLineStyle(3);
   myline->Draw();
   c1->Print("jesseplot_"+TString::Format("%0.1f",zcut)+"_"+TString::Format("%0.3f",zcut)+".pdf");

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
   chi2plot_combined_copy->Draw();
   chi2plot_combinedg_copy->SetLineColor(2);
   chi2plot_combinedq_copy->SetLineColor(4);
   chi2plot_combinedg_copy->Draw("same");
   chi2plot_combinedq_copy->Draw("same");
   l.DrawLatex(0.2,0.95,"#bf{#scale[0.5]{Softdrop @ MLL, 1402.2657}}"); 
   l.DrawLatex(0.2,0.92,"#bf{#scale[0.5]{p_{T} = 500 GeV, z_{cut} = "+TString::Format("%0.1f",zcut)+", #alpha_{s} = "+TString::Format("%0.1f",myalphas)+", f_{g,1} = "+TString::Format("%0.0f%%",100*gfrac1)+", f_{g,2} = "+TString::Format("%0.0f%%",100*gfrac2)+", 100k events}}");     
   
   TLegend* leg2 = new TLegend(.5,.6,0.85,.93);
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

   c1->Print("jesseplot2"+TString::Format("%0.1f",zcut)+"_"+TString::Format("%0.3f",zcut)+".pdf");
 }
}

TH1F* chi2plot_combinedq = new TH1F("","",100,0.05,0.2);
TH1F* chi2plot_combinedg = new TH1F("","",100,0.05,0.2);
TH1F* chi2plot_combined = new TH1F("","",100,0.05,0.2);

for (int ii=1; ii<=chi2plot_combined->GetNbinsX(); ii++){
  chi2plot_combinedq->SetBinContent(ii,1);
  chi2plot_combinedg->SetBinContent(ii,1);
  chi2plot_combined->SetBinContent(ii,1);
}

for (int i=0; i<all_histosq.size(); i++){
  for (int j=0; j<all_histosg.size(); j++){
    for (int ii=1; ii<=all_histosq[i]->GetNbinsX(); ii++){
        double maxprobg = 0;
        double maxprobq = 0;
        for (int jj=1; jj<=all_histosq[i]->GetNbinsY(); jj++){
          if (jj==all_histosq[i]->GetYaxis()->FindBin(gfrac2)) chi2plot_combinedq->SetBinContent(ii,chi2plot_combinedq->GetBinContent(ii)*pow(all_histosq[i]->GetBinContent(ii,jj),2));
          if (jj==all_histosg[j]->GetYaxis()->FindBin(gfrac1)) chi2plot_combinedg->SetBinContent(ii,chi2plot_combinedg->GetBinContent(ii)*pow(all_histosg[j]->GetBinContent(ii,jj),2));
          if (all_histosg[j]->GetBinContent(ii,jj) > maxprobg) maxprobg = all_histosg[j]->GetBinContent(ii,jj);
          if (all_histosq[i]->GetBinContent(ii,jj) > maxprobq) maxprobq = all_histosq[i]->GetBinContent(ii,jj);
        }
        chi2plot_combined->SetBinContent(ii,maxprobg*maxprobq*chi2plot_combined->GetBinContent(ii)); 
       }
  }
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
   chi2plot_combined_copy->Draw();
   chi2plot_combinedg_copy->SetLineColor(2);
   chi2plot_combinedq_copy->SetLineColor(4);
   chi2plot_combinedg_copy->Draw("same");
   chi2plot_combinedq_copy->Draw("same");
   l.DrawLatex(0.2,0.95,"#bf{#scale[0.5]{Softdrop @ MLL, 1402.2657}}"); 
   l.DrawLatex(0.2,0.92,"#bf{#scale[0.5]{p_{T} = 500 GeV, sum over z_{cut}, #beta, f_{g,1} = "+TString::Format("%0.0f%%",100*gfrac1)+", f_{g,2} = "+TString::Format("%0.0f%%",100*gfrac2)+", 100k events}}");     
   
   TLegend* leg2 = new TLegend(.5,.6,0.85,.93);
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

   c1->Print("BLAH.pdf");  

}



