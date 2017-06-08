/*

Formulae are from 1402.2657, coded up by Ian Moult.

*/

#include <stdio.h>
#include <stdarg.h>
#include <math.h>

//Plotting stuff
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TString.h"
#include "TGraphErrors.h"
#include "TLine.h"


//define constants

//Z mass
double mz=91.18;

//pi
double pi=3.141592654;

//Color charges
double CF=4.0/3.0;
double CA=3.0;

//definitions associated with alpha_s
double alpha_s_mZ=0.118;
double beta0=23.0/(12.0*pi);
double Bg=-11.0/12.0+5.0/18.0;
double Bq=-3.0/4.0;

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

//Derivative of Sudakov (here this is just a simple numerical derivative with step delta)
double diffsudakov(const double e2, const double Ci, const double Bi, const double zcut, const double sdbeta, const double as_mz, const double pt, const double delta){

	return (1/delta)*(sudakov(e2+delta,Ci,Bi,zcut,sdbeta,as_mz,pt)-sudakov(e2,Ci,Bi,zcut,sdbeta,as_mz,pt));

}


//cutoff
double cutoff(const double e2, const double zcut, const double pt){
        
        return step(e2-1/(zcut*pt*pt));

}



//result to plot
double result(const double e2, const double Ci, const double Bi, const double zcut, const double sdbeta, const double as_mz, const double pt, const double delta){

	return e2*diffsudakov(e2,Ci,Bi,zcut,sdbeta,as_mz,pt,delta)*cutoff(e2,zcut,pt)/(diffsudakov(zcut+0.0001,Ci,Bi,zcut,sdbeta,as_mz,pt,delta));

}

//-------------------------------------------------------------------------------------------------------------
//test evaluating functions at several points

void mymain() {

   //Let's plot the soft-drop mass distribution
   vector<double> bins;
   bins.push_back(1);
   while (bins[bins.size()-1] > 1e-4){
    bins.push_back(bins[bins.size()-1]/1.1);
   }
   std::reverse(bins.begin(), bins.end());

   TH1F* gluons_beta0 = new TH1F("","",bins.size()-1,&bins[0]);
   TH1F* quarks_beta0 = new TH1F("","",bins.size()-1,&bins[0]);

   TH1F* gluons_beta1 = new TH1F("","",bins.size()-1,&bins[0]);
   TH1F* quarks_beta1 = new TH1F("","",bins.size()-1,&bins[0]);

   for (int i=1; i<=gluons_beta0->GetNbinsX(); i++){
      double e2 = gluons_beta0->GetXaxis()->GetBinCenter(i);
      double zcut = 0.1;

      gluons_beta0->SetBinContent(i, result(e2,CA,Bg,zcut,0.,0.10,500.0,0.000001));
      quarks_beta0->SetBinContent(i, result(e2,CF,Bq,zcut,0.,0.10,500.0,0.000001)); 

      gluons_beta1->SetBinContent(i, result(e2,CA,Bg,zcut,1.,0.10,500.0,0.000001));
      quarks_beta1->SetBinContent(i, result(e2,CF,Bq,zcut,1.,0.10,500.0,0.000001));         
   }

   TCanvas *c1 = new TCanvas("","",500,500);
   gPad->SetLogx();
   gPad->SetBottomMargin(0.15);
   gPad->SetLeftMargin(0.15);
   gStyle->SetOptStat(0);
   gluons_beta0->GetYaxis()->SetNdivisions(505);
   gluons_beta0->GetXaxis()->SetTitleOffset(1.4);
   gluons_beta0->GetYaxis()->SetTitleOffset(1.6);
   gluons_beta0->GetYaxis()->SetRangeUser(0,0.2);
   gluons_beta0->GetXaxis()->SetTitle("m/p_{T}");
   gluons_beta0->GetYaxis()->SetTitle("e_{2} d#sigma / de_{2}");
   gluons_beta0->SetLineColor(2);
   gluons_beta0->Draw();
   quarks_beta0->SetLineColor(4);
   quarks_beta0->Draw("same");

   gluons_beta1->SetLineStyle(3);
   quarks_beta1->SetLineStyle(3);
   gluons_beta1->SetLineColor(2);
   gluons_beta1->Draw("same");
   quarks_beta1->SetLineColor(4);
   quarks_beta1->Draw("same");

   TLatex l  = TLatex(); 
   l.SetNDC();
   l.SetTextColor(1);
   l.DrawLatex(0.2,0.85,"#bf{#scale[0.5]{Softdrop @ MLL, 1402.2657}}"); 
   l.DrawLatex(0.2,0.8,"#bf{#scale[0.5]{p_{T} = 500 GeV, z_{cut} = 0.1}}");    
   TLegend* leg = new TLegend(.7,.73,0.85,.93);
   leg->SetTextFont(42);
   leg->SetHeader("");
   leg->SetNColumns(1);
   leg->AddEntry(gluons_beta0,"gluons (#beta = 0)","l");
   leg->AddEntry(quarks_beta0,"quarks (#beta = 0)","l");
   leg->AddEntry(gluons_beta1,"gluons (#beta = 1)","l");
   leg->AddEntry(quarks_beta1,"quarks (#beta = 1)","l");   
   leg->SetFillStyle(0);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->Draw();    
   c1->Print("sf_gluons.pdf");























   //test alpha_s
   printf("alpha_s = %f\n", alpha_s(150,0.118));  

   //test log
   printf("log = %f\n", xlog(0.1));

   //test radiator for a quark (use Ci=CF, and Bi=Bq)
   printf("Test Radiator Quark = %f\n", Rad(0.0001,CF,Bq,beta0,0.12,0.1,0.0));

   //test Sudakov for a quark (use Ci=CF, and Bi=Bq)
   printf("Test Sudakov Quark = %f\n", sudakov(0.0001,CF,Bq,0.1,0.0,0.12,500));

   //test diff Sudakov for a quark (use Ci=CF, and Bi=Bq)
   printf("Test diff Sudakov Quark = %f\n", diffsudakov(0.0001,CF,Bq,0.1,0.0,0.12,500,0.000001));

   //call it for a quark (use Ci=CF, and Bi=Bq)
   printf("Test Value Quark = %f\n", result(0.0001,CF,Bq,0.1,0.0,0.12,500.0,0.000001));

   //test diff Sudakov for a gluon (use Ci=CA, and Bi=Bg)
   printf("Test diff Sudakov Gluon = %f\n", diffsudakov(0.0001,CA,Bg,0.1,0.0,0.12,500,0.000001));

   //call it for a gluon (use CA, and Bg)
   printf("Test Value Gluon = %f\n", result(0.0001,CA,Bg,0.1,0.0,0.12,500.0,0.000001));


}



