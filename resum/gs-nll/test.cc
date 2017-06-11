//NLL
#include "resum-SD.hh"
#include "resum-mMDT.hh"
#include <iostream>

using namespace std;


int main() {
  double pt = 500;
  double R = 0.8;
  double zcut = 0.2;
  double beta = 0.0;
  double alphasMZ = 0.118;
  double muR = 1.;
  double muC = 1.;
  double mufr = 1.; //GeV
  double endpoint = 0.279303; //1/4 for g -> qq; this is the maximal mass / pTR.  Small correction from that for not small angle approx.
  
  //Let's plot the soft-drop mass distribution
  unsigned int nbins=50;
  vector<double> bins(nbins+1, 0.0);
  for (unsigned int i=0;i<=nbins;++i){
    bins[i]=0.2*i;
  }
 

  std::vector<double> mmdt_weights_q_resum;
  std::vector<double> mmdt_weights_g_resum;
  std::vector<double> sd_weights_q_resum;
  std::vector<double> sd_weights_g_resum;
  std::vector<double> weights_q_lo;
  std::vector<double> weights_g_lo;
  std::vector<double> weights_q_nlo; 
  std::vector<double> weights_g_nlo;


  const std::vector<double> logrhos = bins;

  weights_mMDT(pt, R, zcut, alphasMZ, muR, muC, mufr,
               logrhos, endpoint, 
               mmdt_weights_q_resum,
               mmdt_weights_g_resum,
               weights_q_lo,
               weights_g_lo,
               weights_q_nlo, 
               weights_g_nlo);

  weights_SD(pt, R, zcut, beta, alphasMZ, muR, muC, mufr,
             logrhos, endpoint, 
             sd_weights_q_resum,
             sd_weights_g_resum,
             weights_q_lo,
             weights_g_lo,
             weights_q_nlo, 
             weights_g_nlo);

  for (unsigned int i=0;i<nbins;++i){
    cout << bins[i] << " " << bins[i+1] << " "
         << mmdt_weights_q_resum[i] << " " << mmdt_weights_g_resum[i] << " "
         << sd_weights_q_resum[i] << " " << sd_weights_g_resum[i] << endl;
  }

  return 0;

}



