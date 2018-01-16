#include <cmath>
#include <cassert>
#include <sstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include "resum-SD.hh"
#include "building-blocks.hh"

using namespace std;

class QGVector{
public:
  QGVector() : q(0.0), g(0.0){}
  QGVector(double q_in, double g_in) : q(q_in), g(g_in){}
  
  double q,g;

  QGVector operator -() const{ return QGVector(-q, -g); }
  QGVector operator +(const QGVector &B) const{
    return QGVector(q+B.q, g+B.g);
  }
  QGVector operator -(const QGVector &B) const{
    return QGVector(q-B.q, g-B.g);
  }
  QGVector operator *(const QGVector &B) const{
    return QGVector(q*B.q, g*B.g);
  }
  QGVector operator *(const double &x) const{
    return QGVector(q*x, g*x);
  }
  QGVector operator /(const double &x) const{
    return QGVector(q/x, g/x);
  }

  // exponentiation
  QGVector exponential() const{
    return QGVector(exp(q), exp(g));
  }

  //--------------------------------------------------
  // apply a given quark fraction
  double apply_qfrac(double qfrac) const{
    return qfrac*q + (1-qfrac)*g;
  }
};

// build the special functions we'll need
QGVector exp(const QGVector &A){ return A.exponential(); }

QGVector lngamma(const QGVector &A){
  return QGVector(gsl_sf_lngamma(A.q), gsl_sf_lngamma(A.g));
}

/// simpe math
QGVector operator+(const double a, const QGVector &A){
  return QGVector(a+A.q, a+A.g);
}

QGVector operator*(const double a, const QGVector &A){
  return QGVector(a*A.q, a*A.g);
}

                  
/// overloaded output
std::ostream & operator<<(std::ostream & ostr, const QGVector &A) {
  ostr << A.q << " " << A.g;
  return ostr;
}

//------------------------------------------------------------------------
// SD resummation formulae
class Resum_SD{
public:
  // empty ctor (not to be used for physics)
  Resum_SD(){}

  // ctor with proper initialisation
  Resum_SD(double zcut, double beta_in, double alpha_in, double ptR_in, double alphasMZ_in, double muR_in, double kt_freeze=1.0)
    : alphasMZ(alphasMZ_in), muR(muR_in){
    ptR = ptR_in;
    alpha = alpha_in;
    assert(alpha>0);
    
    // get the alphas value (use the 2-loop expansion independently of
    // the choice for the resummation
    double alphas = alphas_rg(ptR*muR, 91.1876, alphasMZ);
    
    cfg_q = ConfigBase(CF, alphas, kt_freeze/ptR, muR);
    cfg_g = ConfigBase(CA, alphas, kt_freeze/ptR, muR);

    cfg_q.set_two_loops(true);
    cfg_g.set_two_loops(true);

    beta = beta_in;
    Lc = log(1.0/zcut);
  }

  //------------------------------------------------------------------------
  // results for the Sudakov exponent
  //------------------------------------------------------------------------
  //
  // for large rho, everything here is an approximation at best. But
  // it goes beyond our accuracy anyway! We'll therefore take a simple
  // plain jet, double log, approximation
  //
  // finitez controls whether or not finite zcut effects are included
  //
  // zinalphas  controls whether or not we include the z factor in alphas
  QGVector R(double rho) const{
    // compute Lrho (depends on rho and a possible endpoint rhomax)
    double Lrho;
    if (!compute_Lrho(rho,Lrho)) return QGVector();
    return QGVector(R_SD(Lrho, cfg_q), R_SD(Lrho, cfg_g));
  }

  // first order expansion in alphas
  QGVector R_LO(double rho) const{
    // compute Lrho (depends on rho and a possible endpoint rhomax)
    double Lrho;
    if (!compute_Lrho(rho,Lrho)) return QGVector();
    return QGVector(R_SD_LO(Lrho, cfg_q), R_SD_LO(Lrho, cfg_g));
  }

  // second order expansion in alphas
  //
  // this is the alphas expansion of the exponent, i.e. the beta0
  // terms.  The other alphas^2 contribution coming from the expansion
  // of the exponential will be taken care when we do the integration
  QGVector R_NLO(double rho) const{
    // compute Lrho (depends on rho and a possible endpoint rhomax)
    double Lrho;
    if (!compute_Lrho(rho,Lrho)) return QGVector();
    return QGVector(R_SD_NLO(Lrho, cfg_q), R_SD_NLO(Lrho, cfg_g));
  }

  QGVector Rp(double rho) const{
    // compute Lrho (depends on rho and a possible endpoint rhomax)
    double Lrho;
    if (!compute_Lrho(rho,Lrho)) return QGVector();
    return QGVector(Rp_SD(Lrho, cfg_q), Rp_SD(Lrho, cfg_g));
  }

  // first order expansion in alphas
  QGVector Rp_LO(double rho) const{
    // compute Lrho (depends on rho and a possible endpoint rhomax)
    double Lrho;
    if (!compute_Lrho(rho,Lrho)) return QGVector();
    return QGVector(Rp_SD_LO(Lrho, cfg_q), Rp_SD_LO(Lrho, cfg_g));
  }

  // second order expansion in alphas
  //
  // this is the alphas expansion of the exponent, i.e. the beta0
  // terms.  The other alphas^2 contribution coming from the expansion
  // of the exponential will be taken care when we do the integration
  QGVector Rp_NLO(double rho) const{
    // compute Lrho (depends on rho and a possible endpoint rhomax)
    double Lrho;
    if (!compute_Lrho(rho,Lrho)) return QGVector();
    return QGVector(Rp_SD_NLO(Lrho, cfg_q), Rp_SD_NLO(Lrho, cfg_g));
  }

  //----------------------------------------------------------------------
  // Physics building blocks
  //
  // individual flavour contributions for SD radiators, their
  // derivatives and their expansions
  //----------------------------------------------------------------------

  // radiator
  double R_SD(double Lrho, const ConfigBase &cfg) const{
    return ((Lrho>-cfg.Bi) ? cfg.triangle(  0.0, alpha, -cfg.Bi, Lrho, 0.0) : 0.0)
         - ((Lrho>Lc)      ? cfg.triangle(-beta, alpha,      Lc, Lrho, 0.0) : 0.0);
  }

  // LO expansion of the radiator
  double R_SD_LO(double Lrho, const ConfigBase &cfg) const{
    if (abs(alpha-2.0)>0.001){ return 0.0; }
    double B = cfg.Bi;
    return cfg.alphas_ref*cfg.CR/M_PI*
      (((Lrho>-B) ? (Lrho+B) *(Lrho+B) *0.5 : 0.0)
      -((Lrho>Lc) ? (Lrho-Lc)*(Lrho-Lc)/(2.0+beta) : 0.0));
  }

  // NLO expansion of the radiator
  double R_SD_NLO(double Lrho, const ConfigBase &cfg) const{
    if (abs(alpha-2.0)>0.001){ return 0.0; }
    double B = cfg.Bi;
    double K = cfg.Keff;
    return cfg.alphas_ref*cfg.alphas_ref*cfg.CR/M_PI*
      (b0*(((Lrho>-B) ? 0.5*(Lrho+B)*(Lrho+B)*(Lrho-B) : 0.0)
           -((Lrho>Lc) ? 2.0/(3.0*(2.0+beta)*(2.0+beta))*(Lrho-Lc)*(Lrho-Lc)
             *((2*beta+3.0)*Lrho+(3.0+beta)*Lc) : 0.0))
       +K/(2.0*M_PI)*(((Lrho>-B) ? 0.5*(Lrho+B)*(Lrho+B) : 0.0)
                     -((Lrho>Lc) ? (Lrho-Lc)*(Lrho-Lc)/(2.0+beta) : 0.0)));
  }

  // emittor (derivative of the radiator)
  double Rp_SD(double Lrho, const ConfigBase &cfg_base) const{
    // we should not include the 2-loop corrections
    ConfigBase cfg = cfg_base;
    cfg.set_two_loops(false);

    if (abs(alpha-1.0)<1e-4){
      // z t = rho
      // z = zc t^beta
      //   => t = (rho/zc)^{1/(a+b)}
      double Ltmax = (Lrho>Lc) ? (Lrho-Lc)/(alpha+beta) : 0.0;
      return (Lrho>-cfg.Bi) ? cfg.line(Ltmax, Lrho, Lrho+cfg.Bi, Lrho, 0.0) : 0.0;
    }

    // the intersection between the mass line and the SD line is at 
    //  log(1/z theta) = (Lc+(1+b) Lrho)/(2+b)
    // for alpha!=2, we have
    //  z theta^a = rho
    //  z = zc theta^b
    //  => t = (rho/zc)^{1/(a+b)}
    //  => kt = zc (rho/zc)^{(1+b)/(a+b)}
    //        = zc^((a-1)/(a+b)) rho^((1+b)/(a+b))
    double Lktmax = (Lrho>Lc) ? ((alpha-1.0)*Lc+(1.0+beta)*Lrho)/(alpha+beta) : Lrho;
    double result = (Lrho>-cfg.Bi) ? cfg.delta_kt_line(alpha, 0.5*(Lrho-cfg.Bi), Lktmax, 0.0) : 0.0;
    
    // restore 2-loop contributions
    //cfg.set_two_loops(true);

    return result;
  }

  // LO expansion of the emittor
  double Rp_SD_LO(double Lrho, const ConfigBase &cfg) const{
    if (abs(alpha-2.0)>0.001){ return 0.0; }
    double B = cfg.Bi;
    return cfg.alphas_ref*cfg.CR/M_PI*
      (((Lrho>-B) ? (Lrho+B) : 0.0) - ((Lrho>Lc) ? (Lrho-Lc)*2.0/(2.0+beta) : 0.0));
  }

  // NLO expansion of the emittor
  //
  // Note: no K required for our use of R' (as a NLL contribution to
  // the Sudakov)
  double Rp_SD_NLO(double Lrho, const ConfigBase &cfg) const{
    if (abs(alpha-2.0)>0.001){ return 0.0; }
    double B = cfg.Bi;
    return cfg.alphas_ref*cfg.alphas_ref*cfg.CR/M_PI*b0*
      (((Lrho>-B) ? 0.5*(Lrho+B)*(3*Lrho-B) : 0.0)
      -((Lrho>Lc) ? 2.0/((2.0+beta)*(2.0+beta))*(Lrho-Lc)*((2*beta+3.0)*Lrho+Lc) : 0.0));
  }

  // help[er to build Lrho including the endpoint
  //......................................................................
  bool compute_Lrho(double rho, double &Lrho) const{
    double rho_natural = exp(cfg_q.Bi);
    double Lrho_natural = log(1.0/rho_natural);
    if (rhomax>0){
      if (rho>=rhomax){ return false;}
      Lrho = log(1.0/rho + (1.0/rho_natural-1.0/rhomax));
    } else {
      Lrho = log(1.0/rho);
    }

    return (Lrho>Lrho_natural);
  }


  double alphasMZ, muR;
  double rhomax;

  double ptR, Lc, beta, alpha;
  ConfigBase cfg_q, cfg_g;
};

//------------------------------------------------------------------------
// compute the weights for a given pt
void weights_SD(double pt, double R, double zcut, double beta, double alpha,
                double alphasMZ, double muR, double muC, double mufr,
                const vector<double> & logrhos,
                double endpoint,
                vector<double> &weights_q_resum,
                vector<double> &weights_g_resum,
                vector<double> &weights_q_lo,
                vector<double> &weights_g_lo,
                vector<double> &weights_q_nlo, 
                vector<double> &weights_g_nlo){
  Resum_SD sd_resum(zcut, beta, alpha, pt*R, alphasMZ, muR, mufr);
  sd_resum.rhomax = endpoint;

  // resummation-scale uncertainty
  //
  // note that rho will be multiplied by muC, so the endpoint on
  // rho*mu8C should be rhomax*muC instead of rhomax.
  // This can simply be done by multiplying rhomax by muC as well
  sd_resum.rhomax *= muC;
  double lc = log(muC);

  // compute the cumulative distributions at all the l10 points
  vector<double> res_q; res_q.reserve(logrhos.size());
  vector<double> lo_q ; lo_q .reserve(logrhos.size());
  vector<double> nlo_q; nlo_q.reserve(logrhos.size());
  vector<double> res_g; res_g.reserve(logrhos.size());
  vector<double> lo_g ; lo_g .reserve(logrhos.size());
  vector<double> nlo_g; nlo_g.reserve(logrhos.size());
  
  //for (double lrho : logrhos){
  for (unsigned int i=0;i<logrhos.size(); ++i){
    double lrho = logrhos[i];
    double rho = exp(-lrho) * muC;

    QGVector R     = sd_resum.R(rho);     // resummed radiator
    QGVector RLO   = sd_resum.R_LO(rho);  // LO  expansion
    QGVector RNLO  = sd_resum.R_NLO(rho); // NLO expansion

    QGVector Rp    = sd_resum.Rp(rho);     // resummed emittor
    QGVector RpLO  = sd_resum.Rp_LO(rho);  // LO  expansion
    QGVector RpNLO = sd_resum.Rp_NLO(rho); // NLO expansion

    // construct the resummed result
    QGVector resum = exp(-R-M_EULER*Rp-lngamma(1+Rp)-lc*Rp);

    // exp(-R-M_EULER*Rp-lngamma(1+Rp)-lc*Rp)-1
    //   = - R-(M_EULER*Rp+lngamma(1+Rp))-lc*Rp
    //   + 0.5*(R+(M_EULER*Rp+lngamma(1+Rp))+lc*Rp)^2
    //   = as * (-RLO-lc RpLO)
    //   + as^2 (-RNLO-lc*RpNLO + 0.5*(RLO+lc*RpLO)^2 - M_PI^2*RpLO*RpLO/12.0)
    QGVector v_lo  = -RLO-lc*RpLO;
    QGVector v_nlo = -RNLO-lc*RpNLO + 0.5*(RLO+lc*RpLO)*(RLO+lc*RpLO)-M_PI*M_PI/12.0*RpLO*RpLO;

    res_q.push_back(resum.q);
    res_g.push_back(resum.g);

    lo_q .push_back(v_lo .q);
    lo_g .push_back(v_lo .g);

    nlo_q.push_back(v_nlo.q);
    nlo_g.push_back(v_nlo.g);
  }

  // compute diffs
  //
  // Note that this just computes the probabilities of being in each
  // bins! I.e. we should not divide by the bin width (this is taken
  // care of when we make the histogram later)
  weights_q_resum.resize(logrhos.size()-1);
  weights_g_resum.resize(logrhos.size()-1);
  weights_q_lo   .resize(logrhos.size()-1);
  weights_g_lo   .resize(logrhos.size()-1);
  weights_q_nlo  .resize(logrhos.size()-1);
  weights_g_nlo  .resize(logrhos.size()-1);
  for (unsigned int i=0; i<logrhos.size()-1; ++i){
    weights_q_resum[i]=(res_q[i]-res_q[i+1]);
    weights_g_resum[i]=(res_g[i]-res_g[i+1]);
    weights_q_lo   [i]=(lo_q [i]-lo_q [i+1]);
    weights_g_lo   [i]=(lo_g [i]-lo_g [i+1]);
    weights_q_nlo  [i]=(nlo_q[i]-nlo_q[i+1]); 
    weights_g_nlo  [i]=(nlo_g[i]-nlo_g[i+1]); 
  }

}

