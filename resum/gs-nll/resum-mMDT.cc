#include <cmath>
#include <cassert>
#include <sstream>
#include "resum-mMDT.hh"
#include "building-blocks.hh"

using namespace std;

inline double W(double x){ return x*log(x);}

// a quick typedef for a matrix
//   qq  qg
//   gq  gg
class FlavourMatrix{
public:
  FlavourMatrix()
    : qq(0.0), qg(0.0), gq(0.0), gg(0.0){}
  
  FlavourMatrix(double qq_in, double qg_in, double gq_in, double gg_in)
    : qq(qq_in), qg(qg_in), gq(gq_in), gg(gg_in){}
  double qq, qg, gq, gg;

  FlavourMatrix operator -() const{
    return FlavourMatrix(-qq, -qg,
                         -gq, -gg);
  }

  FlavourMatrix operator +(const FlavourMatrix &B) const{
    return FlavourMatrix(qq+B.qq, qg+B.qg,
                         gq+B.gq, gg+B.gg);
  }

  FlavourMatrix operator -(const FlavourMatrix &B) const{
    return FlavourMatrix(qq-B.qq, qg-B.qg,
                         gq-B.gq, gg-B.gg);
  }

  FlavourMatrix operator *(const FlavourMatrix &B) const{
    return FlavourMatrix(qq*B.qq+qg*B.gq, qq*B.qg+qg*B.gg,
                         gq*B.qq+gg*B.gq, gq*B.qg+gg*B.gg);
  }

  FlavourMatrix operator *(const double &x) const{
    return FlavourMatrix(qq*x, qg*x,
                         gq*x, gg*x);
  }

  FlavourMatrix operator /(const double &x) const{
    return FlavourMatrix(qq/x, qg/x,
                         gq/x, gg/x);
  }


  // exponentiation
  // https://en.wikipedia.org/wiki/Matrix_exponential
  FlavourMatrix exponential() const{
    // define (half) the trace of the matrix
    double s = 0.5*(qq+gg);
    
    // A-sI = ((a-s,b)(c,d-s))
    // det = (a-s)*(d-s)-b*c
    assert(qg*gq-(qq-s)*(gg-s) >= 0);
    double q = sqrt(qg*gq-(qq-s)*(gg-s));

    double Afactor=(q==0) ? 1.0 : (exp(q)-exp(-q))/(2*q);  
    double Ifactor=(exp(q)+exp(-q))/2-s*Afactor;
    double prefactor = exp(s);

    return FlavourMatrix(prefactor*(qq*Afactor+Ifactor),
                         prefactor*(qg*Afactor),
                         prefactor*(gq*Afactor),
                         prefactor*(gg*Afactor+Ifactor));
  }

};
  
// build the exponential of a 2x2 matrix
FlavourMatrix exp(const FlavourMatrix &A){ return A.exponential(); }

/// simpe math
FlavourMatrix operator*(const double a, const FlavourMatrix &A){
  return FlavourMatrix(a*A.qq, a*A.qg, a*A.gq, a*A.gg);
}

/// overloaded output for PDFSets
std::ostream & operator<<(std::ostream & ostr, const FlavourMatrix &A) {
  ostr << A.qq << " " << A.qg << " " << A.gq << " " << A.gg;
  return ostr;
}

//------------------------------------------------------------------------
// mMDT resummation formulae
class Resum_mMDT{
public:
  // empty ctor (not to be used for physics)
  Resum_mMDT(){}

  // ctor with proper initialisation
  Resum_mMDT(double zcut, double ptR_in, double alphasMZ_in, double muR_in, double kt_freeze=1.0)
    : alphasMZ(alphasMZ_in), muR(muR_in){
    ptR = ptR_in;
    
    // get the alphas value (use the 2-loop expansion independently of
    // the choice for the resummation
    double alphas = alphas_rg(ptR*muR, 91.1876, alphasMZ);
    
    cfg_q = ConfigBase(CF, alphas, kt_freeze/ptR, muR);
    cfg_g = ConfigBase(CA, alphas, kt_freeze/ptR, muR);

    cfg_q.set_two_loops(false);
    cfg_g.set_two_loops(false);

    zc = zcut;
    Lc = log(1.0/zcut);
    Lm = log(1.0/(1.0-zcut));
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
  FlavourMatrix R(double rho) const{
    // compute Lrho (depends on rho and a possible endpoint rhomax)
    double Lrho;

    // the "natural" endpoint is either at Lrho=-cfg_q.Bi (no finite
    // z) or at 1/2 and we shioft it we shift to log(1/rhomax) [Note
    // that we should take the smallest of Bq and Bg]
    double rho_natural = exp(cfg_q.Bi);
    double Lrho_natural = log(1.0/rho_natural);

    if (rhomax>0){
      if (rho>=rhomax){ return FlavourMatrix();}
      Lrho = log(1.0/rho + (1.0/rho_natural-1.0/rhomax));
    } else {
      Lrho = log(1.0/rho);
    }

    if (Lrho<=Lrho_natural){ return FlavourMatrix();}


    // different treatment depending on the inclusion of finitez effects
    double Sq=0.0, Sg=0.0, Sqtog=0.0, Sgtoq=0.0;
      
    // we include the z factor in alphas
    // 1/z piece up to z=1/2
    Sq = Rplain(Lrho, -cfg_q.Bi, cfg_q) - ((Lrho>Lc) ? Rplain(Lrho, Lc, cfg_q) : 0.0);
    Sg = ((Lrho>-cfg_g.Bi) ? Rplain(Lrho, -cfg_g.Bi, cfg_g) : 0.0)
       - ((Lrho>Lc)        ? Rplain(Lrho, Lc,        cfg_g) : 0.0);
 
    // we include finite z corrections   
    if (Lrho>Lc){
      double as_int = int_alphas(Lrho, Lc);
      Sq += as_int*int_Pqtoq_nosing();
      Sg += as_int*int_Pgtog_nosing();
      
      Sqtog = as_int*int_Pqtog();
      Sgtoq = as_int*int_Pgtoq();
    }

    return FlavourMatrix(Sq+Sqtog,-Sgtoq,-Sqtog,Sg+Sgtoq);
  }

  // first order expansion in alphas
  FlavourMatrix R_LO(double rho) const{
    double rho_natural = exp(cfg_q.Bi);
    double Lrho_natural = log(1.0/rho_natural);
    double Lrho;
    if (rhomax>0){
      if (rho>=rhomax){ return FlavourMatrix();}
      Lrho = log(1.0/rho + (1.0/rho_natural-1.0/rhomax));
    } else {
      Lrho = log(1.0/rho);
    }

    if (Lrho<=Lrho_natural){ return FlavourMatrix();}

    // different treatment depending on the inclusion of finitez effects
    double Sq=0.0, Sg=0.0, Sqtog=0.0, Sgtoq=0.0;
      
    // we include the z factor in alphas
    // 1/z piece up to z=1/2
    Sq = Rplain_LO(Lrho, -cfg_q.Bi, cfg_q) - ((Lrho>Lc) ? Rplain_LO(Lrho, Lc, cfg_q) : 0.0);
    Sg = ((Lrho>-cfg_g.Bi) ? Rplain_LO(Lrho, -cfg_g.Bi, cfg_g) : 0.0)
       - ((Lrho>Lc)        ? Rplain_LO(Lrho, Lc,        cfg_g) : 0.0);
    
    // we include finite z corrections   
    if (Lrho>Lc){
      double as_int = int_alphas_LO(Lrho, Lc);
      Sq += as_int*int_Pqtoq_nosing();
      Sg += as_int*int_Pgtog_nosing();
      
      Sqtog = as_int*int_Pqtog();
      Sgtoq = as_int*int_Pgtoq();
    }

    return FlavourMatrix(Sq+Sqtog,-Sgtoq,-Sqtog,Sg+Sgtoq);
  }


  // second order expansion in alphas
  //
  // this is the alphas expansion of the exponent, i.e. the beta0
  // terms.  The other alphas^2 contribution coming from the expansion
  // of the exponential will be taken care when we do the integration
  FlavourMatrix R_NLO(double rho) const{
    double rho_natural = exp(cfg_q.Bi);
    double Lrho_natural = log(1.0/rho_natural);
    double Lrho;
    if (rhomax>0){
      if (rho>=rhomax){ return FlavourMatrix();}
      Lrho = log(1.0/rho + (1.0/rho_natural-1.0/rhomax));
    } else {
      Lrho = log(1.0/rho);
    }

    if (Lrho<=Lrho_natural){ return FlavourMatrix();}

    // different treatment depending on the inclusion of finitez effects
    double Sq=0.0, Sg=0.0, Sqtog=0.0, Sgtoq=0.0;
      
    // we include the z factor in alphas
    // 1/z piece up to z=1/2
    Sq = Rplain_NLO(Lrho, -cfg_q.Bi, cfg_q) - ((Lrho>Lc) ? Rplain_NLO(Lrho, Lc, cfg_q) : 0.0);
    Sg = ((Lrho>-cfg_g.Bi) ? Rplain_NLO(Lrho, -cfg_g.Bi, cfg_g) : 0.0)
       - ((Lrho>Lc)        ? Rplain_NLO(Lrho, Lc,        cfg_g) : 0.0);
    
    // we include the z factor in alphas
    if (Lrho>Lc){
      double as_int = int_alphas_NLO(Lrho, Lc);
      Sq += as_int*int_Pqtoq_nosing();
      Sg += as_int*int_Pgtog_nosing();
        
      Sqtog = as_int*int_Pqtog();
      Sgtoq = as_int*int_Pgtoq();
    }

    return FlavourMatrix(Sq+Sqtog,-Sgtoq,-Sqtog,Sg+Sgtoq);
  }

  //----------------------------------------------------------------------
  // Physics building blocks
  //----------------------------------------------------------------------

  // the integral of alphas(sqrt(rho)) between rho and rho0
  //......................................................................
  // We actually give Lrho = log(1/rho)
  //
  //                  (1)         (2)                   (3)
  // as(sqrt(rho)) = as/D - as^2 b1/b0 log(D)/D^2 + Keff/(2pi) as^2/D^2
  // D=1-as bo Lrho
  // \int_Lrho0^Lrho dLrho = 1/(as b0) \int_D1^D0 dD
  // Li=log(Di)
  //
  // (1) = 1/b0 log(D0/D1) = 1/b0 (L0-L1)
  // (2) = as b1/b0^2 [(L1+1)/D1-(L0+1)/D0]
  //
  // -int_1/b^1/a du log(u) =  (log(a)+1)/a - b...
  //
  // Resummed result
  double int_alphas(double Lrho, double Lrho0) const{
    double as = cfg_q.alphas_ref;
    double L0 = log(1-as*b0*Lrho0);

    if (0.5*Lrho>cfg_q.Lfr){
      double L1 = log(1-2.0*as*b0*cfg_q.Lfr);
      return 1.0/b0*(L0-L1)+cfg_q.alphas(cfg_q.Lfr)*(Lrho-2*cfg_q.Lfr);
    }
    
    double L1 = log(1-as*b0*Lrho);
    return 1.0/b0*(L0-L1);
  }

  // LO expansion
  double int_alphas_LO(double Lrho, double Lrho0) const{
    double as = cfg_q.alphas_ref;
    return as*(Lrho-Lrho0);
  }

  // NLO expansion
  double int_alphas_NLO(double Lrho, double Lrho0) const{
    double as = cfg_q.alphas_ref;
    // Note that we should shift K to include the possible feed-down
    // from muR!=1 in the leading term
    return 0.5*b0*(Lrho+Lrho0)*(Lrho-Lrho0)*as*as;
  }

  // integrals over the splitting functions between the appropriate
  // limits in the flavour matrix
  //......................................................................
  //
  // The "nosing" variants do not include the 1/z part of the
  // splitting function (neither for the log(1/zc) nor for the
  // log(1/(1-zc)))
  inline double int_Pqtoq() const{ return CF/M_PI*(Lc-Lm+(-0.75)*(1-2*zc)); }
  inline double int_Pgtog() const{ return CA/M_PI*(Lc-Lm-(1-2*zc)*(0.75+(1-2*nf*TR/CA)/6.0*(1-zc+zc*zc))); }
  inline double int_Pqtog() const{ return CF/M_PI*(Lm-0.5*zc*(1+0.5*zc)); }
  inline double int_Pgtoq() const{ return nf*TR/M_PI*zc*(1.0-zc+2.0*zc*zc/3.0); }
  inline double int_Pqtoq_nosing() const{ return CF/M_PI*(-Lm+1.5*zc); }
  inline double int_Pgtog_nosing() const{ return CA/M_PI*(-Lm+zc*(2-0.5*zc+zc*zc/3.0)-nf*TR/CA*zc*(1.0-zc+2.0*zc*zc/3.0)); }
  
  // plain jet mass Sudakov
  //......................................................................
  double Rplain(double Lrho, double Lrho0, const ConfigBase &cfg) const{
    return Lrho<=Lrho0 ? 0.0 : cfg.triangle(0, 2, Lrho0, Lrho, 0.0);
  }
  double Rplain_LO(double Lrho, double Lrho0, const ConfigBase &cfg) const{
    double as = cfg.alphas_ref;
    double CR = cfg.CR;
    double tmp=Lrho-Lrho0;
    return Lrho<=Lrho0 ? 0.0 : as*CR/(2.0*M_PI)*tmp*tmp;
  }
  double Rplain_NLO(double Lrho, double Lrho0, const ConfigBase &cfg) const{
    if (Lrho<=Lrho0) return 0.0;
    double as = cfg.alphas_ref;
    double CR = cfg.CR;
    double tmp=Lrho-Lrho0;
    double twopi = 2.0*M_PI;
    return as*as*CR/twopi * tmp*tmp * b0*(Lrho+Lrho0);
  }

  double alphasMZ, muR;
  double rhomax;

  double ptR, zc, Lc, Lm;
  ConfigBase cfg_q, cfg_g;
};

//------------------------------------------------------------------------
// compute the weights for a given pt
void weights_mMDT(double pt, double R, double zcut, double beta,
                  double alphasMZ, double muR, double muC, double mufr,
                  const vector<double> & logrhos,
                  double endpoint,
                  vector<double> &weights_q_resum,
                  vector<double> &weights_g_resum,
                  vector<double> &weights_q_lo,
                  vector<double> &weights_g_lo,
                  vector<double> &weights_q_nlo, 
                  vector<double> &weights_g_nlo){
  Resum_mMDT mmdt_resum(zcut, pt*R, alphasMZ, muR, mufr);
  mmdt_resum.rhomax = endpoint;

  // resummation-scale uncertainty
  //
  // note that rho will be multiplied by muC, so the endpoint on
  // rho*mu8C should be rhomax*muC instead of rhomax.
  // This can simply be done by multiplying rhomax by muC as well
  mmdt_resum.rhomax *= muC;

  // compute the cumulative distributions at all the l10 points
  vector<double> res_q; res_q.reserve(logrhos.size());
  vector<double> lo_q ; lo_q .reserve(logrhos.size());
  vector<double> nlo_q; nlo_q.reserve(logrhos.size());
  vector<double> res_g; res_g.reserve(logrhos.size());
  vector<double> lo_g ; lo_g .reserve(logrhos.size());
  vector<double> nlo_g; nlo_g.reserve(logrhos.size());

  for (double lrho : logrhos){
    double rho = exp(-lrho) * muC;

    FlavourMatrix R    = exp(-mmdt_resum.R(rho)); // resummed
    FlavourMatrix RLO  = -mmdt_resum.R_LO(rho);   // LO expansion
    FlavourMatrix RNLO0= mmdt_resum.R_NLO(rho);   // R at NLO
    FlavourMatrix RNLO = 0.5*(RLO*RLO)-RNLO0;     // full NLO expansion

    res_q.push_back(R   .qq+R   .gq);
    res_g.push_back(R   .qg+R   .gg);
    
    lo_q .push_back(RLO .qq+RLO .gq);
    lo_g .push_back(RLO .qg+RLO .gg);

    nlo_q.push_back(RNLO.qq+RNLO.gq);
    nlo_g.push_back(RNLO.qg+RNLO.gg); 
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

