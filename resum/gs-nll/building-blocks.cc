#include "building-blocks.hh"

#include <cmath>
#include <cassert>
#include <limits>
#include <iostream>

using namespace std;

//----------------------------------------------------------------------
// internal flags and helpers
//----------------------------------------------------------------------
const double go_below_freeze = 1.0; // 1 for continuous freezing, 0 for a sharp cut

#define lam(X) (2.0*alphas_ref*b0*(X))


//========================================================================
// sudakov frequent bits
//
// Note: to obtain fixed coupling (at scale kt), set Lfr=0
// i.e. mufr=1.0 and the appropriate alphas value in the ConfigBase.
//========================================================================


// this is the expression for alphas suited for NLL resummation
double ConfigBase::alphas(double log_kt, bool allow_freezing) const{
  double lkt_eff = (allow_freezing && (log_kt>Lfr)) ? Lfr : log_kt;
  double den = 1.0-lam(lkt_eff);

  // if we're sitting below the lamdau pole, return "infinity"
  if (den<=0) return numeric_limits<double>::max();

  return alphas_ref/den
    * (two_loops
       ? 1.0+alphas_ref*(-b1/b0*log(den) + 0.5*Keff/M_PI)/den
       : 1.0);
}


// this is the expression for alphas at 1 loop suited for NLL resummation
double ConfigBase::alphas_1loop(double log_kt) const{
  double den = 1.0-lam(log_kt);

  // if we're sitting below the lamdau pole, return "infinity"
  return (den<=0) ? numeric_limits<double>::max() : alphas_ref/den;
}


double ConfigBase::delta_kt_line(double alpha, double Ltop, double Lbot, double B) const{
  assert(alpha!=1);

  // if teh line is upside-down, just revert it
  if (Lbot<Ltop) return delta_kt_line(alpha, Lbot, Ltop, B);
  assert(Lbot-Ltop>=-1e-8);

  // for the Bi term, we need alpha>0
  if (std::abs(B)>1e-4) assert(alpha>0);
  
  // handle the Landau pole if we do not freeze the coupling
  if ((!freeze) && (2*alphas_ref*b0*Lbot>=1.0)) return 0.0;

  // determine the scale that optionally enters in the B term
  double LB = ((alpha>0) && (alpha<1)) ? Lbot : Ltop;
  if (LB>Lfr) LB=Lfr;

  // now the results
  if (Lbot<Lfr){
    double lbot = lam(Lbot);
    double ltop = lam(Ltop);
    return CR/(M_PI*b0*std::abs(alpha-1))*
      (log((1-ltop)/(1-lbot))
       +(two_loops
         ? alphas_ref*Keff/(2*M_PI)*(lbot-ltop)/((1-ltop)*(1-lbot))
         - alphas_ref*b1/b0*( log(1-lbot)/(1-lbot) - log(1-ltop)/(1-ltop)
                          + (lbot-ltop)/((1-ltop)*(1-lbot)) )
         : 0.0))
      +2*alphas_1loop(LB)*CR/(alpha*M_PI)*B;
  }

  if (Ltop<Lfr){
    // split into a contribution with RC from Lfr to Ltop and a frozen
    // contribution between Lbot and Lfr
    //
    // Note that we need to be careful to discard te B term if it
    // falls below the freezing region and we're actually not going
    // there!
    double lbot = lam(Lfr);  // trick to get the same RC expression as above
    double ltop = lam(Ltop);
    return CR/(M_PI*b0*std::abs(alpha-1))*
      (log((1-ltop)/(1-lbot))
       +(two_loops
         ? alphas_ref*Keff/(2*M_PI)*(lbot-ltop)/((1-ltop)*(1-lbot))
         - alphas_ref*b1/b0*( log(1-lbot)/(1-lbot) - log(1-ltop)/(1-ltop)
                        + (lbot-ltop)/((1-ltop)*(1-lbot)) )
         : 0.0))
      +go_below_freeze*2*alphas(Lfr)*CR/(std::abs(alpha-1)*M_PI)*(Lbot-Lfr)
      +(LB<Lfr || go_below_freeze ? 1.0 : 0.0)*2*alphas_1loop(LB)*CR/(alpha*M_PI)*B;
  }
  
  return go_below_freeze*2*CR/M_PI
    *( alphas(Lfr)*(Lbot-Ltop)/std::abs(alpha-1) + alphas_1loop(Lfr)*B/alpha );
}


// same as above, specified on a basis of a line given by the
// coordinates of its endpoints. This also supports a constant kt
// line.
//
// The line starts from (theta0, kt0) to (theta1, kt1)
// Note that this is theta, not theta^2
//
// Note also that there is no check performed if z=1 is included,
// i.e. the B term is included blindly and it it the end-user's
// responsibility to set B to 0 when needed.
double ConfigBase::line(double Ltheta0, double Lkt0, double Ltheta1, double Lkt1, double B) const{
  
  // figure out the slope
  double alpha = 1-(Lkt1-Lkt0)/(Ltheta1-Ltheta0);
  if (std::abs(alpha-1)>1e-4)
    return delta_kt_line(alpha, Lkt0, Lkt1, B);

  double Las = Lkt0 < Lfr ? Lkt0 : Lfr;
  return 2.0*CR/M_PI*(alphas(Las)*std::abs(Ltheta0-Ltheta1)+alphas_1loop(Las)*B);
}

// a triangle bound by an upper theta limit, a kt limit and a z
// theta^alpha limit
//
// For alpha>1, the kt limit is an upper one, the z theta^alpha limit
// is a lower one
// For alpha<1, the kt limit is a lower one, the z theta^alpha limit
// is an upper one
//
// For alpha=0 and B!=0, a B term is included. It is the user's
// responsibility to included it only when needed
//
// constraints: - alpha != 1
//              - Ltop<Lbot
//
// assumptions: the full triangle is in the kinematic acceptance
// [note however that we can trick it if we want to take
// differences of two such triangles where th epart outside of the
// soundary cancels]
double ConfigBase::kt_triangle(double alpha, double Ltop, double Lbot, double B) const{
  // consistency checks
  assert(alpha!=1.0);
  if (Lbot-Ltop<=-1e-8) cerr << "---> " << Ltop << " " << Lbot << endl;
  assert(Lbot-Ltop>=-1e-8);
  assert((std::abs(alpha)<1e-5) || (std::abs(B)<1e-5));

  // handle the Landau pole if we do not freeze the coupling
  if ((!freeze) && (2*alphas_ref*b0*Lbot>=1.0)) return numeric_limits<double>::max();

  if (alpha>1){
    // both limits above the freezing scale
    if (Lbot<Lfr){
      double lbot = lam(Lbot);
      double ltop = lam(Ltop);
      double lgbot = log(1-lbot);
      double lgtop = log(1-ltop);
      double lg = lgbot-lgtop;
      return CR/(2*M_PI*alphas_ref*b0*b0*(alpha-1))
        *((1-lbot)*lg + (lbot-ltop)
          + (two_loops
             ? 0.5*alphas_ref*Keff/M_PI*((ltop-lbot)/(1-ltop)-lg)
             - alphas_ref*b1/b0*(0.5*lgtop*lgtop-0.5*lgbot*lgbot-lgbot+(1-lbot)/(1-ltop)*lgtop+(ltop-lbot)/(1-ltop))
             : 0.0));
    }
    
    // sitting over the freezing scale
    if (Ltop<Lfr){
      double lbot = lam(Lbot);
      double lfr  = lam(Lfr);
      double ltop = lam(Ltop);
      double lgfr  = log(1-lfr);
      double lgtop = log(1-ltop);
      return CR/(2*M_PI*alphas_ref*b0*b0*(alpha-1))
        *((1-lbot)*(lgfr-lgtop) + (lfr-ltop)
        + (two_loops
           ? 0.5*alphas_ref*Keff/M_PI*((ltop-lfr)*(1-lbot)/((1-ltop)*(1-lfr))+(lgtop-lgfr))
           - alphas_ref*b1/b0*(0.5*lgtop*lgtop-0.5*lgfr*lgfr-(1-lbot)/(1-lfr)*lgfr+(1-lbot)/(1-ltop)*lgtop+(ltop-lfr)*(1-lbot)/((1-ltop)*(1-lfr)))
           : 0.0))
        +go_below_freeze*alphas(Lfr)*CR/(M_PI*(alpha-1))*(Lbot-Lfr)*(Lbot-Lfr);
    }
    
    // both limits below the freezing scale
    return go_below_freeze*alphas(Lfr)*CR/(M_PI*(alpha-1))*(Lbot-Ltop)*(Lbot-Ltop);
  }

  // now triangles pointing upwards
  if (Lbot<Lfr){
    double lbot = lam(Lbot);
    double ltop = lam(Ltop);
    double lgbot = log(1-lbot);
    double lgtop = log(1-ltop);
    double lg = lgtop-lgbot;
    double Bfact = 1.0+lam(B)/(1-ltop);
    return CR/(2*M_PI*alphas_ref*b0*b0*(1-alpha))
      *((1-ltop)*Bfact*lg + (ltop-lbot)
        + (two_loops
           ? 0.5*alphas_ref*Keff/M_PI*((lbot-ltop)/(1-lbot)-lg)
           - alphas_ref*b1/b0*(0.5*lgbot*lgbot-0.5*lgtop*lgtop
                         -lgtop+(1-ltop)/(1-lbot)*lgbot
                         +(lbot-ltop)/(1-lbot))
           : 0.0));
  }
    
  // sitting over the freezing scale
  if (Ltop<Lfr){
    double lbot = lam(Lfr);  // misnamed so we can reuse the same final expression as above
    double ltop = lam(Ltop);
    double lgbot = log(1-lbot);
    double lgtop = log(1-ltop);
    double lg = lgtop-lgbot;
    double Bfact = 1.0+lam(B)/(1-ltop);
    return CR/(2*M_PI*alphas_ref*b0*b0*(1-alpha))
      *((1-ltop)*Bfact*lg + (ltop-lbot)
        + (two_loops
           ? 0.5*alphas_ref*Keff/M_PI*((lbot-ltop)/(1-lbot)-lg)
           - alphas_ref*b1/b0*(0.5*lgbot*lgbot-0.5*lgtop*lgtop
                         -lgtop+(1-ltop)/(1-lbot)*lgbot
                         +(lbot-ltop)/(1-lbot))
           : 0.0))
      +go_below_freeze*CR/(M_PI*(1-alpha))*(alphas(Lfr)*(Lbot+Lfr-2*Ltop)+2*alphas_1loop(Lfr)*B)*(Lbot-Lfr);
  }


  // both limits below the freezing scale
  return go_below_freeze*CR/(M_PI*(1-alpha))*(Lbot-Ltop)*(alphas(Lfr)*(Lbot-Ltop)+2*alphas_1loop(Lfr)*B);
}

//------------------------------------------------------------------------
double ConfigBase::triangle(double alpha, double beta, double Ltop, double Lbot, double B) const{
  assert(beta>alpha);

  // start dealing with cases with simplifications
  if (std::abs(alpha-1)<1e-4) return kt_triangle(beta ,Ltop,Lbot,0.0);
  if (std::abs(beta -1)<1e-4) return kt_triangle(alpha,Ltop,Lbot,  B);

  // the intermediate point
  double Lmed = ((beta-1)*Ltop + (1-alpha)*Lbot)/(beta-alpha);

  if (alpha<1){
    if (beta>1)
      return kt_triangle(alpha,Ltop,Lmed,B)+kt_triangle(beta,Lmed,Lbot,0.0);

    // note that we force the 2nd B to 0. This could probably be done a bit more elegantly
    return kt_triangle(alpha,Ltop,Lmed,B)-kt_triangle(beta,Lbot,Lmed,0.0);
  }
  
  return kt_triangle(beta,Lmed,Lbot,0.0)-kt_triangle(alpha,Lmed,Ltop,0.0);
}


//------------------------------------------------------------------------
// the R' corresponding to the bottom-line of the triangle
double ConfigBase::triangle_prime(double alpha, double beta, double Ltop, double Lbot, double B) const{
  assert(beta>alpha);

  // the kt of the intermediate point
  double Lmed = ((beta-1)*Ltop + (1-alpha)*Lbot)/(beta-alpha);
  if (std::abs(beta-1)>1e-4)
    return delta_kt_line(beta, Lmed, Lbot, B);
  
  // angle of the intermediate point
  //
  // The top line has z theta^alpha = cst
  // The top can be taken to have Ltheta=0, Lz=Ltop
  //  => the theta of the intermediate point is
  //     (alpha-1) Ltheta + Lmed = Ltop
  double Ltheta = (Ltop-Lmed)/(alpha-1);
  return line(0.0, Lbot, Ltheta, Lmed, B);
}

