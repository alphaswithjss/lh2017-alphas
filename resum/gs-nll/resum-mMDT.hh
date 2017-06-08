#ifndef __RESUM_MMDT_HH__
#define __RESUM_MMDT_HH__

#include <vector>

/// compute the LL distribution (including finite zcut corrections)
/// using the modified Mass-Drop Tagger (mMDT) for a given kinematic
/// configuration. We also compute the expansion of the LL results at
/// LO and NLO
///
/// Arguments are as follows
///  Kinematic parameters:
///    \param pt  :  jet pt
///    \param R   : the jet radius
///    \param zcut: Soft-Drop zcut parameter
///  Scale choices:
///    \param alphasMZ: value of alphas at the Z mass (evolved w 2-loop running)
///    \param muR     : choice of renormalisation scale (factor of pt R)
///    \param muC     : choice of resummation scale (factor of pt R)
///    \param mufr    : choice of freezing scale (in GeV) (see notes bbelow)
///  Binning information:
///    \param logrhos : the bin edges in log(pt^2 R^2/m^2)
///    \param endpoint: the desired eendpoint of the distribution (see notes)
///  Results:
///    \param weight_q_resum: NLL distribution for quark-initiated jets
///    \param weight_g_resum: NLL distribution for gluon-initiated jets
///    \param weight_q_lo   : NLL@LO distribution for quark-initiated jets
///    \param weight_g_lo   : NLL@LO distribution for gluon-initiated jets
///    \param weight_q_nlo  : NLL@NLO distribution for quark-initiated jets
///    \param weight_g_nlo  : NLL@NLO distribution for gluon-initiated jets
///
/// Notes:
///
///  - freezing scale: it is just meant to guarantee a smooth
///    transition dow nto very small masses. They can be viewed as a
///    double-counting w MC-based NP hadronisation effects. Our final
///    answer should not depend on this
/// 
///  - endpoint: for R=0.8, this is 0.279303 for LO
///                                 0.44974  for NLO
///
///  - At this stage, the results are weights in each bin, i.e. not
///    normalised by the bin width
///
/// TODO:
///  - add support for angularities
///  - add support for single-emission variants
void weights_mMDT(double pt, double R, double zcut,
                  double alphasMZ, double muR, double muC, double mufr,
                  const std::vector<double> & logrhos,
                  double endpoint, 
                  std::vector<double> &weights_q_resum,
                  std::vector<double> &weights_g_resum,
                  std::vector<double> &weights_q_lo,
                  std::vector<double> &weights_g_lo,
                  std::vector<double> &weights_q_nlo, 
                  std::vector<double> &weights_g_nlo);

#endif // __RESUM_MMDT_HH__
