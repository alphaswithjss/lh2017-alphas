#include <cmath>

#include "fastjet/PseudoJet.hh"
//#include "root/Math/ProbFunc.h"
#include "Math/ProbFunc.h"

#include "FJUI.hh"

#ifndef __PUPPI_HH__
#define __PUPPI_HH__

#define P_R0             0.3
#define P_RMIN           0.02
#define P_WCUT           0.1
#define P_CENT_MAX       5.0
#define P_PT_CUT(npu)   (0.1 + npu * 0.007)

using namespace std;

class PUPPI {
public:
  PUPPI() {}
  PUPPI(int npu, double R0 = P_R0, double Rmin = P_RMIN, double wcut = P_WCUT) :
  _R0(R0), _Rmin(Rmin), _wcut(wcut), _pTcut(P_PT_CUT(npu)) {}
  void compute_median_left_RMS(vector<fastjet::PseudoJet> & tot_particles);
  void compute_weight(fastjet::PseudoJet & pj) const;

private:
  double _R0, _Rmin, _wcut, _pTcut, _median_alpha, _left_RMS;
  bool _has_median_left_RMS;
};

#endif // __PUPPI_HH__
