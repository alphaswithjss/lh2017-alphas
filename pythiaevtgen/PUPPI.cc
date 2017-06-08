#include "PUPPI.hh"

/** 
 * A function called on all the particles in the event to compute the median and left-RMS according 
 * to (2.5) and (2.6) in 1407.6013. the alpha value for each partilce is computed and stored in the 
 * user info associated with that particle via the .puppi() method.
 */
void PUPPI::compute_median_left_RMS(vector<fastjet::PseudoJet> & tot_particles) {

  const int M = tot_particles.size();
  vector<double> alphas;

  // loop over particles to compute alphas
  for (int i = 0; i < M; i++) {

    fastjet::PseudoJet & pj_i = tot_particles[i];
    double alphai(0);

    // restrict to central region of detector
    // default value is 5.0 to basically just avoid unnecessary computation with the beam remnants. 
    // we do not care about the puppi weight for these particles
    if (abs(pj_i.rap()) > P_CENT_MAX) continue;

    // loop over particles to do sum for each alpha
    for (int j = 0; j < M; j++) {

      fastjet::PseudoJet & pj_j = tot_particles[j];

      // skip if not charged LV particle
      if (!pj_j.user_info<MyUserInfo>().isCharged() || 
           pj_j.user_info<MyUserInfo>().vertex() > -1) continue;

      // skip if outside of annulus 
      double dist_ij = pj_i.delta_R(pj_j);
      if (dist_ij < _Rmin || dist_ij > _R0) continue;

      alphai += pj_j.pt() / dist_ij;
    }

    alphai = log(alphai);

    // update .puppi property of user info (requires copying user info)
    MyUserInfo* myui = new MyUserInfo(pj_i.user_info<MyUserInfo>());
    myui->puppi(alphai);
    pj_i.set_user_info(myui);

    // if this is charged pileup with a finite alpha value, record for later stats
    if (pj_i.user_info<MyUserInfo>().isCharged() && 
        pj_i.user_info<MyUserInfo>().vertex() > -1 &&
        !isinf(alphai)) alphas.push_back(alphai);
  }

  // compute median
  sort(alphas.begin(), alphas.end());
  double median_alpha(0), left_RMS(0);
  if (alphas.size() > 0) {
    int midi = floor(alphas.size()/2);
    if (alphas.size() % 2 == 1)
      median_alpha = alphas[midi];
    else {
      median_alpha = (alphas[midi] + alphas[midi-1]) / 2.0;
    }

    // compute left_RMS
    int left_RMS_counter(0);
    for (int i = 0; alphas[i] < median_alpha; i++, left_RMS_counter++) 
      left_RMS += pow(alphas[i] - median_alpha, 2);
    if (left_RMS_counter > 0) left_RMS = sqrt(left_RMS/left_RMS_counter);
    else left_RMS = 0.1;
  }
  else {
    // to deal with the case of no pileup
    median_alpha = -100000000000000000000000.;
    left_RMS = 1;
  }

  // store values
  _median_alpha = median_alpha;
  _left_RMS = left_RMS;
  _has_median_left_RMS = true;
}

/**
 * A function to calculate the PUPPI weight for a single particle, 
 * given the already computed mediand and left-RMS
 */
void PUPPI::compute_weight(fastjet::PseudoJet & pj) const {

  // require the median and left RMS to have been computed 
  assert(_has_median_left_RMS);

  // get stored alpha value for particle
  double alpha = pj.user_info<MyUserInfo>().puppi();

  MyUserInfo* myui = new MyUserInfo(pj.user_info<MyUserInfo>());
  double w;

  // if alpha is negative infinity or less than median, set weight to zero
  if (pj.user_info<MyUserInfo>().isCharged()) {
    if (pj.user_info<MyUserInfo>().vertex() == -1) w = 1;
      else w = 0;
  }
  else {
    if (isinf(alpha) || alpha < _median_alpha) w = 0;
    else {

        // compute weight according to (2.7) and (2.8) of 1407.6013
        w = ROOT::Math::chisquared_cdf(pow((alpha - _median_alpha)/_left_RMS, 2), 1);
        if (w < _wcut || (w*pj).pt() < _pTcut) w = 0;
    }
  }
  myui->puppi(w);
  pj.set_user_info(myui);
}
