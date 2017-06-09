  //Based on Andy Buckley's Boost 2014 Analysis

  // -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Recluster.hh"
#include "fastjet/contrib/ModifiedMassDropTagger.hh"
#include "fastjet/contrib/SoftDrop.hh"
using namespace fastjet;
#include <string>



namespace Rivet {
  
  using namespace Cuts;
  
    // a helper that holds
    //  - a jet radius
    //  - a kappa value
    //  - a beta value
    //  - an associated histogram with lin binning
    //  - an associated histogram with log binning
  class HistogramHolder{
  public:
    HistogramHolder(Histo1DPtr h_in, Histo1DPtr h_log_in)
    : h(h_in), h_log(h_log_in){}
    Histo1DPtr h, h_log;
  };
  
  typedef pair<double, double> pair_double; // will hold a (kappa,beta) pair
  
  
    /// Standard jet radius used in this analysis (for both kT and anti-kT)
  
  class MC_LHJSS_Zjet : public Analysis {
    
  private:
      /// parameters
    const double BOSON_PTMIN=500.;     ///< minimal boson pt
    const double JET_PTMIN=500.; ///< min value of jet_pt/boson_pt
    const double DELTA_RAP_MAX_ZJET=1.0; ///< max rapidity difference between Z and jet
    const double PARTICLE_RAPMAX=2.5;    ///< maximal rapidity allowed for particles
    const double JET_RAPMAX=2.5;         ///< maximal rapidity allowed for jets
    const vector<double> zcuts={0.05,0.1,0.2};
    const vector<double> betas={0.,1.,2.};
    const vector<double> alphas={0.5,1.,2.};
      // Fix the anti-kt radius for now to 0.8.
    double const antiKT_R=0.8;
    
    const map<string,int> obsmap = {{"Angularity",1},{"theta_g",2}};
    Recluster ca_wta_recluster;
    SharedPtr<contrib::SoftDrop> softdrop;
    
  public:
    
      /// Constructor
    MC_LHJSS_Zjet()
    : Analysis("MC_LHJSS_Zjet")
    {}
    
      /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs(-PARTICLE_RAPMAX, PARTICLE_RAPMAX, 0.0*GeV);
        // for the Z boson (-> mumu)
      Cut lepton_cut =  Rivet::Cuts::Quantity::pT >= 5.0*GeV; // we already know they have |y|<2.5
      ZFinder zfinder_mm_dressed(fs, lepton_cut,
                                 PID::MUON, 66.0*GeV, 116.0*GeV, 0.1,
                                 ZFinder::CLUSTERNODECAY, ZFinder::NOTRACK);
      addProjection(zfinder_mm_dressed, "ZFinder_mm_dressed");
        // for the jets
      VetoedFinalState jet_input(fs);
      jet_input.vetoNeutrinos();
      jet_input.addVetoPairId(PID::MUON);
      addProjection(jet_input, "JET_INPUT");
        // mMDT
      ca_wta_recluster = Recluster(JetDefinition(cambridge_algorithm,
                                                 JetDefinition::max_allowable_R,
                                                 WTA_pt_scheme),
                                   false,
                                   Recluster::keep_only_hardest);
        // histogram bookings
      for ( auto Omap : obsmap ){
        for ( double zcut : zcuts ){
          for ( double beta : betas ){
            for ( double alpha : alphas ){
                // plain jet quantities
              _obs.push_back(HistogramHolder(bookHisto1D(Omap.first
                                                          +"_zcut_"+std::to_string((float)zcut)
                                                          +"_beta_"+std::to_string((float)beta)
                                                          +"_alpha_"+std::to_string((float)alpha),
                                                          100,0.,1),
                                              bookHisto1D("log"+Omap.first
                                                          +std::to_string((float)zcut)
                                                          +"_beta_"+std::to_string((float)beta)+"_alpha_"
                                                          +std::to_string((float)alpha),
                                                          logspace(100,0.0001,1.))));
              
              
                // control plots
                //h_delta_phi_Zjet.push_back(bookHisto1D("deltaphi_Zjet"+Qlab+Rlab, 100, 0.0, pi));
            }
          }
        }
      }
    };
    
    
      /// Perform the per-event analysis
    void analyze(const Event& e) {
        // see if we find a Z boson
      const ZFinder& zfinder  = applyProjection<ZFinder>(e, "ZFinder_mm_dressed"   );
      if (zfinder.bosons().size() != 1) return;
      
        // deduce the cut to apply on jets
      const PseudoJet zmom = zfinder.bosons()[0].pseudojet();
      double zpt = zmom.pt();
      if (zpt < BOSON_PTMIN) return;
      
        // by default we impose the Z pt cut at the lowest scale under consideration.
        // We check here if we pass more stringent constraints
      
        // a few shortcuts
      const double weight = e.weight();
      const VetoedFinalState &fs = applyProjection<VetoedFinalState>(e, "JET_INPUT");
      vector<PseudoJet> particles;
      particles.reserve(fs.particles().size());
      foreach (const Particle &p, fs.particles()){
        particles.push_back(p.pseudojet());
      }
      
      
        // do the anti-kt clustering
      double ptmin_jet = JET_PTMIN;
      JetDefinition jet_def(antikt_algorithm, antiKT_R);
      vector<PseudoJet> jets = (SelectorAbsRapMax(JET_RAPMAX) * SelectorPtMin(ptmin_jet))(jet_def(particles));
      if(jets.empty()) return;
      
        // select only the hardest jet
      PseudoJet orig_jet = (SelectorNHardest(1)(jets))[0];
      
        // require that the jet is within 1 unit in rapidity of the Z boson
      if (std::abs(zmom.rap()-orig_jet.rap())>DELTA_RAP_MAX_ZJET) return;///TODO ???
      
      for ( auto Omap : obsmap ){
        for ( double zcut : zcuts ){
          for ( double beta : betas ){
            for ( double alpha : alphas ){
              
              softdrop.reset(new contrib::SoftDrop(beta,zcut));
              softdrop->set_grooming_mode();
              softdrop->set_reclustering(false); // we'll recluster ourselves with C/A and the WTA axis
              
                // recluster the jet to get a broadening-free axis and apply grooming
              PseudoJet jet = ca_wta_recluster(orig_jet);
              PseudoJet softdrop_jet = (*softdrop)(jet);
              
                // now compute the angularities for the plain and groomed jets
              compute_and_record(jet,Omap.first,      zcut, beta, alpha, weight);
              compute_and_record(softdrop_jet,Omap.first, zcut, beta, alpha, weight);
              
              
            }
          }
        }
        
        
        
        
      }
    };
    
      /// Normalise histograms etc., after the run
    void finalize() {
      foreach (const HistogramHolder &ho, _obs){
        normalize(ho.h);
        normalize(ho.h_log);
      }

    };
    
    
  private:
    
    vector<HistogramHolder> _obs;
    
    

    
    
    
    void compute_and_record(const PseudoJet &jet,
                            string obs,
                            double zcut,
                            double beta,
                            double alpha,
                            double weight){
      
      
      assert(obsmap.find(obs)!=obsmap.end());
      
      
        //{{"Angularity",1},{"theta_g",2}};
      int tmp=obsmap.at(obs);
      switch (tmp) {
        case 1:
          assert(1==obsmap.at("Angularity"));
            //do something
          break;
          
          
        case 2:
          assert(1==obsmap.at("theta_g"));
            //do something
          break;
          
        default:
          assert(false);
          break;
      }
      
      /*
       
       // init angularity sums to 0
       for(unsigned int i=0; i<ngas;i++) gas[i]=0.0;
       
       // sum over the constituents
       foreach (const PseudoJet& p, jet.constituents()){
       double pt = p.pt();
       double theta = p.squared_distance(jet);
       for(unsigned int i=0; i<ngas;i++)
       gas[i]+=pow(pt, _kappa_betas[i].first) * pow(theta, 0.5*_kappa_betas[i].second);
       }
       
       // normalise and fill histograms
       for(unsigned int i=0; i<ngas;i++){
       gas[i]/=(pow(jet.pt(), _kappa_betas[i].first) * pow(R, _kappa_betas[i].second));
       
       // loop over all ptcuts we're satisfying
       for (unsigned int iQ=0;iQ<nQ;++iQ){
       unsigned int offset = ngas*(iQ*2*nRADII+R_offset);
       _gas[offset+i].h->fill(gas[i], weight);
       // for the log binning make sure we avoid taking log(0) and put
       // that in the most -ve bin.
       _gas[offset+i].h_log->fill(gas[i]>0 ? log(gas[i]) : 1e-5-LOG_SCALE_MAX, weight);
       }
       }
       }
       */
    };
    
    
    
      // Hook for the plugin system
    DECLARE_RIVET_PLUGIN(MC_LHJSS_Zjet);
    
  };
}
