/*

  Code setup by E. Metodiev, B. Nachman, P. Komiske
 
 */

#include <fstream>
#include <iostream>
#include <string>
#include <ctime>
#include <cstdlib>

#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include "CmdLine.hh"

#define LV_JET_MIN_PT     100
#define RECO_MIN_PT       25
#define MIN_QUARK_PT      95
#define RAP_MAX           2.5
#define DEFAULT_JET_R     1.0
#define MWIDTH            1.0
#define REF_JET_R         0.4

#define PRINT_FREQ        500
#define SMUSH             1.e-60
#define OUTPRECISION      8
#define NEVENT            10000
#define DEFAULT_SEED      -1
#define DEFAULT_NPU       150

using namespace std;
using namespace fastjet;

int getSeed(int seed, int npu, string bf){
  if (seed > 0) return seed;

  // extract job id from filename
  reverse(bf.begin(), bf.end());
  string js = bf.substr(bf.find(".") + 1);
  js = js.substr(0,js.find("_"));
  stringstream ss(js);
  int job_id;
  ss >> job_id;
  //int job_id = stoi(js);

  // combine time seeding, npu, and job_id to ensure a good seed
  srand(time(NULL));
  return ((rand() + npu * 10000 + job_id) % 900000000);
}

template <typename T>
string To_String(T val)
{
  //this is to_string, but for older versions of C++
  stringstream stream;
  stream << val;
  return stream.str();
}

int main(int argc, char** argv) {
  
  CmdLine cmdline(argc, argv);
  
  // Settings
  int     nEvent       =   cmdline.value("-nev", NEVENT);
  int     seed         =   cmdline.value("-seed", DEFAULT_SEED);
  double  jet_R        =   cmdline.value("-R", DEFAULT_JET_R);
  int     NPU          =   cmdline.value("-npu", DEFAULT_NPU);
  double  rap_max      =   cmdline.value("-rap-max", RAP_MAX);
  double  min_pt    =   cmdline.value("-lv-min-pt", LV_JET_MIN_PT);
  double  reco_min_pt  =   cmdline.value("-reco-min-pt", RECO_MIN_PT);
  double  min_q_pt     =   cmdline.value("-min-q-pt", MIN_QUARK_PT);
  int     print_freq   =   cmdline.value("-print", PRINT_FREQ);
  bool    do_UE        =  !cmdline.present("-no-ue");
  bool    jet_only     =   cmdline.present("-jet-only");
  string  filename     =   cmdline.value<string>("-out", "test.txt");

  // output setup
  ofstream outstream(filename.c_str());
  outstream << "# " << cmdline.command_line() << endl;
  outstream << "# date: " << cmdline.time_stamp() << endl;
   
  cmdline.assert_all_options_used();
  
  // PYTHIA SETUP
  Pythia8::Pythia pythia8b;

  // Specify processes
  pythia8b.readString("HardQCD:all = on");                                                                                             
  pythia8b.readString("PhaseSpace:pTHatMin  = 200");
  pythia8b.readString("PhaseSpace:pTHatMax  = 220");

  pythia8b.readString("Beams:idA = 2212");
  pythia8b.readString("Beams:idB = 2212");
  pythia8b.readString("Beams:eCM = 13000");
  pythia8b.init();

  // Random seed
  seed = getSeed(seed, NPU, filename);
  outstream << "# seed = " << seed << endl;
  pythia8b.readString("Random:setSeed = on");
  pythia8b.readString("Random:seed = " + To_String(seed));

  // event listing
  pythia8b.readString("Next:numberShowEvent = 0");
  pythia8b.readString("Next:numberShowProcess = 0");
  pythia8b.readString("Next:numberShowInfo = 0");

  pythia8b.init();

  // EVENT LOOP
  for (int iEvent = 0; iEvent < nEvent;) {

    // JET SETUP
    vector<PseudoJet> particles, jetsAkT, jetsCA;
    JetDefinition jet_defAkT = JetDefinition(antikt_algorithm, jet_R);
    JetDefinition jet_defCA = JetDefinition(cambridge_algorithm, jet_R);
    Selector selector = SelectorAbsRapMax(rap_max) * SelectorNHardest(2);

    // LV particle loop
    if (!pythia8b.next()) continue;

    int id = 0;
    for (int i = 0; i < pythia8b.event.size(); i++) {

      if (!pythia8b.event[i].isFinal()     || pythia8b.event[i].idAbs() == 12 ||
	  pythia8b.event[i].idAbs() == 13 || pythia8b.event[i].idAbs() == 14 ||
	  pythia8b.event[i].idAbs() == 16) continue;
      
      PseudoJet p(pythia8b.event[i].px(), 
		  pythia8b.event[i].py(), 
		  pythia8b.event[i].pz(), 
		  pythia8b.event[i].e());
      p.reset_PtYPhiM(p.pt(), p.rap(), p.phi(), 0);
      particles.push_back(p);
      id++;
    }

    //Now, we just want to save event-level information.
    ClusterSequence clsq(particles, jet_defAkT);
    jetsAkT = sorted_by_pt(selector(clsq.inclusive_jets(min_pt)));
    if (jetsAkT.size() < 1) continue;

    ClusterSequence clsqCA(jetsAkT[0].constituents(), jet_defCA);
    jetsCA = sorted_by_pt(selector(clsqCA.inclusive_jets(min_pt)));

    double beta = 0.;
    double zcut = 0.1;
    
    fastjet::contrib::SoftDrop sd(beta, zcut);
    fastjet::PseudoJet sdj1 = sd(jetsCA[0]);

    iEvent++;
    if (iEvent % print_freq == 0) cout << "Generated " << iEvent << " events so far..." << endl;
    
  } // end event loop
	
  outstream << "# End of file confirmation\n";

  // Statistics
  pythia8b.stat();

  return 0;
}
