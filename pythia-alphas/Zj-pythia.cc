/// 
/// Code to generate events for Jet Information project.
/// Frédéric Dreyer, Patrick Komiske, Eric Metodiev, Jesse Thaler
/// MIT, 2017
///
/// Events.cc generates pythia events for a given input process and
/// saves ECFG outputs and jet images.
///
/// usage:
///   ./events [-nev D] [-seed I] [-pthatmin D] [-noUE,noISR,noFSR,parton]
///            [-R D] [-rapmax D] [-ptjetmin D] [-out S,out-img S,out-ecfg S]
///            [-nmax I] [-npixel I] [-process I,uds,cb,gluon,W,Z,top,higgs]
///            


#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#include "CmdLine.hh"
#include <string>

using namespace Pythia8;
using namespace std;

// default definitions
#define PRINT_FREQ       1000
#define NEVS             1000
#define DEFAULT_SEED     -1
#define NPIX             33

#define JET_R            0.8
#define RAP_MAX          1.7
#define PT_HAT_MIN       150.0
#define PT_JET_MIN       200.0
#define PT_POW           -1.0

int main(int argc, char** argv) {
  CmdLine cmdline(argc, argv);

  double alphas   =  cmdline.value("-alphas",0.1365);
  double nEvents  =  cmdline.value("-nev",      NEVS);
  int seed        =  cmdline.value("-seed",     DEFAULT_SEED);
  double pthatmin =  cmdline.value("-pthatmin", PT_HAT_MIN);
  double ptpow    =  cmdline.value("-ptpow",    PT_POW);
  bool do_UE      = !cmdline.present("-noUE");
  bool do_hadr    = !cmdline.present("-parton");
  bool do_FSR     = !cmdline.present("-noFSR");
  bool do_ISR     = !cmdline.present("-noISR");
  string fn_out   =  cmdline.value<string>("-out","Zj-pythia.hepmc");
  string fn_in    =  cmdline.value<string>("-in","Zj.lhe");
  
  // Interface for conversion from Pythia8::Event to HepMC event.
  HepMC::Pythia8ToHepMC ToHepMC;
  
  // Specify file where HepMC events will be stored.
  HepMC::IO_GenEvent ascii_io(fn_out, std::ios::out);
  
  // Generator
  Pythia pythia;

  // start by setting up process settings
  bool read_from_lhe = false;
  if (cmdline.present("-Zj")) {
    // Take the Z+jet process
    pythia.readString("WeakZ0:gmZmode = 2");  
    // set up the quark or gluon channels
    if (!cmdline.present("-gluon")){
      pythia.readString("WeakBosonAndParton:qg2gmZq = on");
    }
    if (!cmdline.present("-quark")){
      pythia.readString("WeakBosonAndParton:qqbar2gmZg = on");
    }
  } else if (cmdline.present("-dijets")) {
    // Alternatively use dijets
    pythia.readString("HardQCD:all = on");
  } else {
    // Otherwise read from LHE file
    read_from_lhe = true;
    pythia.readString("Beams:frameType = 4");
    pythia.readString("Beams:LHEF = "+fn_in);
  }
  if (!read_from_lhe) {
    pythia.settings.flag("Random:setSeed", true);
    pythia.settings.mode("Random:seed", seed);
    pythia.settings.parm("PhaseSpace:pTHatMin", pthatmin);
    pythia.settings.parm("PhaseSpace:bias2SelectionPow", ptpow);
    pythia.settings.flag("PhaseSpace:bias2Selection", ptpow >= 0 ? true : false);
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eCM = 13000.");
  }

  // set the alphas coupling constant
  pythia.settings.parm("SigmaProcess:alphaSvalue", alphas);
  pythia.settings.parm("MultipartonInteractions:alphaSvalue", alphas);
  pythia.settings.parm("SpaceShower:alphaSvalue", alphas);
  pythia.settings.parm("TimeShower:alphaSvalue", alphas);
  
  // now set up shower settings
  pythia.settings.flag("PartonLevel:MPI", do_UE);
  pythia.settings.flag("PartonLevel:ISR", do_ISR);
  pythia.settings.flag("PartonLevel:FSR", do_FSR);
  pythia.settings.flag("HadronLevel:Hadronize", do_hadr);
  // turn off QED shower
  pythia.settings.flag("TimeShower:QEDshowerByQ", false);
  pythia.settings.flag("SpaceShower:QEDshowerByQ", false);
  pythia.settings.flag("TimeShower:QEDshowerByL", false);
  pythia.settings.flag("SpaceShower:QEDshowerByL", false);

  pythia.init();
  if (read_from_lhe)
    nEvents = pythia.mode("Main:numberOfEvents");

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nEvents; ++ iEvent) {
    if (!pythia.next()) continue;

    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
    ToHepMC.fill_next_event( pythia, hepmcevt );

    ascii_io << hepmcevt;
    delete hepmcevt;
  } // End of event loop.

  // Statistics
  pythia.stat();

  return 0;
}
