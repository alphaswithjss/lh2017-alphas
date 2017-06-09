#include "LHEF.h"
#include <cstdlib>
#include <string>
#include <iostream>

//ROOT stuff
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TCanvas.h"

void mymain(){

  double R = 0.8;

  // Init readers and writers
  LHEF::Reader* reader = 0;
  reader = new LHEF::Reader("cut.lhe"); //unweighted_events.lhe");

  vector<double> bins;
  bins.push_back(1);
  while (bins[bins.size()-1] > 1e-4){
    bins.push_back(bins[bins.size()-1]/1.1);
  }
  std::reverse(bins.begin(), bins.end());

  TH1F* myDRmin = new TH1F("","",20,0,1);
  TH1F* sdmass = new TH1F("","",bins.size()-1,&bins[0]);

  while (reader->readEvent()) {
    
    std::vector<TLorentzVector> particles;
    for (int i = 0; i < reader->hepeup.NUP; ++i) {
      TLorentzVector myv = TLorentzVector(reader->hepeup.PUP[i][0],reader->hepeup.PUP[i][1],reader->hepeup.PUP[i][2],reader->hepeup.PUP[i][3]);
      if (fabs(myv.Eta()) <  5) particles.push_back(myv); //can get the type by reader->hepeup.IDUP[i]
    }
    
    double mindR = 999.;
    TLorentzVector j1;
    TLorentzVector j2;
    for (int i=0; i<particles.size(); i++){
      for (int j=i+1; j<particles.size(); j++){
	if (particles[i].DeltaR(particles[j]) < mindR){
	  mindR = particles[i].DeltaR(particles[j]);
	  j1 = particles[i];
	  j2 = particles[j];
	}
      }
    }
    
    myDRmin->Fill(mindR);
    if (mindR < 0.8){
      sdmass->Fill(pow((j1+j2).M()/(j1+j2).Pt(),2));
      std::cout << (j1+j2).M()/(j1+j2).Pt() << std::endl;
    }
  }

  TCanvas* c1 = new TCanvas("","",500,500);
  myDRmin->Draw();
  c1->Print("DR.pdf");
  gPad->SetLogx();
  sdmass->Draw();
  c1->Print("sdmass_LO.pdf");

}
