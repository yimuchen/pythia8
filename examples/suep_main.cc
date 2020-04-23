// This is an example main pythia script to generate to generate a dark sector shower in a strongly coupled, quasi-conformal
// hidden valley, often referred to as "soft unclustered energy patterns (SUEP)" or "softbomb" events.
// The code for the dark shower itself is in suep_shower.cc

// The algorithm relies on arXiv:1305.5226. See arXiv:1612.00850 for a description of the model. 
// Please cite both papers when using this code.

// Written by Simon Knapen on 12/22/2019

// pythia headers
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"

// libraries needed for the softbomb part
#include <iostream>
#include <math.h>
#include <boost/math/tools/roots.hpp>
#include <boost/bind.hpp>
#include<string>
#include<stdio.h>
#include<stdlib.h>

// suep_shower header
#include "SuepHook.h"

using namespace Pythia8;

template <typename T> string tostr(const T& t) { 
   ostringstream os; 
   os<<t; 
   return os.str(); 
}

int main(int argc, char *argv[]) {
     
   // read model parameters from the command line
  if(!(argc==6)){
    std::cout << "I need the following arguments: M m T decaycard outputfilename randomseed\n";
    std::cout << "with\n";
    std::cout << "M: mass of heavy scalar\n";
    std::cout << "m: mass of dark mesons\n";
    std::cout << "T: Temperature parameter\n";
    std::cout << "outputfilename: filename where events will be written\n";
    std::cout << "randomseed: an integer, specifying the random seed\n";
    return 0;
  }
     
  
  // model parameters and settings
  float mh, mX,T;
  string seed, filename, cardfilename;    
  mh=atof(argv[1]);
  mX=atof(argv[2]);
  T=atof(argv[3]);
  filename=tostr(argv[4]);
  seed=tostr(argv[5]);    
  
  // number of events    
  int Nevents=100;    
    
  // Interface for conversion from Pythia8::Event to HepMC event.
  HepMC::Pythia8ToHepMC ToHepMC;

  // Specify file where HepMC events will be stored.
  HepMC::IO_GenEvent ascii_io(filename, std::ios::out);

  // We will run pythia twice: Once for Higgs production, and once to decay the final state mesons of softbomb.
  // We therefore need two different pythia objects, each with different settings
  Pythia pythia;
  
  //Settings for the production Pythia object, before softbomb shower
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("HiggsSM:all = on");
  pythia.readString("25:mayDecay = off"); // turn off SM decays for Higgs
  pythia.readString("25:m0 = "+tostr(mh)); // set the mass of "Higgs" scalar
  pythia.readString("25:m0 = "+tostr(mh)); // set the mass of "Higgs" scalar
  pythia.readString("999999:all = GeneralResonance void 0 0 0 "+tostr(mX)+" 0.001 0.0 0.0 0.0");
  pythia.readString("999999:oneChannel = 1 1.0 101 1 -1");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = "+seed); 
  pythia.readString("Next:numberShowEvent = 0");

  //imported from cmssw
  pythia.readString("Tune:pp 14"); 
  pythia.readString("Tune:ee 7"); 
  pythia.readString("MultipartonInteractions:ecmPow=0.1391"); 
  pythia.readString("PDF:pSet=17"); 
  pythia.readString("MultipartonInteractions:bProfile=2"); 
  pythia.readString("MultipartonInteractions:pT0Ref=2.306"); 
  pythia.readString("MultipartonInteractions:coreRadius=0.3755"); 
  pythia.readString("MultipartonInteractions:coreFraction=0.3269"); 
  pythia.readString("ColourReconnection:range=2.323"); 
  pythia.readString("SigmaTotal:zeroAXB=off"); 
  pythia.readString("SpaceShower:rapidityOrder=off"); 
  pythia.readString("SpaceShower:alphaSvalue=0.13"); 
  pythia.readString("TimeShower:alphaSvalue=0.13");
  pythia.readString("Tune:preferLHAPDF = 2"); 
  pythia.readString("Main:timesAllowErrors = 10000"); 
  pythia.readString("Check:epTolErr = 0.01"); 
  pythia.readString("Beams:setProductionScalesFromLHEF = off"); 
  pythia.readString("SLHA:keepSM = on"); 
  pythia.readString("SLHA:minMassSM = 1000."); 
  pythia.readString("ParticleDecays:limitTau0 = on"); 
  pythia.readString("ParticleDecays:tau0Max = 10"); 
  pythia.readString("ParticleDecays:allowPhotonRadiation = on");

  //temporary
  pythia.readString("PartonLevel:MPI = off");

  SuepHook* hook = new SuepHook(25,999999,T);
  pythia.setUserHooksPtr(hook);

  pythia.init();
 
  // Shortcuts
  Event& event = pythia.event;
  
  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < Nevents; ++iEvent) {
    if (!pythia.next()) continue; //generate an event with the production pythia kernel, higgs is stable.
    
    // Print out a few events
    if (iEvent<1){  
        event.list();
    }
    
    // Construct new empty HepMC event and fill it.
    // Units will be as chosen for HepMC build; but can be changed
    // by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)
    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::MM);
    ToHepMC.fill_next_event( pythia, hepmcevt );

    // Write the HepMC event to file. Done with it.
    ascii_io << hepmcevt;
    delete hepmcevt;

  // End of event loop.
  }
  // print the cross sections etc    
  pythia.stat();

  // Done.
  return 0;
}
