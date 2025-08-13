// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include <cmath>
namespace Rivet {

  /// @brief This analysis compares NPE Invarient cross section between event generator and data, as well as charmed and bottomed electrons
  class PHENIX_2009_I816469 : public Analysis {
  public:
    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PHENIX_2009_I816469);

    /// @name Analysis methods
    /// @{
    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 0.35);
      declare(fs, "fs");

      // Book histograms
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      
      book(_h["All"], 3, 1, 1);
      book(_h["Charmed"], 4, 1, 1);
      book(_h["Bottom"], 5, 1, 1);

    }
    /// Perform the per-event analysis
    void analyze(const Event& event) {

      Particles fsParticles = applyProjection<FinalState>(event,"fs").particles();
      // Fill histogram
      for(const Particle& p : fsParticles){
        //Electrons and positrons
        if(p.pid()==-11 || p.pid() == 11){
          if(p.fromBottom() || p.fromCharm()){
            if(p.pT()>=1.8 && p.pT()<= 9.0){
              double binCenter;
              double binCenter2=5.63; 

              if(p.pT()<=2.0){
                binCenter = 1.9; 
              }else if(p.pT()<=2.2){
                binCenter = 2.1; 
                binCenter2 = 2.32; 
              }else if (p.pT()<=2.4){
                binCenter = 2.3; 
                binCenter2 = 2.32; 
              }else if (p.pT()<=2.6){
                binCenter = 2.5; 
                binCenter2 = 2.32; 
              }else if (p.pT()<=2.8){
                binCenter = 2.7; 
                binCenter2 = 2.32; 
              }else if (p.pT()<=3.0){
                binCenter = 2.9; 
                binCenter2 = 2.32; 
              }else if (p.pT()<=3.2){
                binCenter = 3.1; 
                binCenter2 = 3.35; 
              }else if (p.pT()<=3.4){
                binCenter = 3.3;
                binCenter2 = 3.35; 
              }else if (p.pT()<=3.6){
                binCenter = 3.5;
                binCenter2 = 3.35; 
              }else if (p.pT()<=3.8){
                binCenter = 3.7;
                binCenter2 = 3.35; 
              }else if (p.pT()<=4.0){
                binCenter = 3.9;
                binCenter2 = 3.35; 
              }else if (p.pT()<=4.2){
                binCenter = 4.1;
                binCenter2 = 4.38; 
              }else if (p.pT()<=4.4){
                binCenter = 4.3;
                binCenter2 = 4.38; 
              }else if (p.pT()<=4.6){
                binCenter = 4.5;
                binCenter2 = 4.38; 
              }else if (p.pT()<=4.8){
                binCenter = 4.7;
                binCenter2 = 4.38; 
              }else if (p.pT()<=5.0){
                binCenter = 4.9; 
                binCenter2 = 4.38; 
              }else if (p.pT()<=6.0){
                binCenter = 5.5;
                binCenter2 = 5.63;
              }else if (p.pT()<=7.0){
                binCenter = 6.5;
                binCenter2 = 5.63;
              }else if (p.pT()<=8.0){
                binCenter = 7.5;
              }else{
                binCenter = 8.5; 
              }
            
              _h["All"]->fill(p.pT()/GeV,1/binCenter);
              if(p.fromBottom()){
                _h["Bottom"]->fill(p.pT()/GeV,1/binCenter2);
              }else{
                _h["Charmed"]->fill(p.pT()/GeV,1/binCenter2);
              }

            } 
          }
        }
      }
    }
    /// Normalise histograms etc., after the run
    void finalize() {

      //Scale histogram to go from N to cross section
      
      scale(_h["All"], 42/(100000000*2*M_PI));
      scale(_h["Charmed"], 42/(100000000*2*M_PI));
      scale(_h["Bottom"], 42/(100000000*2*M_PI));


//        divide(_h["AAAA"], _s["AAAA_Scaled"], _h["AAAA"]);
//      _h["AAAA"] = _h["AAAA_scaled"];

//      normalize(_h["XXXX"]); // normalize to unity
  //    normalize(_h["AAAA"], crossSection()/millibarn); // normalize to generated cross-section in mb (no cuts)
//      scale(_h["ZZZZ"], crossSection()/picobarn/sumW()); // norm to generated cross-section in pb (after cuts)
    }
    /// @}
    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
//    map<string, Profile1DPtr> _p;
//    map<string, CounterPtr> _c;
    map<string, Scatter2DPtr> _s;

    /// @}
  };
  RIVET_DECLARE_PLUGIN(PHENIX_2009_I816469);
}