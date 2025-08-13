// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include <cmath>
namespace Rivet {

  /// @brief This analysis compares NPE Invarient cross section between event generator and data
  class PHENIX_2005_I689883 : public Analysis {
  public:
    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PHENIX_2005_I689883);

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
      
      book(_h["AAAA"], 2, 1, 1);

    }
    /// Perform the per-event analysis
    void analyze(const Event& event) {

      Particles fsParticles = applyProjection<FinalState>(event,"fs").particles();
      // Fill histogram
      for(const Particle& p : fsParticles){
        //Filter for electrons and positrons
        if(p.pid()==-11 || p.pid() == 11){
          if(p.fromBottom() || p.fromCharm()){
            if(p.pT()>0.4 && p.pT()<=5){

              double binCenter; 
              //Dividing by the bin center
              if(p.pT()<=0.5){
                binCenter = 0.45; 
              }else if(p.pT()<=0.6){
                binCenter = 0.55; 
              }else if (p.pT()<=0.8){
                binCenter = 0.7; 
              }else if (p.pT()<=1.0){
                binCenter = 0.9; 
              }else if (p.pT()<=1.2){
                binCenter = 1.1; 
              }else if (p.pT()<=1.4){
                binCenter = 1.3; 
              }else if (p.pT()<=1.6){
                binCenter = 1.5; 
              }else if (p.pT()<=2.0){
                binCenter = 1.8;
              }else if (p.pT()<=2.5){
                binCenter = 2.25;
              }else if (p.pT()<=3.0){
                binCenter = 2.75;
              }else if (p.pT()<=4.0){
                binCenter = 3.5;
              }else{
                binCenter = 4.5;
              }


              _h["AAAA"]->fill(p.pT()/GeV,1/binCenter);

            }
          }
        }
      }
    }
    /// Normalise histograms etc., after the run
    void finalize() {

      //Scale histogram to go from N to cross section
      
      scale(_h["AAAA"], 42/(100000000*2*M_PI));

      
    }
    /// @}
    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;

    map<string, Scatter2DPtr> _s;

    /// @}
  };
  RIVET_DECLARE_PLUGIN(PHENIX_2005_I689883);
}