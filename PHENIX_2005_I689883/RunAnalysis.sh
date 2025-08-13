#!/bin/bash
rivet-build RivetPHENIX_2005_I689883.so PHENIX_2005_I689883.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -a PHENIX_2005_I689883 -o Output.yoda ../testfiles/events1.dat
