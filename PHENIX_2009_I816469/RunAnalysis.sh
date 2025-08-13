#!/bin/bash
rivet-build RivetPHENIX_2009_I816469.so PHENIX_2009_I816469.cc
export RIVET_ANALYSIS_PATH=$PWD
rivet --pwd -a PHENIX_2009_I816469 -o Output.yoda ../testfiles/events1.dat
