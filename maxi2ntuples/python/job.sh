#!/bin/sh

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/m/molszews/CMSSW/tautau/M2N/CMSSW747allpairs/src/m2n/maxi2ntuples
cmsenv
cmsRun ./python/maxi2ntuples.py

#bsub -o log -N  job.sh
