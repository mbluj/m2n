import os
from subprocess import call

#Directory with input file(s). Do not put ".root" files there that are not maent to be processed.
directory = '/afs/cern.ch/work/m/molszews/CMSSW/Data/EmAOD_VBF/'
#directory = '/afs/cern.ch/work/m/molszews/CMSSW/tautau/LLRHiggsTauTau/CMSSW_7_4_5/src/LLRHiggsTauTau/NtupleProducer/test/test_F/crab_m2n/results/'


#Directory with outputfile(s). Can be of course the same as the above one, but remember to remove ouput files before another run.
outputdir = '/afs/cern.ch/work/m/molszews/CMSSW/Data/ntuple_VBF/' 

infiles = ['root://xrootd.unl.edu//store/user/molszews/EnrichedMiniAOD/DYJetsToLL_M-50_13TeV-madgraph-pythia8/DYJetsToLL_Phys14_CMSSW_7_2/150803_081314/0000/Enriched_miniAOD_1.root'];
outfile = outputdir+"test.root";

'''
for dirname, dirnames, filenames in os.walk(directory):
    for filename in filenames:
#        if '.root' in filename:
        if '10.root' in filename:
            infiles.append("file:"+directory+filename);
'''

call(["cmsRun", "maxi2ntuples.py", outfile] + infiles)



'''
#        if 'Enriched_miniAOD_16.root' in filename:
            print "running: " + filename;
            call(["cmsRun", "python/maxi2ntuples.py", directory, outputdir, filename])
'''
