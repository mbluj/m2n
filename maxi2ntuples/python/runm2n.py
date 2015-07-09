import os
from subprocess import call

#Directory with input file(s). Do not put ".root" files there that are not maent to be processed.
directory = '/afs/cern.ch/work/m/molszews/CMSSW/Data/EmAOD_VBF/'
#Directory with outputfile(s). Can be of course the same as the above one, but remember to remove ouput files before another run.
outputdir = '/afs/cern.ch/work/m/molszews/CMSSW/Data/ntuple_VBF/' 

for dirname, dirnames, filenames in os.walk(directory):
    for filename in filenames:
        if '.root' in filename:
            print "running: " + filename;
            call(["cmsRun", "python/maxi2ntuples.py", directory, outputdir, filename])

