import os
from subprocess import call

for dirname, dirnames, filenames in os.walk('/opt/CMMSW/Data/m2m/WZJetsTo3LNu/'):
    for filename in filenames:
        if '.root' in filename:
            print "running: " + filename;
            call(["cmsRun", "maxi2ntuples.py", filename])

