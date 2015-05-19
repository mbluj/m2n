from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'm2n'
config.General.workArea = '/afs/cern.ch/work/m/molszews/CMSSW/tautau/CMSSW_7_2_3_patch1/src/m2n/maxi2ntuples/test/'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'maxi2ntuples.py'

config.Data.inputDataset = '/DYJetsToLL_M-50_13TeV-madgraph-pythia8/bluj-EnrMiniAOD_72X_synch_v2-00fe95d9ba807a4943e3f7b9ec0d4e02/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 30
config.Data.totalUnits = -1
config.Data.outLFNDirBase = '/store/user/molszews/' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.publishDataName = 'synch_test'

config.Site.storageSite = 'T2_PL_Warsaw'
config.Site.blacklist = ['T2_US_Florida']
