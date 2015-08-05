from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'm2n'
config.General.workArea = './test/'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'maxi2ntuples.py'
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = '/DYJetsToLL_M-50_13TeV-madgraph-pythia8/molszews-DYJetsToLL_Phys14_CMSSW_7_2-14021cf92b609bd2cdfd7d04edfc7ab0/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 30
config.Data.totalUnits = -1
config.Data.outLFNDirBase = '/store/user/molszews/' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.publishDataName = 'synch_test'

config.Site.storageSite = 'T2_PL_Warsaw'
config.Site.blacklist = ['T2_US_Florida']

