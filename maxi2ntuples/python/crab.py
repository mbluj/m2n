from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'm2n'
config.General.workArea = './susy_synch_02/'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'maxi2ntuples.py'
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = '/SUSYGluGluToHToTauTau_M-160_TuneCUETP8M1_13TeV-pythia8/molszews-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9_EmAOD-cd18adc7e66073d148bac345dd59a4ab/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 300
config.Data.totalUnits = -1
config.Data.outLFNDirBase = '/store/user/molszews/susy/' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.publishDataName = 'susy'

#config.Site.storageSite = 'T2_PL_Warsaw'
config.Site.storageSite = 'T2_PL_Swierk'
config.Site.blacklist = ['T2_US_Florida']

