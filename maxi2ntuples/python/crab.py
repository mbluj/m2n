from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'm2n'
config.General.workArea = './vbf_synch_o1/'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'maxi2ntuples.py'
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = '/VBF_HToTauTau_M-125_13TeV-powheg-pythia6/molszews-Phys14DR-PU20bx25_tsg_PHYS14_EmAOD-258495b13287afbd0ec47b699869be96/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 300
config.Data.totalUnits = -1
config.Data.outLFNDirBase = '/store/user/molszews/' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.publishDataName = 'vbf_synch_o1'

#config.Site.storageSite = 'T2_PL_Warsaw'
config.Site.storageSite = 'T2_PL_Swierk'
config.Site.blacklist = ['T2_US_Florida']

