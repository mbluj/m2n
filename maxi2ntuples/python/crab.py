from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'm2n'
config.General.workArea = './WJetsToLNu_01/'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'maxi2ntuples.py'
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/molszews-20082015_RunIISpring15DR74-Asympt50ns_MCRUN2_EmAOD-f66e06a4fb4f3dd49e2c4e93f2dbebe2/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 205
config.Data.totalUnits = -1
config.Data.outLFNDirBase = '/store/user/molszews/' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.publishDataName = 'w'

#config.Site.storageSite = 'T2_PL_Warsaw'
config.Site.storageSite = 'T2_PL_Swierk'
config.Site.blacklist = ['T2_US_Florida']

