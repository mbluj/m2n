from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'dy_02'
#config.General.workArea = './wjets_04/'
#config.General.workArea = './data_04/'
config.General.workArea = './crab_02/'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'maxi2ntuples.py'
config.JobType.allowUndistributedCMSSW = True

#config.Data.inputDataset = '/SingleMuon/molszews-15082015_Run2015B-PromptReco-v1_EmAOD-8e2075ef5c9b4e6dfd1c693c6f35f159/USER'
#config.Data.inputDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/molszews-20082015_RunIISpring15DR74-Asympt50ns_MCRUN2_EmAOD-f66e06a4fb4f3dd49e2c4e93f2dbebe2/USER'
#config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/molszews-16082015_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2_EmAOD-f66e06a4fb4f3dd49e2c4e93f2dbebe2/USER'
config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/tstreble-DYJets_HTauTau_05102015-2fbdb6aac8792f6fa940b34e6e638022/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 400
config.Data.totalUnits = 1200
config.Data.outLFNDirBase = '/store/user/molszews/' # or '/store/group/<subdir>'
config.Data.publication = False
config.Data.publishDataName = 'w'

#config.Site.storageSite = 'T2_PL_Warsaw'
config.Site.storageSite = 'T2_PL_Swierk'
config.Site.blacklist = ['T2_US_Florida']

