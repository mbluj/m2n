from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.workArea = '/afs/cern.ch/work/m/molszews/CMSSW/Data/m2n/03'
config.General.requestName = 'DYJetsToLL_M_5to50_TuneCUETP8M1_13TeV_madgraphMLM_pythia8_v3'
config.section_('JobType')
config.JobType.psetName = 'maxi2ntuples.py'
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.section_('Data')
config.Data.inputDataset = '/DYJetsToLL_M-5to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
config.Data.outputDatasetTag = 'DYJetsToLL_M_5to50_TuneCUETP8M1_13TeV_madgraphMLM_pythia8_v3'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 25
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/user/molszews/m2n/DYJetsToLL_M_5to50_TuneCUETP8M1_13TeV_madgraphMLM_pythia8_v3'
config.Data.publication = False
config.section_('Site')
config.Site.storageSite = 'T2_PL_Swierk'
config.section_('User')
config.section_('Debug')
