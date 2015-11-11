from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.workArea = '/afs/cern.ch/work/m/molszews/CMSSW/Data/m2n/01'
config.General.requestName = 'TTJets_TuneCUETP8M1_13TeV_amcatnloFXFX_pythia8_v1'
config.section_('JobType')
config.JobType.psetName = 'maxi2ntuples.py'
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.section_('Data')
config.Data.inputDataset = '/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/molszews-TTJets_TuneCUETP8M1_13TeV_amcatnloFXFX_pythia8_v1-de0a0f5c6787d3ec8a72b569e41467a2/USER'
config.Data.outputDatasetTag = 'TTJets_TuneCUETP8M1_13TeV_amcatnloFXFX_pythia8_v1'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 100
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/user/molszews/m2n/TTJets_TuneCUETP8M1_13TeV_amcatnloFXFX_pythia8_v1'
config.Data.publication = False
config.section_('Site')
config.Site.storageSite = 'T2_PL_Swierk'
config.section_('User')
config.section_('Debug')
