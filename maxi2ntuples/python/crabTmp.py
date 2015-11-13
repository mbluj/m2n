from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.workArea = '/afs/cern.ch/work/m/molszews/CMSSW/Data/m2n/03'
config.General.requestName = 'SingleMuon_Run2015D_PromptReco_v3_v1_25ns'
config.section_('JobType')
config.JobType.psetName = 'maxi2ntuples.py'
config.JobType.pluginName = 'Analysis'
config.JobType.allowUndistributedCMSSW = True
config.section_('Data')
config.Data.inputDataset = '/SingleMuon/Run2015D-PromptReco-v3/MINIAOD'
config.Data.outputDatasetTag = 'SingleMuon_Run2015D_PromptReco_v3_v1_25ns'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 100
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/user/molszews/m2n/SingleMuon_Run2015D_PromptReco_v3_v1_25ns'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-258159_13TeV_PromptReco_Collisions15_25ns_JSON_v3.txt'
config.Data.publication = False
config.section_('Site')
config.Site.storageSite = 'T2_PL_Swierk'
config.section_('User')
config.section_('Debug')
