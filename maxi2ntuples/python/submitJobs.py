#!/usr/bin/env python

import os, re
import commands
import math

#from crab3 import *
from CRABClient.UserUtilities import config 
config = config()
#########################################
#########################################
def prepareCrabCfg(dataset,
                   crabCfgName,
                   eventsPerJob,
                   jsonFile,
                   storage_element, 
                   publish_data_suffix,
                   splitmethod):

    workdir = publish_data_suffix
    shortName = dataset.split("/")[1]
    if dataset.split("/")[2].find("Run2015")!=-1:
        shortName += "_"+dataset.split("/")[2]
    shortName = shortName.replace("-","_")
    shortName+="_"+publish_data_suffix
    ##Modify CRAB3 configuration
    config.General.requestName = shortName
    config.General.workArea = '/afs/cern.ch/work/m/molszews/CMSSW/Data/m2n/03'
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = 'maxi2ntuples.py'
    config.JobType.allowUndistributedCMSSW = True
    config.Data.inputDBS = 'global'
    config.Data.splitting = splitmethod
    config.Data.unitsPerJob = eventsPerJob
    config.Data.totalUnits = -1
    config.Data.inputDataset = dataset
    config.Data.outLFNDirBase = '/store/user/molszews/m2n/'+shortName
    config.Data.publication = False
    config.Data.outputDatasetTag = shortName
    config.Site.storageSite = storage_element
    if len(jsonFile):
        config.Data.lumiMask=jsonFile
    out = open('crabTmp.py','w')
    out.write(config.pythonise_())
    out.close()
    os.system("crab submit -c crabTmp.py")
#########################################
#########################################
units = 25

datasetsMC = [
    ##'/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/tstreble-WJets_HTauTau_05102015-2fbdb6aac8792f6fa940b34e6e638022/USER' #LLR WJets campaig1.0
    ##'/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/tstreble-DYJets_HTauTau_05102015-2fbdb6aac8792f6fa940b34e6e638022/USER' #LLR WY campaig1.0
    ##'/SingleMuon/molszews-15082015_Run2015B-PromptReco-v1_EmAOD-8e2075ef5c9b4e6dfd1c693c6f35f159/USER'
    ##'/SingleMuon/molszews-SingleMuon_Run2015C_PromptReco_v1_v1_25ns-42a7e4340eeaef6d4f87c5c983a49b31/USER'
    ##'/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/molszews-20082015_RunIISpring15DR74-Asympt50ns_MCRUN2_EmAOD-f66e06a4fb4f3dd49e2c4e93f2dbebe2/USER'
    ##'/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/molszews-16082015_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2_EmAOD-f66e06a4fb4f3dd49e2c4e93f2dbebe2/USER'
    ##'/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/molszews-TTJets_TuneCUETP8M1_13TeV_amcatnloFXFX_pythia8_v1-de0a0f5c6787d3ec8a72b569e41467a2/USER'#TTbezSVfita
    ## minioadV2 ----------------------
    ##'/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
    ##'/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
    ##'/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v3/MINIAODSIM'
    ##'/DYJetsToLL_M-5to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
    '/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'    
]

########################################################
for dataset in datasetsMC:
    prepareCrabCfg(crabCfgName="crab3.py",
                   dataset=dataset,
                   eventsPerJob=units,
                   jsonFile="",
                   storage_element="T2_PL_Swierk",
                   publish_data_suffix = "v3",
                   splitmethod = 'FileBased')
########################################################


datasetsDATA = [
    ##"/MuonEG/Run2015B-05Aug2015-v1/MINIAOD",
    ##"/MuonEG/Run2015B-17Jul2015-v1/MINIAOD",
    ##"/MuonEG/Run2015B-PromptReco-v1/MINIAOD",
    ##"/MuonEG/Run2015C-PromptReco-v1/MINIAOD"
    ##'/SingleMuon/Run2015D-05Oct2015-v1/MINIAOD'
]

jsonFile50ns = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-255031_13TeV_PromptReco_Collisions15_50ns_JSON_v2.txt"
#jsonFile25ns = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-255031_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt"
jsonFile25ns = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260426_13TeV_PromptReco_Collisions15_25ns_JSON.txt"

jsonFile = jsonFile25ns
########################################################
for dataset in datasetsDATA:
    prepareCrabCfg(crabCfgName="crab3.py",
                   dataset=dataset,
                   eventsPerJob=units,
                   jsonFile=jsonFile,
                   storage_element="T2_PL_Swierk",
                   publish_data_suffix = "v1_25ns",
                   splitmethod = 'FileBased')
########################################################

