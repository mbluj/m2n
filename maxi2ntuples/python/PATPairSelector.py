import FWCore.ParameterSet.Config as cms

process = cms.Process("m2n")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

# Other statements
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:/opt/CMMSW/Data/Input/Phys14MiniAOD/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6.root'
        'file:/afs/cern.ch/work/m/molszews/CMSSW/Data/enrichedminiAOD.root'
    ),
#    dropDescendantsOfDroppedBranches=cms.untracked.bool(False), 
    inputCommands=cms.untracked.vstring(
      'keep *',
#      'drop *_*_*_HLT',
    )
)


process.out = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
    fileName = cms.untracked.string('/afs/cern.ch/work/m/molszews/CMSSW/Data/PATPairSelector.root'),
)


process.channelselector = cms.EDProducer("ChannelSelector",
    pairs = cms.InputTag("svfit"),
    channel = cms.string("mutau"),
)


process.pairselector = cms.EDProducer('PATPairSelector',
    leptonPair = cms.InputTag("channelselector"), 
    muCut = cms.string(''),
    elCut = cms.string(''),
    tauCut = cms.string(''),
    deltaR_ = cms.double(0.),
)


process.paircheck = cms.EDFilter("PatPairExistenceFilter",
    pairs = cms.InputTag("pairselector"),
)


process.p = cms.Path(
        process.channelselector*
        process.pairselector*
        process.paircheck
)

process.e = cms.EndPath(process.out)
