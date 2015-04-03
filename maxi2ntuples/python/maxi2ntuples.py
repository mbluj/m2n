import FWCore.ParameterSet.Config as cms
#import sys

process = cms.Process("maxi2ntuples")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load('Configuration.StandardSequences.Services_cff')                                                                                                   
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('JetMETCorrections.Configuration.JetCorrectionProducers_cff')
process.load('RecoMET.METPUSubtraction.mvaPFMET_cff')

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


#inoutfilename = sys.argv[2];
#direc = 'WZJetsTo3LNu/'
#print inoutfilename;

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/opt/CMMSW/Data/enritchedGluGluToHToTauTau.root'
    )
)

outputFile = "/opt/CMMSW/Data/ntuples.root";
process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile))


############### JETS ##############################

process.jetsSelected = cms.EDFilter("PATJetSelector",
        src = cms.InputTag("slimmedJets"),
        cut = cms.string('abs(eta) < 4.7  & pt > 20.'
           #     'et > 30.'
            ),
        filter = cms.bool(False)
        )


############### PAIRS #############################
process.mutauPairs = cms.EDProducer("ChannelSelector",
    pairs = cms.InputTag("pairs"),
    channel = cms.string("mutau"),
)

process.baselineselected = cms.EDProducer("PairBaselineSelection",
    pairs = cms.InputTag("mutauPairs"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
)


process.mutauClean = cms.EDProducer('PATPairSelector',
    pairs = cms.InputTag("baselineselected"), 
    muCut = cms.string(''
#        'pt > 18. & abs(eta) < 2.1'
        ),
    elCut = cms.string(''),
    tauCut = cms.string(''
#        "tauID('byCombinedIsolationDeltaBetaCorrRaw3Hits') < 1.5 & "
#        "tauID('againstMuonTight3')"
        ),
    deltaR_ = cms.double(0.),
)

process.paircheck = cms.EDFilter("PatPairExistenceFilter",
    pairs = cms.InputTag("mutauClean"),
)


process.m2n = cms.EDAnalyzer('maxi2ntuples',

    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    taus = cms.InputTag("slimmedTaus"),
    photons = cms.InputTag("slimmedPhotons"),
    jets = cms.InputTag("jetsSelected"),
    fatjets = cms.InputTag("slimmedJetsAK8"),
    mets = cms.InputTag("slimmedMETs"),
    pairs = cms.InputTag("mutauClean"),
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("selectedPatTrigger"),

)

process.p = cms.Path(
        process.jetsSelected* 
        process.mutauPairs* process.baselineselected* process.mutauClean* 
        process.paircheck*
        process.m2n
)
  
