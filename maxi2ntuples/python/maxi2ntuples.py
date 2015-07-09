import FWCore.ParameterSet.Config as cms
import sys

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

#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag.globaltag = 'PHYS14_25_V1::All'
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#################################### FILES #####################################################
#dopisuje na koniec pliku
def pisanie(plik, tekst):
    with open(plik,"a") as f:
        f.write(tekst)

def nadpisanie(plik, tekst):
    f = open(plik,"w")
    f.write(tekst)
    f.close()


#dir = '/afs/cern.ch/work/m/molszews/CMSSW/Data/EmAOD_VBF/'
#inputFile = "Enriched_miniAOD_100";
#outputFile = "ntuples.root";
dir = sys.argv[2]
odir = sys.argv[3]
inputFile = sys.argv[4];

#Number_of_events = 464755;
#xSection = 3.748;

#nadpisanie(dir+inputFile+'.json', "Number_of_events = "+ str(Number_of_events)+'\n')
#pisanie(dir+inputFile+'.json', "xSection = "+ str(xSection)+'\n')

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:'+dir+inputFile
#        'file:/afs/cern.ch/work/m/molszews/CMSSW/Data/mbluj/Enriched_miniAOD_100_1_qrj.root'
    )
)

'''
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('miniAODMVAMET.root'),
)
'''

process.TFileService = cms.Service("TFileService", fileName = cms.string(odir+'ntuple_'+inputFile))
######################################################################################################

############### JETS ##############################

process.jetsSelected = cms.EDFilter("PATJetSelector",
        src = cms.InputTag("slimmedJets"),
        cut = cms.string('abs(eta) < 4.7  & pt > 20.'
           #     'et > 30.'
            ),
        filter = cms.bool(False)
        )

############### PAIRS #############################

MVAPairMET = ();
for index in range(100):
  MVAMETName = "pfMETMVA%i" % index
  MVAPairMET += (cms.InputTag(MVAMETName),)


process.pairswithmet = cms.EDProducer("AddMVAMET",
    pairs = cms.InputTag("SVllCand"),
    mets = cms.InputTag("slimmedMETs"),
    pairsmets = cms.VInputTag(MVAPairMET),
    mvamet = cms.InputTag("pfMETMVA0"),
    useMVAMET  = cms.untracked.bool(False),
    usePairMET = cms.untracked.bool(False),
)

process.pairs = cms.EDProducer("ChannelSelector",
    pairs = cms.InputTag("pairswithmet"),
    channel = cms.string("tautau"),
)


process.clean = cms.EDProducer('PATPairSelector',
    pairs = cms.InputTag("pairs"), 
#    pairs = cms.InputTag('pairs'), 
    muCut = cms.string(''
#        'pt > 18. & abs(eta) < 2.1'
        ),
    elCut = cms.string(''),
    tauCut = cms.string(''
#        "tauID('byCombinedIsolationDeltaBetaCorrRaw3Hits') < 1.5 & "
#        "tauID('againstMuonTight3')"
        ),
    deltaR_ = cms.double(0.5),
)

process.selected = cms.EDProducer("PairBaselineSelection",
    pairs = cms.InputTag("clean"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
)

process.paircheck = cms.EDFilter("PatPairExistenceFilter",
    pairs = cms.InputTag("selected"),
)


#process.m2n = cms.EDAnalyzer('maxi2ntuples',
process.m2n = cms.EDAnalyzer('tautau',

    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    taus = cms.InputTag("slimmedTaus"),
    photons = cms.InputTag("slimmedPhotons"),
    jets = cms.InputTag("jetsSelected"),
    fatjets = cms.InputTag("slimmedJetsAK8"),
    mets = cms.InputTag("slimmedMETs"),
    pairs = cms.InputTag("selected"),
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("selectedPatTrigger"),
    prunedGenParticles = cms.InputTag("prunedGenParticles"),
    packedGenParticles = cms.InputTag("packedGenParticles"),
    lheprod = cms.InputTag("source"),
#    lheprod = cms.InputTag("externalLHEProducer"),
    pileupinfo = cms.InputTag("addPileupInfo"),
)

process.p = cms.Path(
        process.jetsSelected
        *process.pairswithmet
        *process.pairs
        *process.clean
        *process.selected
        *process.paircheck
        *process.m2n
)
 
#process.e = cms.EndPath(process.out)
