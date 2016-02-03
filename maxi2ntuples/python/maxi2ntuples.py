import FWCore.ParameterSet.Config as cms
import sys
import os

process = cms.Process("maxi2ntuples")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.load('Configuration.StandardSequences.Services_cff')                                                                                                   
process.load('JetMETCorrections.Configuration.JetCorrectionProducers_cff')
#process.load('RecoMET.METPUSubtraction.mvaPFMET_cff')

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

#####################################################################################

#####################################################################################

#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff') #MC
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff') #data

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

#process.GlobalTag.globaltag = 'PHYS14_25_V1::All'  #phys14 MC;
process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v2'
#process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v1' #25ns data
#process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v2' #25ns data
#process.GlobalTag.globaltag = '74X_dataRun2_reMiniAOD_v0' #25ns data october
#process.GlobalTag.globaltag = 'MCRUN2_74_V9A'  #spring15 50ns MC;
#process.GlobalTag.globaltag = 'MCRUN2_74_V9'  #spring15 25ns MC;
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

mc=True; #if MC then true; if data then  false
sample = 1; #0 -data; 1-DY; 2-WJets; 3-TTbar; 4-QCD
outfile = "DYv3.root";
vbf=False
grid=False
aod = True
minioadv2 = True
#####################################################################################

#####################################################################################


#Directory with input file(s). Do not put ".root" files there that are not maent to be processed.
directory = '/afs/cern.ch/work/m/molszews/CMSSW/Data/EmAOD/'
files = ['DYJetsToLLv3.root'];


#Directory with outputfile(s). Can be of course the same as the above one, but remember to remove ouput files before another run.
outputdir = '/afs/cern.ch/work/m/molszews/CMSSW/Data/ntuple_VBF/' 

def getfiles(directory, files = []):
    infiles =[];
    for dirname, dirnames, filenames in os.walk(directory):
        if not files:
            for filename in filenames:
                if '.root' in filename:
                    infiles.append("file:"+directory+filename);
        else:
            for filename in filenames:
                if any(filename in s for s in files):
                    infiles.append("file:"+directory+filename);
    print infiles;
    return infiles;           



process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )
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
#myfilelist = cms.untracked.vstring()
#myfilelist.extend(sys.argv[3:]);

#print myfilelist;
#Number_of_events = 464755;
#xSection = 3.748;

#nadpisanie(dir+inputFile+'.json', "Number_of_events = "+ str(Number_of_events)+'\n')
#pisanie(dir+inputFile+'.json', "xSection = "+ str(xSection)+'\n')
if grid:
    process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
else:
    process.source = cms.Source("PoolSource",
        # replace 'myfile.root' with the source file you want to use
        fileNames = cms.untracked.vstring(
            getfiles(directory, files)        
    #        'file:/afs/cern.ch/work/m/molszews/CMSSW/Data/mbluj/Enriched_miniAOD_100_1_qrj.root'
        )
    )

'''
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('miniAODMVAMET.root'),
)
'''
if grid:
    process.TFileService = cms.Service("TFileService", fileName = cms.string(outfile))
else:
    process.TFileService = cms.Service("TFileService", fileName = cms.string(outputdir+outfile))

######################################################################################################
process.eventnumberfilter = cms.EDFilter("EventNumberFilter",
)

process.ininfo = cms.EDAnalyzer("ininfo",
    mc = cms.bool(mc),
#    pairs = cms.InputTag("SVllCand"),
    pairs = cms.InputTag("SVbypass"),
)

############### ELECTRON MVA ID ###################
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff']
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)



############### JETS ##############################
process.jetsSelected = cms.EDFilter("PATJetSelector",
        src = cms.InputTag("slimmedJets"),
        cut = cms.string('abs(eta) < 4.7  & pt > 20.'
           #     'et > 30.'
            ),
        filter = cms.bool(False)
        )

process.jetsIDSelected = cms.EDProducer("JetsSelector",
    jets = cms.InputTag("jetsSelected"),
        )

############### PAIRS #############################


##
## Build ll candidates (here OS)
##
process.pairmaker = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("slimmedMuons slimmedTaus"),
    cut = cms.string("mass>0"),
    checkCharge = cms.bool(False)
)


###################################################
MVAPairMET = ();
for index in range(100):
  MVAMETName = "pfMETMVA%i" % index
  MVAPairMET += (cms.InputTag(MVAMETName),)


process.pairswithmet = cms.EDProducer("AddMVAMET",
#    pairs = cms.InputTag("SVllCand"),
#    pairs = cms.InputTag("SVbypass"),
    mets = cms.InputTag("slimmedMETs"),
    pairsmets = cms.VInputTag(MVAPairMET),
    mvamet = cms.InputTag("pfMETMVA0"),
    useMVAMET  = cms.untracked.bool(False),
    usePairMET = cms.untracked.bool(False),
)
if aod:
    process.pairswithmet.pairs =  cms.InputTag("pairmaker") #reco::CompositeCandidateCollection
else:
    process.pairswithmet.pairs =  cms.InputTag("SVllCand") #pat::CompositeCandidateCollection
    #pairs = cms.InputTag("SVbypass"),

process.channel = cms.EDProducer("ChannelSelector",
    pairs = cms.InputTag("pairswithmet"),
    channel = cms.string("mutau"), #mutau, etau, tautau...
)

process.pairchecka = cms.EDFilter("PatPairExistenceFilter",
    pairs = cms.InputTag("channel"),
)

process.clean = cms.EDProducer('PATPairSelector',
    pairs = cms.InputTag("channel"), 
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

process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")

#BASELINE SELECTION
process.selected = cms.EDProducer("PairBaselineSelection",
    pairs = cms.InputTag("clean"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    mc = cms.bool(mc),
    eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
    sample = cms.string('spring15')  #options: "phys14", "spring15"
)

process.hlt= cms.EDProducer("HLTforPair",
    pairs = cms.InputTag("selected"),
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("selectedPatTrigger"),
)



#VETO SELECTION
process.vetoed = cms.EDProducer("PostSynchSelection",
    pairs = cms.InputTag("hlt"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("selectedPatTrigger"),
    mc = cms.bool(mc),
    eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
    sample = cms.string('spring15')  #options: "phys14", "spring15"
)



process.paircheckb = cms.EDFilter("PatPairExistenceFilter",
#    pairs = cms.InputTag("selected"),
    pairs = cms.InputTag("vetoed"),
)

process.eventskimmer = cms.EDProducer("EventsSkimmer",
    mc = cms.bool(mc),
    pairs = cms.InputTag("vetoed"),
)



process.bestpair = cms.EDProducer("BestPairSelector",
    pairs = cms.InputTag("eventskimmer"),
)

process.load('L1Trigger.Skimmer.l1Filter_cfi')
process.l1Filter.algorithms = cms.vstring('L1_Mu16er_TauJet20er', 'L1_SingleMu20er')

process.load('HLTrigger.HLTfilters.hltLevel1GTSeed_cfi')
process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(False)
process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('(L1_SingleMuBeamHalo OR L1_SingleMuOpen)  AND NOT L1_SingleJet6U')

#process.m2n = cms.EDAnalyzer('maxi2ntuples',
#process.m2n = cms.EDAnalyzer('tautau',
process.m2n = cms.EDAnalyzer('ntuple',

    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    taus = cms.InputTag("slimmedTaus"),
    photons = cms.InputTag("slimmedPhotons"),
    jets = cms.InputTag("jetsIDSelected"),
    fatjets = cms.InputTag("slimmedJetsAK8"),
    mets = cms.InputTag("slimmedMETs"),
    pairs = cms.InputTag("bestpair"),
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("selectedPatTrigger"),
    prunedGenParticles = cms.InputTag("prunedGenParticles"),
    packedGenParticles = cms.InputTag("packedGenParticles"),
    pileupinfo = cms.InputTag("addPileupInfo"),
    vertexScores = cms.InputTag("offlineSlimmedPrimaryVertices"),
    src = cms.InputTag("packedPFCandidates"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    mc = cms.bool(mc),
    sample = cms.int32(sample)
)

process.synchtree = cms.EDAnalyzer('synchronization',

    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    taus = cms.InputTag("slimmedTaus"),
    photons = cms.InputTag("slimmedPhotons"),
    jets = cms.InputTag("jetsIDSelected"),
    fatjets = cms.InputTag("slimmedJetsAK8"),
    mets = cms.InputTag("slimmedMETs"),
    pairs = cms.InputTag("bestpair"),
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("selectedPatTrigger"),
    prunedGenParticles = cms.InputTag("prunedGenParticles"),
    packedGenParticles = cms.InputTag("packedGenParticles"),
    mc = cms.bool(mc),
)

if minioadv2:
    process.m2n.pileupinfo = cms.InputTag("slimmedAddPileupInfo")
else:
    process.m2n.pileupinfo = cms.InputTag("addPileupInfo")

if mc and not vbf:
    process.m2n.lheprod =  cms.InputTag("externalLHEProducer");
    process.synchtree.lheprod =  cms.InputTag("externalLHEProducer");
else:
    process.m2n.lheprod =  cms.InputTag("source");
    process.synchtree.lheprod =  cms.InputTag("source");

process.maxi2ntuples = cms.EDAnalyzer('maxi2ntuples',

    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    taus = cms.InputTag("slimmedTaus"),
    photons = cms.InputTag("slimmedPhotons"),
    jets = cms.InputTag("jetsSelected"),
    fatjets = cms.InputTag("slimmedJetsAK8"),
    mets = cms.InputTag("slimmedMETs"),
    pairs = cms.InputTag("channel"),
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("selectedPatTrigger"),
)

process.p = cms.Path(
        process.ininfo
        *process.pairmaker
        *process.egmGsfElectronIDSequence
        *process.jetsSelected
        *process.jetsIDSelected
        *process.pairswithmet
        *process.channel
#        *process.pairchecka
#        *process.hltLevel1GTSeed
#        *process.l1Filter
        *process.clean
        *process.electronMVAValueMapProducer
        *process.selected
        *process.hlt
        *process.vetoed
        *process.eventskimmer
#        *process.paircheckb
        *process.bestpair
        *process.m2n
#        *process.synchtree
)
process.o = cms.Path(
        process.ininfo
        *process.egmGsfElectronIDSequence
        *process.jetsSelected
        *process.jetsIDSelected
        *process.pairswithmet
        *process.channel
#        *process.pairchecka
#        *process.hltLevel1GTSeed
#        *process.l1Filter
        *process.clean
        *process.electronMVAValueMapProducer
        *process.selected
        *process.hlt
        *process.vetoed
        *process.eventskimmer
#        *process.paircheckb
        *process.bestpair
        *process.m2n
#        *process.synchtree
)
'''
process.p = cms.Path(
        process.eventnumberfilter
        *process.pairswithmet
        *process.pairchecka
        *process.channel
        *process.maxi2ntuples
)
'''
#process.e = cms.EndPath(process.out)
#process.o = cms.Path(process.pairmaker)
if aod:
    process.schedule = cms.Schedule(process.p)
else:
    process.schedule = cms.Schedule(process.o)

#process.schedule = cms.Schedule(process.p)

