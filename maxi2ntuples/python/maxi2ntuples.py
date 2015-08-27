import FWCore.ParameterSet.Config as cms
import sys
import os

process = cms.Process("maxi2ntuples")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load('Configuration.StandardSequences.Services_cff')                                                                                                   
process.load('JetMETCorrections.Configuration.JetCorrectionProducers_cff')
process.load('RecoMET.METPUSubtraction.mvaPFMET_cff')

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')


mc=True; #if MC then true; if data then  false
vbf=True
grid=True;

#Directory with input file(s). Do not put ".root" files there that are not maent to be processed.
directory = '/afs/cern.ch/work/m/molszews/CMSSW/Data/EmAOD_VBF/'
files = ['VBF.root'];


#Directory with outputfile(s). Can be of course the same as the above one, but remember to remove ouput files before another run.
outputdir = '/afs/cern.ch/work/m/molszews/CMSSW/Data/ntuple_VBF/' 
outfile = files[0];

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


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag.globaltag = 'PHYS14_25_V1::All'  #phys14 MC;
#process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v0' #50ns data
#process.GlobalTag.globaltag = 'MCRUN2_74_V9A'  #spring15 50ns MC;
#process.GlobalTag.globaltag = 'MCRUN2_74_V9'  #spring15 25ns MC;
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
# overwrite a default parameter: for miniAOD, the collection name is a slimmed one
process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons')
from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)


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


process.ininfo = cms.EDAnalyzer("ininfo",
    mc = cms.bool(mc),
    pairs = cms.InputTag("SVllCand"),
)


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
#    pairs = cms.InputTag("SVbypass"),
    mets = cms.InputTag("slimmedMETs"),
    pairsmets = cms.VInputTag(MVAPairMET),
    mvamet = cms.InputTag("pfMETMVA0"),
    useMVAMET  = cms.untracked.bool(False),
    usePairMET = cms.untracked.bool(False),
)

process.pairs = cms.EDProducer("ChannelSelector",
    pairs = cms.InputTag("pairswithmet"),
    channel = cms.string("mutau"),
)


process.clean = cms.EDProducer('PATPairSelector',
    pairs = cms.InputTag("pairs"), 
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
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("selectedPatTrigger"),
    mc = cms.bool(mc),
    sample = cms.string('phys14')  #options: "phys14", "spring15"
)

process.bestpair = cms.EDProducer("BestPairSelector",
    pairs = cms.InputTag("selected"),
    
)

process.paircheck = cms.EDFilter("PatPairExistenceFilter",
    pairs = cms.InputTag("selected"),
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
    jets = cms.InputTag("jetsSelected"),
    fatjets = cms.InputTag("slimmedJetsAK8"),
    mets = cms.InputTag("slimmedMETs"),
    pairs = cms.InputTag("bestpair"),
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("selectedPatTrigger"),
    prunedGenParticles = cms.InputTag("prunedGenParticles"),
    packedGenParticles = cms.InputTag("packedGenParticles"),
    pileupinfo = cms.InputTag("addPileupInfo"),
    mc = cms.bool(mc),
)

process.synchtree = cms.EDAnalyzer('synchronization',

    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    taus = cms.InputTag("slimmedTaus"),
    photons = cms.InputTag("slimmedPhotons"),
    jets = cms.InputTag("jetsSelected"),
    fatjets = cms.InputTag("slimmedJetsAK8"),
    mets = cms.InputTag("slimmedMETs"),
    pairs = cms.InputTag("bestpair"),
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("selectedPatTrigger"),
    prunedGenParticles = cms.InputTag("prunedGenParticles"),
    packedGenParticles = cms.InputTag("packedGenParticles"),
    pileupinfo = cms.InputTag("addPileupInfo"),
    mc = cms.bool(mc),
)

if mc and not vbf:
    process.m2n.lheprod =  cms.InputTag("externalLHEProducer");
    process.synchtree.lheprod =  cms.InputTag("externalLHEProducer");
else:
    process.m2n.lheprod =  cms.InputTag("source");
    process.synchtree.lheprod =  cms.InputTag("source");

process.p = cms.Path(
        process.ininfo
        *process.jetsSelected
        *process.pairswithmet
        *process.pairs
#        *process.hltLevel1GTSeed
#        *process.l1Filter
        *process.clean
        *process.selected
        *process.paircheck
        *process.bestpair
        *process.m2n
        *process.synchtree
)
 
#process.e = cms.EndPath(process.out)
