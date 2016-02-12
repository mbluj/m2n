//
// Package:    m2n/maxi2ntuples
// Class:      ntuple
// 
/**\class ntuple ntuple.cc m2n/maxi2ntuples/plugins/ntuple.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Wed, 17 Dec 2014 09:06:30 GMT
//
//


// system include files
#include <memory>
#include <vector>
#include <string>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include <DataFormats/Candidate/interface/ShallowCloneCandidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"


#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include <DataFormats/METReco/interface/CommonMETData.h>
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalTrajectoryExtrapolatorToLine.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"

#include "m2n/maxi2ntuples/interface/utilities.h"
#include "m2n/HTTDataFormats/interface/HTTEvent.h"

#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TMVA/Reader.h"
#include "TLorentzVector.h"
#include "RConfig.h" 
#include "TObject.h"
#include "TH1F.h"

class ntuple : public edm::EDAnalyzer {
 public:
  explicit ntuple(const edm::ParameterSet&);
  ~ntuple();
  
  //static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  ///Find PV using different discriminators. Return false
  ///if no vertices were found in the event
  ///pfPV  - using PF particles for score calulation
  ///pt2PV - using sum pt^2 os all tracks assigned to vertex.
  ///        as tracks pt<1 do not have errors, we take
  ///        estimate based on TRK-11-001 note    
  bool findPrimaryVertices(const edm::Event & iEvent, const edm::EventSetup & iSetup);

  //Refit PV using track information stored in miniAOD
  ///This has to be run AFTER finding tau candidates. WARING:
  ///first verion will not exclude tau tracks from fit.
  bool refitPV(const edm::Event & iEvent, const edm::EventSetup & iSetup,
	       const TLorentzVector & leadingTkTau1,
	       const TLorentzVector & leadingTkTau2);

  ///Calculate vector from PV to Point of Closest Approach (PCA).
  ///Using leading track extracted from pat::Tau,
  ///and PV position passed as aPoint
  TVector3 getPCA(const edm::Event & iEvent, const edm::EventSetup & iSetup,
		  const reco::Track *aTrack, const GlobalPoint & aPoint);

  ///Set PCA vector for reco taus. This must be run AFTER
  ///vertex refitting. 
  ///WARNING: The T object is modified in this method.
  template<typename T> void setPCAVectorsOnObject(T & aObject, 			   
						  const reco::Track *aTrack,
						  const edm::Event & iEvent, const edm::EventSetup & iSetup);
    
 private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob(){};

  void clearNtuple();
  void getCollections(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void getMCCollections(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  void fillEventHeaderData(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillBosonDecayMode();
  void fillGenTausAndPV();
  void fillGenData();
  void fillTauPairData(const edm::Event&, const edm::EventSetup&);
  void fillMuonLeg(const reco::Candidate *, const reco::Candidate *, const edm::Event&, const edm::EventSetup&);
  void fillElectronLeg(const reco::Candidate *, const reco::Candidate *, const edm::Event&, const edm::EventSetup&);
  void fillTauLeg(const reco::Candidate *, const reco::Candidate *, const edm::Event&, const edm::EventSetup&);
  void fillMET(const reco::Candidate *, const pat::CompositeCandidate &);
  void fillJetsData();

  void fillGenTauData(const reco::GenParticleRef &, const reco::GenParticleRefVector &);

  typedef std::vector<PileupSummaryInfo> PileupSummaryInfoCollection;
  
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken;
  edm::EDGetTokenT<edm::ValueMap<float> > scoresToken;  
  edm::EDGetTokenT<pat::MuonCollection> muonToken;
  edm::EDGetTokenT<pat::ElectronCollection> electronToken;
  edm::EDGetTokenT<pat::TauCollection> tauToken;
  edm::EDGetTokenT<pat::PhotonCollection> photonToken;
  edm::EDGetTokenT<pat::JetCollection> jetToken;
  edm::EDGetTokenT<pat::JetCollection> fatjetToken;
  edm::EDGetTokenT<pat::METCollection> metToken;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> pairToken;
  edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsToken;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesToken;
  edm::EDGetTokenT<reco::GenParticleCollection> prunedGenToken;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGenToken;
  edm::EDGetTokenT<LHEEventProduct> lheprodToken;
  edm::EDGetTokenT<PileupSummaryInfoCollection> PileupSummaryInfoToken;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate> > candsToken;
  
  edm::Handle<PileupSummaryInfoCollection> genPileUpInfos;
  edm::Handle<LHEEventProduct> lheInfo;
  edm::Handle<GenEventInfoProduct> genInfo;
  edm::Handle<reco::GenParticleCollection> genParticlesPruned;
  edm::Handle<pat::PackedGenParticleCollection> genParticlesPacked;

  edm::Handle<reco::BeamSpot> beamSpot;
  edm::Handle<reco::VertexCollection> vertices;
  edm::Handle<edm::ValueMap<float> > scores;
  edm::Handle<edm::View<pat::PackedCandidate> >  cands;
  edm::Handle<pat::JetCollection> jets; 
  edm::Handle<pat::METCollection> mets;
  edm::Handle<pat::CompositeCandidateCollection> pairs;
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

  const int sample;

  ///Helper variables used for VX refitting.
  unsigned int pfPVIndex;

    TTree * eventTree; 
    TH1F* events;

    Wevent *wevent;
    WtauCollection wtaucollection;
    WtauCollection wtauGencollection;
    WmuCollection wmucollection;
    WelectronCollection welectroncollection;
    WpairCollection wpaircollection;
    WmetCollection wmetcollection;
    WjetCollection wjetcollection;


};




