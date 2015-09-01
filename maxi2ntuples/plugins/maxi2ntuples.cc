//
// Package:    m2n/maxi2ntuples
// Class:      maxi2ntuples
// 
/**\class maxi2ntuples maxi2ntuples.cc m2n/maxi2ntuples/plugins/maxi2ntuples.cc

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

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include <DataFormats/METReco/interface/CommonMETData.h>
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "m2n/maxi2ntuples/interface/utilities.h"

#include "TTree.h"
#include "TFile.h"
#include "TMVA/Reader.h"
#include "TLorentzVector.h"
#include "RConfig.h" 
#include "TObject.h"

//#include "clasadict.h"

//
// class declaration
//

//gSystem->Load("libPhysics.so");
class maxi2ntuples : public edm::EDAnalyzer {
   public:
      explicit maxi2ntuples(const edm::ParameterSet&);
      ~maxi2ntuples();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      typedef std::vector<PileupSummaryInfo> PileupSummaryInfoCollection;


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      void fillmuon(const pat::Muon*, const reco::Vertex&, const reco::Candidate*);
      void filltau(const pat::Tau*, const reco::Candidate*);

  // ----------member data ---------------------------
    edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
    edm::EDGetTokenT<pat::MuonCollection> muonToken_;
    edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
    edm::EDGetTokenT<pat::TauCollection> tauToken_;
    edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
    edm::EDGetTokenT<pat::JetCollection> jetToken_;
    edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
    edm::EDGetTokenT<pat::METCollection> metToken_;
    edm::EDGetTokenT<pat::CompositeCandidateCollection> pairToken_;
    edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
//    edm::EDGetTokenT<reco::GenParticleCollection> prunedGenToken_;
//    edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGenToken_;
//    edm::EDGetTokenT<LHEEventProduct> lheprodToken_;
//    edm::EDGetTokenT<PileupSummaryInfoCollection> PileupSummaryInfoToken_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//


// -*- C++ -*-
maxi2ntuples::maxi2ntuples(const edm::ParameterSet& iConfig):
    vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
    electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
    tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
    photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
    jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
    fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
    metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
    pairToken_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("pairs"))),
    triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
    triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
    triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales")))
//    prunedGenToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
//    packedGenToken_(consumes<pat::PackedGenParticleCollection>(iConfig.getParameter<edm::InputTag>("packedGenParticles"))),
//    lheprodToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheprod"))),
//    PileupSummaryInfoToken_(consumes<PileupSummaryInfoCollection>(iConfig.getParameter<edm::InputTag>("pileupinfo")))
{
   //now do what ever initialization is needed
}


maxi2ntuples::~maxi2ntuples()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
maxi2ntuples::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);
    if (vertices->empty()) return; // skip the event if no PV found
    const reco::Vertex &PV = vertices->front();

    edm::Handle<pat::METCollection> mets;
    iEvent.getByToken(metToken_, mets);
 //   const pat::MET &smet = mets->front();


    edm::Handle<pat::CompositeCandidateCollection> pairs;
    iEvent.getByToken(pairToken_, pairs);

    edm::Handle<edm::TriggerResults> triggerBits;
    iEvent.getByToken(triggerBits_, triggerBits);
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    iEvent.getByToken(triggerObjects_, triggerObjects);
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    iEvent.getByToken(triggerPrescales_, triggerPrescales);
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

    edm::Handle<pat::JetCollection> jets; 
    iEvent.getByToken(jetToken_, jets);

    std::cout << "################################################################\n";
    std::cout << "paircount = " << pairs->size() << 
        " run = " <<iEvent.id().run() << 
        " lumi = " << iEvent.luminosityBlock() <<
        " evt =  " << iEvent.id().event() <<
        " npv = " <<  vertices->size()
        << std::endl;



    for (const pat::CompositeCandidate &lP : *pairs){
             
        //const reco::Candidate * l1, *l2; 
        //l1 = lP.daughter(0); l2 = lP.daughter(1);
        const reco::Candidate *met = lP.daughter(2);
        if (lP.daughter(2)->hasMasterClone())
            const edm::Ref<reco::PFMETCollection> met = lP.daughter(2)->masterClone().castTo<edm::Ref<reco::PFMETCollection>>();
        else
            met = dynamic_cast<const pat::MET*>(lP.daughter(2));

        std::cout << "svfit: " << lP.userFloat("SVfitMass") << "; metpt: " << met->pt();

        const pat::Muon *mu  = dynamic_cast<const pat::Muon*>(lP.daughter(0)->masterClone().get());
    //    const pat::Electron *el = dynamic_cast<const pat::Electron*>(lP.daughter(1)->masterClone().get());
        const pat::Tau *tau  = dynamic_cast<const pat::Tau*>(lP.daughter(1)->masterClone().get());
        std::cout << "; mupt: " << mu->pt() <<"; taupt: "<< tau->pt()
             << "; mueta: " << mu->eta() <<"; taueta: "<< tau->eta()
             << "; mudxy: " << mu->innerTrack()->dxy( PV.position()) <<"; mudz: "<< mu->innerTrack()->dxy( PV.position())
             << "; isMediumMuon: " << mu->isMediumMuon() 
             << "; decayModeFindingNewDMs: " << tau->tauID("decayModeFindingNewDMs") 
             << "; taucharge: " << tau->charge() << "; deltaR: " << deltaR(mu->p4(), tau->p4()) 
             << std::endl;


    std::cout << "................................................................\n";
    }
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
        if(triggerBits->accept(i)){
            if(names.triggerName(i) == "HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1"
                    || names.triggerName(i) == "HLT_IsoMu24_eta2p1_v1")
                std::cout <<  names.triggerName(i) << std::endl;
        }
    }
    std::cout << "################################################################\n";


}


// ------------ method called once each job just before starting event loop  ------------
void maxi2ntuples::beginJob(){

}

// ------------ method called once each job just after ending the event loop  ------------
void 
maxi2ntuples::endJob() 
{
}


// ------------ method called when starting to processes a run  ------------
/*
void 
maxi2ntuples::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
maxi2ntuples::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
maxi2ntuples::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
maxi2ntuples::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
maxi2ntuples::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(maxi2ntuples);


