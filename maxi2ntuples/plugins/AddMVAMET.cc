// -*- C++ -*-
//
// Package:    m2n/AddMVAMET
// Class:      AddMVAMET
// 
/**\class AddMVAMET AddMVAMET.cc m2n/AddMVAMET/plugins/AddMVAMET.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Fri, 03 Apr 2015 09:13:55 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// user include files
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/Candidate/interface/ShallowCloneCandidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include <DataFormats/METReco/interface/CommonMETData.h>

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"

#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

#include <vector>
#include <string>

//
// class declaration
//

class AddMVAMET : public edm::EDProducer {
   public:
      explicit AddMVAMET(const edm::ParameterSet&);
      ~AddMVAMET();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
    //  edm::EDGetTokenT<reco::PFMETCollection> mvametToken_;
      edm::EDGetTokenT<pat::METCollection> metToken_;
      edm::EDGetTokenT<pat::CompositeCandidateCollection> pairsToken_;
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
AddMVAMET::AddMVAMET(const edm::ParameterSet& iConfig):
//     mvametToken_(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("mvamet"))),
     metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
     pairsToken_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("pairs")))
{
   //register your products
   //now do what ever other initialization is needed
   produces<pat::CompositeCandidateCollection>();
  
}


AddMVAMET::~AddMVAMET()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
AddMVAMET::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 

    std::unique_ptr<pat::CompositeCandidateCollection> pairsCollection(new pat::CompositeCandidateCollection());

//    edm::Handle<reco::PFMETCollection> mvamets;
//    iEvent.getByToken(mvametToken_, mvamets);
//    const pat::MET &met = mvamets->front();

    edm::Handle<pat::METCollection> mets;
    iEvent.getByToken(metToken_, mets);
    const pat::MET &met = mets->front();

    edm::Handle<pat::CompositeCandidateCollection> pairs;
    iEvent.getByToken(pairsToken_, pairs);

    if (!pairs.isValid()) return;

    for (const pat::CompositeCandidate &lP : *pairs) {
    
        pat::CompositeCandidate pair(lP);
        pair.addDaughter(met);
        AddFourMomenta addP4;
        addP4.set(pair);
        //std::cout << "PT: " << pair.pt() << std::endl;
        pairsCollection->push_back(pair);
        pair.clearDaughters();
    
    }

    iEvent.put(std::move(pairsCollection));
}

// ------------ method called once each job just before starting event loop  ------------
void 
AddMVAMET::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AddMVAMET::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
AddMVAMET::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
AddMVAMET::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
AddMVAMET::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
AddMVAMET::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AddMVAMET::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AddMVAMET);
