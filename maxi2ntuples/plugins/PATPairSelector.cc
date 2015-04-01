// -*- C++ -*-
//
// Package:    m2n/PATPairSelector
// Class:      PATPairSelector
// 
/**\class PATPairSelector PATPairSelector.cc m2n/PATPairSelector/plugins/PATPairSelector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Tue, 10 Mar 2015 14:55:08 GMT
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


#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/Utilities/interface/InputTag.h>
#include <DataFormats/Candidate/interface/ShallowCloneCandidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include <DataFormats/METReco/interface/CommonMETData.h>
#include <TauAnalysis/SVfitStandalone/interface/SVfitStandaloneAlgorithm.h>

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/SingleElementCollectionSelector.h"

#include "DataFormats/Math/interface/deltaR.h"


#include <vector>
#include <string>
//
// class declaration
//

class PATPairSelector : public edm::EDProducer {
   public:
      explicit PATPairSelector(const edm::ParameterSet&);
      ~PATPairSelector();

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
      edm::EDGetTokenT<pat::CompositeCandidateCollection> PairToken_;

      const StringCutObjectSelector<pat::Muon, true> muCut;
      const StringCutObjectSelector<pat::Electron, true> elCut;
      const StringCutObjectSelector<pat::Tau, true> tauCut;
      const double deltaR_;
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
PATPairSelector::PATPairSelector(const edm::ParameterSet& iConfig):
    PairToken_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("pairs"))),
    muCut(iConfig.getParameter<std::string>("muCut")),
    elCut(iConfig.getParameter<std::string>("elCut")),
    tauCut(iConfig.getParameter<std::string>("tauCut")),
    deltaR_(iConfig.getParameter<double>("deltaR_"))
{
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is neede
   produces<pat::CompositeCandidateCollection>();
}


PATPairSelector::~PATPairSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
PATPairSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    edm::Handle<pat::CompositeCandidateCollection> leptonPair;
    iEvent.getByToken(PairToken_, leptonPair);

    std::unique_ptr<pat::CompositeCandidateCollection> selectedPair(new pat::CompositeCandidateCollection());

    for (const pat::CompositeCandidate &lP : *leptonPair){

        const reco::Candidate * l1, *l2; 
        l1 = lP.daughter("leptonOne"); l2 = lP.daughter("leptonTwo");

        if ( deltaR(l1->eta(), l1->phi(), l2->eta(), l2->phi()) < deltaR_ )
            continue;

        if (l1->isMuon()){
            const pat::Muon *mu = dynamic_cast<const pat::Muon*>(lP.daughter("leptonOne"));    
            if (!muCut(*mu)) 
                continue;
        } else if(l1->isElectron()){
            const pat::Electron* el = dynamic_cast<const pat::Electron*>(lP.daughter("leptonOne"));
            if (!elCut(*el)) 
                continue;
        } else{
            const pat::Tau* tau = dynamic_cast<const pat::Tau*>(lP.daughter("leptonOne"));
            if (!tauCut(*tau)) 
                continue;
        }

        if (l2->isMuon()){
            const pat::Muon *mu = dynamic_cast<const pat::Muon*>(lP.daughter("leptonTwo"));    
            if (!muCut(*mu)) 
                continue;
        } else if(l2->isElectron()){
            const pat::Electron* el = dynamic_cast<const pat::Electron*>(lP.daughter("leptonTwo"));
            if (!elCut(*el)) 
                continue;
        } else{
            const pat::Tau* tau = dynamic_cast<const pat::Tau*>(lP.daughter("leptonTwo"));
            if (!tauCut(*tau))
                continue;
        }

        selectedPair->push_back(lP);
    }
 
   iEvent.put(std::move(selectedPair));
}

// ------------ method called once each job just before starting event loop  ------------
void 
PATPairSelector::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PATPairSelector::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
PATPairSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
PATPairSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
PATPairSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
PATPairSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PATPairSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PATPairSelector);
