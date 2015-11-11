// -*- C++ -*-
//
// Package:    m2n/BestPairSelector
// Class:      BestPairSelector
// 
/**\class BestPairSelector BestPairSelector.cc m2n/maxi2ntuples/plugins/BestPairSelector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Wed, 11 Mar 2015 17:34:48 GMT
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

#include "m2n/maxi2ntuples/interface/utilities.h"

#include <vector>
#include <string>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"

//
// class declaration
//

class BestPairSelector : public edm::EDProducer {
   public:
      explicit BestPairSelector(const edm::ParameterSet&);
      ~BestPairSelector();

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
//      std::string channel;
//
      TH1F* reject_power;
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
BestPairSelector::BestPairSelector(const edm::ParameterSet& iConfig):
    PairToken_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("pairs")))
//    channel(iConfig.getParameter<std::string>("channel"))
{
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
   produces<pat::CompositeCandidateCollection>();
  
}


BestPairSelector::~BestPairSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
BestPairSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    //USE THIS CLASS AFTER CHANNEL SELECTION!!!
    using namespace edm;
    edm::Handle<pat::CompositeCandidateCollection> leptonPair;
    iEvent.getByToken(PairToken_, leptonPair);

    std::unique_ptr<pat::CompositeCandidateCollection> selectedPair(new pat::CompositeCandidateCollection());


    if (!leptonPair.isValid() || leptonPair->size()==0){
        iEvent.put(std::move(selectedPair));
        return;
    }

    
    if(leptonPair->size() == 1){
        reject_power->Fill(1.,1.);
        selectedPair->push_back(leptonPair->front());
        iEvent.put(std::move(selectedPair));
        return;   
    }

    pat::CompositeCandidateCollection buf;
    buf.push_back(leptonPair->front());
    for (const pat::CompositeCandidate &lP : *leptonPair){
        const reco::Candidate * l = lP.daughter(0);
        static float iso_l = 0, iso_max = 0;
        if(l->isMuon()){
            iso_l = utilities::relIso(*(dynamic_cast<const pat::Muon*>(l->masterClone().get())), 0.5);
            iso_max = utilities::relIso(*(dynamic_cast<const pat::Muon*>((buf.front()).daughter(0)->masterClone().get())), 0.5);
        }
        else if(l->isElectron()){
            iso_l = utilities::relIso(*(dynamic_cast<const pat::Electron*>(l->masterClone().get())), 0.5);
            iso_max =  utilities::relIso(*(dynamic_cast<const pat::Electron*>((buf.front()).daughter(0)->masterClone().get())), 0.5);
        }else{
            iso_l = (dynamic_cast<const pat::Tau*>(l->masterClone().get()))->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
            iso_max = (dynamic_cast<const pat::Tau*>((buf.front()).daughter(0)->masterClone().get()))->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
        }
        if(iso_l < iso_max){buf.clear(); buf.push_back(lP);}
        if(iso_l == iso_max){buf.push_back(lP);}
    }
    if(buf.size() == 1){
        reject_power->Fill(2.,1.);
        selectedPair->push_back(buf.front());
        iEvent.put(std::move(selectedPair));
        return;
    }
    pat::CompositeCandidateCollection buf_;
    buf_.push_back(buf.front());
    for (const pat::CompositeCandidate &lP : buf){
        static float pt_l = lP.daughter(0)->pt();
        static float pt_max = buf_.front().daughter(0)->pt();
        if(pt_l > pt_max){ buf_.clear(); buf_.push_back(lP);}
        if(pt_l == pt_max){ buf_.push_back(lP);}

    }
    if(buf_.size() == 1){
        reject_power->Fill(3.,1.);
        selectedPair->push_back(buf_.front());
        iEvent.put(std::move(selectedPair));
        return;   
    }
    buf.clear();
    buf.push_back(buf_.front());
    for (const pat::CompositeCandidate &lP : buf_){
        const reco::Candidate * l = lP.daughter(1);
        static float iso_l = 0, iso_max = 0;
        if(l->isMuon()){
            iso_l = utilities::relIso(*(dynamic_cast<const pat::Muon*>(l->masterClone().get())), 0.5);
            iso_max = utilities::relIso(*(dynamic_cast<const pat::Muon*>((buf.front()).daughter(1)->masterClone().get())), 0.5);
        }
        else if(l->isElectron()){
            iso_l = utilities::relIso(*(dynamic_cast<const pat::Electron*>(l->masterClone().get())), 0.5);
            iso_max =  utilities::relIso(*(dynamic_cast<const pat::Electron*>((buf.front()).daughter(1)->masterClone().get())), 0.5);
        }else{
            iso_l = (dynamic_cast<const pat::Tau*>(l->masterClone().get()))->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
            iso_max = (dynamic_cast<const pat::Tau*>((buf.front()).daughter(1)->masterClone().get()))->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
        }
        if(iso_l < iso_max){ buf.clear(); buf.push_back(lP);}
        if(iso_l == iso_max){ buf.push_back(lP);}
    }
    if(buf.size() == 1){
        reject_power->Fill(4.,1.);
        selectedPair->push_back(buf.front());
        iEvent.put(std::move(selectedPair));
        return;   
    }
    buf_.clear();
    buf_.push_back(buf.front());
    for (const pat::CompositeCandidate &lP : buf){
        static float pt_l = lP.daughter(1)->pt();
        static float pt_max = buf_.front().daughter(1)->pt();
        if(pt_l > pt_max){ buf_.clear(); buf_.push_back(lP);}
        if(pt_l == pt_max){ buf_.push_back(lP);}
    }
    if(buf_.size() == 1){
        reject_power->Fill(5.,1.);
        selectedPair->push_back(buf_.front());
        iEvent.put(std::move(selectedPair));
        return;   
    }

    reject_power->Fill(6.,1.);
    //std::cout << "Still more then one pair - do sth";
    selectedPair->push_back(buf_.front());
    iEvent.put(std::move(selectedPair));
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
BestPairSelector::beginJob()
{
   edm::Service<TFileService> theFileService;
   reject_power = theFileService->make<TH1F>("reject_power","reject_power",10,0,10);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BestPairSelector::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
BestPairSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
BestPairSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
BestPairSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
BestPairSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BestPairSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BestPairSelector);
