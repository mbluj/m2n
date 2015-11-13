// -*- C++ -*-
//
// Package:    m2n/HLTforPair
// Class:      HLTforPair
// 
/**\class HLTforPair HLTforPair.cc m2n/HLTforPair/plugins/HLTforPair.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Tue, 31 Mar 2015 13:41:38 GMT
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
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "m2n/maxi2ntuples/interface/utilities.h"

#include <vector>
#include <string>
#include <math.h>


#include "CondFormats/L1TObjects/interface/L1GtTriggerMenuFwd.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "TMath.h"
#include "TMVA/MethodBDT.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//
// class declaration
//

class HLTforPair : public edm::EDProducer {
   public:
      explicit HLTforPair(const edm::ParameterSet&);
      ~HLTforPair();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      bool getpaths(std::string hlt, edm::Handle<edm::TriggerResults>& triggerBits , edm::Handle<pat::PackedTriggerPrescales>& triggerPrescales, const edm::TriggerNames &names);
            
      bool getfilters(const reco::Candidate* lepton, std::vector<std::string> filters, edm::Handle<pat::TriggerObjectStandAloneCollection>& triggerObjects, edm::Handle<pat::PackedTriggerPrescales>& triggerPrescales, const edm::TriggerNames &names);


      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<pat::CompositeCandidateCollection> PairToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      //edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
      edm::EDGetToken electronToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
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
HLTforPair::HLTforPair(const edm::ParameterSet& iConfig):
    PairToken_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("pairs"))),
    triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
    triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
    triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales")))
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


HLTforPair::~HLTforPair()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
HLTforPair::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
 

    edm::Handle<edm::TriggerResults> triggerBits;
    iEvent.getByToken(triggerBits_, triggerBits);
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    iEvent.getByToken(triggerObjects_, triggerObjects);
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    iEvent.getByToken(triggerPrescales_, triggerPrescales);
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

/*
    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByToken(muonToken_, muons);

    edm::Handle<pat::ElectronCollection> electrons;
    iEvent.getByToken(electronToken_, electrons);
*/

    edm::Handle<pat::CompositeCandidateCollection> leptonPair;
    iEvent.getByToken(PairToken_, leptonPair);

    std::unique_ptr<pat::CompositeCandidateCollection> selectedPair(new pat::CompositeCandidateCollection());


    if (!leptonPair.isValid()){
        iEvent.put(std::move(selectedPair));
        return;
    }

    

    for (const pat::CompositeCandidate &lP : *leptonPair){

        const reco::Candidate * l1, *l2; 
        l1 = lP.daughter(0); l2 = lP.daughter(1);
//        short unsigned int channelcase; 

        //mutau
        float HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1 = 0,
              HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v2 = 0,
              HLT_IsoMu24_eta2p1_v1 = 0,
              HLT_IsoMu24_eta2p1_v2 = 0,
              HLT_IsoMu27_v1 = 0,
              HLT_IsoMu17_eta2p1 = 0,
              HLT_IsoMu18_v1 = 0,
              HLT_IsoMu22_v1 = 0;


 //       float passed = 0;


        //if (!l1->isMuon())
        //    std::cout << "HLTforPair module ERROR: l1 is not muon!!!\n";

        if(
           getpaths("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1", triggerBits, triggerPrescales, names) 
           && getfilters(l1, {"hltOverlapFilterIsoMu17LooseIsoPFTau20", "hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09"}, triggerObjects, triggerPrescales, names) 
           && getfilters(l2, {"hltPFTau20TrackLooseIsoAgainstMuon", "hltOverlapFilterIsoMu17LooseIsoPFTau20"}, triggerObjects, triggerPrescales, names) 
        )  HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1 = 1;

        if(
           getpaths("HLT_IsoMu24_eta2p1_v1", triggerBits, triggerPrescales, names) 
           && getfilters(l1, {"hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09"}, triggerObjects, triggerPrescales, names) 
           && l1->pt() > 25.
        ) HLT_IsoMu24_eta2p1_v1 = 1;

        if(
            getpaths("HLT_IsoMu27_v1",  triggerBits, triggerPrescales, names)
            && getfilters(l1, {"hltL3crIsoL1sMu25L1f0L2f10QL3f27QL3trkIsoFiltered0p09"}, triggerObjects, triggerPrescales, names)
        ) HLT_IsoMu27_v1 = 1;

        if(
            getpaths("HLT_IsoMu17_eta2p1",  triggerBits, triggerPrescales, names)
            && getfilters(l1, {"hltL3crIsoL1sSingleMu16erL1f0L2f10QL3f17QL3trkIsoFiltered0p09"}, triggerObjects, triggerPrescales, names)
        ) HLT_IsoMu17_eta2p1 = 1;

        if(
           getpaths("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v2", triggerBits, triggerPrescales, names) 
           && getfilters(l1, {"hltOverlapFilterIsoMu17LooseIsoPFTau20", "hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09"}, triggerObjects, triggerPrescales, names) 
           && getfilters(l2, {"hltPFTau20TrackLooseIsoAgainstMuon", "hltOverlapFilterIsoMu17LooseIsoPFTau20"}, triggerObjects, triggerPrescales, names) 
        )  HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v2 = 1;

        if(
           getpaths("HLT_IsoMu24_eta2p1_v2", triggerBits, triggerPrescales, names) 
           && getfilters(l1, {"hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09"}, triggerObjects, triggerPrescales, names) 
           && l1->pt() > 25.
        ) HLT_IsoMu24_eta2p1_v2 = 1;

        if(
            getpaths("HLT_IsoMu18_v1",  triggerBits, triggerPrescales, names)
            && getfilters(l1, {"hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09"}, triggerObjects, triggerPrescales, names)
        ) HLT_IsoMu18_v1 = 1;
        
        if(
            getpaths("HLT_IsoMu22_v1",  triggerBits, triggerPrescales, names)
            && getfilters(l1, {"hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09"}, triggerObjects, triggerPrescales, names)
        ) HLT_IsoMu22_v1 = 1;



        pat::CompositeCandidate pair(lP);
        //pair.addUserFloat("HLTforPair", passed);
        pair.addUserFloat("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1",HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1);
        pair.addUserFloat("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v2",HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v2);
        pair.addUserFloat("HLT_IsoMu24_eta2p1_v1",HLT_IsoMu24_eta2p1_v1);
        pair.addUserFloat("HLT_IsoMu24_eta2p1_v2",HLT_IsoMu24_eta2p1_v2);
        pair.addUserFloat("HLT_IsoMu27_v1",HLT_IsoMu27_v1);
        pair.addUserFloat("HLT_IsoMu17_eta2p1",HLT_IsoMu17_eta2p1);
        pair.addUserFloat("HLT_IsoMu18_v1",HLT_IsoMu18_v1);
        pair.addUserFloat("HLT_IsoMu22_v1",HLT_IsoMu22_v1);

        selectedPair->push_back(pair);
        //if (pass)
        //    selectedPair->push_back(lP);
    //        std::cout << "Passeddddddddddddddddd!\n"; 
    }

    iEvent.put(std::move(selectedPair));
}

// ------------ method called once each job just before starting event loop  ------------
void 
HLTforPair::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HLTforPair::endJob() {
}




bool  HLTforPair::getpaths(std::string hlt, edm::Handle<edm::TriggerResults>& triggerBits , edm::Handle<pat::PackedTriggerPrescales>& triggerPrescales, const edm::TriggerNames &names){
    
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
 //       std::cout << names.triggerName(i) << std::endl;
        if(triggerBits->accept(i))
            if (names.triggerName(i) == hlt)
                return true;
    }

    return false;

}
bool  HLTforPair::getfilters(const reco::Candidate* lepton, std::vector<std::string> filters, edm::Handle<pat::TriggerObjectStandAloneCollection>& triggerObjects, edm::Handle<pat::PackedTriggerPrescales>& triggerPrescales, const edm::TriggerNames &names){
 
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
        obj.unpackPathNames(names);
        // for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << obj.filterLabels()[h]  << std::endl;
        if(deltaR(obj.triggerObject().p4(), lepton->p4()) < 0.5) 
            for(auto&& i : filters)
                if (obj.hasFilterLabel(i))
                    return true;
    } 

    return false;
}




// ------------ method called when starting to processes a run  ------------
/*
void
HLTforPair::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
HLTforPair::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
HLTforPair::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
HLTforPair::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HLTforPair::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HLTforPair);
