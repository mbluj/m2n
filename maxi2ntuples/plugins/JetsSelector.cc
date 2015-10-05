// -*- C++ -*-
//
// Package:    m2n/JetsSelector
// Class:      JetsSelector
// 
/**\class JetsSelector JetsSelector.cc m2n/JetsSelector/plugins/JetsSelector.cc

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

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"


#include <vector>
#include <string>
//
// class declaration
//

class JetsSelector : public edm::EDProducer {
   public:
      explicit JetsSelector(const edm::ParameterSet&);
      ~JetsSelector();

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
      edm::EDGetTokenT<pat::JetCollection> jetToken_;

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
JetsSelector::JetsSelector(const edm::ParameterSet& iConfig):
    jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets")))
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
   produces<pat::JetCollection>();
  
}


JetsSelector::~JetsSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
JetsSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    edm::Handle<pat::JetCollection> jets; 
    iEvent.getByToken(jetToken_, jets);
    std::unique_ptr<pat::JetCollection> selectedjets(new pat::JetCollection());
    

    for (const pat::Jet &j : *jets){


        float NHF = j.neutralHadronEnergyFraction();
        float NEMF = j.neutralEmEnergyFraction();
        float NumConst = j.chargedMultiplicity() +j.neutralMultiplicity();
        //float NumNeutralParticles = j.neutralMultiplicity(); 

        float CHF = j.chargedHadronEnergyFraction();
        float CEMF = j.chargedEmEnergyFraction();
        float CHM = j.chargedMultiplicity();
        float MUF = j.muonEnergyFraction();
        
        float AJE = fabs(j.eta());
        /*
        if (
                ((NHF<0.99 && NEMF<0.99 && NumConst>1) && ((fabs(j.eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || fabs(j.eta())>2.4) && fabs(j.eta())<=3.0) 
             || ((NEMF<0.90 && NumNeutralParticles >10 && fabs(j.eta())>3.0 ))
           )
            selectedjets->push_back(j);
        
        */
        /*
        if ( AJE <= 3. && !(NHF < 0.99 && NEMF < 0.99 && NumConst> 1))
            continue;
        if ( AJE <= 2.4 && !(CHF > 0. && CHM > 0. && CEMF < 0.99))
            continue;
        if ( AJE > 3. && !(NEMF < 0.9 && NumNeutralParticles > 10))
            continue;
 
        */       
        if(!(NHF<0.99 && NEMF<0.99 && NumConst>1 && MUF<0.8))
            continue;
        if (AJE <= 2.4 && !( CHF>0 && CHM>0 && CEMF<0.99))
            continue;


        selectedjets->push_back(j);
        
    }


    iEvent.put(std::move(selectedjets));
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
JetsSelector::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetsSelector::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
JetsSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
JetsSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
JetsSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
JetsSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JetsSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetsSelector);
