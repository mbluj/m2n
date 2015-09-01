// -*- C++ -*-
//
// Package:    m2n/EventNumberFilter
// Class:      EventNumberFilter
// 
/**\class EventNumberFilter EventNumberFilter.cc m2n/EventNumberFilter/plugins/EventNumberFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Wed, 11 Feb 2015 16:30:28 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <FWCore/Framework/interface/ESHandle.h>
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




//
// class declaration
//

class EventNumberFilter : public edm::EDFilter {
   public:
      explicit EventNumberFilter(const edm::ParameterSet&);
      ~EventNumberFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
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
EventNumberFilter::EventNumberFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

}


EventNumberFilter::~EventNumberFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
EventNumberFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    if(
 iEvent.id().event() ==            40501
|| iEvent.id().event() ==            66443
|| iEvent.id().event() ==            90387
|| iEvent.id().event() ==            93091
|| iEvent.id().event() ==            103917
|| iEvent.id().event() ==            107309
|| iEvent.id().event() ==            118085
|| iEvent.id().event() ==            130529
|| iEvent.id().event() ==            167363
|| iEvent.id().event() ==            187115
|| iEvent.id().event() ==            249638
|| iEvent.id().event() ==            255111
|| iEvent.id().event() ==            297149
|| iEvent.id().event() ==            342767
|| iEvent.id().event() ==            396690
)
        return true;
    else
        return false;



/*
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
   return true;
*/
}

// ------------ method called once each job just before starting event loop  ------------
void 
EventNumberFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EventNumberFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
EventNumberFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
EventNumberFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
EventNumberFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
EventNumberFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EventNumberFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(EventNumberFilter);
