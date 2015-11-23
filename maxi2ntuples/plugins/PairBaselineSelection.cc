// -*- C++ -*-
//
// Package:    m2n/PairBaselineSelection
// Class:      PairBaselineSelection
// 
/**\class PairBaselineSelection PairBaselineSelection.cc m2n/PairBaselineSelection/plugins/PairBaselineSelection.cc

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

class PairBaselineSelection : public edm::EDProducer {
   public:
      explicit PairBaselineSelection(const edm::ParameterSet&);
      ~PairBaselineSelection();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      float tautau(const reco::Vertex &, const pat::Tau*, const pat::Tau*);
      float mutau(const reco::Vertex &, const reco::Candidate*, const reco::Candidate*);
      float etau(const reco::Vertex &, const pat::Electron*, const pat::Tau*);
      float mumu(const pat::Muon*, const pat::Muon*);
      float ee(const pat::Electron*, const pat::Electron*);
      float emu(const reco::Vertex &, const pat::Electron*,const pat::Muon* );

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
      const bool mc;
      const std::string sample;
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
PairBaselineSelection::PairBaselineSelection(const edm::ParameterSet& iConfig):
    PairToken_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("pairs"))),
    vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
//    electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
    mc(iConfig.getParameter<bool>("mc")),
    sample(iConfig.getParameter<std::string>("sample")),
    eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap")))
    
{
    electronToken_ = mayConsume<edm::View<reco::GsfElectron> >
        (iConfig.getParameter<edm::InputTag>
         ("electrons"));
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


PairBaselineSelection::~PairBaselineSelection()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
PairBaselineSelection::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
 
 //   bool goodvertex = false;
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);
    const reco::Vertex &PV = vertices->front();


/*
    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByToken(muonToken_, muons);

    edm::Handle<pat::ElectronCollection> electrons;
    iEvent.getByToken(electronToken_, electrons);
*/
    edm::Handle<edm::View<reco::GsfElectron> > electrons;
    iEvent.getByToken(electronToken_,electrons);
    edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
    iEvent.getByToken(eleTightIdMapToken_,tight_id_decisions);


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
    /*
            if ( lP.daughter(0)->hasMasterClone() ) {   
               // reco::CandidateBaseRef master = lP.daughter(0)->masterClone();
                pat::MuonRef master = lP.daughter(0)->masterClone().castTo<pat::MuonRef>();
                std::cout << master->isMediumMuon();
            }
            */
        float passed = 1;

        if (l1->isMuon()){
            if (l2->isMuon()){
                //std::cout << "muon-muon channel" <<std::endl;
                const pat::Muon *mu  = dynamic_cast<const pat::Muon*>(l1->masterClone().get());
                const pat::Muon *mu_ = dynamic_cast<const pat::Muon*>(l2->masterClone().get());
                passed = mumu(mu, mu_);
            }
            else if(l2->isElectron()){
                //std::cout << "muon-electron channel" <<std::endl;
                const pat::Muon *mu  = dynamic_cast<const pat::Muon*>(l1->masterClone().get());
                const pat::Electron *e = dynamic_cast<const pat::Electron*>(l2->masterClone().get());
                for (size_t i = 0; i < electrons->size(); ++i){
                     const auto el = electrons->ptrAt(i);
                     if(el->p4() == e->p4())
                        if((*tight_id_decisions)[el]){
                            passed =  emu(PV, e, mu);
                            if (passed)
                                break;
                        }
                }
            }
            else {
                //std::cout << "muon-tau channell" <<std::endl;
                passed = mutau(PV, l1,  l2);
            }
        }
        else if(l1->isElectron()){
            if (l2->isMuon()){
                //std::cout << "muon-electron channel" <<std::endl;
                const pat::Electron *e = dynamic_cast<const pat::Electron*>(l1->masterClone().get());
                const pat::Muon *mu = dynamic_cast<const pat::Muon*>(l2->masterClone().get());
                for (size_t i = 0; i < electrons->size(); ++i){
                     const auto el = electrons->ptrAt(i);
                     if(el->p4() == e->p4())
                        if((*tight_id_decisions)[el]){
                            passed =  emu(PV, e, mu);
                            if (passed)
                                break;
                        }
                }
            }
            else if(l2->isElectron()){
                //std::cout << "electron-electron channel" <<std::endl;
                const pat::Electron *el = dynamic_cast<const pat::Electron*>(l1->masterClone().get());
                const pat::Electron *el_ = dynamic_cast<const pat::Electron*>(l2->masterClone().get());
                passed = ee(el,el_);
            }
            else {
                //std::cout << "electron-tau channel" <<std::endl;
                const pat::Electron *e = dynamic_cast<const pat::Electron*>(l1->masterClone().get());
                const pat::Tau *tau  = dynamic_cast<const pat::Tau*>(l2->masterClone().get());
                for (size_t i = 0; i < electrons->size(); ++i){
                     const auto el = electrons->ptrAt(i);
                     if(el->p4() == e->p4())
                        if((*tight_id_decisions)[el]){
                            passed = etau(PV, e,tau);
                            if (passed)
                                break;
                        }
                }
                //std::cout << *el << std::endl;
            }
        }
        else {
            if (l2->isMuon()){
                //std::cout << "muon-tau channel" <<std::endl;
                const pat::Tau *tau  = dynamic_cast<const pat::Tau*>(l1->masterClone().get());
                const pat::Muon *mu = dynamic_cast<const pat::Muon*>(l2->masterClone().get());
                passed = mutau(PV, mu,  tau);
            }
            else if(l2->isElectron()){
                //std::cout << "electron-tau channel" <<std::endl;
                const pat::Tau *tau  = dynamic_cast<const pat::Tau*>(l1->masterClone().get());
                const pat::Electron *e = dynamic_cast<const pat::Electron*>(l2->masterClone().get());
                for (size_t i = 0; i < electrons->size(); ++i){
                     const auto el = electrons->ptrAt(i);
                     if(el->p4() == e->p4())
                        if((*tight_id_decisions)[el])
                            passed = etau(PV, e,tau);
                }
            }
            else {
                //std::cout << "tau-tau channel" <<std::endl;
                const pat::Tau *tau  = dynamic_cast<const pat::Tau*>(l1->masterClone().get());
                const pat::Tau *tau_  = dynamic_cast<const pat::Tau*>(l2->masterClone().get());
                passed = tautau(PV, tau, tau_);
            }
        }

        pat::CompositeCandidate pair(lP);
        pair.addUserFloat("PairBaselineSelection", passed);
        selectedPair->push_back(pair);
        //if (pass)
        //    selectedPair->push_back(lP);
    //        std::cout << "Passeddddddddddddddddd!\n"; 
    }

    iEvent.put(std::move(selectedPair));
}

// ------------ method called once each job just before starting event loop  ------------
void 
PairBaselineSelection::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PairBaselineSelection::endJob() {
}


float PairBaselineSelection::tautau(const reco::Vertex &PV, const pat::Tau* tau, const pat::Tau* tau_){


    pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau->leadChargedHadrCand().get());
    pat::PackedCandidate const* packedTrailTauCand = dynamic_cast<pat::PackedCandidate const*>(tau_->leadChargedHadrCand().get());
    if(sample.find("spring15") != std::string::npos){
        if( 
            tau->tauID("decayModeFindingNewDMs") <= 0.5 || 
            fabs(packedLeadTauCand->dz()) >= 0.2 ||
            tau->pt() <= 45. ||
            fabs(tau->eta()) >= 2.1  ||
            tau_->tauID("decayModeFindingNewDMs") <= 0.5 || 
            fabs(packedTrailTauCand->dz()) >= 0.2 ||
            tau_->pt() <= 45. || 
            fabs(tau_->eta()) >= 2.1 

            || deltaR(tau->p4(), tau_->p4()) <= 0.5

          ) return false;


    }

    return true;
}



float PairBaselineSelection::mutau(const reco::Vertex &PV, const reco::Candidate* l1, const reco::Candidate* l2){


    if (deltaR(l1->p4(), l2->p4()) < 0.5) return false;
    
    const pat::Muon *mu  = dynamic_cast<const pat::Muon*>(l1->masterClone().get());
    const pat::Tau *tau  = dynamic_cast<const pat::Tau*>(l2->masterClone().get());

    //Spring15 MC samples
    if(sample.find("spring15") != std::string::npos){
        if(
            //mu
            mu->isMediumMuon()
            && mu->pt() > 18.
            && fabs(mu->eta()) < 2.1 
            && fabs(mu->muonBestTrack()->dxy( PV.position()) )  < 0.045 
            && fabs(mu->muonBestTrack()->dz(PV.position())) < 0.2  
            //!utilities::heppymuonID(*mu, "POG_ID_Medium") && 
            //utilities::relIso(*mu, 0.5) >= 0.1 && 

            //tau
            && tau->pt() > 20. 
            && tau->tauID("decayModeFindingNewDMs") > 0.5 
            && fabs(tau->eta()) < 2.3 
            && abs(tau->charge()) == 1){
                pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau->leadChargedHadrCand().get());
                if(fabs(packedLeadTauCand->dz())< 0.2)
                    return true;
        }
    }
    return false;

}

float PairBaselineSelection::etau(const reco::Vertex &PV, const pat::Electron* e, const pat::Tau* tau){

    pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau->leadChargedHadrCand().get());

    if(sample.find("spring15") != std::string::npos){
        if(
            e->pt() <= 23. || fabs(e->eta()) >= 2.1       
            || fabs(e->gsfTrack()->dxy(PV.position()))  >= 0.045 || fabs(e->gsfTrack()->dz(PV.position())) >= 0.2
            || e->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 1
            || !(e->passConversionVeto())

            || tau->pt() <= 20. 
            || fabs(tau->eta()) >= 2.3
            || tau->tauID("decayModeFindingNewDMs") <= 0.5      
            || fabs(packedLeadTauCand->dz()) >= 0.2 

            || deltaR(e->p4(), tau->p4()) <= 0.5
        ) return false;

    }
    return false;   
}
float PairBaselineSelection::mumu(const pat::Muon* mu, const pat::Muon* mu_){
    return false;   
}
float PairBaselineSelection::ee(const pat::Electron* e, const pat::Electron* e_){
    return false;
}
float PairBaselineSelection::emu(const reco::Vertex &PV, const pat::Electron* e,const pat::Muon* mu){
    if(sample.find("spring15") != std::string::npos){

        if(
                e->pt() <= 13. || fabs(e->eta() >=2.5) 
                || fabs(e->gsfTrack()->dxy(PV.position()))  >= 0.045 || fabs(e->gsfTrack()->dz(PV.position()) >= 0.2)
                || e->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 1
                || !(e->passConversionVeto())

                || mu->pt() <= 10. || fabs(mu->eta()) >= 2.4
                || fabs(mu->innerTrack()->dxy( PV.position() ))  >= 0.045 || fabs(mu->innerTrack()->dz(PV.position())) >= 0.2
                || !(mu->isMediumMuon())
                //|| !utilities::heppymuonID(*mu, "POG_ID_Medium") 

                || deltaR(e->p4(), mu->p4()) <= 0.3

                ) return false;
        
    }


    return false;
}


// ------------ method called when starting to processes a run  ------------
/*
void
PairBaselineSelection::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
PairBaselineSelection::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
PairBaselineSelection::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
PairBaselineSelection::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PairBaselineSelection::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PairBaselineSelection);
