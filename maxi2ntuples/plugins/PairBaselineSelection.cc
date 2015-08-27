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

      bool tautau(const reco::Vertex &, edm::Handle<edm::TriggerResults>&, edm::Handle<pat::TriggerObjectStandAloneCollection>&, edm::Handle<pat::PackedTriggerPrescales>&,const edm::TriggerNames &  ,const pat::Tau*, const pat::Tau*);
      bool mutau(const reco::Vertex &,const L1GtTriggerMenu*, edm::Handle<edm::TriggerResults>&, edm::Handle<pat::TriggerObjectStandAloneCollection>&, edm::Handle<pat::PackedTriggerPrescales>&, const edm::TriggerNames &, const pat::Muon*, const pat::Tau*);
      bool etau(const reco::Vertex &, edm::Handle<edm::TriggerResults>&, edm::Handle<pat::TriggerObjectStandAloneCollection>&, edm::Handle<pat::PackedTriggerPrescales>&,const edm::TriggerNames &,const pat::Electron*, const pat::Tau*);
      bool mumu(const pat::Muon*, const pat::Muon*);
      bool ee(const pat::Electron*, const pat::Electron*);
      bool emu(const reco::Vertex &, edm::Handle<edm::TriggerResults>&, edm::Handle<pat::TriggerObjectStandAloneCollection>&, edm::Handle<pat::PackedTriggerPrescales>&, const edm::TriggerNames &, const pat::Electron*,const pat::Muon* );

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<pat::CompositeCandidateCollection> PairToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
      const bool mc;
      const std::string sample;
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
    electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
    triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
    triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
    triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
    mc(iConfig.getParameter<bool>("mc")),
    sample(iConfig.getParameter<std::string>("sample"))
    
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
    if (!leptonPair.isValid()) return;


      edm::ESHandle<L1GtTriggerMenu> menuRcd;
      iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
      const L1GtTriggerMenu* menu = menuRcd.product();




    std::unique_ptr<pat::CompositeCandidateCollection> selectedPair(new pat::CompositeCandidateCollection());
    

    for (const pat::CompositeCandidate &lP : *leptonPair){

        const reco::Candidate * l1, *l2; 
        l1 = lP.daughter(0); l2 = lP.daughter(1);
//        short unsigned int channelcase; 

    bool pass = false; 

    if (l1->isMuon()){
        if (l2->isMuon()){
            //std::cout << "muon-muon channel" <<std::endl;
            const pat::Muon *mu  = dynamic_cast<const pat::Muon*>(l1->masterClone().get());
            const pat::Muon *mu_ = dynamic_cast<const pat::Muon*>(l2->masterClone().get());
            pass = mumu(mu, mu_);
        }
        else if(l2->isElectron()){
            //std::cout << "muon-electron channel" <<std::endl;
            const pat::Muon *mu  = dynamic_cast<const pat::Muon*>(l1->masterClone().get());
            const pat::Electron *el = dynamic_cast<const pat::Electron*>(l2->masterClone().get());
            pass =  emu(PV, triggerBits,triggerObjects,  triggerPrescales, names, el, mu);
        }
        else {
            //std::cout << "muon-tau channell" <<std::endl;
            const pat::Muon *mu  = dynamic_cast<const pat::Muon*>(l1->masterClone().get());
            const pat::Tau *tau  = dynamic_cast<const pat::Tau*>(l2->masterClone().get());
            pass = mutau(PV, menu ,triggerBits,triggerObjects,  triggerPrescales, names,  mu,  tau);
        }
    }
    else if(l1->isElectron()){
        if (l2->isMuon()){
            //std::cout << "muon-electron channel" <<std::endl;
            const pat::Electron *el = dynamic_cast<const pat::Electron*>(l1->masterClone().get());
            const pat::Muon *mu = dynamic_cast<const pat::Muon*>(l2->masterClone().get());
            pass = emu(PV, triggerBits, triggerObjects, triggerPrescales, names, el, mu);
        }
        else if(l2->isElectron()){
            //std::cout << "electron-electron channel" <<std::endl;
            const pat::Electron *el = dynamic_cast<const pat::Electron*>(l1->masterClone().get());
            const pat::Electron *el_ = dynamic_cast<const pat::Electron*>(l2->masterClone().get());
            pass = ee(el,el_);
        }
        else {
            //std::cout << "electron-tau channel" <<std::endl;
            const pat::Electron *el = dynamic_cast<const pat::Electron*>(l1->masterClone().get());
            const pat::Tau *tau  = dynamic_cast<const pat::Tau*>(l2->masterClone().get());
            pass = etau(PV, triggerBits, triggerObjects, triggerPrescales,names ,el,tau);
        }
    }
    else {
        if (l2->isMuon()){
            //std::cout << "muon-tau channel" <<std::endl;
            const pat::Tau *tau  = dynamic_cast<const pat::Tau*>(l1->masterClone().get());
            const pat::Muon *mu = dynamic_cast<const pat::Muon*>(l2->masterClone().get());
            pass = mutau(PV, menu,triggerBits,triggerObjects,  triggerPrescales, names,  mu,  tau);
        }
        else if(l2->isElectron()){
            //std::cout << "electron-tau channel" <<std::endl;
            const pat::Tau *tau  = dynamic_cast<const pat::Tau*>(l1->masterClone().get());
            const pat::Electron *el = dynamic_cast<const pat::Electron*>(l2->masterClone().get());
            pass = etau(PV, triggerBits, triggerObjects, triggerPrescales,names ,el,tau);
        }
        else {
            //std::cout << "tau-tau channel" <<std::endl;
            const pat::Tau *tau  = dynamic_cast<const pat::Tau*>(l1->masterClone().get());
            const pat::Tau *tau_  = dynamic_cast<const pat::Tau*>(l2->masterClone().get());
            pass = tautau(PV, triggerBits, triggerObjects, triggerPrescales,names , tau, tau_);
        }
    }

    if (pass)
        selectedPair->push_back(lP);
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


bool PairBaselineSelection::tautau(const reco::Vertex &PV, edm::Handle<edm::TriggerResults>& triggerBits, edm::Handle<pat::TriggerObjectStandAloneCollection>& triggerObjects, edm::Handle<pat::PackedTriggerPrescales>& triggerPrescales,const edm::TriggerNames &names , const pat::Tau* tau, const pat::Tau* tau_){


    pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau->leadChargedHadrCand().get());
    pat::PackedCandidate const* packedTrailTauCand = dynamic_cast<pat::PackedCandidate const*>(tau_->leadChargedHadrCand().get());
    if(sample.find("spring15") != std::string::npos){
        if( 
            tau_->tauID("decayModeFindingNewDMs") <= 0.5 || 
            abs(packedLeadTauCand->dz()) >= 0.2 ||
            //tau->vertex().z() != PV.z() ||
            tau->pt() <= 45. ||
            abs(tau->eta()) >= 2.1  ||
            tau_->tauID("decayModeFindingNewDMs") <= 0.5 || 
            abs(packedTrailTauCand->dz()) >= 0.2 ||
            //tau_->vertex().z() != PV.z() ||
            tau_->pt() <= 45. || 
            abs(tau_->eta()) >= 3.1 

            || deltaR(tau->p4(), tau_->p4()) <= 0.5

          ) return false;

        std::string path = "";
        for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
            if(triggerBits->accept(i)){
                path += names.triggerName(i);
            }
        }

        if (path.find("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v1") != std::string::npos){ 
            for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
                obj.unpackPathNames(names);
                // for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << obj.filterLabels()[h]  << std::endl;
                if( 
                    deltaR( obj.triggerObject().p4(), tau->p4()) < 0.5 
                    && (std::find(obj.filterLabels().begin(), obj.filterLabels().end(), "hltDoublePFTau40TrackPt1MediumIsolationDz02Reg")!=obj.filterLabels().end())
                ){
                    for (pat::TriggerObjectStandAlone obj_ : *triggerObjects) {
                        obj_.unpackPathNames(names);
                        if(
                            deltaR(obj_.triggerObject().p4(), tau_->p4()) < 0.5 
                            && (std::find(obj.filterLabels().begin(), obj.filterLabels().end(), "hltDoublePFTau40TrackPt1MediumIsolationDz02Reg")!=obj.filterLabels().end())
                        ) return true;
                    }
                    break;
                }
            } 
        }
    }

    if(sample.find("phys14") != std::string::npos){
        if( 
            (tau->tauID("decayModeFindingNewDMs") <= 0.5 && tau->tauID("decayModeFinding") <= 0.5)
            //|| tau->vertex().z() != PV.z()
            || abs(tau->vertex().z() - PV.z()) < 0.2 
            || tau->pt() <= 45. 
            || abs(tau->eta()) >= 2.1 

            || (tau_->tauID("decayModeFindingNewDMs") <= 0.5 && tau_->tauID("decayModeFinding") <= 0.5)
            || tau_->vertex().z() != PV.z()
            ||tau_->pt() <= 45. 
            || abs(tau_->eta()) >= 2.1 

            || deltaR(tau->p4(), tau_->p4()) <= 0.5

          ) return false;

        std::string path = "";
        for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
            if(triggerBits->accept(i)){
                path += names.triggerName(i);
            }
        }

        if (path.find("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v1") != std::string::npos){ 
            for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
                obj.unpackPathNames(names);
                // for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << obj.filterLabels()[h]  << std::endl;
                if( 
                    deltaR( obj.triggerObject().p4(), tau->p4()) < 0.5 
                    && (std::find(obj.filterLabels().begin(), obj.filterLabels().end(), "hltDoublePFTau40TrackPt1MediumIsolationDz02Reg")!=obj.filterLabels().end())
                ){
                    for (pat::TriggerObjectStandAlone obj_ : *triggerObjects) {
                        obj_.unpackPathNames(names);
                        if(
                            deltaR(obj_.triggerObject().p4(), tau_->p4()) < 0.5 
                            && (std::find(obj.filterLabels().begin(), obj.filterLabels().end(), "hltDoublePFTau40TrackPt1MediumIsolationDz02Reg")!=obj.filterLabels().end())
                        ) return true;
                    }
                    break;
                }
            } 
        }
    }

    return false;
}



bool PairBaselineSelection::mutau(const reco::Vertex &PV, const L1GtTriggerMenu* menu, edm::Handle<edm::TriggerResults>& triggerBits, edm::Handle<pat::TriggerObjectStandAloneCollection>& triggerObjects, edm::Handle<pat::PackedTriggerPrescales>& triggerPrescales,const edm::TriggerNames &names, const pat::Muon* mu, const pat::Tau* tau){

    pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau->leadChargedHadrCand().get());
    //Phys'14 MC samples
    if(sample.find("phys14") != std::string::npos){
        if(
            //mu
            mu->pt() <= 18. ||
            abs(mu->eta()) >= 2.1 ||
            abs(mu->innerTrack()->dxy( PV.position()) )  >= 0.045 || 
            abs(mu->innerTrack()->dz(PV.position())) >= 0.2  
            || !mu->isMediumMuon()

            //!utilities::heppymuonID(*mu, "POG_ID_Medium") || 
            //utilities::relIso(*mu, 0.5) >= 0.1 || 
            
            //tau
            || tau->pt() <= 20. 
            || abs(tau->eta()) >= 2.3 
            || (tau->tauID("decayModeFindingNewDMs") <= 0.5 && tau->tauID("decayModeFinding") <= 0.5)
 //           || abs(tau->vertex().z() - PV.z()) > 0.2 

            //abs(packedLeadTauCand->dz()) >= 0.2
            
            || deltaR(mu->p4(), tau->p4()) <=0.5
        ) return false;

        std::string path = "";
        for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
            if(triggerBits->accept(i)){
                path += names.triggerName(i);
                //std::cout <<  names.triggerName(i) << std::endl;
            }
        }

        std::string hlta = mc ? "HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1" : "HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v2";
        if (path.find(hlta) != std::string::npos){ 
            for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
                obj.unpackPathNames(names);
                // for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << obj.filterLabels()[h]  << std::endl;
                if( 
                    deltaR( obj.triggerObject().p4(), mu->p4()) < 0.5 && 
                    std::find(obj.filterLabels().begin(), obj.filterLabels().end(), "hltOverlapFilterIsoMu17LooseIsoPFTau20")!=obj.filterLabels().end() 
                ){
                    for (pat::TriggerObjectStandAlone obj_ : *triggerObjects) {
                        obj_.unpackPathNames(names);
                        if(
                            deltaR(obj_.triggerObject().p4(), tau->p4()) < 0.5 && 
                            (std::find(obj_.filterLabels().begin(), obj_.filterLabels().end(), "hltL1sMu16erTauJet20er ")!=obj_.filterLabels().end() ||
                            std::find(obj_.filterLabels().begin(), obj_.filterLabels().end(), "hltOverlapFilterIsoMu17LooseIsoPFTau20")!=obj_.filterLabels().end() )
                        ) return true;
                    }
                    break;
                }
            } 
        }
        if (path.find("HLT_IsoMu24_eta2p1_IterTrk02_v1") != std::string::npos  ){
            for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
                obj.unpackPathNames(names);
                if( 
                    deltaR( obj.triggerObject().p4(), mu->p4()) < 0.5 && 
                    std::find(obj.filterLabels().begin(), obj.filterLabels().end(), "hltL3crIsoL1sMu20Eta2p1L1f0L2f20QL3f24QL3crIsoRhoFiltered0p15IterTrk02")!=obj.filterLabels().end() 
                ) return true;
            } 
        }
    }

    //Spring15 MC samples
    if(sample.find("spring15") != std::string::npos){
        if(
            //mu
            mu->pt() <= 18.
            || abs(mu->eta()) >= 2.1 
            || abs(mu->muonBestTrack()->dxy( PV.position()) )  >= 0.045 
            || abs(mu->muonBestTrack()->dz(PV.position())) >= 0.2  
            || !mu->isMediumMuon()
            //!utilities::heppymuonID(*mu, "POG_ID_Medium") || 
            //utilities::relIso(*mu, 0.5) >= 0.1 || 
            
            //tau
            || tau->pt() <= 20. 
            || abs(tau->eta()) >= 2.3 
            || tau->tauID("decayModeFindingNewDMs") <= 0.5 
            || abs(packedLeadTauCand->dz()) >= 0.2
            //abs(tau->vertex().z() - PV.z()) < 0.2 
            || abs(tau->charge()) != 1

            || deltaR(mu->p4(), tau->p4()) <= 0.5 
        ) return false;

        std::string path = "";
        for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
            if(triggerBits->accept(i)){
                path += names.triggerName(i);
                //std::cout <<  names.triggerName(i) << std::endl;
            }
        }

        std::string hlta = mc ? "HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1" : "HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v2";
        std::string hltb = mc ? "HLT_IsoMu27_v1" : "HLT_IsoMu24_eta2p1_v2";
        if (path.find(hlta) != std::string::npos){
            for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
                obj.unpackPathNames(names);
                // for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << obj.filterLabels()[h]  << std::endl;
                if( 
                    deltaR( obj.triggerObject().p4(), mu->p4()) < 0.5 
                    && std::find(obj.filterLabels().begin(), obj.filterLabels().end(), "hltOverlapFilterIsoMu17LooseIsoPFTau20")!=obj.filterLabels().end() 
                    && std::find(obj.filterLabels().begin(), obj.filterLabels().end(), "hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09 ")!=obj.filterLabels().end() 
                ){
                    for (pat::TriggerObjectStandAlone obj_ : *triggerObjects) {
                        obj_.unpackPathNames(names);
                        if(
                            deltaR(obj_.triggerObject().p4(), tau->p4()) < 0.5 
                            && std::find(obj_.filterLabels().begin(), obj_.filterLabels().end(), "hltPFTau20TrackLooseIsoAgainstMuon")!=obj_.filterLabels().end() 
                            && std::find(obj_.filterLabels().begin(), obj_.filterLabels().end(), "hltOverlapFilterIsoMu17LooseIsoPFTau20")!=obj_.filterLabels().end() 
                        ) return true;
                    }
                    break;
                }
            } 
        }
        if (path.find(hltb) != std::string::npos  ){
            for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
                obj.unpackPathNames(names);
                if(
                    deltaR( obj.triggerObject().p4(), mu->p4()) < 0.5 
                    && std::find(obj.filterLabels().begin(), obj.filterLabels().end(), "hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09")!=obj.filterLabels().end() 
                    && mu->pt() > 25.
                ) return true;
            }
        }
    }

    return false;


    // VETO
    // Second lepton veto
    /*
    for (const pat::Muon &vetomu : *muons) {
        if (abs(vetomu.px() - mu->px()) +  abs(vetomu.py() - mu->py()) + abs(vetomu.pz() - mu->pz()) > 0.01 && !vetomu.innerTrack().isNull()){
            if (vetomu.charge() != mu->charge()) {
                if (vetomu.pt() > 15 && abs(vetomu.eta()) < 2.4  && 
                    vetomu.isGlobalMuon() && vetomu.isTrackerMuon() && 
                    vetomu.isPFMuon() &&
                    abs(vetomu.innerTrack()->dz( PV.position()) ) < 0.2 &&
                    abs(vetomu.innerTrack()->dxy( PV.position()) ) < 0.045  &&
                    utilities::relIso(vetomu, 0.5) < 0.3)
                    continue;
            }
            // Third lepton veto : next muon veto
            if (vetomu.pt() > 10           &&
                abs(vetomu.eta()) < 2.4    &&
                utilities::heppymuonID(vetomu, "POG_ID_Medium")   &&
                abs(vetomu.innerTrack()->dxy( PV.position())) < 0.045  &&
                abs(vetomu.innerTrack()->dz( PV.position())) < 0.2   &&
                utilities::relIso(vetomu, 0.5) < 0.3)
                continue;
        }
    }
    // Third lepton veto: Next electron veto
    for (const pat::Electron &vetoel : *electrons) {
        if (!vetoel.gsfTrack().isNull() &&  abs(vetoel.px() - mu->px()) +  abs(vetoel.py() - mu->py()) + abs(vetoel.pz() - mu->pz()) > 0.01){
            if (vetoel.pt() > 10                                     &&
                abs(vetoel.eta()) < 2.5                              &&
                abs(vetoel.gsfTrack()->dxy(PV.position())) < 0.045                          &&
                abs(vetoel.gsfTrack()->dz(PV.position())) < 0.2                             &&
           //  FIXME    vetoel.cutBasedId('POG_PHYS14_25ns_v1_Veto')       &&
                utilities::relIso(vetoel, 0.5) < 0.3 )
                    continue;
        }
    }
    */
}
bool PairBaselineSelection::etau(const reco::Vertex &PV, edm::Handle<edm::TriggerResults>& triggerBits, edm::Handle<pat::TriggerObjectStandAloneCollection>& triggerObjects, edm::Handle<pat::PackedTriggerPrescales>& triggerPrescales,const edm::TriggerNames &names, const pat::Electron* e, const pat::Tau* tau){

    pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau->leadChargedHadrCand().get());

    if(sample.find("spring15") != std::string::npos){
        if(
            e->pt() <= 23. || abs(e->eta()) >= 2.1       
            || abs(e->gsfTrack()->dxy(PV.position()))  >= 0.045 || abs(e->gsfTrack()->dz(PV.position())) >= 0.2
            || e->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 1
            || e->passConversionVeto()

            || tau->pt() <= 20 
            || abs(tau->eta()) >= 2.3
            || tau->tauID("decayModeFindingNewDMs") <= 0.5      
            || abs(packedLeadTauCand->dz()) >= 0.2 

            || deltaR(e->p4(), tau->p4()) <= 0.5
        ) return false;

        std::string path = "";
        for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
            if(triggerBits->accept(i)){
                path += names.triggerName(i);
            }
        }
        
        /*
        //std::string seed = "";
        for (CItAlgo algo = menu->gtAlgorithmMap().begin(); algo!=menu->gtAlgorithmMap().end(); ++algo) {
           // std::cout << "Name: " << (algo->second).algoName() << " Alias: " << (algo->second).algoAlias() << std::endl;
           //seed += (algo->second).algoName();
        }
        */
        if(mc){
            if (path.find("HLT_Ele22_eta2p1_WP75_Gsf_LooseIsoPFTau20_v1") != std::string::npos){
                for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
                    obj.unpackPathNames(names);
                    // for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << obj.filterLabels()[h]  << std::endl;
                    if( 
                        deltaR( obj.triggerObject().p4(), e->p4()) < 0.5  
                        && std::find(obj.filterLabels().begin(), obj.filterLabels().end(), "hltEle22WP75L1IsoEG20erTau20erGsfTrackIsoFilter ")!=obj.filterLabels().end() 
                        && std::find(obj.filterLabels().begin(), obj.filterLabels().end(), "hltOverlapFilterIsoEle22WP75GsfLooseIsoPFTau20")!=obj.filterLabels().end() 
                    ){
                        for (pat::TriggerObjectStandAlone obj_ : *triggerObjects) {
                            obj_.unpackPathNames(names);
                            if(
                                deltaR(obj_.triggerObject().p4(), tau->p4()) < 0.5 && 
                                (std::find(obj_.filterLabels().begin(), obj_.filterLabels().end(), "hltPFTau20TrackLooseIso")!=obj_.filterLabels().end() 
                                && std::find(obj_.filterLabels().begin(), obj_.filterLabels().end(), "hltOverlapFilterIsoEle22WP75GsfLooseIsoPFTau20")!=obj_.filterLabels().end() )
                            ) return true;
                        }
                        break;
                    }
                } 
            }
            if (path.find("HLT_Ele32_eta2p1_WP75_Gsf_v1") != std::string::npos  ){
                for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
                    obj.unpackPathNames(names);
                    if( 
                        deltaR( obj.triggerObject().p4(), e->p4()) < 0.5 
                        && std::find(obj.filterLabels().begin(), obj.filterLabels().end(), "hltEle32WP75GsfTrackIsoFilter")!=obj.filterLabels().end() 
                        && e->pt() > 33.
                    ) return true;
                } 
            }
        }
        else {
            if (path.find("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v1") != std::string::npos){
                for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
                    obj.unpackPathNames(names);
                    // for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << obj.filterLabels()[h]  << std::endl;
                    if( 
                        deltaR( obj.triggerObject().p4(), e->p4()) < 0.5  
                        && std::find(obj.filterLabels().begin(), obj.filterLabels().end(), "hltSingleEle22WPLooseGsfTrackIsoFilter")!=obj.filterLabels().end() 
                        && std::find(obj.filterLabels().begin(), obj.filterLabels().end(), "hltOverlapFilterIsoEle22WPLooseGsfLooseIsoPFTau20")!=obj.filterLabels().end() 
                    ){
                        for (pat::TriggerObjectStandAlone obj_ : *triggerObjects) {
                            obj_.unpackPathNames(names);
                            if(
                                deltaR(obj_.triggerObject().p4(), tau->p4()) < 0.5 && 
                                (std::find(obj_.filterLabels().begin(), obj_.filterLabels().end(), "hltPFTau20TrackLooseIso")!=obj_.filterLabels().end() 
                                && std::find(obj_.filterLabels().begin(), obj_.filterLabels().end(), "hltOverlapFilterIsoEle22WPLooseGsfLooseIsoPFTau20")!=obj_.filterLabels().end() )
                            ) return true;
                        }
                        break;
                    }
                }
            }
            if (path.find("HLT_Ele32_eta2p1_WPTight_Gsf_v1") != std::string::npos  ){
                for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
                    obj.unpackPathNames(names);
                    if( 
                        deltaR( obj.triggerObject().p4(), e->p4()) < 0.5 
                        && std::find(obj.filterLabels().begin(), obj.filterLabels().end(), "hltEle32WPTightGsfTrackIsoFilter")!=obj.filterLabels().end() 
                        && e->pt() > 33.
                    ) return true;
                } 
            }
        }
    }
    if(sample.find("phys14") != std::string::npos){
        if(
            e->pt() <= 23. || abs(e->eta()) >= 2.5 
            || abs(e->gsfTrack()->dxy(PV.position()))  >= 0.045 || abs(e->gsfTrack()->dz(PV.position())) >= 0.2
            || (abs(e->eta()) < 0.8 && e->electronID("POG_MVA_ID_Run2_NonTrig_Tight") <= 0.73)
            || (abs(e->eta()) < 1.479 && abs(e->eta()) >= 0.8 && e->electronID("POG_MVA_ID_Run2_NonTrig_Tight") <= 0.57)
            || (abs(e->eta()) >= 1.479  && e->electronID("POG_MVA_ID_Run2_NonTrig_Tight") <= 0.05)

            || tau->pt() <= 20 
            || abs(tau->eta()) >= 2.3
            || (tau->tauID("decayModeFindingNewDMs") <= 0.5 && tau->tauID("decayModeFinding") <= 0.5)
            || abs(tau->vertex().z() - PV.z()) < 0.2 
            //|| (tau->zImpact() <= 0.5 && tau->zImpact() >= -1.5)
            //|| tau->vertex().z() != PV.z()

            || deltaR(e->p4(), tau->p4()) <= 0.5
        ) return false;

        std::string path = "";
        for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
            if(triggerBits->accept(i)){
                path += names.triggerName(i);
            }
        }
        
        /*
        //std::string seed = "";
        for (CItAlgo algo = menu->gtAlgorithmMap().begin(); algo!=menu->gtAlgorithmMap().end(); ++algo) {
           // std::cout << "Name: " << (algo->second).algoName() << " Alias: " << (algo->second).algoAlias() << std::endl;
           //seed += (algo->second).algoName();
        }
        */
        if(mc){
            if (path.find("HLT_Ele22_eta2p1_WP85_Gsf_LooseIsoPFTau20_v1") != std::string::npos){
                for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
                    obj.unpackPathNames(names);
                    // for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << obj.filterLabels()[h]  << std::endl;
                    if( 
                        deltaR( obj.triggerObject().p4(), e->p4()) < 0.5  
                        && std::find(obj.filterLabels().begin(), obj.filterLabels().end(), "hltOverlapFilterIsoEle22WP85GsfLooseIsoPFTau20")!=obj.filterLabels().end() 
                    ){
                        for (pat::TriggerObjectStandAlone obj_ : *triggerObjects) {
                            obj_.unpackPathNames(names);
                            if(
                                deltaR(obj_.triggerObject().p4(), tau->p4()) < 0.5 && 
                                (std::find(obj_.filterLabels().begin(), obj_.filterLabels().end(), "hltL1sL1IsoEG20erTauJet20er")!=obj_.filterLabels().end() 
                                && std::find(obj_.filterLabels().begin(), obj_.filterLabels().end(), "hltOverlapFilterIsoEle22WP85GsfLooseIsoPFTau20")!=obj_.filterLabels().end() )
                            ) return true;
                        }
                        break;
                    }
                } 
            }
            if (path.find("HLT_Ele27_eta2p1_WP85_Gsf_v1") != std::string::npos  ){
                for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
                    obj.unpackPathNames(names);
                    if( 
                        deltaR( obj.triggerObject().p4(), e->p4()) < 0.5 
                        && std::find(obj.filterLabels().begin(), obj.filterLabels().end(), "hltEle27WP85GsfTrackIsoFilter")!=obj.filterLabels().end() 
                    ) return true;
                } 
            }
        }
        else {
            if (path.find("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v1") != std::string::npos){
                for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
                    obj.unpackPathNames(names);
                    // for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << obj.filterLabels()[h]  << std::endl;
                    if( 
                        deltaR( obj.triggerObject().p4(), e->p4()) < 0.5  
                        && std::find(obj.filterLabels().begin(), obj.filterLabels().end(), "hltSingleEle22WPLooseGsfTrackIsoFilter")!=obj.filterLabels().end() 
                        && std::find(obj.filterLabels().begin(), obj.filterLabels().end(), "hltOverlapFilterIsoEle22WPLooseGsfLooseIsoPFTau20")!=obj.filterLabels().end() 
                    ){
                        for (pat::TriggerObjectStandAlone obj_ : *triggerObjects) {
                            obj_.unpackPathNames(names);
                            if(
                                deltaR(obj_.triggerObject().p4(), tau->p4()) < 0.5 && 
                                (std::find(obj_.filterLabels().begin(), obj_.filterLabels().end(), "hltPFTau20TrackLooseIso")!=obj_.filterLabels().end() 
                                && std::find(obj_.filterLabels().begin(), obj_.filterLabels().end(), "hltOverlapFilterIsoEle22WPLooseGsfLooseIsoPFTau20")!=obj_.filterLabels().end() )
                            ) return true;
                        }
                        break;
                    }
                }
            }
            if (path.find("HLT_Ele32_eta2p1_WPTight_Gsf_v1") != std::string::npos  ){
                for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
                    obj.unpackPathNames(names);
                    if( 
                        deltaR( obj.triggerObject().p4(), e->p4()) < 0.5 
                        && std::find(obj.filterLabels().begin(), obj.filterLabels().end(), "hltEle32WPTightGsfTrackIsoFilter")!=obj.filterLabels().end() 
                        && e->pt() > 33.
                    ) return true;
                } 
            }
        }
    }



    /*
    if(abs(e->gsfTrack()->dxy(PV.position()))  >= 0.045 || abs(e->gsfTrack()->dz(PV.position())) >= 0.2 )
        continue;
    if (!e->electronID("POG_MVA_ID_Run2_NonTrig_Tight"))
        continue;
    if(utilities::relIso(*e, 0.5) >= 0.1)       
        continue;
    if (e->pt() <= 23. || abs(e->eta()) >= 2.5)
        continue;

    if( tau->tauID("decayModeFinding") <= 0.5 || 
        tau->tauID("decayModeFindingNewDMs") <= 0.5 ||  
        tau->tauID("againstElectronTightMVA5") <= 0.5 ||
        tau->tauID("againstMuonLoose3") <= 0.5 ||
        (tau->vertex().z() + 130./tan(tau->theta()) <= 0.5 && tau->vertex().z() + 130./tan(tau->theta()) >= -1.5) || //  tau.zImpact() in heppy
        tau->vertex() != PV.position()) 
        continue;
    if(tau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") >= 1.5)
        continue;
    if (tau->pt() <= 20 and abs(tau->eta()) >= 2.3)
        continue;

    //V E T O
    for (const pat::Electron &vetoel : *electrons) {
        if (!vetoel.gsfTrack().isNull() && abs(vetoel.px() - e->px()) +  abs(vetoel.py() - e->py()) + abs(vetoel.pz() - e->pz()) > 0.01){
            if (vetoel.pt() > 10                                     &&
                abs(vetoel.eta()) < 2.5                              &&
                abs(vetoel.gsfTrack()->dxy(PV.position())) < 0.045                          &&
                abs(vetoel.gsfTrack()->dz(PV.position())) < 0.2                             &&
           //  FIXME    vetoel.cutBasedId('POG_PHYS14_25ns_v1_Veto')       &&
                utilities::relIso(vetoel, 0.5) < 0.3 )
                    continue;
        }
    }

    for (const pat::Muon &vetomu : *muons) {
        if (!vetomu.innerTrack().isNull() && abs(vetomu.px() - e->px()) +  abs(vetomu.py() - e->py()) + abs(vetomu.pz() - e->pz()) > 0.01){
            if (vetomu.pt() > 10           &&
                abs(vetomu.eta()) < 2.4    &&
                utilities::heppymuonID(vetomu, "POG_ID_Medium")   &&
                abs(vetomu.innerTrack()->dxy( PV.position())) < 0.045  &&
                abs(vetomu.innerTrack()->dz( PV.position())) < 0.2   &&
                utilities::relIso(vetomu, 0.5) < 0.3)
                continue;
        }
    }
    */
    return false;   

}
bool PairBaselineSelection::mumu(const pat::Muon* mu, const pat::Muon* mu_){
    return false;   
}
bool PairBaselineSelection::ee(const pat::Electron* e, const pat::Electron* e_){
    return false;
}
bool PairBaselineSelection::emu(const reco::Vertex &PV, edm::Handle<edm::TriggerResults>& triggerBits, edm::Handle<pat::TriggerObjectStandAloneCollection>& triggerObjects, edm::Handle<pat::PackedTriggerPrescales>& triggerPrescales,const edm::TriggerNames &names, const pat::Electron* e,const pat::Muon* mu){
    if(sample.find("spring15") != std::string::npos){

        if(
                e->pt() <= 13. || abs(e->eta() >=2.5) 
                || abs(e->gsfTrack()->dxy(PV.position()))  >= 0.045 || abs(e->gsfTrack()->dz(PV.position()) >= 0.2)
                || e->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 1
                || e->passConversionVeto()

                || mu->pt() <= 10. || abs(mu->eta()) >= 2.4
                || abs(mu->innerTrack()->dxy( PV.position() ))  >= 0.045 || abs(mu->innerTrack()->dz(PV.position())) >= 0.2
                || !utilities::heppymuonID(*mu, "POG_ID_Medium") 

                || deltaR(e->p4(), mu->p4()) <= 0.3

                ) return false;

        std::string path = "";
        for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
            if(triggerBits->accept(i)){
                path += names.triggerName(i);
                //std::cout <<  names.triggerName(i) << std::endl;
            }
        }
        if(mc){
            if (path.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1") != std::string::npos){
                for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
                    obj.unpackPathNames(names);
                    // for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << obj.filterLabels()[h]  << std::endl;
                    if( 
                        deltaR( obj.triggerObject().p4(), mu->p4()) < 0.5 
                        && std::find(obj.filterLabels().begin(), obj.filterLabels().end(), "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23")!=obj.filterLabels().end() 
                    ){
                        for (pat::TriggerObjectStandAlone obj_ : *triggerObjects) {
                            obj_.unpackPathNames(names);
                            if(
                                deltaR(obj_.triggerObject().p4(), e->p4()) < 0.5 
                                && std::find(obj_.filterLabels().begin(), obj_.filterLabels().end(), "hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter")!=obj_.filterLabels().end() 
                            ) return true;
                        }
                        break;
                    }
                }
            }
            if (path.find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v1") != std::string::npos){
                for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
                    obj.unpackPathNames(names);
                    // for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << obj.filterLabels()[h]  << std::endl;
                    if( 
                        deltaR( obj.triggerObject().p4(), mu->p4()) < 0.5 
                        && std::find(obj.filterLabels().begin(), obj.filterLabels().end(), "hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8")!=obj.filterLabels().end() 
                    ){
                        for (pat::TriggerObjectStandAlone obj_ : *triggerObjects) {
                            obj_.unpackPathNames(names);
                            if(
                                deltaR(obj_.triggerObject().p4(), e->p4()) < 0.5 
                                && std::find(obj_.filterLabels().begin(), obj_.filterLabels().end(), "hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter")!=obj_.filterLabels().end() 
                            ) return true;
                        }
                        break;
                    }
                }
            }
        }
    }

    if(sample.find("phys14") != std::string::npos){
        if(
            e->pt() <= 13. || abs(e->eta() >=2.5) 
            || abs(e->gsfTrack()->dxy(PV.position()))  >= 0.045 || abs(e->gsfTrack()->dz(PV.position()) >= 0.2)

            || mu->pt() <= 9. || abs(mu->eta()) >= 2.4
            || abs(mu->innerTrack()->dxy( PV.position() ))  >= 0.045 || abs(mu->innerTrack()->dz(PV.position())) >= 0.2
            || !utilities::heppymuonID(*mu, "POG_ID_Medium") 

            || deltaR(e->p4(), mu->p4()) <= 0.3

        ) return false;

        std::string path = "";
        for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
            if(triggerBits->accept(i)){
                path += names.triggerName(i);
                //std::cout <<  names.triggerName(i) << std::endl;
            }
        }
        if(mc){
            if (path.find("HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v1") != std::string::npos){
                for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
                    obj.unpackPathNames(names);
                    // for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << obj.filterLabels()[h]  << std::endl;
                    if( 
                        deltaR( obj.triggerObject().p4(), mu->p4()) < 0.5 
                        && std::find(obj.filterLabels().begin(), obj.filterLabels().end(), "hltL1Mu12EG7L3IsoMuFiltered23")!=obj.filterLabels().end() 
                    ){
                        for (pat::TriggerObjectStandAlone obj_ : *triggerObjects) {
                            obj_.unpackPathNames(names);
                            if(
                                deltaR(obj_.triggerObject().p4(), e->p4()) < 0.5 
                                && std::find(obj_.filterLabels().begin(), obj_.filterLabels().end(), "hltMu23Ele12GsfTrackIsoLegEle12GsfCaloIdTrackIdIsoMediumWPFilter")!=obj_.filterLabels().end() 
                            ) return true;
                        }
                        break;
                    }
                }
            }
            if (path.find("HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v1") != std::string::npos){
                for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
                    obj.unpackPathNames(names);
                    // for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << obj.filterLabels()[h]  << std::endl;
                    if( 
                        deltaR( obj.triggerObject().p4(), mu->p4()) < 0.5 
                        && std::find(obj.filterLabels().begin(), obj.filterLabels().end(), "hltL1sL1Mu5EG20ORL1Mu5IsoEG18L3IsoFiltered8")!=obj.filterLabels().end() 
                    ){
                        for (pat::TriggerObjectStandAlone obj_ : *triggerObjects) {
                            obj_.unpackPathNames(names);
                            if(
                                deltaR(obj_.triggerObject().p4(), e->p4()) < 0.5 
                                && std::find(obj_.filterLabels().begin(), obj_.filterLabels().end(), "hltMu8Ele23GsfTrackIsoLegEle23GsfCaloIdTrackIdIsoMediumWPFilter")!=obj_.filterLabels().end() 
                            ) return true;
                        }
                        break;
                    }
                }
            }
        }
    }

    /*
    if(
        (abs(el->gsfTrack()->dxy(PV.position()))  >= 0.045 || abs(el->gsfTrack()->dz(PV.position())) >= 0.2)  ||
        !el->electronID("POG_MVA_ID_Run2_NonTrig_Tight") ||
        utilities::relIso(*el, 0.5) >= 0.15 ||
        el->pt() <= 13 || abs(el->eta()) >= 2.5 ||
        (abs(mu->innerTrack()->dxy( PV.position() ))  >= 0.045 || abs(mu->innerTrack()->dz(PV.position())) >= 0.2 ) ||
        !utilities::heppymuonID(*mu, "POG_ID_Medium") ||
        utilities::relIso(*mu, 0.5) >= 0.15 ||
        (mu->pt() <= 9 || abs(mu->eta()) >= 2.4)
    ) return true;
    */
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
