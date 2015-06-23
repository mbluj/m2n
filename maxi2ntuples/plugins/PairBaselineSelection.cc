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

#include "DataFormats/Math/interface/deltaR.h"
#include "m2n/maxi2ntuples/interface/utilities.h"

#include <vector>
#include <string>
#include <math.h>
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
      
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<pat::CompositeCandidateCollection> PairToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
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
    electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons")))
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

    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByToken(muonToken_, muons);

    edm::Handle<pat::ElectronCollection> electrons;
    iEvent.getByToken(electronToken_, electrons);

    edm::Handle<pat::CompositeCandidateCollection> leptonPair;
    iEvent.getByToken(PairToken_, leptonPair);
    if (!leptonPair.isValid()) return;


    std::unique_ptr<pat::CompositeCandidateCollection> selectedPair(new pat::CompositeCandidateCollection());
    

    for (const pat::CompositeCandidate &lP : *leptonPair){

        const reco::Candidate * l1, *l2; 
        l1 = lP.daughter(0); l2 = lP.daughter(1);
        if (l1->isMuon()){
            if (l2->isMuon()){
                //std::cout << "muon-muon channel" <<std::endl;
//                const pat::Muon *mu  = dynamic_cast<const pat::Muon*>(l1);
//                const pat::Muon *mu_  = dynamic_cast<const pat::Muon*>(l2);
            }
            else if(l2->isElectron()){
            /*
                //std::cout << "muon-electron channel" <<std::endl;
                const pat::Muon *mu  = dynamic_cast<const pat::Muon*>(l1);
                const pat::Electron* el = dynamic_cast<const pat::Electron*>(l2);
                if(abs(el->gsfTrack()->dxy(PV.position()))  >= 0.045 || abs(el->gsfTrack()->dz(PV.position())) >= 0.2 )
                    continue;
                if (!el->electronID("POG_MVA_ID_Run2_NonTrig_Tight"))
                    continue;
                if(utilities::relIso(*el, 0.5) >= 0.15)       
                    continue;
                if (el->pt() <= 13 || abs(el->eta()) >= 2.5)
                    continue;
                
                if(abs(mu->innerTrack()->dxy( PV.position() ))  >= 0.045 || abs(mu->innerTrack()->dz(PV.position())) >= 0.2 )
                    continue;
                if (!utilities::heppymuonID(*mu, "POG_ID_Medium")){
                    continue;
                }
                if(utilities::relIso(*mu, 0.5) >= 0.15)       
                    continue;
                if (mu->pt() <= 9 || abs(mu->eta()) >= 2.4)
                    continue;
            */
            }
            else {
                //std::cout << "muon-tau channell" <<std::endl;
               const pat::Muon* mu = dynamic_cast<const pat::Muon*>(l1->masterClone().get());
               const pat::Tau* tau = dynamic_cast<const pat::Tau*>(l2->masterClone().get());
                
                if(abs(mu->innerTrack()->dxy( PV.position()) )  >= 0.045 || abs(mu->innerTrack()->dz(PV.position())) >= 0.2 )
                    continue;
                if (!utilities::heppymuonID(*mu, "POG_ID_Medium"))
                    continue;
                if(utilities::relIso(*mu, 0.5) >= 0.1)       
                    continue;
                if (mu->pt() <= 18 && abs(mu->eta()) >= 2.1)
                    continue;
                /*
                if( tau->tauID("decayModeFinding") <= 0.5 || 
                    tau->tauID("decayModeFindingNewDMs") <= 0.5 ||  
                    tau->tauID("againstElectronVLooseMVA5") <= 0.5 ||
                    tau->tauID("againstMuonTight3") <= 0.5 ||
                    tau->vertex() != PV.position()) 
                    continue;
               */
                if (tau->pt() <= 20 and abs(tau->eta()) >= 2.3)
                    continue;
                // VETO
                // Second lepton veto
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
                
            }
        }/*
        else if(l1->isElectron()){
            if (l2->isMuon()){
                //std::cout << "muon-electron channel" <<std::endl;
                const pat::Electron* el = dynamic_cast<const pat::Electron*>(l1);
                const pat::Muon *mu  = dynamic_cast<const pat::Muon*>(l2);
                if(abs(el->gsfTrack()->dxy(PV.position()))  >= 0.045 || abs(el->gsfTrack()->dz(PV.position())) >= 0.2 )
                    continue;
                if (!el->electronID("POG_MVA_ID_Run2_NonTrig_Tight"))
                    continue;
                if(utilities::relIso(*el, 0.5) >= 0.15)       
                    continue;
                if (el->pt() <= 13 || abs(el->eta()) >= 2.5)
                    continue;
                
                if(abs(mu->innerTrack()->dxy( PV.position() ))  >= 0.045 || abs(mu->innerTrack()->dz(PV.position())) >= 0.2 )
                    continue;
                if (!utilities::heppymuonID(*mu, "POG_ID_Medium")){
                    continue;
                }
                if(utilities::relIso(*mu, 0.5) >= 0.15)       
                    continue;
                if (mu->pt() <= 9 || abs(mu->eta()) >= 2.4)
                    continue;
            }
            else if(l2->isElectron()){
                //std::cout << "electron-electron channel" <<std::endl;
  //              const pat::Electron* el = dynamic_cast<const pat::Electron*>(l1);
  //              const pat::Electron* el_ = dynamic_cast<const pat::Electron*>(l2);
            }
            else {
                //std::cout << "electron-tau channel" <<std::endl;
                const pat::Electron* el = dynamic_cast<const pat::Electron*>(l1);
                const pat::Tau* tau = dynamic_cast<const pat::Tau*>(l2);

                if(abs(el->gsfTrack()->dxy(PV.position()))  >= 0.045 || abs(el->gsfTrack()->dz(PV.position())) >= 0.2 )
                    continue;
                if (!el->electronID("POG_MVA_ID_Run2_NonTrig_Tight"))
                    continue;
                if(utilities::relIso(*el, 0.5) >= 0.1)       
                    continue;
                if (el->pt() <= 23. || abs(el->eta()) >= 2.5)
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
                    if (!vetoel.gsfTrack().isNull() && abs(vetoel.px() - el->px()) +  abs(vetoel.py() - el->py()) + abs(vetoel.pz() - el->pz()) > 0.01){
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
                    if (!vetomu.innerTrack().isNull() && abs(vetomu.px() - el->px()) +  abs(vetomu.py() - el->py()) + abs(vetomu.pz() - el->pz()) > 0.01){
                        if (vetomu.pt() > 10           &&
                            abs(vetomu.eta()) < 2.4    &&
                            utilities::heppymuonID(vetomu, "POG_ID_Medium")   &&
                            abs(vetomu.innerTrack()->dxy( PV.position())) < 0.045  &&
                            abs(vetomu.innerTrack()->dz( PV.position())) < 0.2   &&
                            utilities::relIso(vetomu, 0.5) < 0.3)
                            continue;
                    }
                }

            }
        
        }
        else {
            if (l2->isMuon()){
                //std::cout << "muon-tau channel" <<std::endl;
                const pat::Tau* tau = dynamic_cast<const pat::Tau*>(l1);
                const pat::Muon *mu  = dynamic_cast<const pat::Muon*>(l2);
                if(abs(mu->innerTrack()->dxy( PV.position()) )  >= 0.045 || abs(mu->innerTrack()->dz(PV.position())) >= 0.2 )
                    continue;
                if (!utilities::heppymuonID(*mu, "POG_ID_Medium"))
                    continue;
                if(utilities::relIso(*mu, 0.5) >= 0.1)       
                    continue;
                if (mu->pt() <= 18 && abs(mu->eta()) >= 2.1)
                    continue;

                if( tau->tauID("decayModeFinding") <= 0.5 || 
                    tau->tauID("decayModeFindingNewDMs") <= 0.5 ||  
                    tau->tauID("againstElectronVLooseMVA5") <= 0.5 ||
                    tau->tauID("againstMuonTight3") <= 0.5 ||
                    tau->vertex() != PV.position()) 
                    continue;
                if(tau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") >= 1.5)
                    continue;
                if (tau->pt() <= 20 and abs(tau->eta()) >= 2.3)
                    continue;
                // VETO
                // Second lepton veto
                for (const pat::Muon &vetomu : *muons) {
                    if (abs(vetomu.px() - mu->px()) +  abs(vetomu.py() - mu->py()) + abs(vetomu.pz() - mu->pz()) > 0.01 && !vetomu.innerTrack().isNull() ){
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
            }
            else if(l2->isElectron()){
                //std::cout << "electron-tau channel" <<std::endl;
                const pat::Tau* tau = dynamic_cast<const pat::Tau*>(l1);
                const pat::Electron* el = dynamic_cast<const pat::Electron*>(l2);

                if(abs(el->gsfTrack()->dxy(PV.position()))  >= 0.045 || abs(el->gsfTrack()->dz(PV.position())) >= 0.2 )
                    continue;
                if (!el->electronID("POG_MVA_ID_Run2_NonTrig_Tight"))
                    continue;
                if(utilities::relIso(*el, 0.5) >= 0.1)       
                    continue;
                if (el->pt() <= 23. || abs(el->eta()) >= 2.5)
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
                    if (!vetoel.gsfTrack().isNull() && abs(vetoel.px() - el->px()) +  abs(vetoel.py() - el->py()) + abs(vetoel.pz() - el->pz()) > 0.01){
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
                    if (!vetomu.innerTrack().isNull() && abs(vetomu.px() - el->px()) +  abs(vetomu.py() - el->py()) + abs(vetomu.pz() - el->pz()) > 0.01){
                        if (vetomu.pt() > 10           &&
                            abs(vetomu.eta()) < 2.4    &&
                            utilities::heppymuonID(vetomu, "POG_ID_Medium")   &&
                            abs(vetomu.innerTrack()->dxy( PV.position())) < 0.045  &&
                            abs(vetomu.innerTrack()->dz( PV.position())) < 0.2   &&
                            utilities::relIso(vetomu, 0.5) < 0.3)
                            continue;
                    }
                }
            }
            else {
                //std::cout << "tau-tau channel" <<std::endl;
                const pat::Tau* tau = dynamic_cast<const pat::Tau*>(l1);
                const pat::Tau* tau_ = dynamic_cast<const pat::Tau*>(l2);
                //leading and trailing taus has the same cuts so far => no if needed
//                if (tau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") < tau_->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") ){
                    if( tau->tauID("decayModeFinding") <= 0.5 || 
                        tau->tauID("decayModeFindingNewDMs") <= 0.5 ||  
                        tau->tauID("againstElectronVLooseMVA5") <= 0.5 ||
                        tau->tauID("againstMuonLoose3") <= 0.5 ||
                        tau->vertex() != PV.position()) 
                        continue;
                    if(tau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") >= 1.)
                        continue;
                    if (tau->pt() <= 45. and abs(tau->eta()) >= 2.1)
                        continue;
                
                    if( tau_->tauID("decayModeFinding") <= 0.5 || 
                        tau_->tauID("decayModeFindingNewDMs") <= 0.5 ||  
                        tau_->tauID("againstElectronVLooseMVA5") <= 0.5 ||
                        tau_->tauID("againstMuonLoose3") <= 0.5 ||
                        tau_->vertex() != PV.position()) 
                        continue;
                    if(tau_->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") >= 1.)
                        continue;
                    if (tau_->pt() <= 45. and abs(tau_->eta()) >= 2.1)
                        continue;
//                }
//                else{
//                }
            }
        
        }*/




          //  std::cout << "passsed" <<std::endl;
            


        selectedPair->push_back(lP);
    
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
