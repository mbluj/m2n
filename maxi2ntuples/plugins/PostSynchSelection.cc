// -*- C++ -*-
//
// Package:    m2n/PostSynchSelection
// Class:      PostSynchSelection
// 
/**\class PostSynchSelection PostSynchSelection.cc m2n/PostSynchSelection/plugins/PostSynchSelection.cc

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

class PostSynchSelection : public edm::EDProducer {
public:
  explicit PostSynchSelection(const edm::ParameterSet&);
  ~PostSynchSelection(){;}
  
private:
  virtual void beginJob(){;}
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob(){;}
  
  bool checkDiMuonVeto(const pat::MuonCollection &);
  bool checkThirdLeptonVeto(const pat::MuonCollection &, const pat::ElectronCollection &);
  
  float tautau(const reco::Vertex &,  const pat::Tau*, const pat::Tau*);
  float mutau(const reco::Vertex &,  const pat::Muon*, const pat::Tau*);
  float etau(const reco::Vertex &, const pat::Electron*, const pat::Tau*);
  float mumu(const pat::Muon*, const pat::Muon*);
  float ee(const pat::Electron*, const pat::Electron*);
  float emu(const reco::Vertex &, const pat::Electron*,const pat::Muon* );
  
  // ----------member data ---------------------------
  edm::EDGetTokenT<pat::CompositeCandidateCollection> PairToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  const bool mc;
  const std::string sample;
  edm::EDGetToken electronToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
};
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
PostSynchSelection::PostSynchSelection(const edm::ParameterSet& iConfig):
    PairToken_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("pairs"))),
    vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
    triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
    triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
    triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
    mc(iConfig.getParameter<bool>("mc")),
    sample(iConfig.getParameter<std::string>("sample")),
    electronToken_(mayConsume<edm::View<reco::GsfElectron> > (iConfig.getParameter<edm::InputTag>("electrons"))),
    eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))){

   produces<pat::CompositeCandidateCollection>();
  
}
/////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////// 
bool PostSynchSelection::checkDiMuonVeto(const pat::MuonCollection &muons){

  if(muons.size()>1){
    for(auto aMuon1: muons){
      bool passDiMuonSelection1 = aMuon1.pt()>15 && aMuon1.isGlobalMuon() && aMuon1.isTrackerMuon() && aMuon1.isPFMuon();
      if(passDiMuonSelection1){
	for(auto aMuon2: muons){
	  bool passDiMuonSelection2 = aMuon2.pt()>15 && aMuon2.isGlobalMuon() && aMuon2.isTrackerMuon() && aMuon2.isPFMuon();
	  if(passDiMuonSelection2 && deltaR(aMuon1.p4(), aMuon2.p4())>0.15) return false; 
	}
      }
    }    
  }
  return true;
}
/////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////// 
bool PostSynchSelection::checkThirdLeptonVeto(const pat::MuonCollection &muons, const pat::ElectronCollection &electrons){

  unsigned int nMuons = 0;
  for(auto aMuon: muons){
    if(aMuon.isMediumMuon()) ++nMuons;
  }
  unsigned int nElectrons = electrons.size();

  return !(nMuons>1 && nElectrons>0);
  
}
/////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////// 
void PostSynchSelection::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

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

    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByToken(muonToken_, muons);

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

    bool veto = true;
    
    pat::MuonCollection vetomuons;
    for (const pat::Muon &mu : *muons) {
        float iso = (mu.pfIsolationR03().sumChargedHadronPt + 
            std::max(mu.pfIsolationR03().sumNeutralHadronEt + mu.pfIsolationR03().sumPhotonEt - 0.5 * mu.pfIsolationR03().sumPUPt, 0.0)) / mu.pt();
        if(
            mu.pt() >10
            && fabs(mu.eta()) < 2.4
            && fabs(mu.muonBestTrack()->dz(PV.position())) < 0.2
            && fabs(mu.muonBestTrack()->dxy( PV.position()) ) < 0.045
            && iso < 0.3
        ) vetomuons.push_back(mu);
    }
    
    pat::ElectronCollection vetoelectrons;
    for (size_t i = 0; i < electrons->size(); ++i){
         const auto el = electrons->ptrAt(i);
         pat::Electron elec(el);
         float iso = (el->pfIsolationVariables().sumChargedHadronPt + std::max(
                 el->pfIsolationVariables().sumNeutralHadronEt + el->pfIsolationVariables().sumPhotonEt - 
                 0.5 * el->pfIsolationVariables().sumPUPt, 0.0)) / el->pt();
         if(
             el->pt() >10
             && fabs(el->eta())  < 2.5
             && fabs(el->gsfTrack()->dxy(PV.position())) < 0.045
             && fabs(el->gsfTrack()->dz(PV.position())) < 0.2
             && elec.passConversionVeto()
             && el->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)  <=1
             && iso < 0.3
             && (*tight_id_decisions)[el]
         ) vetoelectrons.push_back(*el);
    }


    veto &= checkDiMuonVeto(vetomuons);
    veto &= checkThirdLeptonVeto(vetomuons, vetoelectrons);
    
    for (const pat::CompositeCandidate &lP : *leptonPair){

        const reco::Candidate * l1, *l2; 
        l1 = lP.daughter(0); l2 = lP.daughter(1);

        float pass = false; 

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
                const pat::Electron *e = dynamic_cast<const pat::Electron*>(l2->masterClone().get());
                pass =  emu(PV, e, mu);
            }
            else {
                //std::cout << "muon-tau channell" <<std::endl;
                const pat::Muon *mu  = dynamic_cast<const pat::Muon*>(l1->masterClone().get());
                const pat::Tau *tau  = dynamic_cast<const pat::Tau*>(l2->masterClone().get());
                pass = mutau(PV,  mu,  tau);
            }
        }
        else if(l1->isElectron()){
            if (l2->isMuon()){
                //std::cout << "muon-electron channel" <<std::endl;
                const pat::Electron *e = dynamic_cast<const pat::Electron*>(l1->masterClone().get());
                const pat::Muon *mu = dynamic_cast<const pat::Muon*>(l2->masterClone().get());
                pass =  emu(PV, e, mu);
            }
            else if(l2->isElectron()){
                //std::cout << "electron-electron channel" <<std::endl;
                const pat::Electron *el = dynamic_cast<const pat::Electron*>(l1->masterClone().get());
                const pat::Electron *el_ = dynamic_cast<const pat::Electron*>(l2->masterClone().get());
                pass = ee(el,el_);
            }
            else {
                //std::cout << "electron-tau channel" <<std::endl;
                const pat::Electron *e = dynamic_cast<const pat::Electron*>(l1->masterClone().get());
                const pat::Tau *tau  = dynamic_cast<const pat::Tau*>(l2->masterClone().get());
                pass = etau(PV, e,tau);
            }
        }
        else {
            if (l2->isMuon()){
                //std::cout << "muon-tau channel" <<std::endl;
                const pat::Tau *tau  = dynamic_cast<const pat::Tau*>(l1->masterClone().get());
                const pat::Muon *mu = dynamic_cast<const pat::Muon*>(l2->masterClone().get());
                pass = mutau(PV, mu,  tau);
            }
            else if(l2->isElectron()){
                //std::cout << "electron-tau channel" <<std::endl;
                const pat::Tau *tau  = dynamic_cast<const pat::Tau*>(l1->masterClone().get());
                const pat::Electron *e = dynamic_cast<const pat::Electron*>(l2->masterClone().get());
                pass = etau(PV, e,tau);
            }
            else {
                //std::cout << "tau-tau channel" <<std::endl;
                const pat::Tau *tau  = dynamic_cast<const pat::Tau*>(l1->masterClone().get());
                const pat::Tau *tau_  = dynamic_cast<const pat::Tau*>(l2->masterClone().get());
                pass = tautau(PV,tau, tau_);
            }
        }
        bool passed = veto && pass;
        pat::CompositeCandidate pair(lP);
        pair.addUserFloat("PostSynchSelection", (float)passed);
        selectedPair->push_back(pair);
    }

    iEvent.put(std::move(selectedPair));
}
/////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////// 
float PostSynchSelection::mutau(const reco::Vertex &PV,  const pat::Muon* mu, const pat::Tau* tau){

    if(
        tau->tauID("againstElectronVLooseMVA5") <= 0.5
        || tau->tauID("againstMuonTight3") <= 0.5
        || tau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits") >= 0.5
    ) return false;

    return true;

}
/////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////// 
float PostSynchSelection::tautau(const reco::Vertex &PV, const pat::Tau* tau, const pat::Tau* tau_){ return true; }
float PostSynchSelection::etau(const reco::Vertex &PV,  const pat::Electron* e, const pat::Tau* tau){ return false; }
float PostSynchSelection::mumu(const pat::Muon* mu, const pat::Muon* mu_){ return false;}
float PostSynchSelection::ee(const pat::Electron* e, const pat::Electron* e_){return false;}
float PostSynchSelection::emu(const reco::Vertex &PV,  const pat::Electron* e,const pat::Muon* mu){return false;}

//define this as a plug-in
DEFINE_FWK_MODULE(PostSynchSelection);
