//
// Package:    m2n/maxi2ntuples
// Class:      ntuple
// 
/**\class ntuple ntuple.cc m2n/maxi2ntuples/plugins/ntuple.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Wed, 17 Dec 2014 09:06:30 GMT
//
//


// system include files
#include <memory>
#include <vector>
#include <string>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include <DataFormats/Candidate/interface/ShallowCloneCandidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include <DataFormats/METReco/interface/CommonMETData.h>
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "m2n/maxi2ntuples/interface/utilities.h"

#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TMVA/Reader.h"
#include "TLorentzVector.h"
#include "RConfig.h" 
#include "TObject.h"
#include "TH1F.h"

//#include "clasadict.h"

//
// class declaration
//

//gSystem->Load("libPhysics.so");
class ntuple : public edm::EDAnalyzer {
   public:
      explicit ntuple(const edm::ParameterSet&);
      ~ntuple();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      typedef std::vector<PileupSummaryInfo> PileupSummaryInfoCollection;


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;


  // ----------member data ---------------------------
    edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
    edm::EDGetTokenT<pat::MuonCollection> muonToken_;
    edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
    edm::EDGetTokenT<pat::TauCollection> tauToken_;
    edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
    edm::EDGetTokenT<pat::JetCollection> jetToken_;
    edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
    edm::EDGetTokenT<pat::METCollection> metToken_;
    edm::EDGetTokenT<pat::CompositeCandidateCollection> pairToken_;
    edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
    edm::EDGetTokenT<reco::GenParticleCollection> prunedGenToken_;
    edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGenToken_;
    edm::EDGetTokenT<LHEEventProduct> lheprodToken_;
    edm::EDGetTokenT<PileupSummaryInfoCollection> PileupSummaryInfoToken_;
    const bool mc;



//    TNtuple* evt;
    TTree * event; 
    TTree * isZ;
    TTree * trigger;
    TNtuple* mu;
    TNtuple* mu_;
    TNtuple* e;
    TNtuple* e_;
    TNtuple* tau;
    TNtuple* tau_;
    TNtuple* pair;
    TNtuple* met;
    TNtuple* leadingjet;
    TNtuple* trailingjet;
    TNtuple* jetpair;
    TH1F* events;


    int  run, lumi, evt, npv, nup, paircount, npu; 
    bool isZtt, isZmt, isZet, isZee, isZmm, isZem, isZEE, isZMM, isZLL;

    std::vector<std::string> hltmatch;
    std::vector<std::string> hltpaths;




/*
        hdijetphi, visjeteta; 
*/




     
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


// -*- C++ -*-
ntuple::ntuple(const edm::ParameterSet& iConfig):
    vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
    electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
    tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
    photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
    jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
    fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
    metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
    pairToken_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("pairs"))),
    triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
    triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
    triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
    prunedGenToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
    packedGenToken_(consumes<pat::PackedGenParticleCollection>(iConfig.getParameter<edm::InputTag>("packedGenParticles"))),
    lheprodToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheprod"))),
    PileupSummaryInfoToken_(consumes<PileupSummaryInfoCollection>(iConfig.getParameter<edm::InputTag>("pileupinfo"))),
    mc(iConfig.getParameter<bool>("mc"))
{
   //now do what ever initialization is needed
}


ntuple::~ntuple()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ntuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
 

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);
    if (vertices->empty()) return; // skip the event if no PV found
    const reco::Vertex &PV = vertices->front();

    edm::Handle<pat::METCollection> mets;
    iEvent.getByToken(metToken_, mets);
 //   const pat::MET &smet = mets->front();


    edm::Handle<pat::CompositeCandidateCollection> pairs;
    iEvent.getByToken(pairToken_, pairs);

    edm::Handle<edm::TriggerResults> triggerBits;
    iEvent.getByToken(triggerBits_, triggerBits);
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    iEvent.getByToken(triggerObjects_, triggerObjects);
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    iEvent.getByToken(triggerPrescales_, triggerPrescales);
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

    edm::Handle<pat::JetCollection> jets; 
    iEvent.getByToken(jetToken_, jets);

    if(mc){

        edm::Handle<PileupSummaryInfoCollection> genPileUpInfos;
        iEvent.getByToken(PileupSummaryInfoToken_, genPileUpInfos);

        npu = -1;
        for(const PileupSummaryInfo &pusi : *genPileUpInfos){
            int bx = pusi.getBunchCrossing();
            int nPU = pusi.getPU_NumInteractions();
            if(bx == 0)
                npu = nPU;
        }
        try{
            edm::Handle<LHEEventProduct> evnt;
            iEvent.getByToken(lheprodToken_, evnt);
            const lhef::HEPEUP hepeup_ = evnt->hepeup();
            nup =  hepeup_.NUP;
        }
        catch(const std::exception& e){
            std::cout << "No lheprod collection"; 
        }

    }

    evt = iEvent.id().event(); 
    paircount = pairs->size();
    run = iEvent.id().run();
    lumi = iEvent.luminosityBlock();
    npv =  vertices->size();
    event->Fill();
    events->Fill(1.,1.);
    if(mc){
        edm::Handle<GenEventInfoProduct> genEvt;
        iEvent.getByLabel("generator",genEvt);
       // event weight
        double weightevt=genEvt->weight(); 
        events->Fill(2., weightevt);
    }

    isZtt = isZmt = isZet = isZee = isZmm = isZem = isZEE = isZMM = isZLL = false;
    if(mc){
        edm::Handle<reco::GenParticleCollection> pruned; 
        iEvent.getByToken(prunedGenToken_, pruned);
    //    Handle<edm::View<pat::PackedGenParticle> > packed;
    //    iEvent.getByToken(packedGenToken_,packed);
        edm::Handle<pat::PackedGenParticleCollection> packed; 
        iEvent.getByToken(packedGenToken_, packed);
        std::vector<reco::GenParticle> gentaus; 
        std::vector<reco::GenParticle> gentauleps; 
        std::vector<reco::GenParticle> genleps; 
    //    for(size_t j=0; j<packed->size();j++){
        for (const reco::GenParticle &p : *pruned){
            
            if(abs(p.pdgId()) == 15)
                gentaus.push_back(p);
            if(abs(p.pdgId()) == 11 || abs(p.pdgId()) == 13){
                if(abs(p.mother(0)->pdgId()) == 15)
                    gentauleps.push_back(p);
                else
                    genleps.push_back(p);
            }
                    
        }
        if (gentaus.size() + gentauleps.size() == 2){
        
            if (gentaus.size() == 2)
                isZtt = true;
            else if(gentaus.size() == 1){
                if(abs(gentauleps[0].pdgId()) == 11)
                    isZet = true;
                if(abs(gentauleps[0].pdgId()) == 13)
                    isZmt = true;
            }
            else if(gentaus.size() == 0){
                if(abs(gentauleps[0].pdgId()) == 11 && abs(gentauleps[1].pdgId()) == 11)
                    isZee = true;
                else if(abs(gentauleps[0].pdgId()) == 13 && abs(gentauleps[1].pdgId()) == 13)
                    isZmm = true;
                else
                    isZem = true;
            }
        
        }
        else if(genleps.size() == 2){
            isZLL = true;
            if(abs(genleps[0].pdgId()) == 11 && abs(genleps[1].pdgId()) == 11)
                isZEE = true;
            else if(abs(genleps[0].pdgId()) == 13 && abs(genleps[1].pdgId()) == 13)
                isZMM = true;
        }
    }

    const reco::Candidate *mety = (pairs->front()).daughter(2);

    const pat::CompositeCandidate& lP = pairs->front();
    float pairarr[] = {lP.userFloat("SVfitMass"),(float)((pairs->front()).daughter(0)->charge() * (pairs->front()).daughter(1)->charge()), 
        (float)((pairs->front()).daughter(0)->p4() + (pairs->front()).daughter(1)->p4() + mety->p4()).pt(), 
        (float)((pairs->front()).daughter(0)->p4() + (pairs->front()).daughter(1)->p4()).pt(), 
        (float)((pairs->front()).daughter(0)->p4() + (pairs->front()).daughter(1)->p4()).mass() };
    pair->Fill(pairarr);


    const reco::Candidate * l1 = pairs->front().daughter(0);;
    if(l1->isMuon()){
        const pat::Muon *muon  = dynamic_cast<const pat::Muon*>((pairs->front()).daughter(0)->masterClone().get());
        float muarr[] = {(float)muon->pt(), (float)muon->eta(),(float) muon->phi(),(float) muon->mass(), (float)muon->charge(),(float)(muon->innerTrack()->dxy( PV.position())),
               (float)(muon->innerTrack()->dz(PV.position())), (float)sqrt(pow((muon->p4()).pt() + (mety->p4()).pt(),2) - pow((muon->p4() + mety->p4()).pt(),2)),
                (float)muon->isLooseMuon(), (float)muon->isTightMuon(PV), (float)muon->isHighPtMuon(PV),(float)utilities::heppymuonID(*muon, "POG_ID_Medium"), 
                (float)utilities::heppymuonID(*muon, "POG_ID_TightNoVtx"),  (float)utilities::relIso(*muon, 0.5)};
        mu->Fill(muarr);
    }
    else if(l1->isElectron()){
        const pat::Electron *electron = dynamic_cast<const pat::Electron*>((pairs->front()).daughter(0)->masterClone().get());
        float earr[] = {(float)electron->pt(), (float)electron->eta(),(float)electron->phi(),(float)electron->mass(), (float)electron->charge()};
        e->Fill(earr);
    }
    else{
        const pat::Tau *taon  = dynamic_cast<const pat::Tau*>((pairs->front()).daughter(1)->masterClone().get());
        float tauarr[] = {(float)taon->pt(), (float)taon->eta(), (float)taon->phi(), (float)taon->mass(), (float)taon->charge(), 
                (float)sqrt(pow((taon->p4()).pt() + (mety->p4()).pt(),2) - pow((taon->p4() + mety->p4()).pt(),2)),
                (float)taon->tauID("decayModeFinding"), taon->tauID("decayModeFindingNewDMs"), 
                (float)taon->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"), taon->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"), 
                (float)taon->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"), taon->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"), 
                (float)taon->tauID("chargedIsoPtSum"), taon->tauID("neutralIsoPtSum"), taon->tauID("puCorrPtSum"),taon->tauID("againstMuonLoose3"),
                (float)taon->tauID("againstMuonTight3"), taon->tauID("againstElectronVLooseMVA5"), taon->tauID("againstElectronLooseMVA5"), 
                (float)taon->tauID("againstElectronMediumMVA5"), taon->tauID("againstElectronTightMVA5"), taon->tauID("againstElectronVTightMVA5"), 
                (float)taon->tauID("byIsolationMVA3newDMwoLTraw"), (float)taon->tauID("byIsolationMVA3oldDMwoLTraw"), taon->tauID("byIsolationMVA3newDMwLTraw"), 
                (float)taon->tauID("byIsolationMVA3oldDMwLTraw")};
        tau_->Fill(tauarr);
    }

    const reco::Candidate * l2 = pairs->front().daughter(1);
    if(l2->isMuon()){
        const pat::Muon *muon  = dynamic_cast<const pat::Muon*>((pairs->front()).daughter(0)->masterClone().get());
        float muarr[] = {(float)muon->pt(), (float)muon->eta(),(float) muon->phi(),(float) muon->mass(), (float)muon->charge(),(float)(muon->innerTrack()->dxy( PV.position())),
               (float)(muon->innerTrack()->dz(PV.position())), (float)sqrt(pow((muon->p4()).pt() + (mety->p4()).pt(),2) - pow((muon->p4() + mety->p4()).pt(),2)),
                (float)muon->isLooseMuon(), (float)muon->isTightMuon(PV), (float)muon->isHighPtMuon(PV),(float)utilities::heppymuonID(*muon, "POG_ID_Medium"), 
                (float)utilities::heppymuonID(*muon, "POG_ID_TightNoVtx"),  (float)utilities::relIso(*muon, 0.5)};
        mu_->Fill(muarr);
    
    }
    else if(l2->isElectron()){
        const pat::Electron *electron = dynamic_cast<const pat::Electron*>((pairs->front()).daughter(1)->masterClone().get());
        float earr[] = {(float)electron->pt(), (float)electron->eta(),(float)electron->phi(),(float)electron->mass(), (float)electron->charge()};
        e_->Fill(earr);
    
    }
    else{
        const pat::Tau *taon  = dynamic_cast<const pat::Tau*>((pairs->front()).daughter(1)->masterClone().get());
        float tauarr[] = {(float)taon->pt(), (float)taon->eta(), (float)taon->phi(), (float)taon->mass(), (float)taon->charge(), 
                (float)sqrt(pow((taon->p4()).pt() + (mety->p4()).pt(),2) - pow((taon->p4() + mety->p4()).pt(),2)),
                (float)taon->tauID("decayModeFinding"), taon->tauID("decayModeFindingNewDMs"), 
                (float)taon->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"), taon->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"), 
                (float)taon->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"), taon->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"), 
                (float)taon->tauID("chargedIsoPtSum"), taon->tauID("neutralIsoPtSum"), taon->tauID("puCorrPtSum"),taon->tauID("againstMuonLoose3"),
                (float)taon->tauID("againstMuonTight3"), taon->tauID("againstElectronVLooseMVA5"), taon->tauID("againstElectronLooseMVA5"), 
                (float)taon->tauID("againstElectronMediumMVA5"), taon->tauID("againstElectronTightMVA5"), taon->tauID("againstElectronVTightMVA5"), 
                (float)taon->tauID("byIsolationMVA3newDMwoLTraw"), (float)taon->tauID("byIsolationMVA3oldDMwoLTraw"), taon->tauID("byIsolationMVA3newDMwLTraw"), 
                (float)taon->tauID("byIsolationMVA3oldDMwLTraw")};
        tau->Fill(tauarr);
    }


    float metarr[] = {(float)mety->px(), (float)mety->pt(), (float)mety->phi(),// mety->sumEt(),
        lP.userFloat("MEt_cov00"), lP.userFloat("MEt_cov01"), lP.userFloat("MEt_cov10"), lP.userFloat("MEt_cov11") };
    met->Fill(metarr);

    hltmatch.clear();  hltpaths.clear();
    std::string temp;
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
        obj.unpackPathNames(names);
        if(deltaR( obj.triggerObject().p4(), (pairs->front()).daughter(0)->p4()) < 0.5)
            for (unsigned h = 0; h < obj.filterLabels().size(); ++h) temp += (obj.filterLabels()[h] + "/");
            
    }
    hltmatch.push_back(temp);

    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
        if(triggerBits->accept(i)){
            std::string path = names.triggerName(i);
            if( path != "generation_step" &&  path !=  "digitisation_step" && path != "L1simulation_step"  &&  path != "digi2raw_step"  
                    &&  path != "HLT_ReducedIterativeTracking_v1"   &&  path != "HLT_ZeroBias_v1"  &&  path != "HLT_Physics_v1"  &&  path != "HLTriggerFinalPath")
            hltpaths.push_back(names.triggerName(i)); 
        }
    }
    trigger->Fill();


//    float* leadingjetarr = new float[leadingjet->GetNvar()]();
    std::vector<float> leadingjetarr(leadingjet->GetNvar(),0.);
    std::vector<float> trailingjetarr(trailingjet->GetNvar(),0.);
    std::vector<float> jetpairarr(jetpair->GetNvar(),0.);
    if(jets->size() > 0){
        leadingjetarr = {(float)jets->at(0).pt(), (float)jets->at(0).eta(), (float)jets->at(0).phi(),  (float)jets->at(0).userFloat("pileupJetId:fullDiscriminant"),
            (float)jets->at(0).bDiscriminator("jetBProbabilityBJetTags"), (float)jets->at(0).bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags"),
            (float)utilities::isbJet(jets->at(0)), jets->at(0).jecFactor("Uncorrected"), 
            (float)utilities::jetID(jets->at(0)), (float)utilities::pujetetaid(jets->at(0))};

        if(jets->size() > 1){
            trailingjetarr = {(float)jets->at(1).pt(), (float)jets->at(1).eta(), (float)jets->at(1).phi(),  (float)jets->at(1).userFloat("pileupJetId:fullDiscriminant"),
                jets->at(1).bDiscriminator("jetBProbabilityBJetTags"), jets->at(1).bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags"),
                (float)utilities::isbJet(jets->at(1)), jets->at(1).jecFactor("Uncorrected"), 
                (float)utilities::jetID(jets->at(1)), (float)utilities::pujetetaid(jets->at(1))};

            TLorentzVector jjp4 = TLorentzVector(jets->at(0).px(), jets->at(0).py(), jets->at(0).pz(), jets->at(0).energy()) +  
                TLorentzVector(jets->at(1).px(), jets->at(1).py(), jets->at(1).pz(), jets->at(1).energy());
            float njetingap=0;
            if(jets->size() > 2){
                for(std::vector<int>::size_type i = 2; i != jets->size(); i++){
                    if ( (*jets)[i].pt() < 30 || ! (*jets)[i].userFloat("pileupJetId:fullDiscriminant") || !utilities::pujetetaid( (*jets)[i]) )
                        continue;
                    if( (*jets)[i].eta() > std::min((*jets)[0].eta(),(*jets)[1].eta() ) &&  (*jets)[i].eta() < std::max((*jets)[0].eta(),(*jets)[1].eta() )   ){
                        njetingap+=1;         
                    }
                }
            }
            float hdijetphi = deltaPhi(jjp4.Phi(), ((pairs->front()).daughter(0)->p4() + (pairs->front()).daughter(1)->p4()).phi() );
            float visjeteta =  std::min( fabs(jets->at(0).eta() - ((pairs->front()).daughter(0)->p4() + (pairs->front()).daughter(1)->p4()).eta()), fabs(jets->at(1).eta() - ((pairs->front()).daughter(0)->p4() + (pairs->front()).daughter(1)->p4()).eta()));
            jetpairarr = {(float)jjp4.M(), (float)jjp4.Pt(), (float)jjp4.Phi(), (float)fabs( jets->at(0).eta() - jets->at(1).eta() ), (float)fabs( jets->at(0).phi() - jets->at(1).phi()), njetingap, hdijetphi, visjeteta};
            
        }
    }
    leadingjet->Fill(&leadingjetarr[0]);
    trailingjet->Fill(&trailingjetarr[0]);
    jetpair->Fill(&jetpairarr[0]);

}


// ------------ method called once each job just before starting event loop  ------------
void ntuple::beginJob(){

   edm::Service<TFileService> theFileService;
   event = theFileService->make<TTree>("event", "event");
   event->Branch("evt", &evt);
   event->Branch("nup",&nup);
   event->Branch("paircount",&paircount); 
   event->Branch("run",&run); 
   event->Branch("lumi",&lumi); 
   event->Branch("npv",&npv); 
   event->Branch("npu",&npu); 
    
   isZ = theFileService->make<TTree>("isZ", "isZ");
    isZ->Branch("isZtt", &isZtt);
    isZ->Branch("isZmt", &isZmt);
    isZ->Branch("isZet", &isZet);
    isZ->Branch("isZee", &isZee);
    isZ->Branch("isZmm", &isZmm);
    isZ->Branch("isZem", &isZem);
    isZ->Branch("isZEE", &isZEE);
    isZ->Branch("isZMM", &isZMM);
    isZ->Branch("isZLL", &isZLL);

    trigger = theFileService->make<TTree>("trigger", "trigger");
    trigger->Branch("hltmatch", &hltmatch);
    trigger->Branch("hltpaths", &hltpaths);


   mu = theFileService->make<TNtuple>("mu", "mu", "pt:eta:phi:mass:charge:d0:dz:mt:isLooseMuon:isTightMuon:isHighPtMuon:isMediumMuon:isTightnovtxMuon:iso");
   mu_ = theFileService->make<TNtuple>("mu_", "mu_", "pt:eta:phi:mass:charge:d0:dz:mt:isLooseMuon:isTightMuon:isHighPtMuon:isMediumMuon:isTightnovtxMuon:iso");

   e = theFileService->make<TNtuple>("e", "e", "pt:eta:phi:mass:charge");
   e_ = theFileService->make<TNtuple>("e_", "e_", "pt:eta:phi:mass:charge");

   tau = theFileService->make<TNtuple>("tau", "tau", "pt:eta:phi:mass:charge:mt:decayModeFinding:decayModeFindingNewDMs:byCombinedIsolationDeltaBetaCorrRaw3Hits:lbyCombinedIsolationDeltaBetaCorrRaw3Hits:mbyCombinedIsolationDeltaBetaCorrRaw3Hits:tbyCombinedIsolationDeltaBetaCorrRaw3Hits:chargedIsoPtSum:neutralIsoPtSum:puCorrPtSum:againstMuonLoose3:againstMuonTight3:againstElectronVLooseMVA5:againstElectronLooseMVA5:againstElectronMediumMVA5:againstElectronTightMVA5:againstElectronVTightMVA5:byIsolationMVA3newDMwoLTraw:byIsolationMVA3oldDMwoLTraw:byIsolationMVA3newDMwLTraw:byIsolationMVA3oldDMwLTraw");
   tau_ = theFileService->make<TNtuple>("tau_", "tau_", "pt:eta:phi:mass:charge:mt:decayModeFinding:decayModeFindingNewDMs:byCombinedIsolationDeltaBetaCorrRaw3Hits:lbyCombinedIsolationDeltaBetaCorrRaw3Hits:mbyCombinedIsolationDeltaBetaCorrRaw3Hits:tbyCombinedIsolationDeltaBetaCorrRaw3Hits:chargedIsoPtSum:neutralIsoPtSum:puCorrPtSum:againstMuonLoose3:againstMuonTight3:againstElectronVLooseMVA5:againstElectronLooseMVA5:againstElectronMediumMVA5:againstElectronTightMVA5:againstElectronVTightMVA5:byIsolationMVA3newDMwoLTraw:byIsolationMVA3oldDMwoLTraw:byIsolationMVA3newDMwLTraw:byIsolationMVA3oldDMwLTraw");

   pair =  theFileService->make<TNtuple>("pair", "pair", "svfit:diq:pth:ptvis:m_vis");

   met =  theFileService->make<TNtuple>("met", "met", "metpx:metpt:metphi:mvacov00:mvacov01:mvacov10:mvacov11"); //metsumEt

   leadingjet =  theFileService->make<TNtuple>("leadingjet", "leadingjet", "pt:eta:phi:id:bptag:csvtag:bjet:jecfactor:jetlooseID:pujetetaid");

   trailingjet =  theFileService->make<TNtuple>("trailingjet", "trailingjet", "pt:eta:phi:id:bptag:csvtag:bjet:jecfactor:jetlooseID:pujetetaid");

   jetpair = theFileService->make<TNtuple>("jetpair", "jetpair", "mass:pt:phi:deta:dphi:njetingap:hdijetphi:visjeteta");

   events = theFileService->make<TH1F>("hvar","hvar title",10,0,10);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
ntuple::endJob() 
{
}



// ------------ method called when starting to processes a run  ------------
/*
void 
ntuple::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
ntuple::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
ntuple::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
ntuple::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ntuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ntuple);


