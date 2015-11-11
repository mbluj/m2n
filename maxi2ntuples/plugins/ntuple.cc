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
#include "WarsawAnalysis/HTTDataFormats/interface/HTTEvent.h"

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
    const int sample;



//    TNtuple* evt;
    TTree * newevent; 
    TTree * isZ;
//    TTree * trigger;
//    TNtuple* leadingjet;
//    TNtuple* trailingjet;
//    TNtuple* jetpair;
    TH1F* events;

    bool isZtt, isZmt, isZet, isZee, isZmm, isZem, isZEE, isZMM, isZLL;

    Wevent *wevent;
    //Wtriggers *wtriggers;
    WtauCollection wtaucollection;
    WmuCollection wmucollection;
    WelectronCollection welectroncollection;
    WpairCollection wpaircollection;
    WmetCollection wmetcollection;
    WjetCollection wjetcollection;


    //std::vector<std::string> hltmatch;
    //std::vector<std::string> hltpaths;




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
    mc(iConfig.getParameter<bool>("mc")),
    sample(iConfig.getParameter<int>("sample"))
{
   //now do what ever initialization is needed
}


ntuple::~ntuple()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
    delete wevent;
    //delete wtriggers;
}


//
// member functions
//

// ------------ method called for each event  ------------
void
ntuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    wtaucollection.clear();
    wmucollection.clear();
    wpaircollection.clear();
    wmetcollection.clear();
    wjetcollection.clear();

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

        for(const PileupSummaryInfo &pusi : *genPileUpInfos){
            int bx = pusi.getBunchCrossing();
            int nPU = pusi.getPU_NumInteractions();
            if(bx == 0){
                wevent->npu(nPU);
            }
        }
        try{
            edm::Handle<LHEEventProduct> evnt;
            iEvent.getByToken(lheprodToken_, evnt);
            const lhef::HEPEUP hepeup_ = evnt->hepeup();
            wevent->nup(hepeup_.NUP);
        }
        catch(const std::exception& e){
            //std::cout << "No lheprod collection"; 
        }

    }

    wevent->run(iEvent.id().run());
    wevent->lumi(iEvent.luminosityBlock());
    wevent->event(iEvent.id().event());
    wevent->npv(vertices->size());
    wevent->paircount(pairs->size());
    wevent->sample(sample);


    events->Fill(1.,1.);
    if(mc){
        edm::Handle<GenEventInfoProduct> genEvt;
        iEvent.getByLabel("generator",genEvt);
       // event weight
        double weightevt=genEvt->weight(); 
        events->Fill(2., weightevt);
        wevent->genevtweight(weightevt);
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
    isZ->Fill();
    if (!pairs.isValid() || pairs->size()==0){
        newevent->Fill();
        return;
    }

    for (const pat::CompositeCandidate &lP : *pairs){
        
        const reco::Candidate * l1, *l2; 
        l1 = lP.daughter(0); l2 = lP.daughter(1);

        const reco::Candidate *mety = lP.daughter(2);

        Wpair wpair;
        wpair.svfit(lP.userFloat("SVfitMass"));
        wpair.diq(l1->charge() * l2->charge());
        wpair.pth((l1->p4() + l2->p4() + mety->p4()).pt());
        wpair.ptvis((l1->p4() + l2->p4()).pt());
        wpair.m_vis((l1->p4() + l2->p4()).mass());
        wpair.ChannelSelector(lP.userFloat("ChannelSelector"));
        wpair.PATPairSelector(lP.userFloat("PATPairSelector"));
        wpair.PairBaselineSelection(lP.userFloat("PairBaselineSelection"));
        wpair.PostSynchSelection(lP.userFloat("PostSynchSelection"));

        wpair.trigger(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1, lP.userFloat("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1"));
        wpair.trigger(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v1, lP.userFloat("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v1"));
        wpair.trigger(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2, lP.userFloat("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2"));
        wpair.trigger(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v2, lP.userFloat("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v2"));
        wpair.trigger(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3, lP.userFloat("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3"));
        wpair.trigger(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3, lP.userFloat("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3"));
        wpair.trigger(HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1, lP.userFloat("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1"));
        wpair.trigger(HLT_IsoMu24_eta2p1_v1, lP.userFloat("HLT_IsoMu24_eta2p1_v1"));
        wpair.trigger(HLT_IsoMu27_v1, lP.userFloat("HLT_IsoMu27_v1"));
        wpair.trigger(HLT_IsoMu18_v, lP.userFloat("HLT_IsoMu18_v"));
        wpair.trigger(HLT_IsoMu22_v1, lP.userFloat("HLT_IsoMu22_v1"));
        wpair.trigger(HLT_IsoMu17_eta2p1, lP.userFloat("HLT_IsoMu17_eta2p1"));
        wpair.trigger(HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v2, lP.userFloat("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v2"));
        wpair.trigger(HLT_IsoMu24_eta2p1_v2, lP.userFloat("HLT_IsoMu24_eta2p1_v2"));

        wpaircollection.push_back(wpair);

        if(l1->isMuon()){
            const pat::Muon *muon  = dynamic_cast<const pat::Muon*>(l1->masterClone().get());

            Wmu wmu;
            wmu.pt(muon->pt());
            wmu.eta(muon->eta());
            wmu.phi(muon->phi());
            wmu.mass(muon->mass());
            wmu.charge(muon->charge());
            wmu.mt(sqrt(pow((muon->p4()).pt() + (mety->p4()).pt(),2) - pow((muon->p4() + mety->p4()).pt(),2)));
            wmu.d0(muon->innerTrack()->dxy( PV.position()));
            wmu.dz(muon->innerTrack()->dz(PV.position()));
            wmu.isLooseMuon(muon->isLooseMuon());
            wmu.isTightMuon(muon->isTightMuon(PV));
            wmu.isHighPtMuon(muon->isHighPtMuon(PV));
            wmu.isMediumMuon(muon->isMediumMuon());
            wmu.isTightnovtxMuon(utilities::heppymuonID(*muon, "POG_ID_TightNoVtx"));
            wmu.iso(((muon->pfIsolationR03().sumChargedHadronPt
                    + std::max( muon->pfIsolationR03().sumNeutralHadronEt + muon->pfIsolationR03().sumPhotonEt - 0.5 * muon->pfIsolationR03().sumPUPt, 0.0)
                    ) / muon->pt())); 

            wmucollection.push_back(wmu);

        }
        else if(l1->isElectron()){
            const pat::Electron *electron = dynamic_cast<const pat::Electron*>(l1->masterClone().get());
            Welectron welectron;
            welectron.pt(electron->pt());
            welectron.eta(electron->eta());
            welectron.phi(electron->phi());
            welectron.mass(electron->mass());
            welectron.charge(electron->charge());
            welectroncollection.push_back(welectron);

        }
        else{
            //const pat::Tau *taon  = dynamic_cast<const pat::Tau*>(l1->masterClone().get());
        }

        if(l2->isMuon()){
            //const pat::Muon *muon  = dynamic_cast<const pat::Muon*>(l2->masterClone().get());
        
        }
        else if(l2->isElectron()){
            //const pat::Electron *electron = dynamic_cast<const pat::Electron*>(l2->masterClone().get());
        
        }
        else{

            const pat::Tau *taon  = dynamic_cast<const pat::Tau*>(l2->masterClone().get());

            Wtau wtau;
            wtau.pt(taon->pt());
            wtau.eta(taon->eta());
            wtau.phi(taon->phi());
            wtau.mass(taon->mass());
            wtau.charge(taon->charge());
            wtau.mt(sqrt(pow((taon->p4()).pt() + (mety->p4()).pt(),2) - pow((taon->p4() + mety->p4()).pt(),2)));
            wtau.tauID(decayModeFinding, taon->tauID("decayModeFinding"));
            wtau.tauID(decayModeFindingNewDMs, taon->tauID("decayModeFindingNewDMs"));
            wtau.tauID(byCombinedIsolationDeltaBetaCorrRaw3Hits, taon->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"));
            wtau.tauID(byLooseCombinedIsolationDeltaBetaCorr3Hits, taon->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"));
            wtau.tauID(byMediumCombinedIsolationDeltaBetaCorr3Hits, taon->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"));
            wtau.tauID(byTightCombinedIsolationDeltaBetaCorr3Hits, taon->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"));
            wtau.tauID(chargedIsoPtSum, taon->tauID("chargedIsoPtSum"));
            wtau.tauID(neutralIsoPtSum, taon->tauID("neutralIsoPtSum"));
            wtau.tauID(puCorrPtSum, taon->tauID("puCorrPtSum"));
            wtau.tauID(againstMuonLoose3, taon->tauID("againstMuonLoose3"));
            wtau.tauID(againstMuonTight3, taon->tauID("againstMuonTight3"));
            wtau.tauID(againstElectronVLooseMVA5, taon->tauID("againstElectronVLooseMVA5"));
            wtau.tauID(againstElectronLooseMVA5, taon->tauID("againstElectronLooseMVA5"));
            wtau.tauID(againstElectronMediumMVA5, taon->tauID("againstElectronMediumMVA5"));
            wtau.tauID(againstElectronTightMVA5, taon->tauID("againstElectronTightMVA5"));
            wtau.tauID(againstElectronVTightMVA5, taon->tauID("againstElectronVTightMVA5"));
            wtau.tauID(byIsolationMVA3newDMwoLTraw, taon->tauID("byIsolationMVA3newDMwoLTraw"));
            wtau.tauID(byIsolationMVA3oldDMwoLTraw, taon->tauID("byIsolationMVA3oldDMwoLTraw"));
            wtau.tauID(byIsolationMVA3newDMwLTraw, taon->tauID("byIsolationMVA3newDMwLTraw"));
            wtau.tauID(byIsolationMVA3oldDMwLTraw, taon->tauID("byIsolationMVA3oldDMwLTraw"));

            wtaucollection.push_back(wtau);

        }

        Wmet wmet;
        wmet.metpx(mety->px());
        wmet.metpt(mety->pt());
        wmet.metphi(mety->phi());
        wmet.mvacov00(lP.userFloat("MEt_cov00"));
        wmet.mvacov01(lP.userFloat("MEt_cov01"));
        wmet.mvacov10(lP.userFloat("MEt_cov10"));
        wmet.mvacov11(lP.userFloat("MEt_cov11"));
        wmetcollection.push_back(wmet);





    }
    /*
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
    */

//    float* leadingjetarr = new float[leadingjet->GetNvar()]();
    /*
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
    */

    int cou = 0;
    for(const pat::Jet &j : *jets){
        if (cou >1)
            break;
        Wjet wjet;
        wjet.pt(j.pt());
        wjet.eta(j.eta());
        wjet.phi(j.phi());
        wjet.id(j.userFloat("pileupJetId:fullDiscriminant"));
        wjet.bptag(j.bDiscriminator("jetBProbabilityBJetTags"));
        wjet.csvtag(j.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags"));
        wjet.bjet(utilities::isbJet(j));
        wjet.jecfactor(j.jecFactor("Uncorrected"));
        wjet.jetlooseID(utilities::jetID(j));
        wjet.pujetetaid(utilities::pujetetaid(j));
        wjetcollection.push_back(wjet);
        cou+=1;
    }
    if (cou ==2){
    
    
    
    }

    newevent->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void ntuple::beginJob(){

   edm::Service<TFileService> theFileService;
    

   newevent = theFileService->make<TTree>("newevent", "newevent");
   wevent = 0;
   //wtriggers = 0;
   newevent->Branch("wevent", &wevent);
   newevent->Branch("wtau", &wtaucollection);
   newevent->Branch("wmu", &wmucollection);
   newevent->Branch("welectron", &welectroncollection);
   newevent->Branch("wpair", &wpaircollection);
   newevent->Branch("wmet", &wmetcollection);
   newevent->Branch("wjet", &wjetcollection);
   //newevent->Branch("wtriggers", &wtriggers);

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

    //trigger = theFileService->make<TTree>("trigger", "trigger");
    //trigger->Branch("hltmatch", &hltmatch);
    //trigger->Branch("hltpaths", &hltpaths);


   //leadingjet =  theFileService->make<TNtuple>("leadingjet", "leadingjet", "pt:eta:phi:id:bptag:csvtag:bjet:jecfactor:jetlooseID:pujetetaid");

   //trailingjet =  theFileService->make<TNtuple>("trailingjet", "trailingjet", "pt:eta:phi:id:bptag:csvtag:bjet:jecfactor:jetlooseID:pujetetaid");

   //jetpair = theFileService->make<TNtuple>("jetpair", "jetpair", "mass:pt:phi:deta:dphi:njetingap:hdijetphi:visjeteta");

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


