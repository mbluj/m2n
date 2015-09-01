//
// Package:    m2n/maxi2ntuples
// Class:      synchronization
// 
/**\class synchronization synchronization.cc m2n/maxi2ntuples/plugins/synchronization.cc

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

//#include "clasadict.h"

//
// class declaration
//

//gSystem->Load("libPhysics.so");
class synchronization : public edm::EDAnalyzer {
   public:
      explicit synchronization(const edm::ParameterSet&);
      ~synchronization();

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

      void fillmuon(const pat::Muon*, const reco::Vertex&, const reco::Candidate*);
      void filltau(const pat::Tau*, const reco::Candidate*);

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
    TNtuple* synchtree;

     
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
synchronization::synchronization(const edm::ParameterSet& iConfig):
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


synchronization::~synchronization()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
synchronization::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

    edm::Handle<reco::GenParticleCollection> pruned; 
    iEvent.getByToken(prunedGenToken_, pruned);
//    Handle<edm::View<pat::PackedGenParticle> > packed;
//    iEvent.getByToken(packedGenToken_,packed);
    edm::Handle<pat::PackedGenParticleCollection> packed; 
    iEvent.getByToken(packedGenToken_, packed);

    float npu =0, nup=0;
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
            //std::cout << "No lheprod (LHEEventProduct) collection"; 
        }

    }

    
    float isZtt = 0, isZmt = 0, isZet = 0, isZee = 0, isZmm = 0, isZem = 0, isZEE = 0, isZMM = 0, isZLL = 0;
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


    const reco::Candidate *l1 = (pairs->front()).daughter(0);
    const reco::Candidate *l2 = (pairs->front()).daughter(1);
    const reco::Candidate *met = (pairs->front()).daughter(2);

    const pat::CompositeCandidate& lP = pairs->front();
    
    std::vector<float> arr(synchtree->GetNvar(),0.);
    
      arr[0] = (float)iEvent.id().run();   //run
      arr[1] = (float)iEvent.luminosityBlock();   //lumi
      arr[2] = (float)iEvent.id().event();   //evt
      arr[3] = (float)isZtt;   //isZtt
      arr[4] = (float)isZmt;   //isZmt
      arr[5] = (float)isZet;   //isZet
      arr[6] = (float)isZee;   //isZee
      arr[7] = (float)isZmm;   //isZmm
      arr[8] = (float)isZem;   //isZem
      arr[9] = (float)isZEE;   //isZEE
      arr[10] = (float)isZMM;   //isZMM
      arr[11] = (float)isZLL;   //isZLL
      arr[12] = (float)0;   //isFake
      arr[13] = (float)nup;   //NUP
      arr[14] = (float)0;   //weight
      arr[15] = (float)0;   //puweight
      arr[16] = (float)vertices->size();   //npv
      arr[17] = (float)npu;   //npu
      arr[18] = (float)0;   //rho
      arr[19] = (float)l1->pt();   //pt_1
      arr[20] = (float)l1->phi();   //phi_1
      arr[21] = (float)l1->eta();   //eta_1
      arr[22] = (float)l1->mass();   //m_1
      arr[23] = (float)l1->charge();   //q_1

      if(l1->isMuon()){
          const pat::Muon *mu = dynamic_cast<const pat::Muon*>(l1->masterClone().get());
          arr[24] = (float) mu->innerTrack()->dxy( PV.position());   //d0_1
          arr[25] = (float) mu->innerTrack()->dz(PV.position());   //dZ_1
          arr[27] = (float) utilities::relIso(*mu, 0.5);   //iso_1
          arr[28] = (float) mu->isLooseMuon();   //id_m_loose_1
          arr[29] = (float) utilities::heppymuonID(*mu, "POG_ID_Medium");   //id_m_medium_1
          arr[30] = (float) mu->isTightMuon(PV);   //id_m_tight_1
          arr[31] = (float) utilities::heppymuonID(*mu, "POG_ID_TightNoVtx");   //id_m_tightnovtx_1
          arr[32] = (float) mu->isHighPtMuon(PV);   //id_m_highpt_1
      }
      arr[26] = (float) sqrt(pow((l1->p4()).pt()+(met->p4()).pt(),2)-pow((l1->p4()+met->p4()).pt(),2));//mt_1

      arr[33] = (float)0;   //id_e_mva_nt_loose_1
      arr[34] = (float)0;   //id_e_cut_veto_1
      arr[35] = (float)0;   //id_e_cut_loose_1
      arr[36] = (float)0;   //id_e_cut_medium_1
      arr[37] = (float)0;   //id_e_cut_tight_1
      arr[38] = (float)0;   //trigweight_1
      arr[39] = (float)0;   //againstElectronLooseMVA5_1
      arr[40] = (float)0;   //againstElectronMediumMVA5_1
      arr[41] = (float)0;   //againstElectronTightMVA5_1
      arr[42] = (float)0;   //againstElectronVLooseMVA5_1
      arr[43] = (float)0;   //againstElectronVTightMVA5_1
      arr[44] = (float)0;   //againstMuonLoose3_1
      arr[45] = (float)0;   //againstMuonTight3_1
      arr[46] = (float)0;   //byCombinedIsolationDeltaBetaCorrRaw3Hits_1
      arr[47] = (float)0;   //byIsolationMVA3newDMwoLTraw_1
      arr[48] = (float)0;   //byIsolationMVA3oldDMwoLTraw_1
      arr[49] = (float)0;   //byIsolationMVA3newDMwLTraw_1
      arr[50] = (float)0;   //byIsolationMVA3oldDMwLTraw_1
      arr[51] = (float)0;   //chargedIsoPtSum_1
      arr[52] = (float)0;   //decayModeFinding_1
      arr[53] = (float)0;   //decayModeFindingNewDMs_1
      arr[54] = (float)0;   //neutralIsoPtSum_1
      arr[55] = (float)0;   //puCorrPtSum_1
      arr[56] = (float)l2->pt();   //pt_2
      arr[57] = (float)l2->phi();   //phi_2
      arr[58] = (float)l2->eta();   //eta_2
      arr[59] = (float)l2->mass();   //m_2
      arr[60] = (float)l2->charge();   //q_2
      arr[61] = (float)0;   //d0_2
      arr[62] = (float)0;   //dZ_2
      arr[63] = (float)sqrt(pow((l2->p4()).pt()+(met->p4()).pt(),2)-pow((l2->p4()+met->p4()).pt(),2)); //mt_2
      arr[64] = (float)0;   //iso_2
      arr[65] = (float)0;   //id_m_loose_2
      arr[66] = (float)0;   //id_m_medium_2
      arr[67] = (float)0;   //id_m_tight_2
      arr[68] = (float)0;   //id_m_tightnovtx_2
      arr[69] = (float)0;   //id_m_highpt_2
      arr[70] = (float)0;   //id_e_mva_nt_loose_2
      arr[71] = (float)0;   //id_e_cut_veto_2
      arr[72] = (float)0;   //id_e_cut_loose_2
      arr[73] = (float)0;   //id_e_cut_medium_2
      arr[74] = (float)0;   //id_e_cut_tight_2
      arr[75] = (float)0;   //trigweight_2
      if(!(l2->isMuon() || l2->isElectron())){
          const pat::Tau *tau = dynamic_cast<const pat::Tau*>(l2->masterClone().get());
          arr[76] = (float)tau->tauID("againstElectronVLooseMVA5");   //againstElectronLooseMVA5_2
          arr[77] = (float)tau->tauID("againstElectronLooseMVA5");   //againstElectronMediumMVA5_2
          arr[78] = (float)tau->tauID("againstElectronMediumMVA5");   //againstElectronTightMVA5_2
          arr[79] = (float)tau->tauID("againstElectronTightMVA5");   //againstElectronVLooseMVA5_2
          arr[80] = (float)tau->tauID("againstElectronVTightMVA5");   //againstElectronVTightMVA5_2
          arr[81] = (float)tau->tauID("againstMuonLoose3");   //againstMuonLoose3_2
          arr[82] = (float)tau->tauID("againstMuonTight3");   //againstMuonTight3_2
          arr[83] = (float)tau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"); //byCombinedIsolationDeltaBetaCorrRaw3Hits_2
          arr[84] = (float)tau->tauID("byIsolationMVA3newDMwoLTraw");   //byIsolationMVA3newDMwoLTraw_2
          arr[85] = (float)tau->tauID("byIsolationMVA3oldDMwoLTraw");   //byIsolationMVA3oldDMwoLTraw_2
          arr[86] = (float)tau->tauID("byIsolationMVA3newDMwLTraw");   //byIsolationMVA3newDMwLTraw_2
          arr[87] = (float)tau->tauID("byIsolationMVA3oldDMwLTraw");   //byIsolationMVA3oldDMwLTraw_2
          arr[88] = (float)0;   //chargedIsoPtSum_2
          arr[89] = (float)tau->tauID("decayModeFinding");   //decayModeFinding_2
          arr[90] = (float)tau->tauID("decayModeFindingNewDMs");   //decayModeFindingNewDMs_2
          arr[91] = (float)0;   //neutralIsoPtSum_2
          arr[92] = (float)0;   //puCorrPtSum_2
      }
      arr[93] = (float)0;   //pth
      arr[94] = (float)((pairs->front()).daughter(0)->p4()+(pairs->front()).daughter(1)->p4()).mass();//m_vis
      arr[95] = (float)lP.userFloat("SVfitMass");   //m_sv
      arr[96] = (float)((pairs->front()).daughter(0)->p4()+(pairs->front()).daughter(1)->p4()).pt();//pt_sv
      arr[97] = (float)0;   //eta_sv
      arr[98] = (float)0;   //phi_sv
      arr[99] = (float)0;   //met_sv
      arr[100] = (float)0;   //met
      arr[101] = (float)0;   //metphi
      arr[102] = (float)0;   //mvamet
      arr[103] = (float)0;   //mvametphi
      arr[104] = (float)0;   //pzetavis
      arr[105] = (float)0;   //pzetamiss
      arr[106] = (float)0;   //mvacov00
      arr[107] = (float)0;   //mvacov01
      arr[108] = (float)0;   //mvacov10
      arr[109] = (float)0;   //mvacov11
      arr[110] = (float)0;   //mjj
      arr[111] = (float)0;   //jdeta
      arr[112] = (float)0;   //njetingap
      arr[113] = (float)0;   //jdphi
      arr[114] = (float)0;   //dijetpt
      arr[115] = (float)0;   //dijetphi
      arr[116] = (float)0;   //hdijetphi
      arr[117] = (float)0;   //visjeteta
      arr[118] = (float)0;   //ptvis
      arr[119] = (float)0;   //nbtag
      arr[120] = (float)0;   //njets
      arr[121] = (float)0;   //njetspt20
      arr[122] = (float)0;   //jpt_1
      arr[123] = (float)0;   //jeta_1
      arr[124] = (float)0;   //jphi_1
      arr[125] = (float)0;   //jrawf_1
      arr[126] = (float)0;   //jmva_1
      arr[127] = (float)0;   //jpfid_1
      arr[128] = (float)0;   //jpuid_1
      arr[129] = (float)0;   //jcsv_1
      arr[130] = (float)0;   //jpt_2
      arr[131] = (float)0;   //jeta_2
      arr[132] = (float)0;   //jphi_2
      arr[133] = (float)0;   //jrawf_2
      arr[134] = (float)0;   //jmva_2
      arr[135] = (float)0;   //jpfid_2
      arr[136] = (float)0;   //jpuid_2
      arr[137] = (float)0;   //jcsv_2
      arr[138] = (float)0;   //bpt_1
      arr[139] = (float)0;   //beta_1
      arr[140] = (float)0;   //bphi_1
      arr[141] = (float)0;   //brawf_1
      arr[142] = (float)0;   //bmva_1
      arr[143] = (float)0;   //bpfid_1
      arr[144] = (float)0;   //bpuid_1
      arr[145] = (float)0;   //bcsv_1
      arr[146] = (float)0;   //bpt_2
      arr[147] = (float)0;   //beta_2
      arr[148] = (float)0;   //bphi_2
      arr[149] = (float)0;   //brawf_2
      arr[150] = (float)0;   //bmva_2
      arr[151] = (float)0;   //bpfid_2
      arr[152] = (float)0;   //bpuid_2
      arr[153] = (float)0;   //bcsv_2



    synchtree->Fill(&arr[0]);

}


// ------------ method called once each job just before starting event loop  ------------
void synchronization::beginJob(){

    edm::Service<TFileService> theFileService;

//   synchtree = theFileService->make<TNtuple>("mu", "mu", "pt:eta:phi:mass:charge:d0:dz:mt:isLooseMuon:isTightMuon:isHighPtMuon:isMediumMuon:isTightnovtxMuon:iso");
    synchtree=theFileService->make<TNtuple>("synchtree","synchtree","run:lumi:evt:isZtt:isZmt:isZet:isZee:isZmm:isZem:isZEE:isZMM:isZLL:isFake:NUP:weight:puweight:npv:npu:rho:pt_1:phi_1:eta_1:m_1:q_1:d0_1:dZ_1:mt_1:iso_1:id_m_loose_1:id_m_medium_1:id_m_tight_1:id_m_tightnovtx_1:id_m_highpt_1:id_e_mva_nt_loose_1:id_e_cut_veto_1:id_e_cut_loose_1:id_e_cut_medium_1:id_e_cut_tight_1:trigweight_1:againstElectronLooseMVA5_1:againstElectronMediumMVA5_1:againstElectronTightMVA5_1:againstElectronVLooseMVA5_1:againstElectronVTightMVA5_1:againstMuonLoose3_1:againstMuonTight3_1:byCombinedIsolationDeltaBetaCorrRaw3Hits_1:byIsolationMVA3newDMwoLTraw_1:byIsolationMVA3oldDMwoLTraw_1:byIsolationMVA3newDMwLTraw_1:byIsolationMVA3oldDMwLTraw_1:chargedIsoPtSum_1:decayModeFinding_1:decayModeFindingNewDMs_1:neutralIsoPtSum_1:puCorrPtSum_1:pt_2:phi_2:eta_2:m_2:q_2:d0_2:dZ_2:mt_2:iso_2:id_m_loose_2:id_m_medium_2:id_m_tight_2:id_m_tightnovtx_2:id_m_highpt_2:id_e_mva_nt_loose_2:id_e_cut_veto_2:id_e_cut_loose_2:id_e_cut_medium_2:id_e_cut_tight_2:trigweight_2:againstElectronLooseMVA5_2:againstElectronMediumMVA5_2:againstElectronTightMVA5_2:againstElectronVLooseMVA5_2:againstElectronVTightMVA5_2:againstMuonLoose3_2:againstMuonTight3_2:byCombinedIsolationDeltaBetaCorrRaw3Hits_2:byIsolationMVA3newDMwoLTraw_2:byIsolationMVA3oldDMwoLTraw_2:byIsolationMVA3newDMwLTraw_2:byIsolationMVA3oldDMwLTraw_2:chargedIsoPtSum_2:decayModeFinding_2:decayModeFindingNewDMs_2:neutralIsoPtSum_2:puCorrPtSum_2:pth:m_vis:m_sv:pt_sv:eta_sv:phi_sv:met_sv:met:metphi:mvamet:mvametphi:pzetavis:pzetamiss:mvacov00:mvacov01:mvacov10:mvacov11:mjj:jdeta:njetingap:jdphi:dijetpt:dijetphi:hdijetphi:visjeteta:ptvis:nbtag:njets:njetspt20:jpt_1:jeta_1:jphi_1:jrawf_1:jmva_1:jpfid_1:jpuid_1:jcsv_1:jpt_2:jeta_2:jphi_2:jrawf_2:jmva_2:jpfid_2:jpuid_2:jcsv_2:bpt_1:beta_1:bphi_1:brawf_1:bmva_1:bpfid_1:bpuid_1:bcsv_1:bpt_2:beta_2:bphi_2:brawf_2:bmva_2:bpfid_2:bpuid_2:bcsv_2");



}

// ------------ method called once each job just after ending the event loop  ------------
void 
synchronization::endJob() 
{
}

void synchronization::fillmuon(const pat::Muon* mu,  const reco::Vertex& PV, const reco::Candidate* met){

}


void synchronization::filltau(const pat::Tau* tau, const reco::Candidate*met ){

}


// ------------ method called when starting to processes a run  ------------
/*
void 
synchronization::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
synchronization::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
synchronization::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
synchronization::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
synchronization::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(synchronization);


