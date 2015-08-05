//
// Package:    m2n/maxi2ntuples
// Class:      maxi2ntuples
// 
/**\class maxi2ntuples maxi2ntuples.cc m2n/maxi2ntuples/plugins/maxi2ntuples.cc

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
class maxi2ntuples : public edm::EDAnalyzer {
   public:
      explicit maxi2ntuples(const edm::ParameterSet&);
      ~maxi2ntuples();

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
//    edm::EDGetTokenT<reco::GenParticleCollection> prunedGenToken_;
//    edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGenToken_;
//    edm::EDGetTokenT<LHEEventProduct> lheprodToken_;
//    edm::EDGetTokenT<PileupSummaryInfoCollection> PileupSummaryInfoToken_;



    TTree * mt;                                                                                                //<-------------------------------------------------
    ///////////////////////////// P A I R     S P E C I F I C: ////////////////////////////////

    std::vector<std::pair<float, float>> muon;
    std::vector<float> 
        mupt, muphi, mueta, mum, muq, mud0, mudz, mumt, isLooseMuon, isTightMuon, isHighPtMuon, isMediumMuon, isTightnovtxMuon,
        taupt, tauphi, taueta, taum, tauq, taumt,
        svfit, metpx, metpt, metphi, metsumEt, muiso, decayModeFinding;
    std::vector<int> diq;
    std::vector<float>  // decayModeFindingOldDMs, 
        decayModeFindingNewDMs, byCombinedIsolationDeltaBetaCorrRaw3Hits, lbyCombinedIsolationDeltaBetaCorrRaw3Hits, mbyCombinedIsolationDeltaBetaCorrRaw3Hits, tbyCombinedIsolationDeltaBetaCorrRaw3Hits, chargedIsoPtSum, neutralIsoPtSum,
        puCorrPtSum, againstMuonLoose3, againstMuonTight3, againstElectronVLooseMVA5, againstElectronLooseMVA5, againstElectronMediumMVA5, againstElectronTightMVA5, againstElectronVTightMVA5,
        byIsolationMVA3newDMwoLTraw, byIsolationMVA3oldDMwoLTraw, byIsolationMVA3newDMwLTraw, byIsolationMVA3oldDMwLTraw,
        pth, ptvis, m_vis,
        mvacov00, mvacov01, mvacov10, mvacov11,
        hdijetphi, visjeteta; 
    std::vector<std::string> hltmatch;

    ///////////////////////////// E V E N T     S P E C I F I C: ////////////////////////////////
    //trigger paths
    std::vector<std::string> hltpaths;
    //jets
    std::vector<float> jetpt, jeteta, jetphi, pujetid, jetbptag, jetcsvtag, bjet, jecfactor, jetlooseID, pujetetaid;
    float  run, lumi, evt, npv, npu, mjj, ptjj, phijj, deta, dphi, njetingap;
    unsigned short paircount;//, mu_dz ; 
    //genParticles
    bool isZtt, isZmt, isZet, isZee, isZmm, isZem, isZEE, isZMM, isZLL;
    int  isFake, nup;


     
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
maxi2ntuples::maxi2ntuples(const edm::ParameterSet& iConfig):
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
    triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales")))
//    prunedGenToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
//    packedGenToken_(consumes<pat::PackedGenParticleCollection>(iConfig.getParameter<edm::InputTag>("packedGenParticles"))),
//    lheprodToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheprod"))),
//    PileupSummaryInfoToken_(consumes<PileupSummaryInfoCollection>(iConfig.getParameter<edm::InputTag>("pileupinfo")))
{
   //now do what ever initialization is needed
    edm::Service<TFileService> theFileService;
    mt = theFileService->make<TTree>("ntuple", "mt");
    mt->Branch("muon", &muon);
    mt->Branch("mupt", &mupt);
    mt->Branch("muphi", &muphi);
    mt->Branch("mueta", &mueta);
    mt->Branch("mum", &mum);
    mt->Branch("muq", &muq);
    mt->Branch("mud0", &mud0);
    mt->Branch("mudz", &mudz);
    mt->Branch("mumt", &mumt);
    mt->Branch("taupt", &taupt);
    mt->Branch("tauphi", &tauphi);
    mt->Branch("taueta", &taueta);
    mt->Branch("taum", &taum);
    mt->Branch("tauq", &tauq);
    mt->Branch("taumt", &taumt);
    mt->Branch("svfit", &svfit);
    mt->Branch("metpx", &metpx);
    mt->Branch("metpt", &metpt);
    mt->Branch("metphi", &metphi);                                                                        //<-----------------------------------------------
    mt->Branch("metsumEt", &metsumEt);
    mt->Branch("hltmatch", &hltmatch);
    mt->Branch("hltpaths", &hltpaths);
    mt->Branch("decayModeFinding", &decayModeFinding);
//    mt->Branch("decayModeFindingOldDMs", &decayModeFindingOldDMs);
    mt->Branch("decayModeFindingNewDMs", &decayModeFindingNewDMs);
    mt->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits", &byCombinedIsolationDeltaBetaCorrRaw3Hits );
    mt->Branch("lbyCombinedIsolationDeltaBetaCorrRaw3Hits", &lbyCombinedIsolationDeltaBetaCorrRaw3Hits );
    mt->Branch("mbyCombinedIsolationDeltaBetaCorrRaw3Hits", &mbyCombinedIsolationDeltaBetaCorrRaw3Hits );
    mt->Branch("tbyCombinedIsolationDeltaBetaCorrRaw3Hits", &tbyCombinedIsolationDeltaBetaCorrRaw3Hits );
    mt->Branch("chargedIsoPtSum", &chargedIsoPtSum );
    mt->Branch("neutralIsoPtSum", &neutralIsoPtSum );
    mt->Branch("puCorrPtSum", &puCorrPtSum );
    mt->Branch("againstMuonLoose3", &againstMuonLoose3 );
    mt->Branch("againstMuonTight3", &againstMuonTight3 );
    mt->Branch("againstElectronVLooseMVA5", &againstElectronVLooseMVA5 );
    mt->Branch("againstElectronLooseMVA5", &againstElectronLooseMVA5 );
    mt->Branch("againstElectronMediumMVA5", &againstElectronMediumMVA5 );
    mt->Branch("againstElectronTightMVA5", &againstElectronTightMVA5);
    mt->Branch("againstElectronVTightMVA5", &againstElectronVTightMVA5);
    mt->Branch("byIsolationMVA3newDMwoLTraw",&byIsolationMVA3newDMwoLTraw); 
    mt->Branch("byIsolationMVA3oldDMwoLTraw",&byIsolationMVA3oldDMwoLTraw); 
    mt->Branch("byIsolationMVA3newDMwLTraw",&byIsolationMVA3newDMwLTraw); 
    mt->Branch("byIsolationMVA3oldDMwLTraw",&byIsolationMVA3oldDMwLTraw); 
    mt->Branch("diq", &diq);

    mt->Branch("pth", &pth);
    mt->Branch("ptvis", &ptvis);
    mt->Branch("m_vis", &m_vis); 

    mt->Branch("mvacov00",&mvacov00); 
    mt->Branch("mvacov01",&mvacov01); 
    mt->Branch("mvacov10",&mvacov10); 
    mt->Branch("mvacov11",&mvacov11); 

    mt->Branch("hdijetphi",&hdijetphi); 
    mt->Branch("visjeteta",&visjeteta); 

    mt->Branch("muiso", &muiso);
    mt->Branch("jetpt", &jetpt);
    mt->Branch("jeteta", &jeteta);
    mt->Branch("jetphi", &jetphi);
    mt->Branch("pujetid", &pujetid);
    mt->Branch("jetbptag", &jetbptag);
    mt->Branch("jetcsvtag", &jetcsvtag);
    mt->Branch("jecfactor", &jecfactor);
    mt->Branch("jetlooseID", &jetlooseID);
    mt->Branch("pujetetaid", &pujetetaid);
    mt->Branch("bjet", &bjet);
    mt->Branch("mjj", &mjj);
    mt->Branch("ptjj", &ptjj);
    mt->Branch("phijj", &phijj);
    mt->Branch("deta", &deta);
    mt->Branch("dphi", &dphi);
    mt->Branch("njetingap", &njetingap);

    mt->Branch("isLooseMuon", &isLooseMuon);
    mt->Branch("isTightMuon", &isTightMuon);
    mt->Branch("isHighPtMuon", &isHighPtMuon);
    mt->Branch("isMediumMuon", &isMediumMuon);
    mt->Branch("isTightnovtxMuon", &isTightnovtxMuon);
    mt->Branch("run", &run);
    mt->Branch("lumi", &lumi);
    mt->Branch("evt", &evt);
    mt->Branch("npv", &npv);
    mt->Branch("npu", &npu);
    mt->Branch("paircount", &paircount);
  
    mt->Branch("isZtt", &isZtt);
    mt->Branch("isZmt", &isZmt);
    mt->Branch("isZet", &isZet);
    mt->Branch("isZee", &isZee);
    mt->Branch("isZmm", &isZmm);
    mt->Branch("isZem", &isZem);
    mt->Branch("isZEE", &isZEE);
    mt->Branch("isZMM", &isZMM);
    mt->Branch("isZLL", &isZLL);
    mt->Branch("isFake", &isFake);
    mt->Branch("nup", &nup);
}


maxi2ntuples::~maxi2ntuples()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
maxi2ntuples::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
 
    mupt.clear(); muphi.clear(); mueta.clear(); mum.clear(); muq.clear(); mud0.clear(); mudz.clear(); mumt.clear();
    isLooseMuon.clear(); isTightMuon.clear(); isHighPtMuon.clear(); isMediumMuon.clear(); isTightnovtxMuon.clear();
    taupt.clear(); tauphi.clear(); taueta.clear(); taum.clear(); tauq.clear(); taumt.clear();
    svfit.clear(); 
    metpx.clear(); metpt.clear(); metphi.clear(); metsumEt.clear();                                    //<----------------------------------------------------
    hltmatch.clear();  hltpaths.clear();
  //  decayModeFindingOldDMs.clear(); 
    decayModeFindingNewDMs.clear();
    diq.clear(); pth.clear(); ptvis.clear(); m_vis.clear();     

    mvacov00.clear();mvacov01.clear();mvacov10.clear();mvacov11.clear();
    hdijetphi.clear(); visjeteta.clear();

    muiso.clear(); decayModeFinding.clear();
    byCombinedIsolationDeltaBetaCorrRaw3Hits.clear(); lbyCombinedIsolationDeltaBetaCorrRaw3Hits.clear();  mbyCombinedIsolationDeltaBetaCorrRaw3Hits.clear(); tbyCombinedIsolationDeltaBetaCorrRaw3Hits.clear(); chargedIsoPtSum.clear();
    neutralIsoPtSum.clear(); puCorrPtSum.clear(); againstMuonLoose3.clear(); againstMuonTight3.clear(); 
    againstElectronVLooseMVA5.clear(); againstElectronLooseMVA5.clear(); againstElectronMediumMVA5.clear(); againstElectronTightMVA5.clear(); againstElectronVTightMVA5.clear();
    byIsolationMVA3newDMwoLTraw.clear();  byIsolationMVA3oldDMwoLTraw.clear();  byIsolationMVA3newDMwLTraw.clear();  byIsolationMVA3oldDMwLTraw.clear() ;
 
    jetpt.clear(); jeteta.clear(); jetphi.clear(); pujetid.clear(); jetbptag.clear(); jetcsvtag.clear(); bjet.clear(); jecfactor.clear();
    jetlooseID.clear(); pujetetaid.clear();
    mjj = -1; ptjj = -1; phijj = -10; deta = -10; dphi = -1; isFake = -1; npu = -1; njetingap = 0;
    lumi  = -1; evt = -1;

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

//    edm::Handle<reco::GenParticleCollection> pruned; 
//    iEvent.getByToken(prunedGenToken_, pruned);
//    edm::Handle<pat::PackedGenParticleCollection> packed; 
//    iEvent.getByToken(packedGenToken_, packed);

//    edm::Handle<LHEEventProduct> evnt;
//    iEvent.getByToken(lheprodToken_, evnt);

    /*
    edm::Handle<PileupSummaryInfoCollection> genPileUpInfos;
    iEvent.getByToken(PileupSummaryInfoToken_, genPileUpInfos);

    for(const PileupSummaryInfo &pusi : *genPileUpInfos){
        int bx = pusi.getBunchCrossing();
        int nPU = pusi.getPU_NumInteractions();
        if(bx == 0)
            npu = nPU;
    }
*/
//    const lhef::HEPEUP hepeup_ = evnt->hepeup();
//    nup =  hepeup_.NUP;
    
    
    paircount = pairs->size();
    run = iEvent.id().run();
    lumi = iEvent.luminosityBlock();
    evt = iEvent.id().event(); 
    npv =  vertices->size();

//    metpx.push_back(smet.px());
//    metpt.push_back(smet.pt());
//    metphi.push_back(smet.phi());
//    metsumEt.push_back(smet.sumEt());

    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
        if(triggerBits->accept(i)){
            std::string path = names.triggerName(i);
            if( path != "generation_step" &&  path !=  "digitisation_step" && path != "L1simulation_step"  &&  path != "digi2raw_step"  
                    &&  path != "HLT_ReducedIterativeTracking_v1"   &&  path != "HLT_ZeroBias_v1"  &&  path != "HLT_Physics_v1"  &&  path != "HLTriggerFinalPath")
            hltpaths.push_back(names.triggerName(i)); 
        }
    }


    //jets
    for (const pat::Jet &j : *jets) {
        jetpt.push_back(j.pt());
        jeteta.push_back(j.eta());
        jetphi.push_back(j.phi());
        pujetid.push_back(j.userFloat("pileupJetId:fullDiscriminant"));
        jetbptag.push_back(j.bDiscriminator("jetBProbabilityBJetTags")); 
        jetcsvtag.push_back(j.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags"));
        bjet.push_back(utilities::isbJet(j));
        jecfactor.push_back(j.jecFactor("Uncorrected"));
        jetlooseID.push_back(utilities::jetID(j));
        pujetetaid.push_back(utilities::pujetetaid(j));
    }

    if(jets->size() > 1){
        TLorentzVector jjp4 = TLorentzVector(jets->at(0).px(), jets->at(0).py(), jets->at(0).pz(), jets->at(0).energy()) +  TLorentzVector(jets->at(1).px(), jets->at(1).py(), jets->at(1).pz(), jets->at(1).energy());
        mjj = jjp4.M();
        ptjj = jjp4.Pt();
        phijj = jjp4.Phi();
        deta = fabs( jets->at(0).eta() - jets->at(1).eta() );
        dphi = fabs( jets->at(0).phi() - jets->at(1).phi() );

    }
    if(jets->size() > 2){
        for(std::vector<int>::size_type i = 2; i != jets->size(); i++){
            if ( (*jets)[i].pt() < 30 || ! (*jets)[i].userFloat("pileupJetId:fullDiscriminant") || !utilities::pujetetaid( (*jets)[i]) )
                continue;
            if( (*jets)[i].eta() > std::min((*jets)[0].eta(),(*jets)[1].eta() ) &&  (*jets)[i].eta() < std::max((*jets)[0].eta(),(*jets)[1].eta() )   ){
                njetingap+=1;         
            }
        
        }
    }


    isZtt = false;
    isZmt = false;
    isZet = false;
    isZee = false;
    isZmm = false;
    isZem = false;
    isZEE = false;
    isZMM = false;
    isZLL = false;

    /*
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
    */
      //  TLorentzVector J1(jet1->px, jet1->py, jet1->pz, jet1->energy);
      //  TLorentzVector J2(jet2->px, jet2->py, jet2->pz, jet2->energy);
      //  TLorentzVector dijet = J1 + J2;
      //  Float_t mjj = dijet.M();

    for (const pat::CompositeCandidate &lP : *pairs){
             
        const reco::Candidate * l1, *l2; 
        l1 = lP.daughter(0); l2 = lP.daughter(1);
        const reco::Candidate *met = lP.daughter(2);
        if (lP.daughter(2)->hasMasterClone())
            const edm::Ref<reco::PFMETCollection> met = lP.daughter(2)->masterClone().castTo<edm::Ref<reco::PFMETCollection>>();
        else
            met = dynamic_cast<const pat::MET*>(lP.daughter(2));



        svfit.push_back(lP.userFloat("SVfitMass"));

        metpx.push_back(met->px());
        //std::cout << "met->px(): " << met->px() << " ; METx: " << lP.userFloat("MEt_px") << " : ";//  << met->electronEtFraction() << std::endl;
        metpt.push_back(met->pt());
        metphi.push_back(met->phi());
//        metsumEt.push_back(met->sumEt());


        mvacov00.push_back(lP.userFloat("MEt_cov00"));
        mvacov01.push_back(lP.userFloat("MEt_cov01"));
        mvacov10.push_back(lP.userFloat("MEt_cov10"));
        mvacov11.push_back(lP.userFloat("MEt_cov11"));


        diq.push_back((int)(l1->charge() * l2->charge()));
        pth.push_back( (l1->p4() + l2->p4() + met->p4()).pt()   );
        ptvis.push_back( (l1->p4() + l2->p4()).pt()   );
        m_vis.push_back((l1->p4() + l2->p4()).mass()  );


        std::string temp;
        for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
            obj.unpackPathNames(names);
            if(deltaR( obj.triggerObject().p4(), l1->p4()) < 0.5)
                for (unsigned h = 0; h < obj.filterLabels().size(); ++h) temp += (obj.filterLabels()[h] + "/");
                
        }
        hltmatch.push_back(temp);

        if(jets->size() > 1){
        
            TLorentzVector jjp4 = TLorentzVector(jets->at(0).px(), jets->at(0).py(), jets->at(0).pz(), jets->at(0).energy()) 
                +  TLorentzVector(jets->at(1).px(), jets->at(1).py(), jets->at(1).pz(), jets->at(1).energy());
            hdijetphi.push_back(deltaPhi(jjp4.Phi(), (l1->p4() + l2->p4()).phi() ));
            visjeteta.push_back(std::min( fabs(jets->at(0).eta() - (l1->p4() + l2->p4()).eta()), 
                        fabs(jets->at(1).eta() - (l1->p4() + l2->p4()).eta())));
        
        }

        if (l1->isMuon()){
            if (l2->isMuon()){
                //std::cout << "muon-muon channel" <<std::endl;
                const pat::Muon *mu  = dynamic_cast<const pat::Muon*>(lP.daughter(0)->masterClone().get());
                const pat::Muon *mu_ = dynamic_cast<const pat::Muon*>(lP.daughter(1)->masterClone().get());
                fillmuon(mu, PV, met);
                fillmuon(mu_, PV, met);
            }
            else if(lP.daughter(1)->isElectron()){
                //std::cout << "muon-electron channel" <<std::endl;
                const pat::Muon *mu  = dynamic_cast<const pat::Muon*>(lP.daughter(0)->masterClone().get());
            //    const pat::Electron *el = dynamic_cast<const pat::Electron*>(lP.daughter(1)->masterClone().get());
                fillmuon(mu, PV, met);
            }
            else {
                //std::cout << "muon-tau channell" <<std::endl;
                const pat::Muon *mu  = dynamic_cast<const pat::Muon*>(lP.daughter(0)->masterClone().get());
                const pat::Tau *tau  = dynamic_cast<const pat::Tau*>(lP.daughter(1)->masterClone().get());
                fillmuon(mu, PV, met);
                filltau(tau, met);
            }
        }
        else if(lP.daughter(0)->isElectron()){
            if (lP.daughter(1)->isMuon()){
                //std::cout << "muon-electron channel" <<std::endl;
             //   const pat::Electron *el = dynamic_cast<const pat::Electron*>(lP.daughter(0)->masterClone().get());
                const pat::Muon *mu = dynamic_cast<const pat::Muon*>(lP.daughter(1)->masterClone().get());
                fillmuon(mu, PV, met);
            }
            else if(lP.daughter(1)->isElectron()){
                //std::cout << "electron-electron channel" <<std::endl;
             //   const pat::Electron *el = dynamic_cast<const pat::Electron*>(lP.daughter(0)->masterClone().get());
              //  const pat::Electron *el_ = dynamic_cast<const pat::Electron*>(lP.daughter(1)->masterClone().get());
            }
            else {
                //std::cout << "electron-tau channel" <<std::endl;
             //   const pat::Electron *el = dynamic_cast<const pat::Electron*>(lP.daughter(0)->masterClone().get());
                const pat::Tau *tau  = dynamic_cast<const pat::Tau*>(lP.daughter(1)->masterClone().get());
                filltau(tau, met);
            }
        }
        else {
            if (lP.daughter(1)->isMuon()){
                //std::cout << "muon-tau channel" <<std::endl;
                const pat::Tau *tau  = dynamic_cast<const pat::Tau*>(lP.daughter(0)->masterClone().get());
                const pat::Muon *mu = dynamic_cast<const pat::Muon*>(lP.daughter(1)->masterClone().get());
                fillmuon(mu, PV, met);
                filltau(tau, met);
            }
            else if(lP.daughter(1)->isElectron()){
                //std::cout << "electron-tau channel" <<std::endl;
                const pat::Tau *tau  = dynamic_cast<const pat::Tau*>(lP.daughter(0)->masterClone().get());
            //    const pat::Electron *el = dynamic_cast<const pat::Electron*>(lP.daughter(1)->masterClone().get());
                filltau(tau, met);
            }
            else {
                //std::cout << "tau-tau channel" <<std::endl;
                const pat::Tau *tau  = dynamic_cast<const pat::Tau*>(lP.daughter(0)->masterClone().get());
                const pat::Tau *tau_  = dynamic_cast<const pat::Tau*>(lP.daughter(1)->masterClone().get());
                filltau(tau, met);
                filltau(tau_, met);
            }
        }


        
    }
    if(mupt.size()>0 && taupt.size()>0)
        mt->Fill();


}


// ------------ method called once each job just before starting event loop  ------------
void maxi2ntuples::beginJob(){

}

// ------------ method called once each job just after ending the event loop  ------------
void 
maxi2ntuples::endJob() 
{
}

void maxi2ntuples::fillmuon(const pat::Muon* mu,  const reco::Vertex& PV, const reco::Candidate* met){

    muon.push_back(std::make_pair(1,1));
    std::cout << "pt: "<< mu->pt() << " \n";
    mupt.push_back(mu->pt());
    muphi.push_back(mu->phi());
    mueta.push_back(mu->eta());
    mum.push_back(mu->mass());
    muq.push_back(mu->charge());
    mud0.push_back(mu->innerTrack()->dxy( PV.position()));
    mudz.push_back(mu->innerTrack()->dz(PV.position()));
    mumt.push_back(sqrt(pow((mu->p4()).pt() + (met->p4()).pt(),2) - pow((mu->p4() + met->p4()).pt(),2)));
    isLooseMuon.push_back(mu->isLooseMuon());
    isTightMuon.push_back(mu->isTightMuon(PV));
    isHighPtMuon.push_back( mu->isHighPtMuon(PV));
    isMediumMuon.push_back(utilities::heppymuonID(*mu, "POG_ID_Medium"));
    isTightnovtxMuon.push_back(utilities::heppymuonID(*mu, "POG_ID_TightNoVtx"));
    muiso.push_back( utilities::relIso(*mu, 0.5));



}


void maxi2ntuples::filltau(const pat::Tau* tau, const reco::Candidate*met ){


    taupt.push_back(tau->pt());
    tauphi.push_back(tau->phi());
    taueta.push_back(tau->eta());
    taum.push_back(tau->mass());
    tauq.push_back(tau->charge());
    taumt.push_back(sqrt(pow((tau->p4()).pt() + (met->p4()).pt(),2) - pow((tau->p4() + met->p4()).pt(),2)));
    decayModeFinding.push_back(tau->tauID("decayModeFinding"));
//      decayModeFindingOldDMs.push_back(tau->tauID("decayModeFindingOldDMs"));
    decayModeFindingNewDMs.push_back(tau->tauID("decayModeFindingNewDMs"));
    byCombinedIsolationDeltaBetaCorrRaw3Hits.push_back(tau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"));
    lbyCombinedIsolationDeltaBetaCorrRaw3Hits.push_back(tau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")); 
    mbyCombinedIsolationDeltaBetaCorrRaw3Hits.push_back(tau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits")); 
    tbyCombinedIsolationDeltaBetaCorrRaw3Hits.push_back(tau->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits")); 
    chargedIsoPtSum.push_back(tau->tauID("chargedIsoPtSum"));
    neutralIsoPtSum.push_back(tau->tauID("neutralIsoPtSum"));
    puCorrPtSum.push_back(tau->tauID("puCorrPtSum"));
    againstMuonLoose3.push_back(tau->tauID("againstMuonLoose3"));
    againstMuonTight3.push_back(tau->tauID("againstMuonTight3"));
    againstElectronVLooseMVA5.push_back(tau->tauID("againstElectronVLooseMVA5"));
    againstElectronLooseMVA5.push_back(tau->tauID("againstElectronLooseMVA5"));
    againstElectronMediumMVA5.push_back(tau->tauID("againstElectronMediumMVA5"));
    againstElectronTightMVA5.push_back(tau->tauID("againstElectronTightMVA5"));
    againstElectronVTightMVA5.push_back(tau->tauID("againstElectronVTightMVA5"));

    byIsolationMVA3newDMwoLTraw.push_back(tau->tauID("byIsolationMVA3newDMwoLTraw"));  
    byIsolationMVA3oldDMwoLTraw.push_back(tau->tauID("byIsolationMVA3oldDMwoLTraw"));  
    byIsolationMVA3newDMwLTraw.push_back(tau->tauID("byIsolationMVA3newDMwLTraw"));  
    byIsolationMVA3oldDMwLTraw.push_back(tau->tauID("byIsolationMVA3oldDMwLTraw"));



}


// ------------ method called when starting to processes a run  ------------
/*
void 
maxi2ntuples::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
maxi2ntuples::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
maxi2ntuples::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
maxi2ntuples::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
maxi2ntuples::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(maxi2ntuples);


