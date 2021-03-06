//
// Package:    m2n/maxi2ntuples
// Class:      tautau
// 
/**\class tautau tautau.cc m2n/maxi2ntuples/plugins/tautau.cc

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
//
// class declaration
//

class tautau : public edm::EDAnalyzer {
   public:
      explicit tautau(const edm::ParameterSet&);
      ~tautau();

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



    TTree * t;                                                                                                //<--------------------------------------------------------------------------------------------------------------
    ///////////////////////////// P A I R     S P E C I F I C: ////////////////////////////////
    std::vector<float> 
        tau_pt, tau_phi, tau_eta, tau_m, tau_q, tau_mt,
        taupt, tauphi, taueta, taum, tauq, taumt,
        svfit, metpx, metpt, metphi, metsumEt, tau_iso, decayModeFinding;
    std::vector<int> diq;
    std::vector<float>  // decayModeFindingOldDMs, 
        decayModeFindingNewDMs, byCombinedIsolationDeltaBetaCorrRaw3Hits, lbyCombinedIsolationDeltaBetaCorrRaw3Hits, mbyCombinedIsolationDeltaBetaCorrRaw3Hits, tbyCombinedIsolationDeltaBetaCorrRaw3Hits, chargedIsoPtSum, neutralIsoPtSum,
        puCorrPtSum, againstMuonLoose3, againstMuonTight3, againstElectronVLooseMVA5, againstElectronLooseMVA5, againstElectronMediumMVA5, againstElectronTightMVA5, againstElectronVTightMVA5,
        byIsolationMVA3newDMwoLTraw, byIsolationMVA3oldDMwoLTraw, byIsolationMVA3newDMwLTraw, byIsolationMVA3oldDMwLTraw,
        pth, ptvis, m_vis,
        mvacov00, mvacov01, mvacov10, mvacov11,
        hdijetphi, visjeteta; 

    ///////////////////////////// E V E N T     S P E C I F I C: ////////////////////////////////
    //trigger paths
    std::vector<std::string> hltpaths;
    //jets
    std::vector<float> jetpt, jeteta, jetphi, pujetid, jetbptag, jetcsvtag, bjet, jecfactor, jetlooseID, pujetetaid;
    float  run, lumi, evt, npv, npu, mjj, ptjj, phijj, deta, dphi, njetingap;
    unsigned short paircount;//, tau__dz ; 
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
tautau::tautau(const edm::ParameterSet& iConfig):
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
    PileupSummaryInfoToken_(consumes<PileupSummaryInfoCollection>(iConfig.getParameter<edm::InputTag>("pileupinfo")))
{
   //now do what ever initialization is needed

    edm::Service<TFileService> theFileService;
    t = theFileService->make<TTree>("ntuple", "tautau");
    t->Branch("tau_pt", &tau_pt);
    t->Branch("tau_phi", &tau_phi);
    t->Branch("tau_eta", &tau_eta);
    t->Branch("tau_m", &tau_m);
    t->Branch("tau_q", &tau_q);
    t->Branch("tau_mt", &tau_mt);
    t->Branch("taupt", &taupt);
    t->Branch("tauphi", &tauphi);
    t->Branch("taueta", &taueta);
    t->Branch("taum", &taum);
    t->Branch("tauq", &tauq);
    t->Branch("taumt", &taumt);
    t->Branch("svfit", &svfit);
    t->Branch("metpx", &metpx);
    t->Branch("metpt", &metpt);
    t->Branch("metphi", &metphi);                                                                        //<-----------------------------------------------
    t->Branch("metsumEt", &metsumEt);
    t->Branch("hltpaths", &hltpaths);
    t->Branch("decayModeFinding", &decayModeFinding);
//    t->Branch("decayModeFindingOldDMs", &decayModeFindingOldDMs);
    t->Branch("decayModeFindingNewDMs", &decayModeFindingNewDMs);
    t->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits", &byCombinedIsolationDeltaBetaCorrRaw3Hits );
    t->Branch("lbyCombinedIsolationDeltaBetaCorrRaw3Hits", &lbyCombinedIsolationDeltaBetaCorrRaw3Hits );
    t->Branch("mbyCombinedIsolationDeltaBetaCorrRaw3Hits", &mbyCombinedIsolationDeltaBetaCorrRaw3Hits );
    t->Branch("tbyCombinedIsolationDeltaBetaCorrRaw3Hits", &tbyCombinedIsolationDeltaBetaCorrRaw3Hits );
    t->Branch("chargedIsoPtSum", &chargedIsoPtSum );
    t->Branch("neutralIsoPtSum", &neutralIsoPtSum );
    t->Branch("puCorrPtSum", &puCorrPtSum );
    t->Branch("againstMuonLoose3", &againstMuonLoose3 );
    t->Branch("againstMuonTight3", &againstMuonTight3 );
    t->Branch("againstElectronVLooseMVA5", &againstElectronVLooseMVA5 );
    t->Branch("againstElectronLooseMVA5", &againstElectronLooseMVA5 );
    t->Branch("againstElectronMediumMVA5", &againstElectronMediumMVA5 );
    t->Branch("againstElectronTightMVA5", &againstElectronTightMVA5);
    t->Branch("againstElectronVTightMVA5", &againstElectronVTightMVA5);
    t->Branch("byIsolationMVA3newDMwoLTraw",&byIsolationMVA3newDMwoLTraw); 
    t->Branch("byIsolationMVA3oldDMwoLTraw",&byIsolationMVA3oldDMwoLTraw); 
    t->Branch("byIsolationMVA3newDMwLTraw",&byIsolationMVA3newDMwLTraw); 
    t->Branch("byIsolationMVA3oldDMwLTraw",&byIsolationMVA3oldDMwLTraw); 
    t->Branch("diq", &diq);

    t->Branch("pth", &pth);
    t->Branch("ptvis", &ptvis);
    t->Branch("m_vis", &m_vis); 

    t->Branch("mvacov00",&mvacov00); 
    t->Branch("mvacov01",&mvacov01); 
    t->Branch("mvacov10",&mvacov10); 
    t->Branch("mvacov11",&mvacov11); 

    t->Branch("hdijetphi",&hdijetphi); 
    t->Branch("visjeteta",&visjeteta); 

    t->Branch("tau_iso", &tau_iso);
    t->Branch("jetpt", &jetpt);
    t->Branch("jeteta", &jeteta);
    t->Branch("jetphi", &jetphi);
    t->Branch("pujetid", &pujetid);
    t->Branch("jetbptag", &jetbptag);
    t->Branch("jetcsvtag", &jetcsvtag);
    t->Branch("jecfactor", &jecfactor);
    t->Branch("jetlooseID", &jetlooseID);
    t->Branch("pujetetaid", &pujetetaid);
    t->Branch("bjet", &bjet);
    t->Branch("mjj", &mjj);
    t->Branch("ptjj", &ptjj);
    t->Branch("phijj", &phijj);
    t->Branch("deta", &deta);
    t->Branch("dphi", &dphi);
    t->Branch("njetingap", &njetingap);

    t->Branch("run", &run);
    t->Branch("lumi", &lumi);
    t->Branch("evt", &evt);
    t->Branch("npv", &npv);
    t->Branch("npu", &npu);
    t->Branch("paircount", &paircount);
  
    t->Branch("isZtt", &isZtt);
    t->Branch("isZmt", &isZmt);
    t->Branch("isZet", &isZet);
    t->Branch("isZee", &isZee);
    t->Branch("isZmm", &isZmm);
    t->Branch("isZem", &isZem);
    t->Branch("isZEE", &isZEE);
    t->Branch("isZMM", &isZMM);
    t->Branch("isZLL", &isZLL);
    t->Branch("isFake", &isFake);
    t->Branch("nup", &nup);
}


tautau::~tautau()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
tautau::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
 
    tau_pt.clear(); tau_phi.clear(); tau_eta.clear(); tau_m.clear(); tau_q.clear(); tau_mt.clear();
    taupt.clear(); tauphi.clear(); taueta.clear(); taum.clear(); tauq.clear(); taumt.clear();
    svfit.clear(); 
    metpx.clear(); metpt.clear(); metphi.clear(); metsumEt.clear();                                    //<----------------------------------------------------
    hltpaths.clear();
  //  decayModeFindingOldDMs.clear(); 
    decayModeFindingNewDMs.clear();
    diq.clear(); pth.clear(); ptvis.clear(); m_vis.clear();     

    mvacov00.clear();mvacov01.clear();mvacov10.clear();mvacov11.clear();
    hdijetphi.clear(); visjeteta.clear();

    tau_iso.clear(); decayModeFinding.clear();
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
//    const reco::Vertex &PV = vertices->front();

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

    edm::Handle<LHEEventProduct> evnt;
    iEvent.getByToken(lheprodToken_, evnt);

    edm::Handle<PileupSummaryInfoCollection> genPileUpInfos;
    iEvent.getByToken(PileupSummaryInfoToken_, genPileUpInfos);

    for(const PileupSummaryInfo &pusi : *genPileUpInfos){
        int bx = pusi.getBunchCrossing();
        int nPU = pusi.getPU_NumInteractions();
        if(bx == 0)
            npu = nPU;
    }

    const lhef::HEPEUP hepeup_ = evnt->hepeup();
    nup =  hepeup_.NUP;
    
    
    paircount = pairs->size();
    run = iEvent.id().run();
    lumi = iEvent.luminosityBlock();
    evt = iEvent.id().event(); 
    npv =  vertices->size();
     for (const pat::CompositeCandidate &lP : *pairs){
       
        const pat::Tau *tau = 0, *tau_ = 0; 
        tau = dynamic_cast<const pat::Tau*>(lP.daughter(0)->masterClone().get());    
        tau_ = dynamic_cast<const pat::Tau*>(lP.daughter(1)->masterClone().get());    

        const reco::Candidate *met = lP.daughter(2);
        if (lP.daughter(2)->hasMasterClone())
            const edm::Ref<reco::PFMETCollection> met = lP.daughter(2)->masterClone().castTo<edm::Ref<reco::PFMETCollection>>();
        else
            met = dynamic_cast<const pat::MET*>(lP.daughter(2));

        taupt.push_back(tau->pt());
        tauphi.push_back(tau->phi());
        taueta.push_back(tau->eta());
        taum.push_back(tau->mass());
        tauq.push_back(tau->charge());
        taumt.push_back(sqrt(pow((tau->p4()).pt() + (met->p4()).pt(),2) - pow((tau->p4() + met->p4()).pt(),2)));

        tau_pt.push_back(tau_->pt());
        tau_phi.push_back(tau_->phi());
        tau_eta.push_back(tau_->eta());
        tau_m.push_back(tau_->mass());
        tau_q.push_back(tau_->charge());
        tau_mt.push_back(sqrt(pow((tau_->p4()).pt() + (met->p4()).pt(),2) - pow((tau_->p4() + met->p4()).pt(),2)));

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

        diq.push_back((int)(tau_->charge() * tau->charge()));
        pth.push_back( (tau_->p4() + tau->p4() + met->p4()).pt()   );
        ptvis.push_back( (tau_->p4() + tau->p4()).pt()   );
        m_vis.push_back((tau_->p4() + tau->p4()).mass()  );



        if(jets->size() > 1){
        
            TLorentzVector jjp4 = TLorentzVector(jets->at(0).px(), jets->at(0).py(), jets->at(0).pz(), jets->at(0).energy()) 
                                + TLorentzVector(jets->at(1).px(), jets->at(1).py(), jets->at(1).pz(), jets->at(1).energy());
            hdijetphi.push_back(deltaPhi(jjp4.Phi(), (tau_->p4() + tau->p4()).phi() ));
            visjeteta.push_back(std::min( fabs(jets->at(0).eta() - (tau_->p4() + tau->p4()).eta()), 
                        fabs(jets->at(1).eta() - (tau_->p4() + tau->p4()).eta())));
        }

     }


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

      //  TLorentzVector J1(jet1->px, jet1->py, jet1->pz, jet1->energy);
      //  TLorentzVector J2(jet2->px, jet2->py, jet2->pz, jet2->energy);
      //  TLorentzVector dijet = J1 + J2;
      //  Float_t mjj = dijet.M();

     t->Fill();


}


// ------------ method called once each job just before starting event loop  ------------
void tautau::beginJob(){

}

// ------------ method called once each job just after ending the event loop  ------------
void 
tautau::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
tautau::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
tautau::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
tautau::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
tautau::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
tautau::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(tautau);
