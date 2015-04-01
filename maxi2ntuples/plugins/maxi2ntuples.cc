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


#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
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

class maxi2ntuples : public edm::EDAnalyzer {
   public:
      explicit maxi2ntuples(const edm::ParameterSet&);
      ~maxi2ntuples();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


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


    TTree * t;                                                                                                //<--------------------------------------------------------------------------------------------------------------
    ///////////////////////////// P A I R     S P E C I F I C: ////////////////////////////////
    std::vector<float> mupt, taupt, svfit, metpx, metpt, metphi, metsumEt, mucombreliso, decayModeFinding;
    std::vector<int> charge;
    std::vector<float>  // decayModeFindingOldDMs, 
        decayModeFindingNewDMs, cidbcr3h, lcidbcr3h, mcidbcr3h, tcidbcr3h, chargedIsoPtSum, neutralIsoPtSum,
        puCorrPtSum, againstMuonLoose3, againstMuonTight3, againstElectronVLooseMVA5, againstElectronLooseMVA5, againstElectronMediumMVA5; 
    std::vector<std::string> hltmatch;

    ///////////////////////////// E V E N T     S P E C I F I C: ////////////////////////////////
    //trigger paths
    std::vector<std::string> hltpaths;
    //jets
    std::vector<float> jetpt, pujetid, jetbptag, jetcsvtag;
    float  NPV, mjj, deta;
    unsigned short paircount;//, mu_dz ; 

//      short isLooseMuon, isTightMuon;
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

{
   //now do what ever initialization is needed

    edm::Service<TFileService> theFileService;
    t = theFileService->make<TTree>("ntuple", "tautau");
    t->Branch("mupt", &mupt);
    t->Branch("taupt", &taupt);
    t->Branch("svfit", &svfit);
    t->Branch("metpx", &metpx);
    t->Branch("metpt", &metpt);
    t->Branch("metphi", &metphi);                                                                        //<-----------------------------------------------
    t->Branch("metsumEt", &metsumEt);
    t->Branch("hltmatch", &hltmatch);
    t->Branch("hltpaths", &hltpaths);
    t->Branch("decayModeFinding", &decayModeFinding);
//    t->Branch("decayModeFindingOldDMs", &decayModeFindingOldDMs);
    t->Branch("decayModeFindingNewDMs", &decayModeFindingNewDMs);
    t->Branch("cidbcr3h", &cidbcr3h );
    t->Branch("lcidbcr3h", &lcidbcr3h );
    t->Branch("mcidbcr3h", &mcidbcr3h );
    t->Branch("tcidbcr3h", &tcidbcr3h );
    t->Branch("chargedIsoPtSum", &chargedIsoPtSum );
    t->Branch("neutralIsoPtSum", &neutralIsoPtSum );
    t->Branch("puCorrPtSum", &puCorrPtSum );
    t->Branch("againstMuonLoose3", &againstMuonLoose3 );
    t->Branch("againstMuonTight3", &againstMuonTight3 );
    t->Branch("againstElectronVLooseMVA5", &againstElectronVLooseMVA5 );
    t->Branch("againstElectronLooseMVA5", &againstElectronLooseMVA5 );
    t->Branch("againstElectronMediumMVA5", &againstElectronMediumMVA5 );
    t->Branch("charge", &charge);
    t->Branch("mucombreliso", &mucombreliso);
    t->Branch("jetpt", &jetpt);
    t->Branch("pujetid", &pujetid);
    t->Branch("jetbptag", &jetbptag);
    t->Branch("jetcsvtag", &jetcsvtag);
    t->Branch("mjj", &mjj);
    t->Branch("deta", &deta);


//    t->Branch("mu_dz", &mu_dz);
//    t->Branch("isLooseMuon,", &isLooseMuon);
//    t->Branch("isTightMuon", &isTightMuon);
    t->Branch("NPV", &NPV);
    t->Branch("paircount", &paircount);

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
 
    mupt.clear(); taupt.clear(); svfit.clear(); 
    metpx.clear(); metpt.clear(); metphi.clear(); metsumEt.clear();                                    //<----------------------------------------------------
    hltmatch.clear();  hltpaths.clear();
  //  decayModeFindingOldDMs.clear(); 
    decayModeFindingNewDMs.clear();
    charge.clear(); mucombreliso.clear(); decayModeFinding.clear();
    cidbcr3h.clear(); lcidbcr3h.clear();  mcidbcr3h.clear(); tcidbcr3h.clear(); chargedIsoPtSum.clear();
    neutralIsoPtSum.clear(); puCorrPtSum.clear(); againstMuonLoose3.clear(); againstMuonTight3.clear(); 
    againstElectronVLooseMVA5.clear(); againstElectronLooseMVA5.clear(); againstElectronMediumMVA5.clear();
    jetpt.clear(); pujetid.clear(); jetbptag.clear(); jetcsvtag.clear();
    mjj = -1; deta = 0;



    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);
    if (vertices->empty()) return; // skip the event if no PV found
    //const reco::Vertex &PV = vertices->front();
    NPV =  vertices->size();

    /*
     edm::Handle<pat::MuonCollection> muons;
     iEvent.getByToken(muonToken_, muons);
     for (const pat::Muon &mu : *muons) {
        if (mu.pt() < 5 || !mu.isLooseMuon()) continue;
//        mu_pt = mu.pt();
        mu_dz = mu.muonBestTrack()->dz(PV.position());
        isLooseMuon =  mu.isLooseMuon();
        isTightMuon =  mu.isTightMuon(PV);
     }
*/


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


     paircount = pairs->size();
     for (const pat::CompositeCandidate &lP : *pairs){
        const pat::Muon *mu = 0; 
        const pat::Tau *tau = 0;    
        const pat::MET *met = 0;

        if(lP.daughter("leptonTwo")->isMuon()){
            mu = dynamic_cast<const pat::Muon*>(lP.daughter("leptonTwo"));    
            tau = dynamic_cast<const pat::Tau*>(lP.daughter("leptonOne"));    
            met =  dynamic_cast<const pat::MET*>(lP.daughter("met"));
        }
        else{
          //  std::cout << lP.daughter("leptonOne")->isMuon() << "; " << lP.daughter("leptonOne")->isElectron() << std::endl;
            mu = dynamic_cast<const pat::Muon*>(lP.daughter("leptonOne"));    
            tau = dynamic_cast<const pat::Tau*>(lP.daughter("leptonTwo"));    
            met =  dynamic_cast<const pat::MET*>(lP.daughter("met"));
        }
        mupt.push_back(mu->pt());
        taupt.push_back(tau->pt());
        svfit.push_back(lP.userFloat("SVfitMass"));
        metpx.push_back(met->px());
        metpt.push_back(met->pt());
        metphi.push_back(met->phi());
        metsumEt.push_back(met->sumEt());


        decayModeFinding.push_back(tau->tauID("decayModeFinding"));
  //      decayModeFindingOldDMs.push_back(tau->tauID("decayModeFindingOldDMs"));
        decayModeFindingNewDMs.push_back(tau->tauID("decayModeFindingNewDMs"));
        cidbcr3h.push_back(tau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"));
        lcidbcr3h.push_back(tau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")); 
        mcidbcr3h.push_back(tau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits")); 
        tcidbcr3h.push_back(tau->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits")); 
        chargedIsoPtSum.push_back(tau->tauID("chargedIsoPtSum"));
        neutralIsoPtSum.push_back(tau->tauID("neutralIsoPtSum"));
        puCorrPtSum.push_back(tau->tauID("puCorrPtSum"));
        againstMuonLoose3.push_back(tau->tauID("againstMuonLoose3"));
        againstMuonTight3.push_back(tau->tauID("againstMuonTight3"));
        againstElectronVLooseMVA5.push_back(tau->tauID("againstElectronVLooseMVA5"));
        againstElectronLooseMVA5.push_back(tau->tauID("againstElectronLooseMVA5"));
        againstElectronMediumMVA5.push_back(tau->tauID("againstElectronMediumMVA5"));

        charge.push_back((int)(mu->charge() * tau->charge()));
        mucombreliso.push_back( utilities::relIso(*mu, 0.5));

        std::string temp;
        for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
            obj.unpackPathNames(names);
            if(deltaR( obj.triggerObject().p4(), mu->p4()) < 0.5)
                for (unsigned h = 0; h < obj.filterLabels().size(); ++h) temp += (obj.filterLabels()[h] + "/");
                
        }
        hltmatch.push_back(temp);

     }


    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
        if(triggerBits->accept(i)){
            std::string path = names.triggerName(i);
            if( path != "generation_step" &&  path !=  "digitisation_step" && path != "L1simulation_step"  &&  path != "digi2raw_step"  &&  path != "HLT_ReducedIterativeTracking_v1"   &&  path != "HLT_ZeroBias_v1"  &&  path != "HLT_Physics_v1"  &&  path != "HLTriggerFinalPath")
            hltpaths.push_back(names.triggerName(i)); 
        }
    }


    //jets
    for (const pat::Jet &j : *jets) {
        jetpt.push_back(j.pt());
        pujetid.push_back(j.userFloat("pileupJetId:fullDiscriminant"));
        jetbptag.push_back(j.bDiscriminator("jetBProbabilityBJetTags")); 
        jetcsvtag.push_back(j.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags"));
    }
    if(jets->size() > 1){
        mjj = (TLorentzVector(jets->at(0).px(), jets->at(0).py(), jets->at(0).pz(), jets->at(0).energy()) +  TLorentzVector(jets->at(1).px(), jets->at(1).py(), jets->at(1).pz(), jets->at(1).energy())).M();
        deta = fabs( jets->at(0).eta() - jets->at(1).eta() );

    }


      //  TLorentzVector J1(jet1->px, jet1->py, jet1->pz, jet1->energy);
      //  TLorentzVector J2(jet2->px, jet2->py, jet2->pz, jet2->energy);
      //  TLorentzVector dijet = J1 + J2;
      //  Float_t mjj = dijet.M();

     t->Fill();


}


// ------------ method called once each job just before starting event loop  ------------
void maxi2ntuples::beginJob(){

}

// ------------ method called once each job just after ending the event loop  ------------
void 
maxi2ntuples::endJob() 
{
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
