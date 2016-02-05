#include "m2n/maxi2ntuples/plugins/ntuple.h"

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
ntuple::ntuple(const edm::ParameterSet& iConfig):
  beamSpotToken(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))), 
  vtxToken(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  scoresToken(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("vertexScores"))),
  muonToken(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronToken(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  tauToken(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  photonToken(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  jetToken(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  fatjetToken(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
  metToken(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
  pairToken(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("pairs"))),
  triggerBitsToken(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjectsToken(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  triggerPrescalesToken(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  prunedGenToken(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
  packedGenToken(consumes<pat::PackedGenParticleCollection>(iConfig.getParameter<edm::InputTag>("packedGenParticles"))),
  lheprodToken(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheprod"))),
  PileupSummaryInfoToken(consumes<PileupSummaryInfoCollection>(iConfig.getParameter<edm::InputTag>("pileupinfo"))),
  candsToken(consumes<edm::View<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("src"))),  
  sample(iConfig.getParameter<int>("sample")){
  ;
}
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
ntuple::~ntuple(){
 
    delete wevent;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void ntuple::beginJob(){

   edm::Service<TFileService> theFileService;
    
   eventTree = theFileService->make<TTree>("eventTree", "eventTree");
   wevent = 0;

   eventTree->Branch("wevent", &wevent);
   eventTree->Branch("wtau", &wtaucollection);
   eventTree->Branch("wmu", &wmucollection);
   eventTree->Branch("welectron", &welectroncollection);
   eventTree->Branch("wpair", &wpaircollection);
   eventTree->Branch("wmet", &wmetcollection);
   eventTree->Branch("wjet", &wjetcollection);

   events = theFileService->make<TH1F>("hvar","hvar title",10,0,10);

}
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
void ntuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  //const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  //events->Fill(1.,1.);
  //double weightevt=genEvt->weight(); 
  //events->Fill(2., weightevt);
  
  using namespace edm;
  
  clearNtuple();  
  getCollections(iEvent,iSetup);
  
  if (vertices->empty()) return; // skip the event if no PV found
  
  fillEventHeaderData(iEvent, iSetup);
  findPrimaryVertices(iEvent, iSetup);
  refitPV(iEvent, iSetup);

  if(!iEvent.isRealData()){
    getMCCollections(iEvent,iSetup);
    fillGenData();
  }

  eventTree->Fill();
  
  ///Stop processing of not tau pair is found in the event.
  if(!pairs.isValid() || pairs->size()==0) return;
  
  fillTauPairData(iEvent, iSetup);
  fillJetsData();

  eventTree->Fill();
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void ntuple::clearNtuple(){

  wtaucollection.clear();
  wmucollection.clear();
  wpaircollection.clear();
  wmetcollection.clear();
  wjetcollection.clear();
  
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void ntuple::getCollections(const edm::Event& iEvent, const edm::EventSetup& iSetup){

    iEvent.getByToken(beamSpotToken, beamSpot);
    iEvent.getByToken(vtxToken, vertices);
    iEvent.getByToken(scoresToken, scores);
    iEvent.getByToken(candsToken, cands);
    iEvent.getByToken(jetToken, jets);
    iEvent.getByToken(metToken, mets);
    iEvent.getByToken(pairToken, pairs);
    iEvent.getByToken(triggerBitsToken, triggerBits);
    iEvent.getByToken(triggerObjectsToken, triggerObjects);
    iEvent.getByToken(triggerPrescalesToken, triggerPrescales);
   
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void ntuple::getMCCollections(const edm::Event& iEvent, const edm::EventSetup& iSetup){

 try {
      iEvent.getByToken(PileupSummaryInfoToken, genPileUpInfos);
    }
    catch (...) {;}

    try{
      iEvent.getByToken(lheprodToken, lheInfo);
    }
    catch(...){;}

    try{
      iEvent.getByLabel("generator",genInfo);
    }
    catch (...) {;}

    ///FIX ME: check is reading MC
    
    iEvent.getByToken(prunedGenToken, genParticlesPruned);
    iEvent.getByToken(packedGenToken, genParticlesPacked);
    
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void ntuple::fillEventHeaderData(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  wevent->run(iEvent.id().run());
  wevent->lumi(iEvent.luminosityBlock());
  wevent->event(iEvent.id().event());
  wevent->npv(vertices->size());
  wevent->paircount(pairs->size());
  wevent->sample(sample);
  
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void ntuple::fillGenData(){
  
  if(genPileUpInfos.isValid()){
    for(const PileupSummaryInfo &pusi : *genPileUpInfos){
      if(pusi.getBunchCrossing() == 0){
	wevent->npu(pusi.getPU_NumInteractions());
	break;
      }
    }
  }
    
  if(lheInfo.isValid()) wevent->nup(lheInfo->hepeup().NUP);
  if(genInfo.isValid()) wevent->genevtweight(genInfo->weight());

  fillGenTausAndDecayMode();  
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void ntuple::fillGenTausAndDecayMode(){

  std::vector<reco::GenParticle> gentaus; 
  std::vector<reco::GenParticle> gentauleps; 
  std::vector<reco::GenParticle> genleps;
  
  for (const reco::GenParticle &aParticle : *genParticlesPruned){
            
    if(abs(aParticle.pdgId()) == 15){
      gentaus.push_back(aParticle);

      std::cout<<"aParticle.mother(): "<<aParticle.mother()<<std::endl;
      std::cout<<"aParticle.mother().pdgid(): "<<aParticle.mother()<<std::endl;
      
    }
    ///FIX me: write correct code for tau tau decay mode
  }
  
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void ntuple::fillTauPairData(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  
  for (const pat::CompositeCandidate &aPair : *pairs){
        
    const reco::Candidate *leg1 = aPair.daughter(0);
    const reco::Candidate *leg2 = aPair.daughter(1);
    const reco::Candidate *metCand = aPair.daughter(2);
    
    Wpair wpair;
    wpair.svfit(aPair.userFloat("SVfitMass"));
    wpair.diq(leg1->charge() * leg2->charge());
    wpair.pth((leg1->p4() + leg2->p4() + metCand->p4()).pt());
    wpair.ptvis((leg1->p4() + leg2->p4()).pt());
    wpair.m_vis((leg1->p4() + leg2->p4()).mass());
    wpair.ChannelSelector(aPair.userFloat("ChannelSelector"));
    wpair.PATPairSelector(aPair.userFloat("PATPairSelector"));
    wpair.PairBaselineSelection(aPair.userFloat("PairBaselineSelection"));
    wpair.PostSynchSelection(aPair.userFloat("PostSynchSelection"));
    wpair.trigger(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, aPair.userFloat("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL"));
    wpair.trigger(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, aPair.userFloat("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL"));
    wpair.trigger(HLT_IsoMu17_eta2p1_LooseIsoPFTau20, aPair.userFloat("HLT_IsoMu17_eta2p1_LooseIsoPFTau20"));
    wpair.trigger(HLT_IsoMu24_eta2p1, aPair.userFloat("HLT_IsoMu24_eta2p1"));
    wpair.trigger(HLT_IsoMu27, aPair.userFloat("HLT_IsoMu27"));
    wpair.trigger(HLT_IsoMu18, aPair.userFloat("HLT_IsoMu18"));
    wpair.trigger(HLT_IsoMu22, aPair.userFloat("HLT_IsoMu22"));
    wpair.trigger(HLT_IsoMu17_eta2p1, aPair.userFloat("HLT_IsoMu17_eta2p1"));
    wpair.trigger(HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1, aPair.userFloat("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1"));

    wpaircollection.push_back(wpair);

    if(leg1->isMuon()) fillMuonLeg(leg1, metCand, iEvent, iSetup);
    else if(leg1->isElectron()) fillElectronLeg(leg1, metCand, iEvent, iSetup);
    else fillTauLeg(leg1, metCand, iEvent, iSetup);

    if(leg2->isMuon()) fillMuonLeg(leg2, metCand, iEvent, iSetup);
    else if(leg2->isElectron()) fillElectronLeg(leg2, metCand, iEvent, iSetup);
    else fillTauLeg(leg2, metCand, iEvent, iSetup);
    
    fillMET(metCand, aPair);        
  }   
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////  
void ntuple::fillMuonLeg(const reco::Candidate *aCandidate, const reco::Candidate *aMETCandidate,
			 const edm::Event& iEvent, const edm::EventSetup& iSetup){

  const pat::Muon *muon  = dynamic_cast<const pat::Muon*>(aCandidate->masterClone().get());
  
  Wmu wmu;
  wmu.pt(muon->pt());
  wmu.eta(muon->eta());
  wmu.phi(muon->phi());
  wmu.mass(muon->mass());
  wmu.charge(muon->charge());
  setPCAVectors<Wmu>(wmu, muon->bestTrack(), iEvent, iSetup);
  wmu.mt(sqrt(pow((muon->p4()).pt() + (aMETCandidate->p4()).pt(),2) - pow((muon->p4() + aMETCandidate->p4()).pt(),2)));  
  wmu.d0(muon->innerTrack()->dxy((*vertices)[0].position()));
  wmu.dz(muon->innerTrack()->dz((*vertices)[0].position()));
  wmu.isLooseMuon(muon->isLooseMuon());
  wmu.isTightMuon(muon->isTightMuon((*vertices)[0]));
  wmu.isHighPtMuon(muon->isHighPtMuon((*vertices)[0]));
  wmu.isMediumMuon(muon->isMediumMuon());
  wmu.isTightnovtxMuon(utilities::heppymuonID(*muon, "POG_ID_TightNoVtx"));
  wmu.iso(((muon->pfIsolationR03().sumChargedHadronPt
	    + std::max( muon->pfIsolationR03().sumNeutralHadronEt + muon->pfIsolationR03().sumPhotonEt - 0.5 * muon->pfIsolationR03().sumPUPt, 0.0)
	    ) / muon->pt())); 
  
  wmucollection.push_back(wmu);
  
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////  
void ntuple::fillElectronLeg(const reco::Candidate *aCandidate, const reco::Candidate *aMETCandidate,
			     const edm::Event& iEvent, const edm::EventSetup& iSetup){

  const pat::Electron *electron = dynamic_cast<const pat::Electron*>(aCandidate->masterClone().get());
  Welectron welectron;
  welectron.pt(electron->pt());
  welectron.eta(electron->eta());
  welectron.phi(electron->phi());
  welectron.mass(electron->mass());
  welectron.charge(electron->charge());
  welectroncollection.push_back(welectron);
  setPCAVectors<Welectron>(welectron, electron->bestTrack(), iEvent, iSetup);

}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void ntuple::fillTauLeg(const reco::Candidate *aCandidate, const reco::Candidate *aMETCandidate,
			const edm::Event& iEvent, const edm::EventSetup& iSetup){

  const pat::Tau *taon  = dynamic_cast<const pat::Tau*>(aCandidate->masterClone().get());
  
  Wtau wtau;
  wtau.pt(taon->pt());
  wtau.eta(taon->eta());
  wtau.phi(taon->phi());
  wtau.mass(taon->mass());
  wtau.charge(taon->charge());
  wtau.decayMode(taon->decayMode());

  TLorentzVector a4v(taon->leadChargedHadrCand()->p4().px(),
		     taon->leadChargedHadrCand()->p4().py(),
		     taon->leadChargedHadrCand()->p4().pz(),
		     taon->leadChargedHadrCand()->p4().e());
  wtau.leadingTk(a4v);
  setPCAVectors<Wtau>(wtau, taon->leadChargedHadrCand()->bestTrack(), iEvent, iSetup);
  wtau.mt(sqrt(pow((taon->p4()).pt() + (aMETCandidate->p4()).pt(),2) - pow((taon->p4() + aMETCandidate->p4()).pt(),2)));
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
  wtau.tauID(byIsolationMVA3newDMwLTraw, taon->tauID("byIsolationMVA3newDMwLTraw"));
  wtau.tauID(byIsolationMVA3oldDMwLTraw, taon->tauID("byIsolationMVA3oldDMwLTraw"));
  
  wtaucollection.push_back(wtau);
  
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void ntuple::fillMET(const reco::Candidate *aCandidate,
		     const pat::CompositeCandidate & aPair){
  
  Wmet wmet;
  wmet.metpx(aCandidate->px());
  wmet.metpt(aCandidate->pt());
  wmet.metphi(aCandidate->phi());
  wmet.mvacov00(aPair.userFloat("MEt_cov00"));
  wmet.mvacov01(aPair.userFloat("MEt_cov01"));
  wmet.mvacov10(aPair.userFloat("MEt_cov10"));
  wmet.mvacov11(aPair.userFloat("MEt_cov11"));
  wmetcollection.push_back(wmet);
  
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
void ntuple::fillJetsData(){

  for(const pat::Jet &aJet : *jets){
    Wjet wjet;
    wjet.pt(aJet.pt());
    wjet.eta(aJet.eta());
    wjet.phi(aJet.phi());
    wjet.id(aJet.userFloat("pileupJetId:fullDiscriminant"));
    wjet.bptag(aJet.bDiscriminator("jetBProbabilityBJetTags"));
    wjet.csvtag(aJet.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags"));
    wjet.bjet(utilities::isbJet(aJet));
    wjet.jecfactor(aJet.jecFactor("Uncorrected"));
    wjet.jetlooseID(utilities::jetID(aJet));
    wjet.pujetetaid(utilities::pujetetaid(aJet));
    wjetcollection.push_back(wjet);
    }  
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
template<typename T> void ntuple::setPCAVectors(T & aObject, 			   
						const reco::Track *aTrack,
						const edm::Event & iEvent, const edm::EventSetup & iSetup){
  GlobalPoint aPoint(wevent->refitPfPV().X(),
		     wevent->refitPfPV().Y(),
		     wevent->refitPfPV().Z());
  aObject.nPCARefitvx(getPCA(iEvent, iSetup, aTrack, aPoint));
    
  aPoint = GlobalPoint(wevent->thePV().X(),
		       wevent->thePV().Y(),
		       wevent->thePV().Z());
  aObject.nPCA(getPCA(iEvent, iSetup, aTrack, aPoint));
  
  aPoint = GlobalPoint(wevent->genPV().X(),
  		       wevent->genPV().Y(),
  		       wevent->genPV().Z());
  aObject.nPCAGenvx(getPCA(iEvent, iSetup, aTrack, aPoint));
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
bool ntuple::findPrimaryVertices(const edm::Event & iEvent, const edm::EventSetup & iSetup){

 
  if(vertices->size()==0) return false;   //at least one vertex

  TVector3 aPV((*vertices)[0].x(),(*vertices)[0].y(),(*vertices)[0].z());
  wevent->thePV(aPV);

  //Find vertex with highest score with PF (miniAOD like)
  //miniAOD uses PF particles instead of tracks
  size_t iPfVtx=0;
  float score=-1;
  for(size_t iVx=0; iVx<vertices->size(); ++iVx){
    reco::VertexRef vtxPrt(vertices,iVx);   
    if( (*scores)[vtxPrt] > score){
      score = (*scores)[vtxPrt];
      iPfVtx=iVx;
    }
  }

  aPV.SetXYZ((*vertices)[iPfVtx].x(),(*vertices)[iPfVtx].y(),(*vertices)[iPfVtx].z());
  wevent->pfPV(aPV);
  pfPVIndex = iPfVtx;

  return true;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
bool ntuple::refitPV(const edm::Event & iEvent, const edm::EventSetup & iSetup){

  edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder);
 
  TransientVertex transVtx, transVtxNoBS;

  //Get tracks associated wiht pfPV
  reco::TrackCollection pvTracks;
  TLorentzVector aTrack;
  for(size_t i=0; i<cands->size(); ++i){
    if((*cands)[i].charge()==0 || (*cands)[i].vertexRef().isNull()) continue;
    if(!(*cands)[i].bestTrack()) continue;
    ///Skip tracks comming from tau decay.
    /*
    aTrack.SetPxPyPzE((*cands)[i].px(),(*cands)[i].py(),(*cands)[i].pz(),(*cands)[i].energy());
    if(myEvent->recoEvent.piMinus.DeltaR(aTrack)<0.01 ||
       myEvent->recoEvent.piPlus.DeltaR(aTrack)<0.01) continue;       
    */
    unsigned int key = (*cands)[i].vertexRef().key();
    int quality = (*cands)[i].pvAssociationQuality();
    if(key!=pfPVIndex ||
       (quality!=pat::PackedCandidate::UsedInFitTight &&
	quality!=pat::PackedCandidate::UsedInFitLoose)) continue;

    pvTracks.push_back(*((*cands)[i].bestTrack()));
  }
  ///Built transient tracks from tracks.
  std::vector<reco::TransientTrack> transTracks;  
  for(auto iter: pvTracks) transTracks.push_back(transTrackBuilder->build(iter));
  wevent->nTracksInRefit(transTracks.size());

  bool fitOk = false;  
  if(transTracks.size() >= 2 ) {
    AdaptiveVertexFitter avf;
    avf.setWeightThreshold(0.1); //weight per track. allow almost every fit, else --> exception    
    try {
      transVtxNoBS = avf.vertex(transTracks);      
      transVtx = avf.vertex(transTracks, *beamSpot);
      fitOk = true; 
    } catch (...) {
      fitOk = false; 
      std::cout<<"Vtx fit failed!"<<std::endl;
    }
  }
  
  if(fitOk && transVtx.isValid() && transVtxNoBS.isValid()) { 
    ///Here we put z position of original vertex, as it gave the best results.
    ///This has to be understood.
    TVector3 aPV(transVtx.position().x(),transVtx.position().y(),(*vertices)[0].z());
    wevent->refitPfPV(aPV);
    aPV.SetXYZ(transVtxNoBS.position().x(),transVtxNoBS.position().y(),(*vertices)[0].z());
    wevent->refitPfPVNoBS(aPV);
    wevent->isRefit(true);
  }
  else {
    TVector3 aPV((*vertices)[0].x(), (*vertices)[0].y(), (*vertices)[0].z());
    wevent->refitPfPV(aPV);
    wevent->refitPfPVNoBS(aPV);
    wevent->isRefit(false);
  }

  return true;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
TVector3 ntuple::getPCA(const edm::Event & iEvent, const edm::EventSetup & iSetup,
			const reco::Track *aTrack,	   
			const GlobalPoint & aPoint){

  TVector3 aPCA;
  if(!aTrack) return aPCA;

  edm::ESHandle<TransientTrackBuilder> transTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transTrackBuilder);  
  reco::TransientTrack transTrk=transTrackBuilder->build(aTrack);
     
  //TransverseImpactPointExtrapolator extrapolator(transTrk.field());
  AnalyticalImpactPointExtrapolator extrapolator(transTrk.field());
  GlobalPoint pos  = extrapolator.extrapolate(transTrk.impactPointState(),aPoint).globalPosition();

  aPCA.SetX(pos.x() - aPoint.x());
  aPCA.SetY(pos.y() - aPoint.y());
  aPCA.SetZ(pos.z() - aPoint.z());

  return aPCA;
}
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

//define this as a plug-in
DEFINE_FWK_MODULE(ntuple);


