#ifndef m2n_maxi2ntuples_utilities_h
#define m2n_maxi2ntuples_utilities_h


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <DataFormats/PatCandidates/interface/Tau.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <vector> 
namespace utilities{


    inline double absIso(const pat::Muon& l, double dBetaFactor = 0){
    
        if(dBetaFactor>0 && l.puChargedHadronIso()<0)
            std::cout << "If you want to use dbeta corrections, you must make sure that the pu charged hadron iso is available. This should never happen" ;
        double neutralIso = l.neutralHadronIso() + l.photonIso();
        double corNeutralIso = neutralIso - dBetaFactor * l.puChargedHadronIso();
        double charged = l.chargedHadronIso();  
        return charged + std::max(0., corNeutralIso);

    }


    inline double absIso(const pat::Electron& l, double dBetaFactor = 0){
    
        if(dBetaFactor>0 && l.puChargedHadronIso()<0)
            std::cout << "If you want to use dbeta corrections, you must make sure that the pu charged hadron iso is available. This should never happen" ;
        double neutralIso = l.neutralHadronIso() + l.photonIso();
        double corNeutralIso = neutralIso - dBetaFactor * l.puChargedHadronIso();
        double charged = l.chargedHadronIso();  
        return charged + std::max(0., corNeutralIso);

    }
    

    inline double relIso(const pat::Muon& l, double dBetaFactor = 0){
        return absIso(l, dBetaFactor)/l.pt();
    }


    inline double relIso(const pat::Electron& l, double dBetaFactor = 0){
        return absIso(l, dBetaFactor)/l.pt();
    }



    inline bool heppymuonID(const pat::Muon& l, std::string name){
 
        if(name == "POG_ID_TightNoVtx")
            return l.isLooseMuon() && 
                l.isGlobalMuon() && 
                l.globalTrack()->normalizedChi2() < 10 && 
                l.globalTrack()->hitPattern().numberOfValidMuonHits() > 0 && 
                l.numberOfMatchedStations()>1 && 
                l.innerTrack()->hitPattern().numberOfValidPixelHits()>0 && 
                l.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5;
        if(name == "POG_ID_Medium"){
            bool goodGlb = l.isGlobalMuon() && l.globalTrack()->normalizedChi2() < 3 && l.combinedQuality().chi2LocalPosition < 12 && l.combinedQuality().trkKink < 20;
            return l.innerTrack()->validFraction() >= 0.8 && l.segmentCompatibility() >= (goodGlb ? 0.303 : 0.451);
        }
        return false;
    }

    inline bool isbJet(const pat::Jet& j){
    
        return j.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags") > 0.814;
    }

    inline bool jetID(const pat::Jet& j){

        float energy = 1;
        if(true)
            energy = (j.p4()*j.jecFactor("Uncorrected")).energy();
        float NHF = j.neutralHadronEnergyFraction()/energy;
        float NEMF = j.neutralEmEnergyFraction()/energy;
        float CHF = j.chargedHadronEnergyFraction()/energy;
        float MUF = j.muonEnergyFraction()/energy;
        float CEMF = j.chargedEmEnergyFraction()/energy;
        float NumConst = j.chargedMultiplicity()+j.neutralMultiplicity();
        float CHM = j.chargedMultiplicity(); 
        return (NHF<0.99 && NEMF<0.99 && NumConst>1 && MUF<0.8) && ((abs(j.eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(j.eta())>2.4); //looseID
        return (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((abs(j.eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || abs(j.eta())>2.4); //tightID
    }

    inline bool pujetetaid(const pat::Jet& j){
        
       if(abs(j.eta()) < 2.5)
          return j.userFloat("pileupJetId:fullDiscriminant") > -0.63;
       else if(abs(j.eta()) < 2.75)
          return j.userFloat("pileupJetId:fullDiscriminant") > -0.60;
       else if(abs(j.eta()) < 3.0)
          return j.userFloat("pileupJetId:fullDiscriminant") > -0.55;
       else if(abs(j.eta()) < 5.2)
          return j.userFloat("pileupJetId:fullDiscriminant") > -0.45;
       else{
           std::cout << "ERROR: utilities:pujetID";
           return -100;
        }
    }

};



#endif 
