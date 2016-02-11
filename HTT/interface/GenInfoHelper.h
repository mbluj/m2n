// M. Bluj

#ifndef WarsawAnalysis_HTT_GenInfoHelper_h
#define WarsawAnalysis_HTT_GenInfoHelper_h

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/TauReco/interface/PFTauDecayMode.h"
#include <utility>

//Root
#include "TLorentzVector.h"
#include "TVector3.h"

namespace WawGenInfoHelper {
  
  typedef reco::PFTauDecayMode::hadronicTauDecayModes tauDecayModes;

  
  enum bosonDecayModes {kUndefined = -1,
			 kMuTau,kETau,kTauTau,kMuMu,kEE,kEMu,
			 kEEPrompt, kMuMuPrompt, kTauTauPrompt};
  
  typedef reco::GenParticleCollection::const_iterator IG;
  typedef reco::GenParticleRefVector::const_iterator IGR;

  bool isBoson(const reco::GenParticleRef& particle, 
	      bool checkLastCopy=false, bool tauDec=true);
  const reco::GenParticleRef getFinalClone(const reco::GenParticleRef& particle, bool isUnstable=true);
  bool isFinalClone(const reco::GenParticleRef& particle, bool isUnstable=true);
  const reco::GenParticleRef getInitialClone(const reco::GenParticleRef& initialClone);
  bool isInitialClone(const reco::GenParticleRef& particle);
  int getBosonDecayMode(const reco::GenParticleRef& particle);
  int getLeptonPairDecayMode(const reco::GenParticleCollection& sourceParticles);
  int getDiTauDecayMode(int tau1DecayMode, int tau2DecayMode);
  int getTauDecayMode(const reco::GenParticleRefVector& products);
  int getTauDirDecayMode(const reco::GenParticleRefVector& products);
  int getTausDecays(const reco::GenParticleRef& tau,
		    reco::GenParticleRefVector& products,
		    bool ignoreNus=true,bool direct=false);
  const reco::GenParticleRef getLeadChParticle(const reco::GenParticleRefVector& products);
  TLorentzVector getP4(const reco::GenParticleRef& part);
  TLorentzVector getP4(const reco::GenParticle&  part);
  void setP4Ptr(const TLorentzVector& p4, TLorentzVector *p4Ptr);
  void setV3Ptr(const TVector3& v3, TVector3 *v3Ptr);
  TLorentzVector getCombinedP4(const reco::GenParticleRefVector& products);
  TLorentzVector getLeadChParticleP4(const reco::GenParticleRefVector& products);
  TVector3 getVertex(const reco::GenParticleRef& part);
  void getVertex(const reco::GenParticleRef& part, TVector3 *vtx);
  TVector3 impactParameter(const TVector3& pv, const TVector3& sv, const TLorentzVector& p4);
  void impactParameter(const TVector3& pv, const TVector3& sv, const TLorentzVector& p4,
		       TVector3 *ip);
  TLorentzVector getGenMet(const reco::GenParticleRefVector& particles);
  TLorentzVector getGenMet(const reco::GenParticleCollection&  particles);
  TLorentzVector getTauNuMet(const reco::GenParticleRefVector& taus);
  std::pair<float,float> angleBetweenPlanes(const TLorentzVector &prod1, 
					    const TLorentzVector &prod12,
					    const TLorentzVector &prod2, 
					    const TLorentzVector &prod22);
  //Check recursively if any ancestor of particle is the given one
  bool isAncestor(const reco::Candidate *ancestor, const reco::Candidate *particle);
  //Check recursively if any descendent of particle is the given one
  bool isDescendent(const reco::Candidate *descendent, const reco::Candidate *particle);

  /// find all particles of a given pdgId and status
  /// copy from PhysicsTools/HepMCCandAlgos/interface/GenParticlesHelper.h" which allows ignore status(<=0) or pdgId(==0)
  void findParticles(const reco::GenParticleCollection& sourceParticles,
                     reco::GenParticleRefVector& particleRefs, 
                     int pdgId, int status);
  /// find all descendents of a given status and pdgId (recursive)
  /// copy from PhysicsTools/HepMCCandAlgos/interface/GenParticlesHelper.h" which allows ignore status(<0)
  void findDescendents(const reco::GenParticleRef& base, 
		       reco::GenParticleRefVector& descendents, 
		       int status, int pdgId=0,
		       bool skipPhotonsPi0AndFSR=false);

  void findAncestors(const reco::GenParticleRef& base, 
		     reco::GenParticleRefVector& ancestors, 
		     int status, int pdgId=0);
}

#endif
