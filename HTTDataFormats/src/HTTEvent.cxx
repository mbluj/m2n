#ifdef PROJECT_NAME
#include "m2n/HTTDataFormats/interface/HTTEvent.h"
#include "m2n/maxi2ntuples/interface/GenInfoHelper.h"
#else
#include "HTTEvent.h"
#endif


Wtau::Wtau(){
  clear();
}

Wtau::~Wtau(){;}

void Wtau::clear(){

  pt_ = -999;

  leadingTk_ = TLorentzVector();

}

Wmu::Wmu(){
  clear();
}
Wmu::~Wmu(){;}

void Wmu::clear(){

  pt_ = -999;
}

Welectron::Welectron(){
  clear();
}
Welectron::~Welectron(){;}

void Welectron::clear(){

  pt_ = -999;
}

Wevent::Wevent(){run_ = 0; lumi_ = 0; event_ = 0;}
Wevent::~Wevent(){;}

void Wevent::clear(){

#ifdef PROJECT_NAME
  decModeMinus_ = WawGenInfoHelper::tauDecayModes::tauDecayOther;
  decModePlus_ = WawGenInfoHelper::tauDecayModes::tauDecayOther;
#else
  decModeMinus_  = 99;
  decModePlus_  = 99;
#endif

  thePV_ = TVector3();
  pfPV_ = TVector3();
  refitPfPV_ = TVector3();
  refitPfPVNoBS_ = TVector3();

  isRefit_ = false;
  nTracksInRefit_ = 0;

}

Wpair::Wpair(){;}
Wpair::~Wpair(){;}
void Wpair::clear(){

  svMinus_ = TVector3();
  svPlus_ = TVector3();
  nPiPlus_ = TVector3();
  nPiMinus_ = TVector3();

}
