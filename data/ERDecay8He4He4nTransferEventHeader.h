// -----       ERDecay8He4He4nTransferEventHeader header file        -----
// -----                  Created 08/24  by V. Chudoba               -----
// -------------------------------------------------------------------------
#ifndef ERDecay8He4He4nTransferEventHeader_H
#define ERDecay8He4He4nTransferEventHeader_H

// #include "TLorentzVector.h"
// #include "TArrayI.h"

#include "ERDecayMCEventHeader.h"

class ERDecay8He4He4nTransferEventHeader : public ERDecayMCEventHeader
{
private:
  TLorentzVector fHe8_beam;
  TLorentzVector fHe4_target;
  TLorentzVector fHe8;
  TLorentzVector fHe4;
  
  Int_t fTrigger = 0;
  Int_t fTriggerPriority = 0;
  Float_t fTime = -1.;
  Float_t fThetaCM = -1.;

public:
  ERDecay8He4He4nTransferEventHeader() : fTrigger(0), fTriggerPriority(0),
                                         fTime(-1.), fThetaCM(-1.) {}
  void SetData(const TVector3 &position, const TLorentzVector &beam,
               const TLorentzVector &target,
               const TLorentzVector &He8, const TLorentzVector &He4,
               const Float_t time, const Float_t thetaCM);

  void SetTrigger(Int_t trigger) { fTrigger = trigger; }

  TLorentzVector GetBeam() const { return fHe8_beam; }
  TLorentzVector GetTarget() const { return fHe4_target; }
  TLorentzVector GetHe8() const { return fHe8; }
  TLorentzVector GetHe4() const { return fHe4; }

  Int_t GetTrigger() const { return fTrigger; }
  Int_t GetTriggerPriority() const { return fTriggerPriority; }
  Float_t GetTime() const { return fTime; }
  Float_t GetThetaCM() const { return fThetaCM; }

  void Clear();

  ClassDef(ERDecay8He4He4nTransferEventHeader, 1)
};

#endif
