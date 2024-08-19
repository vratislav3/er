// -----                      EREXP1811EventHeader source file              -----
// -----                  Created 03/16  by V. Schetinin               -----
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
  //   TLorentzVector fp2;
  //   TLorentzVector fp3;
  //   TLorentzVector fp4;

  Int_t fTrigger = 0;
  Int_t fTriggerPriority = 0;
  float fTime = -1.;

public:
  ERDecay8He4He4nTransferEventHeader() : fTrigger(0), fTriggerPriority(0) {}
  void SetData(const TVector3 &position, const TLorentzVector &beam,
               const TLorentzVector &target,
               const TLorentzVector &He8, const TLorentzVector &He4,
               //    const TLorentzVector& p2, const TLorentzVector& p3,
               //    const TLorentzVector& p4,
               float time);

  void SetTrigger(Int_t trigger) { fTrigger = trigger; }

  Int_t GetTrigger() const { return fTrigger; }
  Int_t GetTriggerPriority() const { return fTriggerPriority; }
  TLorentzVector GetBeam() const { return fHe8_beam; }
  TLorentzVector GetTarget() const { return fHe4_target; }
  TLorentzVector GetHe8() const { return fHe8; }
  TLorentzVector GetHe4() const { return fHe4; }
  //   TLorentzVector Getp2() const { return fp2; }
  //   TLorentzVector Getp3() const { return fp3; }
  //   TLorentzVector Getp4() const { return fp4; }

  void Clear();

  ClassDef(ERDecay8He4He4nTransferEventHeader, 1)
};

#endif
