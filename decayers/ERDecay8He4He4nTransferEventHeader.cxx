/********************************************************************************
 *              Copyright (C) Joint Institute for Nuclear Research              *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *         GNU Lesser General Public Licence version 3 (LGPL) version 3,        *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/

#include "ERDecay8He4He4nTransferEventHeader.h"

#include "FairLogger.h"

void ERDecay8He4He4nTransferEventHeader::SetData(const TVector3 &position, const TLorentzVector &beam,
                                                 const TLorentzVector &target,
                                                 const TLorentzVector &He8, const TLorentzVector &He4,
                                                 // const TLorentzVector& p2,
                                                 // const TLorentzVector& p3, const TLorentzVector& p4,
                                                 const float time)
{
  fReactionPos = position;
  fHe8_beam = beam;
  fHe4_target = target;
  fHe8 = He8;
  fHe4 = He4;
  //   fp1 = p1;
  //   fp2 = p2;
  //   fp3 = p3;
  //   fp4 = p4;
  fTime = time;
}
// -------------------------------------------------------------------------
void ERDecay8He4He4nTransferEventHeader::Clear()
{
  ERDecayMCEventHeader::Clear();
  fHe8_beam.SetXYZM(0, 0, 0, 0);
  fHe4_target.SetXYZM(0, 0, 0, 0);
  fHe8.SetXYZM(0, 0, 0, 0);
  fHe4.SetXYZM(0, 0, 0, 0);
  // fp2.SetXYZM(0, 0, 0, 0);
  // fp3.SetXYZM(0, 0, 0, 0);
  // fp4.SetXYZM(0, 0, 0, 0);
  fTrigger = 0;
  fTriggerPriority = 0;
  fTime = -1.;
}
// -------------------------------------------------------------------------

ClassImp(ERDecay8He4He4nTransferEventHeader)