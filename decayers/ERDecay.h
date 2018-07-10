/********************************************************************************
 *              Copyright (C) Joint Institute for Nuclear Research              *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *         GNU Lesser General Public Licence version 3 (LGPL) version 3,        *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/
#ifndef ERDECAYER_H
#define ERDECAYER_H

#include "TRandom3.h"
#include "TString.h"
#include "TGenPhaseSpace.h"
#include "TDatabasePDG.h"               //for TDatabasePDG
#include "TVector3.h"

class ERDecay{
public:
  ERDecay(TString name);
  ~ERDecay();
  /*Accessors*/
  TString GetName() const {return fName;}
  
  /*Modifiers*/
  
  /** @brief Defines name of volume where occures interaction of ion.
   ** @param name  volume name
  **/   
  void SetInteractionVolumeName(TString name) {fInteractionVolumeName = name;}
  
  /** @brief Defines nuclear interaction length.
   ** @param lambda  nuclear interaction length
  **/ 
  void SetNuclearInteractionLength(Double_t lambda) {fNuclearInteractionLength = lambda;}
  
  /** @brief Defines maximum path length for particles in the volume that is defined 
             in SetInteractionVolumeName().
             In macro, method must be called after SetNuclearInteractionLength().
   ** @param pathLength  maximum path length in target volume
  **/  
  void SetMaxPathLength(Double_t pathLength);

  /** @brief Method returns interaction volume name.
  **/
  TString GetInteractionVolumeName() {return fInteractionVolumeName;}

  virtual Bool_t Stepping() = 0;
  virtual void   BeginEvent() = 0;
  virtual void   FinishEvent() = 0;
  virtual Bool_t Init() = 0;

protected:
  /** @brief  
   ** @param Method must be called in Stepping() when particle makes first step in
             interaction volume;
  **/  
  Bool_t FindInteractionPoint();
  // void   CalculateTargetParameters();
  
	TString  fName;
  Bool_t   fIsInterationPointFound;
  Double_t fDistanceToInteractPoint;
  Double_t fDistanceFromEntrance;

private:
  TString  fInteractionVolumeName;
  Double_t fNuclearInteractionLength;
  Double_t fInteractionProbability;
  // Double_t fTargetBoundBoxDiagonal;
  Double_t fNormalizingProbability;
  TRandom3 *fRnd1;
  TRandom3 *fRnd2;

	ClassDef(ERDecay,1)
};

#endif
