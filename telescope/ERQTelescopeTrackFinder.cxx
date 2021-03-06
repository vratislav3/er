/********************************************************************************
 *              Copyright (C) Joint Institute for Nuclear Research              *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *         GNU Lesser General Public Licence version 3 (LGPL) version 3,        *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/

#include "ERQTelescopeTrackFinder.h"

#include "TVector3.h"
#include "TMath.h"

#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "FairLogger.h"

#include "ERBeamDetTrack.h"

using namespace std;

//--------------------------------------------------------------------------------------------------
ERQTelescopeTrackFinder::ERQTelescopeTrackFinder()
  : FairTask("ER qtelescope track finding scheme"),
    fUserTargetPointIsSet(kFALSE),
    fBeamDetTrack(NULL)
{
}
//--------------------------------------------------------------------------------------------------
ERQTelescopeTrackFinder::ERQTelescopeTrackFinder(Int_t verbose)
  : FairTask("ER qtelescope track finding scheme ", verbose),
    fUserTargetPointIsSet(kFALSE),
    fBeamDetTrack(NULL)
{
}
//--------------------------------------------------------------------------------------------------
ERQTelescopeTrackFinder::~ERQTelescopeTrackFinder() {
}
//--------------------------------------------------------------------------------------------------
void ERQTelescopeTrackFinder::SetHitStation(TString stationID) {
  TString xStripArrayName = stationID;
  TString yStripArrayName = stationID;
  fSiHitStationsPair.emplace(make_pair(stationID, 
                                       pair<TString, TString>(xStripArrayName.Append("_X"), 
                                                              yStripArrayName.Append("_Y"))));
}
//--------------------------------------------------------------------------------------------------
void ERQTelescopeTrackFinder::SetHitStation(TString xStationID, TString yStationID) {
  fSiHitStationsPair.emplace(make_pair( xStationID + yStationID,
                                       pair<TString, TString>(xStationID, yStationID)));
}
void ERQTelescopeTrackFinder::SetStripEdepRange(Double_t edepMin, Double_t edepMax) {
  fSiDigiEdepMin = edepMin; 
  fSiDigiEdepMax = edepMax;
}
//--------------------------------------------------------------------------------------------------
void ERQTelescopeTrackFinder::SetTargetPoint(Double_t x, Double_t y, Double_t z) {
  fUserTargetPointIsSet = kTRUE;
  fTargetX = x;
  fTargetY = y;
  fTargetZ = z;
}
//--------------------------------------------------------------------------------------------------
InitStatus ERQTelescopeTrackFinder::Init() {
  FairRootManager* ioman = FairRootManager::Instance();
  if ( ! ioman ) Fatal("Init", "No FairRootManager");

  TList* allbrNames = ioman->GetBranchNameList();
  TIter nextBranch(allbrNames);
  TObjString* bName;
  vector<TString> pointBranches;
  while (bName = (TObjString*)nextBranch()) {
    TString bFullName = bName->GetString();
    if (bFullName.Contains("Digi") && bFullName.Contains("QTelescope")) {
      Int_t bPrefixNameLength = bFullName.First('_'); 
      TString brName(bFullName(bPrefixNameLength + 1, bFullName.Length()));
      fQTelescopeDigi[brName] = (TClonesArray*) ioman->GetObject(bFullName);
    }
  }
  for (const auto itHitStations : fSiHitStationsPair) {
    fQTelescopeTrack[itHitStations.first] = new TClonesArray("ERQTelescopeTrack");
    ioman->Register("ERQTelescopeTrack_" + itHitStations.first, "QTelescope", 
                    fQTelescopeTrack[itHitStations.first], kTRUE);
  }

  fBeamDetTrack = (TClonesArray*) ioman->GetObject("BeamDetTrack");   
  if (!fUserTargetPointIsSet) {
    if (!fBeamDetTrack) {
      LOG(DEBUG) << "ERQTelescopeTrackFinder: target point not initialized by user " 
                 << "(by means of SetTargetPoint()) and there is no ERBeamDetTrack branch" 
                 <<FairLogger::endl;
      return kFATAL;
    } 
  }

  fQTelescopeSetup = ERQTelescopeSetup::Instance();
  fQTelescopeSetup->ReadGeoParamsFromParContainer();
  return kSUCCESS;
}
//--------------------------------------------------------------------------------------------------
void ERQTelescopeTrackFinder::Exec(Option_t* opt) { 
  Reset();
  for (const auto &itSiStationPair : fSiHitStationsPair) {
    vector<pair<Int_t, Int_t>>  hitTelescopePoint;
    vector<Int_t>               correctStripsX;
    vector<Int_t>               correctStripsY;
    TString xDigiBranchName = itSiStationPair.second.first;
    TString yDigiBranchName = itSiStationPair.second.second;  
    TClonesArray *xDigi = fQTelescopeDigi[xDigiBranchName];
    TClonesArray *yDigi = fQTelescopeDigi[yDigiBranchName];
    if ( !xDigi || !yDigi) {
      continue;
    } 
    for (Int_t iXDigi  = 0; iXDigi < xDigi->GetEntriesFast(); iXDigi++) {
      Double_t xStripEdep = ((ERQTelescopeSiDigi*)xDigi->At(iXDigi))->GetEdep();
      if (xStripEdep > fSiDigiEdepMin && xStripEdep < fSiDigiEdepMax) {
        correctStripsX.push_back(iXDigi);
      }
    }
    for (Int_t iYDigi  = 0; iYDigi < yDigi->GetEntriesFast(); iYDigi++) {
      Double_t yStripEdep = ((ERQTelescopeSiDigi*)yDigi->At(iYDigi))->GetEdep();
      if (yStripEdep > fSiDigiEdepMin && yStripEdep < fSiDigiEdepMax) {
        correctStripsY.push_back(iYDigi);
      }
    }
    for (const auto itCorrectStripsX : correctStripsX) {
      Double_t xStripEdep = ((ERQTelescopeSiDigi*)xDigi->At(itCorrectStripsX))->GetEdep();
      for (const auto itCorrectStripsY : correctStripsY) {
        Double_t yStripEdep = ((ERQTelescopeSiDigi*)yDigi->At(itCorrectStripsY))->GetEdep();
        if (TMath::Abs(xStripEdep - yStripEdep) < fEdepDiffXY) {
          hitTelescopePoint.push_back(pair<Int_t, Int_t>(itCorrectStripsX, itCorrectStripsY));
        }
      }
    }
    LOG(DEBUG) << "Strips array pair " << itSiStationPair.second.first << " " 
                                       << itSiStationPair.second.second << FairLogger::endl; 
    for (auto &itHitPoint : hitTelescopePoint) {
      Int_t xStripNb = ((ERQTelescopeSiDigi*)xDigi->At(itHitPoint.first))->GetStripNb();
      Int_t yStripNb = ((ERQTelescopeSiDigi*)yDigi->At(itHitPoint.second))->GetStripNb();
      LOG(DEBUG) << "  Strips pair " << itHitPoint.first << " " << itHitPoint.second << FairLogger::endl;
      Double_t xQtelescopeHit = fQTelescopeSetup->GetStripX(xDigiBranchName, xStripNb);
      Double_t yQtelescopeHit = fQTelescopeSetup->GetStripY(yDigiBranchName, yStripNb);
      Double_t zQtelescopeHit = (fQTelescopeSetup->GetStripZ(xDigiBranchName, xStripNb) 
                              +  fQTelescopeSetup->GetStripZ(yDigiBranchName, yStripNb)) / 2;
      if (!fUserTargetPointIsSet) {
        ERBeamDetTrack* trackFromMWPC = (ERBeamDetTrack*)fBeamDetTrack->At(0);
        if (!trackFromMWPC) {
          FairRun* run = FairRun::Instance();
          run->MarkFill(kFALSE);
          return ;
        }
        fTargetX = trackFromMWPC->GetTargetX();
        fTargetY = trackFromMWPC->GetTargetY();
        fTargetZ = trackFromMWPC->GetTargetZ();
      }
      if (!fBeamDetTrack) {
        FairRun* run = FairRun::Instance();
          run->MarkFill(kFALSE);
          return ;
      }
      AddTrack(fTargetX, fTargetY, fTargetZ, xQtelescopeHit, yQtelescopeHit, zQtelescopeHit,
               itSiStationPair.first);
    }
  }
}
//--------------------------------------------------------------------------------------------------
void ERQTelescopeTrackFinder::Reset() {
  for (const auto itTrackBranches : fQTelescopeTrack) {
    if (itTrackBranches.second) {
      itTrackBranches.second->Delete();
    }
  }
}
//--------------------------------------------------------------------------------------------------
void ERQTelescopeTrackFinder::Finish() {   
}
//--------------------------------------------------------------------------------------------------
ERQTelescopeTrack* ERQTelescopeTrackFinder::AddTrack(Double_t targetX, Double_t targetY,   
                                                     Double_t targetZ, 
                                                     Double_t telescopeX, Double_t telescopeY, 
                                                     Double_t telescopeZ,
                                                     TString digiBranchName) 
{
  ERQTelescopeTrack *track = new((*fQTelescopeTrack[digiBranchName])
                                                   [fQTelescopeTrack[digiBranchName]->GetEntriesFast()])
                             ERQTelescopeTrack( targetX,    targetY,    targetZ, 
                                                telescopeX, telescopeY, telescopeZ);
  return track;
}
//--------------------------------------------------------------------------------------------------
void ERQTelescopeTrackFinder::SetParContainers() {
  // Get run and runtime database
  FairRun* run = FairRun::Instance();
  if ( ! run ) Fatal("SetParContainers", "No analysis run");

  FairRuntimeDb* rtdb = run->GetRuntimeDb();
  if ( ! rtdb ) Fatal("SetParContainers", "No runtime database");
}
ClassImp(ERQTelescopeTrackFinder)