// -------------------------------------------------------------------------
// -----                   ERIonMixGenerator source file                 -----
// -----                   Created by            -----
// -------------------------------------------------------------------------
#include "ERIonMixGenerator.h"

#include "FairIon.h"                    // for FairIon
#include "FairParticle.h"               // for FairParticle
#include "FairPrimaryGenerator.h"       // for FairPrimaryGenerator
#include "FairRunSim.h"                 // for FairRunSim
#include "FairLogger.h"                 // for logging

#include "Riosfwd.h"                    // for ostream
#include "TDatabasePDG.h"               // for TDatabasePDG
#include "TObjArray.h"                  // for TObjArray
#include "TParticle.h"                  // for TParticle
#include "TParticlePDG.h"               // for TParticlePDG
#include "TRandom.h"

#include <stdio.h>                      // for NULL, sprintf

// -----   Default constructor   ------------------------------------------
ERIonMixGenerator::ERIonMixGenerator()
  : ERIonGenerator()
{
//  LOG(WARNING) << "ERIonMixGenerator: "
//               << " Please do not use the default constructor! " 
//               << FairLogger::endl;
}
// ------------------------------------------------------------------------

// -----   Default constructor   ------------------------------------------
ERIonMixGenerator::ERIonMixGenerator(TString name, Int_t z, Int_t a, Int_t q, Int_t mult)
  : ERIonGenerator(name, z, a, q, mult)
{
  fBgIons.clear();
  fBgIons.insert(std::make_pair(name, 0));
}
//_________________________________________________________________________

// -----   Destructor   ---------------------------------------------------
ERIonMixGenerator::~ERIonMixGenerator()
{
// if (fIon) delete fIon;
}

void ERIonMixGenerator::AddBackgroundIon(TString name, Int_t z, Int_t a, Int_t q, Double_t probability)
{
  SetPhiRange();
  static Double_t sumProbability;

  if(fBgIons.size() == 0)
  {
    sumProbability = 0;
  }

  sumProbability += probability;

  if((sumProbability) >= 1)
  {
    LOG(DEBUG) << "Summary probability of appearing background ions more then 1"
       << FairLogger::endl;
    return ;
  }

  fBgIons.insert(std::make_pair(name, sumProbability));

  //fBgIons.at(fIon->GetName()) = 1;

  FairRunSim* run = FairRunSim::Instance();
  if ( ! run ) {
    LOG(ERROR) << "No FairRun instantised!" 
         << FairLogger::endl;
  } else {           
    run->AddNewIon(new FairIon(name, z, a, q));
  }
}

Bool_t ERIonMixGenerator::ReadEvent(FairPrimaryGenerator* primGen)
{
  Double_t randResult;
  TString  ionName;
  // Generate particles
  for (Int_t k = 0; k < fMult; k++) {
    spreadingParameters();

    randResult = gRandom->Uniform(0, 1);

    auto it = std::find_if(fBgIons.begin(), fBgIons.end(),
                            [randResult](const std::pair<TString, Double_t> &t)->bool
                            {
                              return (t.second > randResult);
                            }
                          ); 
    (it == fBgIons.end()) ? ionName = fIon->GetName() : ionName = it->first;


    TParticlePDG* thisPart =
    TDatabasePDG::Instance()->GetParticle(ionName);
    if ( ! thisPart ) {
      LOG(WARNING) << "ERIonGenerator: Ion " << ionName
      << " not found in database!" << FairLogger::endl;
      return kFALSE;
    }

    int pdgType = thisPart->PdgCode();

    std::cout << "ERIonGenerator: Generating " << fMult << " ions of type "
        << ionName << " (PDG code " << pdgType << ")" 
        << std::endl;
    std::cout << "    Momentum (" << fPx << ", " << fPy << ", " << fPz
        << ") Gev from vertex (" << fX << ", " << fY
        << ", " << fZ << ") cm" << std::endl;
    primGen->AddTrack(pdgType, fPx, fPy, fPz, fX, fY, fZ);
  }

  return kTRUE;

}

ClassImp(ERIonMixGenerator)