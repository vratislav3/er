// -----       			  script for vacuum target used to study of kinematics        -----
// -----                        Created 08/24  by V. Chudoba                      -----
// ------------------------------------------------------------------------------------

#if !defined(__CLING__)

// standard ROOT includes
#include "TString.h"
#include "TSystem.h"
#include "TGeoManager.h"
#include "TROOT.h"
#include "TGeoTube.h"
#include "TFile.h"

// FairRoot includes
#include "FairGeoLoader.h"
#include "FairGeoInterface.h"
#include "FairGeoMedia.h"
#include "FairGeoBuilder.h"
#include "FairGeoMedium.h"

#endif

void create_target_vacuum()
{
  TString erPath = gSystem->Getenv("VMCWORKDIR");

  // Output paths
  TString outGeoFilenameRoot = erPath + "/geometry/target.Vacuum.geo.root";

  // Input paths
  TString medFile = erPath + "/geometry/media.geo";

  // Materials and media
  FairGeoLoader *geoLoad = new FairGeoLoader("TGeo", "FairGeoLoader");
  FairGeoInterface *geoFace = geoLoad->getGeoInterface();
  geoFace->setMediaFile(medFile);
  geoFace->readMedia();
  FairGeoMedia *geoMedia = geoFace->getMedia();
  FairGeoBuilder *geoBuild = geoLoad->getGeoBuilder();

  // Geometry manager
  TGeoManager *geoM = (TGeoManager *)gROOT->FindObject("FAIRGeom");

  TString mediumName;

  mediumName = "vacuum";
  FairGeoMedium *mVacuum = geoMedia->getMedium(mediumName);
  if (!mVacuum)
    Fatal("create_target_vacuum", "FairMedium %s not found", mediumName.Data());
  geoBuild->createMedium(mVacuum);
  TGeoMedium *pVacuum = geoM->GetMedium(mediumName);
  if (!pVacuum)
    Fatal("create_target_vacuum", "Medium %s not found", mediumName.Data());

  // General dimensions
  Double_t transX = 0.; // cm
  Double_t transY = 0.; // cm
  // Double_t transZ = 1.; // cm
  Double_t transZ = 0.;        // cm
  Double_t targetRadius = 1.5; // cm
  Double_t targetZ = 0.6;      // cm

  // Shapes
  TGeoTube *targetShape = new TGeoTube("targetShape",
                                       0.,
                                       targetRadius,
                                       targetZ / 2.);

  // Volumes
  TGeoVolume *targetVol = new TGeoVolume("targetSensVol", targetShape, pVacuum);

  // Matrices
  TGeoRotation *rotNoRot = new TGeoRotation("rotNoRot", 0., 0., 0.);
  rotNoRot->RegisterYourself();

  // This is the one but last level in the hierarchy
  // This volume-assembly is the only volume to be inserted into TOP
  TGeoVolumeAssembly *subdetectorVolAss = new TGeoVolumeAssembly("target_Vacuum");
  subdetectorVolAss->AddNode(targetVol, 1,
                             new TGeoCombiTrans("mTargetVolInTarget", transX, transY, transZ, rotNoRot));

  // World ------------------------------------
  TGeoVolumeAssembly *topVolAss = new TGeoVolumeAssembly("TOP");
  topVolAss->AddNode(subdetectorVolAss, 1);

  // Finalize
  geoM->SetTopVolume(topVolAss);
  geoM->CloseGeometry();
  geoM->CheckOverlaps();
  geoM->PrintOverlaps();
  // geoM->CheckGeometry();
  // geoM->CheckGeometryFull();
  // geoM->Test();

  // Export
  TFile *outGeoFileRoot = new TFile(outGeoFilenameRoot, "RECREATE");
  geoM->GetTopVolume()->Write();
  outGeoFileRoot->Close();

  // Draw
  geoM->GetTopVolume()->Draw("");
}
