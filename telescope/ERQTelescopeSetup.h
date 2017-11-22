// -------------------------------------------------------------------------

// -----                        ERQTelescopeSetup header file              -----

// -----                              -----

// -------------------------------------------------------------------------
#ifndef ERQTelescopeSetup_H
#define ERQTelescopeSetup_H

#include <map>
#include <vector>

#include "TString.h"
#include <TXMLNode.h>
#include "Rtypes.h"

using namespace std;

class ERQTelescopeSetup {
public:
  ERQTelescopeSetup();
  ~ERQTelescopeSetup();
  static ERQTelescopeSetup* Instance();

  /* Accessors */
  static void PrintDetectorParameters(void);
  static void PrintDetectorParametersToFile(TString fileName);

  /* Modifiers */
  static void SetXmlParametersFile(TString xmlFileName) {fParamsXmlFileName = xmlFileName;}
  static void AddSi(TString type, Double_t position, TString orientAroundZ); 
  static void AddCsI(TString type, Double_t position);

  static void ConstructGeometry();
  static Int_t SetParContainers();

private:
  static void ParseXmlParameters();
  static void GetCsIParameters(TXMLNode *node);
  static void GetSingleSiParameters(TXMLNode *node);
  static void GetDoubleSiParameters(TXMLNode *node);

  // ----- Sensetive volumes parameters parameters ----------------------------
  static ERQTelescopeSetup* fInstance;
  static TString          fParamsXmlFileName;
  static vector<TString>  fSensVolumes;

  // ----- SingleSi parameters --------------------------------------------------
  static Int_t            fDoubleSiCount;
  static vector<TString>  fDoubleSiType;
  static vector<Double_t> fDoubleSiPosZ;
  static vector<TString>  fDoubleSiOrientAroundZ;
  static vector<Double_t> fDoubleSiX;
  static vector<Double_t> fDoubleSiY;
  static vector<Double_t> fDoubleSiZ;
  static vector<Double_t> fDoubleSiSensX;
  static vector<Double_t> fDoubleSiSensY;
  static vector<Double_t> fDoubleSiSensZ;
  static vector<Int_t>    fDoubleSiStripCountX;
  static vector<Int_t>    fDoubleSiStripCountY;
  static vector<TString>  fDoubleSiMedia;
  // ----- DoubleSi parameters --------------------------------------------------
  static Int_t            fSingleSiCount;
  static vector<Double_t> fSingleSiPosZ;
  static vector<TString>  fSingleSiType;  
  static vector<TString>  fSingleSiOrientAroundZ;
  static vector<Double_t> fSingleSiX;
  static vector<Double_t> fSingleSiY;
  static vector<Double_t> fSingleSiZ;
  static vector<Double_t> fSingleSiSensX;
  static vector<Double_t> fSingleSiSensY;
  static vector<Double_t> fSingleSiSensZ;
  static vector<Int_t>    fSingleSiStripCount;
  static vector<TString>  fSingleSiMedia;
  // ----- CsI parameters -------------------------------------------------------
  static Int_t            fCsICount;
  static vector<TString>  fCsIType;
  static vector<Double_t> fCsIPosZ;  
  static vector<TString>  fCsIOrientAroundZ;
  static vector<Double_t> fCsIX;
  static vector<Double_t> fCsIY;
  static vector<Double_t> fCsIZ;
  static vector<Int_t>    fCsICubesCountX;
  static vector<Int_t>    fCsICubesCountY;
  static vector<TString>  fCsIMedia;

  ClassDef(ERQTelescopeSetup,1)
};
#endif

