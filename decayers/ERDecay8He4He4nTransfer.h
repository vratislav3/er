// -----       			  ERDecay8He4He4nTransfer header file        -----
// -----                  Created 08/24  by V. Chudoba               -----
// -------------------------------------------------------------------------
#ifndef ERDecay8He4He4nTransfer_H
#define ERDecay8He4He4nTransfer_H

// reaction 4He(8He,4He)8He

#include "ERDecay.h"

#include "FairIon.h"

#include <fstream>
#include "TGraph.h"
#include "TF1.h"
#include "TVectorD.h"

class ERDecay8He4He4nTransfer : public ERDecay
{

public:
	ERDecay8He4He4nTransfer();
	~ERDecay8He4He4nTransfer();

	/*Modifiers*/
	void SetMinStep(Double_t minStep) { fMinStep = minStep; }
	void SetTargetThickness(Double_t targetThickness) { fTargetThickness = targetThickness; }
	// void Set4nMass(Double_t mass)
	// {
	// 	f4nMass = mass;
	// 	fIs4nUserMassSet = true;
	// }
	void Set8HeExcitation(Double_t excMean, Double_t fwhm, Double_t distibWeight);
	// void SetDecayFile(const TString &filePath, Double_t excitationEnergyInFile /*[GeV]*/) { fDecayFilePath = filePath; }

	/** @brief Sets distribution is contained in file.
	 ** @param ADfile  file with angular distribution.
	 **/
	void SetAngularDistribution(TString ADfile);

	void Print8HeExcitation();

public:
	Bool_t Init();
	Bool_t Stepping();

	void BeginEvent();
	void FinishEvent();

private:
	/** @brief Body reaction in phase space approach.
	 ** @param Ecm     Total energy in CM.
	 ** @oaram h7Mass  H7 ion mass.
	 **/
	void ReactionPhaseGenerator(Double_t Ecm, Double_t massOf8HeProduct);

	Bool_t DecayPhaseGenerator(Double_t excitation);

	std::vector<TLorentzVector> ReadDecayEvent();

private:
	TRandom3 *fRnd;
	// TRandom3 *fRnd2;

	// TParticlePDG *f8He;
	// TParticlePDG *f4He;
	// TParticlePDG *f4n;
	
	TLorentzVector *fLv8He; //!		product of 4n transfer
	TLorentzVector *fLv4He;	//!		recoil
	Float_t fTheta;			//!		theta of reaction (in CM or LAB?)
	
	TLorentzVector *fLv6He;	//!		6He from 8He decay
	TLorentzVector *fLvn1;	//!		neutron from 8He decay
	TLorentzVector *fLvn2;	//!		neutron from 8He decay

	// FairIon *fIon8He;
	// FairIon *fIon4He;

	// TGenPhaseSpace *fReactionPhaseSpace;	
	
	//TODO:  check its functionality, it should be changed for our analytical formulae
	TGenPhaseSpace *fDecayPhaseSpace;	
	Double_t fTargetReactZ;
	Double_t fMinStep;
	Double_t fTargetThickness;
	Bool_t fDecayFinish;

	std::vector<Double_t> f8HeExcitationStateMeanEnergy;
	std::vector<Double_t> f8HeExcitationStateSigma;
	std::vector<Double_t> f8HeExcitationStateWeight;
	
	Bool_t fIs8HeExcitationSet;

	// Double_t f4nMass;
	// Bool_t fIs4nUserMassSet;

	// TString fDecayFilePath;
	// Double_t fDecayFileExcitation = 1. /*[GeV]*/;
	// Bool_t fDecayFileFinished;
	// Int_t fDecayFileCurrentEvent;
	// std::ifstream fDecayFile;

	TGraph *fADInput = nullptr; //!   distribution (angular distribution) graph containing AD input
	TF1 *fADFunction = nullptr; //!   function describing AD (angular distribution) of binary reaction
	Double_t fThetaMin = 0.;
	Double_t fThetaMax = 0.;

	// ADEvaluate function is necessary for TF1 constructor
	Double_t ADEvaluate(Double_t *x, Double_t *p);

	ClassDef(ERDecay8He4He4nTransfer, 1)
};

#endif
