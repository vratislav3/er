
/********************************************************************************
 *              Copyright (C) Joint Institute for Nuclear Research              *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *         GNU Lesser General Public Licence version 3 (LGPL) version 3,        *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/

#include "ERDecay8He4He4nTransfer.h"

// #include <iostream>
// #include <string>
// #include <sstream>
// using namespace std;

#include "TVirtualMC.h"
#include "Math/Integrator.h"
// #include "TLorentzVector.h"
// #include "TMCProcess.h"
// #include "TRandom.h"
// #include "TVector.h"

#include "FairRunSim.h"
#include "FairLogger.h"

#include "ERDecayMCEventHeader.h"
#include "ERDecay8He4He4nTransferEventHeader.h"
#include "ERMCEventHeader.h"

#include "G4IonTable.hh"
#include "G4Neutron.hh"

ERDecay8He4He4nTransfer::ERDecay8He4He4nTransfer() : ERDecay("8He4He4nTransfer"),
													 fDecayFinish(kFALSE),
													 fMinStep(0.01),
													 fADInput(NULL),
													 fADFunction(NULL)
//   fDecayFilePath(""),
//   fDecayFileFinished(kFALSE),
//   fDecayFileCurrentEvent(0)
{
	fRnd = new TRandom3();
	// fRnd->SetSeed();
	//   fRnd2 = new TRandom3();
	//   fRnd2->SetSeed();
	// fReactionPhaseSpace = new TGenPhaseSpace();
	fDecayPhaseSpace = new TGenPhaseSpace();

	fLv8He = new TLorentzVector();
	fLv4He = new TLorentzVector();

	fLv6He = new TLorentzVector();
	fLvn1 = new TLorentzVector();
	fLvn2 = new TLorentzVector();

	fLv6HeCMdecay = new TLorentzVector();
	fLvn1CMdecay = new TLorentzVector();
	fLvn2CMdecay = new TLorentzVector();

	fLvNN = new TLorentzVector();

	// It is necessary to create FairIon objects for all
	// isotopes which will be used in (this?) class. It seems
	// that only if they are added to the current FairRunSim
	// instance, the are added to the TDatabasePDG iontable.
	// This procedure need to be done in constructor of this
	// class. If done in this->Init(), segmentation violation
	// will appear.
	//
	// The masses of ions added to TDatabasePDG were checked
	// and correspond to Geant4 masses.
	//
	// The beam ion added in sim.C using the ERIonMixGenerator
	// is added to the TDatabasePDG automatically in other
	// place and it may be not treated here.
	//
	// TODO: solve this possible unexpected behaviour

	FairRunSim *run = FairRunSim::Instance();

	FairIon *fIon6He = new FairIon("6He", 2, 6, 2);
	FairIon *fIon8He = new FairIon("8He", 2, 8, 2);

	run->AddNewIon(fIon6He);
	run->AddNewIon(fIon8He);

	LOG(INFO) << "[ERDecay8He4He4nTransfer::ERDecay8He4He4nTransfer()] object " << GetName() << " created" << FairLogger::endl;
}

//-------------------------------------------------------------------------------------------------
ERDecay8He4He4nTransfer::~ERDecay8He4He4nTransfer()
{
	//   if (fDecayFile.is_open())
	//     fDecayFile.close();
	//   if (fDecayFilePath == ""){ // LV from TGenPhaseSpace will be deleted in TGenPhaseSpace

	delete fDecayPhaseSpace;

	delete fLv8He;
	delete fLv4He;

	delete fLv6He;
	delete fLvn1;
	delete fLvn2;

	delete fLv6HeCMdecay;
	delete fLvn1CMdecay;
	delete fLvn2CMdecay;

	delete fLvNN;
}

//-------------------------------------------------------------------------------------------------
void ERDecay8He4He4nTransfer::Set8HeExcitation(Double_t excMean, Double_t fwhm, Double_t distibWeight)
{
	f8HeExcitationStateMeanEnergy.push_back(excMean);
	f8HeExcitationStateSigma.push_back(fwhm / 2.355);
	if (!fIs8HeExcitationSet)
	{
		f8HeExcitationStateWeight.push_back(distibWeight);
		fIs8HeExcitationSet = true;
		return;
	}
	f8HeExcitationStateWeight.push_back(f8HeExcitationStateWeight.back() + distibWeight);
}

//-------------------------------------------------------------------------------------------------
void ERDecay8He4He4nTransfer::Print8HeExcitation()
{

	LOG(INFO) << "ERDecay8He4He4nTransfer::Print8HeExcitation" << FairLogger::endl;
	std::cout << "\tfIs8HeExcitationSet = " << fIs8HeExcitationSet << std::endl;

	for (Int_t i; i < f8HeExcitationStateMeanEnergy.size(); i++)
	{
		std::cout << "\t" << i << "\t" << f8HeExcitationStateMeanEnergy[i]
				  << "\t" << f8HeExcitationStateSigma[i]
				  << "\t" << f8HeExcitationStateWeight[i]
				  << std::endl;
	}
}

//-------------------------------------------------------------------------------------------------
Bool_t ERDecay8He4He4nTransfer::Init()
{

	LOG(INFO) << "[ERDecay8He4He4nTransfer::Init] started" << FairLogger::endl;

	// check if all particles participating in reaction (and decay) are in TDatabasePDG
	TParticlePDG *p8He = TDatabasePDG::Instance()->GetParticle("8He");
	p8He->Print();
	if (!p8He)
	{
		LOG(ERROR) << "-W- [ERDecay8He4He4nTransfer::Init]: Ion 8He not found in TDatabasePDG!" << FairLogger::endl;
		return kFALSE;
	}

	TParticlePDG *p6He = TDatabasePDG::Instance()->GetParticle("6He");
	p6He->Print();
	if (!p6He)
	{
		LOG(ERROR) << "-W- [ERDecay8He4He4nTransfer::Init]: Ion 6He not found in TDatabasePDG!" << FairLogger::endl;
		return kFALSE;
	}

	TParticlePDG *p4He = TDatabasePDG::Instance()->GetParticle(1000020040);
	// TParticlePDG *p4He = TDatabasePDG::Instance()->GetParticle("Alpha");
	p4He->Print();
	if (!p4He)
	{
		LOG(ERROR) << "-W- [ERDecay8He4He4nTransfer::Init]: Ion 4He not found in TDatabasePDG!" << FairLogger::endl;
		return kFALSE;
	}

	TParticlePDG *pNeutron = TDatabasePDG::Instance()->GetParticle("neutron");
	pNeutron->Print();
	if (!pNeutron)
	{
		LOG(ERROR) << "-W- [ERDecay8He4He4nTransfer::Init]: Neutron not found in TDatabasePDG!" << FairLogger::endl;
		return kFALSE;
	}

	CalculateTargetParameters();

	return kTRUE;
}

//-------------------------------------------------------------------------------------------------
Bool_t ERDecay8He4He4nTransfer::Stepping()
{
	if (!fDecayFinish && gMC->TrackPid() == 1000020080 && TString(gMC->CurrentVolName()).Contains(GetInteractionVolumeName()))
	{
		if (!fIsInterationPointFound)
		{
			if (!FindInteractionPoint())
			{
				fDecayFinish = kTRUE;
				return kTRUE;
			}
			else
			{
				fDistanceFromEntrance = 0;
			}
		}
		gMC->SetMaxStep(fMinStep);
		TLorentzVector curPos;
		gMC->TrackPosition(curPos);
		Double_t trackStep = gMC->TrackStep();
		fDistanceFromEntrance += trackStep;

		if (fDistanceFromEntrance > fDistanceToInteractPoint)
		{
			// 8He + 4He → 4He + 8He

			// beam
			TLorentzVector beam;
			gMC->TrackMomentum(beam);
			LOG(DEBUG2) << "[ERDecay8He4He4nTransfer::Stepping] Mass of beam: " << beam.M() << " GeV" << FairLogger::endl;

			if (beam.P() == 0) // TODO: solve this temporary fix of bug
			{				   // temporary fix of bug with zero kinetic energy
				return kTRUE;
			}

			// target
			Double_t mTarget = G4IonTable::GetIonTable()->GetIon(2, 4)->GetPDGMass() * 1e-3;
			TLorentzVector target(0., 0., 0., mTarget);
			TLorentzVector lvReaction;
			lvReaction = beam + target;

			const TVector3 boost = lvReaction.BoostVector(); // Get Pcm 3 vector
			Double_t ECM = 0;
			TLorentzVector beamCM, targetCM;
			beamCM = beam;
			targetCM = target;
			beamCM.Boost(-boost);
			targetCM.Boost(-boost);
			ECM = beamCM.E() + targetCM.E();

			Int_t reactionHappen = kFALSE;

			Double_t excited8HeMass;
			Int_t reactionAttempsCounter = 0;
			Double_t excitation = 0; // excitation energy
			while (reactionHappen == kFALSE)
			{ // while reaction condition is not fullfilled
				if (fIs8HeExcitationSet)
				{
					Double_t randWeight = gRandom->Uniform(0., f8HeExcitationStateWeight.back());
					Int_t distribNum = 0;
					//       // choose distribution by weight
					for (; distribNum < f8HeExcitationStateWeight.size(); distribNum++)
					{
						if (randWeight < f8HeExcitationStateWeight[distribNum])
						{
							break;
						}
					}
					excitation = gRandom->Gaus(f8HeExcitationStateMeanEnergy[distribNum], f8HeExcitationStateSigma[distribNum]);
				}
				excited8HeMass = G4IonTable::GetIonTable()->GetIonMass(2, 8) / 1000. + excitation;
				Double_t alhpaRecoilMass = G4IonTable::GetIonTable()->GetIonMass(2, 4) / 1000.;
				LOG(DEBUG) << "[ERDecay8He4He4nTransfer::Stepping] Mass of 8He* product is " << excited8HeMass << " GeV" << FairLogger::endl;
				if ((ECM - alhpaRecoilMass - excited8HeMass) > 0)
				{ // выход из цикла while для PhaseGenerator
					reactionHappen = kTRUE;
					LOG(DEBUG) << "[ERDecay8He4He4nTransfer::Stepping] Reaction happened." << FairLogger::endl;
				}
				reactionAttempsCounter++;
				if (reactionAttempsCounter > 1000)
				{
					LOG(DEBUG) << "[ERDecay8He4He4nTransfer::Stepping] Reaction is forbidden for this CM energy." << FairLogger::endl;
					fDecayFinish = kTRUE;
					return kTRUE;
				}
			} // while

			ReactionPhaseGenerator(ECM, excited8HeMass);
			fLv8He->Boost(boost);
			fLv4He->Boost(boost);

			// 8He → 6He + n + n
			if (!DecayPhaseGenerator(excitation))
			{
				LOG(WARNING) << "[ERDecay8He4He4nTransfer::Stepping] 3-body decay did not occur!!" << FairLogger::endl;
			}

			PushTracks(curPos);

			gMC->StopTrack();
			fDecayFinish = kTRUE;
			// TODO: check why MaxStep is set to 100. here
			gMC->SetMaxStep(100.);

			FairRunSim *run = FairRunSim::Instance();
			if (TString(run->GetMCEventHeader()->ClassName()).Contains("ERDecay8He4He4nTransferEventHeader"))
			{
				ERDecay8He4He4nTransferEventHeader *header = (ERDecay8He4He4nTransferEventHeader *)run->GetMCEventHeader();
				// TODO: check meaning of time here
				header->SetData(curPos.Vect(),
								beam, target, *fLv8He, *fLv4He,
								*fLv6He, *fLvn1, *fLvn2,
								*fLv6HeCMdecay, *fLvn1CMdecay, *fLvn2CMdecay,
								*fLvNN, fE_T,
								0., fTheta);
				header->SetTrigger(1);
			}
		}
	}
	return kTRUE;
}

//-------------------------------------------------------------------------------------------------
void ERDecay8He4He4nTransfer::BeginEvent()
{
	fDecayFinish = kFALSE;
	fIsInterationPointFound = kFALSE;

	fLv8He->SetXYZM(0., 0., 0., 0.);
	fLv4He->SetXYZM(0., 0., 0., 0.);

	fLv6He->SetPxPyPzE(0., 0., 0., 0.);
	fLvn1->SetPxPyPzE(0., 0., 0., 0.);
	fLvn2->SetPxPyPzE(0., 0., 0., 0.);

	fLv6HeCMdecay->SetPxPyPzE(0., 0., 0., 0.);
	fLvn1CMdecay->SetPxPyPzE(0., 0., 0., 0.);
	fLvn2CMdecay->SetPxPyPzE(0., 0., 0., 0.);

	fLvNN->SetPxPyPzE(0., 0., 0., 0);
}

//-------------------------------------------------------------------------------------------------
void ERDecay8He4He4nTransfer::FinishEvent()
{
	FairRunSim *run = FairRunSim::Instance();
	if (TString(run->GetMCEventHeader()->ClassName()).Contains("ERDecayMCEventHeader"))
	{
		ERDecayMCEventHeader *header = (ERDecayMCEventHeader *)run->GetMCEventHeader();
		header->Clear();
	}
	if (TString(run->GetMCEventHeader()->ClassName()).Contains("ERDecay8He4He4nTransferEventHeader"))
	{
		ERDecay8He4He4nTransferEventHeader *header = (ERDecay8He4He4nTransferEventHeader *)run->GetMCEventHeader();
		header->Clear();
	}
}

//-------------------------------------------------------------------------------------------------
void ERDecay8He4He4nTransfer::ReactionPhaseGenerator(Double_t Ecm, Double_t massOf8HeProduct)
{
	// particle 1: 4He (m1, E1), recoil
	// particle 2: 8He (m2, E2), product

	Double_t m1 = G4IonTable::GetIonTable()->GetIonMass(2, 4) * 1e-3;
	// Double_t m2 = G4IonTable::GetIonTable()->GetIonMass(2, 8) * 1e-3;
	Double_t m2 = massOf8HeProduct;

	LOG(DEBUG2) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] Mass of the recoil 4He: " << m1 << " GeV" << FairLogger::endl;
	LOG(DEBUG2) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] Mass of the product 8He*: " << m2 << " GeV" << FairLogger::endl;
	LOG(DEBUG3) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] Energy in CM system of the reaction: " << Ecm << " GeV" << FairLogger::endl;

	// Energy of 1-st particle in cm.
	// total energy of the first particle is calculated as
	Double_t E1 = 0.5 * (Ecm * Ecm + m1 * m1 - m2 * m2) / Ecm;
	LOG(DEBUG2) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] Total energy of 4He: " << E1 << " GeV" << FairLogger::endl;

	// Impulse in CM
	Double_t Pcm = TMath::Sqrt(E1 * E1 - m1 * m1);
	// Generate angles of particles in CM
	Double_t thetaCM;
	if (!fADInput)
	{ // if file with angular distribution isn't setted than isotropic distribution is generated
		LOG(DEBUG) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] default isotropic distribution" << FairLogger::endl;
		thetaCM = TMath::ACos(gRandom->Uniform(-1, 1));
	}
	else
	{
		LOG(DEBUG) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] distribution from defined fADFunction" << FairLogger::endl;
		thetaCM = fADFunction->GetRandom(fThetaMin, fThetaMax) * TMath::DegToRad();
	}
	fTheta = thetaCM;
	Double_t phi = gRandom->Uniform(0., 2. * TMath::Pi());
	TVector3 Pcmv;
	Pcmv.SetMagThetaPhi(Pcm, thetaCM, phi);

	fLv8He->SetXYZM(Pcmv(0), Pcmv(1), Pcmv(2), m2);
	fLv4He->SetXYZM(-Pcmv(0), -Pcmv(1), -Pcmv(2), m1);

	LOG(DEBUG2) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] Pcm(8He): " << fLv8He->P() << FairLogger::endl;
	LOG(DEBUG2) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] thetaCM(8He): " << fLv8He->Theta() << FairLogger::endl;
	LOG(DEBUG2) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] phi(8He): " << fLv8He->Phi() << FairLogger::endl;

	LOG(DEBUG2) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] Pcm(4He): " << fLv4He->P() << FairLogger::endl;
	LOG(DEBUG2) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] thetaCM(4He): " << fLv4He->Theta() << FairLogger::endl;
	LOG(DEBUG2) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] phi(4He): " << fLv4He->Phi() << FairLogger::endl;
}

//-------------------------------------------------------------------------------------------------
Bool_t ERDecay8He4He4nTransfer::DecayPhaseGenerator(const Double_t excitation)
{
	// if (fDecayFilePath == "")
	// { // if decay file not defined, per morm decay using phase space
	LOG(DEBUG) << "[ERDecay8He4He4nTransfer::DecayPhaseGenerator] excitation of 8He is set to "
			   << excitation * 1000. << " MeV" << FairLogger::endl;

	Double_t decayMasses[3];
	// mass of 6He
	decayMasses[0] = G4IonTable::GetIonTable()->GetIonMass(2, 6) / 1000.;
	LOG(DEBUG2) << "[ERDecay8He4He4nTransfer::DecayPhaseGenerator] mass of 6He is "
				<< decayMasses[0] << " GeV" << FairLogger::endl;
	// mass of neutron
	decayMasses[1] = G4IonTable::GetIonTable()->GetIonMass(0, 1) / 1000.;
	decayMasses[2] = G4IonTable::GetIonTable()->GetIonMass(0, 1) / 1000.;
	LOG(DEBUG2) << "[ERDecay8He4He4nTransfer::DecayPhaseGenerator] mass of neutron is " << decayMasses[1] << " GeV" << FairLogger::endl;

	Double_t decayThreshold = decayMasses[0] + decayMasses[1] + decayMasses[2];
	if (excitation + G4IonTable::GetIonTable()->GetIonMass(2, 8) / 1000. < decayThreshold)
	{
		LOG(WARNING) << "[ERDecay8He4He4nTransfer::DecayPhaseGenerator] excitation of 8He is below 3-body decay threshold: "
					 << excitation * 1000. << " MeV < "
					 << decayThreshold * 1000. - G4IonTable::GetIonTable()->GetIonMass(2, 8) << " MeV"
					 << FairLogger::endl;
		return kFALSE;
	}

	fDecayPhaseSpace->SetDecay(*fLv8He, 3, decayMasses);
	// TODO: Save the outpupt of Generate into eventHeader
	// and investigate its distribution;
	// Pay attention that TGenPhaseSpace is not perfectly
	// working
	fDecayPhaseSpace->Generate();

	*fLv6He = *(fDecayPhaseSpace->GetDecay(0));
	*fLvn1 = *(fDecayPhaseSpace->GetDecay(1));
	*fLvn2 = *(fDecayPhaseSpace->GetDecay(2));

	const TVector3 boostCMdecay = fLv8He->BoostVector();

	*fLv6HeCMdecay = *fLv6He;
	*fLvn1CMdecay = *fLvn1;
	*fLvn2CMdecay = *fLvn2;

	fLv6HeCMdecay->Boost(-boostCMdecay);
	fLvn1CMdecay->Boost(-boostCMdecay);
	fLvn2CMdecay->Boost(-boostCMdecay);

	// epsilon = E_{nn}/E_T
	fE_T = fLv6HeCMdecay->E() - fLv6HeCMdecay->M() + fLvn1CMdecay->E() - fLvn1CMdecay->M() + fLvn2CMdecay->E() - fLvn2CMdecay->M();
	*fLvNN = *fLvn1CMdecay + *fLvn2CMdecay;
	// Double_t epsilon = (fLvNN->E() - fLvNN->M())/fE_T;

	return kTRUE;
	// }
	// else {
	// 	return kFALSE;
	// }

	// if (fDecayFile.eof())
	// {
	// 	LOG(ERROR) << "Decay file finished! There are no more events in file " << fDecayFilePath
	// 			   << " to be processed." << FairLogger::endl;
	// 	return kFALSE;
	// }
	// std::string event_line;
	// std::getline(fDecayFile, event_line);
	// std::istringstream iss(event_line);
	// std::vector<std::string> outputs_components((std::istream_iterator<std::string>(iss)),
	// 											std::istream_iterator<std::string>());
	// if (outputs_components.size() < 4 * 3)
	// {
	// 	LOG(ERROR) << "Wrong components number in raw in decay file!" << FairLogger::endl;
	// 	return kFALSE;
	// }

	// // Fill momentum vectors in CM.
	// TVector3 pn1(std::stod(outputs_components[0]), std::stod(outputs_components[1]),
	// 			 std::stod(outputs_components[2]));
	// TVector3 pn2(std::stod(outputs_components[3]), std::stod(outputs_components[4]),
	// 			 std::stod(outputs_components[5]));
	// TVector3 pn3(std::stod(outputs_components[6]), std::stod(outputs_components[7]),
	// 			 std::stod(outputs_components[8]));
	// TVector3 pn4(std::stod(outputs_components[9]), std::stod(outputs_components[10]),
	// 			 std::stod(outputs_components[11]));
	// // Apply scale factor
	// const auto excitationScale = excitation > 0. ? sqrt(excitation / fDecayFileExcitation) : 1.;
	// const auto MeV2GeV = 1. / 1000.;
	// const auto scale = excitationScale * MeV2GeV;
	// pn1 *= scale;
	// pn2 *= scale;
	// pn3 *= scale;
	// pn4 *= scale;
	// const auto fill_output_lorentz_vectors_in_lab =
	// 	[this](TLorentzVector *lv, const TVector3 &p, const Double_t mass)
	// {
	// 	lv->SetXYZM(p.X(), p.Y(), p.Z(), mass);
	// 	lv->Boost(fLv4n->BoostVector());
	// };
	// //   fill_output_lorentz_vectors_in_lab(fLvn1, pn1, fn->Mass());
	// //   fill_output_lorentz_vectors_in_lab(fLvn2, pn2, fn->Mass());
	// //   fill_output_lorentz_vectors_in_lab(fLvn3, pn3, fn->Mass());
	// //   fill_output_lorentz_vectors_in_lab(fLvn4, pn4, fn->Mass());
	// return kTRUE;
}

//-------------------------------------------------------------------------------------------------
Double_t ERDecay8He4He4nTransfer::ADEvaluate(Double_t *x, Double_t *p)
{
	if (fADInput->IsZombie())
	{
		Error("ERDecay8He4He4nTransfer::ADEvaluate", "AD input was not loaded");
		return -1;
	}
	// on each step of creating distribution function returns interpolated value of input data
	return fADInput->Eval(x[0]);
}

//-------------------------------------------------------------------------------------------------
void ERDecay8He4He4nTransfer::SetAngularDistribution(TString ADFile)
{
	// set the fADFunction basing on points from the text file
	//
	// in case of DEBUG2, resulting function will be drawn

	TString ADFilePath = gSystem->Getenv("VMCWORKDIR");
	ADFilePath = ADFile;
	std::ifstream f;
	LOG(INFO) << "[ERDecay8He4He4nTransfer::SetAngularDistribution] loading of file: "
			  << ADFilePath << FairLogger::endl;
	f.open(ADFilePath.Data());
	if (!f.is_open())
	{
		LOG(FATAL) << "Can't open file " << ADFilePath << FairLogger::endl;
	}

	std::vector<Double_t> theta;
	std::vector<Double_t> sigma;
	LOG(DEBUG3) << "[ERDecay8He4He4nTransfer::SetAngularDistribution] angular points before reading the file: "
				<< theta.size() << FairLogger::endl;
	Int_t i = 0;
	std::string line;
	Double_t currTheta, currSigma; // theta and sigma values for current line
	while (std::getline(f, line))
	{
		if (line.empty() || line[0] == '#' || line.substr(0, 2) == "//")
		{
			continue;
		}

		// Create a stringstream to parse the line
		std::stringstream ss(line);

		// read two columns with angle and cross section
		if (ss >> currTheta >> currSigma) // check for demanded two column format
		{
			theta.push_back(currTheta);
			sigma.push_back(currSigma);
		}
		else // handling error
		{
			LOG(ERROR) << "[ERDecay8He4He4nTransfer::SetAngularDistribution] "
					   << "Invalid format in line: " << line << FairLogger::endl;
			theta.push_back(0.);
			sigma.push_back(0.);
		}
		LOG(DEBUG4) << "\tread from file: " << i << ": " << theta[i] << "\t" << sigma[i] << FairLogger::endl;

		i++;
	}
	LOG(DEBUG3) << "[ERDecay8He4He4nTransfer::SetAngularDistribution] angular points after reading the file: "
				<< theta.size() << FairLogger::endl;

	if (theta.size() != sigma.size())
	{
		LOG(ERROR) << "[ERDecay8He4He4nTransfer::SetAngularDistribution] "
				   << "Different numbers of angle and crosssection values in the file "
				   << ADFilePath << FairLogger::endl;
	}
	Int_t graphSize = static_cast<Int_t>(theta.size());
	fADInput = new TGraph(graphSize, theta.data(), sigma.data());
	if (fADInput->GetN() <= 0)
	{ // if there are no points in input file
		LOG(INFO) << "ERDecay8He4He4nTransfer::SetAngularDistribution: "
				  << "Too few inputs for creation of AD function!" << FairLogger::endl;
		return;
	}
	Double_t *angle = fADInput->GetX(); // get first column variables that contains angle

	// Creation of angular distribution function using class member function.
	// Constructor divides interval (0; fADInput->GetN()-1) into grid.
	// On each step of grid it calls ADEvaluate() to get interpolated values of input data.
	fThetaMin = angle[0];
	fThetaMax = angle[fADInput->GetN() - 1];
	LOG(DEBUG2) << "[ERDecay8He4He4nTransfer::SetAngularDistribution] ThetaMin: " << fThetaMin << "\tThetaMax: " << fThetaMax << FairLogger::endl;

	if (FairLogger::GetLogger()->IsLogNeeded(DEBUG3))
	{
		for (size_t j = 0; j < fADInput->GetN(); j++)
		{
			if (j < 15 || j > (fADInput->GetN() - 15))
				LOG(DEBUG3) << "\tread from file: " << j << ": " << theta[j] << "\t" << sigma[j] << FairLogger::endl;
		}
	}

	fADFunction = new TF1("angDistr", this, &ERDecay8He4He4nTransfer::ADEvaluate,
						  fThetaMin, fThetaMax, 0, "ERDecay8He4He4nTransfer", "ADEvaluate");

	// temporary workaround of errors
	ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");
	fADFunction->GetRandom(fThetaMin, fThetaMax);
	// ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("AdaptiveSingular");

	if (FairLogger::GetLogger()->IsLogNeeded(DEBUG2))
		fADFunction->Draw();
}

//-------------------------------------------------------------------------------------------------
void ERDecay8He4He4nTransfer::PushTracks(TLorentzVector currentPosition)
{

	Int_t He8BeamTrackNb, He8TrackNb, He4TrackNb;
	Int_t He6TrackNb, n1TrackNb, n2TrackNb;

	He8BeamTrackNb = gMC->GetStack()->GetCurrentTrackNumber();
	LOG(DEBUG2) << "[ERDecay8He4He4nTransfer::Stepping] He8BeamTrackNb " << He8BeamTrackNb << FairLogger::endl;

	//TODO: treat the different tracks for different scenarios

	// binary reaction
	PushTrack(1, He8BeamTrackNb,
			  G4IonTable::GetIonTable()->GetNucleusEncoding(2, 8),
			  fLv8He, currentPosition, gMC->TrackTime(), He8TrackNb);

	PushTrack(1, He8BeamTrackNb,
			  G4IonTable::GetIonTable()->GetNucleusEncoding(2, 4),
			  fLv4He, currentPosition, gMC->TrackTime(), He4TrackNb);

	// decay
	PushTrack(1, He8TrackNb,
			  G4IonTable::GetIonTable()->GetNucleusEncoding(2, 6),
			  //   G4IonTable::GetIonTable()->GetIon(2,6)->GetPDGEncoding(),
			  //   1000020060,
			  fLv6He, currentPosition, gMC->TrackTime(), He6TrackNb);

	PushTrack(1, He8TrackNb,
			  TDatabasePDG::Instance()->GetParticle("neutron")->PdgCode(),
			  fLvn1, currentPosition, gMC->TrackTime(), n1TrackNb);

	PushTrack(1, He8TrackNb,
			  G4Neutron::Neutron()->GetPDGEncoding(),
			  //   TDatabasePDG::Instance()->GetParticle("neutron")->PdgCode(),
			  fLvn2, currentPosition, gMC->TrackTime(), n2TrackNb);

}

//-------------------------------------------------------------------------------------------------
void ERDecay8He4He4nTransfer::PushTrack(Int_t toBeDone, Int_t parentId, Int_t pdgCode,
										TLorentzVector *particle, TLorentzVector currentPosition,
										Double_t time, Int_t &trackNumber)
{

	gMC->GetStack()->PushTrack(1, parentId, pdgCode,
							   particle->Px(), particle->Py(), particle->Pz(), particle->E(),
							   currentPosition.X(), currentPosition.Y(), currentPosition.Z(),
							   time, 0., 0., 0.,
							   kPDecay, trackNumber, 1, 0);
}

//-------------------------------------------------------------------------------------------------
ClassImp(ERDecay8He4He4nTransfer)