
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

ERDecay8He4He4nTransfer::ERDecay8He4He4nTransfer() : ERDecay("8He4He4nTransfer"),
													 fDecayFinish(kFALSE),
													 fTargetReactZ(0.),
													 fMinStep(0.01),
													//  f8He(NULL),
													//  f4He(NULL),
													 //   f6Li(NULL),
													 //   f4n (NULL),
													 //   fn  (NULL),
													 //  fIon8He(NULL),
													 //  fIon4He(NULL),
													 //   f4nMass(0.),
													 //   fIs4nUserMassSet(false),
													 //   fIs4nExcitationSet(false),
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
	fReactionPhaseSpace = new TGenPhaseSpace();
	//   fDecayPhaseSpace = new TGenPhaseSpace();
	FairRunSim *run = FairRunSim::Instance();
	//   fUnstable4n = new FairIon("4n",  0, 4, 0);
	//   fIon6Li     = new FairIon("6Li", 3, 6, 3);
	// fIon8He = new FairIon("8He", 2, 8, 2);
	// fIon4He = new FairIon("4He", 2, 4, 2);
	// std::cout << fIon4He->GetName() << std::endl;
	//   std::cout << "akdjfas" << std::endl;
	// run->AddNewIon(fIon8He);
	// run->AddNewIon(fIon4He);
	// std::cout << fIon4He->GetName() << std::endl;

	fLv8He = new TLorentzVector();
	fLv4He = new TLorentzVector();
}

//-------------------------------------------------------------------------------------------------
ERDecay8He4He4nTransfer::~ERDecay8He4He4nTransfer()
{
	//   if (fDecayFile.is_open())
	//     fDecayFile.close();
	//   if (fDecayFilePath == ""){ // LV from TGenPhaseSpace will be deleted in TGenPhaseSpace
	//     //   delete fLvn1;
	//     //   delete fLvn2;
	//     //   delete fLvn3;
	//     //   delete fLvn4;
	//   }
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

	// std::cout << "Decayer Init." << std::endl;
	LOG(INFO) << "[ERDecay8He4He4nTransfer::Init] started" << FairLogger::endl;

	if (!G4IonTable::GetIonTable()->GetIon(2,8))
	// f8He = TDatabasePDG::Instance()->GetParticle("8He");
	// f8He->Print();
	// if (!f8He)
	{
		std::cerr << "-W- [ERDecay8He4He4nTransfer::Init]: Ion 8He not found in database!" << std::endl;
		return kFALSE;
	}

	// f4He = TDatabasePDG::Instance()->GetParticle("4He");
	// std::cout << fIon4He->GetName() << std::endl;
	// f4He = TDatabasePDG::Instance()->GetParticle(fIon4He->GetName());
	// f4He = TDatabasePDG::Instance()->GetParticle("Alpha");
	// f4He->Print();
	// if (!f4He)
	// {
	// 	std::cerr << "-W- ERDecay8He4He4nTransfer: Ion Alpha not found in database!" << std::endl;
	// 	return kFALSE;
	// }

	// f4n = TDatabasePDG::Instance()->GetParticle("4n");
	// if (!f4n)
	// {
	// 	std::cerr << "-W- ERDecay2H_6Li: Ion 4n not found in database!" << std::endl;
	// 	return kFALSE;
	// }

	//   f6Li = TDatabasePDG::Instance()->GetParticle(fIon6Li->GetName());
	//   if ( ! f6Li ) {
	//     std::cerr  << "-W- ERDecay2H_6Li: Ion 6Li not found in database!" << std::endl;
	//     return kFALSE;
	//   }

	//   fn = TDatabasePDG::Instance()->GetParticle("neutron");
	//   if ( ! fn ) {
	//     std::cerr  << "-W- ERDecay2H_6Li: Particle neutron not found in database!" << std::endl;
	//     return kFALSE;
	//   }
	//   if (fIs4nUserMassSet) {
	//     fUnstable4n->SetMass(f4nMass / .931494028);
	//   } else {
	//     f4nMass = f4n->Mass(); // if user mass is not defined in ERDecay2H_6Li::SetH7Mass() than get a GEANT mass
	//   }

	// std::cout << "adfkjasdf1" << std::endl;
	CalculateTargetParameters();
	// std::cout << "adfkjasdf" << std::endl;

	//   if (fDecayFilePath != ""){
	//     LOG(INFO) << "Use decay kinematics from external text file" << FairLogger::endl;
	//     fDecayFile.open(fDecayFilePath.Data());
	//     if (!fDecayFile.is_open())
	//       LOG(FATAL) << "Can`t open decay file " << fDecayFilePath << FairLogger::endl;
	//     //Пропускаем шапку файла
	//     std::string header;
	//     std::getline(fDecayFile,header);

	//     // fLvn1 = new TLorentzVector();
	//     // fLvn2 = new TLorentzVector();
	//     // fLvn3 = new TLorentzVector();
	//     // fLvn4 = new TLorentzVector();
	//   }

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
			LOG(INFO) << "[ERDecay8He4He4nTransfer::Stepping] beam mass: " << beam.M() << FairLogger::endl;

			if (beam.P() == 0)
			{ // temporary fix of bug with zero kinetic energy
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
			// ECM = beamCM(3) + targetCM(3);
			// std::cout << ECM << std::endl;
			ECM = beamCM.E() + targetCM.E();
			// std::cout << ECM << std::endl;

			Int_t reactionHappen = kFALSE;

			//   Double_t decay4nMass;
			Double_t excited8HeMass;
			Int_t reactionAttempsCounter = 0;
			Double_t excitation = 0; // excitation energy
			while (reactionHappen == kFALSE)
			{ // while reaction condition is not fullfilled
				//     decay4nMass = f4nMass;
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
					// fUnstable4n->SetExcEnergy(excitation);
				}
				// Double_t Excited8HeMass = f8He->Mass() + excitation;
				// 
				excited8HeMass = G4IonTable::GetIonTable()->GetIonMass(2,8)/1000. + excitation;
				Double_t alhpaRecoilMass = G4IonTable::GetIonTable()->GetIonMass(2,4)/1000.;
				LOG(DEBUG) << "[ERDecay8He4He4nTransfer::Stepping] Excited8HeMass is " << excited8HeMass << FairLogger::endl;
				//     const float li6_mass = G4IonTable::GetIonTable()->GetIon(3,6)->GetPDGMass() * 1e-3;
				if ((ECM - alhpaRecoilMass - excited8HeMass) > 0)
				{ // выход из цикла while для PhaseGenerator
					reactionHappen = kTRUE;
					LOG(DEBUG) << "[ERDecay8He4He4nTransfer::Stepping] Reaction is happen" << FairLogger::endl;
				}
				reactionAttempsCounter++;
				if (reactionAttempsCounter > 1000)
				{
					LOG(DEBUG) << "[ERDecay8He4He4nTransfer::Stepping] Reaction is forbidden for this CM energy" << FairLogger::endl;
					fDecayFinish = kTRUE;
					return kTRUE;
				}
			}//while

			ReactionPhaseGenerator(ECM, excited8HeMass);
			fLv8He->Boost(boost);
			fLv4He->Boost(boost);

			//   //4n → n +n +n +n
			//   if (!DecayPhaseGenerator(excitation)){
			//     fDecayFinish = kTRUE;
			//     return kTRUE;
			//   }

			//   Int_t He8TrackNb, tetraNTrackNb, Li6TrackNb, n1TrackNb, n2TrackNb, n3TrackNb, n4TrackNb;
			Int_t He8BeamTrackNb, He8TrackNb, He4TrackNb;

			He8BeamTrackNb = gMC->GetStack()->GetCurrentTrackNumber();
			LOG(DEBUG2) << "[ERDecay8He4He4nTransfer::Stepping] He8BeamTrackNb " << He8BeamTrackNb << FairLogger::endl;

			// gMC->GetStack()->PushTrack(1, He8BeamTrackNb, f8He->PdgCode(),
			// 						   fLv8He->Px(), fLv8He->Py(), fLv8He->Pz(),
			// 						   fLv8He->E(), curPos.X(), curPos.Y(), curPos.Z(),
			// 						   gMC->TrackTime(), 0., 0., 0.,
			// 						   kPDecay, He8TrackNb, 1, 0);

			gMC->GetStack()->PushTrack(1, He8BeamTrackNb, G4IonTable::GetIonTable()->GetNucleusEncoding(2,8),
									   fLv8He->Px(), fLv8He->Py(), fLv8He->Pz(),
									   fLv8He->E(), curPos.X(), curPos.Y(), curPos.Z(),
									   gMC->TrackTime(), 0., 0., 0.,
									   kPDecay, He8TrackNb, 1, 0);

			gMC->GetStack()->PushTrack(1, He8BeamTrackNb, G4IonTable::GetIonTable()->GetNucleusEncoding(2,4),
									   fLv4He->Px(), fLv4He->Py(), fLv4He->Pz(),
									   fLv4He->E(), curPos.X(), curPos.Y(), curPos.Z(),
									   gMC->TrackTime(), 0., 0., 0.,
									   kPDecay, He4TrackNb, 1, 0);

			//   gMC->GetStack()->PushTrack(1, He8TrackNb, f6Li->PdgCode(),
			//                              fLv6Li->Px(), fLv6Li->Py(), fLv6Li->Pz(),
			//                              fLv6Li->E(), curPos.X(), curPos.Y(), curPos.Z(),
			//                              gMC->TrackTime(), 0., 0., 0.,
			//                              kPDecay, Li6TrackNb, f6Li->Mass(), 0);
			//   gMC->GetStack()->PushTrack(1, He8TrackNb, fn->PdgCode(),
			//                              fLvn1->Px(),fLvn1->Py(),fLvn1->Pz(),
			//                              fLvn1->E(), curPos.X(), curPos.Y(), curPos.Z(),
			//                              gMC->TrackTime(), 0., 0., 0.,
			//                              kPDecay, n1TrackNb, fn->Mass(), 0);
			//   gMC->GetStack()->PushTrack(1, He8TrackNb, fn->PdgCode(),
			//                              fLvn2->Px(),fLvn2->Py(),fLvn2->Pz(),
			//                              fLvn2->E(), curPos.X(), curPos.Y(), curPos.Z(),
			//                              gMC->TrackTime(), 0., 0., 0.,
			//                              kPDecay, n2TrackNb, fn->Mass(), 0);
			//   gMC->GetStack()->PushTrack(1, He8TrackNb, fn->PdgCode(),
			//                              fLvn3->Px(),fLvn3->Py(),fLvn3->Pz(),
			//                              fLvn3->E(), curPos.X(), curPos.Y(), curPos.Z(),
			//                              gMC->TrackTime(), 0., 0., 0.,
			//                              kPDecay, n3TrackNb, fn->Mass(), 0);
			//   gMC->GetStack()->PushTrack(1, He8TrackNb, fn->PdgCode(),
			//                              fLvn4->Px(),fLvn4->Py(),fLvn4->Pz(),
			//                              fLvn4->E(), curPos.X(), curPos.Y(), curPos.Z(),
			//                              gMC->TrackTime(), 0., 0., 0.,
			//                              kPDecay, n4TrackNb, fn->Mass(), 0);
			gMC->StopTrack();
			fDecayFinish = kTRUE;
			gMC->SetMaxStep(100.);

			FairRunSim *run = FairRunSim::Instance();
			if (TString(run->GetMCEventHeader()->ClassName()).Contains("ERDecayMCEventHeader"))
			{
				ERDecayMCEventHeader *header = (ERDecayMCEventHeader *)run->GetMCEventHeader();
				header->SetReactionPos(curPos.Vect());
				header->SetInputIon(He8BeamTrackNb);
				header->AddOutputParticle(He8TrackNb);
				header->AddOutputParticle(He4TrackNb);
				// header->AddOutputParticle(n1TrackNb);
				// header->AddOutputParticle(n2TrackNb);
				// header->AddOutputParticle(n3TrackNb);
				// header->AddOutputParticle(n4TrackNb);
			}
			if (TString(run->GetMCEventHeader()->ClassName()).Contains("ERDecay8He4He4nTransferEventHeader"))
			{
				ERDecay8He4He4nTransferEventHeader *header = (ERDecay8He4He4nTransferEventHeader *)run->GetMCEventHeader();
				header->SetData(curPos.Vect(), beam, target, *fLv8He, *fLv4He,
								// *fLvn1, *fLvn2, *fLvn3, *fLvn4,
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
	fTargetReactZ = fRnd->Uniform(-fTargetThickness / 2, fTargetThickness / 2);
	FairRunSim *run = FairRunSim::Instance();
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
	// particle 1: 4He (m1, E1)
	// particle 2: 8He (m2, E2)

	Double_t m1 = G4IonTable::GetIonTable()->GetIonMass(2, 4) * 1e-3;
	// Double_t m2 = G4IonTable::GetIonTable()->GetIonMass(2, 8) * 1e-3;
	Double_t m2 = massOf8HeProduct;

	LOG(DEBUG) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] m1(G4GetIonMass):\t " << m1 << FairLogger::endl;
	// LOG(DEBUG) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] m1(FairIon):\t " << fIon4He->GetMass() << FairLogger::endl;
	// LOG(DEBUG) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] m1(GetIonMass):\t " << G4IonTable::GetIonTable()->GetIonMass(2,4) << FairLogger::endl;
	// LOG(DEBUG) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] GetIonName:\t " << G4IonTable::GetIonTable()->GetIonName(2,8) << FairLogger::endl;
	LOG(DEBUG) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] m2: " << m2 << FairLogger::endl;
	// LOG(DEBUG) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] m2: " << fIon8He->GetMass() << FairLogger::endl;
	LOG(DEBUG) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] Ecm: " << Ecm << FairLogger::endl;

	// Energy of 1-st particle in cm.
	// total energy of the first particle is calculated as
	Double_t E1 = 0.5 * (Ecm * Ecm + m1 * m1 - m2 * m2) / Ecm;
	LOG(DEBUG) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] E1: " << E1 << FairLogger::endl;

	// Impulse in CM
	Double_t Pcm = TMath::Sqrt(E1 * E1 - m1 * m1);
	// Generate angles of particles in CM
	Double_t thetaCM;
	if (!fADInput)
	{ // if file with angular distribution isn't setted than isotropic distribution is generated
		// LOG(DEBUG) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] default isotropic distrubution" << FairLogger::endl;
		LOG(INFO) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] default isotropic distrubution" << FairLogger::endl;
		thetaCM = TMath::ACos(gRandom->Uniform(-1, 1));
	}
	else
	{
		LOG(INFO) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] distrubution from defined fADFunction" << FairLogger::endl;
		thetaCM = fADFunction->GetRandom(fThetaMin, fThetaMax) * TMath::DegToRad();
	}
	fTheta = thetaCM;
	Double_t phi = gRandom->Uniform(0., 2. * TMath::Pi());
	LOG(DEBUG) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] Pcm: " << Pcm << FairLogger::endl;
	LOG(DEBUG) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] theta: " << thetaCM << FairLogger::endl;
	LOG(DEBUG) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] phi: " << phi << FairLogger::endl;
	TVector3 Pcmv;
	Pcmv.SetMagThetaPhi(Pcm, thetaCM, phi);
	LOG(DEBUG) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] Pcm: " << Pcmv.Mag() << FairLogger::endl;
	LOG(DEBUG) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] thetaCM: " << Pcmv.Theta() << FairLogger::endl;
	LOG(DEBUG) << "[ERDecay8He4He4nTransfer::ReactionPhaseGenerator] phi: " << Pcmv.Phi() << FairLogger::endl;

	fLv8He->SetXYZM(0., 0., 0., 0.);
	fLv4He->SetXYZM(0., 0., 0., 0.);
	fLv8He->SetXYZM(Pcmv(0), Pcmv(1), Pcmv(2), m2);
	fLv4He->SetXYZM(-Pcmv(0), -Pcmv(1), -Pcmv(2), m1);
}

//-------------------------------------------------------------------------------------------------
// Bool_t ERDecay8He4He4nTransfer::DecayPhaseGenerator(const Double_t excitation) {
//   if (fDecayFilePath == ""){ // if decay file not defined, per morm decay using phase space
//     Double_t decayMasses[4];
//     // decayMasses[0] = fn->Mass();
//     // decayMasses[1] = fn->Mass();
//     // decayMasses[2] = fn->Mass();
//     // decayMasses[3] = fn->Mass();
//     fDecayPhaseSpace->SetDecay(*fLv4n, 4, decayMasses);
//     fDecayPhaseSpace->Generate();
//     // fLvn1 = fDecayPhaseSpace->GetDecay(0);
//     // fLvn2 = fDecayPhaseSpace->GetDecay(1);
//     // fLvn3 = fDecayPhaseSpace->GetDecay(2);
//     // fLvn4 = fDecayPhaseSpace->GetDecay(3);
//     return kTRUE;
//   }
//   if (fDecayFile.eof()){
//     LOG(ERROR) << "Decay file finished! There are no more events in file " << fDecayFilePath
//                << " to be processed." << FairLogger::endl;
//     return kFALSE;
//   }
//   std::string event_line;
//   std::getline(fDecayFile,event_line);
//   std::istringstream iss(event_line);
//   std::vector<std::string> outputs_components((std::istream_iterator<std::string>(iss)),
//                                                std::istream_iterator<std::string>());
//   if (outputs_components.size() < 4*3){
//     LOG(ERROR) << "Wrong components number in raw in decay file!" << FairLogger::endl;
//     return kFALSE;
//   }
//   // Fill momentum vectors in CM.
//   TVector3 pn1(std::stod(outputs_components[0]),std::stod(outputs_components[1]),
//                std::stod(outputs_components[2]));
//   TVector3 pn2(std::stod(outputs_components[3]),std::stod(outputs_components[4]),
//                std::stod(outputs_components[5]));
//   TVector3 pn3(std::stod(outputs_components[6]),std::stod(outputs_components[7]),
//                std::stod(outputs_components[8]));
//   TVector3 pn4(std::stod(outputs_components[9]),std::stod(outputs_components[10]),
//                std::stod(outputs_components[11]));
//   // Apply scale factor
//   const auto excitationScale = excitation > 0. ? sqrt(excitation / fDecayFileExcitation) : 1.;
//   const auto MeV2GeV = 1./1000.;
//   const auto scale = excitationScale * MeV2GeV;
//   pn1 *= scale;
//   pn2 *= scale;
//   pn3 *= scale;
//   pn4 *= scale;
//   const auto fill_output_lorentz_vectors_in_lab =
//       [this](TLorentzVector* lv, const TVector3& p, const Double_t mass) {
//         lv->SetXYZM(p.X(), p.Y(), p.Z(), mass);
//         lv->Boost(fLv4n->BoostVector());
//       };
// //   fill_output_lorentz_vectors_in_lab(fLvn1, pn1, fn->Mass());
// //   fill_output_lorentz_vectors_in_lab(fLvn2, pn2, fn->Mass());
// //   fill_output_lorentz_vectors_in_lab(fLvn3, pn3, fn->Mass());
// //   fill_output_lorentz_vectors_in_lab(fLvn4, pn4, fn->Mass());
//   return kTRUE;
// }

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
		// std::cout << line << std::endl;

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
ClassImp(ERDecay8He4He4nTransfer)