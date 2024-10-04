#if !defined(__CLING__)

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TLegend.h"

#include <iostream>

using std::cout;
using std::endl;

#endif

TTree *tr;

TCut cReaction = "MCEventHeader.fThetaCM>0.";

void InitTree(TString inputFile);

void DrawBeam();
void DrawKin();
void SetGraphColor(TString particle, Color_t col);
void SetAxesTitles(TString histTitle, TString xAxisTitle, TString yAxisTitle = "");

void drawKinematics()
{

	InitTree("rootfiles/sim.root");

	const Bool_t canvasBeam = 0;
	const Bool_t canvasKin = 1;

	Double_t w = 1800;
	Double_t h = 1200;

	DrawBeam();
	DrawKin();
}

void InitTree(TString inputFile)
{

	TFile *fr = new TFile(inputFile, "READ");
	if (!fr->IsOpen())
	{
		Error("InitTree", "File \"%s\" was not open.", inputFile.Data());
		return;
	}
	tr = (TTree *)fr->Get("er");
	if (tr->IsZombie())
	{
		Error("InitTree", "Tree \"kin\" was not open.");
		return;
	}

	tr->SetLineWidth(3);
	tr->SetLineColor(kBlue + 2);
}

void DrawBeam()
{

	Double_t w = 1800;
	Double_t h = 1200;

	// TH1F *currentHist;

	auto cBeam = new TCanvas("cBeam", "Beam properties", w, h);
	cBeam->Divide(3, 2);

	cBeam->cd(1);
	// tr->Draw("MCEventHeader.fHe8_beam.Theta()","","");
	// tr->Draw("MCEventHeader.fHe8_beam.Theta()", "MCEventHeader.fHe8_beam.Theta()>0.", "");
	tr->Draw("MCEventHeader.fHe8_beam.Theta()*TMath::RadToDeg()", cReaction, "");
	// cout << (TH1F*)gPad->GetPrimitive("htemp")->GetName() << endl;
	SetAxesTitles("Beam Theta", "Theta [deg]");
	// currentHist = (TH1F *)gPad->GetPrimitive("htemp");
	// currentHist->SetTitle("Beam Theta");
	// currentHist->GetXaxis()->SetTitle("Theta [deg]");

	cBeam->cd(2);
	tr->Draw("(MCEventHeader.fHe8_beam.E()-MCEventHeader.fHe8_beam.Mag())*1000", "MCEventHeader.fHe8_beam.Mag()>0.", "");
	SetAxesTitles("Beam Kinetic Energy", "T [MeV]");
	// currentHist = (TH1F *)gPad->GetPrimitive("htemp");
	// currentHist->SetTitle("Beam Kinetic Energy");
	// currentHist->GetXaxis()->SetTitle("T [MeV]");

	cBeam->cd(3);
	tr->Draw("(MCEventHeader.fHe8_beam.E()-MCEventHeader.fHe8_beam.Mag())*1000./8.", "MCEventHeader.fHe8_beam.Mag()>0.", "");
	SetAxesTitles("Beam Kinetic Energy", "T [MeV/A]");
	// currentHist = (TH1F *)gPad->GetPrimitive("htemp");
	// currentHist->SetTitle("Beam Kinetic Energy");
	// currentHist->GetXaxis()->SetTitle("T [MeV/A]");

	cBeam->cd(4);
	tr->Draw("fReactionPos.x()", cReaction, "");
	SetAxesTitles("Reaction Position", "x [cm]");
	// currentHist = (TH1F *)gPad->GetPrimitive("htemp");
	// currentHist->SetTitle("Reaction Position");
	// currentHist->GetXaxis()->SetTitle("x [cm]");

	cBeam->cd(5);
	tr->Draw("fReactionPos.y()", cReaction, "");
	SetAxesTitles("Reaction Position", "y [cm]");
	// currentHist = (TH1F *)gPad->GetPrimitive("htemp");
	// currentHist->SetTitle("Reaction Position");
	// currentHist->GetXaxis()->SetTitle("y [cm]");

	cBeam->cd(6);
	tr->Draw("fReactionPos.z()", cReaction, "");
	SetAxesTitles("Reaction Position", "z [cm]");
	// currentHist = (TH1F *)gPad->GetPrimitive("htemp");
	// currentHist->SetTitle("Reaction Position");
	// currentHist->GetXaxis()->SetTitle("z [cm]");
}

void DrawKin()
{

	// TH1F *currentHist;
	// TGraph *currentGraph;

	Double_t w = 1800;
	Double_t h = 1200;

	auto cKin = new TCanvas("cKin", "Kinematics", w, h);
	cKin->Divide(3, 2);

	cKin->cd(1);
	// {
	{
		gPad->Divide(2, 2);
		cKin->cd(1)->cd(1);
		tr->Draw("MCEventHeader.fThetaCM*TMath::RadToDeg()", "", "");
		SetAxesTitles("Reaction Theta", "", "Theta [deg]");
		// tr->Draw(Form("vBeam.E()-%f", massBeam), "");

		cKin->cd(1)->cd(2);
		tr->Draw("((MCEventHeader.fHe4.E())+(MCEventHeader.fHe8.E())-(MCEventHeader.fHe8_beam.E()+MCEventHeader.fHe4_target.E()))*1000.", cReaction, "");
		SetAxesTitles("Reaction Q", "Q [MeV]");
		// tr->Draw(Form("vRProductDcm.E()-%f", massRProduct), "");

		cKin->cd(1)->cd(3);
		tr->Draw("MCEventHeader.fHe8.M()", "");
		// tr->SetLineColor(kRed);
		// tr->Draw("vRProduct.Pz()+vRecoil.Pz()", "", "same");
		// tr->SetLineColor(kBlue + 2);

		cKin->cd(1)->cd(4);
		tr->Draw("MCEventHeader.fHe8.M()-MCEventHeader.fHe8_beam.M()", "");
		// tr->Draw("MCEventHeader.fHe8_beam.M()", "");

		cKin->Update();
	}

	// currentHist = (TH1F *)gPad->GetPrimitive("htemp");
	// currentHist->SetTitle("Reaction Theta");
	// currentHist->GetXaxis()->SetTitle("Theta [deg]");
	// }

	cKin->cd(4);
	{
		// currentHist = (TH1F *)gPad->GetPrimitive("htemp");
		// currentHist->SetTitle("Reaction Q");
	}

	cKin->cd(2);
	{
		tr->Draw("MCEventHeader.fHe4.Theta()*TMath::RadToDeg():MCEventHeader.fThetaCM*TMath::RadToDeg()", cReaction, "");
		SetGraphColor("4He", kRed);

		tr->Draw("MCEventHeader.fHe8.Theta()*TMath::RadToDeg():MCEventHeader.fThetaCM*TMath::RadToDeg()", cReaction, "same");
		SetGraphColor("8He", kBlue);

		SetAxesTitles("Binary Kinematics", "ThetaLab [deg]", "ThetaCM [deg]");
		gPad->BuildLegend();
	}

	cKin->cd(3);
	{
		tr->Draw("(MCEventHeader.fHe4.E()-MCEventHeader.fHe4.Mag())*1000./4.:MCEventHeader.fThetaCM*TMath::RadToDeg()", cReaction, "");
		SetGraphColor("4He", kRed);

		tr->Draw("(MCEventHeader.fHe8.E()-MCEventHeader.fHe8.Mag())*1000./8.:MCEventHeader.fThetaCM*TMath::RadToDeg()", cReaction, "same");
		SetGraphColor("8He", kBlue);

		SetAxesTitles("Binary Kinematics", "T [MeV]", "ThetaCM [deg]");
		gPad->BuildLegend();
	}

	cKin->cd(5);
	{
		tr->Draw("MCEventHeader.fHe4.Theta()*TMath::RadToDeg():MCEventHeader.fHe8.Theta()*TMath::RadToDeg()", cReaction, "");

		SetAxesTitles("Binary Kinematics", "ThetaLab (8He) [deg]", "ThetaLab (4He) [deg]");
	}

	cKin->cd(6);
	{
		tr->Draw("(MCEventHeader.fHe4.E()-MCEventHeader.fHe4.Mag())*1000./4.:MCEventHeader.fHe4.Theta()*TMath::RadToDeg()", cReaction, "");
		SetGraphColor("4He", kRed);

		tr->Draw("(MCEventHeader.fHe8.E()-MCEventHeader.fHe8.Mag())*1000./8.:MCEventHeader.fHe8.Theta()*TMath::RadToDeg()", cReaction, "same");
		SetGraphColor("8He", kBlue);

		SetAxesTitles("Binary Kinematics", "ThetaLab [deg]", "T_lab [MeV/A]");
		gPad->BuildLegend();
	}
}

void SetGraphColor(TString particle, Color_t col)
{

	auto currentGraph = (TGraph *)gPad->GetPrimitive("Graph");
	currentGraph->SetMarkerColor(col);
	currentGraph->SetLineColor(col);
	currentGraph->SetTitle(particle);
	currentGraph->SetName("g" + particle);
}

void SetAxesTitles(TString histTitle, TString xAxisTitle, TString yAxisTitle)
{

	auto currentHist = (TH1F *)gPad->GetPrimitive("htemp");
	currentHist->SetTitle(histTitle);
	currentHist->GetYaxis()->SetTitle(xAxisTitle);
	currentHist->GetXaxis()->SetTitle(yAxisTitle);
}