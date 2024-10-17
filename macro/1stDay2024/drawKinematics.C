#if !defined(__CLING__)

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TF1.h"
#include "TMath.h"

#include <iostream>

using std::cout;
using std::endl;

#endif

TTree *tr;

TCut cReaction = "MCEventHeader.fThetaCM>0.";

void InitTree(TString inputFile);

void DrawBeam();
void DrawKin();
void DrawDecay();
void DrawDecayCM();
void SetGraphColor(TString particle, Color_t col);
void SetAxesTitles(TString histTitle, TString xAxisTitle, TString yAxisTitle = "");

void drawKinematics()
{

	InitTree("rootfiles/sim.root");

	const Bool_t canvasBeam = 0;
	const Bool_t canvasKin = 1;

	Double_t w = 1800;
	Double_t h = 1200;

	// DrawBeam();
	DrawKin();
	// DrawDecay();
	// DrawDecayCM();
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

void DrawDecay()
{

	Double_t w = 1800;
	Double_t h = 1200;

	auto cDecay = new TCanvas("cDecay", "3-body decay", w, h);
	cDecay->Divide(3, 2);

	tr->SetAlias("He8", "MCEventHeader.fHe8");
	tr->SetAlias("He6", "MCEventHeader.fHe6");
	tr->SetAlias("n1", "MCEventHeader.fn1");
	tr->SetAlias("n2", "MCEventHeader.fn2");

	cDecay->cd(1);
	// gPad->Divide(2, 2);
	gPad->Divide(1, 3);
	cDecay->cd(1)->cd(1);
	tr->Draw("He8.Px() - (He6.Px() + n1.Px() + n2.Px())", "");

	cDecay->cd(1)->cd(2);
	tr->Draw("He8.Py() - (He6.Py() + n1.Py() + n2.Py())", "");

	cDecay->cd(1)->cd(3);
	tr->Draw("He8.Pz() - (He6.Pz() + n1.Pz() + n2.Pz())", "");

	cDecay->cd(2);
	gPad->Divide(2, 2);
	cDecay->cd(2)->cd(1);
	tr->Draw("He8.Phi()", "");

	cDecay->cd(2)->cd(2);
	tr->Draw("n1.Phi()", "");

	cDecay->cd(2)->cd(3);
	tr->Draw("n2.Phi()", "", "");

	cDecay->cd(3);
	gPad->Divide(3, 1);
	cDecay->cd(3)->cd(1);
	tr->Draw("He8.Theta()", "");

	cDecay->cd(3)->cd(2);
	tr->Draw("n1.Theta()", "");

	cDecay->cd(3)->cd(3);
	tr->Draw("n2.Theta()", "");

	cDecay->cd(4);
	tr->Draw("He6.E()-He6.Mag():He6.Theta()", "");

	cDecay->cd(5);
	tr->Draw("n1.E()-n1.Mag():n1.Theta()", "");

	cDecay->cd(6);
	tr->Draw("n2.E()-n2.Mag():n2.Theta()", "", "surf");
}

void DrawDecayCM()
{

	Double_t w = 1800;
	Double_t h = 1200;

	auto cDecayCM = new TCanvas("cDecayCM", "3-body decay", w, h);
	cDecayCM->Divide(3, 2);

	// tr->SetAlias("He8", "MCEventHeader.fHe6DecayCM");
	tr->SetAlias("He6cm", "MCEventHeader.fHe6DecayCM");
	tr->SetAlias("n1cm", "MCEventHeader.fn1DecayCM");
	tr->SetAlias("n2cm", "MCEventHeader.fn2DecayCM");

	TF1 *fitThetaCMD6He = new TF1("fitFunc6He", "[0] * sin(x)", 0., TMath::Pi());
	auto hThetaCMD6He = new TH1F("histTheta6He", "Theta_CM(6He)", 100, 0., TMath::Pi());
	hThetaCMD6He->SetLineWidth(3);
	hThetaCMD6He->SetXTitle("ThetaCM (rad)");

	TF1 *fitThetaCMDn1 = new TF1("fitFuncN1", "[0] * sin(x)", 0., TMath::Pi());
	auto hThetaCMDn1 = new TH1F("histThetaN1", "Theta_CM(n1)", 100, 0., TMath::Pi());
	hThetaCMDn1->SetLineWidth(3);
	hThetaCMDn1->SetXTitle("ThetaCM (rad)");

	TF1 *fitThetaCMDn2 = new TF1("fitFuncN2", "[0] * sin(x)", 0., TMath::Pi());
	auto hThetaCMDn2 = new TH1F("histThetaN2", "Theta_CM(n2)", 100, 0., TMath::Pi());
	hThetaCMDn2->SetLineWidth(3);
	hThetaCMDn2->SetXTitle("ThetaCM (rad)");

	cDecayCM->cd(1);
	tr->Draw(Form("He6cm.Theta() >> %s", hThetaCMD6He->GetName()), "");
	// fitThetaCMD6He->SetParameter(0, hThetaCMD6He->GetMaximum());
	hThetaCMD6He->Fit(fitThetaCMD6He);

	cDecayCM->cd(2);
	tr->Draw(Form("n1cm.Theta() >> %s", hThetaCMDn1->GetName()), "");
	// fitThetaCMD6He->SetParameter(0, hThetaCMDn1->GetMaximum());
	hThetaCMDn1->Fit(fitThetaCMDn1);

	cDecayCM->cd(3);
	tr->Draw(Form("n2cm.Theta() >> %s", hThetaCMDn2->GetName()), "");
	// tr->Draw("n2cm.Theta()", "");
	hThetaCMDn2->Fit(fitThetaCMDn2);

	cDecayCM->cd(4);
	// tr->SetAlias("He6_kin", "(He6cm.E() - He6cm.Mag())*1000.");
	// tr->SetAlias("n1_kin", "(n1cm.E() - n1cm.Mag())*1000.");
	// tr->SetAlias("n2_kin", "(n2cm.E() - n2cm.Mag())*1000.");
	tr->SetAlias("E_T", "MCEventHeader.fE_T*1000.");
	
	tr->Draw("E_T>>histET", "");
	auto hET = (TH1F *)gPad->GetPrimitive("histET");
	hET->Fit("gaus");

	// tr->Draw("MCEventHeader.fE_T*1000.", "", "same");
	// tr->Draw("E_T - MCEventHeader.fE_T*1000.", "", "");
	// tr->Draw("He6cm.E() - He6cm.Mag()", "");

	cDecayCM->cd(5);
	// tr->SetAlias("Pnn", "TMath::Sqrt( TMath::Power(n1cm.Px()+n2cm.Px(),2)"
	// 								 "+ TMath::Power(n1cm.Py()+n2cm.Py(),2)"
	// 								 "+ TMath::Power(n1cm.Pz()+n2cm.Pz(),2) )");

	// (fLvNN->E() - fLvNN->M())/fE_T;

	TF1 *epsilon = new TF1("epsilon", "[0]*TMath::Sqrt(x*(1-x))", 0., 1.);
	// epsilon->Draw();
	// tr->Draw("(n1cm.E() - n1cm.Mag())*1000.", "");
	// tr->Draw("(MCEventHeader.fNNdecayCM.E() - MCEventHeader.fNNdecayCM.M())*1000./E_T", "");

	tr->Draw("(MCEventHeader.fNNdecayCM.M()-2*n2cm.M())*1000./E_T >> histEpsilon", "");
	auto hEpsilon = (TH1F *)gPad->GetPrimitive("histEpsilon");
	hEpsilon->SetLineColor(kGreen);
	hEpsilon->Fit(epsilon, "R");

	// tr->Draw("(MCEventHeader.fNNdecayCM.E() - MCEventHeader.fNNdecayCM.M())*1000./E_T", "");

	// tr->Draw("Pnn", "");
	// tr->Draw("((Pnn-TMath::Sqrt(n1cm.E()*n2cm.E()))/(E_T/1000.))", "");
	// tr->Draw("(Pnn-(n1cm.M()+n2cm.M()))/(E_T/1000.)", "");
	// return;
	cDecayCM->cd(6);
	tr->Draw("(n2cm.E() - n2cm.Mag())*1000.", "");
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