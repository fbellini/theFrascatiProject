/* akalweit@cern.ch, fbellini@cern.ch
   17.01.2018 - The Frascati plot
*/
#include "TMath.h"
#include "./generateBWpredictionsB2.C" //ADAPT ME

TGraphErrors * MakeB2TheoryGraphCoalescence(Double_t mT = 1.0, Double_t objSize = 3.2);
TGraphErrors * MakeB3TheoryGraphCoalescence(Double_t mT = 1.0, Double_t objSize = 1.75);

TGraphErrors * getB2_pp7TeV(Bool_t plotSys = 0, Double_t pToA = 0.75);
TGraphErrors * getB2_pPb5TeV(Bool_t plotSys = 0, Double_t pToA = 0.75);
TGraphErrors * getB2_PbPb5TeV(Bool_t plotSys = 0, Double_t pToA = 0.75);
TGraphErrors * getB2_PbPb276TeV(Bool_t plotSys = 0, Double_t pToA = 0.75);

TGraphErrors      * getB3_PbPb5TeV(Bool_t plotSys = 0, Double_t pToAb3 = 0.733);
TGraphErrors      * getB3_PbPb276TeV(Bool_t plotSys = 0, Double_t pToAb3 = 0.733);
TGraphAsymmErrors * getB3_pp7TeV(Bool_t plotSys = 0, Double_t pToAb3pp = 0.800);
TGraphAsymmErrors * getB3Lambda_PbPb276TeV(Bool_t plotSys = 0, Double_t pToAb3Lambda = 1.);

TGraphAsymmErrors * getBlastB2_PbPb276TeV(Bool_t plotSys = 0, Double_t pToA = 0.75);
TGraphAsymmErrors * getBlastB2_PbPb502TeV(Bool_t plotSys = 0, Double_t pToA = 0.75);
TGraphAsymmErrors * getBlastB2_pPb502TeV(Bool_t plotSys = 0, Double_t pToA = 0.75);
TGraphAsymmErrors * getBlastB2_pp7TeV(Bool_t plotSys = 0, Double_t pToA = 0.75);

TGraphAsymmErrors * getBlastB3_PbPb276TeV(Bool_t plotSys = 0, Double_t pToAb3 = 0.733);
TGraphAsymmErrors * getBlastB3Lambda_PbPb276TeV(Bool_t plotSys = 0, Double_t pToAb3Lambda = 1.0);

void MakeUp(TGraphErrors* obj, Color_t color, Color_t Fill_Color, Int_t Fill_Style, Int_t Line_Style, Int_t Line_Width, Int_t Marker_Style, Float_t Marker_Size);
void MakeUp(TGraphAsymmErrors* obj, Color_t color, Color_t Fill_Color, Int_t Fill_Style, Int_t Line_Style, Int_t Line_Width, Int_t Marker_Style, Float_t Marker_Size);


  
Int_t B2vsVolume(Bool_t plotLinX = 1, Double_t pToA = 0.75, Double_t pToAb3 = 0.733, Double_t pToAb3pp = 0.800, Double_t pToAb3Lambda = 1.)
{
  
  //--------------------
  //data
  //--------------------
  TGraphErrors* gB2vsR_pp7TeV = (TGraphErrors *) getB2_pp7TeV(kFALSE, pToA);
  TGraphErrors* gB2vsR_pp7TeV_sys = (TGraphErrors *) getB2_pp7TeV(kTRUE, pToA);

  TGraphErrors* gB2vsR_pPb5TeV = (TGraphErrors *) getB2_pPb5TeV(kFALSE, pToA);
  TGraphErrors* gB2vsR_pPb5TeV_sys = (TGraphErrors *) getB2_pPb5TeV(kTRUE, pToA);
  
  TGraphErrors* gB2vsR_PbPb5TeV = (TGraphErrors *) getB2_PbPb5TeV(kFALSE, pToA);
  TGraphErrors* gB2vsR_PbPb5TeV_sys = (TGraphErrors *) getB2_PbPb5TeV(kTRUE, pToA);
  
  TGraphErrors* gB2vsR_PbPb276TeV = (TGraphErrors *) getB2_PbPb276TeV(kFALSE, pToA);
  TGraphErrors* gB2vsR_PbPb276TeV_sys = (TGraphErrors *) getB2_PbPb276TeV(kTRUE, pToA);

  TGraphErrors* gB3vsR_PbPb5TeV = (TGraphErrors *) getB3_PbPb5TeV(kFALSE, pToAb3);
  TGraphErrors* gB3vsR_PbPb5TeV_sys = (TGraphErrors *) getB3_PbPb5TeV(kTRUE, pToAb3);
  
  TGraphErrors* gB3vsR_PbPb276TeV = (TGraphErrors *) getB3_PbPb276TeV(kFALSE, pToAb3);
  TGraphErrors* gB3vsR_PbPb276TeV_sys = (TGraphErrors *) getB3_PbPb276TeV(kTRUE, pToAb3);

  TGraphAsymmErrors* gB3vsR_pp7TeV = (TGraphAsymmErrors *) getB3_pp7TeV(kFALSE, pToAb3pp);
  TGraphAsymmErrors* gB3vsR_pp7TeV_sys = (TGraphAsymmErrors *) getB3_pp7TeV(kTRUE, pToAb3pp);

  TGraphAsymmErrors* gB3LambdavsR_PbPb276TeV = (TGraphAsymmErrors *) getB3Lambda_PbPb276TeV(kFALSE, pToAb3Lambda);
  TGraphAsymmErrors* gB3LambdavsR_PbPb276TeV_sys = (TGraphAsymmErrors *) getB3Lambda_PbPb276TeV(kTRUE, pToAb3Lambda);

  //-----------------------------
  //theory - Blast Wave + thermal
  //-----------------------------
  TGraphAsymmErrors* gBlastB2vsR_PbPb276TeV = (TGraphAsymmErrors *)  getBlastB2_PbPb276TeV(kFALSE, pToA);
  TGraphAsymmErrors* gBlastB2vsR_PbPb502TeV = (TGraphAsymmErrors *)  getBlastB2_PbPb502TeV(kFALSE, pToA);
  TGraphAsymmErrors* gBlastB2vsR_pPb502TeV = (TGraphAsymmErrors *)  getBlastB2_pPb502TeV(kFALSE, pToA);
  TGraphAsymmErrors* gBlastB2vsR_pp7TeV = (TGraphAsymmErrors *)  getBlastB2_pp7TeV(kFALSE, pToA);

  TGraphAsymmErrors* gBlastB3vsR_PbPb276TeV = (TGraphAsymmErrors *)  getBlastB3_PbPb276TeV(kFALSE, pToAb3);
  TGraphAsymmErrors* gBlastB3LambdavsR_PbPb276TeV = (TGraphAsymmErrors *)  getBlastB3Lambda_PbPb276TeV(kFALSE, pToAb3);
  
  //--------------------
  //theory - coalescence
  //--------------------
  //mT is the mass of the particle relative to which the HBT radius is calculated.
  //We use now the mass of the proton and the pT per nucleon
  Double_t mT = TMath::Sqrt(pToA * pToA + 0.938 * 0.938);
  // objSize = 3.2; //fm for the deuteron
  // objSize = 1.75; //fm for the 3^He
  
  TGraphErrors* hB2_coalescence = (TGraphErrors*) MakeB2TheoryGraphCoalescence(mT);
  hB2_coalescence->SetMarkerStyle(20);
  hB2_coalescence->SetLineWidth(3);

  TGraphErrors* hB2_coalescence_pointlike = (TGraphErrors*) MakeB2TheoryGraphCoalescence(mT, 0.0);
  hB2_coalescence_pointlike->SetMarkerStyle(24);
  hB2_coalescence_pointlike->SetMarkerSize(0.4);
  hB2_coalescence_pointlike->SetLineWidth(3);
  hB2_coalescence_pointlike->SetLineColor(kGray);

  TGraphErrors* hB3_coalescence = (TGraphErrors*) MakeB3TheoryGraphCoalescence(mT);
  hB3_coalescence->SetMarkerStyle(20);
  hB3_coalescence->SetMarkerColor(kOrange+3);
  hB3_coalescence->SetLineColor(kOrange+3);
  hB3_coalescence->SetLineWidth(3);

  TGraphErrors* hB3_coalescence_pointlike = (TGraphErrors*) MakeB3TheoryGraphCoalescence(mT, 0.0);
  hB3_coalescence_pointlike->SetMarkerStyle(24);
  hB3_coalescence_pointlike->SetMarkerColor(kOrange);
  hB3_coalescence_pointlike->SetLineColor(kOrange);
  hB3_coalescence_pointlike->SetMarkerSize(0.4);
  hB3_coalescence_pointlike->SetLineWidth(3);

  TGraphErrors* hB3L_coalescence = (TGraphErrors*) MakeB3TheoryGraphCoalescence(mT, 10.6);
  hB3L_coalescence->SetMarkerStyle(20);
  hB3L_coalescence->SetMarkerColor(kAzure-7);
  hB3L_coalescence->SetLineColor(kAzure-7);
  hB3L_coalescence->SetLineWidth(3);
  hB3L_coalescence->SetLineStyle(9);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.02); 

  //------------------------------
  // PLOT
  //------------------------------

  //make up options
  Int_t Fill_Style = 1001;
  Int_t Line_Style = 1;
  Int_t Line_Style_Blast = 2;
  Int_t Line_Width = 3;
  Float_t Marker_Size = 1.3;

  enum EPlotEntries { kPP7, kPPB502, kPBPB276, kPBPB502,
		      kPP7blast, kPPB502blast, kPBPB276blast, kPBPB502blast,
		      kB3_PP7, kB3_PPB502, kB3_PBPB276, kB3_PBPB502,
		      kB3_PP7blast, kB3_PPB502blast, kB3_PBPB276blast, kB3_PBPB502blast};
  
  Color_t color[]      = {kRed-2, kRed-7, kRed+1, kRed+2,
			   kBlue-5, kBlue-7, kBlue, kBlue+2,
			   kGreen+4, kSpring+4, kSpring-1, kGreen+2,
			   kBlue-5, kBlue-7, kBlue, kBlue+2};
  
  Int_t Marker_Style[] = { 34, 22, 21, 20,
			   28, 26, 25, 24,
			   34, 22, 21, 20,
			   28, 26, 25, 24};
  
  MakeUp(gB2vsR_pp7TeV_sys, color[EPlotEntries::kPP7], color[EPlotEntries::kPP7], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kPP7], Marker_Size);
  MakeUp(gB2vsR_pp7TeV    , color[EPlotEntries::kPP7], color[EPlotEntries::kPP7], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kPP7], Marker_Size);

  MakeUp(gB2vsR_pPb5TeV_sys, color[EPlotEntries::kPPB502], color[EPlotEntries::kPPB502], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kPPB502], Marker_Size);
  MakeUp(gB2vsR_pPb5TeV    , color[EPlotEntries::kPPB502], color[EPlotEntries::kPPB502], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kPPB502], Marker_Size);

  MakeUp(gB2vsR_PbPb276TeV_sys, color[EPlotEntries::kPBPB276], color[EPlotEntries::kPBPB276], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kPBPB276], Marker_Size);
  MakeUp(gB2vsR_PbPb276TeV    , color[EPlotEntries::kPBPB276], color[EPlotEntries::kPBPB276], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kPBPB276], Marker_Size);

  MakeUp(gB2vsR_PbPb5TeV_sys, color[EPlotEntries::kPBPB502], color[EPlotEntries::kPBPB502], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kPBPB502], Marker_Size);
  MakeUp(gB2vsR_PbPb5TeV    , color[EPlotEntries::kPBPB502], color[EPlotEntries::kPBPB502], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kPBPB502], Marker_Size);

  MakeUp(gBlastB2vsR_pp7TeV    , color[EPlotEntries::kPP7blast], color[EPlotEntries::kPP7blast], Fill_Style, Line_Style_Blast, Line_Width, Marker_Style[EPlotEntries::kPP7blast], Marker_Size);
  MakeUp(gBlastB2vsR_pPb502TeV    , color[EPlotEntries::kPPB502blast], color[EPlotEntries::kPPB502blast], Fill_Style, Line_Style_Blast, Line_Width, Marker_Style[EPlotEntries::kPPB502blast], Marker_Size);
  MakeUp(gBlastB2vsR_PbPb276TeV    , color[EPlotEntries::kPBPB276blast], color[EPlotEntries::kPBPB276blast], Fill_Style, Line_Style_Blast, Line_Width, Marker_Style[EPlotEntries::kPBPB276blast], Marker_Size);
  MakeUp(gBlastB2vsR_PbPb502TeV    , color[EPlotEntries::kPBPB502blast], color[EPlotEntries::kPBPB502blast], Fill_Style, Line_Style_Blast, Line_Width, Marker_Style[EPlotEntries::kPBPB502blast], Marker_Size);

  MakeUp(gB3vsR_pp7TeV_sys, color[EPlotEntries::kB3_PP7], color[EPlotEntries::kB3_PP7], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kB3_PP7], Marker_Size);
  MakeUp(gB3vsR_pp7TeV    , color[EPlotEntries::kB3_PP7], color[EPlotEntries::kB3_PP7], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kB3_PP7], Marker_Size);

  MakeUp(gB3vsR_PbPb276TeV_sys, color[EPlotEntries::kB3_PBPB276], color[EPlotEntries::kB3_PBPB276], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kB3_PBPB276], Marker_Size);
  MakeUp(gB3vsR_PbPb276TeV    , color[EPlotEntries::kB3_PBPB276], color[EPlotEntries::kB3_PBPB276], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kB3_PBPB276], Marker_Size);

  MakeUp(gB3vsR_PbPb5TeV_sys, color[EPlotEntries::kB3_PBPB502], color[EPlotEntries::kB3_PBPB502], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kB3_PBPB502], Marker_Size);
  MakeUp(gB3vsR_PbPb5TeV    , color[EPlotEntries::kB3_PBPB502], color[EPlotEntries::kB3_PBPB502], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kB3_PBPB502], Marker_Size);

  MakeUp(gBlastB3vsR_PbPb276TeV    , color[EPlotEntries::kB3_PBPB276blast], color[EPlotEntries::kB3_PBPB276blast], Fill_Style, Line_Style_Blast, Line_Width, Marker_Style[EPlotEntries::kB3_PBPB276blast], Marker_Size);
  
  //display
  TCanvas * cb2 = new TCanvas("cb2", "Frascati plot", 1600, 1000);
  cb2->SetBottomMargin(0.02);
  cb2->SetTopMargin(0.02);
  cb2->SetLeftMargin(0.15);
  cb2->SetRightMargin(0.02);
  TH2D * hframe = new TH2D("hframe", "B_{2} vs radius; radius (fm); #it{B}_{2} (GeV^{2}/#it{c}^{3})", 1000, 0.01, 6.0, 2000, 1.e-4, 0.1);
  hframe->GetXaxis()->SetTitleSize(0.05);
  hframe->GetYaxis()->SetTitleSize(0.05);
  if (plotLinX) hframe->GetXaxis()->SetRangeUser(0.01, 8.5);
  else  hframe->GetXaxis()->SetRangeUser(0.1, 10.5);
  TH2D * hframe3 = new TH2D("hframe3", "B_{3} vs radius; radius (fm); #it{B}_{3} (GeV^{4}/#it{c}^{6})", 1000, 0.01, 6.0, 2000, 1.e-9, 1.e-1);
  hframe3->GetXaxis()->SetTitleSize(0.05);
  hframe3->GetYaxis()->SetTitleSize(0.05);
  if (plotLinX) hframe3->GetXaxis()->SetRangeUser(0.01, 8.5);
  else  hframe3->GetXaxis()->SetRangeUser(0.1, 10.5);
  
  //Legends
  int nl = 10;
  TLegend * legB2;
  if (plotLinX) legB2 = new TLegend(0.35, 0.95-nl*0.03, 0.6, 0.95);
  else legB2 = new TLegend(0.2, 0.15, 0.55, 0.15+nl*0.03);
  legB2->SetFillStyle(0);
  legB2->SetTextSize(0.035);
  legB2->SetBorderSize(0);
  legB2->SetTextSize(0.025);
  
  legB2->AddEntry(hB2_coalescence, "#it{B}_{2} coalesc., #it{r}(d) = 3.2 fm", "l");
  legB2->AddEntry(hB2_coalescence_pointlike, "#it{B}_{2} coalesc., #it{r}(d) = 0 (point-like)", "l");
  legB2->AddEntry(gB2vsR_PbPb5TeV_sys, "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, prelim.", "pf");
  legB2->AddEntry(gB2vsR_PbPb276TeV_sys, "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV [PRC 93, 0249717 (2016)]", "pf");
  legB2->AddEntry(gB2vsR_pPb5TeV_sys, "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, prelim.", "pf");
  legB2->AddEntry(gB2vsR_pp7TeV_sys, "pp #sqrt{#it{s}} = 7 TeV, prelim.", "pf");
  legB2->AddEntry(gBlastB2vsR_PbPb502TeV, "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, BW + GSI (T = 156 MeV)", "l");
  legB2->AddEntry(gBlastB2vsR_PbPb276TeV, "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV, BW + GSI (T = 156 MeV)", "l");
  legB2->AddEntry(gBlastB2vsR_pPb502TeV, "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, BW + GSI (T = 156 MeV)", "l");
  legB2->AddEntry(gBlastB2vsR_pp7TeV, "pp #sqrt{#it{s}} = 7 TeV, BW + GSI (T = 156 MeV)", "l");
  
  nl = 8;
  TLegend * legB3;
  if (plotLinX) legB3 = new TLegend(0.35, 0.95-nl*0.03, 0.6, 0.95);
  else legB3 = new TLegend(0.2, 0.15, 0.55, 0.15+nl*0.03);
  legB3->SetFillStyle(0);
  legB3->SetTextSize(0.035);
  legB3->SetBorderSize(0);
  legB3->SetTextSize(0.025);
  legB3->AddEntry(hB3_coalescence, "#it{B}_{3} coalesc., #it{r}(^{3}He) = 1.75 fm", "l");
  legB3->AddEntry(hB3_coalescence_pointlike, "#it{B}_{3} coalesc., #it{r}(^{3}He) = 0 (point-like)", "l");
  legB3->AddEntry(hB3L_coalescence, "#it{B}_{3,#Lambda} coalesc., #it{r}(^{3}_{#Lambda}H) = 10.6 fm", "l");
  legB3->AddEntry(gB3vsR_PbPb5TeV_sys, "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, prelim.", "pf");
  legB3->AddEntry(gB3vsR_PbPb276TeV_sys, "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV [PRC 93, 0249717 (2016)]", "pf");
  legB3->AddEntry(gB3vsR_pp7TeV_sys, "pp #sqrt{#it{s}} = 7 TeV [arXiv:1709.08522]", "pf");
  legB3->AddEntry(gB3LambdavsR_PbPb276TeV_sys, "#it{B}_{3,#Lambda}, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV [PLB 754, 360-372 (2016)]", "pf");
  legB3->AddEntry(gBlastB3vsR_PbPb276TeV, "#it{B}_{3}, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV, BW + GSI (T = 156 MeV)", "l");
  legB3->AddEntry(gBlastB3LambdavsR_PbPb276TeV, "#it{B}_{3,#Lambda}, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV, BW + GSI (T = 156 MeV)", "l");

  cb2->Divide(2,1);

  //plot B2
  cb2->cd(1);
  gPad->SetLogy();
  if (!plotLinX) gPad->SetLogx();
  hframe->Draw();
  hB2_coalescence->Draw("l");
  hB2_coalescence_pointlike->Draw("lsame");
  gB2vsR_pp7TeV_sys->Draw("p3");
  gB2vsR_pp7TeV->Draw("samep");
  gB2vsR_pPb5TeV_sys->Draw("p3");
  gB2vsR_pPb5TeV->Draw("samep");
  gB2vsR_PbPb5TeV_sys->Draw("p3");
  gB2vsR_PbPb5TeV->Draw("samep");
  gB2vsR_PbPb276TeV_sys->Draw("p3");
  gB2vsR_PbPb276TeV->Draw("samep");
  gBlastB2vsR_PbPb276TeV->Draw("samel");
  gBlastB2vsR_PbPb502TeV->Draw("samel");
  gBlastB2vsR_pPb502TeV->Draw("samel");
  gBlastB2vsR_pp7TeV->Draw("samel");
  legB2->Draw();

  //plot B3
  cb2->cd(2);
  gPad->SetLogy();
  if (!plotLinX) gPad->SetLogx();
  hframe3->Draw();
  hB3_coalescence->Draw("l");
  hB3_coalescence_pointlike->Draw("lsame");
  hB3L_coalescence->Draw("lsame");
  gB3vsR_PbPb5TeV_sys->Draw("p3");
  gB3vsR_PbPb5TeV->Draw("samep");
  gB3vsR_PbPb276TeV_sys->Draw("samep3");
  gB3vsR_PbPb276TeV->Draw("samep");
  gB3vsR_pp7TeV_sys->Draw("samep2");
  gB3vsR_pp7TeV->Draw("samep");
  gB3LambdavsR_PbPb276TeV_sys->Draw("samep2");
  gB3LambdavsR_PbPb276TeV->Draw("samep");
  gBlastB3vsR_PbPb276TeV->Draw("samel");
  gBlastB3LambdavsR_PbPb276TeV->Draw("samel");
  legB3->Draw();

  //--------------------
  //Alternative plotting with legends on the side -- B2
  //--------------------
  nl = 5;
  TLegend * legB2data = new TLegend(0.1, 0.95-nl*0.04, 0.45, 0.95, "ALICE");
  legB2data->SetFillStyle(0);
  legB2data->SetTextSize(0.035);
  legB2data->SetBorderSize(0);
  legB2data->AddEntry(gB2vsR_PbPb5TeV_sys, "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, prelim.", "pf");
  legB2data->AddEntry(gB2vsR_PbPb276TeV_sys, "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV [PRC 93, 0249717 (2016)]", "pf");
  legB2data->AddEntry(gB2vsR_pPb5TeV_sys, "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, prelim.", "pf");
  legB2data->AddEntry(gB2vsR_pp7TeV_sys, "pp #sqrt{#it{s}} = 7 TeV, prelim.", "pf");
 
  nl = 5;
  TLegend * legB2blast = new TLegend(0.1, 0.65-nl*0.04, 0.45, 0.65, "Blast-Wave (#pi,K,p) + GSI-Heid. (T = 156 MeV)");
  legB2blast->SetFillStyle(0);
  legB2blast->SetTextSize(0.035);
  legB2blast->SetBorderSize(0);
  legB2blast->AddEntry(gBlastB2vsR_PbPb502TeV, "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV", "l");
  legB2blast->AddEntry(gBlastB2vsR_PbPb276TeV, "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "l");
  legB2blast->AddEntry(gBlastB2vsR_pPb502TeV, "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV", "l");
  legB2blast->AddEntry(gBlastB2vsR_pp7TeV, "pp #sqrt{#it{s}} = 7 TeV", "l");
  
  nl = 3;
  TLegend * legB2coal = new TLegend(0.1, 0.40-nl*0.04, 0.45, 0.40, "Coalescence");
  legB2coal->SetFillStyle(0);
  legB2coal->SetTextSize(0.035);
  legB2coal->SetBorderSize(0);
  legB2coal->AddEntry(hB2_coalescence, "#it{B}_{2} coalesc., #it{r}(d) = 3.2 fm", "l");
  legB2coal->AddEntry(hB2_coalescence_pointlike, "#it{B}_{2} coalesc., #it{r}(d) = 0 (point-like)", "l");
  
  TCanvas * cb2opta = new TCanvas("cb2opta", "Frascati plot B2", 1000, 600);
  cb2opta->SetBottomMargin(0.02);
  cb2opta->SetTopMargin(0.05);
  cb2opta->SetLeftMargin(0.15);
  cb2opta->SetRightMargin(0.02);
  
  cb2opta->Divide(2,1);
  cb2opta->cd(1);
  gPad->SetLogy();
  if (!plotLinX) gPad->SetLogx();
  hframe->Draw();
  hB2_coalescence->Draw("l");
  hB2_coalescence_pointlike->Draw("lsame");
  gB2vsR_pp7TeV_sys->Draw("p3");
  gB2vsR_pp7TeV->Draw("samep");
  gB2vsR_pPb5TeV_sys->Draw("p3");
  gB2vsR_pPb5TeV->Draw("samep");
  gB2vsR_PbPb5TeV_sys->Draw("p3");
  gB2vsR_PbPb5TeV->Draw("samep");
  gB2vsR_PbPb276TeV_sys->Draw("p3");
  gB2vsR_PbPb276TeV->Draw("samep");
  gBlastB2vsR_PbPb276TeV->Draw("samel");
  gBlastB2vsR_PbPb502TeV->Draw("samel");
  gBlastB2vsR_pPb502TeV->Draw("samel");
  gBlastB2vsR_pp7TeV->Draw("samel");

  cb2opta->cd(2);
  legB2data->Draw();
  legB2coal->Draw();
  legB2blast->Draw();

  //--------------------
  //Alternative plotting with legends on the side -- B3
  //--------------------
  nl = 8;
  TLegend * legB3data = new TLegend(0.1, 0.95-nl*0.04, 0.45, 0.95, "ALICE");
  legB3data->SetFillStyle(0);
  legB3data->SetTextSize(0.035);
  legB3data->SetBorderSize(0);
  legB3data->AddEntry(gB3vsR_PbPb5TeV_sys, "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, prelim.", "pf");
  legB3data->AddEntry(gB3vsR_PbPb276TeV_sys, "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV [PRC 93, 0249717 (2016)]", "pf");
  legB3data->AddEntry(gB3vsR_pp7TeV_sys, "pp #sqrt{#it{s}} = 7 TeV [arXiv:1709.08522]", "pf");
  legB3data->AddEntry(gB3LambdavsR_PbPb276TeV_sys, "#it{B}_{3,#Lambda}, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV [PLB 754, 360-372 (2016)]", "pf");

  nl = 2;
  TLegend * legB3blast = new TLegend(0.1, 0.55-nl*0.04, 0.45, 0.55, "Blast-Wave (#pi,K,p) + GSI-Heid. (T = 156 MeV)");
  legB3blast->SetFillStyle(0);
  legB3blast->SetTextSize(0.035);
  legB3blast->SetBorderSize(0);
  legB3blast->AddEntry(gBlastB3vsR_PbPb276TeV, "#it{B}_{3}, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "l");
  legB3blast->AddEntry(gBlastB3LambdavsR_PbPb276TeV, "#it{B}_{3,#Lambda}, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "l");

  nl = 4;
  TLegend * legB3coal = new TLegend(0.1, 0.4-nl*0.04, 0.45, 0.4, "Coalescence");
  legB3coal->SetFillStyle(0);
  legB3coal->SetTextSize(0.035);
  legB3coal->SetBorderSize(0);
  legB3coal->AddEntry(hB3_coalescence, "#it{B}_{3} coalesc., #it{r}(^{3}He) = 1.75 fm", "l");
  legB3coal->AddEntry(hB3_coalescence_pointlike, "#it{B}_{3} coalesc., #it{r}(^{3}He) = 0 (point-like)", "l");
  legB3coal->AddEntry(hB3L_coalescence, "#it{B}_{3,#Lambda} coalesc., #it{r}(^{3}_{#Lambda}H) = 10.6 fm", "l");

  TCanvas * cb3opta = new TCanvas("cb3opta", "Frascati plot B3", 1000, 600);
  cb3opta->SetBottomMargin(0.02);
  cb3opta->SetTopMargin(0.05);
  cb3opta->SetLeftMargin(0.15);
  cb3opta->SetRightMargin(0.02);
  
  cb3opta->Divide(2,1);
  cb3opta->cd(1);
  gPad->SetLogy();
  if (!plotLinX)  gPad->SetLogx();
  hframe3->Draw();
  hB3_coalescence->Draw("l");
  hB3_coalescence_pointlike->Draw("lsame");
  hB3L_coalescence->Draw("lsame");
  gB3vsR_PbPb5TeV_sys->Draw("p3");
  gB3vsR_PbPb5TeV->Draw("samep");
  gB3vsR_PbPb276TeV_sys->Draw("samep3");
  gB3vsR_PbPb276TeV->Draw("samep");
  gB3vsR_pp7TeV_sys->Draw("samep2");
  gB3vsR_pp7TeV->Draw("samep");
  gB3LambdavsR_PbPb276TeV_sys->Draw("samep2");
  gB3LambdavsR_PbPb276TeV->Draw("samep");
  gBlastB3vsR_PbPb276TeV->Draw("samel");
  gBlastB3LambdavsR_PbPb276TeV->Draw("samel");


  cb3opta->cd(2);
  legB3data->Draw();
  legB3coal->Draw();
  legB3blast->Draw();

  
  return 0;  
}


void getRadiusFromParameterisation(Double_t * multi, Double_t * radius)
{
  // parameterisation to convert multiplicity into Rside
  // this is the parameterisation for the pion radius
  //
  if (!multi || !radius) return;
  Double_t  multi3 = TMath::Power(multi[0], 1./3.);  
  //
  // Here is the crucial mapping between HBT radii and multi^(1/3)
  //
  // NEW VERSION (5th February 2018):
  // We assume 0.85fm for pp as Kfir (at dNdeta 6.01)
  // We assume R=4.5fm based on arXiv:1012.4035 (figure 2)
  // We interpolate linearly and take the highest kT bin, because
  // a pT/A of 0.8 GeV/c corresponds to mT=1.23GeV which is a kT of pions of 1.2GeV
  //
  Double_t radiusVal = 0.177825 + 0.36733 * multi3; 
  //
  //
  // OLD versions:
  // we consider Rside as a proxy for Rinv
  // the Rside vs dN/detaˆ1/3 is take from figure 9 of http://aliceinfo.cern.ch/ArtSubmission/node/1183
  // input has to be dNch/deta at midrapidity in VO* bins
  //
  //  Double_t radiusVal = -0.750 + 0.625 * multi3; //0.18 to fall on the b3 in pp
  //
  //
  Double_t radiusErr = radiusVal / 3.0 * multi[1] / multi[0];
  radius[0] = radiusVal;
  radius[1] = radiusErr;
  return; 
}



Double_t getB2fromRadius(Double_t homogR, Double_t mT, Double_t objSize)
{
 
  // formula 6.3, 4.12 from Scheibl, Heinz 1999 paper arXiv:nucl-th/9809092v2
  Double_t Rpar = homogR;
  Double_t Rperp = homogR;
  Double_t convFactor_fm2InvGeV = 0.197;

  Double_t invCd = (1.0 + TMath::Power( objSize / (2.0 * Rperp ), 2.0) ) *
		     TMath::Power((1.0 + TMath::Power( objSize / (2.0 * Rpar ), 2.0) ), 0.5);

  Double_t invCdApprox = TMath::Power( (1.0 + TMath::Power( objSize / (2.0 * Rperp ), 2.0)), 1.5);
  
  Double_t Cd = 1.0/invCd;
  Double_t B2 = (3.0 * TMath::Power(TMath::Pi(), 1.5) * Cd ) /
    (2 * mT * Rpar /convFactor_fm2InvGeV * Rperp /convFactor_fm2InvGeV * Rperp / convFactor_fm2InvGeV) ; 

  return B2;

}


Double_t getB3fromRadius(Double_t homogR, Double_t mT, Double_t objSize)
{
  // formula 9 of K. Blum, PRD 96 (2017) 103021
  // https://journals.aps.org/prd/pdf/10.1103/PhysRevD.96.103021
  Double_t R1 = homogR;
  Double_t R2 = homogR;
  Double_t R3 = homogR;

  Double_t convFactor_fm2InvGeV = 0.197;

  Double_t invC3Approx = (1.0 + TMath::Power( objSize / (2.0 * R1), 2.0)) *
    (1.0 + TMath::Power( objSize / (2.0 * R2), 2.0)) *
    (1.0 + TMath::Power( objSize / (2.0 * R3), 2.0));
  
  Double_t C3 = 1.0/invC3Approx;
  Double_t B3 = (TMath::Power(2* TMath::Pi(), 3.0) * C3) *
    TMath::Power(mT * R1 /convFactor_fm2InvGeV * R2 /convFactor_fm2InvGeV * R3 / convFactor_fm2InvGeV, -2.0) /
    ( 4 * TMath::Sqrt(3.0));
  
  return B3;
}


TGraphErrors * MakeB2TheoryGraphCoalescence(Double_t mT, Double_t objSize)
{
  const Int_t nPoints = 1000;
  Double_t gY[nPoints];
  Double_t gR[nPoints];
  
  TGraphErrors * graphOut = new TGraphErrors(nPoints);
  graphOut->SetName("B2_th_coalescence");
  graphOut->SetTitle("B_{2} from coalescence");
  
  
  for (int i = 0; i<nPoints; i++){
    gR[i] = 10.0 * i / nPoints; 
    gY[i] = getB2fromRadius(gR[i], mT, objSize);

    graphOut->SetPoint(i, gR[i], gY[i]);
    graphOut->SetPointError(i, 0.0, 0.0);
  }

  return graphOut;

}


TGraphErrors * MakeB3TheoryGraphCoalescence(Double_t mT, Double_t objSize)
{
  const Int_t nPoints = 1000;
  Double_t gY[nPoints];
  Double_t gR[nPoints];
  
  TGraphErrors * graphOut = new TGraphErrors(nPoints);
  graphOut->SetName("B3_th_coalescence");
  graphOut->SetTitle("B_{3} from coalescence");
  
  
  for (int i = 0; i<nPoints; i++){
    gR[i] = 10.0 * i / nPoints; 
    gY[i] = getB3fromRadius(gR[i], mT, objSize);

    graphOut->SetPoint(i, gR[i], gY[i]);
    graphOut->SetPointError(i, 0.0, 0.0);
  }

  return graphOut;

}

void convertMultiToRadius(TGraphErrors * graph)
{
  //convert multiplicity dN/deta into radius with parameterisation
  if (!graph) return;
  
  for (Int_t ip = 0; ip < graph->GetN(); ip++){
    Double_t xold[2];
    xold[0] = graph->GetX()[ip];
    xold[1] = graph->GetEX()[ip];
    Double_t xnew[2] = {-1.0, -1.0};
    getRadiusFromParameterisation(xold, xnew);
    graph->GetX()[ip] = xnew[0];
    graph->GetEX()[ip] = xnew[1];
  }

  return;
}


void convertMultiToRadius(TGraphAsymmErrors * graph)
{
  //convert multiplicity dN/deta into radius with parameterisation
  if (!graph) return;
  
  for (Int_t ip = 0; ip < graph->GetN(); ip++){
    Double_t xold[2];
    xold[0] = graph->GetX()[ip];
    xold[1] = graph->GetEXlow()[ip];
    Double_t xnew[2] = {-1.0, -1.0};
    getRadiusFromParameterisation(xold, xnew);
    graph->GetX()[ip] = xnew[0];
    graph->GetEXlow()[ip] = xnew[1];
    graph->GetEXhigh()[ip] = xnew[1];
	
  }

  return;
}

//---------------------------------------------------------
//------------------------------ ALICE data B2 --- pp 7 TeV

TGraphErrors * getB2_pp7TeV(Bool_t plotSys, Double_t pToA)
{
  //from Manuel - pp 7 TeV EXA 2017 preliminary
  TFile * f0 = TFile::Open("oB2vsNch.root");
  if (!f0) return NULL;

  TString gName = Form("B2_pToA=%.3f%s", pToA, (plotSys? "_sys" : ""));
  TGraphErrors * graph = (TGraphErrors*) f0->Get(gName.Data());
  if (!graph) {
    Printf("Error: cannot retrieve graph. Check pt bin requested.");
    return NULL;
  }

  convertMultiToRadius(graph);

  graph->SetMarkerColor(kGreen+2);
  graph->SetLineColor(kGreen+2);
  graph->SetFillColorAlpha(kGreen+2, 0.3);
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.2);
  graph->SetMarkerStyle(20);
  return graph;
  
}

//---------------------------------------------------------
//------------------------------ ALICE data B2 --- pPb 5 TeV
TGraphErrors * getB2_pPb5TeV(Bool_t plotSys, Double_t pToA)
{

  TFile * f0 = TFile::Open("oB2vsNch.root");
  if (!f0) return NULL;

  TString gName = Form("B2_pPb_pToA=%4.3f%s", pToA, (plotSys? "_sys" : ""));
  TGraphErrors * graph = (TGraphErrors*) f0->Get(gName.Data());
  if (!graph) {
    Printf("Error: cannot retrieve graph. Check pt bin requested.");
    return NULL;
  }

  convertMultiToRadius(graph);

  graph->SetMarkerColor(kBlue);
  graph->SetLineColor(kBlue);
  graph->SetFillColorAlpha(kBlue, 0.3);
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.2);
  graph->SetMarkerStyle(20);
  return graph;
  
}


//---------------------------------------------------------
//------------------------------ ALICE data B2 --- PbPb 5 TeV
TGraphErrors * getB2_PbPb5TeV(Bool_t plotSys, Double_t pToA)
{
  //Preliminary from Max - QM 2017
  TFile * f0 = TFile::Open("oB2vsNch.root");
  if (!f0) return NULL;

  TString gName = Form("B2_PbPb15_pToA=%4.3f%s", pToA, (plotSys? "_sys" : ""));
  TGraphErrors * graph = (TGraphErrors*) f0->Get(gName.Data());
  if (!graph) {
    Printf("Error: cannot retrieve graph. Check pt bin requested.");
    return NULL;
  }

  convertMultiToRadius(graph);

  graph->SetMarkerColor(kRed);
  graph->SetLineColor(kRed);
  graph->SetFillColorAlpha(kRed, 0.3);
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.2);
  graph->SetMarkerStyle(20);
  return graph;
  
}


//---------------------------------------------------------
//------------------------------ ALICE data B2 --- PbPb 2.76 TeV
TGraphErrors * getB2_PbPb276TeV(Bool_t plotSys, Double_t pToA)
{
  //Published PRC 93, 0249717 (2016)
  TFile * f0 = TFile::Open("oB2vsNch.root");
  if (!f0) return NULL;

  TString gName = Form("B2_PbPb10_pToA=%4.3f%s", pToA, (plotSys? "_sys" : ""));
  TGraphErrors * graph = (TGraphErrors*) f0->Get(gName.Data());
  if (!graph) {
    Printf("Error: cannot retrieve graph. Check pt bin requested.");
    return NULL;
  }

  convertMultiToRadius(graph);

  graph->SetMarkerColor(kRed+2);
  graph->SetLineColor(kRed+2);
  graph->SetFillColorAlpha(kRed+2, 0.3);  
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.2);
  graph->SetMarkerStyle(20);
  return graph;
  
}


//---------------------------------------------------------
//------------------------------ ALICE data B3 --- pp 7 TeV

TGraphAsymmErrors * getB3_pp7TeV(Bool_t plotSys, Double_t pToAb3pp)
{
  //from 1709.08522
  TFile * f0 = TFile::Open("B3pToA_pp7TeV.root");
  if (!f0) return NULL;

  TString gName = Form("B3_pp7TeV_pToA=%4.3f%s", pToAb3pp, (plotSys? "_sys" : ""));
  TGraphAsymmErrors * graph = (TGraphAsymmErrors*) f0->Get(gName.Data());
  if (!graph) {
    Printf("Error: cannot retrieve graph. Check pt bin requested.");
    return NULL;
  }

  convertMultiToRadius(graph);

  graph->SetMarkerColor(kTeal-5);
  graph->SetLineColor(kTeal-5);
  graph->SetFillColorAlpha(kTeal-5, 0.3);
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.2);
  graph->SetMarkerStyle(20);
  return graph;
  
}


//---------------------------------------------------------
//------------------------------ ALICE data B3 --- PbPb 5 TeV
TGraphErrors * getB3_PbPb5TeV(Bool_t plotSys, Double_t pToAb3)
{
  //Preliminary from Max - QM 2017
  TFile * f0 = TFile::Open("B3pToA_PbPb5TeV_preliminarySQM17.root");
  if (!f0) return NULL;

  TString gName = Form("B3_PbPb15_pToA=%4.3f%s", pToAb3, (plotSys? "_sys" : ""));
  TGraphErrors * graph = (TGraphErrors*) f0->Get(gName.Data());
  if (!graph) {
    Printf("Error: cannot retrieve graph. Check pt bin requested.");
    return NULL;
  }

  convertMultiToRadius(graph);

  graph->SetMarkerColor(kRed);
  graph->SetLineColor(kRed);
  graph->SetFillColorAlpha(kRed, 0.3);
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.2);
  graph->SetMarkerStyle(20);
  return graph;
  
}


//---------------------------------------------------------
//------------------------------ ALICE data B3 --- PbPb 2.76 TeV
TGraphErrors * getB3_PbPb276TeV(Bool_t plotSys, Double_t pToAb3)
{
  //Published PRC 93, 0249717 (2016)
  TFile * f0 = TFile::Open("B3pToA_PbPb276TeV.root");
  if (!f0) return NULL;

  TString gName = Form("B3_PbPb10_pToA=%4.3f%s", pToAb3, (plotSys? "_sys" : ""));
  TGraphErrors * graph = (TGraphErrors*) f0->Get(gName.Data());
  if (!graph) {
    Printf("Error: cannot retrieve graph. Check pt bin requested.");
    return NULL;
  }

  convertMultiToRadius(graph);
  graph->SetMarkerColor(kRed+2);
  graph->SetLineColor(kRed+2);
  graph->SetFillColorAlpha(kRed+2, 0.3);  
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.2);
  graph->SetMarkerStyle(20);
  return graph;
  
}

//---------------------------------------------------------
//------------------------------ ALICE data B3Lambda --- PbPb 2.76 TeV
TGraphAsymmErrors * getB3Lambda_PbPb276TeV(Bool_t plotSys, Double_t pToAb3Lambda)
{
  //Published Physics Letters B 754 (2016) 360–372
  TFile * f0 = TFile::Open("B3LambdapToA_PbPb276TeV.root");
  if (!f0) return NULL;
  TString gName = Form("B3Lambda_PbPb276TeV_pToA=%4.3f%s", pToAb3Lambda, (plotSys? "_sys" : ""));
  TGraphAsymmErrors * graph = (TGraphAsymmErrors*) f0->Get(gName.Data());
  if (!graph) {
    Printf("Error: cannot retrieve graph. Check pt bin requested.");
    return NULL;
  }

  convertMultiToRadius(graph);

  graph->SetMarkerColor(kAzure-7);
  graph->SetLineColor(kAzure-7);
  graph->SetFillColorAlpha(kAzure-7, 0.3);  
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.2);
  graph->SetMarkerStyle(34);
  return graph;
  
}


//---------------------------------------------------------
//---------------------- Blast wave + thermal PbPb 2.76 TeV
TGraphAsymmErrors * getBlastB2_PbPb276TeV(Bool_t plotSys, Double_t pToA)
{
  // Published proton yield from Phys. Rev. C 88 (2013) 044910
  // Blast wave params from π,K,p published
  // d/p from thermal model T = 156 MeV
 
  TGraphAsymmErrors* graph = (TGraphAsymmErrors *) generateBWpredictionsB2("PbPb276TeV", "rms", "deuteron", pToA);
  convertMultiToRadius(graph);
  
  graph->SetMarkerColor(kMagenta+2);
  graph->SetLineColor(kMagenta+2);
  graph->SetFillColorAlpha(kMagenta+2, 0.3);  
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.5);
  graph->SetMarkerStyle(33);
  graph->SetLineWidth(3);
  graph->SetLineStyle(2);
  return graph;
}

//---------------------------------------------------------
//---------------------- Blast wave + thermal PbPb 5.02 TeV
TGraphAsymmErrors * getBlastB2_PbPb502TeV(Bool_t plotSys, Double_t pToA)
{
  // Preliminary (~final) proton yields and BW params from:
  //https://gitlab.cern.ch/njacazio/SpectraAnalysisRun2/tree/master/results/spectra/spectra-pag/Preliminaries/QM2017
  // d/p from thermal model T = 156 MeV
 
  TGraphAsymmErrors* graph = (TGraphAsymmErrors *) generateBWpredictionsB2("PbPb502TeV", "rms", "deuteron", pToA);
  convertMultiToRadius(graph);
  
  graph->SetMarkerColor(kMagenta+2);
  graph->SetLineColor(kMagenta+2);
  graph->SetFillColorAlpha(kMagenta+2, 0.3);  
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.5);
  graph->SetMarkerStyle(33);
  graph->SetLineWidth(3);
  graph->SetLineStyle(2);
  return graph;
  
}

//---------------------------------------------------------
//---------------------- Blast wave + thermal pPb 5 TeV
TGraphAsymmErrors * getBlastB2_pPb502TeV(Bool_t plotSys, Double_t pToA)
{
  // Published proton yield  
  // Blast wave params from π,K,p published
  // d/p from thermal model T = 156 MeV
 
  TGraphAsymmErrors* graph = (TGraphAsymmErrors *) generateBWpredictionsB2("pPb502TeV", "rms", "deuteron", pToA);
  convertMultiToRadius(graph);
  
  graph->SetMarkerColor(kMagenta);
  graph->SetLineColor(kMagenta);
  graph->SetFillColorAlpha(kMagenta, 0.3);  
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.5);
  graph->SetMarkerStyle(33);
  graph->SetLineWidth(3);
  graph->SetLineStyle(2);
  return graph;
  
}

//---------------------------------------------------------
//---------------------- Blast wave + thermal pp 7 TeV
TGraphAsymmErrors * getBlastB2_pp7TeV(Bool_t plotSys, Double_t pToA)
{
  // Final proton yield from long paper pp7 TeV
  // Blast wave params from π,K,p published
  // d/p from thermal model T = 156 MeV
 
  TGraphAsymmErrors* graph = (TGraphAsymmErrors *) generateBWpredictionsB2("pp7TeV", "rms", "deuteron", pToA);
  convertMultiToRadius(graph);
  
  graph->SetMarkerColor(kMagenta-2);
  graph->SetLineColor(kMagenta-2);
  graph->SetFillColorAlpha(kMagenta-2, 0.3);  
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.5);
  graph->SetMarkerStyle(33);
  graph->SetLineWidth(3);
  graph->SetLineStyle(2);
  return graph;
  
}


//---------------------------------------------------------
//---------------------- Blast wave + thermal PbPb 2.76 TeV
TGraphAsymmErrors * getBlastB3_PbPb276TeV(Bool_t plotSys, Double_t pToAb3)
{
  // Published proton yield from Phys. Rev. C 88 (2013) 044910
  // Blast wave params from π,K,p published
  // 3He/p from thermal model T = 156 MeV
 
  TGraphAsymmErrors* graph = (TGraphAsymmErrors *) generateBWpredictionsB2("PbPb276TeV", "rms", "He3", pToAb3);
  convertMultiToRadius(graph);
  
  graph->SetMarkerColor(kMagenta+1);
  graph->SetLineColor(kMagenta+1);
  graph->SetFillColorAlpha(kMagenta+1, 0.3);  
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.5);
  graph->SetMarkerStyle(33);
  graph->SetLineWidth(3);
  graph->SetLineStyle(2);
  return graph;
  
}

//---------------------- Blast wave + thermal PbPb 2.76 TeV -- Hypertriton
TGraphAsymmErrors * getBlastB3Lambda_PbPb276TeV(Bool_t plotSys, Double_t pToAb3)
{
  // Published proton yield from Phys. Rev. C 88 (2013) 044910
  // Blast wave params from π,K,p published
  // 3He/p from thermal model T = 156 MeV
  // s3 from hyper-triton paper
 
  TGraphAsymmErrors* graph = (TGraphAsymmErrors *) generateBWpredictionsB2("PbPb276TeV", "rms", "hyper-triton", pToAb3);
  convertMultiToRadius(graph);
  
  graph->SetMarkerColor(kAzure+1);
  graph->SetLineColor(kAzure+1);
  graph->SetFillColorAlpha(kAzure+1, 0.3);  
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.5);
  graph->SetMarkerStyle(24);
  graph->SetLineWidth(3);
  graph->SetLineStyle(2);
  return graph;
  
}




//----------------------------------------------
void MakeUp(TGraphAsymmErrors* obj, Color_t color, Color_t Fill_Color, Int_t Fill_Style, Int_t Line_Style, Int_t Line_Width, Int_t Marker_Style, Float_t Marker_Size)
{
  if (!obj) return;
  obj->SetLineColor(color);
  obj->SetFillColor(Fill_Color);
  obj->SetFillStyle(Fill_Style);
  obj->SetFillColorAlpha(Fill_Color, 0.3);
  obj->SetLineStyle(Line_Style);
  obj->SetLineWidth(Line_Width);
  obj->SetMarkerColor(color);
  obj->SetMarkerStyle(Marker_Style);
  obj->SetMarkerSize(Marker_Size);
  return;
}

//----------------------------------------------
void MakeUp(TGraphErrors* obj, Color_t color, Color_t Fill_Color, Int_t Fill_Style, Int_t Line_Style, Int_t Line_Width, Int_t Marker_Style, Float_t Marker_Size)
{
  if (!obj) return;
  obj->SetLineColor(color);
  obj->SetFillColor(Fill_Color);
  obj->SetFillStyle(Fill_Style);
  obj->SetFillColorAlpha(Fill_Color, 0.3);
  obj->SetLineStyle(Line_Style);
  obj->SetLineWidth(Line_Width);
  obj->SetMarkerColor(color);
  obj->SetMarkerStyle(Marker_Style);
  obj->SetMarkerSize(Marker_Size);
  return;
}
