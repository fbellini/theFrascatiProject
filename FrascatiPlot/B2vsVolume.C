/*
   akalweit@cern.ch, fbellini@cern.ch
   17.01.2018 - The Frascati plot
   last update: 24.07.2020 by fbellini
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//The Frascati project - main figures
/////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TH2D.h"
#include "./GetCoalescence.cpp"
#include "./GetThermalBlastWave.cpp"
//#include "./GetALICEdata.cpp"

// beautification of plots
void MakeUp(TGraphErrors* obj, Color_t color, Color_t Fill_Color, Int_t Fill_Style, Int_t Line_Style, Int_t Line_Width, Int_t Marker_Style, Float_t Marker_Size);
void MakeUp(TGraphAsymmErrors* obj, Color_t color, Color_t Fill_Color, Int_t Fill_Style, Int_t Line_Style, Int_t Line_Width, Int_t Marker_Style, Float_t Marker_Size);

//figure making
void MakePaperFigure1(Bool_t plotLinX, Double_t pToA,
		      TF1 *Cd_coalescence, TF1* Cd_coalescence_pointlike, TF1* Cd_coalescence_radius1third,  TF1* Cd_coalescence_largeradius,
		      TGraphErrors * hB2_coalescence, TGraphErrors * hB2_coalescence_pointlike, TGraphErrors * hB2_coalescence_radius1third, TGraphErrors * hB2_coalescence_largeradius, Bool_t plotvert = 1);

void MakePaperFigure3(Bool_t plotLinX, Double_t pToA, Double_t pToAb3,
		      TGraphErrors * hB2_coalescence, TGraphErrors * hB3_coalescence, 
		      TGraphErrors ** gB2vsR_PbPb276TeV_sys,  TGraphErrors ** gB2vsR_pp7TeVINELg0_sys, TGraphErrors **gB2vsR_PbPb276TeV, TGraphErrors **  gB2vsR_pp7TeVINELg0,
		      TGraphErrors ** gB3vsR_PbPb276TeV_sys,  TGraphAsymmErrors ** gB3vsR_pp7TeV_sys, TGraphErrors ** gB3vsR_PbPb276TeV, TGraphAsymmErrors **  gB3vsR_pp7TeV);

void MakePaperFigure4(Bool_t plotLinX, Double_t pToA, Double_t pToAb3, Double_t pToAb3Lambda, 
		      TGraphErrors * hB2_coalescence, TGraphErrors * hB3_coalescence, TGraphErrors* hB3L_coalescence, TGraphErrors* hB3L_coalescence_largeradius, 
		      TGraphErrors ** gB2vsR_PbPb276TeV_sys,  TGraphErrors ** gB2vsR_pp7TeVINELg0_sys, TGraphErrors **gB2vsR_PbPb276TeV, TGraphErrors **  gB2vsR_pp7TeVINELg0,
		      TGraphErrors ** gB3vsR_PbPb276TeV_sys,  TGraphAsymmErrors ** gB3vsR_pp7TeV_sys, TGraphErrors ** gB3vsR_PbPb276TeV, TGraphAsymmErrors **  gB3vsR_pp7TeV,
		      TGraphAsymmErrors** gBlastB2vsR_PbPb276TeV,  TGraphAsymmErrors** gBlastB3vsR_PbPb276TeV,
		      TGraphAsymmErrors** gB3LambdavsR_PbPb276TeV, TGraphAsymmErrors** gB3LambdavsR_PbPb276TeV_sys, TGraphAsymmErrors** gBlastB3LambdavsR_PbPb276TeV);

//main plotting for Frascati plot
Int_t B2vsVolume(Bool_t plotLinX = 1, Double_t pToA = 0.75, Double_t pToAb3 = 0.733, Double_t pToAb3pp = 0.800, Double_t pToAb3Lambda = 1., Double_t pToAb4 = 0.75,
		 Double_t pToAb4Lambda = 0.75, Bool_t plotOnlyCoalescence = kFALSE, Bool_t plotPaperFigures = 0, Bool_t plotYRFigure = 0, Bool_t plotPseudoData = 0)
{
  //
  // main function which generates the plots of the Frascati project
  //

  //--------------------
  //data
  //--------------------
  const Int_t nParamSet = 4;
  TGraphErrors* gB2vsR_pp7TeV[nParamSet];
  TGraphErrors* gB2vsR_pp7TeV_sys[nParamSet];
  TGraphErrors* gB2vsR_pp7TeVINELg0[nParamSet];
  TGraphErrors* gB2vsR_pp7TeVINELg0_sys[nParamSet];
  TGraphErrors* gB2vsR_pp13TeV[nParamSet];
  TGraphErrors* gB2vsR_pp13TeV_sys[nParamSet];
  TGraphErrors* gB2vsR_pPb5TeV[nParamSet];
  TGraphErrors* gB2vsR_pPb5TeV_sys[nParamSet];
  TGraphErrors* gB2vsR_PbPb5TeV[nParamSet];
  TGraphErrors* gB2vsR_PbPb5TeV_sys[nParamSet];
  TGraphErrors* gB2vsR_PbPb276TeV[nParamSet];
  TGraphErrors* gB2vsR_PbPb276TeV_sys[nParamSet];
  //
  TGraphErrors* gB3vsR_PbPb5TeV[nParamSet];
  TGraphErrors* gB3vsR_PbPb5TeV_sys[nParamSet];
  TGraphErrors* gB3vsR_PbPb276TeV[nParamSet];
  TGraphErrors* gB3vsR_PbPb276TeV_sys[nParamSet];
  TGraphErrors* gB3vsR_pPb5TeV[nParamSet];
  TGraphErrors* gB3vsR_pPb5TeV_sys[nParamSet];
  TGraphAsymmErrors* gB3vsR_pp7TeV[nParamSet];
  TGraphAsymmErrors* gB3vsR_pp7TeV_sys[nParamSet];
  TGraphAsymmErrors* gB3LambdavsR_PbPb276TeV[nParamSet];
  TGraphAsymmErrors* gB3LambdavsR_PbPb276TeV_sys[nParamSet];
  //
  TGraphAsymmErrors* gBlastB2vsR_PbPb276TeV[nParamSet];
  TGraphAsymmErrors* gBlastB2vsR_PbPb502TeV[nParamSet];
  TGraphAsymmErrors* gBlastB2vsR_pPb502TeV[nParamSet];
  TGraphAsymmErrors* gBlastB3vsR_PbPb276TeV[nParamSet];
  TGraphAsymmErrors* gBlastB3LambdavsR_PbPb276TeV[nParamSet];
  TGraphAsymmErrors* gBlastB4vsR_PbPb276TeV[nParamSet];
  TGraphAsymmErrors* gBlastB4LambdavsR_PbPb276TeV[nParamSet];
  //
  TGraphAsymmErrors* gBlastB2vsR_pp7TeV[nParamSet];
  TGraphAsymmErrors* gBlastB3vsR_pp7TeV[nParamSet];
  TGraphAsymmErrors* gBlastB4vsR_pp7TeV[nParamSet];
  TGraphAsymmErrors* gBlastB3LambdavsR_pp7TeV[nParamSet];
  TGraphAsymmErrors* gBlastB4LambdavsR_pp7TeV[nParamSet];

  for (Int_t ip = 0; ip < nParamSet; ip++){
    gB2vsR_pp7TeVINELg0[ip] = (TGraphErrors *) getB2_pp7TeVINELg0(kFALSE, pToA, ip);
    gB2vsR_pp7TeVINELg0_sys[ip] = (TGraphErrors *) getB2_pp7TeVINELg0(kTRUE, pToA, ip);

    gB2vsR_pp7TeV[ip] = (TGraphErrors *) getB2_pp7TeV(kFALSE, pToA, ip);
    gB2vsR_pp7TeV_sys[ip] = (TGraphErrors *) getB2_pp7TeV(kTRUE, pToA, ip);

    gB2vsR_pp13TeV[ip] = (TGraphErrors *) getB2_pp13TeV(kFALSE, pToA, ip);
    gB2vsR_pp13TeV_sys[ip] = (TGraphErrors *) getB2_pp13TeV(kTRUE, pToA, ip);

    gB2vsR_pPb5TeV[ip] = (TGraphErrors *) getB2_pPb5TeV(kFALSE, pToA, ip);
    gB2vsR_pPb5TeV_sys[ip] = (TGraphErrors *) getB2_pPb5TeV(kTRUE, pToA, ip);
    gB2vsR_PbPb5TeV[ip] = (TGraphErrors *) getB2_PbPb5TeV(kFALSE, pToA, ip);
    gB2vsR_PbPb5TeV_sys[ip] = (TGraphErrors *) getB2_PbPb5TeV(kTRUE, pToA, ip);
    gB2vsR_PbPb276TeV[ip] = (TGraphErrors *) getB2_PbPb276TeV(kFALSE, pToA, ip);
    gB2vsR_PbPb276TeV_sys[ip] = (TGraphErrors *) getB2_PbPb276TeV(kTRUE, pToA, ip);
    //
    gB3vsR_PbPb5TeV[ip] = (TGraphErrors *) getB3_PbPb5TeV(kFALSE, pToAb3, ip);
    gB3vsR_PbPb5TeV_sys[ip] = (TGraphErrors *) getB3_PbPb5TeV(kTRUE, pToAb3, ip);
    gB3vsR_PbPb276TeV[ip] = (TGraphErrors *) getB3_PbPb276TeV(kFALSE, pToAb3, ip);
    gB3vsR_PbPb276TeV_sys[ip] = (TGraphErrors *) getB3_PbPb276TeV(kTRUE, pToAb3, ip);
    gB3vsR_pPb5TeV[ip] = (TGraphErrors *) getB3_pPb5TeV(kFALSE, pToAb3, ip);
    gB3vsR_pPb5TeV_sys[ip] = (TGraphErrors *) getB3_pPb5TeV(kTRUE, pToAb3, ip);
    gB3vsR_pp7TeV[ip] = (TGraphAsymmErrors *) getB3_pp7TeVINELg0(kFALSE, pToAb3pp, ip);
    gB3vsR_pp7TeV_sys[ip] = (TGraphAsymmErrors *) getB3_pp7TeVINELg0(kTRUE, pToAb3pp, ip);
    gB3LambdavsR_PbPb276TeV[ip] = (TGraphAsymmErrors *) getB3Lambda_PbPb276TeV(kFALSE, pToAb3Lambda, ip);
    gB3LambdavsR_PbPb276TeV_sys[ip] = (TGraphAsymmErrors *) getB3Lambda_PbPb276TeV(kTRUE, pToAb3Lambda, ip);

    //-----------------------------
    //theory - Blast Wave + thermal
    //-----------------------------
    gBlastB2vsR_PbPb276TeV[ip] = (TGraphAsymmErrors *)  getBlastB2_PbPb276TeV(kFALSE, pToA, ip);
    gBlastB2vsR_PbPb502TeV[ip] = (TGraphAsymmErrors *)  getBlastB2_PbPb502TeV(kFALSE, pToA, ip);
    gBlastB2vsR_pPb502TeV[ip] = (TGraphAsymmErrors *)  getBlastB2_pPb502TeV(kFALSE, pToA, ip);
    gBlastB2vsR_pp7TeV[ip] = (TGraphAsymmErrors *)  getBlastB2_pp7TeV(kFALSE, pToA, ip);
    //
    gBlastB3vsR_PbPb276TeV[ip] = (TGraphAsymmErrors *)  getBlastB3_PbPb276TeV(kFALSE, pToAb3, ip);
    gBlastB3vsR_pp7TeV[ip] = (TGraphAsymmErrors *)  getBlastB3_pp7TeV(kFALSE, pToAb3pp, ip);
    //
    gBlastB3LambdavsR_PbPb276TeV[ip] = (TGraphAsymmErrors *)  getBlastB3Lambda_PbPb276TeV(kFALSE, pToAb3Lambda, ip);
    gBlastB3LambdavsR_pp7TeV[ip] = (TGraphAsymmErrors *)  getBlastB3Lambda_pp7TeV(kFALSE, pToAb3Lambda, ip);
    //
    gBlastB4vsR_PbPb276TeV[ip] = (TGraphAsymmErrors *)  getBlastB4_PbPb276TeV(kFALSE, pToAb4, ip);
    gBlastB4vsR_pp7TeV[ip] = (TGraphAsymmErrors *)  getBlastB4_pp7TeV(kFALSE, pToAb4, ip);
    //
    gBlastB4LambdavsR_PbPb276TeV[ip] = (TGraphAsymmErrors *)  getBlastB4Lambda_PbPb276TeV(kFALSE, pToAb4Lambda, ip);
    gBlastB4LambdavsR_pp7TeV[ip] = (TGraphAsymmErrors *)  getBlastB4Lambda_pp7TeV(kFALSE, pToAb4Lambda, ip);
  }
  
  //--------------------
  //theory - coalescence
  //--------------------
  //mT is the mass of the particle relative to which the HBT radius is calculated.
  //We use now the mass of the proton and the pT per nucleon
  Double_t mT = TMath::Sqrt(pToA * pToA + 0.938 * 0.938);
  Double_t mT3 = TMath::Sqrt(pToAb3 * pToAb3 + 0.938 * 0.938);
  Double_t mT4 = TMath::Sqrt(pToAb4 * pToAb4 + 0.938 * 0.938);
  Double_t mT3L = TMath::Sqrt(pToAb3Lambda * pToAb3Lambda + 0.938 * 0.938);
  Double_t mT4L = TMath::Sqrt(pToAb4Lambda * pToAb4Lambda + 0.938 * 0.938);
  
  // objSize = 3.2; //fm for the deuteron
  // objSize = 2.48; //fm for the 3^He

  TF1 * Cd_coalescence = (TF1*) MakeB2TheoryGraphQMfactor();
  Cd_coalescence->SetLineWidth(3);
  Cd_coalescence->SetLineStyle(1);
  Cd_coalescence->SetLineColor(kBlack);

  TF1 * Cd_coalescence_pointlike = (TF1*) MakeB2TheoryGraphQMfactor(0.0);
  Cd_coalescence_pointlike->SetLineWidth(4);
  Cd_coalescence_pointlike->SetLineStyle(6);
  Cd_coalescence_pointlike->SetLineColor(kBlue);

  TF1 * Cd_coalescence_radius1third = (TF1*) MakeB2TheoryGraphQMfactor(0.3);
  Cd_coalescence_radius1third->SetLineWidth(2);
  Cd_coalescence_radius1third->SetLineStyle(7);
  Cd_coalescence_radius1third->SetLineColor(kRed);

  TF1 * Cd_coalescence_largeradius = (TF1*) MakeB2TheoryGraphQMfactor(10.0);
  Cd_coalescence_largeradius->SetLineWidth(2);
  Cd_coalescence_largeradius->SetLineStyle(9);
  Cd_coalescence_largeradius->SetLineColor(kGreen+1);
  
  TGraphErrors* hB2_coalescence = (TGraphErrors*) MakeB2TheoryGraphCoalescence(mT);
  hB2_coalescence->SetMarkerStyle(20);
  hB2_coalescence->SetLineWidth(3);
  hB2_coalescence->SetLineStyle(1);

  TGraphErrors* hB2_coalescence_pointlike = (TGraphErrors*) MakeB2TheoryGraphCoalescence(mT, 0.0);
  hB2_coalescence_pointlike->SetMarkerStyle(24);
  hB2_coalescence_pointlike->SetMarkerSize(0.4);
  hB2_coalescence_pointlike->SetLineWidth(4);
  hB2_coalescence_pointlike->SetLineStyle(6);
  hB2_coalescence_pointlike->SetLineColor(kBlue);

  TGraphErrors* hB2_coalescence_radius1third = (TGraphErrors*) MakeB2TheoryGraphCoalescence(mT, 0.3);
  hB2_coalescence_radius1third->SetMarkerStyle(24);
  hB2_coalescence_radius1third->SetMarkerSize(0.4);
  hB2_coalescence_radius1third->SetLineWidth(2);
  hB2_coalescence_radius1third->SetLineStyle(7);
  hB2_coalescence_radius1third->SetLineColor(kRed);

  TGraphErrors* hB2_coalescence_largeradius = (TGraphErrors*) MakeB2TheoryGraphCoalescence(mT, 10.0);
  hB2_coalescence_largeradius->SetMarkerStyle(24);
  hB2_coalescence_largeradius->SetMarkerSize(0.4);
  hB2_coalescence_largeradius->SetLineWidth(2);
  hB2_coalescence_largeradius->SetLineStyle(9);
  hB2_coalescence_largeradius->SetLineColor(kGreen+1);

  TGraphErrors* hB3_coalescence = (TGraphErrors*) MakeB3TheoryGraphCoalescence(mT3);
  hB3_coalescence->SetMarkerStyle(20);
  hB3_coalescence->SetMarkerColor(kBlack);
  hB3_coalescence->SetLineColor(kBlack);
  hB3_coalescence->SetLineWidth(3);

  TGraphErrors* hB3_coalescence_pointlike = (TGraphErrors*) MakeB3TheoryGraphCoalescence(mT3, 0.0);
  hB3_coalescence_pointlike->SetMarkerStyle(24);
  hB3_coalescence_pointlike->SetMarkerColor(kGray);
  hB3_coalescence_pointlike->SetLineColor(kGray);
  hB3_coalescence_pointlike->SetMarkerSize(0.4);
  hB3_coalescence_pointlike->SetLineWidth(2);

  TGraphErrors* hB3L_coalescence = (TGraphErrors*) MakeB3TheoryGraphCoalescence(mT3L, 6.8);
  hB3L_coalescence->SetMarkerStyle(20);
  hB3L_coalescence->SetMarkerColor(kBlack);
  hB3L_coalescence->SetLineColor(kBlack);
  hB3L_coalescence->SetLineWidth(3);
  hB3L_coalescence->SetLineStyle(1);

  TGraphErrors* hB3L_coalescence_largeradius = (TGraphErrors*) MakeB3TheoryGraphCoalescence(mT3L, 14.1);
  hB3L_coalescence_largeradius->SetMarkerStyle(20);
  hB3L_coalescence_largeradius->SetMarkerColor(kBlack);
  hB3L_coalescence_largeradius->SetLineColor(kBlack);
  hB3L_coalescence_largeradius->SetLineWidth(2);
  hB3L_coalescence_largeradius->SetLineStyle(7);

  TGraphErrors* hB4_coalescence = (TGraphErrors*) MakeB4TheoryGraphCoalescence(mT4, 1.9); //He4
  hB4_coalescence->SetMarkerStyle(1);
  hB4_coalescence->SetMarkerColor(kBlack);
  hB4_coalescence->SetLineColor(kBlack);
  hB4_coalescence->SetLineWidth(3);

  TGraphErrors* hB4L_coalescence = (TGraphErrors*) MakeB4TheoryGraphCoalescence(mT4L, 2.4); //4LH
  hB4L_coalescence->SetMarkerStyle(20);
  hB4L_coalescence->SetMarkerColor(kBlack);
  hB4L_coalescence->SetLineColor(kBlack);
  hB4L_coalescence->SetLineWidth(3);
  hB4L_coalescence->SetLineStyle(1);

  TGraphErrors* hB4L_coalescence_largeradius = (TGraphErrors*) MakeB4TheoryGraphCoalescence(mT4L, 4.9); //4LH
  hB4L_coalescence_largeradius->SetMarkerStyle(20);
  hB4L_coalescence_largeradius->SetMarkerColor(kBlack);
  hB4L_coalescence_largeradius->SetLineColor(kBlack);
  hB4L_coalescence_largeradius->SetLineWidth(2);
  hB4L_coalescence_largeradius->SetLineStyle(7);
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTopMargin(0.03);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05); 

  // ------------------
  // make the figure with coalescence only
  // ------------------
  MakePaperFigure1(plotLinX, pToA,
		   Cd_coalescence, Cd_coalescence_pointlike, Cd_coalescence_radius1third, Cd_coalescence_largeradius,
		   hB2_coalescence, hB2_coalescence_pointlike, hB2_coalescence_radius1third, hB2_coalescence_largeradius, kTRUE);
  if (plotOnlyCoalescence) return 0;

  //------------------------------
  // PLOT FRASCATI PLOT(S)
  //------------------------------
  //make up options
  Int_t Fill_Style = 1001;
  Int_t Line_Style = 1;
  Int_t Line_Style_Blast = 2;
  Int_t Line_Width = 1; 
  Int_t Line_Width_Blast = 3;
  Float_t Marker_Size = 1.3;

  enum EPlotEntries { kPP7, kPP13, kPPB502, kPBPB276, kPBPB502,
		      kPP7blast, kPPB502blast, kPBPB276blast, kPBPB502blast,
		      kB3_PP7, kB3_PPB502, kB3_PBPB276, kB3_PBPB502, kB3L_PBPB276,
		      kB3_PP7blast, kB3_PPB502blast, kB3_PBPB276blast, kB3_PBPB502blast, kB3L_PBPB276blast};
  
  Color_t color[]      = {kGreen+2, kOrange-3, kBlue+2, kRed, kRed+2,
			  kGreen+2, kBlue+2, kBlue+1, kRed+2,
			  kGreen+2, kBlue+2, kRed, kRed+2, kRed,
			  kBlue-5, kBlue-7, kBlue+1, kRed+2, kAzure+1};
  
  Int_t Marker_Style[] = { 21, 33, 22, 20, 23,
			   21, 22, 20, 23,
			   21, 22, 20, 23, 33,
			   21, 22, 20, 23, 33};
  
  for (Int_t ip = 0; ip < nParamSet; ip++) {
    MakeUp(gB2vsR_pp7TeV_sys[ip], color[EPlotEntries::kPP7], color[EPlotEntries::kPP7], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kPP7], Marker_Size);
    MakeUp(gB2vsR_pp7TeV[ip]    , color[EPlotEntries::kPP7], color[EPlotEntries::kPP7], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kPP7], Marker_Size);
    MakeUp(gB2vsR_pp7TeVINELg0_sys[ip], color[EPlotEntries::kPP7], color[EPlotEntries::kPP7], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kPP7], Marker_Size);
    MakeUp(gB2vsR_pp7TeVINELg0[ip]    , color[EPlotEntries::kPP7], color[EPlotEntries::kPP7], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kPP7], Marker_Size);
    MakeUp(gB2vsR_pp13TeV_sys[ip], color[EPlotEntries::kPP13], color[EPlotEntries::kPP13], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kPP13], Marker_Size);
    MakeUp(gB2vsR_pp13TeV[ip]    , color[EPlotEntries::kPP13], color[EPlotEntries::kPP13], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kPP13], Marker_Size);
    gB2vsR_pp13TeV_sys[ip]->SetMarkerSize(1.7);
    gB2vsR_pp13TeV[ip]->SetMarkerSize(1.7);
    MakeUp(gB2vsR_pPb5TeV_sys[ip], color[EPlotEntries::kPPB502], color[EPlotEntries::kPPB502], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kPPB502], Marker_Size);
    MakeUp(gB2vsR_pPb5TeV[ip]    , color[EPlotEntries::kPPB502], color[EPlotEntries::kPPB502], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kPPB502], Marker_Size);
    MakeUp(gB2vsR_PbPb276TeV_sys[ip], color[EPlotEntries::kPBPB276], color[EPlotEntries::kPBPB276], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kPBPB276], Marker_Size);
    MakeUp(gB2vsR_PbPb276TeV[ip]    , color[EPlotEntries::kPBPB276], color[EPlotEntries::kPBPB276], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kPBPB276], Marker_Size);
    MakeUp(gB2vsR_PbPb5TeV_sys[ip], color[EPlotEntries::kPBPB502], color[EPlotEntries::kPBPB502], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kPBPB502], Marker_Size);
    MakeUp(gB2vsR_PbPb5TeV[ip]    , color[EPlotEntries::kPBPB502], color[EPlotEntries::kPBPB502], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kPBPB502], Marker_Size);
    //blast 
    MakeUp(gBlastB2vsR_pp7TeV[ip]    , color[EPlotEntries::kPP7blast], color[EPlotEntries::kPP7blast], Fill_Style, Line_Style_Blast, Line_Width_Blast, Marker_Style[EPlotEntries::kPP7blast], Marker_Size);
    MakeUp(gBlastB2vsR_pPb502TeV[ip]    , color[EPlotEntries::kPPB502blast], color[EPlotEntries::kPPB502blast], Fill_Style, Line_Style_Blast, Line_Width_Blast, Marker_Style[EPlotEntries::kPPB502blast], Marker_Size);
    MakeUp(gBlastB2vsR_PbPb276TeV[ip]    , color[EPlotEntries::kPBPB276blast], color[EPlotEntries::kPBPB276blast], Fill_Style, Line_Style_Blast, Line_Width_Blast, Marker_Style[EPlotEntries::kPBPB276blast], Marker_Size);
    MakeUp(gBlastB2vsR_PbPb502TeV[ip]    , color[EPlotEntries::kPBPB502blast], color[EPlotEntries::kPBPB502blast], Fill_Style, Line_Style_Blast, Line_Width_Blast, Marker_Style[EPlotEntries::kPBPB502blast], Marker_Size);
    //B3
    MakeUp(gB3vsR_pp7TeV_sys[ip], color[EPlotEntries::kB3_PP7], color[EPlotEntries::kB3_PP7], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kB3_PP7], Marker_Size);
    MakeUp(gB3vsR_pp7TeV[ip]    , color[EPlotEntries::kB3_PP7], color[EPlotEntries::kB3_PP7], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kB3_PP7], Marker_Size);
    MakeUp(gB3vsR_PbPb276TeV_sys[ip], color[EPlotEntries::kB3_PBPB276], color[EPlotEntries::kB3_PBPB276], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kB3_PBPB276], Marker_Size);
    MakeUp(gB3vsR_PbPb276TeV[ip]    , color[EPlotEntries::kB3_PBPB276], color[EPlotEntries::kB3_PBPB276], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kB3_PBPB276], Marker_Size);
    MakeUp(gB3vsR_PbPb5TeV_sys[ip], color[EPlotEntries::kB3_PBPB502], color[EPlotEntries::kB3_PBPB502], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kB3_PBPB502], Marker_Size);
    MakeUp(gB3vsR_PbPb5TeV[ip]    , color[EPlotEntries::kB3_PBPB502], color[EPlotEntries::kB3_PBPB502], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kB3_PBPB502], Marker_Size);
    MakeUp(gB3vsR_pPb5TeV_sys[ip], color[EPlotEntries::kB3_PPB502], color[EPlotEntries::kB3_PPB502], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kB3_PPB502], Marker_Size);
    MakeUp(gB3vsR_pPb5TeV[ip]    , color[EPlotEntries::kB3_PPB502], color[EPlotEntries::kB3_PPB502], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kB3_PPB502], Marker_Size);
    MakeUp(gB3LambdavsR_PbPb276TeV_sys[ip], color[EPlotEntries::kB3L_PBPB276], color[EPlotEntries::kB3L_PBPB276], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kB3_PBPB276], Marker_Size);
    MakeUp(gB3LambdavsR_PbPb276TeV[ip]    , color[EPlotEntries::kB3L_PBPB276], color[EPlotEntries::kB3L_PBPB276], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kB3L_PBPB276], Marker_Size);
    //BLAST B3
    MakeUp(gBlastB3vsR_PbPb276TeV[ip]    , color[EPlotEntries::kB3_PBPB276blast], color[EPlotEntries::kB3_PBPB276blast], Fill_Style, Line_Style_Blast, Line_Width_Blast, Marker_Style[EPlotEntries::kB3_PBPB276blast], Marker_Size);
    MakeUp(gBlastB3LambdavsR_PbPb276TeV[ip] , color[EPlotEntries::kB3L_PBPB276blast], color[EPlotEntries::kB3L_PBPB276blast], Fill_Style, Line_Style_Blast, Line_Width_Blast, Marker_Style[EPlotEntries::kB3L_PBPB276blast], Marker_Size);
    MakeUp(gBlastB3vsR_pp7TeV[ip]    , color[EPlotEntries::kPP7blast], color[EPlotEntries::kPP7blast], Fill_Style, Line_Style_Blast, Line_Width_Blast, Marker_Style[EPlotEntries::kPP7blast], Marker_Size);
    MakeUp(gBlastB3LambdavsR_pp7TeV[ip]    , color[EPlotEntries::kPP7blast], color[EPlotEntries::kPP7blast], Fill_Style, Line_Style_Blast, Line_Width_Blast, Marker_Style[EPlotEntries::kPP7blast], Marker_Size);

    //blast B4
    MakeUp(gBlastB4vsR_PbPb276TeV[ip]    , color[EPlotEntries::kB3_PBPB276blast], color[EPlotEntries::kB3_PBPB276blast], Fill_Style, Line_Style_Blast, Line_Width_Blast, Marker_Style[EPlotEntries::kB3_PBPB276blast], Marker_Size);
    MakeUp(gBlastB4vsR_pp7TeV[ip]    , color[EPlotEntries::kPP7blast], color[EPlotEntries::kPP7blast], Fill_Style, Line_Style_Blast, Line_Width_Blast, Marker_Style[EPlotEntries::kPP7blast], Marker_Size);
    MakeUp(gBlastB4LambdavsR_PbPb276TeV[ip] , color[EPlotEntries::kB3L_PBPB276blast], color[EPlotEntries::kB3L_PBPB276blast], Fill_Style, Line_Style_Blast, Line_Width_Blast, Marker_Style[EPlotEntries::kB3L_PBPB276blast], Marker_Size);
    MakeUp(gBlastB4LambdavsR_pp7TeV[ip]    , color[EPlotEntries::kPP7blast], color[EPlotEntries::kPP7blast], Fill_Style, Line_Style_Blast, Line_Width_Blast, Marker_Style[EPlotEntries::kPP7blast], Marker_Size);

  }

  //---------------------------------------
  // PLOT FRASCATI PLOTS FOR PAPER
  //---------------------------------------  
  if (plotPaperFigures) {

    // MakePaperFigure3(plotLinX, pToA, pToAb3,
    // 		     hB2_coalescence, hB3_coalescence, 
    // 		     gB2vsR_PbPb276TeV_sys,  gB2vsR_pp7TeVINELg0_sys, gB2vsR_PbPb276TeV,  gB2vsR_pp7TeVINELg0,
    // 		     gB3vsR_PbPb276TeV_sys,  gB3vsR_pp7TeV_sys, gB3vsR_PbPb276TeV,  gB3vsR_pp7TeV);

    MakePaperFigure4(plotLinX, pToA, pToAb3, pToAb3Lambda,
		     hB2_coalescence, hB3_coalescence, hB3L_coalescence, hB3L_coalescence_largeradius,
		     gB2vsR_PbPb276TeV_sys,  gB2vsR_pp7TeVINELg0_sys, gB2vsR_PbPb276TeV,  gB2vsR_pp7TeVINELg0,
		     gB3vsR_PbPb276TeV_sys,  gB3vsR_pp7TeV_sys, gB3vsR_PbPb276TeV,  gB3vsR_pp7TeV,
		     gBlastB2vsR_PbPb276TeV, gBlastB3vsR_PbPb276TeV,
		     gB3LambdavsR_PbPb276TeV, gB3LambdavsR_PbPb276TeV_sys, gBlastB3LambdavsR_PbPb276TeV);

    return 0;
  }
    
  //---------------------------------------
  // PLOT FRASCATI PLOTS FOR SLIDES
  //---------------------------------------   
  TH2D * hframe = new TH2D("hframe", "B_{2} vs radius; #it{R} (fm); #it{B}_{2} (GeV^{2}/#it{c}^{3})", 1000, 0.01, 10.0, 2000, 1.e-4, 0.1);
  hframe->GetXaxis()->SetTitleSize(0.06);
  hframe->GetYaxis()->SetTitleSize(0.06);
  hframe->GetXaxis()->SetTitleOffset(0.8);
  hframe->GetXaxis()->SetLabelSize(0.05);
  hframe->GetYaxis()->SetLabelSize(0.05);
  if (plotLinX) hframe->GetXaxis()->SetRangeUser(0.01, 8.5);
  else  hframe->GetXaxis()->SetRangeUser(0.1, 10.5);

  TH2D * hframe3 = new TH2D("hframe3", "B_{3} vs radius; #it{R} (fm); #it{B}_{3} (GeV^{4}/#it{c}^{6})", 1000, 0.01, 10.0, 2000, 1.e-9, 1.e-1);
  hframe3->GetXaxis()->SetTitleSize(0.06);
  hframe3->GetYaxis()->SetTitleSize(0.06);
  hframe3->GetXaxis()->SetTitleOffset(0.8);
  hframe3->GetXaxis()->SetLabelSize(0.05);
  hframe3->GetYaxis()->SetLabelSize(0.05);
  if (plotLinX) hframe3->GetXaxis()->SetRangeUser(0.01, 8.5);
  else  hframe3->GetXaxis()->SetRangeUser(0.1, 10.5);

  //Define pT/A labels only once
  TPaveText * pavept = new TPaveText(0.17, 0.17, 0.7, 0.23, "NDC");
  pavept->SetFillStyle(0);
  pavept->SetTextFont(42);
  pavept->SetBorderSize(0);
  pavept->SetTextSize(0.05);
  pavept->SetTextAlign(12);
  pavept->AddText(Form("#it{B}_{2}: #it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToA));

  //TPaveText * paveptB3 = new TPaveText(0.55, 0.62, 0.95, 0.67, "NDC");
  TPaveText * paveptB3 = new TPaveText(0.17, 0.17, 0.7, 0.23, "NDC");
  paveptB3->SetFillStyle(0);
  paveptB3->SetTextFont(42);
  paveptB3->SetBorderSize(0);
  paveptB3->SetTextSize(0.05);
  paveptB3->SetTextAlign(12);
  paveptB3->AddText(Form("#it{B}_{3}: #it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb3));

  TPaveText * paveptB3L = new TPaveText(0.55, 0.55, 0.95, 0.6, "NDC");
  paveptB3L->SetFillStyle(0);
  paveptB3L->SetBorderSize(0);
  paveptB3L->SetTextFont(42);
  paveptB3L->SetTextSize(0.05);
  paveptB3L->SetTextAlign(12);
  paveptB3L->AddText(Form("#it{B}_{3,#Lambda}: #it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb3Lambda));


  //Chose the set of parameterisation for the radii
  Short_t ip = 3;

//Pave for parameteriszatino of R --> multi 
  TPaveText * paveparam = new TPaveText(0.1, 0.17, 0.7, 0.23, "NDC");
  paveparam->SetFillStyle(0);
  paveparam->SetTextFont(42);
  paveparam->SetBorderSize(0);
  paveparam->SetTextSize(0.05);
  paveparam->SetTextAlign(12);
  if (ip==1) paveparam->AddText(Form("R#rightarrow mult. constrained with ALICE B_{2}"));
  if (ip==0) paveparam->AddText(Form("R#rightarrow mult. fit to ALICE HBT"));
  if (ip==3) paveparam->AddText(Form("R#rightarrow mult. as in Ko et al., arXiv:1812.05175"));

 //display
  TCanvas * cb2 = new TCanvas("cb2", "Frascati plot", 1600, 1000);
  cb2->SetBottomMargin(0.02);
  cb2->SetTopMargin(0.02);
  cb2->SetLeftMargin(0.15);
  cb2->SetRightMargin(0.02);
  
  //Legends
  int nl = 10;
  TLegend * legB2;
  if (plotLinX) legB2 = new TLegend(0.42, 0.95-nl*0.03, 0.6, 0.95);
  else legB2 = new TLegend(0.2, 0.15, 0.55, 0.15+nl*0.03);
  legB2->SetFillStyle(0);
  legB2->SetTextSize(0.035);
  legB2->SetBorderSize(0);
  legB2->SetTextSize(0.025);
  
  legB2->AddEntry(hB2_coalescence, "#it{B}_{2} coalesc., #it{r}(d) = 3.2 fm", "l");
  legB2->AddEntry(hB2_coalescence_pointlike, "#it{B}_{2} coalesc., #it{r}(d) = 0 (point-like)", "l");
  legB2->AddEntry(gB2vsR_PbPb5TeV_sys[ip], "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, prelim.", "pf");
  legB2->AddEntry(gB2vsR_PbPb276TeV_sys[ip], "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV [PRC 93, 0249717 (2016)]", "pf");
  legB2->AddEntry(gB2vsR_pPb5TeV_sys[ip], "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, prelim.", "pf");
  legB2->AddEntry(gB2vsR_pp7TeV_sys[ip], "pp #sqrt{#it{s}} = 7 TeV [PRC 97, 024615 (2018)]", "pf");
  legB2->AddEntry(gBlastB2vsR_PbPb502TeV[ip], "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, BW + GSI (T = 156 MeV)", "l");
  legB2->AddEntry(gBlastB2vsR_PbPb276TeV[ip], "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV, BW + GSI (T = 156 MeV)", "l");
  legB2->AddEntry(gBlastB2vsR_pPb502TeV[ip], "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, BW + GSI (T = 156 MeV)", "l");
  legB2->AddEntry(gBlastB2vsR_pp7TeV[ip], "pp #sqrt{#it{s}} = 7 TeV, BW + GSI (T = 156 MeV)", "l");
  
  nl = 8;
  TLegend * legB3;
  if (plotLinX) legB3 = new TLegend(0.35, 0.95-nl*0.03, 0.6, 0.95);
  else legB3 = new TLegend(0.2, 0.15, 0.55, 0.15+nl*0.03);
  legB3->SetFillStyle(0);
  legB3->SetTextSize(0.035);
  legB3->SetBorderSize(0);
  legB3->SetTextSize(0.025);
  legB3->AddEntry(hB3_coalescence, "#it{B}_{3} coalesc., #it{r}(^{3}He) = 2.48 fm", "l");
  legB3->AddEntry(hB3_coalescence_pointlike, "#it{B}_{3} coalesc., #it{r}(^{3}He) = 0 (point-like)", "l");
  legB3->AddEntry(gB3vsR_PbPb5TeV_sys[ip], "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, prelim.", "pf");
  legB3->AddEntry(gB3vsR_PbPb276TeV_sys[ip], "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV [PRC 93, 0249717 (2016)]", "pf");
  legB3->AddEntry(gB3vsR_pPb5TeV_sys[ip], "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, 21.09.2018", "pf");
  legB3->AddEntry(gB3vsR_pp7TeV_sys[ip], "pp #sqrt{#it{s}} = 7 TeV [PRC 97, 024615 (2018)]", "pf");
  legB3->AddEntry(gBlastB3vsR_PbPb276TeV[ip], "#it{B}_{3}, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV, BW + GSI (T = 156 MeV)", "l");
  legB3->AddEntry(hB3L_coalescence, "#it{B}_{3,#Lambda} coalesc., #it{r}(^{3}_{#Lambda}H) = 6.8 fm", "l");
  legB3->AddEntry(gB3LambdavsR_PbPb276TeV_sys[ip], "#it{B}_{3,#Lambda}, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV [PLB 754, 360-372 (2016)]", "pf");
  legB3->AddEntry(gBlastB3LambdavsR_PbPb276TeV[ip], "#it{B}_{3,#Lambda}, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV, BW + GSI (T = 156 MeV)", "l");

  cb2->Divide(2,1);

  //plot B2
  cb2->cd(1);
  gPad->SetLogy();
  if (!plotLinX) gPad->SetLogx();
  hframe->Draw();
  hB2_coalescence->Draw("l");
  hB2_coalescence_pointlike->Draw("lsame");
  gBlastB2vsR_PbPb276TeV[ip]->Draw("samel");
  gBlastB2vsR_PbPb502TeV[ip]->Draw("samel");
  gBlastB2vsR_pPb502TeV[ip]->Draw("samel");
  gBlastB2vsR_pp7TeV[ip]->Draw("samel");
  gB2vsR_pp7TeV_sys[ip]->Draw("p3");
  gB2vsR_pp7TeV[ip]->Draw("samep");
  gB2vsR_pPb5TeV_sys[ip]->Draw("p3");
  gB2vsR_pPb5TeV[ip]->Draw("samep");
  gB2vsR_PbPb5TeV_sys[ip]->Draw("p3");
  gB2vsR_PbPb5TeV[ip]->Draw("samep");
  gB2vsR_PbPb276TeV_sys[ip]->Draw("p3");
  gB2vsR_PbPb276TeV[ip]->Draw("samep");
  legB2->Draw();
  pavept->Draw();

  //plot B3
  cb2->cd(2);
  gPad->SetLogy();
  if (!plotLinX) gPad->SetLogx();
  hframe3->Draw();
  hB3_coalescence->Draw("l");
  hB3_coalescence_pointlike->Draw("lsame");
  hB3L_coalescence->Draw("lsame");
  gB3vsR_PbPb5TeV_sys[ip]->Draw("p3");
  gB3vsR_PbPb5TeV[ip]->Draw("samep");
  gB3vsR_PbPb276TeV_sys[ip]->Draw("samep3");
  gB3vsR_PbPb276TeV[ip]->Draw("samep");
  gB3vsR_pPb5TeV_sys[ip]->Draw("p3");
  gB3vsR_pPb5TeV[ip]->Draw("samep");
  gB3vsR_pp7TeV_sys[ip]->Draw("samep2");
  gB3vsR_pp7TeV[ip]->Draw("samep");
  gB3LambdavsR_PbPb276TeV_sys[ip]->Draw("samep2");
  gB3LambdavsR_PbPb276TeV[ip]->Draw("samep");
  gBlastB3vsR_PbPb276TeV[ip]->Draw("samel");
  gBlastB3LambdavsR_PbPb276TeV[ip]->Draw("samel");
  legB3->Draw();
  paveptB3->Draw();
  paveptB3L->Draw();

  //--------------------
  //Alternative plotting with legends on the side -- B2
  //--------------------
  nl = 5;
  TLegend * legB2data = new TLegend(0.1, 0.95-nl*0.04, 0.45, 0.95, "ALICE");
  legB2data->SetFillStyle(0);
  legB2data->SetTextSize(0.035);
  legB2data->SetBorderSize(0);
  legB2data->AddEntry(gB2vsR_PbPb5TeV_sys[ip], "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, prelim.", "pf");
  legB2data->AddEntry(gB2vsR_PbPb276TeV_sys[ip], "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV [PRC 93, 0249717 (2016)]", "pf");
  legB2data->AddEntry(gB2vsR_pPb5TeV_sys[ip], "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, prelim.", "pf");
  legB2data->AddEntry(gB2vsR_pp7TeV_sys[ip], "pp #sqrt{#it{s}} = 7 TeV, paper in prep.", "pf");
  legB2data->AddEntry(gB2vsR_pp13TeV_sys[ip], "pp #sqrt{#it{s}} = 13 TeV, prelim.", "pf");
  //legB2data->AddEntry(gB2vsR_pp7TeV_inelg0_sys[ip], "pp #sqrt{#it{s}} = 7 TeV [PRC 97, 024615 (2018)]", "pf");

  nl = 5;
  TLegend * legB2blast = new TLegend(0.1, 0.65-nl*0.04, 0.45, 0.65, "Blast-Wave (#pi,K,p) + GSI-Heid. (T = 156 MeV)");
  legB2blast->SetFillStyle(0);
  legB2blast->SetTextSize(0.035);
  legB2blast->SetBorderSize(0);
  legB2blast->AddEntry(gBlastB2vsR_PbPb502TeV[ip], "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV", "l");
  legB2blast->AddEntry(gBlastB2vsR_PbPb276TeV[ip], "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "l");
  legB2blast->AddEntry(gBlastB2vsR_pPb502TeV[ip], "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV", "l");
  legB2blast->AddEntry(gBlastB2vsR_pp7TeV[ip], "pp #sqrt{#it{s}} = 7 TeV", "l");
  
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
  pavept->Draw();
  if (!plotLinX) gPad->SetLogx();
  hframe->Draw();
  hB2_coalescence->Draw("l");
  hB2_coalescence_pointlike->Draw("lsame");
  
  gBlastB2vsR_PbPb276TeV[ip]->Draw("samel");
  gBlastB2vsR_PbPb502TeV[ip]->Draw("samel");
  gBlastB2vsR_pPb502TeV[ip]->Draw("samel");
  gBlastB2vsR_pp7TeV[ip]->Draw("samel");
  gB2vsR_pp7TeV_sys[ip]->Draw("p3");
  gB2vsR_pp7TeV[ip]->Draw("samep");
  gB2vsR_pPb5TeV_sys[ip]->Draw("p3");
  gB2vsR_pPb5TeV[ip]->Draw("samep");
  gB2vsR_PbPb5TeV_sys[ip]->Draw("p3");
  gB2vsR_PbPb5TeV[ip]->Draw("samep");
  gB2vsR_PbPb276TeV_sys[ip]->Draw("p3");
  gB2vsR_PbPb276TeV[ip]->Draw("samep");
  gB2vsR_pp13TeV_sys[ip]->Draw("p3");
  gB2vsR_pp13TeV[ip]->Draw("samep");
  pavept->Draw();
  cb2opta->cd(2);
  legB2data->Draw();
  legB2coal->Draw();
  legB2blast->Draw();
  paveparam->Draw();  
  //--------------------
  //Alternative plotting with legends on the side -- B3
  //--------------------
  nl = 8;
  TLegend * legB3data = new TLegend(0.1, 0.95-nl*0.04, 0.45, 0.95, "ALICE");
  legB3data->SetFillStyle(0);
  legB3data->SetTextSize(0.04);
  legB3data->SetBorderSize(0);
  legB3data->AddEntry(gB3vsR_PbPb5TeV_sys[ip], "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, prelim.", "pf");
  legB3data->AddEntry(gB3vsR_PbPb276TeV_sys[ip], "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV [PRC 93, 0249717 (2016)]", "pf");
  legB3data->AddEntry(gB3vsR_pPb5TeV_sys[ip], "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, 21.09.2018", "pf");
  legB3data->AddEntry(gB3vsR_pp7TeV_sys[ip], "pp #sqrt{#it{s}} = 7 TeV [arXiv:1709.08522]", "pf");
  //legB3data->AddEntry(gB3LambdavsR_PbPb276TeV_sys[ip], "#it{B}_{3,#Lambda}, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV [PLB 754, 360-372 (2016)]", "pf");

  nl = 3;
  TLegend * legB3blast = new TLegend(0.1, 0.55-nl*0.04, 0.45, 0.55, "Blast-Wave (#pi,K,p) + GSI-Heid. (T = 156 MeV)");
  legB3blast->SetFillStyle(0);
  legB3blast->SetTextSize(0.04);
  legB3blast->SetBorderSize(0);
  legB3blast->AddEntry(gBlastB3vsR_PbPb276TeV[ip], "#it{B}_{3}, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "l");
  //legB3blast->AddEntry(gBlastB3LambdavsR_PbPb276TeV[ip], "#it{B}_{3,#Lambda}, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "l");

  nl = 4;
  TLegend * legB3coal = new TLegend(0.1, 0.4-nl*0.04, 0.45, 0.4, "Coalescence");
  legB3coal->SetFillStyle(0);
  legB3coal->SetTextSize(0.04);
  legB3coal->SetBorderSize(0);
  legB3coal->AddEntry(hB3_coalescence, "#it{B}_{3} coalesc., #it{r}(^{3}He) = 2.48 fm", "l");
  //legB3coal->AddEntry(hB3_coalescence_pointlike, "#it{B}_{3} coalesc., #it{r}(^{3}He) = 0 (point-like)", "l");
  //legB3coal->AddEntry(hB3L_coalescence, "#it{B}_{3,#Lambda} coalesc., #it{r}(^{3}_{#Lambda}H) = 6.8 fm", "l");

  TCanvas * cb3opta = new TCanvas("cb3opta", "Frascati plot B3", 1000, 600);
  cb3opta->SetBottomMargin(0.02);
  cb3opta->SetTopMargin(0.05);
  cb3opta->SetLeftMargin(0.15);
  cb3opta->SetRightMargin(0.02);
  
  cb3opta->Divide(2,1);
  cb3opta->cd(1);
  paveptB3->Draw();

  gPad->SetLogy();
  if (!plotLinX)  gPad->SetLogx();
  hframe3->Draw();
  paveptB3->Draw();
  hB3_coalescence->Draw("l");
  //hB3_coalescence_pointlike->Draw("lsame");
  //hB3L_coalescence->Draw("lsame");
  gBlastB3vsR_PbPb276TeV[ip]->Draw("samel");
  //gBlastB3LambdavsR_PbPb276TeV->Draw("samel");
  gB3vsR_PbPb5TeV_sys[ip]->Draw("p3");
  gB3vsR_PbPb5TeV[ip]->Draw("samep");
  gB3vsR_PbPb276TeV_sys[ip]->Draw("samep3");
  gB3vsR_PbPb276TeV[ip]->Draw("samep");
  gB3vsR_pp7TeV_sys[ip]->Draw("samep2");
  gB3vsR_pp7TeV[ip]->Draw("samep");
  gB3vsR_pPb5TeV_sys[ip]->Draw("p3");
  gB3vsR_pPb5TeV[ip]->Draw("samep");
  // gB3LambdavsR_PbPb276TeV_sys->Draw("samep2");
  // gB3LambdavsR_PbPb276TeV->Draw("samep");

  cb3opta->cd(2);
  legB3data->Draw();
  legB3coal->Draw();
  legB3blast->Draw();
  paveparam->Draw();  

  TCanvas * cb3optaLambda = new TCanvas("cb3optaLambda", "Frascati plot B3 Lambda", 1000, 600);
  cb3optaLambda->SetBottomMargin(0.02);
  cb3optaLambda->SetTopMargin(0.05);
  cb3optaLambda->SetLeftMargin(0.15);
  cb3optaLambda->SetRightMargin(0.02);
  
  cb3optaLambda->Divide(2,1);
  cb3optaLambda->cd(1);
  gPad->SetLogy();
  if (!plotLinX)  gPad->SetLogx();
  hframe3->DrawCopy();
  hB3L_coalescence->Draw("lsame");
  gBlastB3LambdavsR_PbPb276TeV[ip]->Draw("samel");
  gB3LambdavsR_PbPb276TeV_sys[ip]->Draw("samep2");
  gB3LambdavsR_PbPb276TeV[ip]->Draw("samep");


  cb3optaLambda->cd(2);
  legB3data->Draw();
  legB3coal->Draw();
  legB3blast->Draw();
  paveparam->Draw();  
  
  
  return 0;  
}



void MakePaperFigure1(Bool_t plotLinX, Double_t pToA,
		      TF1 *Cd_coalescence, TF1* Cd_coalescence_pointlike, TF1* Cd_coalescence_radius1third,  TF1* Cd_coalescence_largeradius,
		      TGraphErrors * hB2_coalescence, TGraphErrors * hB2_coalescence_pointlike, TGraphErrors * hB2_coalescence_radius1third, TGraphErrors * hB2_coalescence_largeradius, Bool_t plotvert) {
  //
  // Create the (pure theory) figure which plots <C_d> and B2 vs R
  // for different radii (PLOT COALESCENCE ONLY)
  //
  Float_t leftmargin = 0.15;
  Float_t rightmargin = 0.05;

  TH2D * frame_cd = new TH2D("frame_cd", "#LTC_{d}#GT vs radius; #it{R} (fm); #LT#it{C}_{d}#GT", 1000, 0.01, 6.0, 120, -0.03, 1.2);
  frame_cd->GetXaxis()->SetTitleSize(0.07);
  frame_cd->GetYaxis()->SetTitleSize(0.07);
  frame_cd->GetXaxis()->SetLabelSize(0.07);
  frame_cd->GetYaxis()->SetLabelSize(0.07);
  frame_cd->GetYaxis()->SetTitleOffset(0.9);
  frame_cd->GetXaxis()->SetRangeUser(0.01, 6.);

  TH2D * frame_coal = new TH2D("frame_coal", "B_{2} vs radius; #it{R} (fm); #it{B}_{2} (GeV^{2}/#it{c}^{3})", 1000, 0.01, 6.0, 2e5, 5.e-5, 0.3);
  frame_coal->GetXaxis()->SetTitleSize(0.07);
  frame_coal->GetYaxis()->SetTitleSize(0.07);
  frame_coal->GetYaxis()->SetTitleOffset(0.9);
  frame_coal->GetXaxis()->SetLabelSize(0.07);
  frame_coal->GetYaxis()->SetLabelSize(0.07);
  frame_coal->GetXaxis()->SetRangeUser(0.01, 6.);
 
  TLine * lInflectionPoint = new TLine(1.3, -0.03, 1.3, 1.2);
  lInflectionPoint->SetLineStyle(3);
  lInflectionPoint->SetLineWidth(2);
  
  TLegend * legB2_coal;
  legB2_coal = new TLegend(0.52, 0.93-5*0.06, 0.85, 0.93, Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToA));
  legB2_coal->SetFillStyle(0);
  legB2_coal->SetTextSize(0.05);
  legB2_coal->SetBorderSize(0);
  legB2_coal->AddEntry(hB2_coalescence_pointlike, "#it{B}_{2} coalesc., #it{r_{d}} = 0 (point-like)", "l");
  legB2_coal->AddEntry(hB2_coalescence_radius1third, "#it{B}_{2} coalesc., #it{r_{d}} = 0.3 fm", "l");
  legB2_coal->AddEntry(hB2_coalescence, "#it{B}_{2} coalesc., #it{r_{d}} = 3.2 fm", "l");
  legB2_coal->AddEntry(hB2_coalescence_largeradius, "#it{B}_{2} coalesc., #it{r_{d}} = 10 fm", "l");
  
  TCanvas * coalcanv = new TCanvas("coalcanv", "coalescence", 800, 1800);
  //top pad
  TPad * pad1 = new TPad("pad1","pad1", 0.01, 0.51, 0.99, 0.99);
  pad1->SetFillColor(0);
  pad1->SetBorderMode(0);
  pad1->SetBorderSize(0);
  pad1->SetMargin(leftmargin, rightmargin, 0.0, 0.05);
  pad1->SetTicky();
  pad1->SetTickx();

  coalcanv->cd();
  pad1->Draw();
  pad1->cd();
  frame_cd->Draw();
  Cd_coalescence->Draw("same");
  Cd_coalescence_pointlike->Draw("same");
  Cd_coalescence_radius1third->Draw("same");
  Cd_coalescence_largeradius->Draw("same");
  lInflectionPoint->Draw();

  //bottom pad
  TPad * pad2 = new TPad("pad2","pad2", 0.01, 0.01, 0.99, 0.51); 
  pad2->SetFillColor(0);
  pad2->SetBorderMode(0);
  pad2->SetBorderSize(0);
  pad2->SetMargin(leftmargin, rightmargin, 0.15, 0.0);
  pad2->SetLogy();
  pad2->SetTicky();
  pad2->SetTickx();
  
  coalcanv->cd();
  pad2->Draw();
  pad2->cd();
  frame_coal->Draw();
  hB2_coalescence->Draw("same");
  hB2_coalescence_pointlike->Draw("same");
  hB2_coalescence_radius1third->Draw("same");
  hB2_coalescence_largeradius->Draw("same");
  legB2_coal->Draw();
  
  if (plotvert) {
    coalcanv->SaveAs("../Paper/theory_coalescence_Cd_B2_vert.eps");
    coalcanv->SaveAs("../Paper/theory_coalescence_Cd_B2_vert.png");
  } else {
    coalcanv->SaveAs("../Paper/theory_coalescence_Cd_B2.eps");
    coalcanv->SaveAs("../Paper/theory_coalescence_Cd_B2.png");
  }
  return;
}


void MakePaperFigure3(Bool_t plotLinX, Double_t pToA, Double_t pToAb3,
		      TGraphErrors * hB2_coalescence, TGraphErrors * hB3_coalescence, 
		      TGraphErrors ** gB2vsR_PbPb276TeV_sys,  TGraphErrors ** gB2vsR_pp7TeVINELg0_sys, TGraphErrors **gB2vsR_PbPb276TeV, TGraphErrors **  gB2vsR_pp7TeVINELg0,
		      TGraphErrors ** gB3vsR_PbPb276TeV_sys,  TGraphAsymmErrors ** gB3vsR_pp7TeV_sys, TGraphErrors **gB3vsR_PbPb276TeV, TGraphAsymmErrors **  gB3vsR_pp7TeV) {
  //
  // Make the 4-panel figure comparing B2 and B3 with data for the 
  // tuned and non-tuned parameterisation.
  // 

  TH2D * hframe = new TH2D("hframe", "B_{2} vs radius; #it{R} (fm); #it{B}_{2} (GeV^{2}/#it{c}^{3})", 1000, 0.01, 6.0, 2000, 1.e-4, 0.1);
  hframe->GetXaxis()->SetTitleSize(0.06);
  hframe->GetYaxis()->SetTitleSize(0.06);
  hframe->GetXaxis()->SetTitleOffset(0.8);
  hframe->GetXaxis()->SetLabelSize(0.05);
  hframe->GetYaxis()->SetLabelSize(0.05);
  if (plotLinX) hframe->GetXaxis()->SetRangeUser(0.01, 8.5);
  else  hframe->GetXaxis()->SetRangeUser(0.1, 10.5);

  TH2D * hframe3 = new TH2D("hframe3", "B_{3} vs radius; #it{R} (fm); #it{B}_{3} (GeV^{4}/#it{c}^{6})", 1000, 0.01, 6.0, 2000, 1.e-9, 1.e-1);
  hframe3->GetXaxis()->SetTitleSize(0.06);
  hframe3->GetYaxis()->SetTitleSize(0.06);
  hframe3->GetXaxis()->SetTitleOffset(0.8);
  hframe3->GetXaxis()->SetLabelSize(0.05);
  hframe3->GetYaxis()->SetLabelSize(0.05);
  if (plotLinX) hframe3->GetXaxis()->SetRangeUser(0.01, 8.5);
  else  hframe3->GetXaxis()->SetRangeUser(0.1, 10.5);

  //Define pT/A labels only once
  TPaveText * pavept = new TPaveText(0.17, 0.17, 0.7, 0.23, "NDC");
  pavept->SetFillStyle(0);
  pavept->SetTextFont(42);
  pavept->SetBorderSize(0);
  pavept->SetTextSize(0.05);
  pavept->SetTextAlign(12);
  pavept->AddText(Form("#it{B}_{2}: #it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToA));

  //TPaveText * paveptB3 = new TPaveText(0.55, 0.62, 0.95, 0.67, "NDC");
  TPaveText * paveptB3 = new TPaveText(0.17, 0.17, 0.7, 0.23, "NDC");
  paveptB3->SetFillStyle(0);
  paveptB3->SetTextFont(42);
  paveptB3->SetBorderSize(0);
  paveptB3->SetTextSize(0.05);
  paveptB3->SetTextAlign(12);
  paveptB3->AddText(Form("#it{B}_{3}: #it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb3));

  TLegend * legA = new TLegend(0.17, 0.13, 0.45, 0.32);
  legA->SetFillStyle(0);
  legA->SetTextFont(42);
  legA->SetTextSize(0.05);
  legA->SetBorderSize(0);
  legA->AddEntry(hB2_coalescence, "#it{B}_{2} coalesc., #it{r}(d) = 3.2 fm", "l");
  legA->AddEntry(gB2vsR_PbPb276TeV_sys[0], "ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "p");
  legA->AddEntry(gB2vsR_pp7TeVINELg0_sys[0], "ALICE, pp #sqrt{#it{s}} = 7 TeV (INEL>0)", "p");

  TLegend * legC = new TLegend(0.17, 0.13, 0.45, 0.32);
  legC->SetFillStyle(0);
  legC->SetTextFont(42);
  legC->SetTextSize(0.05);
  legC->SetBorderSize(0);
  legC->AddEntry(hB3_coalescence, "#it{B}_{3} coalesc., #it{r}(^{3}He) = 2.48 fm", "l");
  legC->AddEntry(gB2vsR_PbPb276TeV_sys[0], "ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "p");
  legC->AddEntry(gB2vsR_pp7TeVINELg0_sys[0], "ALICE, pp #sqrt{#it{s}} = 7 TeV (INEL>0)", "p");
  
  TPaveText * paveA = new TPaveText(0.17, 0.82, 0.80, 0.95, "NDC");
  paveA->SetFillStyle(0);
  paveA->SetTextFont(42);
  paveA->SetBorderSize(0);
  paveA->SetTextSize(0.05);
  paveA->SetTextAlign(12);
  paveA->InsertText("Linear fit to HBT radii: ");
  paveA->InsertText("#it{R} = #it{a} #LT d#it{N}_{ch}/d#it{#eta}#GT^{1/3} + #it{b}, #it{a} = 0.339, #it{b} = 0.128");

  TPaveText * paveB = new TPaveText(0.17, 0.82, 0.80, 0.95, "NDC");
  paveB->SetFillStyle(0);
  paveB->SetTextFont(42);
  paveB->SetBorderSize(0);
  paveB->SetTextSize(0.05);
  paveB->SetTextAlign(12);
  paveB->InsertText("Fixed to #it{B}_{2} in 0-10% Pb-Pb:");
  //paveB->InsertText("#it{R} = #it{a} #LT d#it{N}_{ch}/d#it{#eta}#GT^{1/3} + #it{b}");
  paveB->InsertText("#it{R} = #it{a} #LT d#it{N}_{ch}/d#it{#eta}#GT^{1/3} + #it{b}, #it{a} = 0.473, #it{b} = 0");

  
  TCanvas * cr1 = new TCanvas("cr1", "compare radii parameterisations", 1000, 1000);
  cr1->SetBottomMargin(0.02);
  cr1->SetTopMargin(0.01);
  cr1->SetLeftMargin(0.12);
  cr1->SetRightMargin(0.02);
  cr1->Divide(2,2);
  
  cr1->cd(1);
  gPad->SetLogy();
  hframe->Draw();
  gB2vsR_pp7TeVINELg0_sys[0]->Draw("samep2");
  gB2vsR_pp7TeVINELg0[0]->Draw("samepz");
  gB2vsR_PbPb276TeV_sys[0]->Draw("samep3");
  gB2vsR_PbPb276TeV[0]->Draw("samepz");
  hB2_coalescence->Draw("l");
  pavept->Draw();
  paveA->Draw();
  
  cr1->cd(2);
  gPad->SetLogy();
  hframe->Draw();
  gB2vsR_pp7TeVINELg0_sys[1]->Draw("samep2");
  gB2vsR_pp7TeVINELg0[1]->Draw("samepz");
  gB2vsR_PbPb276TeV_sys[1]->Draw("samep3");
  gB2vsR_PbPb276TeV[1]->Draw("samepz");
  hB2_coalescence->Draw("l");
  paveB->Draw();
  legA->Draw();

  cr1->cd(3);
  gPad->SetLogy();
  hframe3->Draw();
  gB3vsR_PbPb276TeV_sys[0]->Draw("samep3");
  gB3vsR_PbPb276TeV[0]->Draw("samepz");
  gB3vsR_pp7TeV_sys[0]->Draw("samep2");
  gB3vsR_pp7TeV[0]->Draw("samepz");
  hB3_coalescence->Draw("l");
  paveptB3->Draw();
  paveA->Draw();

  cr1->cd(4);
  gPad->SetLogy();
  hframe3->Draw();
  gB3vsR_PbPb276TeV_sys[1]->Draw("samep3");
  gB3vsR_PbPb276TeV[1]->Draw("samepz");
  gB3vsR_pp7TeV_sys[1]->Draw("samep2");
  gB3vsR_pp7TeV[1]->Draw("samepz");
  hB3_coalescence->Draw("l");
  paveB->Draw();
  legC->Draw();

  // cr1->SaveAs("Paper/radiiParamCompareData.eps");
  // cr1->SaveAs("Paper/radiiParamCompareData.png");
  

}

void MakePaperFigure4(Bool_t plotLinX, Double_t pToA, Double_t pToAb3,  Double_t pToAb3Lambda, 
		      TGraphErrors * hB2_coalescence, TGraphErrors * hB3_coalescence, TGraphErrors* hB3L_coalescence, TGraphErrors* hB3L_coalescence_largeradius, 
		      TGraphErrors ** gB2vsR_PbPb276TeV_sys,  TGraphErrors ** gB2vsR_pp7TeVINELg0_sys, TGraphErrors **gB2vsR_PbPb276TeV, TGraphErrors **  gB2vsR_pp7TeVINELg0,
		      TGraphErrors ** gB3vsR_PbPb276TeV_sys,  TGraphAsymmErrors ** gB3vsR_pp7TeV_sys, TGraphErrors ** gB3vsR_PbPb276TeV, TGraphAsymmErrors **  gB3vsR_pp7TeV,
		      TGraphAsymmErrors** gBlastB2vsR_PbPb276TeV,  TGraphAsymmErrors** gBlastB3vsR_PbPb276TeV,
		      TGraphAsymmErrors** gB3LambdavsR_PbPb276TeV, TGraphAsymmErrors** gB3LambdavsR_PbPb276TeV_sys, TGraphAsymmErrors** gBlastB3LambdavsR_PbPb276TeV) {
  //
  // make figure 4 of the current paper
  //
  TH1D * hframe = new TH1D("hframeFig4", "B_{2} vs radius; #it{R} (fm); #it{B}_{2} (GeV^{2}/#it{c}^{3})", 1000, 0.01, 6.0);
  hframe->GetXaxis()->SetTitleSize(0.07);
  hframe->GetYaxis()->SetTitleSize(0.07);
  hframe->GetYaxis()->SetTitleOffset(1.);
  hframe->GetXaxis()->SetTitleOffset(0.8);
  hframe->GetXaxis()->SetLabelSize(0.06);
  hframe->GetYaxis()->SetLabelSize(0.06);
  if (plotLinX) hframe->GetXaxis()->SetRangeUser(0.01, 8.5);
  else  hframe->GetXaxis()->SetRangeUser(0.1, 10.5);

  TH1D * hframe3 = new TH1D("hframe3Fig4", "B_{3} vs radius; #it{R} (fm); #it{B}_{3} (GeV^{4}/#it{c}^{6})", 1000, 0.01, 6.0);
  hframe3->GetXaxis()->SetTitleSize(0.07);
  hframe3->GetYaxis()->SetTitleSize(0.07);
  hframe3->GetYaxis()->SetTitleOffset(1.);
  hframe3->GetXaxis()->SetTitleOffset(0.8);
  hframe3->GetXaxis()->SetLabelSize(0.06);
  hframe3->GetYaxis()->SetLabelSize(0.06);
  if (plotLinX) hframe3->GetXaxis()->SetRangeUser(0.01, 8.5);
  else  hframe3->GetXaxis()->SetRangeUser(0.1, 10.5);

  TH1D * hframe3L = new TH1D("hframe3LFig4", "B_{3,#Lambda} vs radius; #it{R} (fm); #it{B}_{3,#Lambda} (GeV^{4}/#it{c}^{6})", 1000, 0.01, 6.0);
  hframe3L->GetXaxis()->SetTitleSize(0.06);
  hframe3L->GetYaxis()->SetTitleSize(0.06);
  hframe3L->GetYaxis()->SetTitleOffset(1.2);
  hframe3L->GetXaxis()->SetTitleOffset(0.8);
  hframe3L->GetXaxis()->SetLabelSize(0.05);
  hframe3L->GetYaxis()->SetLabelSize(0.05);
  if (plotLinX) hframe3L->GetXaxis()->SetRangeUser(0.01, 8.5);
  else  hframe3L->GetXaxis()->SetRangeUser(0.1, 10.5);

  //define particle label
  TPaveText * paveLab2 = new TPaveText(0.22, 0.82, 0.3, 0.92, "NDC");
  paveLab2->SetFillStyle(0);
  paveLab2->SetTextFont(42);
  paveLab2->SetBorderSize(0);
  paveLab2->SetTextSize(0.12);
  paveLab2->SetTextAlign(12);
  paveLab2->AddText("#bf{d}");

  TPaveText * paveLab3 = new TPaveText(0.19, 0.85, 0.3, 0.95, "NDC");
  paveLab3->SetFillStyle(0);
  paveLab3->SetTextFont(42);
  paveLab3->SetBorderSize(0);
  paveLab3->SetTextSize(0.11);
  paveLab3->SetTextAlign(12);
  paveLab3->AddText("#bf{^{3}He}");

  TPaveText * paveLab3L = new TPaveText(0.17, 0.85, 0.3, 0.95, "NDC");
  paveLab3L->SetFillStyle(0);
  paveLab3L->SetTextFont(42);
  paveLab3L->SetBorderSize(0);
  paveLab3L->SetTextSize(0.1);
  paveLab3L->SetTextAlign(12);
  paveLab3L->AddText("#bf{ ^{3}_{#Lambda}H}");
  
  //Define pT/A labels only once
  TPaveText * pavept = new TPaveText(0.17, 0.87, 0.7, 0.92, "NDC");
  pavept->SetFillStyle(0);
  pavept->SetTextFont(42);
  pavept->SetBorderSize(0);
  pavept->SetTextSize(0.05);
  pavept->SetTextAlign(12);
  pavept->AddText(Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToA));

  //TPaveText * paveptB3 = new TPaveText(0.55, 0.62, 0.95, 0.67, "NDC");
  TPaveText * paveptB3 = new TPaveText(0.17, 0.87, 0.7, 0.92, "NDC");
  paveptB3->SetFillStyle(0);
  paveptB3->SetTextFont(42);
  paveptB3->SetBorderSize(0);
  paveptB3->SetTextSize(0.05);
  paveptB3->SetTextAlign(12);
  paveptB3->AddText(Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb3));

  TPaveText * paveptB3L = new TPaveText(0.17, 0.87, 0.7, 0.92, "NDC");
  paveptB3L->SetFillStyle(0);
  paveptB3L->SetBorderSize(0);
  paveptB3L->SetTextFont(42);
  paveptB3L->SetTextSize(0.05);
  paveptB3L->SetTextAlign(12);
  paveptB3L->AddText(Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb3Lambda));
  
  /*TCanvas * cr4 = new TCanvas("cr4", "compare coalesence with thermal", 1600, 600);
  cr4->SetBottomMargin(0.02);
  cr4->SetTopMargin(0.01);
  cr4->SetLeftMargin(0.12);
  cr4->SetRightMargin(0.01);
  cr4->Divide(3,1);
  */
   
  TCanvas * cr4 = new TCanvas("cr4", "compare thermal with coalescence", 700, 1400);
  Float_t lowx[2] = {0.001, 0.99};
  Float_t lowy[4] = {0.996, 0.68, 0.38, 0.01};
  gStyle->SetPadBottomMargin(0.12);

  TPad * pad[3][1];
  for (int c = 0; c<1; c++) {
    for (int r = 0; r<3; r++) {
      pad[c][r] = new TPad(Form("pad%i%i",c,r),"pad", lowx[c], lowy[r], lowx[c+1], lowy[r+1]);
    }
  }

  Float_t leftmargin = 0.15;
  Float_t rightmargin = 0.02;
  //top pad
  pad[0][0]->SetFillColor(0);
  pad[0][0]->SetBorderMode(0);
  pad[0][0]->SetBorderSize(0);
  pad[0][0]->SetMargin(leftmargin, rightmargin, 0.001, 0.05);
  pad[0][0]->SetLogy();
  pad[0][0]->SetTicky();
  pad[0][0]->SetTickx();

  //middle pad
  pad[1][0]->SetFillColor(0);
  pad[1][0]->SetBorderMode(0);
  pad[1][0]->SetBorderSize(0);
  pad[1][0]->SetMargin(leftmargin, rightmargin, 0.001, 0.001);
  pad[1][0]->SetTicky();
  pad[1][0]->SetTickx();

  //bottom pad
  pad[2][0]->SetFillColor(0);
  pad[2][0]->SetBorderMode(0);
  pad[2][0]->SetBorderSize(0);
  pad[2][0]->SetMargin(leftmargin, rightmargin, 0.2, 0.001);
  pad[2][0]->SetTicky();
  pad[2][0]->SetTickx();
  
  cr4->cd();
  pad[0][0]->Draw();
  pad[0][0]->cd();
  pad[0][0]->SetLogy();
  hframe->GetYaxis()->SetRangeUser(3.E-5, 8E-2);
  hframe->Draw();
  gBlastB2vsR_PbPb276TeV[1]->Draw("samel");
  gB2vsR_pp7TeVINELg0_sys[1]->Draw("samep2");
  gB2vsR_pp7TeVINELg0[1]->Draw("samepz");
  gB2vsR_PbPb276TeV_sys[1]->Draw("samep3");
  gB2vsR_PbPb276TeV[1]->Draw("samepz");
  //
  TLegend * legB2 = new TLegend(0.2, 0.05, 0.6, 0.45, Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToA));
  legB2->SetFillStyle(0);
  legB2->SetTextSize(0.045);
  legB2->SetBorderSize(0);
  legB2->AddEntry(gB2vsR_pp7TeVINELg0_sys[1], "ALICE, pp #sqrt{#it{s}} = 7 TeV", "pf"); //[PRC 97, 024615 (2018)]
  legB2->AddEntry(gB2vsR_PbPb276TeV_sys[1], "ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "pf"); //[PRC 93, 0249717 (2016)]
  legB2->AddEntry(hB2_coalescence, "Coal., #it{r}(d) = 3.2 fm", "l");
  legB2->AddEntry(gBlastB2vsR_PbPb276TeV[1], "GSI-Heidelberg (#it{T_{chem}} = 156 MeV) + ", "l");
  legB2->AddEntry(gBlastB2vsR_PbPb276TeV[1], "blast-wave (#piKp, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV)", "");

  hB2_coalescence->Draw("l");
  paveLab2->Draw();
  legB2->Draw();
  
  cr4->cd();
  pad[1][0]->Draw();
  pad[1][0]->cd();
  pad[1][0]->SetLogy();
  hframe3->GetYaxis()->SetRangeUser(6.E-10, 8E-2);
  hframe3->Draw();
  hB3_coalescence->Draw("l");
  gBlastB3vsR_PbPb276TeV[1]->Draw("samel");
  gB3vsR_PbPb276TeV_sys[1]->Draw("samep3");
  gB3vsR_PbPb276TeV[1]->Draw("samepz");
  gB3vsR_pp7TeV_sys[1]->Draw("samep2");
  gB3vsR_pp7TeV[1]->Draw("samepz");
  //
  TLegend * legB3 = new TLegend(0.2, 0.05, 0.6, 0.2, Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb3));
  legB3->SetFillStyle(0);
  legB3->SetTextSize(0.055);
  legB3->SetBorderSize(0);
  legB3->AddEntry(hB3_coalescence, "Coal., #it{r}(^{3}He) = 2.48 fm", "l");
  // legB3->AddEntry(gB3vsR_PbPb276TeV_sys[1], "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV [PRC 93, 0249717 (2016)]", "pf");
  // legB3->AddEntry(gB3vsR_pp7TeV_sys[1], "pp #sqrt{#it{s}} = 7 TeV [arXiv:1709.08522]", "pf");
  // legB3->AddEntry(gBlastB3vsR_PbPb276TeV[1], "#it{B}_{3}, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV, BW + GSI (T = 156 MeV)", "l");
  //legB3->Draw();
  //
  //paveptB3
  legB3->Draw();
  paveLab3->Draw();

  cr4->cd();
  pad[2][0]->Draw();
  pad[2][0]->cd();
  pad[2][0]->SetLogy();
  hframe3L->GetYaxis()->SetRangeUser(6.E-10, 8E-2);
  hframe3L->Draw();
  hB3L_coalescence_largeradius->Draw("samel");
  hB3L_coalescence->Draw("l");
  gBlastB3LambdavsR_PbPb276TeV[1]->Draw("samel");
  gB3LambdavsR_PbPb276TeV_sys[1]->Draw("samep2");
  gB3LambdavsR_PbPb276TeV[1]->Draw("samep");
  //paveptB3L->Draw();
  paveLab3L->Draw();
  //
  TLegend * legB3Lambda = new TLegend(0.5, 0.78, 0.8, 0.97, Form("#it{p}_{T}/#it{A} = %1.0f GeV/#it{c}", pToAb3Lambda));
  legB3Lambda->SetFillStyle(0);
  legB3Lambda->SetTextSize(0.045);
  legB3Lambda->SetBorderSize(0);
  legB3Lambda->AddEntry(hB3L_coalescence, "Coal., #it{r}(^{3}_{#Lambda}H) = 6.8 fm", "l");
  legB3Lambda->AddEntry(hB3L_coalescence_largeradius, "Coal., #it{r} (^{3}_{#Lambda}H) = 14.1 fm", "l");

  // legB3Lambda->AddEntry(gB3LambdavsR_PbPb276TeV_sys[1], "#it{B}_{3,#Lambda}, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV [PLB 754, 360-372 (2016)]", "pf");
  // legB3Lambda->AddEntry(gBlastB3LambdavsR_PbPb276TeV[1], "#it{B}_{3,#Lambda}, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV, BW + GSI (T = 156 MeV)", "l");

  legB3Lambda->Draw();
  //

  // TLegend * masterLeg = new TLegend(0.1, 0.3, 0.5, 0.9, "");
  // masterLeg->SetFillStyle(0);
  // masterLeg->SetTextSize(0.05);
  // masterLeg->SetBorderSize(0);
  // masterLeg->AddEntry(gB2vsR_PbPb276TeV_sys[1], "ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "pf");
  // masterLeg->AddEntry(gB2vsR_pp7TeVINELg0_sys[1], "ALICE, pp #sqrt{#it{s}} = 7 TeV (INEL>0)", "pf");
  // masterLeg->AddEntry(gBlastB2vsR_PbPb276TeV[1], "BW + GSI-Heidelberg (#it{T}_{chem} = 156 MeV)", "l");
  // //masterLeg->AddEntry(gBlastB2vsR_PbPb276TeV[1], "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "");
  // masterLeg->AddEntry(hB2_coalescence, "#it{B}_{#it{A}} coalescence", "l");
  // masterLeg->AddEntry(hB2_coalescence, "#it{r} (d) = 3.2 fm", "");
  // masterLeg->AddEntry(hB3_coalescence, "#it{r} (^{3}He) = 2.48 fm", "");
  // masterLeg->AddEntry(hB3L_coalescence, "#it{r} (^{3}_{#Lambda}H) = 6.8 fm", "");
  // masterLeg->AddEntry(hB3L_coalescence_largeradius, "#it{r} (^{3}_{#Lambda}H) = 14.1 fm", "l");

  //cr4->cd(2);
  //masterLeg->Draw();
  
  // cr4->SaveAs("Paper/compareThermalAndCoalescence.eps");
  // cr4->SaveAs("Paper/compareThermalAndCoalescence.png");

}

//----------------------------------------------
void MakeUp(TGraphAsymmErrors* obj, Color_t color, Color_t Fill_Color, Int_t Fill_Style, Int_t Line_Style, Int_t Line_Width, Int_t Marker_Style, Float_t Marker_Size)
{
  if (!obj) return;
  obj->SetLineColor(color);
  obj->SetFillColor(Fill_Color);
  obj->SetFillStyle(Fill_Style);
  obj->SetFillColorAlpha(Fill_Color, 0.2);
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
  obj->SetFillColorAlpha(Fill_Color, 0.2);
  obj->SetLineStyle(Line_Style);
  obj->SetLineWidth(Line_Width);
  obj->SetMarkerColor(color);
  obj->SetMarkerStyle(Marker_Style);
  obj->SetMarkerSize(Marker_Size);
  return;
}

