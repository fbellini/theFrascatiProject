/*
   akalweit@cern.ch, fbellini@cern.ch
   17.01.2018 - The Frascati plot
*/
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TH2D.h"
#include "./generateBWpredictionsB2.C" //ADAPT ME

void convertMultiToRadius(TGraphErrors * graph = 0x0, Int_t paramSet = 0);
void convertMultiToRadius(TGraphAsymmErrors * graph = 0x0, Int_t paramSet = 0);

void getRadiusFromParameterisation(Double_t * multi = 0x0, Double_t * radius = 0x0, Int_t paramSet = 0);
TF1 * MakeB2TheoryGraphQMfactor(Double_t objSize = 3.2);
TGraphErrors * MakeB2TheoryGraphCoalescence(Double_t mT = 1.0, Double_t objSize = 3.2);
TGraphErrors * MakeB3TheoryGraphCoalescence(Double_t mT = 1.0, Double_t objSize = 2.48);

TGraphErrors * getB2_pp7TeV(Bool_t plotSys = 0, Double_t pToA = 0.75, Int_t paramSet = 0);
TGraphErrors * getB2_pp7TeVINELg0(Bool_t plotSys = 0, Double_t pToA = 0.75, Int_t paramSet = 0);
TGraphErrors * getB2_pPb5TeV(Bool_t plotSys = 0, Double_t pToA = 0.75, Int_t paramSet = 0);
TGraphErrors * getB2_PbPb5TeV(Bool_t plotSys = 0, Double_t pToA = 0.75, Int_t paramSet = 0);
TGraphErrors * getB2_PbPb276TeV(Bool_t plotSys = 0, Double_t pToA = 0.75, Int_t paramSet = 0);

TGraphErrors      * getB3_PbPb5TeV(Bool_t plotSys = 0, Double_t pToAb3 = 0.733, Int_t paramSet = 0);
TGraphErrors      * getB3_PbPb276TeV(Bool_t plotSys = 0, Double_t pToAb3 = 0.733, Int_t paramSet = 0);
TGraphAsymmErrors * getB3_pp7TeVINELg0(Bool_t plotSys = 0, Double_t pToAb3pp = 0.800, Int_t paramSet = 0);
TGraphAsymmErrors * getB3Lambda_PbPb276TeV(Bool_t plotSys = 0, Double_t pToAb3Lambda = 1., Int_t paramSet = 0);

TGraphAsymmErrors * getBlastB2_PbPb276TeV(Bool_t plotSys = 0, Double_t pToA = 0.75, Int_t paramSet = 0);
TGraphAsymmErrors * getBlastB2_PbPb502TeV(Bool_t plotSys = 0, Double_t pToA = 0.75, Int_t paramSet = 0);
TGraphAsymmErrors * getBlastB2_pPb502TeV(Bool_t plotSys = 0, Double_t pToA = 0.75, Int_t paramSet = 0);
TGraphAsymmErrors * getBlastB2_pp7TeV(Bool_t plotSys = 0, Double_t pToA = 0.75, Int_t paramSet = 0);

TGraphAsymmErrors * getBlastB3_PbPb276TeV(Bool_t plotSys = 0, Double_t pToAb3 = 0.733, Int_t paramSet = 0);
TGraphAsymmErrors * getBlastB3Lambda_PbPb276TeV(Bool_t plotSys = 0, Double_t pToAb3Lambda = 1.0, Int_t paramSet = 0);

void MakeUp(TGraphErrors* obj, Color_t color, Color_t Fill_Color, Int_t Fill_Style, Int_t Line_Style, Int_t Line_Width, Int_t Marker_Style, Float_t Marker_Size);
void MakeUp(TGraphAsymmErrors* obj, Color_t color, Color_t Fill_Color, Int_t Fill_Style, Int_t Line_Style, Int_t Line_Width, Int_t Marker_Style, Float_t Marker_Size);

void MakePaperFigure2(Bool_t plotLinX, Double_t pToA,
		      TF1 *Cd_coalescence, TF1* Cd_coalescence_pointlike, TF1* Cd_coalescence_radius1third,  TF1* Cd_coalescence_largeradius,
		      TGraphErrors * hB2_coalescence, TGraphErrors * hB2_coalescence_pointlike, TGraphErrors * hB2_coalescence_radius1third, TGraphErrors * hB2_coalescence_largeradius);

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

Int_t B2vsVolume(Bool_t plotLinX = 1, Double_t pToA = 0.75, Double_t pToAb3 = 0.733, Double_t pToAb3pp = 0.800, Double_t pToAb3Lambda = 1.,
		 Bool_t plotOnlyCoalescence = kFALSE, Bool_t plotPaperFigures = kTRUE)
{
  //
  // main function which generates the plots of the Frascati project
  //
  
  //--------------------
  //data
  //--------------------
  const Int_t nParamSet = 2;
  TGraphErrors* gB2vsR_pp7TeV[nParamSet];
  TGraphErrors* gB2vsR_pp7TeV_sys[nParamSet];
  TGraphErrors* gB2vsR_pp7TeVINELg0[nParamSet];
  TGraphErrors* gB2vsR_pp7TeVINELg0_sys[nParamSet];
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
  TGraphAsymmErrors* gB3vsR_pp7TeV[nParamSet];
  TGraphAsymmErrors* gB3vsR_pp7TeV_sys[nParamSet];
  TGraphAsymmErrors* gB3LambdavsR_PbPb276TeV[nParamSet];
  TGraphAsymmErrors* gB3LambdavsR_PbPb276TeV_sys[nParamSet];
  //
  TGraphAsymmErrors* gBlastB2vsR_PbPb276TeV[nParamSet];
  TGraphAsymmErrors* gBlastB2vsR_PbPb502TeV[nParamSet];
  TGraphAsymmErrors* gBlastB2vsR_pPb502TeV[nParamSet];
  TGraphAsymmErrors* gBlastB2vsR_pp7TeV[nParamSet];
  TGraphAsymmErrors* gBlastB3vsR_PbPb276TeV[nParamSet];
  TGraphAsymmErrors* gBlastB3LambdavsR_PbPb276TeV[nParamSet];

  for (Int_t ip = 0; ip < nParamSet; ip++){
    gB2vsR_pp7TeVINELg0[ip] = (TGraphErrors *) getB2_pp7TeVINELg0(kFALSE, pToA, ip);
    gB2vsR_pp7TeVINELg0_sys[ip] = (TGraphErrors *) getB2_pp7TeVINELg0(kTRUE, pToA, ip);

    gB2vsR_pp7TeV[ip] = (TGraphErrors *) getB2_pp7TeVINELg0(kFALSE, pToA, ip);
    gB2vsR_pp7TeV_sys[ip] = (TGraphErrors *) getB2_pp7TeVINELg0(kTRUE, pToA, ip);
      
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
    gBlastB3LambdavsR_PbPb276TeV[ip] = (TGraphAsymmErrors *)  getBlastB3Lambda_PbPb276TeV(kFALSE, pToAb3, ip);
  }
  
  //--------------------
  //theory - coalescence
  //--------------------
  //mT is the mass of the particle relative to which the HBT radius is calculated.
  //We use now the mass of the proton and the pT per nucleon
  Double_t mT = TMath::Sqrt(pToA * pToA + 0.938 * 0.938);
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
  Cd_coalescence_radius1third->SetLineWidth(3);
  Cd_coalescence_radius1third->SetLineStyle(7);
  Cd_coalescence_radius1third->SetLineColor(kRed);

  TF1 * Cd_coalescence_largeradius = (TF1*) MakeB2TheoryGraphQMfactor(10.0);
  Cd_coalescence_largeradius->SetLineWidth(3);
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
  hB2_coalescence_radius1third->SetLineWidth(3);
  hB2_coalescence_radius1third->SetLineStyle(7);
  hB2_coalescence_radius1third->SetLineColor(kRed);

  TGraphErrors* hB2_coalescence_largeradius = (TGraphErrors*) MakeB2TheoryGraphCoalescence(mT, 10.0);
  hB2_coalescence_largeradius->SetMarkerStyle(24);
  hB2_coalescence_largeradius->SetMarkerSize(0.4);
  hB2_coalescence_largeradius->SetLineWidth(3);
  hB2_coalescence_largeradius->SetLineStyle(9);
  hB2_coalescence_largeradius->SetLineColor(kGreen+1);

  TGraphErrors* hB3_coalescence = (TGraphErrors*) MakeB3TheoryGraphCoalescence(mT);
  hB3_coalescence->SetMarkerStyle(20);
  hB3_coalescence->SetMarkerColor(kBlack);
  hB3_coalescence->SetLineColor(kBlack);
  hB3_coalescence->SetLineWidth(3);

  TGraphErrors* hB3_coalescence_pointlike = (TGraphErrors*) MakeB3TheoryGraphCoalescence(mT, 0.0);
  hB3_coalescence_pointlike->SetMarkerStyle(24);
  hB3_coalescence_pointlike->SetMarkerColor(kGray);
  hB3_coalescence_pointlike->SetLineColor(kGray);
  hB3_coalescence_pointlike->SetMarkerSize(0.4);
  hB3_coalescence_pointlike->SetLineWidth(3);

  TGraphErrors* hB3L_coalescence = (TGraphErrors*) MakeB3TheoryGraphCoalescence(mT, 6.8);
  hB3L_coalescence->SetMarkerStyle(20);
  hB3L_coalescence->SetMarkerColor(kBlack);
  hB3L_coalescence->SetLineColor(kBlack);
  hB3L_coalescence->SetLineWidth(3);
  hB3L_coalescence->SetLineStyle(1);

  TGraphErrors* hB3L_coalescence_largeradius = (TGraphErrors*) MakeB3TheoryGraphCoalescence(mT, 14.1);
  hB3L_coalescence_largeradius->SetMarkerStyle(20);
  hB3L_coalescence_largeradius->SetMarkerColor(kBlack);
  hB3L_coalescence_largeradius->SetLineColor(kBlack);
  hB3L_coalescence_largeradius->SetLineWidth(3);
  hB3L_coalescence_largeradius->SetLineStyle(7);
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.02); 



  // ------------------
  // make the figure with coalescence only
  // ------------------
  MakePaperFigure2(plotLinX, pToA,
		   Cd_coalescence, Cd_coalescence_pointlike, Cd_coalescence_radius1third, Cd_coalescence_largeradius,
		   hB2_coalescence, hB2_coalescence_pointlike, hB2_coalescence_radius1third, hB2_coalescence_largeradius);
  if (plotOnlyCoalescence) return 0;

  //  if (plotPaperFigures) return 0;

  //------------------------------
  // PLOT FRASCATI PLOT(S)
  //------------------------------
  //make up options
  Int_t Fill_Style = 1001;
  Int_t Line_Style = 1;
  Int_t Line_Style_Blast = 2;
  Int_t Line_Width = 1;
  Int_t Line_Width_Blast = 5;
  Float_t Marker_Size = 1.3;

  enum EPlotEntries { kPP7, kPPB502, kPBPB276, kPBPB502,
		      kPP7blast, kPPB502blast, kPBPB276blast, kPBPB502blast,
		      kB3_PP7, kB3_PPB502, kB3_PBPB276, kB3_PBPB502, kB3L_PBPB276,
		      kB3_PP7blast, kB3_PPB502blast, kB3_PBPB276blast, kB3_PBPB502blast, kB3L_PBPB276blast};
  
  Color_t color[]      = {kGreen+2, kBlue+2, kRed, kRed+2,
			  kGreen+2, kBlue+2, kBlue+1, kRed+2,
			  kGreen+2, kSpring+4, kRed, kRed+2, kRed,
			  kBlue-5, kBlue-7, kBlue+1, kRed+2, kBlue+1};
  
  Int_t Marker_Style[] = { 21, 22, 20, 23,
			   21, 22, 20, 23,
			   21, 22, 20, 23, 33,
			   21, 22, 20, 23, 33};
  
  for (Int_t ip = 0; ip < nParamSet; ip++) {
    MakeUp(gB2vsR_pp7TeV_sys[ip], color[EPlotEntries::kPP7], color[EPlotEntries::kPP7], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kPP7], Marker_Size);
    MakeUp(gB2vsR_pp7TeV[ip]    , color[EPlotEntries::kPP7], color[EPlotEntries::kPP7], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kPP7], Marker_Size);
    MakeUp(gB2vsR_pp7TeVINELg0_sys[ip], color[EPlotEntries::kPP7], color[EPlotEntries::kPP7], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kPP7], Marker_Size);
    MakeUp(gB2vsR_pp7TeVINELg0[ip]    , color[EPlotEntries::kPP7], color[EPlotEntries::kPP7], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kPP7], Marker_Size);
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
    MakeUp(gB3LambdavsR_PbPb276TeV_sys[ip], color[EPlotEntries::kB3L_PBPB276], color[EPlotEntries::kB3L_PBPB276], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kB3_PBPB276], Marker_Size);
    MakeUp(gB3LambdavsR_PbPb276TeV[ip]    , color[EPlotEntries::kB3L_PBPB276], color[EPlotEntries::kB3L_PBPB276], Fill_Style, Line_Style, Line_Width, Marker_Style[EPlotEntries::kB3L_PBPB276], Marker_Size);
    //BLAST B3
    MakeUp(gBlastB3vsR_PbPb276TeV[ip]    , color[EPlotEntries::kB3_PBPB276blast], color[EPlotEntries::kB3_PBPB276blast], Fill_Style, Line_Style_Blast, Line_Width_Blast, Marker_Style[EPlotEntries::kB3_PBPB276blast], Marker_Size);
    MakeUp(gBlastB3LambdavsR_PbPb276TeV[ip] , color[EPlotEntries::kB3L_PBPB276blast], color[EPlotEntries::kB3L_PBPB276blast], Fill_Style, Line_Style_Blast, Line_Width_Blast, Marker_Style[EPlotEntries::kB3L_PBPB276blast], Marker_Size);
  }


  

  MakePaperFigure3(plotLinX, pToA, pToAb3,
		   hB2_coalescence, hB3_coalescence, 
		   gB2vsR_PbPb276TeV_sys,  gB2vsR_pp7TeVINELg0_sys, gB2vsR_PbPb276TeV,  gB2vsR_pp7TeVINELg0,
		   gB3vsR_PbPb276TeV_sys,  gB3vsR_pp7TeV_sys, gB3vsR_PbPb276TeV,  gB3vsR_pp7TeV);

  MakePaperFigure4(plotLinX, pToA, pToAb3, pToAb3Lambda,
		   hB2_coalescence, hB3_coalescence, hB3L_coalescence, hB3L_coalescence_largeradius,
		   gB2vsR_PbPb276TeV_sys,  gB2vsR_pp7TeVINELg0_sys, gB2vsR_PbPb276TeV,  gB2vsR_pp7TeVINELg0,
		   gB3vsR_PbPb276TeV_sys,  gB3vsR_pp7TeV_sys, gB3vsR_PbPb276TeV,  gB3vsR_pp7TeV,
		   gBlastB2vsR_PbPb276TeV, gBlastB3vsR_PbPb276TeV,
		   gB3LambdavsR_PbPb276TeV, gB3LambdavsR_PbPb276TeV_sys, gBlastB3LambdavsR_PbPb276TeV);

  if (plotPaperFigures) return 0;

  //---------------------------------------
  // PLOT FRASCATI PLOTS FOR SLIDES
  //---------------------------------------  

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

  TPaveText * paveptB3L = new TPaveText(0.55, 0.55, 0.95, 0.6, "NDC");
  paveptB3L->SetFillStyle(0);
  paveptB3L->SetBorderSize(0);
  paveptB3L->SetTextFont(42);
  paveptB3L->SetTextSize(0.05);
  paveptB3L->SetTextAlign(12);
  paveptB3L->AddText(Form("#it{B}_{3,#Lambda}: #it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb3Lambda));


  //Chose the set of parameterisation for the radii
  Short_t ip = 1;

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
  legB2data->AddEntry(gB2vsR_pp7TeV_sys[ip], "pp #sqrt{#it{s}} = 7 TeV [PRC 97, 024615 (2018)]", "pf");
 
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
  pavept->Draw();
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
  legB3data->AddEntry(gB3vsR_PbPb5TeV_sys[ip], "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, prelim.", "pf");
  legB3data->AddEntry(gB3vsR_PbPb276TeV_sys[ip], "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV [PRC 93, 0249717 (2016)]", "pf");
  legB3data->AddEntry(gB3vsR_pp7TeV_sys[ip], "pp #sqrt{#it{s}} = 7 TeV [arXiv:1709.08522]", "pf");
  legB3data->AddEntry(gB3LambdavsR_PbPb276TeV_sys[ip], "#it{B}_{3,#Lambda}, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV [PLB 754, 360-372 (2016)]", "pf");

  nl = 3;
  TLegend * legB3blast = new TLegend(0.1, 0.55-nl*0.04, 0.45, 0.55, "Blast-Wave (#pi,K,p) + GSI-Heid. (T = 156 MeV)");
  legB3blast->SetFillStyle(0);
  legB3blast->SetTextSize(0.035);
  legB3blast->SetBorderSize(0);
  legB3blast->AddEntry(gBlastB3vsR_PbPb276TeV[ip], "#it{B}_{3}, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "l");
  legB3blast->AddEntry(gBlastB3LambdavsR_PbPb276TeV[ip], "#it{B}_{3,#Lambda}, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "l");

  nl = 4;
  TLegend * legB3coal = new TLegend(0.1, 0.4-nl*0.04, 0.45, 0.4, "Coalescence");
  legB3coal->SetFillStyle(0);
  legB3coal->SetTextSize(0.035);
  legB3coal->SetBorderSize(0);
  legB3coal->AddEntry(hB3_coalescence, "#it{B}_{3} coalesc., #it{r}(^{3}He) = 2.48 fm", "l");
  legB3coal->AddEntry(hB3_coalescence_pointlike, "#it{B}_{3} coalesc., #it{r}(^{3}He) = 0 (point-like)", "l");
  legB3coal->AddEntry(hB3L_coalescence, "#it{B}_{3,#Lambda} coalesc., #it{r}(^{3}_{#Lambda}H) = 6.8 fm", "l");

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
  hB3_coalescence_pointlike->Draw("lsame");
  //hB3L_coalescence->Draw("lsame");
  gBlastB3vsR_PbPb276TeV[ip]->Draw("samel");
  //gBlastB3LambdavsR_PbPb276TeV->Draw("samel");
  gB3vsR_PbPb5TeV_sys[ip]->Draw("p3");
  gB3vsR_PbPb5TeV[ip]->Draw("samep");
  gB3vsR_PbPb276TeV_sys[ip]->Draw("samep3");
  gB3vsR_PbPb276TeV[ip]->Draw("samep");
  gB3vsR_pp7TeV_sys[ip]->Draw("samep2");
  gB3vsR_pp7TeV[ip]->Draw("samep");
  // gB3LambdavsR_PbPb276TeV_sys->Draw("samep2");
  // gB3LambdavsR_PbPb276TeV->Draw("samep");


  cb3opta->cd(2);
  legB3data->Draw();
  legB3coal->Draw();
  legB3blast->Draw();

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

  
  return 0;  
}



void MakePaperFigure2(Bool_t plotLinX, Double_t pToA,
		      TF1 *Cd_coalescence, TF1* Cd_coalescence_pointlike, TF1* Cd_coalescence_radius1third,  TF1* Cd_coalescence_largeradius,
		      TGraphErrors * hB2_coalescence, TGraphErrors * hB2_coalescence_pointlike, TGraphErrors * hB2_coalescence_radius1third, TGraphErrors * hB2_coalescence_largeradius) {
  //
  // Create the (pure theory) figure which plots <C_d> and B2 vs R
  // for different radii (PLOT COALESCENCE ONLY)
  //
  TCanvas * coalcanv = new TCanvas("coalcanv", "coalescence", 1600, 1000);
  coalcanv->SetBottomMargin(0.02);
  coalcanv->SetTopMargin(0.02);
  coalcanv->SetLeftMargin(0.15);
  coalcanv->SetRightMargin(0.02);
  coalcanv->Divide(2,1);

  TH2D * frame_cd = new TH2D("frame_cd", "#LTC_{d}#GT vs radius; #it{R} (fm); #LT#it{C}_{d}#GT", 1000, 0.01, 6.0, 2e5, 0, 1.2);
  frame_cd->GetXaxis()->SetTitleSize(0.05);
  frame_cd->GetYaxis()->SetTitleSize(0.05);
  if (plotLinX) frame_cd->GetXaxis()->SetRangeUser(0.01, 8.5);
  else  frame_cd->GetXaxis()->SetRangeUser(0.1, 10.5);
  
  TH2D * frame_coal = new TH2D("frame_coal", "B_{2} vs radius; #it{R} (fm); #it{B}_{2} (GeV^{2}/#it{c}^{3})", 1000, 0.01, 6.0, 2e5, 1.e-4, 0.1);
  frame_coal->GetXaxis()->SetTitleSize(0.05);
  frame_coal->GetYaxis()->SetTitleSize(0.05);
  if (plotLinX) frame_coal->GetXaxis()->SetRangeUser(0.01, 8.5);
  else  frame_coal->GetXaxis()->SetRangeUser(0.1, 10.5);
  
  TLegend * legB2_coal;
  if (plotLinX) legB2_coal = new TLegend(0.4, 0.95-4*0.04, 0.8, 0.95);
  else legB2_coal = new TLegend(0.2, 0.15, 0.75, 0.15+4*0.03);
  legB2_coal->SetFillStyle(0);
  legB2_coal->SetTextSize(0.035);
  legB2_coal->SetBorderSize(0);
  
  legB2_coal->AddEntry(hB2_coalescence_pointlike, "#it{B}_{2} coalesc., #it{r_{d}} = 0 (point-like)", "l");
  legB2_coal->AddEntry(hB2_coalescence_radius1third, "#it{B}_{2} coalesc., #it{r_{d}} = 0.3 fm", "l");
  legB2_coal->AddEntry(hB2_coalescence, "#it{B}_{2} coalesc., #it{r_{d}} = 3.2 fm", "l");
  legB2_coal->AddEntry(hB2_coalescence_largeradius, "#it{B}_{2} coalesc., #it{r_{d}} = 10 fm", "l");
  // add pT over A here as well

  coalcanv->cd(1);
  gPad->SetTickx();
  gPad->SetTicky();
  frame_cd->Draw();
  Cd_coalescence->Draw("same");
  Cd_coalescence_pointlike->Draw("same");
  Cd_coalescence_radius1third->Draw("same");
  Cd_coalescence_largeradius->Draw("same");
  
  coalcanv->cd(2);
  gPad->SetLogy();
  gPad->SetTickx();
  gPad->SetTicky();
  frame_coal->Draw();
  hB2_coalescence->Draw("same");
  hB2_coalescence_pointlike->Draw("same");
  hB2_coalescence_radius1third->Draw("same");
  hB2_coalescence_largeradius->Draw("same");
  legB2_coal->Draw();

  TPaveText * paveptCoalCanv = new TPaveText(0.5, 0.6, 0.9, 0.7, "NDC");
  paveptCoalCanv->SetFillStyle(0);
  paveptCoalCanv->SetTextFont(42);
  paveptCoalCanv->SetBorderSize(0);
  paveptCoalCanv->SetTextSize(0.04);
  paveptCoalCanv->SetTextAlign(12);
  paveptCoalCanv->AddText(Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToA));
  paveptCoalCanv->Draw();


  coalcanv->SaveAs("Paper/theory_coalescence_Cd_B2.eps");
  coalcanv->SaveAs("Paper/theory_coalescence_Cd_B2.png");


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
  TH2D * hframe = new TH2D("hframeFig4", "B_{2} vs radius; #it{R} (fm); #it{B}_{2} (GeV^{2}/#it{c}^{3})", 1000, 0.01, 6.0, 2000, 1.e-4, 0.1);
  hframe->GetXaxis()->SetTitleSize(0.06);
  hframe->GetYaxis()->SetTitleSize(0.06);
  hframe->GetYaxis()->SetTitleOffset(1.3);
  hframe->GetXaxis()->SetTitleOffset(0.8);
  hframe->GetXaxis()->SetLabelSize(0.05);
  hframe->GetYaxis()->SetLabelSize(0.05);
  if (plotLinX) hframe->GetXaxis()->SetRangeUser(0.01, 8.5);
  else  hframe->GetXaxis()->SetRangeUser(0.1, 10.5);

  TH2D * hframe3 = new TH2D("hframe3Fig4", "B_{3} vs radius; #it{R} (fm); #it{B}_{3} (GeV^{4}/#it{c}^{6})", 1000, 0.01, 6.0, 2000, 1.e-9, 1.e-1);
  hframe3->GetXaxis()->SetTitleSize(0.06);
  hframe3->GetYaxis()->SetTitleSize(0.06);
  hframe3->GetYaxis()->SetTitleOffset(1.3);
  hframe3->GetXaxis()->SetTitleOffset(0.8);
  hframe3->GetXaxis()->SetLabelSize(0.05);
  hframe3->GetYaxis()->SetLabelSize(0.05);
  if (plotLinX) hframe3->GetXaxis()->SetRangeUser(0.01, 8.5);
  else  hframe3->GetXaxis()->SetRangeUser(0.1, 10.5);

  TH2D * hframe3L = new TH2D("hframe3LFig4", "B_{3,#Lambda} vs radius; #it{R} (fm); #it{B}_{3,#Lambda} (GeV^{4}/#it{c}^{6})", 1000, 0.01, 6.0, 2000, 1.e-9, 1.e-1);
  hframe3L->GetXaxis()->SetTitleSize(0.06);
  hframe3L->GetYaxis()->SetTitleSize(0.06);
  hframe3L->GetYaxis()->SetTitleOffset(1.3);
  hframe3L->GetXaxis()->SetTitleOffset(0.8);
  hframe3L->GetXaxis()->SetLabelSize(0.05);
  hframe3L->GetYaxis()->SetLabelSize(0.05);
  if (plotLinX) hframe3L->GetXaxis()->SetRangeUser(0.01, 8.5);
  else  hframe3L->GetXaxis()->SetRangeUser(0.1, 10.5);


  //define particle label
  TPaveText * paveLab2 = new TPaveText(0.8, 0.8, 0.9, 0.9, "NDC");
  paveLab2->SetFillStyle(0);
  paveLab2->SetTextFont(42);
  paveLab2->SetBorderSize(0);
  paveLab2->SetTextSize(0.1);
  paveLab2->SetTextAlign(12);
  paveLab2->AddText("#bf{d}");

  TPaveText * paveLab3 = new TPaveText(0.77, 0.8, 0.9, 0.9, "NDC");
  paveLab3->SetFillStyle(0);
  paveLab3->SetTextFont(42);
  paveLab3->SetBorderSize(0);
  paveLab3->SetTextSize(0.1);
  paveLab3->SetTextAlign(12);
  paveLab3->AddText("#bf{^{3}He}");

  TPaveText * paveLab3L = new TPaveText(0.75, 0.8, 0.9, 0.9, "NDC");
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

   
  TCanvas * cr4 = new TCanvas("cr4", "compare thermal with coalescence", 1000, 1000);
  cr4->SetBottomMargin(0.02);
  cr4->SetTopMargin(0.01);
  cr4->SetLeftMargin(0.12);
  cr4->SetRightMargin(0.02);
  cr4->Divide(2,2);
  
  cr4->cd(1);
  gPad->SetLogy();
  hframe->Draw();
  gBlastB2vsR_PbPb276TeV[1]->Draw("samel");
  gB2vsR_pp7TeVINELg0_sys[1]->Draw("samep2");
  gB2vsR_pp7TeVINELg0[1]->Draw("samepz");
  gB2vsR_PbPb276TeV_sys[1]->Draw("samep3");
  gB2vsR_PbPb276TeV[1]->Draw("samepz");
  //
  TLegend * legB2 = new TLegend(0.2, 0.95, 0.6, 0.8, "");
  legB2->SetFillStyle(0);
  legB2->SetTextSize(0.03);
  legB2->SetBorderSize(0);
  legB2->AddEntry(hB2_coalescence, "#it{B}_{2} coalesc., #it{r}(d) = 3.2 fm", "l");
  legB2->AddEntry(gB2vsR_PbPb276TeV_sys[1], "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV [PRC 93, 0249717 (2016)]", "pf");
  legB2->AddEntry(gB2vsR_pp7TeVINELg0_sys[1], "pp #sqrt{#it{s}} = 7 TeV [PRC 97, 024615 (2018)]", "pf");
  legB2->AddEntry(gBlastB2vsR_PbPb276TeV[1], "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV, BW + GSI (T = 156 MeV)", "l");
  //legB2->Draw();
  //
  hB2_coalescence->Draw("l");
  pavept->Draw();
  paveLab2->Draw();
  
  cr4->cd(3);
  gPad->SetLogy();
  hframe3->Draw();
  hB3_coalescence->Draw("l");
  gBlastB3vsR_PbPb276TeV[1]->Draw("samel");
  gB3vsR_PbPb276TeV_sys[1]->Draw("samep3");
  gB3vsR_PbPb276TeV[1]->Draw("samepz");
  gB3vsR_pp7TeV_sys[1]->Draw("samep2");
  gB3vsR_pp7TeV[1]->Draw("samepz");
  //
  TLegend * legB3 = new TLegend(0.2, 0.95, 0.6, 0.8, "");
  legB3->SetFillStyle(0);
  legB3->SetTextSize(0.03);
  legB3->SetBorderSize(0);
  legB3->AddEntry(hB3_coalescence, "#it{B}_{3} coalesc., #it{r}(^{3}He) = 2.48 fm", "l");
  legB3->AddEntry(gB3vsR_PbPb276TeV_sys[1], "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV [PRC 93, 0249717 (2016)]", "pf");
  legB3->AddEntry(gB3vsR_pp7TeV_sys[1], "pp #sqrt{#it{s}} = 7 TeV [arXiv:1709.08522]", "pf");
  legB3->AddEntry(gBlastB3vsR_PbPb276TeV[1], "#it{B}_{3}, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV, BW + GSI (T = 156 MeV)", "l");
  //legB3->Draw();
  //
  paveptB3->Draw();
  paveLab3->Draw();

  cr4->cd(4);
  gPad->SetLogy();
  hframe3L->Draw();
  hB3L_coalescence_largeradius->Draw("samel");
  hB3L_coalescence->Draw("l");
  gBlastB3LambdavsR_PbPb276TeV[1]->Draw("samel");
  gB3LambdavsR_PbPb276TeV_sys[1]->Draw("samep2");
  gB3LambdavsR_PbPb276TeV[1]->Draw("samep");
  paveptB3L->Draw();
  paveLab3L->Draw();
  //
  TLegend * legB3Lambda = new TLegend(0.2, 0.95, 0.6, 0.8, "");
  legB3Lambda->SetFillStyle(0);
  legB3Lambda->SetTextSize(0.03);
  legB3Lambda->SetBorderSize(0);
  legB3Lambda->AddEntry(hB3L_coalescence, "#it{B}_{3,#Lambda} coalesc., #it{r}(^{3}_{#Lambda}H) = 6.8 fm", "l");
  legB3Lambda->AddEntry(gB3LambdavsR_PbPb276TeV_sys[1], "#it{B}_{3,#Lambda}, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV [PLB 754, 360-372 (2016)]", "pf");
  legB3Lambda->AddEntry(gBlastB3LambdavsR_PbPb276TeV[1], "#it{B}_{3,#Lambda}, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV, BW + GSI (T = 156 MeV)", "l");
  //legB3Lambda->Draw();
  //

  TLegend * masterLeg = new TLegend(0.1, 0.3, 0.5, 0.9, "");
  masterLeg->SetFillStyle(0);
  masterLeg->SetTextSize(0.05);
  masterLeg->SetBorderSize(0);
  masterLeg->AddEntry(gB2vsR_PbPb276TeV_sys[1], "ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "pf");
  masterLeg->AddEntry(gB2vsR_pp7TeVINELg0_sys[1], "ALICE, pp #sqrt{#it{s}} = 7 TeV (INEL>0)", "pf");
  masterLeg->AddEntry(gBlastB2vsR_PbPb276TeV[1], "BW + GSI-Heidelberg (#it{T}_{chem} = 156 MeV)", "l");
  //masterLeg->AddEntry(gBlastB2vsR_PbPb276TeV[1], "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "");
  masterLeg->AddEntry(hB2_coalescence, "#it{B}_{#it{A}} coalescence", "l");
  masterLeg->AddEntry(hB2_coalescence, "#it{r} (d) = 3.2 fm", "");
  masterLeg->AddEntry(hB3_coalescence, "#it{r} (^{3}He) = 2.48 fm", "");
  masterLeg->AddEntry(hB3L_coalescence, "#it{r} (^{3}_{#Lambda}H) = 6.8 fm", "");
  masterLeg->AddEntry(hB3L_coalescence_largeradius, "#it{r} (^{3}_{#Lambda}H) = 14.1 fm", "l");

  cr4->cd(2);
  masterLeg->Draw();
  
  // cr4->SaveAs("Paper/compareThermalAndCoalescence.eps");
  // cr4->SaveAs("Paper/compareThermalAndCoalescence.png");

}


void getRadiusFromParameterisation(Double_t * multi, Double_t * radius, Int_t paramSet)
{
  // parameterisation to convert multiplicity into Rside
  // this is the parameterisation for the pion radius
  //
  if (!multi || !radius) return;
  Double_t  multi3 = TMath::Power(multi[0], 1./3.);  
  //
  // Here is the crucial mapping between HBT radii and multi^(1/3)
  //
  // VERSION (5th February 2018):
  // We assume 0.85fm for pp as Kfir (at dNdeta 6.01)
  // We assume R=4.5fm for central Pb-Pb based on arXiv:1012.4035 (figure 2)
  // We interpolate linearly and take the highest kT bin, because
  // a pT/A of 0.8 GeV/c corresponds to mT=1.23GeV which is a kT of pions of 1.2GeV
  // radiusVal = 0.177825 + 0.36733 * multi3;

  // VERSION (17th May 2018):
  // We fit linearly the ALICE data at the kT = 0.887 
  Double_t radiusVal = 0.0;
  if (paramSet==2) {
    //manual hack to have the data points fall onto the U. Heinz curve for 3He
    radiusVal = 0.190 + 0.380 * multi3;
  } else  if (paramSet==1) {
    //manual hack to have the data points fall onto the U. Heinz curve for deuteron
    radiusVal = 0.0 + 0.472949 * multi3; //radius = 0 at 0 dN/deta
    //radiusVal = 0.07412 + 0.46637 * multi3; //radius = 0.85 fm for pp, dN/deta = 4.60 (INELg0)
    //radiusVal = -0.009 + 0.4738 * multi3; //radius = 0.85 fm for pp, dN/deta = 5.98 (INELg0>0)
    //radiusVal = -0.3949 + 0.507865 * multi3; //most central and peripheral PbPb 2.76 TeV on the curve
  } else {
    //fit to the HBT data, kT = 0.887
    radiusVal = 0.128 + 0.339 * multi3; 
  }
  
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

TF1 * MakeB2TheoryGraphQMfactor(Double_t objSize)
{
  TF1 * funcCd = new TF1(Form("funcCd_%i", TMath::Nint(objSize*10)), "1 / TMath::Power(1 + [0]*[0]/(4*x*x), 1.5)", 0., 15.);
  funcCd->SetParameter(0, objSize);
  return funcCd;
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

void convertMultiToRadius(TGraphErrors * graph, Int_t paramSet)
{
  //convert multiplicity dN/deta into radius with parameterisation
  if (!graph) return;
  
  for (Int_t ip = 0; ip < graph->GetN(); ip++){
    Double_t xold[2];
    xold[0] = graph->GetX()[ip];
    xold[1] = graph->GetEX()[ip];
    Double_t xnew[2] = {-1.0, -1.0};
    getRadiusFromParameterisation(xold, xnew, paramSet);
    graph->GetX()[ip] = xnew[0];
    graph->GetEX()[ip] = xnew[1];
  }

  return;
}


void convertMultiToRadius(TGraphAsymmErrors * graph, Int_t paramSet)
{
  //convert multiplicity dN/deta into radius with parameterisation
  if (!graph) return;
  
  for (Int_t ip = 0; ip < graph->GetN(); ip++){
    Double_t xold[2];
    xold[0] = graph->GetX()[ip];
    xold[1] = graph->GetEXlow()[ip];
    Double_t xnew[2] = {-1.0, -1.0};
    getRadiusFromParameterisation(xold, xnew, paramSet);
    graph->GetX()[ip] = xnew[0];
    graph->GetEXlow()[ip] = xnew[1];
    graph->GetEXhigh()[ip] = xnew[1];
	
  }

  return;
}

//---------------------------------------------------------
//------------------------------ ALICE data B2 --- pp 7 TeV

TGraphErrors * getB2_pp7TeV(Bool_t plotSys, Double_t pToA, Int_t paramSet)
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

  convertMultiToRadius(graph, paramSet);

  graph->SetMarkerColor(kGreen+2);
  graph->SetLineColor(kGreen+2);
  graph->SetFillColorAlpha(kGreen+2, 0.1);
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.2);
  graph->SetMarkerStyle(20);
  return graph;
  
}

TGraphErrors * getB2_pp7TeVINELg0(Bool_t plotSys, Double_t pToA, Int_t paramSet)
{
  if (pToA<0.74 || pToA>0.76 ) return 0x0;
  
  //for multi INEL, INEL>0 Eur. Phys. J. C 77 (2017) 33, Link: https://link.springer.com/article/10.1140/epjc/s10052-016-4571-1
  const Int_t nP = 1;

  
  //INEL   Double_t dndeta[nP] = {4.60}; Double_t dndetaErr[nP] = {0.34};
  //INEL >0
  Double_t dndeta[nP] = {5.98};
  Double_t dndetaErr[nP] = {0.09};

  //from ALICE, Phys. Rev. C 97, 024615 (2018) -- arXiv:1709.08522
  //the value is for INEL, so we recale to go to INEL>0
  Float_t factorINELtoINELg0 = 4.60 / 5.98;
  Double_t y[nP] = {0.01842*factorINELtoINELg0};
  Double_t ey[nP] = {0.00074*factorINELtoINELg0};
  Double_t yy[nP] = {0.00275*factorINELtoINELg0};
  TGraphErrors * graph = 0x0;
  if (plotSys) graph = new TGraphErrors(nP, dndeta, y, dndetaErr, yy);
  else graph = new TGraphErrors(nP, dndeta, y, dndetaErr, ey);
  graph->SetName("B2_pp7TeVINELg0");
  if (!graph) {
    Printf("Error: cannot retrieve graph. Check pt bin requested.");
    return NULL;
  }

  convertMultiToRadius(graph, paramSet);

  graph->SetMarkerColor(kGreen+2);
  graph->SetLineColor(kGreen+2);
  graph->SetFillColorAlpha(kGreen+2, 0.1);
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.2);
  graph->SetMarkerStyle(20);
  return graph;
  
}

//---------------------------------------------------------
//------------------------------ ALICE data B2 --- pPb 5 TeV
TGraphErrors * getB2_pPb5TeV(Bool_t plotSys, Double_t pToA, Int_t paramSet)
{

  TFile * f0 = TFile::Open("oB2vsNch.root");
  if (!f0) return NULL;

  TString gName = Form("B2_pPb_pToA=%4.3f%s", pToA, (plotSys? "_sys" : ""));
  TGraphErrors * graph = (TGraphErrors*) f0->Get(gName.Data());
  if (!graph) {
    Printf("Error: cannot retrieve graph. Check pt bin requested.");
    return NULL;
  }

  convertMultiToRadius(graph, paramSet);

  graph->SetMarkerColor(kBlue);
  graph->SetLineColor(kBlue);
  graph->SetFillColorAlpha(kBlue, 0.1);
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.2);
  graph->SetMarkerStyle(20);
  return graph;
  
}


//---------------------------------------------------------
//------------------------------ ALICE data B2 --- PbPb 5 TeV
TGraphErrors * getB2_PbPb5TeV(Bool_t plotSys, Double_t pToA, Int_t paramSet)
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

  convertMultiToRadius(graph, paramSet);

  graph->SetMarkerColor(kRed);
  graph->SetLineColor(kRed);
  graph->SetFillColorAlpha(kRed, 0.1);
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.2);
  graph->SetMarkerStyle(20);
  return graph;
  
}


//---------------------------------------------------------
//------------------------------ ALICE data B2 --- PbPb 2.76 TeV
TGraphErrors * getB2_PbPb276TeV(Bool_t plotSys, Double_t pToA, Int_t paramSet)
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

  convertMultiToRadius(graph, paramSet);

  graph->SetMarkerColor(kRed+2);
  graph->SetLineColor(kRed+2);
  graph->SetFillColorAlpha(kRed+2, 0.1);  
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.2);
  graph->SetMarkerStyle(20);
  return graph;
  
}


//---------------------------------------------------------
//------------------------------ ALICE data B3 --- pp 7 TeV

TGraphAsymmErrors * getB3_pp7TeVINELg0(Bool_t plotSys, Double_t pToAb3pp, Int_t paramSet)
{
  //from ALICE, Phys. Rev. C 97, 024615 (2018) -- arXiv:1709.08522
  // x value for multi is INEL >0
  TFile * f0 = TFile::Open("B3pToA_pp7TeV.root");
  if (!f0) return NULL;

  TString gName = Form("B3_pp7TeV_pToA=%4.3f%s", pToAb3pp, (plotSys? "_sys" : ""));
  TGraphAsymmErrors * graph = (TGraphAsymmErrors*) f0->Get(gName.Data());
  if (!graph) {
    Printf("Error: cannot retrieve graph. Check pt bin requested.");
    return NULL;
  }
  //the value is for INEL, so we rescale to go to INEL>0
  Float_t factorINELtoINELg0 = 4.60 / 5.98;

  graph->GetY()[0] *= (factorINELtoINELg0*factorINELtoINELg0);
  graph->GetEYlow()[0] *= (factorINELtoINELg0*factorINELtoINELg0);
  graph->GetEYhigh()[0] *= (factorINELtoINELg0*factorINELtoINELg0);
    
  convertMultiToRadius(graph, paramSet);

  graph->SetMarkerColor(kTeal-5);
  graph->SetLineColor(kTeal-5);
  graph->SetFillColorAlpha(kTeal-5, 0.1);
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.2);
  graph->SetMarkerStyle(20);
  return graph;
  
}


//---------------------------------------------------------
//------------------------------ ALICE data B3 --- PbPb 5 TeV
TGraphErrors * getB3_PbPb5TeV(Bool_t plotSys, Double_t pToAb3, Int_t paramSet)
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

  convertMultiToRadius(graph, paramSet);

  graph->SetMarkerColor(kRed);
  graph->SetLineColor(kRed);
  graph->SetFillColorAlpha(kRed, 0.1);
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.2);
  graph->SetMarkerStyle(20);
  return graph;
  
}


//---------------------------------------------------------
//------------------------------ ALICE data B3 --- PbPb 2.76 TeV
TGraphErrors * getB3_PbPb276TeV(Bool_t plotSys, Double_t pToAb3, Int_t paramSet)
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

  convertMultiToRadius(graph, paramSet);
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
TGraphAsymmErrors * getB3Lambda_PbPb276TeV(Bool_t plotSys, Double_t pToAb3Lambda, Int_t paramSet)
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

  convertMultiToRadius(graph, paramSet);

  graph->SetMarkerColor(kAzure-7);
  graph->SetLineColor(kAzure-7);
  graph->SetFillColorAlpha(kAzure-7, 0.1);  
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.2);
  graph->SetMarkerStyle(21);
  return graph;
  
}


//---------------------------------------------------------
//---------------------- Blast wave + thermal PbPb 2.76 TeV
TGraphAsymmErrors * getBlastB2_PbPb276TeV(Bool_t plotSys, Double_t pToA, Int_t paramSet)
{
  // Published proton yield from Phys. Rev. C 88 (2013) 044910
  // Blast wave params from π,K,p published
  // d/p from thermal model T = 156 MeV
 
  TGraphAsymmErrors* graph = (TGraphAsymmErrors *) generateBWpredictionsB2("PbPb276TeV", "rms", "deuteron", pToA);
  convertMultiToRadius(graph, paramSet);
  
  graph->SetMarkerColor(kMagenta+2);
  graph->SetLineColor(kMagenta+2);
  graph->SetFillColorAlpha(kMagenta+2, 0.1);  
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.5);
  graph->SetMarkerStyle(33);
  graph->SetLineWidth(3);
  graph->SetLineStyle(2);
  return graph;
}

//---------------------------------------------------------
//---------------------- Blast wave + thermal PbPb 5.02 TeV
TGraphAsymmErrors * getBlastB2_PbPb502TeV(Bool_t plotSys, Double_t pToA, Int_t paramSet)
{
  // Preliminary (~final) proton yields and BW params from:
  //https://gitlab.cern.ch/njacazio/SpectraAnalysisRun2/tree/master/results/spectra/spectra-pag/Preliminaries/QM2017
  // d/p from thermal model T = 156 MeV
 
  TGraphAsymmErrors* graph = (TGraphAsymmErrors *) generateBWpredictionsB2("PbPb502TeV", "rms", "deuteron", pToA);
  convertMultiToRadius(graph, paramSet);
  
  graph->SetMarkerColor(kMagenta+2);
  graph->SetLineColor(kMagenta+2);
  graph->SetFillColorAlpha(kMagenta+2, 0.1);  
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.5);
  graph->SetMarkerStyle(33);
  graph->SetLineWidth(3);
  graph->SetLineStyle(2);
  return graph;
  
}

//---------------------------------------------------------
//---------------------- Blast wave + thermal pPb 5 TeV
TGraphAsymmErrors * getBlastB2_pPb502TeV(Bool_t plotSys, Double_t pToA, Int_t paramSet)
{
  // Published proton yield  
  // Blast wave params from π,K,p published
  // d/p from thermal model T = 156 MeV
 
  TGraphAsymmErrors* graph = (TGraphAsymmErrors *) generateBWpredictionsB2("pPb502TeV", "rms", "deuteron", pToA);
  convertMultiToRadius(graph, paramSet);
  
  graph->SetMarkerColor(kMagenta);
  graph->SetLineColor(kMagenta);
  graph->SetFillColorAlpha(kMagenta, 0.1);  
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.5);
  graph->SetMarkerStyle(33);
  graph->SetLineWidth(3);
  graph->SetLineStyle(2);
  return graph;
  
}

//---------------------------------------------------------
//---------------------- Blast wave + thermal pp 7 TeV
TGraphAsymmErrors * getBlastB2_pp7TeV(Bool_t plotSys, Double_t pToA, Int_t paramSet)
{
  // Final proton yield from long paper pp7 TeV
  // Blast wave params from π,K,p published
  // d/p from thermal model T = 156 MeV
 
  TGraphAsymmErrors* graph = (TGraphAsymmErrors *) generateBWpredictionsB2("pp7TeV", "rms", "deuteron", pToA);
  convertMultiToRadius(graph, paramSet);
  
  graph->SetMarkerColor(kMagenta-2);
  graph->SetLineColor(kMagenta-2);
  graph->SetFillColorAlpha(kMagenta-2, 0.1);  
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.5);
  graph->SetMarkerStyle(33);
  graph->SetLineWidth(3);
  graph->SetLineStyle(2);
  return graph;
  
}


//---------------------------------------------------------
//---------------------- Blast wave + thermal PbPb 2.76 TeV
TGraphAsymmErrors * getBlastB3_PbPb276TeV(Bool_t plotSys, Double_t pToAb3, Int_t paramSet)
{
  // Published proton yield from Phys. Rev. C 88 (2013) 044910
  // Blast wave params from π,K,p published
  // 3He/p from thermal model T = 156 MeV
 
  TGraphAsymmErrors* graph = (TGraphAsymmErrors *) generateBWpredictionsB2("PbPb276TeV", "rms", "He3", pToAb3);
  convertMultiToRadius(graph, paramSet);
  
  graph->SetMarkerColor(kMagenta+1);
  graph->SetLineColor(kMagenta+1);
  graph->SetFillColorAlpha(kMagenta+1, 0.1);  
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.5);
  graph->SetMarkerStyle(33);
  graph->SetLineWidth(5);
  graph->SetLineStyle(2);
  return graph;
  
}

//---------------------- Blast wave + thermal PbPb 2.76 TeV -- Hypertriton
TGraphAsymmErrors * getBlastB3Lambda_PbPb276TeV(Bool_t plotSys, Double_t pToAb3, Int_t paramSet)
{
  // Published proton yield from Phys. Rev. C 88 (2013) 044910
  // Blast wave params from π,K,p published
  // 3He/p from thermal model T = 156 MeV
  // s3 from hyper-triton paper
 
  TGraphAsymmErrors* graph = (TGraphAsymmErrors *) generateBWpredictionsB2("PbPb276TeV", "rms", "hyper-triton", pToAb3);
  convertMultiToRadius(graph, paramSet);
  
  graph->SetMarkerColor(kAzure-7);
  graph->SetLineColor(kAzure-7);
  graph->SetFillColorAlpha(kAzure-7, 0.1);  
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.5);
  graph->SetMarkerStyle(24);
  graph->SetLineWidth(5);
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
