/////////////////////////////////////////////////////////////////////////////////////////////////////////
//The Frascati project - thermal + blastwave
/////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "TGraphErrors.h"
#include "./GetALICEdata.cpp"
#include "./generateBWpredictionsB2.C"

//---------------------------------------------------------------------------------------------------------------------
//thermal model + blast-wave predictions up to A=4
//---------------------------------------------------------------------------------------------------------------------
TGraphAsymmErrors * getBlastB2_PbPb276TeV(Bool_t plotSys = 0, Double_t pToA = 0.75, Int_t paramSet = 0);
TGraphAsymmErrors * getBlastB2_PbPb502TeV(Bool_t plotSys = 0, Double_t pToA = 0.75, Int_t paramSet = 0);
TGraphAsymmErrors * getBlastB2_pPb502TeV(Bool_t plotSys = 0, Double_t pToA = 0.75, Int_t paramSet = 0);
TGraphAsymmErrors * getBlastB2_pp7TeV(Bool_t plotSys = 0, Double_t pToA = 0.75, Int_t paramSet = 0);

TGraphAsymmErrors * getBlastB3_PbPb276TeV(Bool_t plotSys = 0, Double_t pToAb3 = 0.733, Int_t paramSet = 0, Bool_t convertToR = kTRUE);
TGraphAsymmErrors * getBlastB3_pPb5TeV(Bool_t plotSys = 0, Double_t pToAb3 = 0.733, Int_t paramSet = 0, Bool_t convertToR = kTRUE, Bool_t plotcSHM = 0);
TGraphAsymmErrors * getBlastB3_pp7TeV(Bool_t plotSys = 0, Double_t pToAb3 = 0.733, Int_t paramSet = 0, Bool_t convertToR = kTRUE);

TGraphAsymmErrors * getBlastB3Lambda_PbPb276TeV(Bool_t plotSys = 0, Double_t pToAb3Lambda = 1.0, Int_t paramSet = 0);
TGraphAsymmErrors * getBlastB3Lambda_pp7TeV(Bool_t plotSys = 0, Double_t pToAb3Lambda = 1.0, Int_t paramSet = 0);

TGraphAsymmErrors * getBlastB4_PbPb276TeV(Bool_t plotSys = 0, Double_t pToAb4 = 0.75, Int_t paramSet = 0);
TGraphAsymmErrors * getBlastB4_pp7TeV(Bool_t plotSys = 0, Double_t pToAb4 = 0.75, Int_t paramSet = 0);

TGraphAsymmErrors * getBlastB4Lambda_PbPb276TeV(Bool_t plotSys = 0, Double_t pToAb3 = 0.75, Int_t paramSet = 0);
TGraphAsymmErrors * getBlastB4Lambda_pp7TeV(Bool_t plotSys = 0, Double_t pToAb3 = 0.75, Int_t paramSet = 0);

TGraphAsymmErrors * getBAthermalBlast(TString system = "PbPb276TeV", TString particle = "deuteron", Double_t pToA = 0.75, Color_t color = kBlack, Int_t paramSet = 1, Bool_t convertToR = kTRUE);

//---------------------------------------------------------------------------------------------------------------------
TGraphAsymmErrors * getBAthermalBlast(TString system, TString particle, Double_t pToA, Color_t color, Int_t paramSet, Bool_t convertToR)
{
  TGraphAsymmErrors* graph = (TGraphAsymmErrors *) generateBWpredictionsB2(system.Data(), "rms", particle.Data(), pToA);
  //parameterisation 1 of the radius to be used in order to have data points fall onto the U. Heinz curve for deuteron
  if (convertToR) convertMultiToRadius(graph, paramSet);
  //set style
  graph->SetMarkerColor(color);
  graph->SetLineColor(color);
  graph->SetFillColorAlpha(color, 0.1);  
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.5);
  graph->SetMarkerStyle(33);
  graph->SetLineWidth(2);
  graph->SetLineStyle(2);

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
  graph->SetLineWidth(2);
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
  graph->SetLineWidth(2);
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
  graph->SetLineWidth(2);
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
  graph->SetLineWidth(2);
  graph->SetLineStyle(2);
  return graph;
  
}

//---------------------------------------------------------
//---------------------- Blast wave + thermal pp 7 TeV
TGraphAsymmErrors * getBlastB3_pp7TeV(Bool_t plotSys, Double_t pToAb3, Int_t paramSet, Bool_t convertToR)
{
  // Final proton yield from long paper pp7 TeV
  // Blast wave params from π,K,p published
  // d/p from thermal model T = 156 MeV
 
  TGraphAsymmErrors* graph = (TGraphAsymmErrors *) generateBWpredictionsB2("pp7TeV", "rms", "He3", pToAb3);
  if (convertToR) convertMultiToRadius(graph, paramSet);
  
  graph->SetMarkerColor(kMagenta-2);
  graph->SetLineColor(kMagenta-2);
  graph->SetFillColorAlpha(kMagenta-2, 0.1);  
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.5);
  graph->SetMarkerStyle(33);
  graph->SetLineWidth(2);
  graph->SetLineStyle(2);
  return graph;
  
}

//---------------------------------------------------------
//---------------------- Blast wave + thermal pPb 5 TeV
TGraphAsymmErrors * getBlastB3_pPb5TeV(Bool_t plotSys, Double_t pToAb3, Int_t paramSet, Bool_t convertToR, Bool_t plotcSHM)
{
  // Final proton yield from paper spectra p-Pb 5 TeV
  // Blast wave params from π,K,p published
  // 3He/p from thermal model 
 
  TGraphAsymmErrors* graph = (TGraphAsymmErrors *) generateBWpredictionsB2("pPb502TeV", "rms", "He3", pToAb3, plotcSHM);
  if (convertToR) convertMultiToRadius(graph, paramSet);
  
  graph->SetMarkerColor(kMagenta-2);
  graph->SetLineColor(kMagenta-2);
  graph->SetFillColorAlpha(kMagenta-2, 0.1);  
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.5);
  graph->SetMarkerStyle(33);
  graph->SetLineWidth(2);
  graph->SetLineStyle(2);
  return graph;
  
}

//---------------------------------------------------------
//---------------------- Blast wave + thermal PbPb 2.76 TeV
TGraphAsymmErrors * getBlastB3_PbPb276TeV(Bool_t plotSys, Double_t pToAb3, Int_t paramSet, Bool_t convertToR)
{
  // Published proton yield from Phys. Rev. C 88 (2013) 044910
  // Blast wave params from π,K,p published
  // 3He/p from thermal model T = 156 MeV
 
  TGraphAsymmErrors* graph = (TGraphAsymmErrors *) generateBWpredictionsB2("PbPb276TeV", "rms", "He3", pToAb3);
  if (convertToR) convertMultiToRadius(graph, paramSet);
  
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


//---------------------- Blast wave + thermal PbPb 2.76 TeV -- Hypertriton
TGraphAsymmErrors * getBlastB3Lambda_pp7TeV(Bool_t plotSys, Double_t pToAb3, Int_t paramSet)
{
  // Published proton yield from Phys. Rev. C 88 (2013) 044910
  // Blast wave params from π,K,p published
  // 3He/p from thermal model T = 156 MeV
  // s3 from hyper-triton paper
 
  TGraphAsymmErrors* graph = (TGraphAsymmErrors *) generateBWpredictionsB2("pp7TeV", "rms", "hyper-triton", pToAb3);
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
//---------------------- Blast wave + thermal PbPb 2.76 TeV
TGraphAsymmErrors * getBlastB4_PbPb276TeV(Bool_t plotSys, Double_t pToAb4, Int_t paramSet)
{
  // Published proton yield from Phys. Rev. C 88 (2013) 044910
  // Blast wave params from π,K,p published
  // 4He from thermal model T = 156 MeV yield = 7e-7
 
  TGraphAsymmErrors* graph = (TGraphAsymmErrors *) generateBWpredictionsB2("PbPb276TeV", "rms", "He4", pToAb4);
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

//---------------------------------------------------------
//---------------------- Blast wave + thermal pp 7 TeV
TGraphAsymmErrors * getBlastB4_pp7TeV(Bool_t plotSys, Double_t pToAb4, Int_t paramSet)
{
  // Final proton yield from long paper pp7 TeV
  // Blast wave params from π,K,p published
  // 4He yield from thermal model T = 156 MeV
 
  TGraphAsymmErrors* graph = (TGraphAsymmErrors *) generateBWpredictionsB2("pp7TeV", "rms", "He4", pToAb4);
  convertMultiToRadius(graph, paramSet);
  
  graph->SetMarkerColor(kMagenta-2);
  graph->SetLineColor(kMagenta-2);
  graph->SetFillColorAlpha(kMagenta-2, 0.1);  
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.5);
  graph->SetMarkerStyle(33);
  graph->SetLineWidth(2);
  graph->SetLineStyle(2);
  return graph;
  
}


//---------------------- Blast wave + thermal PbPb 2.76 TeV -- 4LH
TGraphAsymmErrors * getBlastB4Lambda_PbPb276TeV(Bool_t plotSys, Double_t pToAb3, Int_t paramSet)
{
  // Published proton yield from Phys. Rev. C 88 (2013) 044910
  // Blast wave params from π,K,p published
  // 4LH yield from thermal model T = 156 MeV
 
  TGraphAsymmErrors* graph = (TGraphAsymmErrors *) generateBWpredictionsB2("PbPb276TeV", "rms", "4LH", pToAb3);
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

//---------------------- Blast wave + thermal PbPb 2.76 TeV -- 4LH
TGraphAsymmErrors * getBlastB4Lambda_pp7TeV(Bool_t plotSys, Double_t pToAb4, Int_t paramSet)
{
  // Published proton yield from Phys. Rev. C 88 (2013) 044910
  // Blast wave params from π,K,p published
  // 4LH yield from thermal model T = 156 MeV
 
  TGraphAsymmErrors* graph = (TGraphAsymmErrors *) generateBWpredictionsB2("pp7TeV", "rms", "4LH", pToAb4);
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