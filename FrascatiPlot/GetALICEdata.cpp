/////////////////////////////////////////////////////////////////////////////////////////////////////////
//The Frascati project - ALICE data - colaescence parameters
/////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "TGraphErrors.h"
#include "./MapMulti2R.C"

//deuteron
TGraphErrors * getB2_pp7TeV(Bool_t plotSys = 0, Double_t pToA = 0.75, Int_t paramSet = 0, Bool_t convertToR = kTRUE);
TGraphErrors * getB2_pp13TeV(Bool_t plotSys = 0, Double_t pToA = 0.75, Int_t paramSet = 0, Bool_t convertToR = kTRUE);
TGraphErrors * getB2_pp7TeVINELg0(Bool_t plotSys = 0, Double_t pToA = 0.75, Int_t paramSet = 0, Bool_t convertToR = kTRUE);
TGraphErrors * getB2_pPb5TeV(Bool_t plotSys = 0, Double_t pToA = 0.75, Int_t paramSet = 0, Bool_t convertToR = kTRUE);
TGraphErrors * getB2_PbPb5TeV(Bool_t plotSys = 0, Double_t pToA = 0.75, Int_t paramSet = 0, Bool_t convertToR = kTRUE);
TGraphErrors * getB2_PbPb276TeV(Bool_t plotSys = 0, Double_t pToA = 0.75, Int_t paramSet = 0, Bool_t convertToR = kTRUE);

//helium-3
TGraphErrors      * getB3_PbPb5TeV(Bool_t plotSys = 0, Double_t pToAb3 = 0.733, Int_t paramSet = 0, Bool_t convertToR = kTRUE);
TGraphErrors      * getB3_PbPb276TeV(Bool_t plotSys = 0, Double_t pToAb3 = 0.733, Int_t paramSet = 0, Bool_t convertToR = kTRUE);
TGraphErrors      * getB3_pPb5TeV(Bool_t plotSys = 0, Double_t pToAb3 = 0.733, Int_t paramSet = 0, Bool_t convertToR = kTRUE);
TGraphAsymmErrors * getB3_pp7TeVINELg0(Bool_t plotSys = 0, Double_t pToAb3pp = 0.800, Int_t paramSet = 0, Bool_t convertToR = kTRUE);

//hypertriton
TGraphAsymmErrors * getB3Lambda_PbPb276TeV(Bool_t plotSys = 0, Double_t pToAb3Lambda = 1., Int_t paramSet = 0);

//---------------------------------------------------------
//------------------------------ ALICE data B2 --- pp 7 TeV
TGraphErrors * getB2_pp7TeV(Bool_t plotSys, Double_t pToA, Int_t paramSet, Bool_t convertToR)
{
  //from Manuel - pp 7 TeV EXA 2017 preliminary
  TFile * f0 = TFile::Open("/Users/fbellini/alice/nucleiB2/data/oB2vsNch.root");
  if (!f0) return NULL;

  TString gName = Form("B2_pToA=%0.3f%s", pToA, (plotSys? "_sys" : ""));
  TGraphErrors * graph = (TGraphErrors*) f0->Get(gName.Data());
  if (!graph) {
    Printf("Error: cannot retrieve graph. Check pt bin requested.");
    return NULL;
  }

  if (convertToR) convertMultiToRadius(graph, paramSet);

  graph->SetMarkerColor(kGreen+2);
  graph->SetLineColor(kGreen+2);
  graph->SetFillColorAlpha(kGreen+2, 0.1);
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.2);
  graph->SetMarkerStyle(20);
  return graph;
  
}

//---------------------------------------------------------
//------------------------------ ALICE data B2 --- pp 13 TeV
TGraphErrors * getB2_pp13TeV(Bool_t plotSys, Double_t pToA, Int_t paramSet, Bool_t convertToR)
{
  //from Luca B. - pp 13 TeV QM 2018 preliminary
  TFile * f0 = TFile::Open("/Users/fbellini/alice/nucleiB2/data/customB2PtFixed75_A.root");
  if (!f0) return NULL;

  TString gName = Form("pp13TeV%s", (plotSys? "Syst" : "Stat"));
  TGraphErrors * graph = (TGraphErrors*) f0->Get(gName.Data());
  if (!graph) {
    Printf("Error: cannot retrieve graph. Check pt bin requested.");
    return NULL;
  }

  if (convertToR) convertMultiToRadius(graph, paramSet);
  graph->SetMarkerColor(kOrange);
  graph->SetLineColor(kOrange);
  graph->SetFillColorAlpha(kOrange, 0.1);
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.6);
  graph->SetMarkerStyle(33);
  return graph;
  
}

//---------------------------------------------------------
//------------------------------ ALICE data B2 --- pp 7 TeV INEL > 0
TGraphErrors * getB2_pp7TeVINELg0(Bool_t plotSys, Double_t pToA, Int_t paramSet, Bool_t convertToR)
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

  if (convertToR) convertMultiToRadius(graph, paramSet);

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
TGraphErrors * getB2_pPb5TeV(Bool_t plotSys, Double_t pToA, Int_t paramSet, Bool_t convertToR)
{

  TFile * f0 = TFile::Open("/Users/fbellini/alice/nucleiB2/data/oB2vsNch.root");
  if (!f0) return NULL;

  TString gName = Form("B2_pPb_pToA=%4.3f%s", pToA, (plotSys? "_sys" : ""));
  TGraphErrors * graph = (TGraphErrors*) f0->Get(gName.Data());
  if (!graph) {
    Printf("Error: cannot retrieve graph. Check pt bin requested.");
    return NULL;
  }

  if (convertToR) convertMultiToRadius(graph, paramSet);

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
TGraphErrors * getB2_PbPb5TeV(Bool_t plotSys, Double_t pToA, Int_t paramSet, Bool_t convertToR)
{
  //Preliminary from Max - QM 2017
  TFile * f0 = TFile::Open("/Users/fbellini/alice/nucleiB2/data/oB2vsNch.root");
  if (!f0) return NULL;

  TString gName = Form("B2_PbPb15_pToA=%4.3f%s", pToA, (plotSys? "_sys" : ""));
  TGraphErrors * graph = (TGraphErrors*) f0->Get(gName.Data());
  if (!graph) {
    Printf("Error: cannot retrieve graph. Check pt bin requested.");
    return NULL;
  }

  if (convertToR) convertMultiToRadius(graph, paramSet);

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
TGraphErrors * getB2_PbPb276TeV(Bool_t plotSys, Double_t pToA, Int_t paramSet, Bool_t convertToR)
{
  //Published PRC 93, 0249717 (2016)
  TFile * f0 = TFile::Open("/Users/fbellini/alice/nucleiB2/data/oB2vsNch.root");
  if (!f0) return NULL;

  TString gName = Form("B2_PbPb10_pToA=%4.3f%s", pToA, (plotSys? "_sys" : ""));
  TGraphErrors * graph = (TGraphErrors*) f0->Get(gName.Data());
  if (!graph) {
    Printf("Error: cannot retrieve graph. Check pt bin requested.");
    return NULL;
  }

  if (convertToR) convertMultiToRadius(graph, paramSet);

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

TGraphAsymmErrors * getB3_pp7TeVINELg0(Bool_t plotSys, Double_t pToAb3pp, Int_t paramSet, Bool_t convertToR)
{
  //from ALICE, Phys. Rev. C 97, 024615 (2018) -- arXiv:1709.08522
  // x value for multi is INEL >0
  TFile * f0 = TFile::Open("/Users/fbellini/alice/nucleiB2/data/B3pToA_pp7TeV.root");
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
    
  if (convertToR) convertMultiToRadius(graph, paramSet);

  graph->SetMarkerColor(kGreen+2);
  graph->SetLineColor(kGreen+2);
  graph->SetFillColorAlpha(kGreen-10, 0.1);
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.2);
  graph->SetMarkerStyle(20);
  return graph;
  
}


//---------------------------------------------------------
//------------------------------ ALICE data B3 --- PbPb 5 TeV
TGraphErrors * getB3_PbPb5TeV(Bool_t plotSys, Double_t pToAb3, Int_t paramSet, Bool_t convertToR)
{
  //Preliminary from Max - QM 2017
  TFile * f0 = TFile::Open("/Users/fbellini/alice/nucleiB2/data/B3pToA_PbPb5TeV_preliminarySQM17.root");
  if (!f0) return NULL;

  TString gName = Form("B3_PbPb15_pToA=%4.3f%s", pToAb3, (plotSys? "_sys" : ""));
  TGraphErrors * graph = (TGraphErrors*) f0->Get(gName.Data());
  if (!graph) {
    Printf("Error: cannot retrieve graph. Check pt bin requested.");
    return NULL;
  }

  if (convertToR) convertMultiToRadius(graph, paramSet);

  graph->SetMarkerColor(kRed);
  graph->SetLineColor(kRed);
  graph->SetFillColorAlpha(kRed, 0.1);
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.2);
  graph->SetMarkerStyle(20);
  return graph;
  
}

//---------------------------------------------------------
//------------------------------ ALICE data B3 p-Pb 5 TeV
TGraphErrors * getB3_pPb5TeV(Bool_t plotSys, Double_t pToAb3, Int_t paramSet, Bool_t convertToR)
{
  //Ongoing Analysis from Sebastian Hornung
  TFile * f0 = TFile::Open("/Users/fbellini/alice/nucleiB2/data/B3pToA_pPb502TeV.root");
  if (!f0) return NULL;
  
  TString gName = Form("B3_pPb15_pToA=%4.3f%s", pToAb3, (plotSys? "_sys" : ""));
  TGraphErrors * graph = (TGraphErrors*) f0->Get(gName.Data());
  if (!graph) {
    Printf("Error: cannot retrieve graph. Check pt bin requested.");
    return NULL;
  }

  if (convertToR) convertMultiToRadius(graph, paramSet);
  
  graph->SetMarkerColor(kBlue+2);
  graph->SetLineColor(kBlue+2);
  graph->SetFillColorAlpha(kBlue+2, 0.1);
  graph->SetFillStyle(1001);
  graph->SetMarkerSize(1.2);
  graph->SetMarkerStyle(20);
  return graph;
  
}


//---------------------------------------------------------
//------------------------------ ALICE data B3 --- PbPb 2.76 TeV
TGraphErrors * getB3_PbPb276TeV(Bool_t plotSys, Double_t pToAb3, Int_t paramSet, Bool_t convertToR)
{
  //Published PRC 93, 0249717 (2016)
  TFile * f0 = TFile::Open("/Users/fbellini/alice/nucleiB2/data/B3pToA_PbPb276TeV.root");
  if (!f0) return NULL;

  TString gName = Form("B3_PbPb10_pToA=%4.3f%s", pToAb3, (plotSys? "_sys" : ""));
  TGraphErrors * graph = (TGraphErrors*) f0->Get(gName.Data());
  if (!graph) {
    Printf("Error: cannot retrieve graph. Check pt bin requested.");
    return NULL;
  }

  if (convertToR) convertMultiToRadius(graph, paramSet);

  graph->SetMarkerColor(kRed+1);
  graph->SetLineColor(kRed+1);
  graph->SetFillColorAlpha(kRed+1, 0.3);  
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
  TFile * f0 = TFile::Open("/Users/fbellini/alice/nucleiB2/data/B3LambdapToA_PbPb276TeV.root");
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