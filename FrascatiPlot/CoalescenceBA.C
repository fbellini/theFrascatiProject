#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TH2D.h"
#include "./B2vsVolume.C"

//coalescence
void    CoalescenceBA(Double_t pToA = 0.75, TString figPath = ".");
Double_t getBAfunc(Double_t *x, Double_t *par);
Double_t getBAfromCoalescence(Double_t A, Double_t JA, Double_t mT, Double_t homogR,  Double_t objSize);
TGraphErrors * MakeBATheoryGraphCoalescence(Double_t A, Double_t JA, Double_t mT, Double_t objSize);
TGraphErrors * MakeBAcoalRelSize(Double_t A, Double_t JA, Double_t pToA, Bool_t isHyper);

//thermal 2 coalescence comparison
void           Thermal2Coalescence(Bool_t draw = 0);
TGraphErrors * Thermal2CoalescenceSeparation(Int_t particle = 0, Double_t relSysUnc = 0.10);

//pseudodata
TH1D * CoalescenceYields010(Int_t iN = 0, Bool_t useBlastWave = kFALSE);
TH1D * GetProtonIn010CentBin(TString file = "~/alice/pwglf-piKp5teV/SpectraAnalysisRun2/results/spectra/spectra-pag/Preliminaries/QM2017/Spectra_PbPbLHC15o_Combined_Histograms.root");
TH1D * GetProtonIn010CentBinBlast();


//frames and plotting
TH1D * CreateFrameBA();
TH1D * CreateFrameSeparation();
TH1D * CreateFrameUncertainty();
void   SetCanvasStyle();

//make coalescence prediction
Double_t nucleiA[8] = {2, 3, 3, 3, 4, 4, 4, 4};
Double_t spin[8] = {1., 0.5, 0.5, 0.5, 0., 0., 1., 0};
Double_t objRadius[8] = {3.2, 2.15, 2.48, 6.8, 1.9, 2.4, 5.5, 2.4};
Color_t color[8] = {kBlack, kBlue+1, kAzure+8, kGreen+2, kRed+1, kOrange+8, kPink-3, kOrange};
Int_t line[8] = {1, 1, 1, 2, 1, 2, 5, 3};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Compare BA from coalescence and thermal model
/////////////////////////////////////////////////////////////////////////////////////////////////////////
TGraphErrors * Thermal2CoalescenceSeparation(Int_t particle, Double_t relSysUnc)
{
  //get separation between thermal and coalescence given pseudo data for nuclei spectra
  // the 1st argument, particle, has to be set as one of the following:
  // 0 = d, 1 = 3H, 2 = 3He, 4 = 3LH, 5 = 4He, 6 = 4LH
  //the 2nd argument, relSysUnc, is the relative sys uncert assumed, with value between 0 and 1.
  //
  Char_t particleName[6][15] = {"deuteron", "triton", "He3", "hyper-triton", "He4", "4LH"};
  Double_t nucleiA[]     = {   2,   3,     3,    3,    4,    4};
  Double_t spin[]        = { 1.0, 0.5,   0.5,  0.5,  0.0,  0.0};
  Double_t objRadius[]   = { 3.2, 2.15, 2.48,  6.8,  1.9,  2.4};
  Double_t pToA[]        = {0.75, 0.77, 0.77, 1.17, 0.75, 0.62};
  Color_t color[]        = {kBlack, kBlue+1, kAzure+8, kGreen+2, kRed+1, kOrange+8};
    
  TF1 * fBAcoalescence;
  TGraphAsymmErrors * gThermalBlast;
  TGraphErrors * gSeparationMaxSys;
  TGraphErrors * gRelativeUncBA;

  Double_t radius[10];
  Double_t radiuserr[10];
  Double_t separationMaxSys[10];
  Double_t separationerr[10];
  Double_t reluncert[10];
  Double_t reluncerterr[10];

  //get pesudodata for PbPb 5 TeV
  TFile * fin = TFile::Open("~/alice/nucleiB2/projectionsYR/ba_300818.root");
  TH1D * hist;
  switch (particle) {
  case 0:
    hist = (TH1D *) fin->Get("badeuteron"); //ptoa = 0.75 GeV/c
    break;
  case 1:
    hist = (TH1D *) fin->Get("batriton");  //ptoa = 0.77 GeV/c
    break;
  case 2:
    hist = (TH1D *) fin->Get("bahe3"); //ptoa = 0.77 GeV/c
    break;
  case 3:
    hist = (TH1D *) fin->Get("bahyperH3"); //ptoa = 1.17 GeV/c
    break;
  case 4:
    hist = (TH1D *) fin->Get("bahe4"); //ptoa = 0.75 GeV/c
    break;
  case 5:
    hist = (TH1D *) fin->Get("bahyperH4"); //ptoa = 0.62 GeV/c
  default:
    Printf("invalid particle selected");
    return 0x0;
  }
  
  if (!hist) Printf("WARNING::: cannot find histogram %i", particle);

  //coalescence predictions - energy indep.
  fBAcoalescence = new TF1("BAcoalescence", getBAfunc, 0., 7.0, 4);
  fBAcoalescence->SetParameters(nucleiA[particle], spin[particle], pToA[particle], objRadius[particle]);
  fBAcoalescence->SetLineColor(color[particle]);
  
  //thermal-blast prediction for 5 TeV
  gThermalBlast = getBAthermalBlast("PbPb502TeV", particleName[particle], pToA[particle], color[particle]);
    
  for (int ip=0; ip<10; ip++) {
    
    radius[ip] = gThermalBlast->GetX()[ip];
    radiuserr[ip] = 0.0;      
    
    //calculate relative uncertainty on BA from pseudodata
    if (hist->GetBinContent(ip+1)>0)
      reluncert[ip] = hist->GetBinError(ip+1)/hist->GetBinContent(ip+1);
    else reluncert[ip] = 0.0;
    reluncerterr[ip] = 0.0;
    
    //calculate separation in nsigma
    if (reluncert[ip]<1.0) {
      Double_t totreluncertMax = TMath::Sqrt(reluncert[ip]*reluncert[ip] + relSysUnc*relSysUnc);
      Double_t totabsuncertMax = totreluncertMax * hist->GetBinContent(ip+1);
      separationMaxSys[ip] = TMath::Abs(gThermalBlast->GetY()[ip] - fBAcoalescence->Eval(radius[ip])) / totabsuncertMax;
    } else {
      separationMaxSys[ip] = 0.0;
      reluncert[ip] = 1.; 
    }
    separationerr[ip] = 0.0;
  }

  //graph separation
  gSeparationMaxSys = new TGraphErrors(10, radius, separationMaxSys, radiuserr, separationerr);
  gSeparationMaxSys->SetMarkerColor(color[particle]);
  gSeparationMaxSys->SetMarkerStyle(20+particle);
  gSeparationMaxSys->SetLineColor(color[particle]);
  gSeparationMaxSys->SetLineStyle(2);
  gSeparationMaxSys->SetLineWidth(2);

  //graph relativeUncBA
  gRelativeUncBA = new TGraphErrors(10, radius, reluncert, radiuserr, reluncerterr);
  gRelativeUncBA->SetLineColor(color[particle]);
  gRelativeUncBA->SetMarkerColor(color[particle]);
  gRelativeUncBA->SetMarkerStyle(20+particle);
  gRelativeUncBA->SetLineWidth(2);
    
  return gSeparationMaxSys;
}

void Thermal2Coalescence(Bool_t draw)
{
  // d, 3H, 3He, 3LH, 4He, 4LH 
  Char_t particle[6][15] = {"deuteron", "triton", "He3", "hyper-triton", "He4", "4LH"};
  Double_t nucleiA[]     = {   2,   3,     3,    3,    4,    4};
  Double_t spin[]        = { 1.0, 0.5,   0.5,  0.5,  0.0,  0.0};
  Double_t objRadius[]   = { 3.2, 2.15, 2.48,  6.8,  1.9,  2.4};
  Double_t pToA[]        = {0.75, 0.77, 0.77, 1.17, 0.75, 0.62};
  Color_t color[]        = {kBlack, kBlue+1, kAzure+8, kGreen+2, kRed+1, kOrange+8};
  
  TCanvas * cr4 = 0x0;  
  if (draw) {
    cr4 = new TCanvas("cr4", "compare thermal with coalescence", 900, 2000);
    SetCanvasStyle();
    TH1D * frame = CreateFrameBA();      
    frame->GetYaxis()->SetTitleOffset(0.8);
    frame->GetYaxis()->SetRangeUser(1e-12, 1);
    frame->GetXaxis()->SetRangeUser(0.1, 6.5);
    
    TH1D * frameSep = CreateFrameSeparation();
    TH1D * frameUnc = CreateFrameUncertainty();
    cr4->Divide(1,3);
    cr4->cd(1); gPad->SetLogy(); gPad->SetTicky();  gPad->SetTickx();
    frame->Draw();
    cr4->cd(2); gPad->SetLogy(); gPad->SetTicky();  gPad->SetTickx();
    frameUnc->Draw();
    cr4->cd(3); gPad->SetLogy(); gPad->SetTicky();  gPad->SetTickx();
    frameSep->Draw();
  }
  
  TF1 * fBAcoalescence[6];
  TGraphAsymmErrors * gThermalBlast[6];
  TGraphErrors * gSeparation[6];
  TGraphErrors * gSeparationMaxSys[6];
  TGraphErrors * gRelativeUncBA[6];

  Double_t radius[6][10];
  Double_t radiuserr[6][10];
  Double_t separation[6][10];
  Double_t separationMaxSys[6][10];
  Double_t separationerr[6][10];
  Double_t reluncert[6][10];
  Double_t reluncerterr[6][10];

  //get pesudodata for PbPb 5 TeV
  TFile * fin = TFile::Open("~/alice/nucleiB2/projectionsYR/ba_300818.root");
  TH1D * hist[6];
  hist[0] = (TH1D *) fin->Get("badeuteron"); //ptoa = 0.75 GeV/c
  hist[1] = (TH1D *) fin->Get("batriton");  //ptoa = 0.77 GeV/c
  hist[2] = (TH1D *) fin->Get("bahe3"); //ptoa = 0.77 GeV/c
  hist[3] = (TH1D *) fin->Get("bahyperH3"); //ptoa = 1.17 GeV/c
  hist[4] = (TH1D *) fin->Get("bahe4"); //ptoa = 0.75 GeV/c
  hist[5] = (TH1D *) fin->Get("bahyperH4"); //ptoa = 0.62 GeV/c
  
  for (int j=0; j<6; j++) {
    if (!hist[j]) Printf("WARNING::: cannot find histogram %i", j);

    //coalescence predictions - energy indep.
    fBAcoalescence[j] = new TF1("BAcoalescence", getBAfunc, 0., 7.0, 4);
    fBAcoalescence[j]->SetParameters(nucleiA[j], spin[j], pToA[j], objRadius[j]);
    fBAcoalescence[j]->SetLineColor(color[j]);

    //thermal-blast prediction for 5 TeV
    gThermalBlast[j] = getBAthermalBlast("PbPb502TeV", particle[j], pToA[j], color[j]);
    
    for (int ip=0; ip<10; ip++) {
      
      radius[j][ip] = gThermalBlast[j]->GetX()[ip];
      radiuserr[j][ip] = 0.0;      
      //      if (j == 0) Printf("ip = %i     radius = %3.2f", ip, radius[j][ip]);
      
      //calculate relative uncertainty on BA from pseudodata
      if (hist[j]->GetBinContent(ip+1)>0)
	reluncert[j][ip] = hist[j]->GetBinError(ip+1)/hist[j]->GetBinContent(ip+1);
      else reluncert[j][ip] = 1.0;
      reluncerterr[j][ip] = 0.0;
      
      //calculate separation in nsigma
      if (reluncert[j][ip]<1.0) {
	Double_t sysreluncertMin = 0.15;
	Double_t sysreluncertMax = 0.20;
	Double_t totreluncertMin = TMath::Sqrt(reluncert[j][ip]*reluncert[j][ip] + sysreluncertMin*sysreluncertMin);
	Double_t totreluncertMax = TMath::Sqrt(reluncert[j][ip]*reluncert[j][ip] + sysreluncertMax*sysreluncertMax);
	Double_t totabsuncertMin = sysreluncertMin * hist[j]->GetBinContent(ip+1);
	Double_t totabsuncertMax = sysreluncertMax * hist[j]->GetBinContent(ip+1);
	separation[j][ip] = TMath::Abs(gThermalBlast[j]->GetY()[ip] - fBAcoalescence[j]->Eval(radius[j][ip])) / totabsuncertMin;
	separationMaxSys[j][ip] = TMath::Abs(gThermalBlast[j]->GetY()[ip] - fBAcoalescence[j]->Eval(radius[j][ip])) / totabsuncertMax;

      } else {
	separation[j][ip] = 0.0;
	separationMaxSys[j][ip] = 0.0;
	reluncert[j][ip] = 1.; 
      }
      separationerr[j][ip] = 0.0;
    }

    //graph separation
    gSeparation[j] = new TGraphErrors(10, radius[j], separation[j], radiuserr[j], separationerr[j]);
    gSeparation[j]->SetMarkerColor(color[j]);
    gSeparation[j]->SetMarkerStyle(20+j);
    gSeparation[j]->SetLineColor(color[j]);
    gSeparation[j]->SetLineWidth(2);

    gSeparationMaxSys[j] = new TGraphErrors(10, radius[j], separationMaxSys[j], radiuserr[j], separationerr[j]);
    gSeparationMaxSys[j]->SetMarkerColor(color[j]);
    gSeparationMaxSys[j]->SetMarkerStyle(20+j);
    gSeparationMaxSys[j]->SetLineColor(color[j]);
    gSeparationMaxSys[j]->SetLineStyle(2);
    gSeparationMaxSys[j]->SetLineWidth(2);

    //graph relativeUncBA
    gRelativeUncBA[j] = new TGraphErrors(10, radius[j], reluncert[j], radiuserr[j], reluncerterr[j]);
    gRelativeUncBA[j]->SetLineColor(color[j]);
    gRelativeUncBA[j]->SetMarkerColor(color[j]);
    gRelativeUncBA[j]->SetMarkerStyle(20+j);
    gRelativeUncBA[j]->SetLineWidth(2);

    //plotting
    if (draw) {
      cr4->cd(1);
      fBAcoalescence[j]->Draw("same");
      gThermalBlast[j]->Draw("samel");
      cr4->cd(2);
      gRelativeUncBA[j]->Draw("samelp");
      cr4->cd(3);
      if (j<5)  gSeparation[j]->Draw("samel");
      // if (j<5)  gSeparationMaxSys[j]->Draw("samel");
    }
  }
  
  //  TLegend * speciesLeg = new TLegend(0.13, 0.13, 0.33, 0.4, "");
  TLegend * speciesLeg = new TLegend(0.13, 0.43, 0.3, 0.93, "");
  speciesLeg->SetFillStyle(0);
  speciesLeg->SetTextSize(0.05);
  speciesLeg->SetBorderSize(0);
  //speciesLeg->SetNColumns(2);
  speciesLeg->AddEntry(gRelativeUncBA[5], " ^{4}_{#Lambda}H", "l");
  speciesLeg->AddEntry(gRelativeUncBA[4], " ^{4}He", "l");
  speciesLeg->AddEntry(gRelativeUncBA[3], " ^{3}_{#Lambda}H", "l");
  speciesLeg->AddEntry(gRelativeUncBA[2], " ^{3}He", "l");
  speciesLeg->AddEntry(gRelativeUncBA[1], " ^{3}H", "l");
  speciesLeg->AddEntry(gRelativeUncBA[0], " d", "l");

  TLegend * modelLeg = new TLegend(0.6, 0.85, 0.95, 0.95, "");
  modelLeg->SetFillStyle(0);
  modelLeg->SetTextSize(0.05);
  modelLeg->SetBorderSize(0);
  modelLeg->AddEntry(fBAcoalescence[0], "coalescence", "l");
  modelLeg->AddEntry(gThermalBlast[0], "thermal + blast-wave", "l");

  if (draw) {
    cr4->cd(1);
    modelLeg->Draw();
    cr4->cd(2);
    speciesLeg->Draw();
  }
  return;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Extract pT spectra for different (hyper-)nuclear species from predictions of BA from coalescence
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void CoalescenceSpectra(Bool_t useBlastWave){

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  TCanvas * cs = new TCanvas("cs", "cs", 800, 800);
  cs->SetBottomMargin(0.13);
  cs->SetTopMargin(0.02);
  cs->SetLeftMargin(0.13);
  cs->SetRightMargin(0.02);

  TH1D * hFrame = new TH1D("frame","frame; #it{p}_{T} (GeV/#it{c}); d^{2}#it{N}/(d#it{y}d#it{p}_{T}) (GeV/#it{c})^{-1}", 150, 0., 15.);
  hFrame->GetYaxis()->SetRangeUser(1E-14, .3);

  TH1D * hProSpectrum = (TH1D *) GetProtonIn010CentBin();//->Clone("pro010");
  hProSpectrum->SetTitle("p");
  hProSpectrum->GetYaxis()->SetRangeUser(1E-14, 300.);
  hProSpectrum->GetYaxis()->SetTitleOffset(1.5);
  hProSpectrum->GetXaxis()->SetTitleOffset(1.2);
  hProSpectrum->GetYaxis()->SetTitle("d^{2}#it{N}/(d#it{y}d#it{p}_{T}) (GeV/#it{c})^{-1}");

  TH1D * hProSpectrumBlast = (TH1D *) GetProtonIn010CentBinBlast();//->Clone("pro010");
  hProSpectrumBlast->SetTitle("p blast-wave");
  cs->cd();
  hFrame->Draw();
  hProSpectrum->Draw("");
  hProSpectrumBlast->Draw("same");
  gPad->SetLogy(); 

  TFile * fout = new TFile(Form("coalescenceYields_%s.root", (useBlastWave?"blast":"data")), "recreate");
  for (int i = 0; i<8 ; i++) {
    TH1D * hASpectrum = (TH1D *) CoalescenceYields010(i, useBlastWave);
    fout->cd();
    hASpectrum->Write();
    cs->cd();
    hASpectrum->Draw("same");
  }
  fout->Close();
  
  TLegend * leg = (TLegend*) gPad->BuildLegend(0.6, 0.75, 0.95, 0.95, "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV (0-10%)");
  leg->SetBorderSize(0);
  leg->SetNColumns(3);
  leg->SetTextSize(0.03);
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Extract yields for different (hyper-)nuclear species from predictions of BA from coalescence
/////////////////////////////////////////////////////////////////////////////////////////////////////////
TH1D * CoalescenceYields010(Int_t iN, Bool_t useBlastWave)
{
  //multiplicity for PbPb 5.02 TeV
  //ALICE, https://arxiv.org/abs/1612.08966
  Double_t multi[2] = {(2035.*2.5 + 1850.*2.5 + 1666.*2.5 + 1505.*2.5)/10.,  (52.*2.5 + 55.*2.5 + 48.*2.5 + 44.*2.5)/10.};
  Double_t radius[2] = {-1., -1.};
  getRadiusFromParameterisation(multi, radius, 1);
  if (iN==0) Printf("Centrality 0-10%% PbPb 5 TeV -- multi = %5.1f    R = %3.2f fm", multi[0], radius[0]);
  
  //GET PROTON SPECTRUM from data
  TH1D * hProData = (TH1D *) GetProtonIn010CentBin();
  const Int_t nBinsIn = hProData->GetNbinsX();
  hProData->SetLineColor(kGray+1);
  hProData->SetMarkerColor(kGray+1);
  hProData->SetMarkerStyle(28);

  //GET PROTON SPECTRUM from blast-wave
  TH1D * hProBlast = (TH1D *) GetProtonIn010CentBinBlast();
  hProBlast->SetLineColor(kGray+1);
  hProBlast->SetMarkerColor(kGray+1);
  hProBlast->SetMarkerStyle(34);

  //chose PROTON SPECTRUM to use
  TH1D * hProSpectrum = 0x0;
  if (!useBlastWave) hProSpectrum = hProData;
  else  hProSpectrum = hProBlast;
    
  //GET ARRAY OF PT FOR NUCLEUS
  Double_t ptBinsOut[nBinsIn];
  for (Int_t ibin = 1; ibin < nBinsIn+1; ibin++){
    ptBinsOut[ibin-1] = nucleiA[iN] * hProSpectrum->GetXaxis()->GetBinLowEdge(ibin);
  }

  //DEFINE NUCLEUS SPECTRUM HISTO
  TString name = "d";
  if (iN == 1) name = "3H";
  if (iN == 2) name = "3He";
  if (iN == 3) name = "3LH";
  if (iN == 4) name = "4He";
  if (iN == 5) name = "4LH";
  if (iN == 6) name = "4LLH";
  if (iN == 7) name = "4LHe";
  
  TH1D * hASpectrum = new TH1D("hASpectrum", "spectrum; #it{p}_{T} (GeV/#it{c}); d#it{N}/(d#it{y}d#it{p}_{T})", nBinsIn-1, ptBinsOut);
  hASpectrum->SetLineColor(color[iN]);
  hASpectrum->SetMarkerColor(color[iN]);
  hASpectrum->SetMarkerStyle(20+iN);
  
  for (Int_t ibin = 1; ibin < hASpectrum->GetNbinsX()+1; ibin++){ 
    Double_t proYield[2];
    Double_t pToA = hProSpectrum->GetXaxis()->GetBinCenter(ibin);
    
    //get yield d^2N/dpTdy and transform to invariant yield
    proYield[0] =  hProSpectrum->GetBinContent(ibin) / (2*TMath::Pi()*pToA);
    proYield[1] =  hProSpectrum->GetBinError(ibin)   / (2*TMath::Pi()*pToA);

    Double_t Ayield[2];
    //Get BA from formula
    Double_t BA = getBAfromCoalescence(nucleiA[iN], spin[iN], pToA, radius[0], objRadius[iN]);

    //get back from invariant yield to d^2N/dpTdy 
    if (proYield[0]>0){
      Ayield[0] = BA * TMath::Power(proYield[0], nucleiA[iN]) * (2*TMath::Pi()*nucleiA[iN]*pToA);
      Ayield[1] = nucleiA[iN] * BA * TMath::Power(proYield[0], nucleiA[iN]-1) * Ayield[0] * (2*TMath::Pi()*nucleiA[iN]*pToA);
      hASpectrum->SetBinContent(ibin, Ayield[0]);
      hASpectrum->SetBinError(ibin, Ayield[1]);
    }
  }

  /*
  TCanvas * c1 = new TCanvas("c1", "c1", 800, 800);
  c1->cd();
  //gPad->SetLogx();
  gPad->SetLogy(); 
  hProSpectrum->Draw("");
  hASpectrum->Draw("same");
  */
  
  hASpectrum->SetName(Form("hSpectrum_%s_%s",name.Data(), (useBlastWave?"blast":"")));
  hASpectrum->SetTitle(Form("%s",name.Data()));
  return hASpectrum;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
TH1D * GetProtonIn010CentBin(TString file)
{
  //proton spectra PbPb 5.02 TeV -- preliminary QM 2017
  //from analysis repository on git
  TFile * fin = TFile::Open(file.Data());
  TList * lin = (TList*) fin->Get("Summed_Proton_Sys");
  TH1D * hin05 = (TH1D*) lin->FindObject("hSpectraSummedProton_PbPb_Combined_0.00to5.00");
  TH1D * hin510 = (TH1D*) lin->FindObject("hSpectraSummedProton_PbPb_Combined_5.00to10.00");
    

  TH1D * hin010 = (TH1D*) hin05->Clone("hSpectraProton_PbPb_010");  
  hin010->Add(hin510, 1.);
  hin010->Scale(0.5); // for the average
  hin010->Scale(0.5); //scale for p+pbar
  /*
  TCanvas * c1 = new TCanvas("c1", "c1", 800, 800);
  c1->cd();
  gPad->SetLogy();
  hin05->Draw();
  hin510->Draw("same");
  hin010->Draw("same");
  */
  return hin010;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
TH1D * GetProtonIn010CentBinBlast()
{
  //use blast-wave from π,K,p to regenerate spectrum for p
  TH1D * hin05 = (TH1D *) generateBWpredictionSpectra(0, "PbPb502TeV", "proton", "rms");
  TH1D * hin510 = (TH1D *) generateBWpredictionSpectra(1, "PbPb502TeV",  "proton", "rms");
  TH1D * hin010 = (TH1D*) hin05->Clone("hSpectraProtonBlast_PbPb_010");  
  hin010->Add(hin510, 1.);
  hin010->Scale(0.5); // for the average
  /*
  TCanvas * c1 = new TCanvas("c1", "c1", 800, 800);
  c1->cd();
  gPad->SetLogy();
  hin05->Draw();
  hin510->Draw("same");
  hin010->Draw("same");
  */
  return hin010;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Make plot with coalescence BA predictions for different (hyper-)nuclear species
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void CoalescenceBA(Double_t pToA, TString figPath)
{
  //set style
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.05); 

  //make coalescence prediction
  // Double_t nucleiA[8] = {2, 3, 3, 3, 4, 4, 4, 4};
  // Double_t spin[8] = {1., 0.5, 0.5, 0.5, 0., 0., 1., 0};
  // Double_t objRadius[8] = {3.2, 2.15, 2.48, 6.8, 1.9, 2.4, 5.5, 2.4};
  // Color_t color[8] = {kBlack, kBlue+1, kAzure+8, kGreen+2, kRed+1, kOrange+8, kPink-3, kOrange};
  // Int_t line[8] = {1, 1, 1, 2, 1, 2, 5, 3};
  TGraphErrors * hBA_coalescence[8];

  for(Int_t j=0; j<8; j++) {
    hBA_coalescence[j] = (TGraphErrors*) MakeBATheoryGraphCoalescence(nucleiA[j], spin[j], pToA, objRadius[j]);
    hBA_coalescence[j]->SetMarkerStyle(1);
    hBA_coalescence[j]->SetLineWidth(3);
    hBA_coalescence[j]->SetLineStyle(line[j]);
    hBA_coalescence[j]->SetLineColor(color[j]);
  }

  TH1D * hframe4 = CreateFrameBA();
  hframe4->GetYaxis()->SetTitleOffset(1.4);
  hframe4->GetYaxis()->SetTitleSize(0.06);
  hframe4->GetYaxis()->SetRangeUser(1e-13, .1);
  //Define pT/A labels only once
  TPaveText * pavept = new TPaveText(0.47, 0.87, 0.87, 0.92, "NDC");
  pavept->SetFillStyle(0);
  pavept->SetTextFont(42);
  pavept->SetBorderSize(0);
  pavept->SetTextSize(0.03);
  pavept->SetTextAlign(12);
  pavept->AddText(Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToA));

  TPaveText * paveA[3];
  for (Int_t j = 2; j<5; j++) {
      paveA[j-2] = new TPaveText(0.8, 0.82-0.21*(j-2), 0.95, 0.82-0.27*(j-2), "NDC");
      paveA[j-2]->SetFillStyle(0);
      paveA[j-2]->SetTextFont(42);
      paveA[j-2]->SetBorderSize(0);
      paveA[j-2]->SetTextSize(0.04);
      paveA[j-2]->SetTextAlign(12);
      paveA[j-2]->AddText(Form("#it{A} = %i", j));
  }
  TLegend * masterLeg = new TLegend(0.02, 0.5, 0.9, 0.95, Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToA));
  masterLeg->SetFillStyle(0);
  masterLeg->SetTextSize(0.1);
  masterLeg->SetBorderSize(0);
  masterLeg->SetNColumns(1);
  masterLeg->AddEntry(hBA_coalescence[0], "d, #it{r} = 3.2 fm", "l");
  masterLeg->AddEntry(hBA_coalescence[1], "^{3}H, #it{r} = 2.15 fm", "l");
  masterLeg->AddEntry(hBA_coalescence[2], "^{3}He, #it{r} = 2.48 fm", "l");
  masterLeg->AddEntry(hBA_coalescence[3], "^{3} _{#Lambda}He, #it{r} = 6.8 fm", "l");
  masterLeg->AddEntry(hBA_coalescence[4], "^{4}He, #it{r} = 1.9 fm", "l");
  masterLeg->AddEntry(hBA_coalescence[5], "^{4} _{#Lambda}H, #it{r} = 2.4 fm", "l");
  masterLeg->AddEntry(hBA_coalescence[6], "^{4} _{#Lambda#Lambda}H, #it{r} = 5.5 fm", "l");
  masterLeg->AddEntry(hBA_coalescence[7], "^{4} _{#Lambda}He, #it{r} = 2.4 fm", "l");



  TCanvas * cr4 = new TCanvas("cr4", "coalescence", 900, 700);
  cr4->cd();
  TPad * pad1 = new TPad("pad1","This is pad1",0.001, 0.001, 0.75,0.999);
  pad1->SetFillColor(0);
  pad1->SetBorderMode(0);
  pad1->SetBorderSize(0);
  pad1->SetMargin(0.2,0.05,0.15,0.05);
  pad1->SetLogy();
  pad1->SetTicky();
  pad1->SetTickx();

  cr4->cd();
  TPad * pad2 = new TPad("pad2","This is pad2",0.75, 0.001, 0.999,0.9999);
  pad2->SetFillColor(0);
  pad2->SetBorderMode(0);
  pad2->SetBorderSize(0);
  pad2->SetMargin(0.05,0.05,0.15,0.05);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  hframe4->Draw();

for(Int_t j=0; j<8; j++) {
   hBA_coalescence[j]->Draw("samel");
   if (j<3) paveA[j]->Draw();
 }
 pad2->cd();
  //pavept->Draw();
 masterLeg->Draw();
 cr4->SaveAs(Form("%s/coalescenceBA%03.0f.pdf", figPath.Data(), pToA*100));
 cr4->SaveAs(Form("%s/coalescenceBA%03.0f.eps", figPath.Data(), pToA*100));
 
 return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
TGraphErrors * MakeBATheoryGraphCoalescence(Double_t A, Double_t JA, Double_t mT, Double_t objSize)
{
  const Int_t nPoints = 1000;
  Double_t gY[nPoints];
  Double_t gR[nPoints];
  
  TGraphErrors * graphOut = new TGraphErrors(nPoints);
  graphOut->SetName("B4_th_coalescence");
  graphOut->SetTitle(Form("B_{%1.0f} from coalescence", A));
  
  
  for (int i = 0; i<nPoints; i++){
    gR[i] = 10.0 * i / nPoints; 
    gY[i] = getBAfromCoalescence(A, JA, mT, gR[i], objSize);

    graphOut->SetPoint(i, gR[i], gY[i]);
    graphOut->SetPointError(i, 0.0, 0.0);
  }

  return graphOut;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t getBAfromCoalescence(Double_t A, Double_t JA, Double_t pToA, Double_t homogR,  Double_t objSize)
{
  // formula 13 of the Frascati paper arXiv:1807.05894
  Double_t convFactor_fm2InvGeV = 0.197;
  Double_t spinFactor = (2*JA+1) / TMath::Power(2., A);
  Double_t expo = 3.*(A-1)/2.;
  Double_t mT = TMath::Sqrt(0.938*0.938 + pToA*pToA);
  Double_t mtfactor = TMath::Power(mT, -(A-1));
  Double_t BA =  spinFactor * 1./TMath::Sqrt(A) * mtfactor *
    TMath::Power(2* TMath::Pi(), expo) *
    TMath::Power((homogR * homogR +  objSize * objSize / 4.)/(convFactor_fm2InvGeV * convFactor_fm2InvGeV), -expo);
  
  return BA;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t getBAfunc(Double_t *x, Double_t *par)
{
  //parameters are
  /*
    par[0] = A;
    par[1] = JA;
    par[2] = pT/A;
    par[3] = object size;
  */
  return getBAfromCoalescence(par[0], par[1], par[2], x[0], par[3]);
    
}

TH1D * CreateFrameSeparation()
{
  TH1D * hframe4 = new TH1D("hframe", "separation; #it{R} (fm); | #it{B}_{A}^{ther}-#it{B}_{A}^{coal} | / #sigma_{stat}^{pseudo}", 3500, 0.01, 7.01);
  hframe4->GetXaxis()->SetTitleSize(0.06);
  hframe4->GetYaxis()->SetTitleSize(0.06);
  hframe4->GetYaxis()->SetTitleOffset(0.8);
  hframe4->GetXaxis()->SetTitleOffset(0.8);
  hframe4->GetXaxis()->SetLabelSize(0.05);
  hframe4->GetYaxis()->SetLabelSize(0.05);
  hframe4->GetXaxis()->SetRangeUser(0.01, 7.);
  hframe4->GetYaxis()->SetRangeUser(1., 2e2);
  return hframe4;
}

TH1D * CreateFrameUncertainty()
{
  TH1D * hframe4 = new TH1D("hframe", "relative uncertainty; #it{R} (fm); #sigma_{stat}^{pseudo}/#it{B}_{A}^{pseudo}", 3500, 0.01, 7.01);
  hframe4->GetXaxis()->SetTitleSize(0.06);
  hframe4->GetYaxis()->SetTitleSize(0.06);
  hframe4->GetYaxis()->SetTitleOffset(0.8);
  hframe4->GetXaxis()->SetTitleOffset(0.8);
  hframe4->GetXaxis()->SetLabelSize(0.05);
  hframe4->GetYaxis()->SetLabelSize(0.05);
  hframe4->GetXaxis()->SetRangeUser(0.01, 7.);
  hframe4->GetYaxis()->SetRangeUser(0.0001, 1.2);
  return hframe4;
}

TH1D * CreateFrameBA()
{
  TH1D * hframeBA = new TH1D("hframeBA", "B_{A} vs radius; #it{R} (fm); #it{B}_{A}   [(GeV^{2}/#it{c}^{3})^{A-1}]", 3500, 0.01, 7.01);
  hframeBA->GetXaxis()->SetTitleSize(0.06);
  hframeBA->GetYaxis()->SetTitleSize(0.06);
  hframeBA->GetYaxis()->SetTitleOffset(1.3);
  hframeBA->GetXaxis()->SetTitleOffset(0.8);
  hframeBA->GetXaxis()->SetLabelSize(0.05);
  hframeBA->GetYaxis()->SetLabelSize(0.05);
  hframeBA->GetYaxis()->SetRangeUser(2.E-12, 3.E-3);
  return hframeBA;
}

void SetCanvasStyle()
{
  //cosmetics
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadRightMargin(0.05);
  return;
}

