#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TH2D.h"
#include "./B2vsVolume.C"

void CoalescenceBA(Double_t pToA = 0.75);
TH1D * CoalescenceYields010(Int_t iN = 0, Bool_t useBlastWave = kFALSE);
Double_t getBAfromCoalescence(Double_t A, Double_t JA, Double_t mT, Double_t homogR,  Double_t objSize);
TGraphErrors * MakeBATheoryGraphCoalescence(Double_t A, Double_t JA, Double_t mT, Double_t objSize);
TGraphErrors * MakeBAcoalRelSize(Double_t A, Double_t JA, Double_t pToA, Bool_t isHyper);
TH2D * CreateFrameBA();
TH1D * GetProtonIn010CentBin(TString file = "~/alice/pwglf-piKp5teV/SpectraAnalysisRun2/results/spectra/spectra-pag/Preliminaries/QM2017/Spectra_PbPbLHC15o_Combined_Histograms.root");
TH1D * GetProtonIn010CentBinBlast();

//make coalescence prediction
Double_t nucleiA[8] = {2, 3, 3, 3, 4, 4, 4, 4};
Double_t spin[8] = {1., 0.5, 0.5, 0.5, 0., 0., 1., 0};
Double_t objRadius[8] = {3.2, 2.15, 2.48, 6.8, 1.9, 2.4, 5.5, 2.4};
Color_t color[8] = {kBlack, kBlue+1, kAzure+8, kGreen+2, kRed+1, kOrange+8, kPink-3, kOrange};
Int_t line[8] = {1, 1, 1, 2, 1, 2, 5, 3};



/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Extract pT spectra for different (hyper-)nuclear species from predictions of BA from coalescence
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void CoalescenceSpectra(Bool_t useBlastWave){

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  TFile * fout = new TFile(Form("coalescenceYields_%s.root", (useBlastWave?"blast":"data")), "recreate");
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

  for (int i = 0; i<8 ; i++) {
    TH1D * hASpectrum = (TH1D *) CoalescenceYields010(i, useBlastWave);
    fout->cd();
    hASpectrum->Write();
    cs->cd();
    hASpectrum->Draw("same");
  }
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
  if (iN==0) Printf("Centrality 0-10% PbPb 5 TeV -- multi = %5.1f    R = %3.2f fm", multi[0], radius[0]);
  
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
  TFile * fin = TFile::Open(file.Data);
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
void CoalescenceBA(Double_t pToA = 0.75)
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
    hBA_coalescence[j]->SetLineWidth(2);
    hBA_coalescence[j]->SetLineStyle(line[j]);
    hBA_coalescence[j]->SetLineColor(color[j]);
  }

  TH2D * hframe4 = CreateFrameBA();
  
  //Define pT/A labels only once
  TPaveText * pavept = new TPaveText(0.47, 0.87, 0.87, 0.92, "NDC");
  pavept->SetFillStyle(0);
  pavept->SetTextFont(42);
  pavept->SetBorderSize(0);
  pavept->SetTextSize(0.04);
  pavept->SetTextAlign(12);
  pavept->AddText(Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToA));

  TLegend * masterLeg = new TLegend(0.1, 0.3, 0.5, 0.9, "");
  masterLeg->SetFillStyle(0);
  masterLeg->SetTextSize(0.05);
  masterLeg->SetBorderSize(0);
  masterLeg->AddEntry(hBA_coalescence[0], "d, #it{r} = 3.2 fm", "l");
  masterLeg->AddEntry(hBA_coalescence[1], "^{3}H, #it{r} = 2.15 fm", "l");
  masterLeg->AddEntry(hBA_coalescence[2], "^{3}He, #it{r} = 2.48 fm", "l");
  masterLeg->AddEntry(hBA_coalescence[3], "^{3} _{#Lambda}He, #it{r} = 6.8 fm", "l");
  masterLeg->AddEntry(hBA_coalescence[4], "^{4}He, #it{r} = 1.9 fm", "l");
  masterLeg->AddEntry(hBA_coalescence[5], "^{4} _{#Lambda}H, #it{r} = 2.4 fm", "l");
  masterLeg->AddEntry(hBA_coalescence[6], "^{4} _{#Lambda#Lambda}H, #it{r} = 5.5 fm", "l");
  masterLeg->AddEntry(hBA_coalescence[7], "^{4} _{#Lambda}He, #it{r} = 2.4 fm", "l");

  TCanvas * cr4 = new TCanvas("cr4", "compare thermal with coalescence", 1600, 900);
  cr4->Divide(2,1);
  cr4->cd(1);
  gPad->SetLogy();
  gPad->SetTicky();
  gPad->SetTickx();
  hframe4->Draw();
 for(Int_t j=0; j<8; j++) {
   hBA_coalescence[j]->Draw("samel");
 }
 pavept->Draw();
 cr4->cd(2);
 masterLeg->Draw();
 cr4->SaveAs(Form("coalescenceBA_%3.2f.pdf", pToA));
 
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
TH2D * CreateFrameBA()
{
  TH2D * hframe4 = new TH2D("hframe4", "B_{A} vs radius; #it{R} (fm); #it{B}_{A}", 2000, 0.01, 7.0, 2000, 1.e-12, 1.e-1);
  hframe4->GetXaxis()->SetTitleSize(0.06);
  hframe4->GetYaxis()->SetTitleSize(0.06);
  hframe4->GetYaxis()->SetTitleOffset(1.3);
  hframe4->GetXaxis()->SetTitleOffset(0.8);
  hframe4->GetXaxis()->SetLabelSize(0.05);
  hframe4->GetYaxis()->SetLabelSize(0.05);
  hframe4->GetXaxis()->SetRangeUser(0.01, 7.);

  return hframe4;
}
