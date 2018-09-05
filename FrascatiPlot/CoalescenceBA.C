#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TH2D.h"


void CoalescenceBA(Double_t pToA = 0.75);
void CoalescenceYields();
Double_t getBAfromCoalescence(Double_t A, Double_t JA, Double_t mT, Double_t homogR,  Double_t objSize);
TGraphErrors * MakeBATheoryGraphCoalescence(Double_t A, Double_t JA, Double_t mT, Double_t objSize);
TGraphErrors * MakeBAcoalRelSize(Double_t A, Double_t JA, Double_t pToA, Bool_t isHyper);
TH2D * CreateFrameBA();
TH1D * GetProtonIn010CentBin();

//make coalescence prediction
Double_t nucleiA[8] = {2, 3, 3, 3, 4, 4, 4, 4};
Double_t spin[8] = {1., 0.5, 0.5, 0.5, 0., 0., 1., 0};
Double_t objRadius[8] = {3.2, 2.15, 2.48, 6.8, 1.9, 2.4, 5.5, 2.4};
Color_t color[8] = {kBlack, kBlue+1, kAzure+8, kGreen+2, kRed+1, kOrange+8, kPink-3, kOrange};
Int_t line[8] = {1, 1, 1, 2, 1, 2, 5, 3};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Extract yields for different (hyper-)nuclear species from predictions of BA from coalescence
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void CoalescenceYields()
{
  Double_t radiusCent = 1.;

  //GetProton Spectrum
  TH1D * hProSpectrum = (TH1D *) GetProtonIn010CentBin();
  const Int_t nBinsIn = hProSpectrum->GetBinsX();
  Double_t ptBinsOut[nBinsIn];
  
  for (Int_t ibin = 1; ibin < nBins; ibin++){
    Double_t proPtBin[4];
    proPtBin[0] = hProSpectrum->GetXaxis()->GetBinCenter(ibin); //center
    proPtBin[1] = hProSpectrum->GetXaxis()->GetBinLowEdge(ibin); //low edge
    proPtBin[2] = hProSpectrum->GetXaxis()->GetBinUpEdge(ibin); //up edge
    proPtBin[3] = proPtBinUp - proPtBinLow; //bin width
    
    proYield[0] =  hProSpectrum->GetXaxis()->GetBinContent(ibin); 
    proYield[1] =  hProSpectrum->GetXaxis()->GetBinError(ibin); 
    
    Double_t APtBin[4];
    for(Int_t j = 0; j < 4; j++) {   
      APtBin[j] = A * proPtBin[j];
    }
    
    Double_t BA[8] = {0,0,0,0,0,0,0,0};
    Double_t Ayield[8];
    Double_t pToA = proPtBin[0];
    
    for(Int_t j=0; j<8; j++) {   
      BA[j] = getBAfromCoalescence(nucleiA[j], spin[j], pToA, radiusCent, objRadius[j]);
    Ayield[j] = BA[j] * TMath::Power(proYield[0], nucleiA[j]);
    }
  }
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
TH1D * GetProtonIn010CentBin()
{
  TFile * fin = TFile::Open("~/alice/pwglf-piKp5teV/SpectraAnalysisRun2/results/spectra/spectra-pag/Preliminaries/QM2017/Spectra_PbPbLHC15o_Combined_Histograms.root");
  TList * lin = (TList*) fin->Get("Summed_Proton_Sys");
  TH1D * hin05 = (TH1D*) lin->FindObject("hSpectraSummedProton_PbPb_Combined_0.00to5.00");
  TH1D * hin510 = (TH1D*) lin->FindObject("hSpectraSummedProton_PbPb_Combined_5.00to10.00");

  TH1D * hin010 = hin05->Clone("hSpectraProton_PbPb_Combined_0.00to10.00");
  hin010->Add("hSpectraSummedProton_PbPb_Combined_5.00to10.00");
  hin010->Scale(0.5);
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
