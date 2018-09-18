/*
   fbellini@cern.ch
   17.09.2018 - Yellow Report plot
*/

#include "./CoalescenceBA.C"

void ConfigTwoPanelsPad(TPad* pad1, TPad* pad2, Float_t leftmargin, Float_t rightmargin);
void ConfigThreePanelsPad(TPad* pad1, TPad* pad2, TPad* pad3, Float_t leftmargin = 0.15, Float_t rightmargin = 0.001);

//pseudodata
TGraphErrors * GetPseudoBAvsMulti(TString particle = "deuteron",  Int_t paramSet = 1);
TGraphErrors * GetPseudoBAvsMultiUnc(TString particle = "deuteron",  Int_t paramSet = 1);

//main
void MakeYRfigure(Double_t pToAb3 = 0.77, Double_t pToAb3Lambda = 1.17, Double_t pToAb4 = 0.75, Double_t pToAb4Lambda = 0.62, Bool_t plotPseudoData = 1)
{
  //
  // main function which generates Figure 1 for the LF chapter of the yellow report
  //
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.02); 

  Double_t Rmin = 0.01; Double_t Rmax = 6.5;
  //------------------------------
  // Pseudodata
  //------------------------------
  TGraphErrors * gPseudo3He = (TGraphErrors *) GetPseudoBAvsMulti("he3", 1);
  TGraphErrors * gPseudo3LambdaH = (TGraphErrors *) GetPseudoBAvsMulti("hyperH3", 1);
  TGraphErrors * gPseudo4He = (TGraphErrors *) GetPseudoBAvsMulti("he4", 1);
  TGraphErrors * gPseudo4LambdaH = (TGraphErrors *) GetPseudoBAvsMulti("hyperH4", 1);
  
  //-----------------------------
  //theory - Blast Wave + thermal
  //-----------------------------
  // pToAb3 = 0.77;
  // pToAb3Lambda = 1.17;
  // pToAb4 = 0.75;
  // pToAb4Lambda = 0.62;
  
  TGraphAsymmErrors* gBlastB3vsR_PbPb502TeV = (TGraphAsymmErrors *)  getBAthermalBlast("PbPb502TeV", "He3", pToAb3, kBlue, 1);
  TGraphAsymmErrors* gBlastB3LambdavsR_PbPb502TeV = (TGraphAsymmErrors *) getBAthermalBlast("PbPb502TeV", "hyper-triton", pToAb3Lambda, kBlue, 1);
  TGraphAsymmErrors* gBlastB4vsR_PbPb502TeV = (TGraphAsymmErrors *) getBAthermalBlast("PbPb502TeV", "He4", pToAb4, kBlue, 1);
  TGraphAsymmErrors* gBlastB4LambdavsR_PbPb502TeV = (TGraphAsymmErrors *) getBAthermalBlast("PbPb502TeV", "4LH", pToAb4Lambda, kBlue, 1);
  
  //--------------------
  //theory - coalescence
  //--------------------
  //mT is the mass of the particle relative to which the HBT radius is calculated.
  //We use now the mass of the proton and the pT per nucleon
  Double_t mT3 = TMath::Sqrt(pToAb3 * pToAb3 + 0.938 * 0.938);
  Double_t mT3L = TMath::Sqrt(pToAb3Lambda * pToAb3Lambda + 0.938 * 0.938);
  Double_t mT4 = TMath::Sqrt(pToAb4 * pToAb4 + 0.938 * 0.938);
  Double_t mT4L = TMath::Sqrt(pToAb4Lambda * pToAb4Lambda + 0.938 * 0.938);
    
  TGraphErrors* hB3_coalescence = (TGraphErrors*) MakeB3TheoryGraphCoalescence(mT3);
  hB3_coalescence->SetMarkerStyle(20);
  hB3_coalescence->SetMarkerColor(kBlack);
  hB3_coalescence->SetLineColor(kBlack);
  hB3_coalescence->SetLineWidth(3);

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
  hB3L_coalescence_largeradius->SetLineStyle(5);

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
  hB4L_coalescence_largeradius->SetLineStyle(5);
 
  //---------------------------------------
  // PLOT Yellow Report figure
  //---------------------------------------  
  TH2D * hframe = new TH2D("hframeFig4", "B_{2} vs radius; #it{R} (fm); #it{B}_{2} (GeV^{2}/#it{c}^{3})", 2000, 0.01, 7.0, 2000, 1.e-4, 0.1);
  hframe->GetXaxis()->SetTitleSize(0.06);
  hframe->GetYaxis()->SetTitleSize(0.06);
  hframe->GetYaxis()->SetTitleOffset(1.3);
  hframe->GetXaxis()->SetTitleOffset(0.8);
  hframe->GetXaxis()->SetLabelSize(0.05);
  hframe->GetYaxis()->SetLabelSize(0.05);
  hframe->GetXaxis()->SetRangeUser(Rmin, Rmax);

  TH2D * hframe3 = new TH2D("hframe3Fig4", "B_{3} vs radius; #it{R} (fm); #it{B}_{3} (GeV^{4}/#it{c}^{6})", 2000, 0.01, 7.0, 2000, 1.e-9, 1.e-1);
  hframe3->GetXaxis()->SetTitleSize(0.06);
  hframe3->GetYaxis()->SetTitleSize(0.06);
  hframe3->GetYaxis()->SetTitleOffset(1.3);
  hframe3->GetXaxis()->SetTitleOffset(0.8);
  hframe3->GetXaxis()->SetLabelSize(0.05);
  hframe3->GetYaxis()->SetLabelSize(0.05);
  hframe3->GetXaxis()->SetRangeUser(Rmin, Rmax);

  TH2D * hframe3L = new TH2D("hframe3LFig4", "B_{3,#Lambda} vs radius; #it{R} (fm); #it{B}_{3,#Lambda} (GeV^{4}/#it{c}^{6})", 2000, 0.01, 7.0, 2000, 1.e-9, 1.e-1);
  hframe3L->GetXaxis()->SetTitleSize(0.06);
  hframe3L->GetYaxis()->SetTitleSize(0.06);
  hframe3L->GetYaxis()->SetTitleOffset(1.3);
  hframe3L->GetXaxis()->SetTitleOffset(0.8);
  hframe3L->GetXaxis()->SetLabelSize(0.05);
  hframe3L->GetYaxis()->SetLabelSize(0.05);
  hframe3L->GetXaxis()->SetRangeUser(Rmin, Rmax);

 TH2D * hframe4 = new TH2D("hframe4Fig4", "B_{4} vs radius; #it{R} (fm); #it{B}_{4} (GeV^{6}/#it{c}^{9})", 2000, 0.01, 7.0, 2000, 1.e-12, 1.e-1);
  hframe4->GetXaxis()->SetTitleSize(0.06);
  hframe4->GetYaxis()->SetTitleSize(0.06);
  hframe4->GetYaxis()->SetTitleOffset(1.3);
  hframe4->GetXaxis()->SetTitleOffset(0.8);
  hframe4->GetXaxis()->SetLabelSize(0.05);
  hframe4->GetYaxis()->SetLabelSize(0.05);
  hframe4->GetXaxis()->SetRangeUser(Rmin, Rmax);

  TH2D * hframe4L = new TH2D("hframe4LFig4", "B_{3,#Lambda} vs radius; #it{R} (fm); #it{B}_{4,#Lambda} (GeV^{6}/#it{c}^{9})", 2000, 0.01, 7.0, 2000, 1.e-12, 1.e-1);
  hframe4L->GetXaxis()->SetTitleSize(0.06);
  hframe4L->GetYaxis()->SetTitleSize(0.06);
  hframe4L->GetYaxis()->SetTitleOffset(1.3);
  hframe4L->GetXaxis()->SetTitleOffset(0.8);
  hframe4L->GetXaxis()->SetLabelSize(0.05);
  hframe4L->GetYaxis()->SetLabelSize(0.05);
  hframe4L->GetXaxis()->SetRangeUser(Rmin, Rmax);
  
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

  TPaveText * paveLab4 = new TPaveText(0.77, 0.8, 0.9, 0.9, "NDC");
  paveLab4->SetFillStyle(0);
  paveLab4->SetTextFont(42);
  paveLab4->SetBorderSize(0);
  paveLab4->SetTextSize(0.1);
  paveLab4->SetTextAlign(12);
  paveLab4->AddText("#bf{^{4}He}");

   TPaveText * paveLab4L = new TPaveText(0.77, 0.8, 0.9, 0.9, "NDC");
  paveLab4L->SetFillStyle(0);
  paveLab4L->SetTextFont(42);
  paveLab4L->SetBorderSize(0);
  paveLab4L->SetTextSize(0.1);
  paveLab4L->SetTextAlign(12);
  paveLab4L->AddText("#bf{ ^{4}_{#Lambda}H}");
  
  //Define pT/A labels only once
  TPaveText * pavept = new TPaveText(0.6, 0.7, 0.9, 0.78, "NDC");
  pavept->SetFillStyle(0);
  pavept->SetTextFont(42);
  pavept->SetBorderSize(0);
  pavept->SetTextSize(0.05);
  pavept->SetTextAlign(12);
  pavept->AddText(Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb4));

  //TPaveText * paveptB3 = new TPaveText(0.55, 0.62, 0.95, 0.67, "NDC");
  TPaveText * paveptB3 = new TPaveText(0.6, 0.7, 0.9, 0.78, "NDC");
  paveptB3->SetFillStyle(0);
  paveptB3->SetTextFont(42);
  paveptB3->SetBorderSize(0);
  paveptB3->SetTextSize(0.05);
  paveptB3->SetTextAlign(12);
  paveptB3->AddText(Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb3));

  TPaveText * paveptB3L = new TPaveText(0.6, 0.7, 0.9, 0.78, "NDC");
  paveptB3L->SetFillStyle(0);
  paveptB3L->SetBorderSize(0);
  paveptB3L->SetTextFont(42);
  paveptB3L->SetTextSize(0.05);
  paveptB3L->SetTextAlign(12);
  paveptB3L->AddText(Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb3Lambda));


  TPaveText * paveptB4L = new TPaveText(0.6, 0.7, 0.9, 0.78, "NDC");
  paveptB4L->SetFillStyle(0);
  paveptB4L->SetBorderSize(0);
  paveptB4L->SetTextFont(42);
  paveptB4L->SetTextSize(0.05);
  paveptB4L->SetTextAlign(12);
  paveptB4L->AddText(Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb4Lambda));

  //-------------------------------
  //-------------------------------
  //-------------------------------
  //   DRAW FIRST OPTION -- 4 panels, pseudodata on top
  //-------------------------------
  //-------------------------------
  //-------------------------------

  TCanvas * cr4 = new TCanvas("cr4", "compare thermal with coalescence", 1000, 1200);
  cr4->SetBottomMargin(0.02);
  cr4->SetTopMargin(0.01);
  cr4->SetLeftMargin(0.12);
  cr4->SetRightMargin(0.02);
  cr4->Divide(2,3);

  //-------------------------------
  //   DRAW 3He
  //-------------------------------
  cr4->cd(1);
  gPad->SetLogy();
  gPad->SetTicky();
  gPad->SetTickx();
  hframe3->Draw();
  hB3_coalescence->Draw("l");
  gBlastB3vsR_PbPb502TeV->Draw("samel");
  paveptB3->Draw();
  paveLab3->Draw();

  //-------------------------------
  //   DRAW 3LH
  //-------------------------------
  cr4->cd(2);
  gPad->SetLogy();
  gPad->SetTicky();
  gPad->SetTickx();
  hframe3L->Draw();
  hB3L_coalescence_largeradius->Draw("samel");
  hB3L_coalescence->Draw("l");
  gBlastB3LambdavsR_PbPb502TeV->Draw("samel");
  paveptB3L->Draw();
  paveLab3L->Draw();

  //-------------------------------
  //   DRAW 4He
  //-------------------------------
  cr4->cd(3);
  gPad->SetLogy();
  gPad->SetTicky();
  gPad->SetTickx();
  hframe4->Draw();
  hB4_coalescence->Draw("samel");
  gBlastB4vsR_PbPb502TeV->Draw("samel");
  paveLab4->Draw();
  pavept->Draw();
  
  //-------------------------------
  //   DRAW 4LH
  //-------------------------------
  cr4->cd(4);
  gPad->SetLogy();
  gPad->SetTicky();
  gPad->SetTickx();
  hframe4L->Draw();
  hB4L_coalescence_largeradius->Draw("samel");
  hB4L_coalescence->Draw("samel");
  gBlastB4LambdavsR_PbPb502TeV->Draw("samel");
  paveptB4L->Draw();
  paveLab4L->Draw();
  
  //-------------------------------
  //   DRAW pseudodata
  //-------------------------------
  if (plotPseudoData) {
    cr4->cd(1);  gPseudo3He->Draw("samep2");
    cr4->cd(2);  gPseudo3LambdaH->Draw("samep2");
    cr4->cd(3);  gPseudo4He->Draw("samep2");
    cr4->cd(4);  gPseudo4LambdaH->Draw("samep2");
  }
  
  //-------------------------------
  //   DRAW Master legend
  //-------------------------------

  TLegend * masterLeg = new TLegend(0.3, 0.4, 0.9, 0.9, "");
  masterLeg->SetFillStyle(0);
  masterLeg->SetTextSize(0.05);
  masterLeg->SetBorderSize(0);
  masterLeg->AddEntry(gBlastB3vsR_PbPb502TeV, "BW + GSI-Heid. (#it{T}_{chem} = 156 MeV), Pb-Pb", "l");
  masterLeg->AddEntry(hB3_coalescence, "#it{B}_{#it{A}} coalescence", "l");
  masterLeg->AddEntry(hB3_coalescence, "#it{r} (^{3}He) = 2.48 fm", "");
  masterLeg->AddEntry(hB3L_coalescence, "#it{r} (^{3}_{#Lambda}H) = 6.8 fm", "");
  masterLeg->AddEntry(hB3L_coalescence, "#it{r} (^{4}_{#Lambda}H) = 2.4 fm", "");
  masterLeg->AddEntry(hB4_coalescence, "#it{r} (^{4}He) = 1.9 fm", "");
  masterLeg->AddEntry(hB3L_coalescence_largeradius, "#it{r} (^{3}_{#Lambda}H) = 14.1 fm", "l");
  masterLeg->AddEntry(hB4L_coalescence_largeradius, "#it{r} (^{4}_{#Lambda}H) = 4.9 fm", "");
  if (plotPseudoData)
    masterLeg->AddEntry(gPseudo3He, "ALICE pseudo-data, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV", "lp");
  cr4->cd(5);
  masterLeg->Draw();  
  cr4->SaveAs("BAmodels.png");

  //-------------------------------
  //-------------------------------
  //-------------------------------
  // DRAW SECOND OPTION -- 12 panels, pseudodata in the middle, separation in the bottom
  //-------------------------------
  //-------------------------------
  //-------------------------------
  TPaveText * paveALICE = new TPaveText(0.1, 0.7, 0.8, 0.9, "NDC");
  paveALICE->SetFillStyle(0);
  paveALICE->SetBorderSize(0);
  paveALICE->SetTextFont(42);
  paveALICE->SetTextSize(0.09);
  paveALICE->SetTextAlign(12);
  paveALICE->AddText(Form("#bf{ALICE Upgrade projection}"));
  paveALICE->AddText("Pb-Pb #sqrt{#it{s}_{NN}} = 5.5 TeV, L = 10 nb^{-1}");
  
  //TGraphErrors * gPseudoDeuteronUnc = (TGraphErrors *) GetPseudoBAvsMultiUnc("deuteron", 1);
  TGraphErrors * gPseudo3HeUnc = (TGraphErrors *) GetPseudoBAvsMultiUnc("he3", 1);
  TGraphErrors * gPseudo3LambdaHUnc = (TGraphErrors *) GetPseudoBAvsMultiUnc("hyperH3", 1);
  TGraphErrors * gPseudo4HeUnc = (TGraphErrors *) GetPseudoBAvsMultiUnc("he4", 1);
  TGraphErrors * gPseudo4LambdaHUnc = (TGraphErrors *) GetPseudoBAvsMultiUnc("hyperH4", 1);
  gPseudo3HeUnc->SetFillStyle(0);
  gPseudo3LambdaHUnc->SetFillStyle(0);
  gPseudo4HeUnc->SetFillStyle(0);
  gPseudo4LambdaHUnc->SetFillStyle(0);

  Double_t assumedRelSysUncMin = 0.10;
  Double_t assumedRelSysUncMax = 0.30;
  Double_t assumedNoSyst = 0.20;

  TGraphErrors * gSepar3HeMin = (TGraphErrors *) Thermal2CoalescenceSeparation(2, assumedRelSysUncMin);
  TGraphErrors * gSepar3LambdaHMin = (TGraphErrors *) Thermal2CoalescenceSeparation(3, assumedRelSysUncMin);
  TGraphErrors * gSepar4HeMin = (TGraphErrors *) Thermal2CoalescenceSeparation(4, assumedRelSysUncMin);
  //TGraphErrors * gSepar4LambdaHMin = (TGraphErrors *) Thermal2CoalescenceSeparation(5, assumedRelSysUncMin);

  gSepar3HeMin->SetLineStyle(1);
  gSepar3LambdaHMin->SetLineStyle(1);
  gSepar4HeMin->SetLineStyle(1);
  //gSepar4LambdaHMin->SetLineStyle(1);
  
  TGraphErrors * gSepar3HeMax = (TGraphErrors *) Thermal2CoalescenceSeparation(2, assumedRelSysUncMax);
  TGraphErrors * gSepar3LambdaHMax = (TGraphErrors *) Thermal2CoalescenceSeparation(3, assumedRelSysUncMax);
  TGraphErrors * gSepar4HeMax = (TGraphErrors *) Thermal2CoalescenceSeparation(4, assumedRelSysUncMax);
  //TGraphErrors * gSepar4LambdaHMax = (TGraphErrors *) Thermal2CoalescenceSeparation(5, assumedRelSysUncMax);

  gSepar3HeMax->SetLineStyle(3);
  gSepar3LambdaHMax->SetLineStyle(3);
  gSepar4HeMax->SetLineStyle(3);
  //gSepar4LambdaHMax->SetLineStyle(3);
  
  TH1D * hframeBA = new TH1D("hframeBA", "B_{A} vs radius; #it{R} (fm); #it{B}_{A}", 2000, 0.01, 6.5);
  hframeBA->GetXaxis()->SetTitleSize(0.06);
  hframeBA->GetYaxis()->SetTitleSize(0.06);
  hframeBA->GetYaxis()->SetTitleOffset(1.3);
  hframeBA->GetXaxis()->SetTitleOffset(0.8);
  hframeBA->GetXaxis()->SetLabelSize(0.05);
  hframeBA->GetYaxis()->SetLabelSize(0.05);
  hframeBA->GetYaxis()->SetRangeUser(2.E-12, 3.E-3);

  TH1D * hframeRStat = new TH1D("hframeRStat", "; #it{R} (fm); #sigma_{stat}/#it{B_{A}}", 2000, 0.01, 6.5);
  hframeRStat->GetXaxis()->SetTitleSize(0.09);
  hframeRStat->GetYaxis()->SetTitleSize(0.09);
  hframeRStat->GetYaxis()->SetTitleOffset(1.);
  hframeRStat->GetXaxis()->SetTitleOffset(1.);
  hframeRStat->GetXaxis()->SetLabelSize(0.07);
  hframeRStat->GetYaxis()->SetLabelSize(0.09);
  hframeRStat->GetYaxis()->SetNdivisions(507);
  hframeRStat->GetYaxis()->SetRangeUser(-1.2, 1.2);
  hframeRStat->GetYaxis()->SetTitle("relative stat. uncert.");

  TH1D * hframeSep = CreateFrameSeparation();
  hframeSep->GetXaxis()->SetRangeUser(0.01, 6.5);
  hframeSep->GetYaxis()->SetRangeUser(0.6, 200);
  hframeSep->GetYaxis()->SetTitleSize(0.09);
  hframeSep->GetYaxis()->SetTitleOffset(1.);
  hframeSep->GetYaxis()->SetLabelSize(0.09);
  hframeSep->GetXaxis()->SetTitleSize(0.09);
  hframeSep->GetXaxis()->SetTitleOffset(1.);
  
  hframeSep->GetYaxis()->SetTitle("separation power (N#sigma)");
  Float_t yl = 0.03; 
  Float_t yu = 0.23; 
  Float_t xl = 0.05; 
  Float_t xu = 0.5; 
  TLegend * masterLeg3 = new TLegend(xl + 0.2, yl, xu, yu, Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb3));
  masterLeg3->SetFillStyle(0);
  masterLeg3->SetTextSize(0.045);
  masterLeg3->SetBorderSize(0);
  masterLeg3->AddEntry(gBlastB3vsR_PbPb502TeV, "BW + GSI-Heid. (#it{T}_{chem} = 156 MeV)", "l");
  masterLeg3->AddEntry(hB3_coalescence, "Coal., #it{r}(^{3}He) = 2.48 fm", "l");

  TLegend * masterLeg3L = new TLegend(xl, yl, xu, yu, Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb3Lambda));
  masterLeg3L->SetFillStyle(0);
  masterLeg3L->SetTextSize(0.05);
  masterLeg3L->SetBorderSize(0);
  masterLeg3L->AddEntry(hB3L_coalescence, "Coal., #it{r} (^{3}_{#Lambda}H) = 6.8 fm", "l");
  masterLeg3L->AddEntry(hB3L_coalescence_largeradius, "Coal., #it{r} (^{3}_{#Lambda}H) = 14.1 fm", "l");
  
  TLegend * masterLeg4 = new TLegend(xl, yl+0.06, xu, yu, Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb4));
  masterLeg4->SetFillStyle(0);
  masterLeg4->SetTextSize(0.05);
  masterLeg4->SetBorderSize(0);
  masterLeg4->AddEntry(hB4_coalescence, "Coal., #it{r} (^{4}He) = 1.9 fm", "l");

  TLegend * masterLeg4L = new TLegend(xl, yl, xu, yu, Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb4Lambda));
  masterLeg4L->SetFillStyle(0);
  masterLeg4L->SetTextSize(0.05);
  masterLeg4L->SetBorderSize(0);
  masterLeg4L->AddEntry(hB3L_coalescence, "Coal., #it{r} (^{4}_{#Lambda}H) = 2.4 fm", "l");
  masterLeg4L->AddEntry(hB4L_coalescence_largeradius, "Coal., #it{r} (^{4}_{#Lambda}H) = 4.9 fm", "l");
  
  TLine * lineat1 = new TLine(0.01, 0., 6.5, 0.);
  lineat1->SetLineStyle(3);
  lineat1->SetLineWidth(1);

  TLine * lineat5 = new TLine(0.01, 5., 6.5, 5.);
  lineat5->SetLineStyle(3);
  lineat5->SetLineWidth(1);

  TCanvas * cr5 = new TCanvas("cr5", "compare thermal with coalescence", 1400, 900);
  cr5->SetMargin(0.15, 0.02, 0.02, 0.02);
  //cr5->Divide(4,1);

  Float_t lowx[5] = {0.001, 0.28, 0.52, 0.76, 0.99};
  Float_t lowy[4] = {0.001, 0.27, 0.52, 0.99};
  
  TPad * pad[4][3];
  for (int c = 0; c<4; c++) {
    for (int r = 0; r<3; r++) {
      pad[c][r] = new TPad(Form("pad%i%i",c,r),"pad", lowx[c], lowy[r], lowx[c+1], lowy[r+1]);
    }
    if (c==0) ConfigThreePanelsPad(pad[c][2], pad[c][1], pad[c][0], 0.2, 0.001);
    else  ConfigThreePanelsPad(pad[c][2], pad[c][1], pad[c][0], 0.001, 0.001);
  }
  
  // TPad * pad2 = new TPad("pad2","This is pad2", lowx[0], 0.001, upx[0], 0.40);
  // TPad * pad3 = new TPad("pad3","This is pad3", lowx[1], 0.400, upx[1], 0.99);
  // TPad * pad4 = new TPad("pad4","This is pad4", lowx[1], 0.001, upx[1], 0.40);
  // TPad * pad5 = new TPad("pad5","This is pad5", lowx[2], 0.400, upx[2], 0.99);
  // TPad * pad6 = new TPad("pad6","This is pad6", lowx[2], 0.001, upx[2], 0.40);
  // TPad * pad7 = new TPad("pad7","This is pad7", lowx[3], 0.400, upx[3], 0.99);
  // TPad * pad8 = new TPad("pad8","This is pad8", lowx[3], 0.001, upx[3], 0.40);
  
  // ConfigTwoPanelsPad(pad1, pad2, 0.2, 0.001);
  // ConfigTwoPanelsPad(pad3, pad4, 0.001, 0.001);
  // ConfigTwoPanelsPad(pad5, pad6, 0.001, 0.001);
  // ConfigTwoPanelsPad(pad7, pad8, 0.001, 0.001);

  //-------------------------------
  //   DRAW 3He
  //-------------------------------
  cr5->cd();
  pad[0][2]->Draw();
  pad[0][2]->cd();
  hframeBA->Draw();
  hB3_coalescence->Draw("l");
  gBlastB3vsR_PbPb502TeV->Draw("samel");
  paveLab3->Draw();
  masterLeg3->Draw();
  
  cr5->cd();
  pad[0][1]->Draw();
  pad[0][1]->cd();
  hframeRStat->Draw();
  lineat1->Draw();
  gPseudo3HeUnc->Draw("sameE2");
  // paveALICE->Draw();

  cr5->cd();
  pad[0][0]->Draw();
  pad[0][0]->cd();
  gPad->SetLogy();
  hframeSep->Draw();
  gSepar3HeMin->Draw("samel");
  gSepar3HeMax->Draw("samel");
  lineat5->Draw();

  //-------------------------------
  //   DRAW 3LH
  //-------------------------------
  cr5->cd();
  pad[1][2]->Draw();
  pad[1][2]->cd();
  hframeBA->Draw();
  hB3L_coalescence_largeradius->Draw("samel");
  hB3L_coalescence->Draw("l");
  gBlastB3LambdavsR_PbPb502TeV->Draw("samel");
  masterLeg3L->Draw();
  paveLab3L->Draw();
  
  
  cr5->cd();
  pad[1][1]->Draw();
  pad[1][1]->cd();
  hframeRStat->Draw();
  lineat1->Draw();
  gPseudo3LambdaHUnc->Draw("sameE2");

  cr5->cd();
  pad[1][0]->Draw();
  pad[1][0]->cd();
  gPad->SetLogy();
  hframeSep->Draw();
  gSepar3LambdaHMin->Draw("samel");
  gSepar3LambdaHMax->Draw("samel");
  lineat5->Draw();

  //-------------------------------
  //   DRAW 4He
  //-------------------------------
  
  cr5->cd();
  pad[2][2]->Draw();
  pad[2][2]->cd();
  hframeBA->Draw();
  hB4_coalescence->Draw("samel");
  gBlastB4vsR_PbPb502TeV->Draw("samel");
  paveLab4->Draw();
  masterLeg4->Draw();
  
  cr5->cd();
  pad[2][1]->Draw();
  pad[2][1]->cd();
  hframeRStat->Draw();
  lineat1->Draw();
  gPseudo4HeUnc->Draw("sameE2");

  cr5->cd();
  pad[2][0]->Draw();
  pad[2][0]->cd();
  gPad->SetLogy();
  hframeSep->Draw();
  gSepar4HeMin->Draw("samel");
  gSepar4HeMax->Draw("samel");
  lineat5->Draw();
  
  //-------------------------------
  //   DRAW 4LH
  //-------------------------------
  cr5->cd();
  pad[3][2]->Draw();
  pad[3][2]->cd();
  hframeBA->Draw();
  hB4L_coalescence_largeradius->Draw("samel");
  hB4L_coalescence->Draw("samel");
  gBlastB4LambdavsR_PbPb502TeV->Draw("samel");
  masterLeg4L->Draw();
  paveLab4L->Draw();

  cr5->cd();
  pad[3][1]->Draw();
  pad[3][1]->cd();
  hframeRStat->Draw();
  lineat1->Draw();
  gPseudo4LambdaHUnc->Draw("sameE2");

  cr5->cd();
  pad[3][0]->Draw();
  pad[3][0]->cd();
  paveALICE->Draw();
  //gSepar4LambdaHMin->Draw("l");
  //gSepar4LambdaHMax->Draw("samel");
 
  cr5->SaveAs("BAmodels_pseudoUnc.png");

  return;
  
}



TGraphErrors * GetPseudoBAvsMulti(TString particle,  Int_t paramSet)
{
  TFile * fin = TFile::Open("~/alice/nucleiB2/projectionsYR/ba_300818.root");
  TH1D * hist = 0x0;
  if (particle.Contains("deuteron")) hist = (TH1D *) fin->Get("badeuteron"); //ptoa = 0.75 GeV/c
  if (particle.Contains("triton")) hist = (TH1D *) fin->Get("batriton");  //ptoa = 0.77 GeV/c
  if (particle.Contains("he3")) hist = (TH1D *) fin->Get("bahe3"); //ptoa = 0.77 GeV/c
  if (particle.Contains("he4")) hist = (TH1D *) fin->Get("bahe4"); //ptoa = 0.75 GeV/c
  if (particle.Contains("hyperH3")) hist = (TH1D *) fin->Get("bahyperH3"); //ptoa = 1.17 GeV/c
  if (particle.Contains("hyperH4")) hist = (TH1D *) fin->Get("bahyperH4"); //ptoa = 0.62 GeV/c

  Double_t multi[10]    = {1942.5, 1585.5, 1180.0, 786.0, 512.0, 318.0, 183.0, 96.3, 44.9, 17.52};
  Double_t multierr[10] = {53.5, 46.0, 31.0, 20.0, 15.0, 12.0, 8.0, 5.8, 3.4, 1.84};
  
  Double_t ba[10], baerr[10];
  for (int i=1;i<11;i++){
    ba[i-1] = hist->GetBinContent(i);
    baerr[i-1] = hist->GetBinError(i);
  }
  
  TGraphErrors * graph = new TGraphErrors(10, multi, ba, multierr, baerr);
  convertMultiToRadius(graph, paramSet);
  graph->SetMarkerStyle(24);
  graph->SetMarkerColor(kBlack);
  graph->SetLineColor(kBlack);
  return graph;
}

TGraphErrors * GetPseudoBAvsMultiUnc(TString particle,  Int_t paramSet)
{
  TFile * fin = TFile::Open("~/alice/nucleiB2/projectionsYR/ba_300818.root");
  TH1D * hist = 0x0;
  if (particle.Contains("deuteron")) hist = (TH1D *) fin->Get("badeuteron"); //ptoa = 0.75 GeV/c
  if (particle.Contains("triton")) hist = (TH1D *) fin->Get("batriton");  //ptoa = 0.77 GeV/c
  if (particle.Contains("he3")) hist = (TH1D *) fin->Get("bahe3"); //ptoa = 0.77 GeV/c
  if (particle.Contains("he4")) hist = (TH1D *) fin->Get("bahe4"); //ptoa = 0.75 GeV/c
  if (particle.Contains("hyperH3")) hist = (TH1D *) fin->Get("bahyperH3"); //ptoa = 1.17 GeV/c
  if (particle.Contains("hyperH4")) hist = (TH1D *) fin->Get("bahyperH4"); //ptoa = 0.62 GeV/c

  Double_t multi[10]    = {1942.5, 1585.5, 1180.0, 786.0, 512.0, 318.0, 183.0, 96.3, 44.9, 17.52};
  Double_t multierr[10] = {53.5, 46.0, 31.0, 20.0, 15.0, 12.0, 8.0, 5.8, 3.4, 1.84};
  
  Double_t ba[10], baerr[10];
  for (int i=1;i<11;i++){
    Double_t relUnc = hist->GetBinError(i)/hist->GetBinContent(i);
    if (relUnc>=1.) ba[i-1] = -10.;
    else {
      ba[i-1] = 0.;
      if (hist->GetBinContent(i)>0) baerr[i-1] = relUnc;
      else baerr[i-1] = 0.0;
    }
  }
  
  TGraphErrors * graph = new TGraphErrors(10, multi, ba, multierr, baerr);
  convertMultiToRadius(graph, paramSet);
  graph->SetMarkerStyle(24);
  graph->SetMarkerColor(kBlack);
  graph->SetLineColor(kBlack);
  return graph;
}

void ConfigTwoPanelsPad(TPad* pad1, TPad* pad2, Float_t leftmargin = 0.15, Float_t rightmargin = 0.001)
{
  //top pad
  pad1->SetFillColor(0);
  pad1->SetBorderMode(0);
  pad1->SetBorderSize(0);
  pad1->SetMargin(leftmargin, rightmargin, 0.001, 0.001);
  pad1->SetLogy();
  pad1->SetTicky();
  pad1->SetTickx();
  //bottom pad
  pad2->SetFillColor(0);
  pad2->SetBorderMode(0);
  pad2->SetBorderSize(0);
  pad2->SetMargin(leftmargin, rightmargin, 0.15, 0.001);
  pad2->SetTicky();
  pad2->SetTickx();
  return;
}

void ConfigThreePanelsPad(TPad* pad1, TPad* pad2, TPad* pad3, Float_t leftmargin = 0.15, Float_t rightmargin = 0.001)
{
  //top pad
  pad1->SetFillColor(0);
  pad1->SetBorderMode(0);
  pad1->SetBorderSize(0);
  pad1->SetMargin(leftmargin, rightmargin, 0.001, 0.05);
  pad1->SetLogy();
  pad1->SetTicky();
  pad1->SetTickx();

  //middle pad
  pad2->SetFillColor(0);
  pad2->SetBorderMode(0);
  pad2->SetBorderSize(0);
  pad2->SetMargin(leftmargin, rightmargin, 0.001, 0.001);
  pad2->SetTicky();
  pad2->SetTickx();

  //bottom pad
  pad3->SetFillColor(0);
  pad3->SetBorderMode(0);
  pad3->SetBorderSize(0);
  pad3->SetMargin(leftmargin, rightmargin, 0.2, 0.001);
  pad3->SetTicky();
  pad3->SetTickx();
  return;
}
