/*
   fbellini@cern.ch
   17.09.2018 - Yellow Report plot
*/

#include "./CoalescenceBA.C"

void ConfigTwoPanelsPad(TPad* pad1, TPad* pad2, Float_t leftmargin, Float_t rightmargin);
void ConfigThreePanelsPad(TPad* pad1, TPad* pad2, TPad* pad3, Float_t leftmargin = 0.15, Float_t rightmargin = 0.001);
void LogPad();
//pseudodata
TGraphErrors * GetPseudoBAvsMulti(TString particle = "deuteron",  Int_t paramSet = 1);
TGraphErrors * GetPseudoBAvsMultiUnc(TString particle = "deuteron",  Int_t paramSet = 1, Color_t color = 0);
TGraphErrors * GetMeasBAvsMultiUnc(TString particle = "he3",  Int_t paramSet = 1, Color_t color = kGreen+2);

//main
void MakePubNoteFigures(Double_t pToAb3 = 0.77, Double_t pToAb3Lambda = 1.17, Double_t pToAb4 = 0.75, Double_t pToAb4Lambda = 0.62, Bool_t plotPseudoData = 1)
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


  Color_t projCol = kRed+1;
  Double_t Rmin = 0.01; Double_t Rmax = 6.5;
  Float_t yl = 0.2; 
  Float_t yu = 0.40; 
  Float_t xl = 0.06; 
  Float_t xu = 0.5; 
  Float_t lowx[5] = {0.001, 0.3, 0.54, 0.76, 0.99};
  Float_t lowy[4] = {0.001, 0.29, 0.52, 0.99};
  //------------------------------
  // Pseudodata
  //------------------------------
  TGraphErrors * gPseudo3He = (TGraphErrors *) GetPseudoBAvsMulti("he3", 1);
  TGraphErrors * gPseudo3LambdaH = (TGraphErrors *) GetPseudoBAvsMulti("hyperH3", 1);
  TGraphErrors * gPseudo4He = (TGraphErrors *) GetPseudoBAvsMulti("he4", 1);
  TGraphErrors * gPseudo4LambdaH = (TGraphErrors *) GetPseudoBAvsMulti("hyperH4", 1);

   //TGraphErrors * gPseudoDeuteronUnc = (TGraphErrors *) GetPseudoBAvsMultiUnc("deuteron", 1);
  TGraphErrors * gPseudo3HeUnc = (TGraphErrors *) GetPseudoBAvsMultiUnc("he3", 1, projCol);
  TGraphErrors * gPseudo3LambdaHUnc = (TGraphErrors *) GetPseudoBAvsMultiUnc("hyperH3", 1, projCol);
  TGraphErrors * gPseudo4HeUnc = (TGraphErrors *) GetPseudoBAvsMultiUnc("he4", 1, projCol);
  TGraphErrors * gPseudo4LambdaHUnc = (TGraphErrors *) GetPseudoBAvsMultiUnc("hyperH4", 1, projCol);
  gPseudo3HeUnc->SetFillStyle(1001);
  gPseudo3LambdaHUnc->SetFillStyle(1001);
  gPseudo4HeUnc->SetFillStyle(1001);
  gPseudo4LambdaHUnc->SetFillStyle(1001);

  //------------------------------
  // Rel unc measured Run 2
  //------------------------------
  TGraphErrors * gMeasHe3Unc = (TGraphErrors*) GetMeasBAvsMultiUnc("he3",  1, kBlack);
  TGraphErrors * gMeas3LambdaHUnc = (TGraphErrors*) GetMeasBAvsMultiUnc("hyperH3",  1, kBlack);
  gMeasHe3Unc->SetFillStyle(0);
  gMeas3LambdaHUnc->SetFillStyle(0);

  Double_t assumedRelSysUncMin = 0.10;
  Double_t assumedRelSysUncMax = 0.20;
  Double_t assumedNoSyst = 0.20;

  TGraphErrors * gSepar3HeMin = (TGraphErrors *) Thermal2CoalescenceSeparation(2, assumedRelSysUncMin);
  TGraphErrors * gSepar3LambdaHMin = (TGraphErrors *) Thermal2CoalescenceSeparation(3, assumedRelSysUncMin);
  TGraphErrors * gSepar4HeMin = (TGraphErrors *) Thermal2CoalescenceSeparation(4, assumedRelSysUncMin);
  //TGraphErrors * gSepar4LambdaHMin = (TGraphErrors *) Thermal2CoalescenceSeparation(5, assumedRelSysUncMin);

  gSepar3HeMin->SetLineStyle(1);
  gSepar3LambdaHMin->SetLineStyle(1);
  gSepar4HeMin->SetLineStyle(1);
  //gSepar4LambdaHMin->SetLineStyle(1);
  gSepar3HeMin->SetLineColor(projCol);
  gSepar3LambdaHMin->SetLineColor(projCol);
  gSepar4HeMin->SetLineColor(projCol);
  
  TGraphErrors * gSepar3HeMax = (TGraphErrors *) Thermal2CoalescenceSeparation(2, assumedRelSysUncMax);
  TGraphErrors * gSepar3LambdaHMax = (TGraphErrors *) Thermal2CoalescenceSeparation(3, assumedRelSysUncMax);
  TGraphErrors * gSepar4HeMax = (TGraphErrors *) Thermal2CoalescenceSeparation(4, assumedRelSysUncMax);
  //TGraphErrors * gSepar4LambdaHMax = (TGraphErrors *) Thermal2CoalescenceSeparation(5, assumedRelSysUncMax);

  gSepar3HeMax->SetLineStyle(3);
  gSepar3LambdaHMax->SetLineStyle(3);
  gSepar4HeMax->SetLineStyle(3);
  gSepar3HeMax->SetLineColor(projCol);
  gSepar3LambdaHMax->SetLineColor(projCol);
  gSepar4HeMax->SetLineColor(projCol);
  //gSepar4LambdaHMax->SetLineStyle(3);

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
  TH1D * hframeBA = new TH1D("hframeBA", "B_{A} vs radius; #it{R} (fm); #it{B}_{A} [(GeV^{2}/#it{c}^{3})^{A-1}]", 2000, 0.01, 6.5);
  hframeBA->GetXaxis()->SetTitleSize(0.09);
  hframeBA->GetXaxis()->SetTitleOffset(0.7);
  hframeBA->GetXaxis()->SetLabelSize(0.07);
  hframeBA->GetXaxis()->SetLabelOffset(-0.004);

  hframeBA->GetYaxis()->SetTitleSize(0.07);
  hframeBA->GetYaxis()->SetTitleOffset(1.5);
  hframeBA->GetYaxis()->SetLabelSize(0.06);
  hframeBA->GetYaxis()->SetRangeUser(2.E-12, 3.E-3);
 
  TH1D * hframeBA2 = new TH1D("hframeBA", "B_{A} vs radius; #it{R} (fm); #it{B}_{A} [(GeV^{2}/#it{c}^{3})^{A-1}]", 2000, 0.01, 6.5);
  hframeBA2->GetXaxis()->SetTitleSize(0.07);
  hframeBA2->GetXaxis()->SetTitleOffset(0.83);
  hframeBA2->GetXaxis()->SetLabelSize(0.05);
  hframeBA2->GetXaxis()->SetLabelOffset(0.008);

  hframeBA2->GetYaxis()->SetTitleSize(0.08);
  hframeBA2->GetYaxis()->SetTitleOffset(1.4);
  hframeBA2->GetYaxis()->SetLabelSize(0.06);
  hframeBA2->GetYaxis()->SetRangeUser(2.E-12, 3.E-3);

  TH1D * hframeRStat = new TH1D("hframeRStat", "; #it{R} (fm); #sigma_{stat}/#it{B_{A}}", 2000, 0.01, 6.5);
  hframeRStat->GetXaxis()->SetTitleSize(0.09);
  hframeRStat->GetXaxis()->SetTitleOffset(0.8);
  hframeRStat->GetXaxis()->SetLabelSize(0.06);
  hframeRStat->GetYaxis()->SetTitleSize(0.09);
  hframeRStat->GetYaxis()->SetTitleOffset(1.);
  hframeRStat->GetYaxis()->SetLabelSize(0.09);
  hframeRStat->GetYaxis()->SetLabelOffset(0.02);
  hframeRStat->GetYaxis()->SetNdivisions(507);
  hframeRStat->GetYaxis()->SetRangeUser(-1.2, 1.2);
  hframeRStat->GetYaxis()->SetTitle("rel. stat. unc.");

  TH1D * hframeSep = CreateFrameSeparation();
  hframeSep->GetXaxis()->SetRangeUser(0.01, 6.5);
  hframeSep->GetXaxis()->SetTitleSize(0.09);
  hframeSep->GetXaxis()->SetLabelSize(0.06);
  hframeSep->GetXaxis()->SetTitleOffset(0.75);
  hframeSep->GetYaxis()->SetRangeUser(0.6, 200);
  hframeSep->GetYaxis()->SetTitleSize(0.09);
  hframeSep->GetYaxis()->SetTitleOffset(1.);
  hframeSep->GetYaxis()->SetLabelSize(0.09);
  hframeSep->GetYaxis()->SetTitle("N#sigma separation");

  //define particle label
  TPaveText * paveLab2 = new TPaveText(0.8, 0.8, 0.9, 0.9, "NDC");
  paveLab2->SetFillStyle(0);
  paveLab2->SetTextFont(42);
  paveLab2->SetBorderSize(0);
  paveLab2->SetTextSize(0.1);
  paveLab2->SetTextAlign(12);
  paveLab2->AddText("#bf{d}");

  TPaveText * paveLab3 = new TPaveText(0.77, 0.83, 0.9, 0.93, "NDC");
  paveLab3->SetFillStyle(0);
  paveLab3->SetBorderSize(0);
  paveLab3->SetTextAlign(12);
  paveLab3->AddText("#bf{^{3}He}");
  paveLab3->SetTextSize(0.09);
  paveLab3->SetTextFont(42);
  // TLatex * paveLab3 = new TLatex(4.5, 2.e-4,"#bf{^{3}He}" );
  // paveLab3->SetTextSize(0.09);
  // paveLab3->SetTextFont(42);

  TPaveText * paveLab3L = new TPaveText(0.75, 0.83, 0.9, 0.93, "NDC");
  // TLegend * paveLab3L = new TLegend(0.75, 0.8, 0.9, 0.9, "#bf{{}^{3}_{#Lambda}H}");
  paveLab3L->SetFillStyle(0);
  paveLab3L->SetTextFont(42);
  paveLab3L->SetBorderSize(0);
  paveLab3L->SetTextSize(0.1);
  paveLab3L->SetTextAlign(12);
  paveLab3L->AddText("#bf{{}^{3}_{#Lambda}H}");
  // TLatex * paveLab3L = new TLatex(4.5, 2.e-4,"#bf{^{3} _{#Lambda}H}");
  // paveLab3L->SetTextSize(0.1);
  // paveLab3L->SetTextFont(42);

  TPaveText * paveLab4 = new TPaveText(0.77, 0.83, 0.9, 0.93, "NDC");
  paveLab4->SetFillStyle(0);
  paveLab4->SetTextFont(42);
  paveLab4->SetBorderSize(0);
  paveLab4->SetTextSize(0.1);
  paveLab4->SetTextAlign(12);
  paveLab4->AddText("#bf{^{4}He}");
  // TLatex * paveLab4 = new TLatex(4.5, 2.e-4,"#bf{^{4}He}" );
  // paveLab4->SetTextSize(0.1);
  // paveLab4->SetTextFont(42);

  TPaveText * paveLab4L = new TPaveText(0.77, 0.83, 0.9, 0.93, "NDC");
  paveLab4L->SetFillStyle(0);
  paveLab4L->SetTextFont(42);
  paveLab4L->SetBorderSize(0);
  paveLab4L->SetTextSize(0.1);
  paveLab4L->SetTextAlign(12);
  paveLab4L->AddText("#bf{{}^{4}_{#Lambda}H}");
  // TLatex * paveLab4L = new TLatex(4.5, 2.e-4,"#bf{{}^{4} _{#Lambda}H}");
  // paveLab4L->SetTextSize(0.1);
  // paveLab4L->SetTextFont(42);
  
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
  //   DRAW ALICE label
  //-------------------------------
  TPaveText * paveALICE = new TPaveText(0.01, 0.6, 0.9, 0.9, "NDC");
  paveALICE->SetFillStyle(0);
  paveALICE->SetBorderSize(0);
  paveALICE->SetTextFont(42);
  paveALICE->SetTextSize(0.08);
  paveALICE->SetTextAlign(12);
  paveALICE->AddText(Form("#bf{ALICE Upgrade projection}"));
  paveALICE->AddText("Pb-Pb #sqrt{#it{s}_{NN}} = 5.5 TeV, #it{L}_{int} = 10 nb^{-1}");

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
    
  TLegend * masterLeg3 = new TLegend(xl + 0.2, yl, xu, yu, Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb3));
  masterLeg3->SetFillStyle(0);
  masterLeg3->SetTextSize(0.045);
  masterLeg3->SetBorderSize(0);
  masterLeg3->AddEntry(gBlastB3vsR_PbPb502TeV, "BW + GSI-Heid. (#it{T}_{chem} = 156 MeV)", "l");
  masterLeg3->AddEntry(hB3_coalescence, "Coal., #it{r}(^{3}He) = 2.48 fm", "l");
  // masterLeg3->AddEntry(gPseudo3He, "Pb-Pb #sqrt{#it{s}_{NN}} = 5.5 TeV, L_{int} = 10 nb^{-1}","lp");


  TLegend * masterLeg3L = new TLegend(xl, yl, xu, yu, Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb3Lambda));
  masterLeg3L->SetFillStyle(0);
  masterLeg3L->SetTextSize(0.055);
  masterLeg3L->SetBorderSize(0);
  masterLeg3L->AddEntry(hB3L_coalescence, "Coal., #it{r} (^{3}_{#Lambda}H) = 6.8 fm", "l");
  masterLeg3L->AddEntry(hB3L_coalescence_largeradius, "Coal., #it{r} (^{3}_{#Lambda}H) = 14.1 fm", "l");
  
  TLegend * masterLeg4 = new TLegend(0.3, 0.67, 0.95, 0.8, Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb4));
  masterLeg4->SetFillStyle(0);
  masterLeg4->SetTextSize(0.055);
  masterLeg4->SetBorderSize(0);
  masterLeg4->AddEntry(hB4_coalescence, "Coal., #it{r} (^{4}He) = 1.9 fm", "l");

  TLegend * masterLeg4L = new TLegend(0.3, 0.60, 0.95, 0.8, Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb4Lambda));
  masterLeg4L->SetFillStyle(0);
  masterLeg4L->SetTextSize(0.055);
  masterLeg4L->SetBorderSize(0);
  masterLeg4L->AddEntry(hB3L_coalescence, "Coal., #it{r} (^{4}_{#Lambda}H) = 2.4 fm", "l");
  masterLeg4L->AddEntry(hB4L_coalescence_largeradius, "Coal., #it{r} (^{4}_{#Lambda}H) = 4.9 fm", "l");

  TLegend * uncLegend = new TLegend(0.24, 0.05, 0.45, 0.35);
  uncLegend->SetFillStyle(0);
  uncLegend->SetTextSize(0.06);
  uncLegend->SetBorderSize(0);
  uncLegend->AddEntry( gPseudo3HeUnc, "Pb-Pb #sqrt{#it{s}_{NN}} = 5.5 TeV (Run 3+4)","f");
  uncLegend->AddEntry( gMeasHe3Unc, "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV (2015 data)","f");
    
  TLegend * sepLegend = new TLegend(0.24, 0.65, 0.5, 0.95, "#sigma^{2} = #sigma_{stat}^{2} + #sigma_{sys}^{2}");
  sepLegend->SetFillStyle(0);
  sepLegend->SetTextSize(0.06);
  sepLegend->SetBorderSize(0);
  sepLegend->AddEntry( gSepar3HeMin, Form("#sigma_{sys}/B_{A} = %2.0f%%",assumedRelSysUncMin*100),"l");
  sepLegend->AddEntry( gSepar3HeMax, Form("#sigma_{sys}/B_{A} = %2.0f%%",assumedRelSysUncMax*100),"l");
  
  TLine * lineat1 = new TLine(0.01, 0., 6.5, 0.);
  lineat1->SetLineStyle(3);
  lineat1->SetLineWidth(1);

  TLine * lineat5 = new TLine(0.01, 5., 6.5, 5.);
  lineat5->SetLineStyle(3);
  lineat5->SetLineWidth(1);


  //-------------------------------
  // DRAW THIRD OPTION -- 12 panels, pseudodata in the middle, separation in the bottom
  //-------------------------------
  TCanvas * cr6 = new TCanvas("cr6", "compare thermal with coalescence", 1400, 550);
  cr6->SetMargin(0.15, 0.02, 0.02, 0.02);
  TPad * pad_models[4];
  for (int c = 0; c<4; c++) {
    pad_models[c]= new TPad(Form("pad_models%i",c),"pad", lowx[c], 0.001, lowx[c+1], 0.99);
    pad_models[c]->SetMargin(0.0, 0.0, 0.15, 0.02);
  }
  pad_models[0]->SetMargin(0.23, 0.0, 0.15, 0.02);
  pad_models[3]->SetMargin(0.0, 0.02, 0.15, 0.02);
  
  cr6->cd();
  pad_models[0]->Draw();
  pad_models[0]->cd();
  LogPad();
  hframeBA2->Draw();
  hB3_coalescence->Draw("l");
  gBlastB3vsR_PbPb502TeV->Draw("samel");
  paveLab3->Draw();
  masterLeg3->Draw();

  cr6->cd();
  pad_models[1]->Draw();
  pad_models[1]->cd();
  LogPad();
  hframeBA->Draw();
  hB3L_coalescence_largeradius->Draw("samel");
  hB3L_coalescence->Draw("l");
  gBlastB3LambdavsR_PbPb502TeV->Draw("samel");
  masterLeg3L->Draw();
  paveLab3L->Draw();
  
  cr6->cd();
  pad_models[2]->Draw();
  pad_models[2]->cd();
  LogPad();
  hframeBA->Draw();
  hB4_coalescence->Draw("samel");
  gBlastB4vsR_PbPb502TeV->Draw("samel");
  paveLab4->Draw();
  masterLeg4->Draw();
  
  cr6->cd();
  pad_models[3]->Draw();
  pad_models[3]->cd();
  LogPad();
  hframeBA->Draw();
  hB4L_coalescence_largeradius->Draw("samel");
  hB4L_coalescence->Draw("samel");
  gBlastB4LambdavsR_PbPb502TeV->Draw("samel");
  masterLeg4L->Draw();
  paveLab4L->Draw();

  cr6->SaveAs("BA_models.png");
  cr6->SaveAs("BA_models.eps");
  cr6->SaveAs("BA_models.pdf");

//-------------------------------------
  TCanvas * cr7 = new TCanvas("cr7", "projections", 1400, 700);
  cr7->SetMargin(0.15, 0.02, 0.02, 0.02);
  
  Float_t lowy_3opt[3] = {0.001, 0.51, 0.99};
  Float_t lowx_3opt[5] = {0.001, 0.3, 0.54, 0.76, 0.99};
  TPad * pad_3opt[4][2];
  pad_3opt[3][1] = new TPad(Form("pad_%i%i",3,1),"pad", lowx_3opt[3], 0.38, lowx_3opt[4], 0.99);
  pad_3opt[3][0] = new TPad(Form("pad_%i%i",3,1),"pad", lowx_3opt[3], 0.1, lowx_3opt[4], 0.37);
  pad_3opt[3][1]->SetRightMargin(0.02);
  pad_3opt[3][0]->SetRightMargin(0.02);
  pad_3opt[3][1]->SetLeftMargin(0.0);
  pad_3opt[3][0]->SetLeftMargin(0.0);
  pad_3opt[3][1]->SetTopMargin(0.05);
  pad_3opt[3][0]->SetTopMargin(0.0);
  pad_3opt[3][1]->SetBottomMargin(0.23);
  pad_3opt[3][0]->SetBottomMargin(0.0);
  
    for (int r=0; r<2;r++){
      pad_3opt[3][r]->SetTickx(); 
      pad_3opt[3][r]->SetTicky(); 
    }
  for (int r=0; r<2;r++){
    for (int c=0; c<3;c++){
        pad_3opt[c][r] = new TPad(Form("pad_%i%i",c,r),"pad", lowx_3opt[c], lowy_3opt[r], lowx_3opt[c+1], lowy_3opt[r+1]);
        pad_3opt[c][r]->SetRightMargin(0.0);
        pad_3opt[c][r]->SetLeftMargin((c==0? 0.2 : 0.0));
        pad_3opt[c][r]->SetTickx(); 
        pad_3opt[c][r]->SetTicky(); 
        if (r==0) {
            pad_3opt[c][r]->SetBottomMargin(0.2);
            pad_3opt[c][r]->SetTopMargin(0.);
        } else {
            pad_3opt[c][r]->SetBottomMargin(0.);
            pad_3opt[c][r]->SetTopMargin(0.02);
        } 
        if (c==3 && r==1) pad_3opt[c][r]->SetBottomMargin(0.2);
    }
    pad_3opt[3][r]->SetRightMargin(0.02);
    pad_3opt[3][r]->SetLeftMargin(0.0);
    pad_3opt[3][r]->SetTopMargin((r==0? 0.0 : 0.0));
    pad_3opt[3][r]->SetBottomMargin(0.2);
    }


  //-------------------------------
  //   DRAW 3He
  //-------------------------------
  cr7->cd();
  pad_3opt[0][1]->Draw();
  pad_3opt[0][1]->cd();
  hframeRStat->Draw("0");
  lineat1->Draw();
  gMeasHe3Unc->Draw("sameE2");
  gPseudo3HeUnc->Draw("sameE2");
  paveLab3->Draw();
  uncLegend->Draw();

  cr7->cd();
  pad_3opt[0][0]->Draw();
  pad_3opt[0][0]->cd();
  gPad->SetLogy();
  hframeSep->Draw();
  gSepar3HeMin->Draw("samel");
  gSepar3HeMax->Draw("samel");
  lineat5->Draw();
  sepLegend->Draw();
  
  //-------------------------------
  //   DRAW 3LH
  //-------------------------------
 
  cr7->cd();
  pad_3opt[1][1]->Draw();
  pad_3opt[1][1]->cd();
  hframeRStat->Draw("0");
  lineat1->Draw();
  paveLab3L->Draw();
  gMeas3LambdaHUnc->Draw("sameE2");
  gPseudo3LambdaHUnc->Draw("sameE2");
  
  cr7->cd();
  pad_3opt[1][0]->Draw();
  pad_3opt[1][0]->cd();
  gPad->SetLogy();
  hframeSep->Draw();
  gSepar3LambdaHMin->Draw("samel");
  gSepar3LambdaHMax->Draw("samel");
  lineat5->Draw();

  //-------------------------------
  //   DRAW 4He
  //-------------------------------

  cr7->cd();
  pad_3opt[2][1]->Draw();
  pad_3opt[2][1]->cd();
  hframeRStat->Draw("0");
  lineat1->Draw();
  paveLab4->Draw();
  gPseudo4HeUnc->Draw("sameE2");

  cr7->cd();
  pad_3opt[2][0]->Draw();
  pad_3opt[2][0]->cd();
  gPad->SetLogy();
  hframeSep->Draw();
  gSepar4HeMin->Draw("samel");
  gSepar4HeMax->Draw("samel");
  lineat5->Draw();
  
  //-------------------------------
  //   DRAW 4LH
  //-------------------------------

  cr7->cd();
  pad_3opt[3][1]->Draw();
  pad_3opt[3][1]->cd();
  hframeRStat->Draw("0");
  lineat1->Draw();
  paveLab4L->Draw();
  gPseudo4LambdaHUnc->Draw("sameE2");

  cr7->cd();
  pad_3opt[3][0]->Draw();
  pad_3opt[3][0]->cd();
  paveALICE->Draw();
 
  cr7->SaveAs("BA_projected_pseudoUnc.png");
  cr7->SaveAs("BA_projected_pseudoUnc.eps");
  cr7->SaveAs("BA_projected_pseudoUnc.pdf");

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
  if (particle.Contains("hyperH4"))
    graph = new TGraphErrors(5, multi, ba, multierr, baerr);

  convertMultiToRadius(graph, paramSet);
  graph->SetMarkerStyle(24);
  graph->SetMarkerColor(kBlack);
  graph->SetLineColor(kBlack);
  return graph;
}

TGraphErrors * GetPseudoBAvsMultiUnc(TString particle,  Int_t paramSet, Color_t color)
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
  graph->SetMarkerColor(color);
  graph->SetLineColor(color);
  graph->SetFillColor(color);
  graph->SetFillStyle(0);
  return graph;
}

TGraphErrors * GetMeasBAvsMultiUnc(TString particle,  Int_t paramSet, Color_t color)
{

  TGraphErrors * graph  = 0x0;

  if (particle.Contains("he3")){
    Double_t multi[]    = {1942.5, (1585.5 + 1180.0 + 786.0) / 3., (512.0 + 318.0 + 183.0 + 96.3 + 44.9 + 17.52 ) / 5.};
    Double_t multierr[] = {53.5, (46.0 + 31.0 + 20.0) / 3., (15.0 + 12.0 + 8.0 + 5.8 + 3.4 + 1.84)/ 5.};
    Double_t ba[3] = { 0., 0., 0.};
    Double_t baerr[3] = { 0.07796, 0.05209, 0.07791};
    graph = new TGraphErrors(3, multi, ba, multierr, baerr);
  }
  if (particle.Contains("hyperH3")){
    Double_t multi[]    = { (1585.5 + 1180.0 + 786.0) / 3.};
    Double_t multierr[] = { (46.0 + 31.0 + 20.0) / 3.};
    Double_t ba[1] = { 0.};
    Double_t baerr[1] = { 0.3106};
    graph = new TGraphErrors(1, multi, ba, multierr, baerr);
  }
  
  convertMultiToRadius(graph, paramSet);
  graph->SetMarkerStyle(24);
  graph->SetMarkerColor(color);
  graph->SetLineColor(color);
  graph->SetFillColor(color);
  graph->SetFillStyle(0);
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

void ConfigThreePanelsPad(TPad* pad1, TPad* pad2, TPad* pad3, Float_t leftmargin = 0.15, Float_t rightmargin = 0.01)
{
  //top pad
  pad1->SetFillColor(0);
  pad1->SetBorderMode(0);
  pad1->SetBorderSize(0);
  pad1->SetMargin(leftmargin, rightmargin, 0.0, 0.05);
  pad1->SetLogy();
  pad1->SetTicky();
  pad1->SetTickx();

  //middle pad
  pad2->SetFillColor(0);
  pad2->SetBorderMode(0);
  pad2->SetBorderSize(0);
  pad2->SetMargin(leftmargin, rightmargin, 0.0, 0.0);
  pad2->SetTicky();
  pad2->SetTickx();

  //bottom pad
  pad3->SetFillColor(0);
  pad3->SetBorderMode(0);
  pad3->SetBorderSize(0);
  pad3->SetMargin(leftmargin, rightmargin, 0.2, 0.00);
  pad3->SetTicky();
  pad3->SetTickx();
  return;
}



void LogPad(){

  gPad->SetLogy();
  gPad->SetTicky();
  gPad->SetTickx();
  return;

}