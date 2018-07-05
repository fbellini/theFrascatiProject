#include "B2vsVolume.C"

TGraphErrors * GetPPbHbtRadius025KT();
TGraphErrors * GetPPbHbtRadius0887KT();
TGraphErrors * GetPbPbHbtRadius0887KT();


void HbtRadii() {
  //
  // small macro to plot HBT radii versus multiplicity to compare
  // them to the parameterisation
  //
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetLineWidth(1);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleOffset(0.8,"y");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TGaxis::SetMaxDigits(4);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin(0.05);

  TCanvas * canvHbtPlot = new TCanvas("canvHbtPlot","canvHbtPlot");
  TH2D * histFrame = new TH2D("histFrame","radius vs mult; #LTd#it{N}_{ch}/d#it{#eta}#GT^{1/3}; #it{R} (fm)",200,0.0,13.0,200,0.0,6.0);
  histFrame->GetXaxis()->SetTitleSize(0.06);
  histFrame->GetXaxis()->SetLabelSize(0.045);
  histFrame->GetXaxis()->SetTitleOffset(1.1);
  histFrame->GetYaxis()->SetTitleOffset(1.2);
  histFrame->GetYaxis()->SetTitleSize(0.06);
  histFrame->GetYaxis()->SetLabelOffset(0.01);
  histFrame->GetYaxis()->SetLabelSize(0.045);

  histFrame->GetYaxis()->SetTitleOffset(1.);
  histFrame->Draw();
  //
  const Long_t nPoints = 500;
  Double_t xP[nPoints];
  Double_t yP[nPoints];
  for (int i = 0; i<nPoints; i++){
    xP[i] = 15.0 * i / nPoints; 
    Double_t multi[2] = {xP[i]*xP[i]*xP[i], 0.0};
    Double_t radius[2] = {0.0, 0.0};
    getRadiusFromParameterisation(multi, radius);
    yP[i] = radius[0];
    //
  }
  //
  // TGraph * grRadiusVsMult = new TGraph(nPoints, xP, yP);
  // grRadiusVsMult->SetMarkerColor(kBlack);
  // grRadiusVsMult->SetLineColor(kBlack);
  // grRadiusVsMult->SetLineWidth(3);
  // grRadiusVsMult->SetLineStyle(9);
  //
  //grRadiusVsMult->Draw("L");
  //

  // add PPB points
  //
  TGraphErrors * grHbtRadiusPPB = GetPPbHbtRadius0887KT();
  grHbtRadiusPPB->SetMarkerStyle(21);
  grHbtRadiusPPB->SetMarkerSize(1.);
  grHbtRadiusPPB->SetLineWidth(1);
  grHbtRadiusPPB->SetLineStyle(1);
  grHbtRadiusPPB->SetLineColor(kBlue);
  //
  // add PbPb points
  //
  TGraphErrors * grHbtRadiusPBPB = GetPbPbHbtRadius0887KT();
  grHbtRadiusPBPB->SetMarkerStyle(20);
  grHbtRadiusPBPB->SetMarkerSize(1.0);
  grHbtRadiusPBPB->SetLineWidth(1);
  grHbtRadiusPBPB->SetLineStyle(1);
  grHbtRadiusPBPB->SetLineColor(kRed);
  //
  // add a pp point
  //http://aliceinfo.cern.ch/ArtSubmission/node/1885, value is from Kfir
  TGraphErrors * grHbtRadiusPP = new TGraphErrors(1);
  grHbtRadiusPP->SetPoint(0, TMath::Power(5.98,1./3.), 0.8);
  grHbtRadiusPP->SetPointError(0, TMath::Power(5.98,-2./3.)*0.09/3., 0.3);
  grHbtRadiusPP->SetLineColor(kGreen+1);
  grHbtRadiusPP->SetMarkerColor(kGreen+1);
  grHbtRadiusPP->SetMarkerStyle(22);

    //canvHbtPlot->BuildLegend()
  TF1 * poly1 = new TF1("poly1", "[0]*x + [1]", 0., 13.);
  poly1->SetMarkerColor(kBlack);
  poly1->SetLineColor(kBlack);
  poly1->SetLineWidth(3);
  poly1->SetLineStyle(9);
  poly1->SetParameter(1, 0.0);

  TMultiGraph * multig = new TMultiGraph("multig", "multigraph");
  multig->Add(grHbtRadiusPP);
  multig->Add(grHbtRadiusPPB);
  multig->Add(grHbtRadiusPBPB);
  multig->Fit(poly1, "RS");
  multig->Draw("PZ");

  grHbtRadiusPP->Draw("PZ");
  grHbtRadiusPPB->Draw("PZ");
  grHbtRadiusPBPB->Draw("PZ");

  TLegend * leg1 = new TLegend(0.2, 0.82, 0.5, 0.92, "Our parameterisation:");
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.04);
  leg1->AddEntry(poly1, "#it{R} = #it{a} #times #LTd#it{#it{N}_{ch}}/d#it{#eta}#GT^{1/3} + #it{b}", "l");
  leg1->Draw();
  TLegend * leg = new TLegend(0.2, 0.6, 0.5, 0.8, "ALICE, #it{k}_{T} = 0.887 GeV/#it{c}");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(grHbtRadiusPP,"pp, #sqrt{s} = 7 TeV", "p");
  leg->AddEntry(grHbtRadiusPPB,"p-Pb, #sqrt{s_{NN}} = 5.02 TeV", "p");
  leg->AddEntry(grHbtRadiusPBPB, "Pb-Pb, #sqrt{s_{NN}} = 2.76 TeV", "p");
  leg->Draw();
  
  TString text="Uncertainties: #sqrt{stat.^{2} + sys.^{2}}";
  TPaveText * pave = new TPaveText(0.2,0.53,0.5,0.58,"brNDC");
  pave->SetBorderSize(0);
  pave->SetFillColor(kWhite);
  pave->SetTextColor(kBlack);
  pave->SetTextFont(42);
  pave->SetTextSize(0.04);
  pave->InsertText(text.Data());
  pave->Draw();

  canvHbtPlot->Print("HbtRadiusParam.eps");
}

TGraphErrors * GetPbPbHbtRadius0887KT() {
  //
  // geometric mean of Rout, Rside, Rlong
  // http://hepdata.cedar.ac.uk/view/ins1384807
  //
  // Rout and dN/deta
  double p8772_d8x1y1_xval[] = { 149., 261., 426., 649., 966., 1294., 1601. };
  double p8772_d8x1y1_xerrminus[] = { 6., 9., 15., 23., 37., 49., 60.};
  double p8772_d8x1y1_xerrplus[] = {6., 9., 15., 23., 37., 49., 60. };
  double p8772_d8x1y1_yval[] = { 1.83, 2.12, 2.39, 2.84, 3.23, 3.5, 3.87};
  double p8772_d8x1y1_yerrminus[] = { 0.21, 0.2, 0.21, 0.23, 0.27, 0.37, 0.35 };
  double p8772_d8x1y1_yerrplus[] = {  0.21, 0.2, 0.21, 0.23, 0.27, 0.37, 0.35};
  double p8772_d8x1y1_ystatminus[] = { 0.03, 0.03, 0.02, 0.03, 0.03, 0.04, 0.04};
  double p8772_d8x1y1_ystatplus[] = { 0.03, 0.03, 0.02, 0.03, 0.03, 0.04, 0.04};
  //
  //
  // Rside
  double p8772_d8x1y2_yval[] = { 2.26, 2.53, 2.76, 3.26, 3.62, 3.76, 4.05};
  double p8772_d8x1y2_yerrminus[] = { 0.22, 0.3, 0.26, 0.28, 0.33, 0.34, 0.46};
  double p8772_d8x1y2_yerrplus[] = {  0.22, 0.3, 0.26, 0.28, 0.33, 0.34, 0.46};
  double p8772_d8x1y2_ystatminus[] = { 0.04, 0.03, 0.02, 0.03, 0.03, 0.04, 0.04};
  double p8772_d8x1y2_ystatplus[] = {  0.04, 0.03, 0.02, 0.03, 0.03, 0.04, 0.04 };
  //
  //
  // Rlong
  double p8772_d8x1y3_yval[] = { 1.94, 2.27, 2.52, 3.06, 3.39, 3.55, 3.83};
  double p8772_d8x1y3_yerrminus[] = { 0.28, 0.28, 0.32, 0.34, 0.34, 0.41, 0.38};
  double p8772_d8x1y3_yerrplus[] = { 0.28, 0.28, 0.32, 0.34, 0.34, 0.41, 0.38};
  double p8772_d8x1y3_ystatminus[] = { 0.04, 0.03, 0.02, 0.04, 0.03, 0.04, 0.04 };
  double p8772_d8x1y3_ystatplus[] = { 0.04, 0.03, 0.02, 0.04, 0.03, 0.04, 0.04 };
  //
  // new TGraphErrors
  //
  Double_t xP[7];
  Double_t yP[7];
  Double_t xPerr[7];
  Double_t yPerr[7];
  //
  for (Int_t iP = 0; iP <7; iP++) {
    xP[iP] = TMath::Power(p8772_d8x1y1_xval[iP], 1./3.);
    xPerr[iP] =  TMath::Power(p8772_d8x1y1_xval[iP], -2./3.) * p8772_d8x1y1_xerrplus[iP]/3.;
    yP[iP] = TMath::Power(p8772_d8x1y2_yval[iP]*p8772_d8x1y2_yval[iP]*p8772_d8x1y1_yval[iP], 1./3.);
    Double_t errSyst = TMath::Power(p8772_d8x1y2_yerrplus[iP]*p8772_d8x1y2_yerrplus[iP]*p8772_d8x1y1_yerrplus[iP], 1./3.);
    Double_t errStat = TMath::Power(p8772_d8x1y2_ystatplus[iP]*p8772_d8x1y2_ystatplus[iP]*p8772_d8x1y1_ystatplus[iP], 1./3.);
    yPerr[iP] = TMath::Sqrt(errSyst*errSyst + errStat*errStat);;
  }
  TGraphErrors * grHbtRadiusPBPB = new TGraphErrors(7, xP, yP, xPerr, yPerr);
  grHbtRadiusPBPB->SetName("grHbtRadiusPBPB");
  grHbtRadiusPBPB->SetTitle("#pi Pb-Pb 2.76 TeV  #it{k}_{T}=0.887 GeV/it{c} #sqrt[3]{R_{out}R_{side}R_{long}}");
  //
  grHbtRadiusPBPB->SetMarkerSize(1.0);
  grHbtRadiusPBPB->SetMarkerStyle(21);
  grHbtRadiusPBPB->SetMarkerColor(kRed);
  grHbtRadiusPBPB->SetLineColor(kRed);
  //
  return grHbtRadiusPBPB;


}


TGraphErrors * GetPPbHbtRadius025KT() {
  //
  // geometric mean of Rout, Rside, Rlong
  // taken from: http://hepdata.cedar.ac.uk/view/ins1342499
  //
  // Pions at kT = 0.25 GeV/c
  //
  // Rout and dN/deta
  double p8772_d8x1y1_xval[] = { 1.88, 2.53, 2.87, 3.27 };
  double p8772_d8x1y1_xerrminus[] = { 0.0, 0.0, 0.0, 0.0 };
  double p8772_d8x1y1_xerrplus[] = { 0.0, 0.0, 0.0, 0.0 };
  double p8772_d8x1y1_yval[] = { 1.17, 1.28, 1.32, 1.53 };
  double p8772_d8x1y1_yerrminus[] = { 0.16278820596099708, 0.15297058540778355, 0.13341664064126335, 0.16278820596099708 };
  double p8772_d8x1y1_yerrplus[] = { 0.16278820596099708, 0.15297058540778355, 0.13341664064126335, 0.16278820596099708 };
  double p8772_d8x1y1_ystatminus[] = { 0.03, 0.03, 0.03, 0.03 };
  double p8772_d8x1y1_ystatplus[] = { 0.03, 0.03, 0.03, 0.03 };
  //
  //
  // Rside
  double p8772_d8x1y2_yval[] = { 1.25, 1.44, 1.56, 1.82 };
  double p8772_d8x1y2_yerrminus[] = { 0.14317821063276354, 0.13341664064126335, 0.16278820596099708, 0.1726267650163207 };
  double p8772_d8x1y2_yerrplus[] = { 0.14317821063276354, 0.13341664064126335, 0.16278820596099708, 0.1726267650163207 };
  double p8772_d8x1y2_ystatminus[] = { 0.03, 0.03, 0.03, 0.03 };
  double p8772_d8x1y2_ystatplus[] = { 0.03, 0.03, 0.03, 0.03 };
  //
  //
  // Rlong
  double p8772_d8x1y3_yval[] = { 1.37, 1.58, 1.67, 2.0 };
  double p8772_d8x1y3_yerrminus[] = { 0.20223748416156687, 0.18248287590894657, 0.1726267650163207, 0.21213203435596423 };
  double p8772_d8x1y3_yerrplus[] = { 0.20223748416156687, 0.18248287590894657, 0.1726267650163207, 0.21213203435596423 };
  double p8772_d8x1y3_ystatminus[] = { 0.03, 0.03, 0.03, 0.03 };
  double p8772_d8x1y3_ystatplus[] = { 0.03, 0.03, 0.03, 0.03 };
  //
  // new TGraphErrors
  //
  Double_t xP[4];
  Double_t yP[4];
  Double_t xPerr[4];
  Double_t yPerr[4];
  //
  for (Int_t iP =0; iP <4; iP++) {
    xP[iP] = p8772_d8x1y1_xval[iP];
    xPerr[iP] = p8772_d8x1y1_xerrplus[iP];
    yP[iP] = TMath::Power(p8772_d8x1y2_yval[iP]*p8772_d8x1y2_yval[iP]*p8772_d8x1y1_yval[iP], 1./3.);
    yPerr[iP] = TMath::Power(p8772_d8x1y2_yerrplus[iP]*p8772_d8x1y2_yerrplus[iP]*p8772_d8x1y1_yerrplus[iP], 1./3.);
  }
  TGraphErrors * grHbtRadiusPPB = new TGraphErrors(4, xP, yP, xPerr, yPerr);
  grHbtRadiusPPB->SetName("grHbtRadiusPPB");

  return grHbtRadiusPPB;
}

TGraphErrors * GetPPbHbtRadius0887KT() {
  //
  // geometric mean of Rout, Rside, Rlong (gaus-gaus-gaus)
  // taken from: http://hepdata.cedar.ac.uk/view/ins1342499 
  //
  // Pions at kT = 0.887 GeV/c
  //
  // Rout and dN/deta
  double p8772_d8x1y1_xval[] = { 2.53, 2.87, 3.27 };
  double p8772_d8x1y1_xerrminus[] = { TMath::Power(2.53, -2./3.)*0.4/3, TMath::Power(2.87,-2./3.)*0.5/3., TMath::Power(3.27,-2./3.)*0.8/3. };
  double p8772_d8x1y1_xerrplus[] = { 0.0, 0.0, 0.0 };
  double p8772_d8x1y1_yval[] = { 0.68, 0.77, 0.94 };
  double p8772_d8x1y1_yerrminus[] = { 0.17, 0.18, 0.23 };
  double p8772_d8x1y1_yerrplus[] = { 0.17, 0.18, 0.23};
  double p8772_d8x1y1_ystatminus[] = { 0.02, 0.03, 0.03 };
  double p8772_d8x1y1_ystatplus[] = { 0.02, 0.03, 0.03 };
  //
  //
  // Rside
  double p8772_d8x1y2_yval[] = { 1.08, 1.14, 1.33};
  double p8772_d8x1y2_yerrminus[] = { 0.25, 0.23, 0.31};
  double p8772_d8x1y2_yerrplus[] = { 0.25, 0.23, 0.31 };
  double p8772_d8x1y2_ystatminus[] = { 0.04, 0.03, 0.03 };
  double p8772_d8x1y2_ystatplus[] = { 0.04, 0.03, 0.03 };
  //
  //
  // Rlong
  double p8772_d8x1y3_yval[] = { 1.10, 1.01, 1.37};
  double p8772_d8x1y3_yerrminus[] = { 0.19, 0.19, 0.19};
  double p8772_d8x1y3_yerrplus[] = { 0.19, 0.19, 0.19 };
  double p8772_d8x1y3_ystatminus[] = { 0.0, 0.0, 0.0 };
  double p8772_d8x1y3_ystatplus[] = { 0.0, 0.0, 0.0 };
  //
  // new TGraphErrors
  //
  Double_t xP[3];
  Double_t yP[3];
  Double_t xPerr[3];
  Double_t yPerr[3];
  //
  for (Int_t iP =0; iP <3; iP++) {
    xP[iP] = p8772_d8x1y1_xval[iP];
    xPerr[iP] = p8772_d8x1y1_xerrminus[iP];
    yP[iP] = TMath::Power(p8772_d8x1y2_yval[iP]*p8772_d8x1y2_yval[iP]*p8772_d8x1y1_yval[iP], 1./3.);
    Double_t errSyst = TMath::Power(p8772_d8x1y2_yerrplus[iP]*p8772_d8x1y2_yerrplus[iP]*p8772_d8x1y1_yerrplus[iP], 1./3.);
    Double_t errStat = TMath::Power(p8772_d8x1y2_ystatplus[iP]*p8772_d8x1y2_ystatplus[iP]*p8772_d8x1y1_ystatplus[iP], 1./3.);
    yPerr[iP] = TMath::Sqrt(errSyst*errSyst + errStat*errStat);;
  }
  TGraphErrors * grHbtRadiusPPB = new TGraphErrors(3, xP, yP, xPerr, yPerr);
  grHbtRadiusPPB->SetName("grHbtRadiusPPB");
  grHbtRadiusPPB->SetTitle("#pi p-Pb 5.02 TeV #it{k}_{T}=0.887 GeV/it{c} #sqrt[3]{R_{out}R_{side}R_{long}}");
  //
  grHbtRadiusPPB->SetMarkerSize(1.0);
  grHbtRadiusPPB->SetMarkerStyle(20);
  grHbtRadiusPPB->SetMarkerColor(kBlue);
  grHbtRadiusPPB->SetLineColor(kBlue);
  //
  return grHbtRadiusPPB;
}

