#include "B2vsVolume.C"

TGraphErrors * GetPPbHbtRadius025KT();
TGraphErrors * GetPPbHbtRadius0887KT();
TGraphErrors * GetPbPbHbtRadius0887KT();


void HbtRadii() {
  //
  // small macro to plot HBT radii versus multiplicity to compare
  // them to the parameterisation
  //
  TCanvas * canvHbtPlot = new TCanvas("canvHbtPlot","canvHbtPlot");
  TH2D * histFrame = new TH2D("histFrame","radius vs mult; #LTd#it{N}/d#it{#eta}#GT^{1/3}; R (fm)",200,0.0,13.0,200,0.0,6.0);
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
  TGraph * grRadiusVsMult = new TGraph(nPoints, xP, yP);
  grRadiusVsMult->SetMarkerStyle(20);
  grRadiusVsMult->SetMarkerColor(kAzure-7);
  grRadiusVsMult->SetLineColor(kAzure-7);
  grRadiusVsMult->SetLineWidth(3);
  grRadiusVsMult->SetLineStyle(9);
  //
  grRadiusVsMult->Draw("L");
  //
  // add PPB points
  //
  TGraphErrors * grHbtRadiusPPB = GetPPbHbtRadius0887KT();
  grHbtRadiusPPB->Draw("P");
  //
  // add PbPb points
  //
  TGraphErrors * grHbtRadiusPBPB = GetPbPbHbtRadius0887KT();
  grHbtRadiusPBPB->Draw("P");
  //
  // add a pp point
  //http://aliceinfo.cern.ch/ArtSubmission/node/1885, value is from Kfir
  TGraphErrors * grHbtRadiusPP = new TGraphErrors(1);
  grHbtRadiusPP->SetPoint(0, TMath::Power(5.98,1./3.), 0.8);
  grHbtRadiusPP->SetPointError(0, TMath::Power(5.98,-2./3.)*0.09/3., 0.3);
  grHbtRadiusPP->SetLineColor(kGreen+1);
  grHbtRadiusPP->SetMarkerColor(kGreen+1);
  grHbtRadiusPP->SetMarkerStyle(20);
  grHbtRadiusPP->Draw("P");
  //
  
  canvHbtPlot->BuildLegend();

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
  for (Int_t iP =0; iP <7; iP++) {
    xP[iP] = TMath::Power(p8772_d8x1y1_xval[iP], 1./3.);
    xPerr[iP] =  TMath::Power(p8772_d8x1y1_xval[iP], -2./3.) * p8772_d8x1y1_xerrplus[iP]/3.;
    yP[iP] = TMath::Power(p8772_d8x1y3_yval[iP]*p8772_d8x1y2_yval[iP]*p8772_d8x1y1_yval[iP], 1./3.);
    Double_t errSyst = TMath::Power(p8772_d8x1y3_yerrplus[iP]*p8772_d8x1y2_yerrplus[iP]*p8772_d8x1y1_yerrplus[iP], 1./3.);
    Double_t errStat = TMath::Power(p8772_d8x1y3_ystatplus[iP]*p8772_d8x1y2_ystatplus[iP]*p8772_d8x1y1_ystatplus[iP], 1./3.);
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
    yP[iP] = TMath::Power(p8772_d8x1y3_yval[iP]*p8772_d8x1y2_yval[iP]*p8772_d8x1y1_yval[iP], 1./3.);
    yPerr[iP] = TMath::Power(p8772_d8x1y3_yerrplus[iP]*p8772_d8x1y2_yerrplus[iP]*p8772_d8x1y1_yerrplus[iP], 1./3.);
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
  double p8772_d8x1y1_xerrminus[] = { 0.0, 0.0, 0.0 };
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
    xPerr[iP] = p8772_d8x1y1_xerrplus[iP];
    yP[iP] = TMath::Power(p8772_d8x1y3_yval[iP]*p8772_d8x1y2_yval[iP]*p8772_d8x1y1_yval[iP], 1./3.);
    Double_t errSyst = TMath::Power(p8772_d8x1y3_yerrplus[iP]*p8772_d8x1y2_yerrplus[iP]*p8772_d8x1y1_yerrplus[iP], 1./3.);
    Double_t errStat = TMath::Power(p8772_d8x1y3_ystatplus[iP]*p8772_d8x1y2_ystatplus[iP]*p8772_d8x1y1_ystatplus[iP], 1./3.);
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

