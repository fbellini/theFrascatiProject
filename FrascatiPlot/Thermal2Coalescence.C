#include "CoalescenceBA.C"
#include "B2vsVolume.C"

TH1D * CreateFrameSeparation();
TH1D * CreateFrameUncertainty();
void SetCanvasStyle();
TLegend * BuildLegendSpecies();
TLegend * BuildLegendModels();


Double_t getBAfunc(Double_t *x, Double_t *par);
TGraphAsymmErrors * getBAthermalBlast(TString particle, Double_t pToA, Color_t color);

void Thermal2Coalescence()
{

  // d, 3H, 3He, 3LH, 4He, 4LH 
  Char_t particle[6][15] = {"deuteron", "triton", "He3", "hyper-triton", "He4", "4LH"};
  Double_t nucleiA[]     = {   2,   3,     3,    3,    4,    4};
  Double_t spin[]        = { 1.0, 0.5,   0.5,  0.5,  0.0,  0.0};
  Double_t objRadius[]   = { 3.2, 2.15, 2.48,  6.8,  1.9,  2.4};
  Double_t pToA[]        = {0.75, 0.77, 0.77, 1.17, 0.75, 0.62};
  Color_t color[]        = {kBlack, kBlue+1, kAzure+8, kGreen+2, kRed+1, kOrange+8};
  
  TCanvas * cr4 = new TCanvas("cr4", "compare thermal with coalescence", 900, 2000);
  SetCanvasStyle();
  TH2D * frame = CreateFrameBA();
  frame->GetYaxis()->SetTitleOffset(0.8);
  frame->GetYaxis()->SetRangeUser(1e-12, 1);
  TH1D * frameSep = CreateFrameSeparation();
  TH1D * frameUnc = CreateFrameUncertainty();
  cr4->Divide(1,3);
  cr4->cd(1); gPad->SetLogy(); gPad->SetTicky();  gPad->SetTickx();
  frame->Draw();
  cr4->cd(3); gPad->SetLogy(); gPad->SetTicky();  gPad->SetTickx();
  frameSep->Draw();
  cr4->cd(2); gPad->SetLogy(); gPad->SetTicky();  gPad->SetTickx();
  frameUnc->Draw();

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
    gThermalBlast[j] = getBAthermalBlast(particle[j], pToA[j], color[j]);
    
    for (int ip=0; ip<10; ip++) {
      
      radius[j][ip] = gThermalBlast[j]->GetX()[ip];
      radiuserr[j][ip] = 0.0;      
      //      if (j == 0) Printf("ip = %i     radius = %3.2f", ip, radius[j][ip]);
      
      //calculate relative uncertainty on BA from pseudodata
      if (hist[j]->GetBinContent(ip+1)>0)
	reluncert[j][ip] = hist[j]->GetBinError(ip+1)/hist[j]->GetBinContent(ip+1);
      else reluncert[j][ip] = 0.0;
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
    cr4->cd(1);
    fBAcoalescence[j]->Draw("same");
    gThermalBlast[j]->Draw("samel");
    cr4->cd(2);
    gRelativeUncBA[j]->Draw("samelp");
    cr4->cd(3);
    if (j<5)  gSeparation[j]->Draw("samel");
    // if (j<5)  gSeparationMaxSys[j]->Draw("samel");
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

  cr4->cd(1);
  modelLeg->Draw();
  cr4->cd(2);
  speciesLeg->Draw();
  return;
}


Double_t getBAfunc(Double_t *x, Double_t *par)
{
  //parameters are
  /*
    par[0] = A;
    par[1] = JA;
    par[2] = pT/A;
    par[3] = object size;
  */
  return getBAfromRadius(par[0], par[1], par[2], x[0], par[3]);
    
}



TGraphAsymmErrors * getBAthermalBlast(TString particle, Double_t pToA, Color_t color)
{
  TGraphAsymmErrors* graph = (TGraphAsymmErrors *) generateBWpredictionsB2("PbPb502TeV", "rms", particle.Data(), pToA);
  //parameterisation of the radius to have data points fall onto the U. Heinz curve for deuteron
  convertMultiToRadius(graph, 1);
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
  hframe4->GetYaxis()->SetRangeUser(1., 1e4);
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
  hframe4->GetYaxis()->SetRangeUser(5E-5, 5.);
  return hframe4;
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
