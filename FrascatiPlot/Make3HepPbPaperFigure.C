#include "./B2vsVolume.C"

void convertRadiusToMulti(TGraphAsymmErrors * graph = 0x0, Int_t paramSet = 1);
void convertRadiusToMulti(TGraphErrors * graph = 0x0, Int_t paramSet = 1);
void getMultiFromR(Double_t * multi = NULL, Double_t * radius = NULL, Int_t paramSet = 1);
void Make3HepPbPaperFigure(Bool_t plotLinX = 0, Double_t pToA = 0.90, Int_t RmappingParam = 1);

  
void Make3HepPbPaperFigure(Bool_t plotLinX, Double_t pToA, Int_t RmappingParam)
{

  Int_t ip = RmappingParam;

  //----------------------------
  //Get data
  //----------------------------
  const Int_t nParamSet = 2;
  TGraphErrors* gB3vsR_pp7TeV[nParamSet];
  TGraphErrors* gB3vsR_pp7TeV_sys[nParamSet];
  // TGraphErrors* gB3vsR_pp13TeV[nParamSet];
  // TGraphErrors* gB3vsR_pp13TeV_sys[nParamSet];
  TGraphErrors* gB3vsR_pPb502TeV[nParamSet];
  TGraphErrors* gB3vsR_pPb502TeV_sys[nParamSet];
  TGraphErrors* gB3vsR_PbPb276TeV[nParamSet];
  TGraphErrors* gB3vsR_PbPb276TeV_sys[nParamSet];
  TGraphErrors* gB3vsR_PbPb502TeV[nParamSet];
  TGraphErrors* gB3vsR_PbPb502TeV_sys[nParamSet];

  //Thermal + blast-wave
  TGraphErrors* gB3blastvsR_PbPb276TeV[nParamSet];
  TGraphErrors* gB3blastvsR_PbPb276TeV_sys[nParamSet];
  
  for (Int_t jj = 0; jj < nParamSet; jj++){
    gB3vsR_pp7TeV[jj] = (TGraphErrors *) getB3_pp7TeVINELg0(kFALSE, 0.8, jj, kFALSE);
    gB3vsR_pp7TeV_sys[jj] = (TGraphErrors *) getB3_pp7TeVINELg0(kTRUE, 0.8, jj, kFALSE);
    // gB3vsR_pp13TeV[jj] = (TGraphErrors *) getB3_pp13TeV(kFALSE, pToA, jj, kFALSE);
    // gB3vsR_pp13TeV_sys[jj] = (TGraphErrors *) getB3_pp13TeV(kTRUE, pToA, jj, kFALSE);
    gB3vsR_pPb502TeV[jj] = (TGraphErrors *) getB3_pPb5TeV(kFALSE, pToA, jj, kFALSE);
    gB3vsR_pPb502TeV_sys[jj] = (TGraphErrors *) getB3_pPb5TeV(kTRUE, pToA, jj, kFALSE);
    gB3vsR_PbPb276TeV[jj] = (TGraphErrors *) getB3_PbPb276TeV(kFALSE, pToA, jj, kFALSE);
    gB3vsR_PbPb276TeV_sys[jj] = (TGraphErrors *) getB3_PbPb276TeV(kTRUE, pToA, jj, kFALSE);
    gB3vsR_PbPb502TeV[jj] = (TGraphErrors *) getB3_PbPb5TeV(kFALSE, pToA, jj, kFALSE);
    gB3vsR_PbPb502TeV_sys[jj] = (TGraphErrors *) getB3_PbPb5TeV(kTRUE, pToA, jj, kFALSE);
  }

  //----------------------------
  //Get thermal + Blast-wave prediction
  //----------------------------
  gB3blastvsR_PbPb276TeV[0] = (TGraphErrors *) getBlastB3_PbPb276TeV(kFALSE, pToA, 0, kFALSE);
  gB3blastvsR_PbPb276TeV[0]->SetLineStyle(5);
  gB3blastvsR_PbPb276TeV[0]->SetLineWidth(4);
  gB3blastvsR_PbPb276TeV[0]->SetLineColor(kRed+2);

  //----------------------------
  //Get coalescence prediction
  //----------------------------
  Double_t mT = TMath::Sqrt(pToA * pToA + 0.938 * 0.938);
  TGraphErrors* hB3_coalescence = (TGraphErrors*) MakeB3TheoryGraphCoalescence(mT);
  hB3_coalescence->SetMarkerStyle(20);
  hB3_coalescence->SetLineWidth(3);
  hB3_coalescence->SetLineStyle(1);

  TGraphErrors * hB3_coalescenceParam0 = (TGraphErrors*) hB3_coalescence->Clone("hB3_coalescenceParam0");
  convertRadiusToMulti(hB3_coalescenceParam0, 0);
  hB3_coalescenceParam0->SetLineStyle(2);
  // hB3_coalescenceParam0->SetLineColor(kRed+2);
  // hB3_coalescenceParam0->SetMarkerColor(kRed+2);

  TGraphErrors * hB3_coalescenceParam1 = (TGraphErrors*) hB3_coalescence->Clone("hB3_coalescenceParam1");
  convertRadiusToMulti(hB3_coalescenceParam1, 1);

  //----------------------------
  //plot
  //----------------------------
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadRightMargin(0.02); 
  
  TH2D * hframe = new TH2D("hframe", "B_{3} vs radius; #it{R} (fm); #it{B}_{3} (GeV^{4}/#it{c}^{6})", 1000, 0.01, 6.0, 2000, 1.e-9, 0.1);
  hframe->GetXaxis()->SetTitleSize(0.06);
  hframe->GetYaxis()->SetTitleSize(0.06);
  hframe->GetXaxis()->SetTitleOffset(0.8);
  hframe->GetXaxis()->SetLabelSize(0.05);
  hframe->GetYaxis()->SetLabelSize(0.05);
  if (plotLinX) hframe->GetXaxis()->SetRangeUser(0.01, 8.5);
  else  hframe->GetXaxis()->SetRangeUser(0.1, 10.5);
  
  //Define pT/A labels only once
  TPaveText * pavept = new TPaveText(0.17, 0.17, 0.7, 0.23, "NDC");
  pavept->SetFillStyle(0);
  pavept->SetTextFont(42);
  pavept->SetBorderSize(0);
  pavept->SetTextSize(0.05);
  pavept->SetTextAlign(12);
  pavept->AddText(Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToA));

  Int_t nl = 6;
  TLegend * legB3data = new TLegend(0.45, 0.9-nl*0.04, 0.85, 0.9);
  legB3data->SetFillStyle(0);
  legB3data->SetTextSize(0.035);
  legB3data->SetBorderSize(0);
  legB3data->AddEntry(gB3vsR_PbPb276TeV_sys[ip], "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "pf");
  legB3data->AddEntry(gB3vsR_pPb502TeV_sys[ip], "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, prelim.", "pf");
  legB3data->AddEntry(gB3vsR_pp7TeV_sys[ip], "pp #sqrt{#it{s}} = 7 TeV, paper in prep.", "pf");
  //  legB3data->AddEntry(gB3vsR_pp13TeV_sys[ip], "pp #sqrt{#it{s}} = 13 TeV, prelim.", "pf");
  legB3data->AddEntry(hB3_coalescence, "#it{B}_{3} coalesc., #it{r}(^{3}He) = 2.48 fm", "l");

  TCanvas * cb2opta = new TCanvas("B3vsR", "B3 vs R", 900, 800);
  cb2opta->SetBottomMargin(0.15);
  cb2opta->SetTopMargin(0.05);
  cb2opta->SetLeftMargin(0.17);
  cb2opta->SetRightMargin(0.05);
  
  cb2opta->cd();
  gPad->SetLogy();
  pavept->Draw();
  if (!plotLinX) gPad->SetLogx();
  hframe->Draw();
  hB3_coalescence->Draw("l");
  
  gB3vsR_PbPb276TeV_sys[ip]->Draw("p3");
  gB3vsR_PbPb276TeV[ip]->Draw("samep");

  gB3vsR_pPb502TeV_sys[ip]->Draw("p3");
  gB3vsR_pPb502TeV[ip]->Draw("samep");

  gB3vsR_pp7TeV_sys[ip]->Draw("p3");
  gB3vsR_pp7TeV[ip]->Draw("samep");
  
  // gB3vsR_pp13TeV_sys[ip]->Draw("p3");
  // gB3vsR_pp13TeV[ip]->Draw("samep");

  pavept->Draw();
  legB3data->Draw();
  
  //Alternative version -- all vs dN/deta
  //Map radius in coalescence to dN/deta
  TH2D * hframeMult = new TH2D("hframeMult", "B_{3} vs mult; #LTd#it{N}_{ch}/d#it{#eta}#GT; #it{B}_{3} (GeV^{4}/#it{c}^{6})", 3000, 0., 3000.0, 2000, 1.e-9, 0.1);
  hframeMult->GetXaxis()->SetTitleSize(0.06);
  hframeMult->GetYaxis()->SetTitleSize(0.06);
  hframeMult->GetXaxis()->SetTitleOffset(1.);
  hframeMult->GetXaxis()->SetLabelSize(0.05);
  hframeMult->GetYaxis()->SetLabelSize(0.05);
  hframeMult->GetXaxis()->SetRangeUser(1., 3.E3);
  
  TCanvas * cb2vsdNdeta = new TCanvas("cb3vsdNdeta", "B3 vs R", 900, 800);
  cb2vsdNdeta->SetBottomMargin(0.15);
  cb2vsdNdeta->SetTopMargin(0.05);
  cb2vsdNdeta->SetLeftMargin(0.17);
  cb2vsdNdeta->SetRightMargin(0.05);

  TLegend * legB3coal = new TLegend(0.45, 0.8, 0.65, 0.95,"#it{B}_{3} coalesc., #it{r}(^{3}He) = 2.48 fm");
  legB3coal->SetFillStyle(0);
  legB3coal->SetTextSize(0.035);
  legB3coal->SetBorderSize(0);
  //legB3coal->AddEntry(hB3_coalescenceParam0, "param. A (fit to HBT radii)", "l");
  legB3coal->AddEntry(hB3_coalescenceParam1, "param. B (constrained to ALICE #it{B}_{2})", "l");

  TLegend * legB3therm = new TLegend(0.45, 0.7, 0.65, 0.8,"GSI-Heid.+ Blast wave (#piKp)");
  legB3therm->SetFillStyle(0);
  legB3therm->SetTextSize(0.035);
  legB3therm->SetBorderSize(0);
  legB3therm->AddEntry(gB3blastvsR_PbPb276TeV[0], "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "l");
 
  TLegend * legB3dataMult = new TLegend(0.2, 0.25, 0.45, 0.25+5*0.05,"ALICE");
  legB3dataMult->SetFillStyle(0);
  legB3dataMult->SetTextSize(0.035);
  legB3dataMult->SetBorderSize(0);
  legB3dataMult->AddEntry(gB3vsR_PbPb276TeV_sys[ip], "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "pf");
  legB3dataMult->AddEntry(gB3vsR_PbPb502TeV_sys[ip], "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, preliminary", "pf");
  legB3dataMult->AddEntry(gB3vsR_pPb502TeV_sys[ip], "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, this analysis", "pf");
  legB3dataMult->AddEntry(gB3vsR_pp7TeV_sys[ip], "pp #sqrt{#it{s}} = 7 TeV (#it{p}_{T}/#it{A} = 0.8 GeV/#it{c})", "pf");
  //  legB3dataMult->AddEntry(gB3vsR_pp13TeV_sys[ip], "pp #sqrt{#it{s}} = 13 TeV, prelim.", "pf");

  
  cb2vsdNdeta->cd();
  gPad->SetLogy();
  gPad->SetLogx();
  pavept->Draw();
  if (!plotLinX) gPad->SetLogx();
  
  hframeMult->Draw();
  hB3_coalescenceParam0->Draw("l");
  hB3_coalescenceParam1->Draw("l");
  gB3blastvsR_PbPb276TeV[0]->Draw("lsame");
  
  gB3vsR_PbPb502TeV_sys[ip]->Draw("p3");
  gB3vsR_PbPb502TeV[ip]->Draw("samep");

  gB3vsR_PbPb276TeV_sys[ip]->Draw("p3");
  gB3vsR_PbPb276TeV[ip]->Draw("samep");

  gB3vsR_pPb502TeV_sys[ip]->Draw("p3");
  gB3vsR_pPb502TeV[ip]->Draw("samep");

  gB3vsR_pp7TeV_sys[ip]->Draw("p3");
  gB3vsR_pp7TeV[ip]->Draw("samep");
  
  // gB3vsR_pp13TeV_sys[ip]->Draw("p3");
  // gB3vsR_pp13TeV[ip]->Draw("samep");

  pavept->Draw();
  legB3dataMult->Draw();
  legB3coal->Draw();
  legB3therm->Draw();
  
  cb2vsdNdeta->SaveAs("B3vsMult.pdf");
  TString foutName = Form("B3vsMult_R%i.root", RmappingParam);
  TFile * fout = new TFile(foutName.Data(), "recreate");
  fout->cd();
  cb2vsdNdeta->Write();
  fout->Close();
  return;
  
}

void getMultiFromR(Double_t * multi, Double_t * radius, Int_t paramSet)
{
  // Here is the crucial mapping between HBT radii and multi^(1/3)
  if (!multi || !radius) return;
  //
  // VERSION (17th May 2018):
  // We fit linearly the ALICE data at the kT = 0.887 
  Double_t radiusVal = radius[0];
  Double_t  multi3 = 0.0;
  if (paramSet==2) {
    //manual hack to have the data points fall onto the U. Heinz curve for 3He
    multi3 = (radiusVal - 0.190) /  0.380;
  } else  if (paramSet==1) {
    //manual hack to have the data points fall onto the U. Heinz curve for deuteron
    multi3 = radiusVal / 0.472949; //radius = 0 at 0 dN/deta
    //radiusVal = 0.07412 + 0.46637 * multi3; //radius = 0.85 fm for pp, dN/deta = 4.60 (INELg0)
    //radiusVal = -0.009 + 0.4738 * multi3; //radius = 0.85 fm for pp, dN/deta = 5.98 (INELg0>0)
    //radiusVal = -0.3949 + 0.507865 * multi3; //most central and peripheral PbPb 2.76 TeV on the curve
  } else {
    //fit to the HBT data, kT = 0.887
    multi3 = (radiusVal - 0.128) / 0.339; 
  }
  multi[0] = TMath::Power(multi3, 3);
  multi[1] = 0.0; //error on multi set to 0
  return; 
}

void convertRadiusToMulti(TGraphErrors * graph, Int_t paramSet)
{
  //convert multiplicity dN/deta into radius with parameterisation
  if (!graph) return;
  
  for (Int_t ip = 0; ip < graph->GetN(); ip++){
    Double_t xold[2];
    xold[0] = graph->GetX()[ip];
    xold[1] = graph->GetEX()[ip];
    Double_t xnew[2] = {-1.0, -1.0};
    getMultiFromR(xnew, xold, paramSet);
    graph->GetX()[ip] = xnew[0];
    graph->GetEX()[ip] = xnew[1];
  }

  return;
}


void convertRadiusToMulti(TGraphAsymmErrors * graph, Int_t paramSet)
{
  //convert multiplicity dN/deta into radius with parameterisation
  if (!graph) return;
  
  for (Int_t ip = 0; ip < graph->GetN(); ip++){
    Double_t xold[2];
    xold[0] = graph->GetX()[ip];
    xold[1] = graph->GetEXlow()[ip];
    Double_t xnew[2] = {-1.0, -1.0};
    getMultiFromR(xnew, xold, paramSet);
    graph->GetX()[ip] = xnew[0];
    graph->GetEXlow()[ip] = xnew[1];
    graph->GetEXhigh()[ip] = xnew[1];
	
  }

  return;
}


Double_t GetB3ratioTritiumToHe(Double_t Rsource = 0.77) // fm pp 7 TeV INEL HBT param A
{

  Double_t R3He = 2.48; //fm
  Double_t Rt = 2.15; //fm
  return TMath::Power( (Rsource*Rsource + R3He*R3He /4.) / (Rsource*Rsource + Rt*Rt /4.) , 3.);

}