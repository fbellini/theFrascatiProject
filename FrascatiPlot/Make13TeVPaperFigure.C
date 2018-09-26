#include "./B2vsVolume.C"

void convertRadiusToMulti(TGraphAsymmErrors * graph = 0x0, Int_t paramSet = 1);
void convertRadiusToMulti(TGraphErrors * graph = 0x0, Int_t paramSet = 1);
void getMultiFromR(Double_t * multi = NULL, Double_t * radius = NULL, Int_t paramSet = 1);
void Make13TeVPaperFigure(Bool_t plotLinX = 0, Double_t pToA = 0.75, Int_t RmappingParam = 1);

  
void Make13TeVPaperFigure(Bool_t plotLinX, Double_t pToA, Int_t RmappingParam)
{

  Int_t ip = RmappingParam;

  //----------------------------
  //Get data
  //----------------------------
  const Int_t nParamSet = 2;
  TGraphErrors* gB2vsR_pp7TeV[nParamSet];
  TGraphErrors* gB2vsR_pp7TeV_sys[nParamSet];
  TGraphErrors* gB2vsR_pp13TeV[nParamSet];
  TGraphErrors* gB2vsR_pp13TeV_sys[nParamSet];
  TGraphErrors* gB2vsR_pPb502TeV[nParamSet];
  TGraphErrors* gB2vsR_pPb502TeV_sys[nParamSet];
  TGraphErrors* gB2vsR_PbPb276TeV[nParamSet];
  TGraphErrors* gB2vsR_PbPb276TeV_sys[nParamSet];

  for (Int_t jj = 0; jj < nParamSet; jj++){
    gB2vsR_pp7TeV[jj] = (TGraphErrors *) getB2_pp7TeV(kFALSE, pToA, jj, kFALSE);
    gB2vsR_pp7TeV_sys[jj] = (TGraphErrors *) getB2_pp7TeV(kTRUE, pToA, jj, kFALSE);
    gB2vsR_pp13TeV[jj] = (TGraphErrors *) getB2_pp13TeV(kFALSE, pToA, jj, kFALSE);
    gB2vsR_pp13TeV_sys[jj] = (TGraphErrors *) getB2_pp13TeV(kTRUE, pToA, jj, kFALSE);
    gB2vsR_pPb502TeV[jj] = (TGraphErrors *) getB2_pPb5TeV(kFALSE, pToA, jj, kFALSE);
    gB2vsR_pPb502TeV_sys[jj] = (TGraphErrors *) getB2_pPb5TeV(kTRUE, pToA, jj, kFALSE);
    gB2vsR_PbPb276TeV[jj] = (TGraphErrors *) getB2_PbPb276TeV(kFALSE, pToA, jj, kFALSE);
    gB2vsR_PbPb276TeV_sys[jj] = (TGraphErrors *) getB2_PbPb276TeV(kTRUE, pToA, jj, kFALSE);
  }

  //----------------------------
  //Get coalescence prediction
  //----------------------------
  Double_t mT = TMath::Sqrt(pToA * pToA + 0.938 * 0.938);
  TGraphErrors* hB2_coalescence = (TGraphErrors*) MakeB2TheoryGraphCoalescence(mT);
  hB2_coalescence->SetMarkerStyle(20);
  hB2_coalescence->SetLineWidth(3);
  hB2_coalescence->SetLineStyle(1);

  TGraphErrors * hB2_coalescenceParam0 = (TGraphErrors*) hB2_coalescence->Clone("hB2_coalescenceParam0");
  convertRadiusToMulti(hB2_coalescenceParam0, 0);
  hB2_coalescenceParam0->SetLineStyle(2);
  hB2_coalescenceParam0->SetLineColor(kRed+2);
  hB2_coalescenceParam0->SetMarkerColor(kRed+2);

  TGraphErrors * hB2_coalescenceParam1 = (TGraphErrors*) hB2_coalescence->Clone("hB2_coalescenceParam1");
  convertRadiusToMulti(hB2_coalescenceParam1, 1);

  //----------------------------
  //plot
  //----------------------------
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadRightMargin(0.02); 
  
  TH2D * hframe = new TH2D("hframe", "B_{2} vs radius; #it{R} (fm); #it{B}_{2} (GeV^{2}/#it{c}^{3})", 1000, 0.01, 6.0, 2000, 1.e-4, 0.1);
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
  TLegend * legB2data = new TLegend(0.45, 0.9-nl*0.04, 0.85, 0.9);
  legB2data->SetFillStyle(0);
  legB2data->SetTextSize(0.035);
  legB2data->SetBorderSize(0);
  legB2data->AddEntry(gB2vsR_PbPb276TeV_sys[ip], "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "pf");
  legB2data->AddEntry(gB2vsR_pPb502TeV_sys[ip], "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, prelim.", "pf");
  legB2data->AddEntry(gB2vsR_pp7TeV_sys[ip], "pp #sqrt{#it{s}} = 7 TeV, paper in prep.", "pf");
  legB2data->AddEntry(gB2vsR_pp13TeV_sys[ip], "pp #sqrt{#it{s}} = 13 TeV, prelim.", "pf");
  legB2data->AddEntry(hB2_coalescence, "#it{B}_{2} coalesc., #it{r}(d) = 3.2 fm", "l");

  TCanvas * cb2opta = new TCanvas("B2vsR", "B2 vs R", 900, 800);
  cb2opta->SetBottomMargin(0.15);
  cb2opta->SetTopMargin(0.05);
  cb2opta->SetLeftMargin(0.17);
  cb2opta->SetRightMargin(0.05);
  
  cb2opta->cd();
  gPad->SetLogy();
  pavept->Draw();
  if (!plotLinX) gPad->SetLogx();
  hframe->Draw();
  hB2_coalescence->Draw("l");
  
  gB2vsR_PbPb276TeV_sys[ip]->Draw("p3");
  gB2vsR_PbPb276TeV[ip]->Draw("samep");

  gB2vsR_pPb502TeV_sys[ip]->Draw("p3");
  gB2vsR_pPb502TeV[ip]->Draw("samep");

  gB2vsR_pp7TeV_sys[ip]->Draw("p3");
  gB2vsR_pp7TeV[ip]->Draw("samep");
  
  gB2vsR_pp13TeV_sys[ip]->Draw("p3");
  gB2vsR_pp13TeV[ip]->Draw("samep");

  pavept->Draw();
  legB2data->Draw();
  
  //Alternative version -- all vs dN/deta
  //Map radius in coalescence to dN/deta
  TH2D * hframeMult = new TH2D("hframeMult", "B_{2} vs mult; #LTd#it{N}_{ch}/d#it{#eta}#GT; #it{B}_{2} (GeV^{2}/#it{c}^{3})", 3000, 0., 3000.0, 2000, 1.e-4, 0.1);
  hframeMult->GetXaxis()->SetTitleSize(0.06);
  hframeMult->GetYaxis()->SetTitleSize(0.06);
  hframeMult->GetXaxis()->SetTitleOffset(1.);
  hframeMult->GetXaxis()->SetLabelSize(0.05);
  hframeMult->GetYaxis()->SetLabelSize(0.05);
  hframeMult->GetXaxis()->SetRangeUser(1., 3.E3);
  
  TCanvas * cb2vsdNdeta = new TCanvas("cb2vsdNdeta", "B2 vs R", 900, 800);
  cb2vsdNdeta->SetBottomMargin(0.15);
  cb2vsdNdeta->SetTopMargin(0.05);
  cb2vsdNdeta->SetLeftMargin(0.17);
  cb2vsdNdeta->SetRightMargin(0.05);

  TLegend * legB2coal = new TLegend(0.45, 0.75, 0.65, 0.9,"#it{B}_{2} coalesc., #it{r}(d) = 3.2 fm");
  legB2coal->SetFillStyle(0);
  legB2coal->SetTextSize(0.035);
  legB2coal->SetBorderSize(0);
  legB2coal->AddEntry(hB2_coalescenceParam0, "param. A (fit to HBT radii)", "l");
  legB2coal->AddEntry(hB2_coalescenceParam1, "param. B (constrained to ALICE #it{B}_{2})", "l");

  TLegend * legB2dataMult = new TLegend(0.2, 0.25, 0.45, 0.25+4*0.05,"ALICE");
  legB2dataMult->SetFillStyle(0);
  legB2dataMult->SetTextSize(0.035);
  legB2dataMult->SetBorderSize(0);
  legB2dataMult->AddEntry(gB2vsR_PbPb276TeV_sys[ip], "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "pf");
  legB2dataMult->AddEntry(gB2vsR_pPb502TeV_sys[ip], "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, prelim.", "pf");
  legB2dataMult->AddEntry(gB2vsR_pp7TeV_sys[ip], "pp #sqrt{#it{s}} = 7 TeV, paper in prep.", "pf");
  legB2dataMult->AddEntry(gB2vsR_pp13TeV_sys[ip], "pp #sqrt{#it{s}} = 13 TeV, prelim.", "pf");

  
  cb2vsdNdeta->cd();
  gPad->SetLogy();
  gPad->SetLogx();
  pavept->Draw();
  if (!plotLinX) gPad->SetLogx();
  
  hframeMult->Draw();
  hB2_coalescenceParam0->Draw("l");
  hB2_coalescenceParam1->Draw("l");
 
  gB2vsR_PbPb276TeV_sys[ip]->Draw("p3");
  gB2vsR_PbPb276TeV[ip]->Draw("samep");

  gB2vsR_pPb502TeV_sys[ip]->Draw("p3");
  gB2vsR_pPb502TeV[ip]->Draw("samep");

  gB2vsR_pp7TeV_sys[ip]->Draw("p3");
  gB2vsR_pp7TeV[ip]->Draw("samep");
  
  gB2vsR_pp13TeV_sys[ip]->Draw("p3");
  gB2vsR_pp13TeV[ip]->Draw("samep");

  pavept->Draw();
  legB2dataMult->Draw();
  legB2coal->Draw();

  cb2vsdNdeta->SaveAs("B2vsMult_w13TeV.pdf");
  //  TString foutName = Form("B2vsR_w13TeV_param%i.root", RmappingParam);
  TString foutName = Form("B2vsMult_w13TeV.root", RmappingParam);
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
