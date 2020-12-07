#include "./B2vsVolume.C"
#include "./GetB3inpToA.C"

TCanvas * MakeLDtalkPlotHe(TString figPath = ".", Double_t pToA = 0.733, Bool_t plotLinX = 0, 
                                Int_t RmappingParam = 1, Bool_t plotcSHM = 0, Bool_t plotPbPb5TeV = 0, 
                                Bool_t plotpPb = 1);

  
TCanvas * MakeLDtalkPlotHe(TString figPath = ".", Double_t pToA, Bool_t plotLinX, Int_t RmappingParam, Bool_t plotcSHM, Bool_t plotPbPb5TeV, Bool_t plotpPb)
{

  //prepare input file
  GetB3inpToA_pPb5TeV(pToA, "B3_Average_pPb_5TeV_30052019.root");

  Int_t ip = RmappingParam;
  Float_t textLabelSize = 0.033;
  Short_t modelCurvesLineWidth = 4;
  //----------------------------
  //Get data
  //----------------------------
  const Int_t nParamSet = 4;
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
  
  for (Int_t jj = 0; jj < nParamSet; jj++){
    gB3vsR_pp7TeV[jj] = (TGraphErrors *) getB3_pp7TeVINELg0(kFALSE, 0.8, jj, kFALSE);
    gB3vsR_pp7TeV_sys[jj] = (TGraphErrors *) getB3_pp7TeVINELg0(kTRUE, 0.8, jj, kFALSE);
    gB3vsR_pp7TeV[jj]->SetMarkerStyle(33);
    gB3vsR_pp7TeV[jj]->SetMarkerSize(1.3);
    gB3vsR_pp7TeV[jj]->SetMarkerColor(kRed+1);
    gB3vsR_pp7TeV[jj]->SetLineColor(kRed+1);
    gB3vsR_pp7TeV_sys[jj]->SetMarkerStyle(33);
    gB3vsR_pp7TeV_sys[jj]->SetMarkerSize(1.3);
    gB3vsR_pp7TeV_sys[jj]->SetFillStyle(1001);    
    gB3vsR_pp7TeV_sys[jj]->SetFillColorAlpha(kRed+1, 0.4);
    // gB3vsR_pp13TeV[jj] = (TGraphErrors *) getB3_pp13TeV(kFALSE, pToA, jj, kFALSE);
    // gB3vsR_pp13TeV_sys[jj] = (TGraphErrors *) getB3_pp13TeV(kTRUE, pToA, jj, kFALSE);
    gB3vsR_pPb502TeV[jj] = (TGraphErrors *) getB3_pPb5TeV(kFALSE, pToA, jj, kFALSE);
    gB3vsR_pPb502TeV_sys[jj] = (TGraphErrors *) getB3_pPb5TeV(kTRUE, pToA, jj, kFALSE);
    gB3vsR_pPb502TeV[jj]->SetMarkerStyle(20);
    gB3vsR_pPb502TeV[jj]->SetMarkerColor(kOrange);
    gB3vsR_pPb502TeV[jj]->SetLineColor(kOrange);
    gB3vsR_pPb502TeV_sys[jj]->SetMarkerStyle(20);
    gB3vsR_pPb502TeV_sys[jj]->SetFillStyle(1001);
    gB3vsR_pPb502TeV_sys[jj]->SetFillColorAlpha(kOrange, 0.4);

    gB3vsR_PbPb276TeV[jj] = (TGraphErrors *) getB3_PbPb276TeV(kFALSE, pToA, jj, kFALSE);
    gB3vsR_PbPb276TeV_sys[jj] = (TGraphErrors *) getB3_PbPb276TeV(kTRUE, pToA, jj, kFALSE);
    gB3vsR_PbPb276TeV[jj]->SetMarkerStyle(21);
    gB3vsR_PbPb276TeV[jj]->SetMarkerColor(kBlack);
    gB3vsR_PbPb276TeV[jj]->SetLineColor(kBlack);
    gB3vsR_PbPb276TeV_sys[jj]->SetMarkerStyle(21);
    gB3vsR_PbPb276TeV_sys[jj]->SetFillStyle(1001);
    gB3vsR_PbPb276TeV_sys[jj]->SetFillColorAlpha(kGray, 0.5);

    gB3vsR_PbPb502TeV[jj] = (TGraphErrors *) getB3_PbPb5TeV(kFALSE, pToA, jj, kFALSE);
    gB3vsR_PbPb502TeV_sys[jj] = (TGraphErrors *) getB3_PbPb5TeV(kTRUE, pToA, jj, kFALSE);
  }

  //----------------------------
  //Get thermal + Blast-wave prediction
  //----------------------------
  Short_t shmBlast_lineStyle = 9;
    //Thermal + blast-wave
  TGraphErrors* gB3blastvsR_PbPb276TeV[nParamSet];
  TGraphErrors* gB3blastvsR_PbPb276TeV_sys[nParamSet];
  TGraphErrors* gB3blastvsR_pPb5TeV[nParamSet];
  for (Int_t jj = 0; jj < nParamSet; jj++){

    gB3blastvsR_PbPb276TeV[jj] = (TGraphErrors *) getBlastB3_PbPb276TeV(kFALSE, pToA, jj, kFALSE);
    gB3blastvsR_PbPb276TeV[jj]->SetLineStyle(shmBlast_lineStyle);
    gB3blastvsR_PbPb276TeV[jj]->SetLineWidth(modelCurvesLineWidth);
    gB3blastvsR_PbPb276TeV[jj]->SetLineColor(kGreen+1);

    gB3blastvsR_pPb5TeV[jj] = (TGraphErrors *) getBlastB3_pPb5TeV(kFALSE, pToA, jj, kFALSE, plotcSHM);
    gB3blastvsR_pPb5TeV[jj]->SetLineStyle(shmBlast_lineStyle);
    gB3blastvsR_pPb5TeV[jj]->SetLineWidth(modelCurvesLineWidth);
    gB3blastvsR_pPb5TeV[jj]->SetLineColor(kBlue);
  }
//----------------------------
  //Get coalescence prediction
  //----------------------------
  Double_t mT = TMath::Sqrt(pToA * pToA + 0.938 * 0.938);
  TGraphErrors* hB3_coalescence = (TGraphErrors*) MakeB3TheoryGraphCoalescence(mT);
  hB3_coalescence->SetMarkerStyle(20);
  hB3_coalescence->SetLineWidth(modelCurvesLineWidth);
  hB3_coalescence->SetLineStyle(1);
  hB3_coalescence->SetLineColor(kBlue);


  TGraphErrors * hB3_coalescenceParam0 = (TGraphErrors*) hB3_coalescence->Clone("hB3_coalescenceParam0");
  convertRadiusToMulti(hB3_coalescenceParam0, 0);
  hB3_coalescenceParam0->SetLineStyle(1);
  hB3_coalescenceParam0->SetLineColor(kBlue);
  hB3_coalescenceParam0->SetLineWidth(modelCurvesLineWidth);
  // hB3_coalescenceParam0->SetLineColor(kRed+2);
  // hB3_coalescenceParam0->SetMarkerColor(kRed+2);

  TGraphErrors * hB3_coalescenceParam1 = (TGraphErrors*) hB3_coalescence->Clone("hB3_coalescenceParam1");
  convertRadiusToMulti(hB3_coalescenceParam1, 1);

  //Plot coalescence curve using same R->multi param. as Doenigus&Ko
  TGraphErrors * hB3_coalescenceParam3 = (TGraphErrors*) hB3_coalescence->Clone("hB3_coalescenceParam3");
  convertRadiusToMulti(hB3_coalescenceParam3, 3);
  hB3_coalescenceParam3->SetLineColor(kOrange);
  hB3_coalescenceParam3->SetLineWidth(3);
  

  Printf("pp 7 TeV  vs multi - B3 3He");
  for (int i = 0; i<gB3vsR_pp7TeV[ip]->GetN(); i++){
    Printf("dN/deta = %5.3f", gB3vsR_pp7TeV[ip]->GetX()[i]);
  } 

  Printf("pPb 5.02 TeV  vs multi - B3 3He");
  for (int i = 0; i<gB3vsR_pPb502TeV[ip]->GetN(); i++){
    Printf("dN/deta = %5.3f", gB3vsR_pPb502TeV[ip]->GetX()[i]);
  } 

    Printf("PbPb 2.76 TeV  vs multi - B3 3He");
  for (int i = 0; i<gB3vsR_PbPb276TeV[ip]->GetN(); i++){
    Printf("dN/deta = %5.3f", gB3vsR_PbPb276TeV[ip]->GetX()[i]);
  } 
  //----------------------------
  //plot
  //----------------------------
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadBottomMargin(0.2);
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadRightMargin(0.02); 
  
  TH2D * hframe = new TH2D("hframe", "B_{3} vs radius; #it{R} (fm); #it{B}_{3} (GeV^{4}/#it{c}^{6})", 1000, 0.01, 6.0, 2000, 1.e-8, 0.1);
  hframe->GetXaxis()->SetTitleSize(0.06);
  hframe->GetYaxis()->SetTitleSize(0.06);
  hframe->GetXaxis()->SetTitleOffset(0.8);
  hframe->GetXaxis()->SetLabelSize(0.05);
  hframe->GetYaxis()->SetLabelSize(0.05);
  if (plotLinX) hframe->GetXaxis()->SetRangeUser(0.01, 8.5);
  else  hframe->GetXaxis()->SetRangeUser(0.1, 10.5);
  
  //Define pT/A labels only once
  TPaveText * pavept = new TPaveText(0.17, 0.17, 0.7, 0.25, "NDC");
  pavept->SetFillStyle(0);
  pavept->SetTextFont(42);
  pavept->SetTextSize(0.045);
  pavept->SetBorderSize(0);
  pavept->SetTextSize(textLabelSize);
  pavept->SetTextAlign(12);
  pavept->AddText(Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToA));

  Int_t nl = 6;
  TLegend * legB3data = new TLegend(0.45, 0.9-nl*0.04, 0.85, 0.9);
  legB3data->SetFillStyle(0);
  legB3data->SetTextSize(0.035);
  legB3data->SetBorderSize(0);
  legB3data->AddEntry(gB3vsR_PbPb276TeV_sys[ip], "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "pf");
  if (plotpPb) legB3data->AddEntry(gB3vsR_pPb502TeV_sys[ip], "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, prelim.", "pf");
  legB3data->AddEntry(gB3vsR_pp7TeV[ip], "pp #sqrt{#it{s}} = 7 TeV, paper in prep.", "pf");
  //  legB3data->AddEntry(gB3vsR_pp13TeV_sys[ip], "pp #sqrt{#it{s}} = 13 TeV, prelim.", "pf");
  legB3data->AddEntry(hB3_coalescence, "#it{B}_{3} coalesc., #it{r}(^{3}He) = 2.48 fm", "l");

  TCanvas * cb2opta = new TCanvas("B3vsR", "B3 vs R", 800, 600);
  cb2opta->SetBottomMargin(0.17);
  cb2opta->SetTopMargin(0.05);
  cb2opta->SetLeftMargin(0.17);
  cb2opta->SetRightMargin(0.05);
  
  cb2opta->cd();
  gPad->SetLogy();
  gPad->SetTickx();
  gPad->SetTicky();
  pavept->Draw();
  if (!plotLinX) gPad->SetLogx();
  hframe->Draw();
  gB3vsR_pp7TeV_sys[ip]->Draw("p3");
  if (plotpPb) gB3vsR_pPb502TeV_sys[ip]->Draw("p3");
  hB3_coalescence->Draw("l");
  gB3vsR_PbPb276TeV_sys[ip]->Draw("p3");
  gB3vsR_PbPb276TeV[ip]->Draw("samep");
  if (plotpPb) gB3vsR_pPb502TeV[ip]->Draw("samep");
  gB3vsR_pp7TeV[ip]->Draw("samep");
  
  // gB3vsR_pp13TeV_sys[ip]->Draw("p3");
  // gB3vsR_pp13TeV[ip]->Draw("samep");
  pavept->Draw();
  legB3data->Draw();
  
  //Alternative version -- all vs dN/deta
  //Map radius in coalescence to dN/deta
  TH2D * hframeMult = new TH2D("hframeMult", "B_{3} vs mult; #LTd#it{N}_{ch}/d#it{#eta}#GT; #it{B}_{3} (GeV^{4}/#it{c}^{6})", 3000, 0., 3000.0, 2000, 1.e-8, 0.01);
  hframeMult->GetXaxis()->SetTitleSize(0.06);
  hframeMult->GetYaxis()->SetTitleSize(0.06);
  hframeMult->GetXaxis()->SetTitleOffset(1.2);
  hframeMult->GetXaxis()->SetLabelSize(0.05);
  hframeMult->GetYaxis()->SetLabelSize(0.05);
  hframeMult->GetYaxis()->SetNdivisions(505);
  hframeMult->GetXaxis()->SetRangeUser(1., 3.E3);
  hframeMult->GetYaxis()->SetRangeUser(7.E-8, 3.E-3);
  
  TCanvas * cb2vsdNdeta = new TCanvas("cb3vsdNdeta", "B3 vs R", 800, 600);
  cb2vsdNdeta->SetBottomMargin(0.17);
  cb2vsdNdeta->SetTopMargin(0.05);
  cb2vsdNdeta->SetLeftMargin(0.17);
  cb2vsdNdeta->SetRightMargin(0.05);

  TLegend * legB3coal = new TLegend(0.5, 0.78, 0.8, 0.85); //"#it{B}_{3} coalesc., #it{r}(^{3}He) = 2.48 fm"
  legB3coal->SetFillStyle(0);
  legB3coal->SetTextSize(textLabelSize);
  legB3coal->SetTextColor(kBlue);
  legB3coal->SetTextSize(0.045);
  legB3coal->SetBorderSize(0);
  if (plotpPb) legB3coal->AddEntry(hB3_coalescenceParam0, "Coalescence, #it{r}_{^{3}He} = 2.48 fm", "l");
  //legB3coal->AddEntry(hB3_coalescenceParam1, "constrained to ALICE #it{B}_{2}", "l");

  TLegend * legB3therm = new TLegend(0.5, 0.85, 0.83, 0.92); //,"SHM + blast-wave (ALICE #piKp)");
  legB3therm->SetFillStyle(0);
  legB3therm->SetTextSize(textLabelSize);
  legB3therm->SetTextColor(kGreen+1);
  legB3therm->SetBorderSize(0);
    legB3therm->SetTextSize(0.045);
  legB3therm->AddEntry(gB3blastvsR_PbPb276TeV[0], "Statistical hadron.", "l"); //"GC GSI-Heid. (T = 156 MeV)", "l");
  //legB3therm->AddEntry(gB3blastvsR_pPb5TeV[0], "CSM (T = 155 MeV)", "l");
 
  TLegend * legB3dataMult = new TLegend(0.2, 0.25, 0.45, 0.25+3*0.08,"ALICE");
  legB3dataMult->SetFillStyle(0);
  legB3dataMult->SetTextSize(0.033);
  legB3dataMult->SetBorderSize(0);
  legB3dataMult->SetTextSize(0.045);

  legB3dataMult->AddEntry(gB3vsR_PbPb276TeV[ip], "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "p");
  if (plotPbPb5TeV)  legB3dataMult->AddEntry(gB3vsR_PbPb502TeV[ip], "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, preliminary", "p");
  if (plotpPb) legB3dataMult->AddEntry(gB3vsR_pPb502TeV[ip], "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV", "p");
  legB3dataMult->AddEntry(gB3vsR_pp7TeV[ip], "pp #sqrt{#it{s}} = 7 TeV", "p");
  //  legB3dataMult->AddEntry(gB3vsR_pp13TeV_sys[ip], "pp #sqrt{#it{s}} = 13 TeV, prelim.", "pf");
  
  cb2vsdNdeta->cd();
  gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetTickx();
  gPad->SetTicky();
  
  pavept->Draw();
  if (!plotLinX) gPad->SetLogx();
  
  hframeMult->Draw();
  gB3vsR_PbPb276TeV_sys[ip]->Draw("p3");
  if (plotpPb) gB3vsR_pPb502TeV_sys[ip]->Draw("p3");

  gB3vsR_pp7TeV_sys[ip]->Draw("2");
  gB3vsR_pp7TeV[ip]->Draw("samep");

  if (plotPbPb5TeV) {
      gB3vsR_PbPb502TeV_sys[ip]->Draw("p3");
      gB3vsR_PbPb502TeV[ip]->Draw("samep");
  }

  gB3blastvsR_PbPb276TeV[0]->Draw("lsame");
  gB3vsR_PbPb276TeV[ip]->Draw("samep");

  if (plotpPb) gB3vsR_pPb502TeV[ip]->Draw("samep");
  if (plotpPb) hB3_coalescenceParam0->Draw("l");
  hB3_coalescenceParam1->Draw("l");

  //gB3blastvsR_pPb5TeV[0]->Draw("lsame");
    
  // gB3vsR_pp13TeV_sys[ip]->Draw("p3");
  // gB3vsR_pp13TeV[ip]->Draw("samep");

  pavept->Draw();
  legB3dataMult->Draw();
  legB3coal->Draw();
  legB3therm->Draw();
  
  if (plotpPb) {
    figPath = "./LDtalkFigure";
    cb2vsdNdeta->SaveAs(Form("%s/B3vsMult_pt%03.0f.pdf", figPath.Data(), pToA*100));
    cb2vsdNdeta->SaveAs(Form("%s/B3vsMult_pt%03.0f.png", figPath.Data(), pToA*100));
    TString foutName = Form("%s/B3vsMult_R%i.root", figPath.Data(), RmappingParam);
    TFile * fout = new TFile(foutName.Data(), "recreate");
    fout->cd();
    cb2vsdNdeta->Write();
    fout->Close();
  } else
  {
    figPath = "../LDtalkFigure/";
    cb2vsdNdeta->SaveAs(Form("%s/B3vsMult%03.0f.pdf", figPath.Data(),pToA*100));
    cb2vsdNdeta->SaveAs(Form("%s/B3vsMult%03.0f.eps", figPath.Data(),pToA*100));
    cb2vsdNdeta->SaveAs(Form("%s/B3vsMult%03.0f.png", figPath.Data(),pToA*100));
  }
  

  return cb2vsdNdeta;
  
}

Double_t GetB3ratioTritiumToHe(Double_t Rsource = 0.77) // fm pp 7 TeV INEL HBT param A
{

  Double_t R3He = 2.48; //fm
  Double_t Rt = 2.15; //fm
  return TMath::Power( (Rsource*Rsource + R3He*R3He /4.) / (Rsource*Rsource + Rt*Rt /4.) , 3.);

}
