
void Make13TeVPaperFigure(Bool_t plotLinX, Double_t pToA, Int_t RmappingParam,
			  TGraphErrors * hB2_coalescence, 
			  TGraphErrors ** gB2vsR_PbPb276TeV_sys,  TGraphErrors ** gB2vsR_pp7TeV_sys, TGraphErrors **gB2vsR_PbPb276TeV, TGraphErrors **  gB2vsR_pp7TeV,
			  TGraphErrors ** gB2vsR_pPb502TeV_sys,  TGraphErrors ** gB2vsR_pp13TeV_sys, TGraphErrors ** gB2vsR_pPb502TeV, TGraphErrors **  gB2vsR_pp13TeV)
{

  Int_t ip = RmappingParam;
  
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

  TString foutName = Form("B2vsR_w13TeV_param%i.root", RmappingParam);
  TFile * fout = new TFile(foutName.Data(), "recreate");
  fout->cd();
  cb2opta->Write();
  fout->Close();
  return;
  
}

