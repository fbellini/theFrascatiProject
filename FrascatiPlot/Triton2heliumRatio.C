#include "./Make3HepPbPaperFigure.C"
Float_t getRho(Float_t radius, Float_t errt, Float_t err3He);
TGraphErrors * MakeRhoGraph();

void Triton2heliumRatio()
{
  //set style
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadBottomMargin(0.2);
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadRightMargin(0.02);
   
  TH2D * hframeMult = new TH2D("hframeMult", "#rho vs mult; #LTd#it{N}_{ch}/d#it{#eta_{lab}}#GT; #rho", 3000, 0., 3000.0, 180, 0.0, 2.2);
  hframeMult->GetXaxis()->SetTitleSize(0.06);
  hframeMult->GetYaxis()->SetTitleSize(0.06);
  hframeMult->GetXaxis()->SetTitleOffset(1.2);
  hframeMult->GetXaxis()->SetLabelSize(0.05);
  hframeMult->GetYaxis()->SetLabelSize(0.05);
  hframeMult->GetYaxis()->SetNdivisions(505);
  hframeMult->GetXaxis()->SetRangeUser(1., 3.E3);

  //Define pT/A labels only once
  TPaveText * pavept = new TPaveText(0.67, 0.85, 0.87, 0.9, "NDC");
  pavept->SetFillStyle(0);
  pavept->SetTextFont(42);
  pavept->SetBorderSize(0);
  pavept->SetTextSize(0.04);
  pavept->SetTextAlign(12);
  pavept->AddText(Form("#it{p}_{T}/#it{A} = 0.75 GeV/#it{c}"));

TGraphErrors * gRhoVsR = (TGraphErrors*) MakeRhoGraph();
gRhoVsR->SetLineStyle(1);
gRhoVsR->SetLineColor(kBlack);
gRhoVsR->SetLineWidth(4);
convertRadiusToMulti(gRhoVsR, 1);

Float_t arrowsize = 0.02;
TArrow * arrpp = new TArrow(2.0, 0.8, 22.0, 0.8, arrowsize, "<|>");
arrpp->SetLineColor(kGreen+2);
arrpp->SetFillColor(kGreen+2);
arrpp->SetFillStyle(1001);
arrpp->SetLineWidth(3);

TPaveText * pavepp = new TPaveText(2.0, 0.70, 20.0, 0.75);
  pavepp->SetFillStyle(0);
  pavepp->SetTextFont(42);
  pavepp->SetBorderSize(0);
  pavepp->SetTextSize(0.04);
  pavepp->SetTextAlign(22);
  pavepp->SetTextColor(kGreen+2);
  pavepp->AddText("ALICE pp 7 TeV");

TArrow * arrpPb = new TArrow(4.0, 0.6, 46.0, 0.6, arrowsize, "<|>");
arrpPb->SetLineColor(kBlue);
arrpPb->SetFillColor(kBlue);
arrpPb->SetFillStyle(1001);
arrpPb->SetLineWidth(3);

TPaveText * paveppb = new TPaveText(4.0, 0.48, 46., 0.53);
  paveppb->SetFillStyle(0);
  paveppb->SetTextFont(42);
  paveppb->SetBorderSize(0);
  paveppb->SetTextSize(0.04);
  paveppb->SetTextAlign(22);
  paveppb->SetTextColor(kBlue);
  paveppb->AddText("ALICE p-Pb 5.02 TeV");

TArrow * arrPbPb = new TArrow(17.0, 0.4, 2050.0, 0.4, arrowsize);
arrPbPb->SetOption("<>");
arrPbPb->SetLineColor(kRed+1);
arrPbPb->SetFillColor(kRed+1);
arrPbPb->SetFillStyle(1001);
arrPbPb->SetLineWidth(3);

TPaveText * pavepbpb = new TPaveText(17.0, 0.30, 2050.0, 0.35);
  pavepbpb->SetFillStyle(0);
  pavepbpb->SetTextFont(42);
  pavepbpb->SetBorderSize(0);
  pavepbpb->SetTextSize(0.04);
  pavepbpb->SetTextAlign(22);
  pavepbpb->SetTextColor(kRed+1);
  pavepbpb->AddText("ALICE Pb-Pb 5.02 TeV");

TLine * thermal = new TLine(1.0, 1.0, 3000., 1.0);
thermal->SetLineColor(kOrange+1);
thermal->SetLineWidth(3);
thermal->SetLineStyle(2);

TLegend * leg = new TLegend(0.4, 0.8, 0.9, 0.9);
leg->SetFillStyle(0);
leg->SetTextSize(0.04);
leg->SetBorderSize(0);
leg->AddEntry(gRhoVsR, "Coalescence, #it{p}_{T}/#it{A} = 0.75 GeV/#it{c}","l");
leg->AddEntry(thermal, "Statistical hadronisation model", "l");

TCanvas * c1 = new TCanvas("c1","c1", 800, 600);
c1->SetTickx();
c1->SetTicky();
c1->cd();
gPad->SetLogx();
hframeMult->Draw();
gRhoVsR->Draw("same l");
arrpp->Draw("same");
arrpPb->Draw("same");
arrPbPb->Draw("same");
pavepp->Draw("same");
paveppb->Draw("same");
pavepbpb->Draw("same");
//pavept->Draw();
leg->Draw();
thermal->Draw();
return;

}

Float_t getRho(Float_t radius, Float_t errt, Float_t err3He)
{
 const Float_t rt = 2.15;
 const Float_t r3He = 2.48;
 return TMath::Power((radius*radius + TMath::Power((r3He+err3He)/2., 2.0))/(radius*radius + TMath::Power((rt+errt)/2., 2.0)), 3.0);
}

TGraphErrors * MakeRhoGraph()
{
  const Int_t nPoints = 1000;
  Double_t gY[nPoints];
  Double_t gR[nPoints];
  
  TGraphErrors * graphOut = new TGraphErrors(nPoints);
  graphOut->SetName("rho_th_coalescence");
  graphOut->SetTitle("#rho from coalescence");
  
  for (int i = 0; i<nPoints; i++){
    gR[i] = 10.0 * i / nPoints; 
    gY[i] = getRho(gR[i], 0.0, 0.0);

    graphOut->SetPoint(i, gR[i], gY[i]);
    graphOut->SetPointError(i, 0.0, 0.0);
  }

  return graphOut;

}
