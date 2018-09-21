void GetB3inpToA_PbPb5TeV();
void GetB3inpToA_PbPb276TeV();
void GetB3inpToA_pp7TeV();
void GetB3LambdainpToA_PbPb276TeV();
void GetB3inpToA_pPb5TeV(Double_t pToA = 0.733);

void GetB3inpToA_pp7TeV()
{
  //for B3 https://arxiv.org/abs/1709.08522
  TFile * fin = TFile::Open("He3-pp-7TeV-INEL-B3.root");
  if (!fin) return;

  TDirectoryFile *c1 = (TDirectoryFile *) fin->Get("pp-7TeV-INEL");
  if (!c1) return;

  TGraphAsymmErrors * gB30100 = (TGraphAsymmErrors *) c1->Get("AntiHe3_B3_Pt");
  TGraphAsymmErrors * gB30100_sys = (TGraphAsymmErrors *) c1->Get("AntiHe3_SystErr_B3_Pt");
  if (!gB30100) return;
  if (!gB30100_sys) return;
  
  //for multi INEL Eur. Phys. J. C 77 (2017) 33, Link: https://link.springer.com/article/10.1140/epjc/s10052-016-4571-1
  const Int_t nP = 3;
  Double_t dndeta[nP] = {4.60};
  Double_t dndetaErr[nP] = {0.34};

  //Following are INEL > 0
  //  Double_t dndeta[nP] = {5.98};
  // Double_t dndetaErr[nP] = {0.09};

  TFile * fout = new TFile("B3pToA_pp7TeV.root","recreate");

  for (Int_t ip = 0; ip<nP; ip++){
    
      Double_t pToA = gB30100->GetX()[ip];
      TGraphAsymmErrors * gout = new TGraphAsymmErrors(1);
      gout->SetName(Form("B3_pp7TeV_pToA=%4.3f", pToA));
      gout->SetPoint(0, dndeta[0], gB30100->GetY()[ip]);
      gout->SetPointEXhigh(0, dndetaErr[0]);
      gout->SetPointEXlow(0, dndetaErr[0]);
      gout->SetPointEYhigh(0, gB30100->GetEYhigh()[ip]);
      gout->SetPointEYlow(0, gB30100->GetEYlow()[ip]);
     
      TGraphAsymmErrors * gout_sys = new TGraphAsymmErrors(1);
      gout_sys->SetName(Form("B3_pp7TeV_pToA=%4.3f_sys", pToA));
      gout_sys->SetPoint(0, dndeta[0], gB30100->GetY()[ip]);
      gout_sys->SetPointEXhigh(0, dndetaErr[0]);
      gout_sys->SetPointEXlow(0, dndetaErr[0]);
      gout_sys->SetPointEYhigh(0, gB30100_sys->GetEYhigh()[ip]);
      gout_sys->SetPointEYlow(0, gB30100_sys->GetEYlow()[ip]);
      
      fout->cd();
      gout->Write();
      gout_sys->Write();

  }
  
  return;

}
  
void GetB3inpToA_PbPb276TeV()
{
  //for B3 https://arxiv.org/abs/1506.08951
  TFile * fin = TFile::Open("B3_PbPb276TeV.root");
  if (!fin) return;

  TCanvas *c1 = (TCanvas *) fin->Get("canvB3vsPtHe");
  if (!c1) return;

  TGraphErrors * gB32080 = (TGraphErrors *) ((TList *)c1->GetListOfPrimitives())->FindObject("Graph_from_hist1stat3He");
  TGraphErrors * gB3020 = (TGraphErrors *) ((TList *)c1->GetListOfPrimitives())->FindObject("Graph_from_hist0stat3He");
  
  TGraphErrors * gB32080_sys = (TGraphErrors *) ((TList *)c1->GetListOfPrimitives())->FindObject("Graph_from_hist1syst3He");
  TGraphErrors * gB3020_sys = (TGraphErrors *) ((TList *)c1->GetListOfPrimitives())->FindObject("Graph_from_hist0syst3He");
  
  //for multi 0-80%  http://prl.aps.org/abstract/PRL/v106/i3/e032301, Phys. Rev. Lett. 106, 032301 (2011)
  Double_t dndeta[2] = {1206.75, 266.17};
  Double_t dndetaErr[2] = {45.75, 9.83};

  TFile * fout = new TFile("B3pToA_PbPb276TeV.root","recreate");

  for (Int_t ip = 0; ip<gB3020->GetN()-1; ip++){
    
      Double_t pToA = gB3020->GetX()[ip];
      TGraphErrors * gout = new TGraphErrors(2);
      gout->SetName(Form("B3_PbPb10_pToA=%4.3f", pToA));
      gout->SetPoint(0, dndeta[0], gB3020->GetY()[ip]);
      gout->SetPointError(0, dndetaErr[0], gB3020->GetEY()[ip]);
      gout->SetPoint(1, dndeta[1], gB32080->GetY()[ip]);
      gout->SetPointError(1, dndetaErr[1], gB32080->GetEY()[ip]);

      TGraphErrors * gout_sys = new TGraphErrors(2);
      gout_sys->SetName(Form("B3_PbPb10_pToA=%4.3f_sys", pToA));
      gout_sys->SetPoint(0, dndeta[0], gB3020_sys->GetY()[ip]);
      gout_sys->SetPointError(0, dndetaErr[0], gB3020_sys->GetEY()[ip]);
      gout_sys->SetPoint(1, dndeta[1], gB32080_sys->GetY()[ip]);
      gout_sys->SetPointError(1, dndetaErr[1], gB32080_sys->GetEY()[ip]);
      
      fout->cd();
      gout->Write();
      gout_sys->Write();

      if (ip == 3){
	TCanvas *c2 = new TCanvas("c2","",700,500);
	c2->SetFillColor(kWhite);
	c2->GetFrame()->SetFillColor(kWhite);
	c2->GetFrame()->SetBorderSize(12);
	c2->cd();
	gout->Draw("APL");
	gout_sys->Draw("APLsame");
      }
  }
  
  return;

}
  

void GetB3inpToA_PbPb5TeV()
{
  TFile * fin = TFile::Open("B3_PbPb5TeV_preliminarySQM17.root");
  if (!fin) return;

  TCanvas *c1 = (TCanvas *) fin->Get("b2canv");
  if (!c1) return;

  TGraphErrors * gB34090 = (TGraphErrors *) ((TList *)c1->GetListOfPrimitives())->FindObject("b22");
  TGraphErrors * gB31040 = (TGraphErrors *) ((TList *)c1->GetListOfPrimitives())->FindObject("b21");
  TGraphErrors * gB3010 = (TGraphErrors *) ((TList *)c1->GetListOfPrimitives())->FindObject("b20");
  
  TGraphErrors * gB34090_sys = (TGraphErrors *) ((TList *)c1->GetListOfPrimitives())->FindObject("b2syst2");
  TGraphErrors * gB31040_sys = (TGraphErrors *) ((TList *)c1->GetListOfPrimitives())->FindObject("b2syst1");
  TGraphErrors * gB3010_sys = (TGraphErrors *) ((TList *)c1->GetListOfPrimitives())->FindObject("b2syst0");
  
    
  Int_t dim[3] = {gB3010->GetN(), gB31040->GetN(), gB34090->GetN()};

  //for 0-80%  http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.116.222302
  //for 80-90% https://doi.org/10.1016/j.physletb.2017.07.017
  Double_t dndeta[3] = {1763.25, 826.00, 131.94};
  Double_t dndetaErr[3] = {49.75, 22.00, 6.2};

  TFile * fout = new TFile("B3pToA.root","recreate");

  for (Int_t ip = 0; ip<gB3010->GetN(); ip++){
    
      Double_t pToA = gB3010->GetX()[ip];
      TGraphErrors * gout = new TGraphErrors(3);
      gout->SetName(Form("B3_PbPb15_pToA=%4.3f", pToA));
      gout->SetPoint(0, dndeta[0], gB3010->GetY()[ip]);
      gout->SetPointError(0, dndetaErr[0], gB3010->GetEY()[ip]);
      gout->SetPoint(1, dndeta[1], gB31040->GetY()[ip]);
      gout->SetPointError(1, dndetaErr[1], gB31040->GetEY()[ip]);

      TGraphErrors * gout_sys = new TGraphErrors(3);
      gout_sys->SetName(Form("B3_PbPb15_pToA=%4.3f_sys", pToA));
      gout_sys->SetPoint(0, dndeta[0], gB3010_sys->GetY()[ip]);
      gout_sys->SetPointError(0, dndetaErr[0], gB3010_sys->GetEY()[ip]);
      gout_sys->SetPoint(1, dndeta[1], gB31040_sys->GetY()[ip]);
      gout_sys->SetPointError(1, dndetaErr[1], gB31040_sys->GetEY()[ip]);

      if (ip<6) {
	gout->SetPoint(2, dndeta[2], gB34090->GetY()[ip]);
	gout->SetPointError(2, dndetaErr[2], gB34090->GetEY()[ip]);
	gout_sys->SetPoint(2, dndeta[2], gB34090_sys->GetY()[ip]);
	gout_sys->SetPointError(2, dndetaErr[2], gB34090_sys->GetEY()[ip]);
      } else {
	gout->RemovePoint(2);
	gout_sys->RemovePoint(2);	
      }
      fout->cd();
      gout->Write();
      gout_sys->Write();
  }
  
  return;

}


void GetB3inpToA_pPb5TeV(Double_t pToA = 0.733)
{
  TFile * fin = TFile::Open("B3_pPb_5TeV_21092018.root");
  if (!fin) return;

  Int_t cent[5] = {0, 10, 20, 40, 100}; 
  TH1D * hin[4];
  TH1D * hin_syst[4];
  //                       0-10    10-20     20-40    40-100
  Double_t dndeta[4] =    {40.6,    30.5,     23.2,     10.1}; 
  Double_t dndetaErr[4] = {0.9,     0.7,     0.5,        0.16};
  
  for (int i = 0; i < 4; i++){
    hin[i] = (TH1D *) fin->Get(Form("B3_3He_%i_%i", cent[i], cent[i+1]));
    hin_syst[i] = (TH1D *) fin->Get(Form("B3_3He_%i_%i_Syst", cent[i], cent[i+1]));
  }

  //selected pt bin  
  TFile * fout = new TFile("B3pToA_pPb502TeV.root","recreate");
  TGraphErrors * gout = new TGraphErrors(4);
  gout->SetName(Form("B3_pPb15_pToA=%4.3f", pToA));
  
  TGraphErrors * gout_sys = new TGraphErrors(4);
  gout_sys->SetName(Form("B3_pPb15_pToA=%4.3f_sys", pToA));
  
  for (Int_t ip = 0; ip<4; ip++){
    Int_t sb = hin[ip]->GetXaxis()->FindBin(pToA);
    
    gout->SetPoint(ip, dndeta[ip], hin[ip]->GetBinContent(sb));
    gout->SetPointError(ip, dndetaErr[ip], hin[ip]->GetBinError(sb));
    
    gout_sys->SetPoint(ip, dndeta[ip], hin_syst[ip]->GetBinContent(sb));
    gout_sys->SetPointError(ip, dndetaErr[ip], hin_syst[ip]->GetBinError(sb));
  } 
  fout->cd();
  gout->Write();
  gout_sys->Write();
  return;
  
}
  
  
void GetB3LambdainpToA_PbPb276TeV()
{

  //from ALICE paper: https://hepdata.net/record/ins1380234
  Double_t xpToA[3] = {1.0, 1.67, 2.67};

  //the values from the paper have to be transformed from B2 to B3 as explained in fig. 7 of
  //https://www.sciencedirect.com/science/article/pii/S0370269316000575?via%3Dihub#ec-research-data
  Double_t B2L[3] = {0.000444, 0.000588, 0.00295};
  Double_t B2Lstat[3] = {0.000111, 0.000202, 0.00134};
  Double_t B2Lsys[3] = {0.0000531, 0.0000117, 0.000553};
  
  Double_t B3L[3];
  Double_t B3Lstat[3];
  Double_t B3Lsys[3];
  
  const Int_t nP = 3;

  const Double_t k2 = TMath::Power(1.875612928 / 0.93827231, 2.0) * 1.115683 / 2.991;
  //transform the B2 into B3L
  for (Int_t ip = 0; ip<nP; ip++){
    B3L[ip] = TMath::Power(B2L[ip], 2.0) / k2;
    B3Lstat[ip] = B3L[ip] * (2.0 * B2Lstat[ip] / B2L[ip] / k2);
    B3Lsys[ip] = B3L[ip] * (2.0 * B2Lsys[ip] / B2L[ip] / k2);
  }
  
  TGraphErrors * gB3L010 =  new TGraphErrors(3, xpToA, B3L, 0, B3Lstat);
  TGraphErrors * gB3L010_sys =  new TGraphErrors(3, xpToA, B3L, 0, B3Lsys);
  
  if (!gB3L010) return;
  if (!gB3L010_sys) return;
  
  Double_t dndeta[1] = {1447.5};
  Double_t dndetaErr[1] = {54.5};

  TFile * fout = new TFile("B3LambdapToA_PbPb276TeV.root","recreate");

  for (Int_t ip = 0; ip<nP; ip++){
    
      Double_t pToA = gB3L010->GetX()[ip];
      TGraphAsymmErrors * gout = new TGraphAsymmErrors(1);
      gout->SetName(Form("B3Lambda_PbPb276TeV_pToA=%4.3f", pToA));
      gout->SetPoint(0, dndeta[0], gB3L010->GetY()[ip]);
      gout->SetPointEXhigh(0, dndetaErr[0]);
      gout->SetPointEXlow(0, dndetaErr[0]);
      gout->SetPointEYhigh(0, gB3L010->GetEY()[ip]);
      gout->SetPointEYlow(0, gB3L010->GetEY()[ip]);
     
      TGraphAsymmErrors * gout_sys = new TGraphAsymmErrors(1);
      gout_sys->SetName(Form("B3Lambda_PbPb276TeV_pToA=%4.3f_sys", pToA));
      gout_sys->SetPoint(0, dndeta[0], gB3L010_sys->GetY()[ip]);
      gout_sys->SetPointEXhigh(0, dndetaErr[0]);
      gout_sys->SetPointEXlow(0, dndetaErr[0]);
      gout_sys->SetPointEYhigh(0, gB3L010_sys->GetEY()[ip]);
      gout_sys->SetPointEYlow(0, gB3L010_sys->GetEY()[ip]);
      
      fout->cd();
      gout->Write();
      gout_sys->Write();

  }
  
  return;


  
}



