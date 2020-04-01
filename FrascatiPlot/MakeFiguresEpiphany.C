#include "./CoalescenceBA.C"
//#include "./MapMulti2R.C"
//#include "./Make3HepPbPaperFigure.Cs"
void MakeFiguresEpiphany()
{
    TString figPath = "../EPIPHANY2019proceedings";

    //Figure 1: BA from coalescence
    CoalescenceBA(0.75, figPath.Data());
   //Make3HepPbPaperFigure(figPath.Data(), 0.733, 0, 1, 1, 0, 0);

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTopMargin(0.03);
    gStyle->SetPadBottomMargin(0.17);
    gStyle->SetPadLeftMargin(0.17);
    gStyle->SetPadRightMargin(0.02); 

    Bool_t convertToR = kFALSE;    //Bool_t plotPseudoData = 0;
    Int_t paramSet = 1;

    Color_t projCol = kRed+1;
    Double_t Rmin = 0.8; 
    Double_t Rmax = 3000;
    if (convertToR) {
        Rmin = 0.01;
        Rmax = 6.5;//fm
    }

    Double_t pToA = 0.75;
    Double_t pToAb3 = 0.733;
    Double_t pToAb3pp = 0.8;
    Double_t pToAb3Lambda = 1.;
    Double_t pToAb4 = 0.75; 
    Double_t pToAb4Lambda = 0.75; 

    //-----------------------------
    //    data
    //-----------------------------
    Int_t Fill_Style = 1001;
    Int_t Line_Style = 1;
    Int_t Line_Style_Blast = 2;
    Int_t Line_Width = 1; 
    Int_t Line_Width_Blast = 3;
    Float_t Marker_Size = 1.3;
    //nb: manual hack to avoid loop over all parameterisations
    const Int_t nParamSet = 1;
    Int_t ip = paramSet-1;

    TGraphErrors* gB2vsR_pp7TeVINELg0[nParamSet];
    TGraphErrors* gB2vsR_pp7TeVINELg0_sys[nParamSet];
    gB2vsR_pp7TeVINELg0[ip] = (TGraphErrors *) getB2_pp7TeVINELg0(kFALSE, pToA, paramSet);
    gB2vsR_pp7TeVINELg0_sys[ip] = (TGraphErrors *) getB2_pp7TeVINELg0(kTRUE, pToA, paramSet);
    MakeUp(gB2vsR_pp7TeVINELg0_sys[ip], kGreen+2, kGreen+1, Fill_Style, Line_Style, Line_Width, 21, Marker_Size);
    MakeUp(gB2vsR_pp7TeVINELg0[ip]    , kGreen+2, kGreen+1, Fill_Style, Line_Style, Line_Width, 21, Marker_Size);

    TGraphErrors* gB2vsR_PbPb276TeV[nParamSet];
    TGraphErrors* gB2vsR_PbPb276TeV_sys[nParamSet];
    gB2vsR_PbPb276TeV[ip] = (TGraphErrors *) getB2_PbPb276TeV(kFALSE, pToA, paramSet);
    gB2vsR_PbPb276TeV_sys[ip] = (TGraphErrors *) getB2_PbPb276TeV(kTRUE, pToA, paramSet);
    MakeUp(gB2vsR_PbPb276TeV_sys[ip], kRed+1, kRed+1, Fill_Style, Line_Style, Line_Width, 20, Marker_Size);
    MakeUp(gB2vsR_PbPb276TeV[ip]    , kRed+1, kRed+1, Fill_Style, Line_Style, Line_Width, 20, Marker_Size);

    TGraphAsymmErrors* gB3vsR_pp7TeV[nParamSet];
    TGraphAsymmErrors* gB3vsR_pp7TeV_sys[nParamSet];
    gB3vsR_pp7TeV[ip] = (TGraphAsymmErrors *) getB3_pp7TeVINELg0(kFALSE, pToAb3pp, paramSet);
    gB3vsR_pp7TeV_sys[ip] = (TGraphAsymmErrors *) getB3_pp7TeVINELg0(kTRUE, pToAb3pp, paramSet);
    MakeUp(gB3vsR_pp7TeV_sys[ip], kGreen+2, kGreen+1, Fill_Style, Line_Style, Line_Width, 21, Marker_Size);
    MakeUp(gB3vsR_pp7TeV[ip]    , kGreen+2, kGreen+1, Fill_Style, Line_Style, Line_Width, 21, Marker_Size);

    TGraphErrors* gB3vsR_PbPb276TeV[nParamSet];
    TGraphErrors* gB3vsR_PbPb276TeV_sys[nParamSet];
    gB3vsR_PbPb276TeV[ip] = (TGraphErrors *) getB3_PbPb276TeV(kFALSE, pToAb3, paramSet);
    gB3vsR_PbPb276TeV_sys[ip] = (TGraphErrors *) getB3_PbPb276TeV(kTRUE, pToAb3, paramSet);
    MakeUp(gB3vsR_PbPb276TeV_sys[ip], kRed+1, kRed+1, Fill_Style, Line_Style, Line_Width, 20, Marker_Size);
    MakeUp(gB3vsR_PbPb276TeV[ip]    , kRed+1, kRed+1, Fill_Style, Line_Style, Line_Width, 20, Marker_Size);


    TGraphAsymmErrors* gB3LambdavsR_PbPb276TeV[nParamSet];
    TGraphAsymmErrors* gB3LambdavsR_PbPb276TeV_sys[nParamSet];
    gB3LambdavsR_PbPb276TeV[ip] = (TGraphAsymmErrors *) getB3Lambda_PbPb276TeV(kFALSE, pToAb3Lambda, paramSet);
    gB3LambdavsR_PbPb276TeV_sys[ip] = (TGraphAsymmErrors *) getB3Lambda_PbPb276TeV(kTRUE, pToAb3Lambda, paramSet);
    MakeUp(gB3LambdavsR_PbPb276TeV_sys[ip], kRed+1, kRed+1, Fill_Style, Line_Style, Line_Width, 20, Marker_Size);
    MakeUp(gB3LambdavsR_PbPb276TeV[ip]    , kRed+1, kRed+1, Fill_Style, Line_Style, Line_Width, 20, Marker_Size);

    if (!convertToR) {
        convertRadiusToMulti(gB2vsR_pp7TeVINELg0[ip], paramSet); 
        convertRadiusToMulti(gB2vsR_pp7TeVINELg0_sys[ip], paramSet); 
        convertRadiusToMulti(gB3vsR_pp7TeV[ip], paramSet); 
        convertRadiusToMulti(gB3vsR_pp7TeV_sys[ip], paramSet); 
        convertRadiusToMulti(gB2vsR_PbPb276TeV[ip], paramSet); 
        convertRadiusToMulti(gB2vsR_PbPb276TeV_sys[ip], paramSet); 
        convertRadiusToMulti(gB3vsR_PbPb276TeV[ip], paramSet); 
        convertRadiusToMulti(gB3vsR_PbPb276TeV_sys[ip], paramSet); 
        convertRadiusToMulti(gB3LambdavsR_PbPb276TeV[ip], paramSet); 
        convertRadiusToMulti(gB3LambdavsR_PbPb276TeV_sys[ip], paramSet); 
    }
    //-----------------------------
    //theory - Blast Wave + thermal
    //-----------------------------
    // pToAb3 = 0.77;
    // pToAb3Lambda = 1.17;
    // pToAb4 = 0.75;
    // pToAb4Lambda = 0.62;
    TGraphAsymmErrors* gBlastB2vsR_PbPb276TeV = (TGraphAsymmErrors *)  getBAthermalBlast("PbPb276TeV", "deuteron", pToA, kBlue, paramSet, convertToR);
    TGraphAsymmErrors* gBlastB3vsR_PbPb276TeV = (TGraphAsymmErrors *)  getBAthermalBlast("PbPb276TeV", "He3", pToAb3, kBlue, paramSet, convertToR);
    TGraphAsymmErrors* gBlastB3LambdavsR_PbPb276TeV = (TGraphAsymmErrors *) getBAthermalBlast("PbPb276TeV", "hyper-triton", pToAb3Lambda, kBlue, paramSet, convertToR);
    TGraphAsymmErrors* gBlastB4vsR_PbPb276TeV = (TGraphAsymmErrors *) getBAthermalBlast("PbPb276TeV", "He4", pToAb4, kBlue, paramSet, convertToR);
    TGraphAsymmErrors* gBlastB4LambdavsR_PbPb276TeV = (TGraphAsymmErrors *) getBAthermalBlast("PbPb276TeV", "4LH", pToAb4Lambda, kBlue, paramSet, convertToR);

    //--------------------
    //theory - coalescence
    //--------------------
    //mT is the mass of the particle relative to which the HBT radius is calculated.
    //We use now the mass of the proton and the pT per nucleon
    Double_t mT2 = TMath::Sqrt(pToA * pToA + 0.938 * 0.938);
    Double_t mT3 = TMath::Sqrt(pToAb3 * pToAb3 + 0.938 * 0.938);
    Double_t mT3L = TMath::Sqrt(pToAb3Lambda * pToAb3Lambda + 0.938 * 0.938);
    Double_t mT4 = TMath::Sqrt(pToAb4 * pToAb4 + 0.938 * 0.938);
    Double_t mT4L = TMath::Sqrt(pToAb4Lambda * pToAb4Lambda + 0.938 * 0.938);
        
    //d
    TGraphErrors* hB2_coalescence = (TGraphErrors*) MakeB2TheoryGraphCoalescence(mT2, 3.2);
    hB2_coalescence->SetMarkerStyle(20);
    hB2_coalescence->SetMarkerColor(kBlack);
    hB2_coalescence->SetLineColor(kBlack);
    hB2_coalescence->SetLineWidth(3);
    //tritium
    TGraphErrors* hB3t_coalescence = (TGraphErrors*) MakeB3TheoryGraphCoalescence(mT3, 2.15);
    hB3t_coalescence->SetMarkerStyle(20);
    hB3t_coalescence->SetMarkerColor(kBlack);
    hB3t_coalescence->SetLineColor(kBlack);
    hB3t_coalescence->SetLineWidth(3);
    //3He
    TGraphErrors* hB3_coalescence = (TGraphErrors*) MakeB3TheoryGraphCoalescence(mT3, 2.48);
    hB3_coalescence->SetMarkerStyle(20);
    hB3_coalescence->SetMarkerColor(kBlack);
    hB3_coalescence->SetLineColor(kBlack);
    hB3_coalescence->SetLineWidth(3);
    //hypertriton
    TGraphErrors* hB3L_coalescence = (TGraphErrors*) MakeB3TheoryGraphCoalescence(mT3L, 6.8);
    hB3L_coalescence->SetMarkerStyle(20);
    hB3L_coalescence->SetMarkerColor(kBlack);
    hB3L_coalescence->SetLineColor(kBlack);
    hB3L_coalescence->SetLineWidth(3);
    hB3L_coalescence->SetLineStyle(1);
    //hypertriton larger r
    TGraphErrors* hB3L_coalescence_largeradius = (TGraphErrors*) MakeB3TheoryGraphCoalescence(mT3L, 14.1);
    hB3L_coalescence_largeradius->SetMarkerStyle(20);
    hB3L_coalescence_largeradius->SetMarkerColor(kBlack);
    hB3L_coalescence_largeradius->SetLineColor(kBlack);
    hB3L_coalescence_largeradius->SetLineWidth(2);
    hB3L_coalescence_largeradius->SetLineStyle(5);
    //4He
    TGraphErrors* hB4_coalescence = (TGraphErrors*) MakeB4TheoryGraphCoalescence(mT4, 1.9); //He4
    hB4_coalescence->SetMarkerStyle(1);
    hB4_coalescence->SetMarkerColor(kBlack);
    hB4_coalescence->SetLineColor(kBlack);
    hB4_coalescence->SetLineWidth(3);
    //4LH
    TGraphErrors* hB4L_coalescence = (TGraphErrors*) MakeB4TheoryGraphCoalescence(mT4L, 2.4); //4LH
    hB4L_coalescence->SetMarkerStyle(20);
    hB4L_coalescence->SetMarkerColor(kBlack);
    hB4L_coalescence->SetLineColor(kBlack);
    hB4L_coalescence->SetLineWidth(3);
    hB4L_coalescence->SetLineStyle(1);
    //4LH larger r
    TGraphErrors* hB4L_coalescence_largeradius = (TGraphErrors*) MakeB4TheoryGraphCoalescence(mT4L, 4.9); //4LH
    hB4L_coalescence_largeradius->SetMarkerStyle(20);
    hB4L_coalescence_largeradius->SetMarkerColor(kBlack);
    hB4L_coalescence_largeradius->SetLineColor(kBlack);
    hB4L_coalescence_largeradius->SetLineWidth(2);
    hB4L_coalescence_largeradius->SetLineStyle(5);

    
    if (!convertToR){
        convertRadiusToMulti(hB2_coalescence, paramSet); 
        convertRadiusToMulti(hB3t_coalescence, paramSet); 
        convertRadiusToMulti(hB3_coalescence, paramSet); 
        convertRadiusToMulti(hB4_coalescence, paramSet); 
        convertRadiusToMulti(hB3L_coalescence, paramSet); 
        convertRadiusToMulti(hB3L_coalescence_largeradius, paramSet); 
        convertRadiusToMulti(hB4L_coalescence, paramSet); 
        convertRadiusToMulti(hB4L_coalescence_largeradius, paramSet); 
    }

    //---------------------------------------
    // PLOT Yellow Report figure
    //---------------------------------------  
    
    TH2D * hframe = new TH2D("hframeFig4", "B_{2} vs radius; #it{R} (fm); #it{B}_{2} (GeV^{2}/#it{c}^{3})", 1000, Rmin, Rmax, 2000, 5.e-5, 0.05);
    hframe->GetXaxis()->SetTitleSize(0.06);
    hframe->GetYaxis()->SetTitleSize(0.06);
    hframe->GetYaxis()->SetTitleOffset(1.3);
    hframe->GetXaxis()->SetTitleOffset(1.1);
    hframe->GetXaxis()->SetLabelSize(0.05);
    hframe->GetYaxis()->SetLabelSize(0.05);
    hframe->GetXaxis()->SetRangeUser(Rmin, Rmax);
    
    TH2D * hframe3t = new TH2D("hframe3Fig4", "B_{3} vs radius; #it{R} (fm); #it{B}_{3} (GeV^{4}/#it{c}^{6})", 1000, Rmin, Rmax, 2000, 5.e-9, 0.05);
    hframe3t->GetXaxis()->SetTitleSize(0.06);
    hframe3t->GetYaxis()->SetTitleSize(0.06);
    hframe3t->GetYaxis()->SetTitleOffset(1.3);
    hframe3t->GetXaxis()->SetTitleOffset(1.1);
    hframe3t->GetXaxis()->SetLabelSize(0.05);
    hframe3t->GetYaxis()->SetLabelSize(0.05);
    hframe3t->GetXaxis()->SetRangeUser(Rmin, Rmax);

    TH2D * hframe3 = new TH2D("hframe3Fig4", "B_{3} vs radius; #it{R} (fm); #it{B}_{3} (GeV^{4}/#it{c}^{6})", 1000, Rmin, Rmax, 2000, 2.e-11, 5.e-3);
    hframe3->GetXaxis()->SetTitleSize(0.06);
    hframe3->GetYaxis()->SetTitleSize(0.06);
    hframe3->GetYaxis()->SetTitleOffset(1.3);
    hframe3->GetXaxis()->SetTitleOffset(1.1);
    hframe3->GetXaxis()->SetLabelSize(0.05);
    hframe3->GetYaxis()->SetLabelSize(0.05);
    hframe3->GetXaxis()->SetRangeUser(Rmin, Rmax);

    TH2D * hframe3L = new TH2D("hframe3LFig4", "B_{3,#Lambda} vs radius; #it{R} (fm); #it{B}_{3,#Lambda} (GeV^{4}/#it{c}^{6})", 1000, Rmin, Rmax, 2000, 2.e-11, 5.e-3);
    hframe3L->GetXaxis()->SetTitleSize(0.06);
    hframe3L->GetYaxis()->SetTitleSize(0.06);
    hframe3L->GetYaxis()->SetTitleOffset(1.3);
    hframe3L->GetXaxis()->SetTitleOffset(1.1);
    hframe3L->GetXaxis()->SetLabelSize(0.05);
    hframe3L->GetYaxis()->SetLabelSize(0.05);
    hframe3L->GetXaxis()->SetRangeUser(Rmin, Rmax);

    TH2D * hframe4 = new TH2D("hframe4Fig4", "B_{4} vs radius; #it{R} (fm); #it{B}_{4} (GeV^{6}/#it{c}^{9})", 1000, Rmin, Rmax, 2000, 2.e-13, 5.e-3);
    hframe4->GetXaxis()->SetTitleSize(0.06);
    hframe4->GetYaxis()->SetTitleSize(0.06);
    hframe4->GetYaxis()->SetTitleOffset(1.3);
    hframe4->GetXaxis()->SetTitleOffset(1.1);
    hframe4->GetXaxis()->SetLabelSize(0.05);
    hframe4->GetYaxis()->SetLabelSize(0.05);
    hframe4->GetXaxis()->SetRangeUser(Rmin, Rmax);

    TH2D * hframe4L = new TH2D("hframe4LFig4", "B_{3,#Lambda} vs radius; #it{R} (fm); #it{B}_{4,#Lambda} (GeV^{6}/#it{c}^{9})", 1000, Rmin, Rmax, 2000, 2.e-13, 5.e-3);
    hframe4L->GetXaxis()->SetTitleSize(0.06);
    hframe4L->GetYaxis()->SetTitleSize(0.06);
    hframe4L->GetYaxis()->SetTitleOffset(1.3);
    hframe4L->GetXaxis()->SetTitleOffset(1.1);
    hframe4L->GetXaxis()->SetLabelSize(0.05);
    hframe4L->GetYaxis()->SetLabelSize(0.05);
    hframe4L->GetXaxis()->SetRangeUser(Rmin, Rmax);
    
    if (!convertToR) {
        hframe->GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta_{lab}}#GT");
        hframe3->GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta_{lab}}#GT");
        hframe3t->GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta_{lab}}#GT");
        hframe3L->GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta_{lab}}#GT");
        hframe4->GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta_{lab}}#GT");
        hframe4L->GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta_{lab}}#GT");
    }

    //define particle label
    TPaveText * paveLab2 = new TPaveText(0.8, 0.8, 0.9, 0.9, "NDC");
    paveLab2->SetFillStyle(0);
    paveLab2->SetTextFont(42);
    paveLab2->SetBorderSize(0);
    paveLab2->SetTextSize(0.1);
    paveLab2->SetTextAlign(12);
    paveLab2->AddText("#bf{d}");

    TPaveText * paveLab3t = new TPaveText(0.77, 0.8, 0.9, 0.9, "NDC");
    paveLab3t->SetFillStyle(0);
    paveLab3t->SetBorderSize(0);
    paveLab3t->SetTextAlign(12);
    paveLab3t->AddText("#bf{^{3}H}");
    paveLab3t->SetTextSize(0.09);
    paveLab3t->SetTextFont(42);
    
    TPaveText * paveLab3 = new TPaveText(0.77, 0.8, 0.9, 0.9, "NDC");
    paveLab3->SetFillStyle(0);
    paveLab3->SetBorderSize(0);
    paveLab3->SetTextAlign(12);
    paveLab3->AddText("#bf{^{3}He}");
    paveLab3->SetTextSize(0.09);
    paveLab3->SetTextFont(42);

    TPaveText * paveLab3L = new TPaveText(0.75, 0.8, 0.9, 0.9, "NDC");
    paveLab3L->SetFillStyle(0);
    paveLab3L->SetTextFont(42);
    paveLab3L->SetBorderSize(0);
    paveLab3L->SetTextSize(0.1);
    paveLab3L->SetTextAlign(12);
    paveLab3L->AddText("#bf{{}^{3}_{#Lambda}H}");

    TPaveText * paveLab4 = new TPaveText(0.77, 0.8, 0.9, 0.9, "NDC");
    paveLab4->SetFillStyle(0);
    paveLab4->SetTextFont(42);
    paveLab4->SetBorderSize(0);
    paveLab4->SetTextSize(0.1);
    paveLab4->SetTextAlign(12);
    paveLab4->AddText("#bf{^{4}He}");

    TPaveText * paveLab4L = new TPaveText(0.77, 0.8, 0.9, 0.9, "NDC");
    paveLab4L->SetFillStyle(0);
    paveLab4L->SetTextFont(42);
    paveLab4L->SetBorderSize(0);
    paveLab4L->SetTextSize(0.1);
    paveLab4L->SetTextAlign(12);
    paveLab4L->AddText("#bf{{}^{4}_{#Lambda}H}");
    
    //Define pT/A labels only once
    TPaveText * pavept = new TPaveText(0.6, 0.7, 0.9, 0.78, "NDC");
    pavept->SetFillStyle(0);
    pavept->SetTextFont(42);
    pavept->SetBorderSize(0);
    pavept->SetTextSize(0.05);
    pavept->SetTextAlign(12);
    pavept->AddText(Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb4));

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


    //legends
    Float_t yl = 0.2; 
    Float_t yu = 0.4; 
    Float_t xl = 0.2; 
    Float_t xu = 0.7; 
    Float_t legTextSize = 0.045;
    TLegend * masterLeg = new TLegend(xl, yl, xu, yu, Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToA));
    masterLeg->SetFillStyle(0);
    masterLeg->SetTextSize(legTextSize);
    masterLeg->SetBorderSize(0);
    masterLeg->AddEntry(gBlastB2vsR_PbPb276TeV, "BW + GSI-Heid. (#it{T}_{chem} = 156 MeV)", "l");
    masterLeg->AddEntry(hB3_coalescence, "Coal., #it{r}(d) = 3.2 fm", "l");

    TLegend * masterLeg3t = new TLegend(xl, yl, xu, yu, Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb3));
    masterLeg3t->SetFillStyle(0);
    masterLeg3t->SetTextSize(legTextSize);
    masterLeg3t->SetBorderSize(0);
    //masterLeg3->AddEntry(gBlastB3vsR_PbPb502TeV, "BW + GSI-Heid. (#it{T}_{chem} = 156 MeV)", "l");
    masterLeg3t->AddEntry(hB3t_coalescence, "Coal., #it{r}(^{3}H) = 2.15 fm", "l");
    
    TLegend * masterLeg3 = new TLegend(xl, yl, xu, yu, Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb3));
    masterLeg3->SetFillStyle(0);
    masterLeg3->SetTextSize(legTextSize);
    masterLeg3->SetBorderSize(0);
    //masterLeg3->AddEntry(gBlastB3vsR_PbPb502TeV, "BW + GSI-Heid. (#it{T}_{chem} = 156 MeV)", "l");
    masterLeg3->AddEntry(hB3_coalescence, "Coal., #it{r}(^{3}He) = 2.48 fm", "l");
    
    TLegend * masterLeg3L = new TLegend(xl, yl, xu, yu, Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb3Lambda));
    masterLeg3L->SetFillStyle(0);
    masterLeg3L->SetTextSize(legTextSize);
    masterLeg3L->SetBorderSize(0);
    masterLeg3L->AddEntry(hB3L_coalescence, "Coal., #it{r} (^{3}_{#Lambda}H) = 6.8 fm", "l");
    masterLeg3L->AddEntry(hB3L_coalescence_largeradius, "Coal., #it{r} (^{3}_{#Lambda}H) = 14.1 fm", "l");
    
    TLegend * masterLeg4 = new TLegend(xl, yl, xu, yu, Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb4));
    masterLeg4->SetFillStyle(0);
    masterLeg4->SetTextSize(legTextSize);
    masterLeg4->SetBorderSize(0);
    masterLeg4->AddEntry(hB4_coalescence, "Coal., #it{r} (^{4}He) = 1.9 fm", "l");

    TLegend * masterLeg4L = new TLegend(xl, yl, xu, yu, Form("#it{p}_{T}/#it{A} = %3.2f GeV/#it{c}", pToAb4Lambda));
    masterLeg4L->SetFillStyle(0);
    masterLeg4L->SetTextSize(legTextSize);
    masterLeg4L->SetBorderSize(0);
    masterLeg4L->AddEntry(hB3L_coalescence, "Coal., #it{r} (^{4}_{#Lambda}H) = 2.4 fm", "l");
    masterLeg4L->AddEntry(hB4L_coalescence_largeradius, "Coal., #it{r} (^{4}_{#Lambda}H) = 4.9 fm", "l");

    TLegend * legB3data = new TLegend(xl, yu, xu, yu+0.15, "ALICE");
    legB3data->SetFillStyle(0);
    legB3data->SetTextSize(legTextSize);
    legB3data->SetBorderSize(0);
    legB3data->AddEntry(gB3vsR_PbPb276TeV_sys[ip], "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "p");
    legB3data->AddEntry(gB3vsR_pp7TeV_sys[ip], "pp #sqrt{#it{s}} = 7 TeV (INEL>0)", "p");

    //-------------------------------
    //-------------------------------
    //-------------------------------
    //   DRAW FIRST OPTION -- 6 panels
    //-------------------------------
    //-------------------------------
    //-------------------------------

    TCanvas * cr4 = new TCanvas("cr4", "compare thermal with coalescence", 1000, 2000);
    cr4->SetBottomMargin(0.02);
    cr4->SetTopMargin(0.01);
    cr4->SetLeftMargin(0.17);
    cr4->SetRightMargin(0.02);
    cr4->Divide(2,3);

    //-------------------------------
    //   DRAW d
    //-------------------------------
    cr4->cd(1);
    gPad->SetLogy();
    gPad->SetLogx();
    gPad->SetTicky();
    gPad->SetTickx();
    hframe->Draw();
    hB2_coalescence->Draw("l");
    gBlastB2vsR_PbPb276TeV->Draw("samel");
    //pavept->Draw();
    gB2vsR_pp7TeVINELg0_sys[ip]->Draw("p3");
    gB2vsR_pp7TeVINELg0[ip]->Draw("samep");
    gB2vsR_PbPb276TeV_sys[ip]->Draw("p3");
    gB2vsR_PbPb276TeV[ip]->Draw("samep");
    paveLab2->Draw();
    masterLeg->Draw();
    legB3data->Draw();  
    //-------------------------------
    //   DRAW 3H
    //-------------------------------
    cr4->cd(2);
    gPad->SetLogy();
    gPad->SetLogx();
    gPad->SetTicky();
    gPad->SetTickx();
    hframe3t->Draw();
    hB3t_coalescence->Draw("l");
    gBlastB3vsR_PbPb276TeV->Draw("samel");
    //paveptB3->Draw();
    paveLab3t->Draw();
    masterLeg3t->Draw();

    //-------------------------------
    //   DRAW 3He
    //-------------------------------
    cr4->cd(3);
    gPad->SetLogy();
    gPad->SetLogx();
    gPad->SetTicky();
    gPad->SetTickx();
    hframe3->Draw();
    hB3_coalescence->Draw("l");
    gBlastB3vsR_PbPb276TeV->Draw("samel");
    gB3vsR_pp7TeV_sys[ip]->Draw("p3");
    gB3vsR_pp7TeV[ip]->Draw("samep");
    gB3vsR_PbPb276TeV_sys[ip]->Draw("p3");
    gB3vsR_PbPb276TeV[ip]->Draw("samep");
    //paveptB3->Draw();
    paveLab3->Draw();
    masterLeg3->Draw();

    //-------------------------------
    //   DRAW 3LH
    //-------------------------------
    cr4->cd(4);
    gPad->SetLogy();
    gPad->SetLogx();
    gPad->SetTicky();
    gPad->SetTickx();
    hframe3L->Draw();
    hB3L_coalescence_largeradius->Draw("samel");
    hB3L_coalescence->Draw("l");
    gBlastB3LambdavsR_PbPb276TeV->Draw("samel");
    gB3LambdavsR_PbPb276TeV_sys[ip]->Draw("p3");
    gB3LambdavsR_PbPb276TeV[ip]->Draw("samep");
    //paveptB3L->Draw();
    paveLab3L->Draw();
    masterLeg3L->Draw();

    //-------------------------------
    //   DRAW 4He
    //-------------------------------
    cr4->cd(5);
    gPad->SetLogy();
    gPad->SetLogx();
    gPad->SetTicky();
    gPad->SetTickx();
    hframe4->Draw();
    hB4_coalescence->Draw("samel");
    gBlastB4vsR_PbPb276TeV->Draw("samel");
    paveLab4->Draw();
    //pavept->Draw();
    masterLeg4->Draw();

    //-------------------------------
    //   DRAW 4LH
    //-------------------------------
    cr4->cd(6);
    gPad->SetLogy();
    gPad->SetLogx();
    gPad->SetTicky();
    gPad->SetTickx();
    hframe4L->Draw();
    hB4L_coalescence_largeradius->Draw("samel");
    hB4L_coalescence->Draw("samel");
    gBlastB4LambdavsR_PbPb276TeV->Draw("samel");
    //paveptB4L->Draw();
    paveLab4L->Draw();
    masterLeg4L->Draw();

    cr4->SaveAs(Form("%s/coal2Thermal2alice.pdf", figPath.Data()));
    cr4->SaveAs(Form("%s/coal2Thermal2alice.eps", figPath.Data()));
    cr4->SaveAs(Form("%s/coal2Thermal2alice.png", figPath.Data()));

    TFile * fout = new TFile(Form("%s/BApredictions_Epiphany.root",figPath.Data()),"RECREATE");
    fout->cd();
    hframe3->Write();
    hB3_coalescence->Write();
    gBlastB3vsR_PbPb276TeV->SetName("gBlastB3vsR_PbPb276TeV");
    gBlastB3vsR_PbPb276TeV->Write();
    hframe4->Write();
    hB4_coalescence->Write();
    gBlastB4vsR_PbPb276TeV->SetName("gBlastB4vsR_PbPb276TeV");
    gBlastB4vsR_PbPb276TeV->Write();
    //cr4->Write();

    return;
}

