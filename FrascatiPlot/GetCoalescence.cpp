/////////////////////////////////////////////////////////////////////////////////////////////////////////
//The Frascati project - coalescence predictions up to A=4
/////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "TGraphErrors.h"
#include "TF1.h"

Double_t getB2fromRadius(Double_t homogR, Double_t mT, Double_t objSize);
Double_t getB3fromRadius(Double_t homogR, Double_t mT, Double_t objSize);
Double_t getB4fromRadius(Double_t homogR, Double_t mT, Double_t objSize);

TF1 *          MakeB2TheoryGraphQMfactor(Double_t objSize = 3.2);
TGraphErrors * MakeB2TheoryGraphCoalescence(Double_t mT = 1.0, Double_t objSize = 3.2);
TGraphErrors * MakeB3TheoryGraphCoalescence(Double_t mT = 1.0, Double_t objSize = 2.48);
TGraphErrors * MakeB4TheoryGraphCoalescence(Double_t mT = 1.0, Double_t objSize = 1.9);


Double_t getB2fromRadius(Double_t homogR, Double_t mT, Double_t objSize)
{
 
  // formula 6.3, 4.12 from Scheibl, Heinz 1999 paper arXiv:nucl-th/9809092v2
  Double_t Rpar = homogR;
  Double_t Rperp = homogR;
  Double_t convFactor_fm2InvGeV = 0.197;

  Double_t invCd = (1.0 + TMath::Power( objSize / (2.0 * Rperp ), 2.0) ) *
		     TMath::Power((1.0 + TMath::Power( objSize / (2.0 * Rpar ), 2.0) ), 0.5);

  Double_t invCdApprox = TMath::Power( (1.0 + TMath::Power( objSize / (2.0 * Rperp ), 2.0)), 1.5);
  
  Double_t Cd = 1.0/invCd;
  Double_t B2 = (3.0 * TMath::Power(TMath::Pi(), 1.5) * Cd ) /
    (2 * mT * Rpar /convFactor_fm2InvGeV * Rperp /convFactor_fm2InvGeV * Rperp / convFactor_fm2InvGeV) ; 

  return B2;

}


Double_t getB3fromRadius(Double_t homogR, Double_t mT, Double_t objSize)
{
  // formula 9 of K. Blum, PRD 96 (2017) 103021
  // https://journals.aps.org/prd/pdf/10.1103/PhysRevD.96.103021
  Double_t R1 = homogR;
  Double_t R2 = homogR;
  Double_t R3 = homogR;

  Double_t convFactor_fm2InvGeV = 0.197;

  Double_t invC3Approx = (1.0 + TMath::Power( objSize / (2.0 * R1), 2.0)) *
    (1.0 + TMath::Power( objSize / (2.0 * R2), 2.0)) *
    (1.0 + TMath::Power( objSize / (2.0 * R3), 2.0));
  
  Double_t C3 = 1.0/invC3Approx;
  Double_t B3 = (TMath::Power(2* TMath::Pi(), 3.0) * C3) *
    TMath::Power(mT * R1 /convFactor_fm2InvGeV * R2 /convFactor_fm2InvGeV * R3 / convFactor_fm2InvGeV, -2.0) /
    ( 4 * TMath::Sqrt(3.0));
  
  return B3;
}

Double_t getB4fromRadius(Double_t homogR, Double_t mT, Double_t objSize)//He4
{
  // formula 13 of the Frascati paper arXiv:1807.05894
  Double_t convFactor_fm2InvGeV = 0.197;

  //calculation for 4He (J_A = 0)
  Double_t B4 = 1./32. * TMath::Power(mT, -3.) *
    TMath::Power(2* TMath::Pi(), 9./2.) *
    TMath::Power((homogR * homogR +  objSize * objSize / 4.)/(convFactor_fm2InvGeV * convFactor_fm2InvGeV), -9./2.);
  
  return B4;
}


TF1 * MakeB2TheoryGraphQMfactor(Double_t objSize)
{
  TF1 * funcCd = new TF1(Form("funcCd_%i", TMath::Nint(objSize*10)), "1 / TMath::Power(1 + [0]*[0]/(4*x*x), 1.5)", 0., 15.);
  funcCd->SetParameter(0, objSize);
  return funcCd;
}


TGraphErrors * MakeB2TheoryGraphCoalescence(Double_t mT, Double_t objSize)
{
  const Int_t nPoints = 1000;
  Double_t gY[nPoints];
  Double_t gR[nPoints];
  
  TGraphErrors * graphOut = new TGraphErrors(nPoints);
  graphOut->SetName("B2_th_coalescence");
  graphOut->SetTitle("B_{2} from coalescence");
  
  
  for (int i = 0; i<nPoints; i++){
    gR[i] = 10.0 * i / nPoints; 
    gY[i] = getB2fromRadius(gR[i], mT, objSize);

    graphOut->SetPoint(i, gR[i], gY[i]);
    graphOut->SetPointError(i, 0.0, 0.0);
  }

  return graphOut;

}


TGraphErrors * MakeB3TheoryGraphCoalescence(Double_t mT, Double_t objSize)
{
  const Int_t nPoints = 1000;
  Double_t gY[nPoints];
  Double_t gR[nPoints];
  
  TGraphErrors * graphOut = new TGraphErrors(nPoints);
  graphOut->SetName("B3_th_coalescence");
  graphOut->SetTitle("B_{3} from coalescence");
  
  
  for (int i = 0; i<nPoints; i++){
    gR[i] = 10.0 * i / nPoints; 
    gY[i] = getB3fromRadius(gR[i], mT, objSize);

    graphOut->SetPoint(i, gR[i], gY[i]);
    graphOut->SetPointError(i, 0.0, 0.0);
  }

  return graphOut;

}

TGraphErrors * MakeB4TheoryGraphCoalescence(Double_t mT, Double_t objSize)
{
  const Int_t nPoints = 1000;
  Double_t gY[nPoints];
  Double_t gR[nPoints];
  
  TGraphErrors * graphOut = new TGraphErrors(nPoints);
  graphOut->SetName("B4_th_coalescence");
  graphOut->SetTitle("B_{4} from coalescence");
  
  
  for (int i = 0; i<nPoints; i++){
    gR[i] = 10.0 * i / nPoints; 
    gY[i] = getB4fromRadius(gR[i], mT, objSize);

    graphOut->SetPoint(i, gR[i], gY[i]);
    graphOut->SetPointError(i, 0.0, 0.0);
  }

  return graphOut;

}
