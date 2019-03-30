//mapping multiplicity into radius
void getRadiusFromParameterisation(Double_t * multi = 0x0, Double_t * radius = 0x0, Int_t paramSet = 1)
{
  // parameterisation to convert multiplicity into Rside
  // this is the parameterisation for the pion radius
  //
  if (!multi || !radius) return;
  Double_t  multi3 = TMath::Power(multi[0], 1./3.);  
  //
  // Here is the crucial mapping between HBT radii and multi^(1/3)
  //
  // VERSION (5th February 2018):
  // We assume 0.85fm for pp as Kfir (at dNdeta 6.01)
  // We assume R=4.5fm for central Pb-Pb based on arXiv:1012.4035 (figure 2)
  // We interpolate linearly and take the highest kT bin, because
  // a pT/A of 0.8 GeV/c corresponds to mT=1.23GeV which is a kT of pions of 1.2GeV
  // radiusVal = 0.177825 + 0.36733 * multi3;

  // VERSION (17th May 2018):
  // We fit linearly the ALICE data at the kT = 0.887 
  Double_t radiusVal = 0.0;
  if (paramSet == 3){
    //manual hack to reproduce Donigus, Ko arXiv:1812.05175 where HBT parameter is based on results for kT = 0.25
    radiusVal = 0.0 + 0.83 * multi3;
  } else if (paramSet==2) {
    //manual hack to have the data points fall onto the U. Heinz curve for 3He
    radiusVal = 0.190 + 0.380 * multi3;
  } else  if (paramSet==1) {
    //manual hack to have the data points fall onto the U. Heinz curve for deuteron
    radiusVal = 0.0 + 0.472949 * multi3; //radius = 0 at 0 dN/deta
    //radiusVal = 0.07412 + 0.46637 * multi3; //radius = 0.85 fm for pp, dN/deta = 4.60 (INELg0)
    //radiusVal = -0.009 + 0.4738 * multi3; //radius = 0.85 fm for pp, dN/deta = 5.98 (INELg0>0)
    //radiusVal = -0.3949 + 0.507865 * multi3; //most central and peripheral PbPb 2.76 TeV on the curve
  } else {
    //fit to the HBT data, kT = 0.887
    radiusVal = 0.128 + 0.339 * multi3; 
  }
  
  //
  //
  // OLD versions:
  // we consider Rside as a proxy for Rinv
  // the Rside vs dN/detaˆ1/3 is take from figure 9 of http://aliceinfo.cern.ch/ArtSubmission/node/1183
  // input has to be dNch/deta at midrapidity in VO* bins
  //
  //  Double_t radiusVal = -0.750 + 0.625 * multi3; //0.18 to fall on the b3 in pp
  //
  //
  Double_t radiusErr = radiusVal / 3.0 * multi[1] / multi[0];
  radius[0] = radiusVal;
  radius[1] = radiusErr;
  return; 
}



void getMultiFromR(Double_t * multi = NULL, Double_t * radius = NULL, Int_t paramSet = 1)
{
  // Here is the crucial mapping between HBT radii and multi^(1/3)
  if (!multi || !radius) return;
  //
  // VERSION (17th May 2018):
  // We fit linearly the ALICE data at the kT = 0.887 
  Double_t radiusVal = radius[0];
  Double_t  multi3 = 0.0;

  if (paramSet==3) {
    //manual hack to reproduce Donigus, Ko arXiv:1812.05175 where HBT parameteris is based on results for kT = 0.25
    multi3 = radiusVal / 0.83;
  } else if (paramSet==2) {
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


void convertMultiToRadius(TGraphErrors * graph, Int_t paramSet)
{
  //convert multiplicity dN/deta into radius with parameterisation
  if (!graph) return;
  
  for (Int_t ip = 0; ip < graph->GetN(); ip++){
    Double_t xold[2];
    xold[0] = graph->GetX()[ip];
    xold[1] = graph->GetEX()[ip];
    Double_t xnew[2] = {-1.0, -1.0};
    getRadiusFromParameterisation(xold, xnew, paramSet);
    graph->GetX()[ip] = xnew[0];
    graph->GetEX()[ip] = xnew[1];
  }

  return;
}


void convertMultiToRadius(TGraphAsymmErrors * graph, Int_t paramSet)
{
  //convert multiplicity dN/deta into radius with parameterisation
  if (!graph) return;
  
  for (Int_t ip = 0; ip < graph->GetN(); ip++){
    Double_t xold[2];
    xold[0] = graph->GetX()[ip];
    xold[1] = graph->GetEXlow()[ip];
    Double_t xnew[2] = {-1.0, -1.0};
    getRadiusFromParameterisation(xold, xnew, paramSet);
    graph->GetX()[ip] = xnew[0];
    graph->GetEXlow()[ip] = xnew[1];
    graph->GetEXhigh()[ip] = xnew[1];
	
  }

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
