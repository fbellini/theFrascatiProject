#include "Rtypes.h"
#include "TF1.h" 
#include "TMath.h"
#include "TCanvas.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//AK implementation of the blast-wave function
/////////////////////////////////////////////////////////////////////////////////////////////////////////
Double_t blast_integrand(const Double_t *x,const Double_t *par){
  Double_t x0=x[0];
  Double_t m=par[0];
  Double_t t=fabs(par[1]);
  Double_t n=par[2];
  Double_t beta_max=par[3];
  Double_t pt=par[4];

  //Keep beta within reasonable limits.
  Double_t beta=beta_max*TMath::Power(x0,n);
  if(beta>0.9999999999999999) beta=0.9999999999999999;

  Double_t mt=TMath::Sqrt(m*m+pt*pt);
  Double_t rho0=TMath::ATanH(beta);
  Double_t a0=pt*TMath::SinH(rho0)/t;
  if(a0>700.) a0=700.;//avoid floating point exception
  Double_t a1=mt*TMath::CosH(rho0)/t;
  return x0*mt*TMath::BesselI0(a0)*TMath::BesselK1(a1);
}

Double_t blast(Double_t *x,Double_t *par){
  static TF1* fint=0;
  if(!fint) fint=new TF1("fint",blast_integrand,0.,1.,5);
  
  fint->SetParameters(par[0],par[2],par[3],par[4],x[0]);
  return fint->Integral(0.,1.)*par[1];
}

Double_t  x_blast(Double_t *x,Double_t *par){
  return x[0]*blast(x,par);
}

void get_average(Double_t &bt, Double_t &bn, Double_t &bb)
{
  bt = (0.143*5 + 0.147*5 + 0.151*10) / 20.0;
  bn = (1.07*5 + 1.14*5 + 1.24*10) / 20.0;
  bb = (0.547*5 + 0.531*5 + 0.511*10) / 20.0;
  return;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Blast-wave parameters from ALICE fits to π,K,p
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void GetParams_pp7TeV(Double_t *parT, Double_t *parTerr, Double_t *parN, Double_t *parNerr, Double_t *parB, Double_t *parBerr, Double_t *dNdeta, Double_t *dNdetaerr, Double_t *protdNdy, Double_t *protdNdyErr, Double_t *piondNdy, Double_t *piondNdyErr)
{    
  //Blast-wave parameters for different multiplicity percentile bins
  //bt is temperature, bn is n power
  //bb is <beta_T>
  //K kaon yield, KT statistical uncertainty, KY systematic uncertainty

 //pp 7 TeV 
  //values received by Livio for pp 7 TeV fit to π,K,p
  //Ranges: pi: 0.5-1, k: 0.2-1.5, p: 0.3-3.0
  Double_t bt[10]    = {0.1683, 0.1676, 0.1696, 0.1722, 0.1737, 0.1740, 0.1744, 0.1758, 0.1735, 0.1650};
  Double_t bterr[10] = {0.0090, 0.0078, 0.0074, 0.0070, 0.0068, 0.0064, 0.0059, 0.0054, 0.0051, 0.0043};

  Double_t bn[10]    = {1.4725, 1.6749, 1.9076, 2.1583, 2.4063, 2.7026, 3.2686, 4.0937, 5.6110, 10.7011};
  Double_t bnerr[10] = {0.1172, 0.1356, 0.1688, 0.2084, 0.2488, 0.2950, 0.3868, 0.5411, 0.8396,  2.0103};
  
  Double_t bb[10]    = {0.4647, 0.4336, 0.4024, 0.3741, 0.3517, 0.3262, 0.2903, 0.2496, 0.2010, 0.1215};
  Double_t bberr[10] = {0.0151, 0.0150, 0.0160, 0.0169, 0.0176, 0.0179, 0.0184, 0.0187, 0.0187, 0.0163};

  Double_t multi[10]     = {21.3, 16.5, 13.5, 11.5, 10.1, 8.45, 6.72, 5.40, 3.90, 2.26};
  Double_t multierr[10]  = { 0.6,  0.5,  0.4,  0.3,  0.3, 0.25, 0.21, 0.17, 0.14, 0.12};

  //(p+pabar)/2
  Double_t dNdyProt[10]    = {0.5488, 0.4369, 0.3599, 0.3106, 0.2741, 0.2316, 0.1860, 0.1491, 0.1070, 0.05856};
  Double_t dNdyErrProt[10] = {0.0393, 0.0301, 0.0246, 0.0210, 0.0185, 0.0156, 0.0125, 0.0101, 0.0080, 0.00532};

  // (pion+ + pion-)/2. --> values from the famous long paper
  Double_t dNdyPion[10]    = {10.035, 7.878, 6.459, 5.554, 4.982, 4.138, 3.326, 2.699, 1.989, 1.210}; //HERE-I-STOP
  Double_t dNdyErrPion[10] = {0.519,  0.377, 0.306, 0.261, 0.228, 0.192, 0.153, 0.125, 0.103, 0.086};

  
  for (Int_t j=0; j<10; j++){
    parT[j] = bt[j];
    parTerr[j] = bterr[j];
    parN[j] = bn[j];
    parNerr[j] = bnerr[j];
    parB[j] = bb[j];
    parBerr[j] = bberr[j];
    dNdeta[j] = multi[j];
    dNdetaerr[j] = multierr[j];
    protdNdy[j] = dNdyProt[j];
    protdNdyErr[j] = dNdyErrProt[j];
    piondNdy[j] = dNdyPion[j];
    piondNdyErr[j] = dNdyErrPion[j];
  }
  return;
}


void GetParams_pPb502TeV(Double_t *parT, Double_t *parTerr, Double_t *parN, Double_t *parNerr, Double_t *parB, Double_t *parBerr, Double_t *dNdeta, Double_t *dNdetaerr, Double_t *protdNdy, Double_t *protdNdyErr, Double_t *piondNdy, Double_t *piondNdyErr)
{    
  //Blast-wave parameters for different multiplicity percentile bins
  //bt is temperature, bn is n power
  //bb is <beta_T>
  //Ref. for parameters and dN/deta: Physics Letters B 728 (2014) 25–38
  // errors on the parameters are only systematics and symmetrised between +sys and -sys using the largest of the two
  Double_t bt[7]    = {0.143, 0.147, 0.151, 0.157, 0.164, 0.169, 0.166};
  Double_t bterr[7] = {0.010, 0.010, 0.020, 0.020, 0.020, 0.020, 0.020};

  Double_t bn[7]    = {1.07, 1.14, 1.24, 1.41, 1.73, 2.4, 3.9};
  Double_t bnerr[7] = {0.09, 0.20, 0.20, 0.20, 0.40, 0.6, 0.7};
  
  Double_t bb[7]    = {0.547, 0.531, 0.511, 0.478, 0.428, 0.360, 0.260};
  Double_t bberr[7] = {0.020, 0.030, 0.030, 0.030, 0.030, 0.040, 0.030};
  
  Double_t multi[7]    = { 45.0, 36.2, 30.5, 23.2, 16.1, 9.8, 4.4};
  Double_t multierr[7] = {  1.0,  0.8,  0.7,  0.5,  0.4, 0.2, 0.1};

  //p+pbar
  Double_t dNdyProt[7]    = {2.280446, 1.853576, 1.577375, 1.221875, 0.8621290, 0.5341445, 0.2307767};
  Double_t dNdyErrProt[7] = {1.585968e-01, 1.284426e-01, 1.084793e-01, 8.190972e-02, 5.608663e-02, 3.408551e-02, 1.450569e-02};

  // pions
  Double_t dNdyPion[7]     = {4.080698e+01, 3.308405e+01, 2.805619e+01, 2.174798e+01, 1.529001e+01, 9.455217e+00, 4.313766e+00};
  Double_t dNdyErrPion[7]  = {1.985452e+00, 1.586469e+00, 1.329896e+00, 1.010231e+00, 7.081661e-01, 4.401432e-01, 1.950208e-01};



   
  for (Int_t j=0; j<7; j++){
    parT[j] = bt[j];
    parTerr[j] = bterr[j];
    parN[j] = bn[j];
    parNerr[j] = bnerr[j];
    parB[j] = bb[j];
    parBerr[j] = bberr[j];
    dNdeta[j] = multi[j];
    dNdetaerr[j] = multierr[j];
    //return (p+pbar)/2
    protdNdy[j] = dNdyProt[j]/2.;
    protdNdyErr[j] = dNdyErrProt[j]/2.;
    //return (pi+ + pi-)/2
    piondNdy[j] = dNdyPion[j]/2.;
    piondNdyErr[j] = dNdyErrPion[j]/2.;
  }
  return;
}

void GetParams_PbPb276TeV(Double_t *parT, Double_t *parTerr, Double_t *parN, Double_t *parNerr, Double_t *parB, Double_t *parBerr, Double_t *dNdeta, Double_t *dNdetaerr, Double_t *protdNdy, Double_t *protdNdyErr, Double_t *piondNdy, Double_t *piondNdyErr, Double_t *lambdadNdy, Double_t *lambdadNdyErr)
{    
  //Blast-wave parameters for different multiplicity percentile bins
  //bt is temperature, bn is n power
  //bb is <beta_T>
  //Ref. for parameters: Phys. Rev. C 88 (2013) 044910 
  //Ref for dN/deta: Phys. Rev. C 88 (2013) 044910 
  // errors on the parameters are only systematics

  Double_t bt[10]    = {0.095, 0.097, 0.099, 0.101, 0.106, 0.112, 0.118, 0.129, 0.139, 0.151};
  Double_t bterr[10] = {0.010, 0.011, 0.011, 0.012, 0.012, 0.013, 0.014, 0.017, 0.027, 0.044};

  Double_t bn[10]    = {0.712, 0.723, 0.738, 0.779, 0.841, 0.944, 1.099, 1.292, 1.578, 2.262};
  Double_t bnerr[10] = {0.086, 0.116, 0.118, 0.133, 0.168, 0.142, 0.187, 0.194, 0.205, 0.498};
  
  Double_t bb[10]    = {0.651, 0.646, 0.639, 0.625, 0.604, 0.574, 0.535, 0.489, 0.438, 0.357};
  Double_t bberr[10] = {0.020, 0.023, 0.022, 0.025, 0.022, 0.016, 0.018, 0.024, 0.039, 0.084};
  
  Double_t multi[10]    = {1601, 1294, 966, 649, 426, 261, 149, 76, 35, 13.4 };
  Double_t multierr[10] = {60, 49, 37, 23, 15, 9, 6, 4, 2, 1.6};

  // (p+pbar)/2.
  Double_t dNdyProt[10]    = {  33,   28,  21.1,  14.5,  9.7,    6.2,   3.7,   2.0,   0.9, 0.36};
  Double_t dNdyErrProt[10] = {   3,    2,   1.8,   1.2,  0.8,    0.5,   0.3,   0.2,  0.08, 0.04};

  // pions
  Double_t dNdyPion[10]     = {733., 606., 455., 307., 201., 124., 71., 37., 17.1, 6.6};
  Double_t dNdyErrPion[10]  = { 54.,  42.,  31.,  20.,  13.,   8.,  5.,  2.,  1.1, 0.4};

  // Lambda only -- table 23 of http://hepdata.cedar.ac.uk/view/ins1243863
  // --> TODO: deal with different binning by interpolation
  // published bins are          (0-5,  5-10, 10-20, 20-40, fixme, 40-60, fixme, 60-80, fixme,  80-90)%
  Double_t dNdyLambda[10]    = {25.61, 21.58, 16.89, 9.79,  9.79,  3.79, 3.79,  0.99,   0.99,   0.214};
  Double_t dNdyErrLambda[10] = { 2.91,  1.9,   1.9,  0.947, 0.947, 0.378,0.378, 0.0937, 0.0937, 0.0166};



  for (Int_t j=0; j<10; j++){
    parT[j] = bt[j];
    parTerr[j] = bterr[j];
    parN[j] = bn[j];
    parNerr[j] = bnerr[j];
    parB[j] = bb[j];
    parBerr[j] = bberr[j];
    dNdeta[j] = multi[j];
    dNdetaerr[j] = multierr[j];
    protdNdy[j] = dNdyProt[j];
    protdNdyErr[j] = dNdyErrProt[j];
    piondNdy[j] = dNdyPion[j];
    piondNdyErr[j] = dNdyErrPion[j];
    lambdadNdy[j] = dNdyLambda[j];
    lambdadNdyErr[j] = dNdyErrLambda[j];
  }
  return;
}


void GetParams_PbPb502TeV(Double_t *parT, Double_t *parTerr, Double_t *parN, Double_t *parNerr, Double_t *parB, Double_t *parBerr, Double_t *dNdeta, Double_t *dNdetaerr , Double_t *protdNdy, Double_t *protdNdyErr,  Double_t *piondNdy, Double_t *piondNdyErr, Double_t *lambdaYield, Double_t *lambdaYieldErr)
{    
  //Blast-wave parameters for different multiplicity percentile bins
  //bt is temperature, bn is n power
  //bb is <beta_T>
  
  //BW Results from Livio (email 27/01/2017): fit of the spectra from 24.01.2017 (approval forum, iteration 1)
  //Parameters and proton yield can be found at:
  //https://gitlab.cern.ch/njacazio/SpectraAnalysisRun2/tree/master/results/spectra/spectra-pag/Preliminaries/QM2017

  Double_t bt[10]    = {0.086, 0.088, 0.091, 0.094, 0.099, 0.106, 0.104, 0.127, 0.147, 0.163};
  Double_t bterr[10] = {0.004, 0.004, 0.004, 0.004, 0.004, 0.005, 0.005, 0.006, 0.007, 0.008};

  Double_t bn[10]    = {0.751, 0.753, 0.752, 0.778, 0.824, 0.890, 1.027, 1.213, 1.627, 2.254};
  Double_t bnerr[10] = {0.017, 0.016, 0.016, 0.017, 0.019, 0.024, 0.031, 0.053, 0.106, 0.278};
  
  Double_t bb[10]    = {0.663, 0.659, 0.655, 0.643, 0.625, 0.600, 0.563, 0.515, 0.440, 0.360};
  Double_t bberr[10] = {0.004, 0.004, 0.004, 0.004, 0.004, 0.005, 0.006, 0.009, 0.014, 0.022};
  
  Double_t multi[10]    = {1942.5, 1585.5, 1180.0, 786.0, 512.0, 318.0, 183.0, 96.3, 44.9, 17.52};
  Double_t multierr[10] = {53.5, 46.0, 31.0, 20.0, 15.0, 12.0, 8.0, 5.8, 3.4, 1.84};

  //p+pbar 
  Double_t dNdyProt[10]    = {74.5852, 61.2207, 47.2026, 33.0090, 22.3932, 14.3862, 8.6260, 4.7102, 2.2640, 0.8936};
  Double_t dNdyErrProt[10] = { 5.0590,  3.8429,  3.0289,  2.0477,  1.3414,  0.8648, 0.5221, 0.3186, 0.1516, 0.0651};

  //pi+  + pi- 
  Double_t dNdyPion[10]    = {1702.3, 1376.59, 1037.78, 710.13, 466.401, 291.383, 170.758, 88.6999, 41.6013, 16.2665};
  Double_t dNdyErrPion[10] = {87.4099, 71.5484, 51.191, 32.7179, 21.4739, 13.3082, 8.10989, 4.42991, 2.05458, 0.840404};

  //Lambda + antiLambda - from David 13.09.2018
  Double_t dNdyLambda[10]    = {64.8, 52.9, 41.30, 28.84, 19.44, 12.28, 7.18, 3.74, 1.75, 0.714};
  Double_t dNdyErrLambda[10] = { 6.0,  5.0,  3.45,  2.38,  1.55,  0.96, 0.56, 0.30, 0.15, 0.076};


  for (Int_t j=0; j<10; j++){
    parT[j] = bt[j];
    parTerr[j] = bterr[j];
    parN[j] = bn[j];
    parNerr[j] = bnerr[j];
    parB[j] = bb[j];
    parBerr[j] = bberr[j];
    dNdeta[j] = multi[j];
    dNdetaerr[j] = multierr[j];
    //return (p+pbar)/2
    protdNdy[j] = dNdyProt[j]/2.;
    protdNdyErr[j] = dNdyErrProt[j]/2.;  
    //return (pi+ + pi-)/2
    piondNdy[j] = dNdyPion[j]/2.;
    piondNdyErr[j] = dNdyErrPion[j]/2.;
    lambdaYield[j] = dNdyLambda[j]/2;
    lambdaYieldErr[j] = dNdyErrLambda[j]/2;
  }
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Generate predictions for BA
////////////////////////////////////////////////////////////////////////////////////////////////////////

TGraphAsymmErrors * generateBWpredictionsB2(TString system = "PbPb276TeV",  TString errType = "rms", TString particle = "deuteron", Double_t pToA = 0.)
{
  //
  // Get the blast-wave function for a given particle type based on the pi/K/p parameters.
  // The normalisation is such that the pT-integrated yields of (anti-nuclei)
  // correspond to the thermal model expectation based on a fixed d/p or 3He/p ratio.
  //
  const Double_t dOverPthermal = 0.00294824; // GSI-Heidelberg at 156 MeV
  const Double_t dOverPiThermal = 0.0001817; // GSI-Heidelberg at 156 MeV
  const Double_t He3OverPiThermal = 5.354e-07; // GSI-Heidelberg at 156 MeV
  const Double_t He3OverPthermal =  8.68733e-06; // GSI-Heidelberg at 156 MeV --- (only 4% difference wrt to old value of 0.00294824/330.)
  const Double_t s3thermal = 0.55; // TODO: this is just read off figure 7 og hyper-triton paper.
  const Double_t mProton = 0.938;
  const Double_t mLambda = 1.115; // FIXME: use proper mass value
  const Double_t He4thermal = 7.e-7; //read off plot 
  const Double_t H4Lthermal = 2.e-7; //read off plot
  //Sed deuteron mass
  Double_t m = 0.0;
  
  //CODATA ref: http://www.codata.org/uploads/RMP.88.035009.pdf
  if (particle.Contains("deuteron")) m = 1.875612928; //CODATA, https://physics.nist.gov/cgi-bin/cuu/Value?mdc2mev
  if (particle.Contains("triton")) m = 2.808921112; //triton mass 
  if (particle.Contains("He3")) m = 2.808921112 - 0.018592017; //using mass difference (m3H - m3He) from CODATA
  if (particle.Contains("hyper-triton")) m = 2.991; //FIXME: use a proper values, but this is correct within two permille 2*0.938+1.115
  if (particle.Contains("He4")) m = 3.72737937823; //using value from CODATA for alpha particle
  if (particle.Contains("4LH")) m = 3.931; 
  //centrality bin
  Int_t cb = -1;
 
  Color_t color[10] = {kRed+2, kRed, kOrange, kYellow+1, kSpring+2, kGreen+1, kAzure+7, kBlue, kViolet-3, kBlue+2};
  Char_t mclass[10][5] = {"I", "II", "III", "IV", "V","VI","VII","VIII","IX","X"};
   
  Double_t bt[10];
  Double_t bterr[10];

  Double_t bn[10];
  Double_t bnerr[10];
  
  Double_t bb[10];
  Double_t bberr[10];

  Double_t multi[10];
  Double_t multiErr[10];

  Double_t protYield[10];
  Double_t protYieldErr[10];

  Double_t pionYield[10];
  Double_t pionYieldErr[10];

  Double_t lambdaYield[10];
  Double_t lambdaYieldErr[10];
  //-------------------------------
  //arrays defined above are filled here usig the proper function and depending on the content of the string 'system'
  //-------------------------------

  if (system.Contains("PbPb276TeV"))
    GetParams_PbPb276TeV(bt, bterr, bn, bnerr, bb, bberr, multi, multiErr, protYield, protYieldErr, pionYield, pionYieldErr, lambdaYield, lambdaYieldErr);
  
  if (system.Contains("PbPb502TeV"))
    GetParams_PbPb502TeV(bt, bterr, bn, bnerr, bb, bberr, multi, multiErr, protYield, protYieldErr, pionYield, pionYieldErr, lambdaYield, lambdaYieldErr);
  
  if (system.Contains("pPb502TeV"))
    GetParams_pPb502TeV(bt, bterr, bn, bnerr, bb, bberr, multi, multiErr, protYield, protYieldErr, pionYield, pionYieldErr);

  if (system.Contains("pp7TeV"))
    GetParams_pp7TeV(bt, bterr, bn, bnerr, bb, bberr, multi, multiErr, protYield, protYieldErr, pionYield, pionYieldErr);

  //-------------------------------
  // coalescence parameter values
  //-------------------------------
  Double_t B2blastW[10]       = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t B2blastWErrLow[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t B2blastWErrUp[10]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  //
  Int_t binCounter = 0;
  //
  for (Int_t j=0; j<10; j++){
    
    if ((cb>0) && (j!=cb)) continue;

    if (system.Contains("pPb502TeV") && j>6) continue; //we have only 7 points for p-Pb
    
    //bb is changed from <beta_T> to beta_S
    bb[j] *= (0.5*(bn[j]+2.));
    
    //-----------------------------------
    // Construct Blast-Wave model protons and Lambda
    //-----------------------------------
    TF1 * fProton = new TF1(Form("BlastWaveProton-%s", mclass[j]), x_blast, 0., 10., 5);
    fProton->SetParameters(mProton, 1., bt[j], bn[j], bb[j]);//m is mass, yield is yield in min bias
    fProton->SetLineColor(color[j]);
    Double_t Iproton = fProton->Integral(0.,30.);
    if (Iproton>0) fProton->SetParameter(1, protYield[j]/Iproton); // normalisation to match yields

    TF1 * fLambda = new TF1(Form("BlastWaveLambda-%s", mclass[j]), x_blast, 0., 10., 5);
    fLambda->SetParameters(mLambda, 1., bt[j], bn[j], bb[j]);//m is mass, yield is yield in min bias
    fLambda->SetLineColor(color[j]);
    Double_t Ilambda = fLambda->Integral(0.,30.);
    if (Ilambda>0) fLambda->SetParameter(1, protYield[j]*(lambdaYield[j]/protYield[j])/Ilambda); // normalisation to match yields
   
    //-----------------------------------
    // Construct Blast-Wave model nuclei
    //-----------------------------------
    TF1 * fNucleus = new TF1(Form("BlastWave%s-%s", particle.Data(), mclass[j]), x_blast, 0., 10., 5);
    fNucleus->SetParameters(m, 1., bt[j], bn[j], bb[j]);//m is mass, yield is yield in min bias
    fNucleus->SetLineColor(color[j]);
    Double_t Inucleus = fNucleus->Integral(0.,30.); //don't we want here only up to 10 GeV/c? effect of tail seen in Manuel's results
    //if (Inucleus>0 && particle.Contains("deuteron")) fNucleus->SetParameter(1, protYield[j]*dOverPthermal/Inucleus); // normalisation to match yields
    //if (Inucleus>0 && particle.Contains("He3")) fNucleus->SetParameter(1, protYield[j]*He3OverPthermal/Inucleus); // normalisation to match yields
    if (Inucleus>0 && particle.Contains("deuteron")) fNucleus->SetParameter(1, pionYield[j]*dOverPiThermal/Inucleus); // normalisation to match yields
    if (Inucleus>0 && particle.Contains("He3")) fNucleus->SetParameter(1, pionYield[j]*He3OverPiThermal/Inucleus); // normalisation to match yields
    if (Inucleus>0 && particle.Contains("hyper-triton")) {
      //Double_t he3Yield = protYield[j]*He3OverPthermal;
      Double_t he3Yield = pionYield[j]*He3OverPiThermal;
      Double_t hyperTritonYield = s3thermal*he3Yield*(lambdaYield[j]/protYield[j]); //now with centrality dependence
      fNucleus->SetParameter(1, hyperTritonYield/Inucleus); // normalisation to match yields
    }
    if (Inucleus>0 && particle.Contains("He4")) fNucleus->SetParameter(1, He4thermal/Inucleus); // normalisation to match yields
    if (Inucleus>0 && particle.Contains("4LH")) fNucleus->SetParameter(1, H4Lthermal/Inucleus); // normalisation to match yields

    
    // Printf(Form("-------------------------------------\nV0M multiplicity class: %s (cb = %i)", mclass[j], j));
    // Printf("Blast-Wave parameters (pi,K,p):\n  T = %6.4f \n  n = %6.4f \n  beta_T = %6.4f", bt[j], bn[j], bb[j]);

    //-----------------------------------
    // do the actual BA calculation
    //-----------------------------------
    Int_t A = 1;
    if (particle.Contains("deuteron")) A = 2;
    else if (particle.Contains("He3") || particle.Contains("hyper-triton")) A = 3;
    else if (particle.Contains("He4") || particle.Contains("4LH")) A = 4;
	
    Double_t invYieldProton = fProton->Eval(pToA)/(2*TMath::Pi()*pToA);
    Double_t invYieldLambda = fLambda->Eval(pToA)/(2*TMath::Pi()*pToA);
    Double_t invYieldNucleus = fNucleus->Eval(pToA*A)/(2*TMath::Pi()*pToA*A);

    if (particle.Contains("deuteron")) B2blastW[j] = invYieldNucleus / (invYieldProton * invYieldProton);
    if (particle.Contains("He3")) B2blastW[j] = invYieldNucleus / (invYieldProton * invYieldProton * invYieldProton);
    if (particle.Contains("hyper-triton")) B2blastW[j] = invYieldNucleus / (invYieldLambda * invYieldProton * invYieldProton);
    if (particle.Contains("He4")) B2blastW[j] = invYieldNucleus / (invYieldProton * invYieldProton * invYieldProton * invYieldProton);
    if (particle.Contains("4LH")) B2blastW[j] = invYieldNucleus / (invYieldLambda * invYieldProton * invYieldProton * invYieldProton);
    binCounter++;

  }
  
  //output graph
  TGraphAsymmErrors * gB2 = new TGraphAsymmErrors(binCounter, multi, B2blastW, 0x0, 0x0, 0x0, 0x0);

  // TCanvas * model = new TCanvas("model","model", 800, 800);
  // model->cd();
  // gB2->Draw("apz");
  
  return gB2;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Generate predictions for nuclei spectra
////////////////////////////////////////////////////////////////////////////////////////////////////////

TH1D * generateBWpredictionSpectra(Int_t selcent = 0, TString system = "PbPb276TeV",  TString particle = "proton", TString errType = "rms")
{
  //
  // Get the blast-wave for proton, based on the pi/K/p parameters.
  // The normalisation is such that the pT-integrated yield of protons corresponds to the measured yield
  const Double_t mProton = 0.938;
  const Double_t mLambda = 1.115; // FIXME: use proper mass value
  //Set mass
  Double_t m = 0.0;
  
  //CODATA ref: http://www.codata.org/uploads/RMP.88.035009.pdf
  if (particle.Contains("deuteron")) m = 1.875612928; //CODATA, https://physics.nist.gov/cgi-bin/cuu/Value?mdc2mev
  if (particle.Contains("triton")) m = 2.808921112; //triton mass 
  if (particle.Contains("He3")) m = 2.808921112 - 0.018592017; //using mass difference (m3H - m3He) from CODATA
  if (particle.Contains("hyper-triton")) m = 2.991; //FIXME: use a proper values, but this is correct within two permille 2*0.938+1.115
  if (particle.Contains("He4")) m = 3.72737937823; //using value from CODATA for alpha particle
  if (particle.Contains("4LH")) m = 3.931; //Ramona uses the same, from STEER/STEERBase/AliPDG.cxx

  //centrality bin
  Int_t cb = -1;
 
  Color_t color[10] = {kRed+2, kRed, kOrange, kYellow+1, kSpring+2, kGreen+1, kAzure+7, kBlue, kViolet-3, kBlue+2};
  Char_t mclass[10][5] = {"I", "II", "III", "IV", "V","VI","VII","VIII","IX","X"};
   
  Double_t bt[10];
  Double_t bterr[10];

  Double_t bn[10];
  Double_t bnerr[10];
  
  Double_t bb[10];
  Double_t bberr[10];

  Double_t multi[10];
  Double_t multiErr[10];

  Double_t protYield[10];
  Double_t protYieldErr[10];

  Double_t pionYield[10];
  Double_t pionYieldErr[10];

  Double_t lambdaYield[10];
  Double_t lambdaYieldErr[10];
  
  //-------------------------------
  //arrays defined above are filled here usig the proper function and depending on the content of the string 'system'
  //-------------------------------

  if (system.Contains("PbPb276TeV"))
    GetParams_PbPb276TeV(bt, bterr, bn, bnerr, bb, bberr, multi, multiErr, protYield, protYieldErr, pionYield, pionYieldErr, lambdaYield, lambdaYieldErr);
  
  if (system.Contains("PbPb502TeV"))
    GetParams_PbPb502TeV(bt, bterr, bn, bnerr, bb, bberr, multi, multiErr, protYield, protYieldErr, pionYield, pionYieldErr, lambdaYield, lambdaYieldErr);
  
  if (system.Contains("pPb502TeV"))
    GetParams_pPb502TeV(bt, bterr, bn, bnerr, bb, bberr, multi, multiErr, protYield, protYieldErr, pionYield, pionYieldErr);

  if (system.Contains("pp7TeV"))
    GetParams_pp7TeV(bt, bterr, bn, bnerr, bb, bberr, multi, multiErr, protYield, protYieldErr, pionYield, pionYieldErr);

  Int_t binCounter = 0;
  TH1D * hPro[10];
  TH1D * hLambda[10];
  
  for (Int_t j=0; j<10; j++){

    if ((cb>0) && (j!=cb)) continue;
    if (system.Contains("pPb502TeV") && j>6) continue; //we have only 7 points for p-Pb
    
    //bb is changed from <beta_T> to beta_S
    bb[j] *= (0.5*(bn[j]+2.));
    
    //-----------------------------------
    // Construct Blast-Wave model protons and Lambda
    //-----------------------------------
    TF1 * fProton = new TF1(Form("BlastWaveProton-%s", mclass[j]), x_blast, 0., 10., 5);
    fProton->SetParameters(mProton, 1., bt[j], bn[j], bb[j]);//m is mass, yield is yield in min bias
    fProton->SetLineColor(color[j]);
    Double_t Iproton = fProton->Integral(0.,10.);
    if (Iproton>0) fProton->SetParameter(1, protYield[j]/Iproton); // normalisation to match yields

    TF1 * fLambda = new TF1(Form("BlastWaveLambda-%s", mclass[j]), x_blast, 0., 10., 5);
    fLambda->SetParameters(mLambda, 1., bt[j], bn[j], bb[j]);//m is mass, yield is yield in min bias
    fLambda->SetLineColor(color[j]);
    Double_t Ilambda = fLambda->Integral(0.,10.);
    if (Ilambda>0) fLambda->SetParameter(1, lambdaYield[j]/Ilambda); // normalisation to match yields // FIXME --> proper LambdaYield

    hPro[j] = new TH1D(Form("hPro%i", j), "p spectrum", 60, 0., 6.);
    hLambda[j] = new TH1D(Form("hLambda%i", j), "#Lambda spectrum", 60, 0., 6.);
    
    for (Int_t ib = 1; ib < 61; ib++){
      Double_t pt = hPro[j]->GetXaxis()->GetBinCenter(ib);
      Double_t dpt = hPro[j]->GetXaxis()->GetBinWidth(ib);

      hPro[j]->SetBinContent(ib, fProton->Eval(pt));
      //assign conervative error of 10%
      hPro[j]->SetBinError(ib, fProton->Eval(pt)*0.1);
      
      hLambda[j]->SetBinContent(ib, fLambda->Eval(pt));
      //assign conervative error of 10%
      hLambda[j]->SetBinError(ib, fLambda->Eval(pt)*0.1);
    }
    
    hPro[j]->SetLineColor(color[j]);
    hPro[j]->SetMarkerColor(color[j]);
    hPro[j]->SetMarkerStyle(20);

    hLambda[j]->SetLineColor(color[j]);
    hLambda[j]->SetMarkerColor(color[j]);
    hLambda[j]->SetMarkerStyle(25);
  }
  /*
  TCanvas * model = new TCanvas("model","model", 800, 800);
  model->cd();
  hPro[0]->Draw();
  hPro[1]->Draw("same");
  hLambda[0]->Draw("same");
  hLambda[1]->Draw("same");
  */
  if (particle.Contains("proton")) return hPro[selcent];
  return hLambda[selcent];
}
