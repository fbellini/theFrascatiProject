void CompareHeinzToSato() {
  //
  // this small macro compares the coalescence models from Heinz and Sato
  // Phys.Lett. 98B (1981) 153-157
  //
  // TF1 * funcB2Sato = new TF1("funcB2Sato","(2*0.938/(0.938*0.938)) * 6.0 * TMath::Power(2*TMath::Pi(),1.5)*TMath::Power(x/0.197*x/0.197 + 3.2/0.197*3.2/0.197/2.0, -1.5)", 0.0, 5.0);
  //  TF1 * funcB2Sato = new TF1("funcB2Sato","(2./0.938) * 6./8.  * TMath::Power(2*TMath::Pi(),1.5)*TMath::Power((TMath::Sqrt(2./3.)*x/0.197)*(TMath::Sqrt(2./3.)*x/0.197) + (3.2/0.197)*(3.2/0.197)/2.0, -1.5)", 0.0, 5.0);
  TF1 * funcB2Sato = new TF1("funcB2Sato","(2./0.938) * 6./(2*2*2.)  * TMath::Power(2*TMath::Pi(),1.5)*TMath::Power((x/0.197)*(x/0.197) + (3.2/0.197)*(3.2/0.197)/2.0, -1.5)", 0.0, 5.0);

  TF1 * funcB2Sato2 = new TF1("funcB2Sato2","(1./0.938) * (3./4)  * TMath::Power(2,1.5)*TMath::Power(4*TMath::Pi() * (0.2/(x*x)/(0.2 + 1./(x*x))), 1.5)*0.197*0.197*0.197*0.197", 0.0, 5.0);
  TF1 * funcB2Heinz = new TF1("funcB2Heinz", " 3./(4*TMath::Sqrt(2)) * TMath::Power(2*TMath::Pi(),1.5) * 1/[0] * TMath::Power(x/0.197*x/0.197 + 3.2/0.197*3.2/0.197/4.0, -1.5)", 0.0, 5.0);
  funcB2Heinz->SetParName(0, "mT");
  //funcB2Heinz->SetParameter(0, TMath::Sqrt(0.73*0.73 + 0.938*0.938));
  funcB2Heinz->SetParameter(0, TMath::Sqrt(0.938*0.938));
  //
  funcB2Heinz->Draw();
  funcB2Sato->SetLineColor(kBlue);
  funcB2Sato->Draw("SAME");

}
