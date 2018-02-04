void GetThermalModelValues() {
  //
  // GSI-Heidelberg values for deuteron-to-proton and 3He-to-proton
  //
  ifstream ifile;
  ifile.open("FromAnton_ratios2pi_T.txt");
  Float_t temp[30], k2pi[30], p2pi[30], lambda2pi[30];
  Float_t ksi2pi[30], omega2pi[30], d2pi[30], he2pi[30];
  Float_t d2p[30], he2p[30];
  int i=0;
  cout << " temperature" << " / " << " d/p " << " / " << " 3He / p " << endl;
  while(ifile) {
    ifile>>temp[i]>>k2pi[i]>>p2pi[i]>>lambda2pi[i]>>ksi2pi[i]>>omega2pi[i]>>d2pi[i]>>he2pi[i];
    if(ifile.eof()) break;
    //
    d2p[i] = d2pi[i]/p2pi[i];
    he2p[i] = he2pi[i]/p2pi[i];
    //
    cout << temp[i] << " / " << d2p[i] << " / " << he2p[i] << endl;
    //
    i++;
  }

}
