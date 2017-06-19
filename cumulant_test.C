void doit(const char*);

void cumulant_test()
{

  doit("testest10_new05.root");
  doit("testest15_new05.root");
  doit("testest20_new05.root");
  doit("testest25_new05.root");
  doit("testest30_new05.root");
  doit("testest400_new05.root");

}


void doit(const char* filename)
{

  cout << filename << endl;
  TFile* file = TFile::Open(filename);

  TProfile* tp1d_two = (TProfile*)file->Get("tp1d_rescursion_0_0");
  TProfile* tp1d_for = (TProfile*)file->Get("tp1d_rescursion_0_2");
  TProfile* tp1d_six = (TProfile*)file->Get("tp1d_rescursion_0_4");
  TProfile* tp1d_eit = (TProfile*)file->Get("tp1d_rescursion_0_6");
  //pow(tp1d_rescursion_0_2->GetBinContent(1),(1.0/4.0))

  double dtwo = tp1d_two->GetBinContent(1);
  double dfor = tp1d_for->GetBinContent(1);
  double dsix = tp1d_six->GetBinContent(1);
  double deit = tp1d_eit->GetBinContent(1);

  double ctwo = dtwo;
  double cfor = dfor - (2*dtwo*dtwo);
  double csix = dsix - (9*dfor*dtwo) + (12*dtwo*dtwo*dtwo);
  double ceit = deit - (16*dsix*dtwo) - (18*dfor*dfor) + (144*dfor*dtwo*dtwo) - (144*dtwo*dtwo*dtwo*dtwo);

  double ftwo = pow(dtwo,(1.0/2.0));
  double ffor = pow(dfor,(1.0/4.0));
  double fsix = pow(dsix,(1.0/6.0));
  double feit = pow(deit,(1.0/8.0));

  double vtwo = pow(ctwo,(1.0/2.0));
  double vfor = pow(-cfor,(1.0/4.0));
  double vsix = pow(0.25*csix,(1.0/6.0));
  double veit = pow(-(1.0/33.0)*ceit,(1.0/8.0));

  cout << ftwo << " " << vtwo << endl;
  cout << ffor << " " << vfor << endl;
  cout << fsix << " " << vsix << endl;
  cout << feit << " " << veit << endl;


}

