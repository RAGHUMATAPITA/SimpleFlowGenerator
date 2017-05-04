void doit(const char*);


void plot()
{

  doit("400tracks_10000events");

}


void doit(const char *flag)
{

  TCanvas *c1 = new TCanvas();

  TFile *file = TFile::Open(Form("out_%s.root",flag));

  TH1D *th1d_v2pT_pure = ((TProfile *)file->Get("tp1d_v2pT_pure"))->ProjectionX();
  TH1D *th1d_v2pT_true = ((TProfile *)file->Get("tp1d_v2pT_true"))->ProjectionX();
  TH1D *th1d_v2pT_reco = ((TProfile *)file->Get("tp1d_v2pT_reco"))->ProjectionX();
  TH1D *th1d_YZ_v2pT_reco = ((TProfile *)file->Get("tp1d_YZ_v2pT_reco"))->ProjectionX();
  TH1D *th1d_c2 = ((TProfile *)file->Get("tp1d_c2"))->ProjectionX();
  TH1D *th1d_d2pT = ((TProfile *)file->Get("tp1d_d2pT"))->ProjectionX();
  TH1D *helper = (TH1D *)th1d_d2pT->Clone();
  TH1D *th1d_YZ_d2pT = ((TProfile *)file->Get("tp1d_YZ_d2pT"))->ProjectionX();
  TH1D *helperYZ = (TH1D *)th1d_YZ_d2pT->Clone();
  TH1D *th1d_YZ_W = ((TProfile *)file->Get("tp1d_YZ_W"))->ProjectionX();
  TH1D *th1d_psi2_tmr = ((TProfile *)file->Get("tp1d_psi2_tmr"))->ProjectionX();

  double resolution = th1d_psi2_tmr->GetBinContent(1);

  th1d_v2pT_true->SetMaximum(0.22);
  th1d_v2pT_true->SetMinimum(0.0);
  th1d_v2pT_true->GetYaxis()->SetTitle("v_{2}, various methods");
  th1d_v2pT_true->GetYaxis()->SetTitleOffset(1.2);
  th1d_v2pT_true->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  th1d_v2pT_true->SetLineColor(kBlack);
  th1d_v2pT_true->Draw();
  th1d_v2pT_pure->SetLineColor(kBlack);
  th1d_v2pT_pure->Draw("same");

  th1d_v2pT_reco->Scale(1.0/resolution);
  th1d_v2pT_reco->SetLineColor(kRed);
  th1d_v2pT_reco->Draw("same");
  th1d_YZ_v2pT_reco->Divide(th1d_YZ_W); // weight
  th1d_YZ_v2pT_reco->SetLineColor(kRed);
  //th1d_YZ_v2pT_reco->Draw("same"); //  weight produces identical result

  float c2 = th1d_c2->GetBinContent(1);
  th1d_d2pT->Scale(1.0/sqrt(c2));
  th1d_d2pT->SetLineColor(kBlue);
  th1d_d2pT->Draw("same");
  th1d_YZ_d2pT->Divide(th1d_YZ_W);
  th1d_YZ_d2pT->Scale(1.0/sqrt(c2));
  th1d_YZ_d2pT->SetLineColor(kBlue);
  //th1d_YZ_d2pT->Draw("same");
  helper->Scale(2*c2);
  helperYZ->Scale(2*c2);

  // ---
  TH1D *th1d_p4pT = ((TProfile *)file->Get("tp1d_p4pT"))->ProjectionX();
  TH1D *th1d_p4 = ((TProfile *)file->Get("tp1d_p4"))->ProjectionX();
  th1d_p4pT->Add(helper,-1.0);
  th1d_p4pT->Scale(-1.0);
  float p4 = th1d_p4->GetBinContent(1);
  float c4 = p4 - 2.0*c2*c2;
  if(c4>0){cout << "LOLOMGWTF" << endl; c4 *= -1;}
  float k4 = pow(-c4,0.75);
  th1d_p4pT->Scale(1.0/k4);
  th1d_p4pT->SetLineColor(kGreen+2);
  th1d_p4pT->Draw("same");
  cout << sqrt(c2) << " " << pow(-c4,0.25) << endl;
  // ---

  // ---
  TH1D *th1d_YZ_p4pT = ((TProfile *)file->Get("tp1d_YZ_p4pT"))->ProjectionX();
  TH1D *th1d_YZ_p4 = ((TProfile *)file->Get("tp1d_YZ_p4"))->ProjectionX();
  th1d_YZ_p4pT->Add(helperYZ,-1.0);
  th1d_YZ_p4pT->Scale(-1.0);
  th1d_YZ_p4pT->Divide(th1d_YZ_W); // weight
  float p4YZ = th1d_YZ_p4->GetBinContent(1);
  float c4YZ = p4YZ - 2.0*c2*c2;
  if(c4>0){cout << "LOLOMGWTF" << endl; c4 *= -1;}
  float k4YZ = pow(-c4YZ,0.75);
  th1d_YZ_p4pT->Scale(1.0/k4YZ);// watch...
  th1d_YZ_p4pT->SetLineColor(kGreen+2);
  //th1d_YZ_p4pT->Draw("same");
  cout << sqrt(c2) << " " << pow(-c4YZ,0.25) << endl; // watch...
  // ---

  TLegend *leg = new TLegend(0.68,0.18,0.88,0.38);
  //TLegend *leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1d_v2pT_true,"truth","elp");
  leg->AddEntry(th1d_v2pT_reco,"v_{2}{#Psi_{2}}","elp");
  leg->AddEntry(th1d_d2pT,"v_{2}{2}","elp");
  leg->AddEntry(th1d_p4pT,"v_{2}{4}","elp");
  leg->SetTextSize(0.05);
  leg->Draw();

  c1->Print(Form("Figures/comparison_%s.png",flag));
  c1->Print(Form("Figures/comparison_%s.pdf",flag));

}
