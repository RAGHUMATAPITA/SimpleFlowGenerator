void doit(const char*);


void plot_prime()
{

  doit("prime_10");
  doit("prime_25");
  doit("prime_50");
  doit("prime_55");
  doit("prime_60");
  doit("prime_65");
  doit("prime_70");
  doit("prime_71");
  doit("prime_72");
  doit("prime_73");
  doit("prime_74");
  doit("prime_75");
  doit("prime_99");
  doit("prime_100");
  doit("prime_200");

}


void doit(const char *flag)
{

  TCanvas *c1 = new TCanvas();

  TFile *file = TFile::Open(Form("out_%s.root",flag));
  if ( file == NULL )
    {
      cout << "file not found" << endl;
      return;
    }

  TH1D *th1d_v2pT_pure = ((TProfile *)file->Get("tp1d_v2pT_pure"))->ProjectionX();
  TH1D *th1d_v2pT_true = ((TProfile *)file->Get("tp1d_v2pT_true"))->ProjectionX();
  TH1D *th1d_v2pT_reco = ((TProfile *)file->Get("tp1d_v2pT_reco"))->ProjectionX();
  TH1D *th1d_v2pT_prim = ((TProfile *)file->Get("tp1d_v2pT_prim"))->ProjectionX(); // psi2 same, particle offset
  TH1D *th1d_psi2_tmr = ((TProfile *)file->Get("tp1d_psi2_tmr"))->ProjectionX();

  TH1D *th1d_v2pT_primereco = ((TProfile *)file->Get("tp1d_v2pT_primereco"))->ProjectionX(); // psi2 offset, particle same
  TH1D *th1d_psi2prime_tmr = ((TProfile *)file->Get("tp1d_psi2prime_tmr"))->ProjectionX();

  TH1D *th1d_v2pT_doubleprimereco = ((TProfile *)file->Get("tp1d_v2pT_doubleprimereco"))->ProjectionX(); // psi2 offset, particle same

  double resolution = th1d_psi2_tmr->GetBinContent(1);
  double primeresolution = th1d_psi2prime_tmr->GetBinContent(1);

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

  th1d_v2pT_prim->Scale(1.0/resolution);
  th1d_v2pT_prim->SetLineColor(kBlue);
  th1d_v2pT_prim->Draw("same");

  th1d_v2pT_primereco->Scale(1.0/primeresolution);
  th1d_v2pT_primereco->SetLineColor(kGreen+2);
  th1d_v2pT_primereco->Draw("same");

  TH1D* th1d_v2pT_primerecoNEW = (TH1D*)th1d_v2pT_primereco->Clone();
  th1d_v2pT_primerecoNEW->Scale(primeresolution/resolution); // use regular resolution...
  th1d_v2pT_primerecoNEW->SetLineColor(kMagenta+2);
  th1d_v2pT_primerecoNEW->Draw("same");

  th1d_v2pT_doubleprimereco->Scale(1.0/resolution);
  th1d_v2pT_doubleprimereco->SetLineColor(kOrange+1);
  th1d_v2pT_doubleprimereco->Draw("same");

  //TLegend *leg = new TLegend(0.48,0.18,0.88,0.38);
  TLegend *leg = new TLegend(0.18,0.68,0.38,0.88);
  leg->AddEntry(th1d_v2pT_true,"truth","elp");
  leg->AddEntry(th1d_v2pT_reco,"v_{2}{#Psi_{2}}","elp");
  leg->AddEntry(th1d_v2pT_prim,"v_{2}{#Psi_{2}}, coord offset","elp");
  leg->AddEntry(th1d_v2pT_primereco,"v_{2}{#Psi_{2}}, plane offset","elp");
  leg->AddEntry(th1d_v2pT_primerecoNEW,"v_{2}{#Psi_{2}}, plane offset, wrong EPRC","elp");
  leg->AddEntry(th1d_v2pT_doubleprimereco,"v_{2}{#Psi_{2}}, coord and plane offset","elp");
  leg->SetTextSize(0.05);
  leg->SetFillStyle(0);
  leg->Draw();

  c1->Print(Form("Figures/comparison_%s.png",flag));
  c1->Print(Form("Figures/comparison_%s.pdf",flag));

  // ---

  TH1D *th1d_psi2_reco = (TH1D *)file->Get("th1d_psi2_reco");
  TH1D *th1d_psi2prime_reco = (TH1D *)file->Get("th1d_psi2prime_reco");

  th1d_psi2_reco->Rebin(10);
  th1d_psi2prime_reco->Rebin(10);

  th1d_psi2prime_reco->SetLineColor(kRed);
  th1d_psi2_reco->SetLineColor(kBlack);

  th1d_psi2prime_reco->Draw();
  th1d_psi2prime_reco->GetYaxis()->SetTitle("Number of events");
  th1d_psi2prime_reco->GetXaxis()->SetTitle("#Psi_{2}");
  th1d_psi2_reco->Draw("same");

  //TLegend *leg2 = new TLegend(0.4,0.2,0.6,0.4);
  TLegend *leg2 = new TLegend(0.18,0.68,0.38,0.88);
  leg2->AddEntry(th1d_psi2_reco,"Original #Psi_{2}","l");
  leg2->AddEntry(th1d_psi2prime_reco,"Offset #Psi_{2}","l");
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.05);
  leg2->Draw();

  c1->Print(Form("Figures/comparison_eventplane_%s.png",flag));
  c1->Print(Form("Figures/comparison_eventplane_%s.pdf",flag));



}
