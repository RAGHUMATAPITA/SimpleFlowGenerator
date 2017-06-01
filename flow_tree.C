#include <iostream>
#include <cstdlib>
#include <vector>


#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TH1.h"
#include "TF1.h"


using namespace std;


int dostuff(const int, TTree*); // function prototype...


const double pi = TMath::Pi();

const double v2prefactor = 0.05;

const int maxmult = 400;
const int nptbins = 100;

int d_mult;
float d_psi2;
int d_charge[maxmult];
float d_v2[maxmult];
float d_pt[maxmult];
float d_phi[maxmult];

TF1 *funpt;
TF1 *funphi;
TF1 *funphiarray[nptbins];


int main()
{

  cout << "setting up functions" << endl;

  funpt = new TF1("funpt","TMath::Exp(-0.35*x)",0.0,2.0); // inverse slope param of 350 MeV
  funphi = new TF1("funphi","1 + 2*[0]*TMath::Cos(2*x)",-pi,pi);
  funphi->SetParameter(0,0.2);

  for(int i=0; i<nptbins; i++)
    {
      float pt = 2.0*i/float(nptbins);
      float v2 = v2prefactor*pt;
      funphiarray[i] = new TF1("funphi","1 + 2*[0]*TMath::Cos(2*x)",-pi,pi);
      funphiarray[i]->SetParameter(0,v2);
    }


  cout << "setting up tree " << endl;


  TFile *file = new TFile("simpletree.root","recreate");

  TH1D *hmult = new TH1D("hmult","",200,0,200);

  TTree *tree = new TTree("simpletree","simpletree");
  tree->Branch("mult",&d_mult,"mult/I");
  tree->Branch("psi2",&d_psi2,"psi2/F");
  tree->Branch("charge",d_charge,"charge[mult]/I");
  tree->Branch("v2",d_v2,"v2[mult]/F");
  tree->Branch("pt",d_pt,"pt[mult]/F");
  tree->Branch("phi",d_phi,"phi[mult]/F");

  cout << "now starting loop" << endl;

  int nevents = 1000000; // way too slow for v2(pT), okay but long for fixed v2
  //int nevents = 100000; // way too slow for v2(pT), fine for fixed v2
  //int nevents = 10000; // good for either
  for(int i=0; i<nevents; i++)
    {
      if(i % 1000 == 0) cout << i << " events processed so far " << endl;
      int ntracks = maxmult; //...
      ntracks = 30;
      int ntrk = dostuff(ntracks,tree);
      hmult->Fill(ntrk);
    }

  file->Write();
  file->Close();
  // delete hmult;
  // delete file;
  // delete tree;

  return 0;

}




int dostuff(const int number, TTree *tree)
{

  float psi2 = gRandom->Uniform(-pi,pi); // phi range for throw of psi2 for event

  for(int i=0; i<number; i++)
    {
      int charge = -1;
      double random = gRandom->Rndm();
      if(random>0.5) charge = 1;

      float pt = funpt->GetRandom();

      float v2 = v2prefactor*pt;
      //funphi->SetParameter(0,v2);
      int ptbin = pt*nptbins/2;
      float phi = funphiarray[ptbin]->GetRandom();
      phi += psi2;
      if(phi>pi) phi -= 2*pi;
      if(phi<-pi) phi += 2*pi;
      // ---------------------------------------------------
      //cout << pt << " " << ptbin << endl;
      //cout << funphiarray[ptbin]->GetParameter(0) << endl;
      // ---
      d_charge[i] = charge;
      d_v2[i] = v2;
      d_pt[i] = pt;
      d_phi[i] = phi;
    }

  d_psi2 = psi2;
  d_mult = number;

  tree->Fill();

  return d_mult;

}
