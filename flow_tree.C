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


double width = 0.5;

int dostuff(const int, TTree*); // function prototype...


const int maxmult = 5000;

int d_mult;
float d_psi2;
int d_charge[maxmult];
float d_pt[maxmult];
float d_phi[maxmult];



int main()
{

  TFile *file = new TFile("simpletree.root","recreate");

  TH1D *hmult = new TH1D("hmult","",200,0,200);

  TTree *tree = new TTree("simpletree","simpletree");
  tree->Branch("mult",&d_mult,"mult/I");
  tree->Branch("psi2",&d_psi2,"psi2/F");
  tree->Branch("charge",d_charge,"charge[mult]/I");
  tree->Branch("pt",d_pt,"pt[mult]/F");
  tree->Branch("phi",d_phi,"phi[mult]/F");

  //int nevents = 100000; // way too slow...
  int nevents = 10000;
  for(int i=0; i<nevents; i++)
    {
      int ntracks = 100;
      int ntrk = dostuff(ntracks,tree);
      hmult->Fill(ntrk);
      if(i % 1000 == 0) cout << i << " events processed so far " << endl;
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

  vector<int> vec_charge;
  vector<double> vec_phi;
  TF1 *funpt = new TF1("funpt","TMath::Exp(-0.35*x)",0.0,2.0); // inverse slope param of 350 MeV
  TF1 *funphi = new TF1("funphi","1 + 2*[0]*TMath::Cos(2*x - [1])",-TMath::Pi(),TMath::Pi());
  float psi2 = gRandom->Uniform(-TMath::Pi(),TMath::Pi()); // phi range for throw of psi2 for event
  funphi->SetParameter(1,psi2);
  for(int i=0; i<number; i++)
    {
      int charge = -1;
      double random = gRandom->Rndm();
      if(random>0.5) charge = 1;
      float pt = funpt->GetRandom(); // slow...
      float v2 = 0.2*pt; // v2 = 0.2 at 1 GeV
      funphi->SetParameter(0,v2);
      float phi = funphi->GetRandom(); // very slow...
      // ---
      d_charge[i] = charge;
      d_pt[i] = pt;
      d_phi[i] = phi;
    }

  d_psi2 = psi2;
  d_mult = number;

  tree->Fill();

  return d_mult;

}
