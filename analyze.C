// Author: Ron Belmont
// Date: 2010-05-17

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include <sys/time.h>




#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

using namespace std;






Long64_t nevents = 0;
Long64_t ntracks = 0;


const int maxn = 5000;

int d_mult;
float d_psi2;
int d_charge[maxn];
float d_pt[maxn];
float d_phi[maxn];

TBranch *b_mult;
TBranch *b_psi2;
TBranch *b_charge;
TBranch *b_pt;
TBranch *b_phi;

// Main part of program
int main(int argc, char *argv[])
{

  //Char_t inFile[100];
  //Char_t outFile[100];
  char inFile[100];
  char outFile[100];

  if(argc==1)
    {
      cout<<"Now beginning program"<<endl;
      cout<<"Program name is "<<argv[0]<<endl;
      cout<<"Please enter input file name"<<endl;
      cin>>inFile;
      cout<<"Input file is "<<inFile<<endl;
      cout<<"Please enter output file name"<<endl;
      cin>>outFile;
      cout<<"Output file is "<<outFile<<endl;
    }
  else if(argc==3)
    {
      strcpy(inFile,argv[1]);
      strcpy(outFile,argv[2]);
      cout<<"Now beginning program"<<endl;
      cout<<"Program name is "<<argv[0]<<endl;
      cout<<"Input file is "<<inFile<<endl;
      cout<<"Output file is "<<outFile<<endl;
    }
  else
    {
      cout<<"Wrong number of input arguments"<<endl;
      cout<<"This program takes 0 or 2 arguments"<<endl;
      cout<<"With 0 arguments it prompts the user for the file list and output root file"<<endl;
      cout<<"With 2 arguments the first is the file list and the second is the output root file"<<endl;
      return 1;
    }

  // ----------------------------

  struct timeval Time;

  gettimeofday(&Time,0);
  int begintime = Time.tv_sec;
  //cout<<"begintime is "<<begintime<<endl;

  // ----------------------------

  TFile *mData = new TFile(outFile,"recreate"); // declare output file

  TH1D *th1d_mult = new TH1D("th1d_mult","",5000,0,5000);
  TH1D *th1d_psi2 = new TH1D("th1d_psi2","",630,-3.2,3.2);
  TH1D *th1d_pt = new TH1D("th1d_pt","",200,0,2);
  TH1D *th1d_phi = new TH1D("th1d_phi","",630,-3.2,3.2);

  TProfile *tp1d_v2pT_true = new TProfile("tp1d_v2pT_true","",200,0,2,-1e10,1e10,"");
  TProfile *tp1d_v2pT_reco = new TProfile("tp1d_v2pT_reco","",200,0,2,-1e10,1e10,"");

  TProfile *tp1d_d2pT = new TProfile("tp1d_d2pT","",200,0,2,-1e10,1e10,"");

  TProfile *tp1d_c2 = new TProfile("tp1d_c2","",1,0,1000,-1e10,1e10,"");
  TProfile *tp1d_c2mult = new TProfile("tp1d_c2mult","",1000,0,1000,-1e10,1e10,"");


  TH1D *th1d_psi2_reco = new TH1D("th1d_psi2_reco","",630,-3.2,3.2);
  TH1D *th1d_psi2_tmr = new TH1D("th1d_psi2_tmr","",630,-3.2,3.2);


  // --- Done with Histograms ---------------------

  //Now read in the pDSTs listed in the input files
  TFile *file = TFile::Open(inFile);
  if(!file)
    {
      cout<<"file input error: file does not exist "<<endl;
      return -2;
    }

  TTree *tree = (TTree *)file->Get("simpletree");
  if(!tree)
    {
      cout<<"file input error: cannot find tree "<<endl;
      return -3;
    }

  b_mult = tree->GetBranch("mult");
  b_mult->SetAddress(&d_mult);
  b_psi2 = tree->GetBranch("psi2");
  b_psi2->SetAddress(&d_psi2);

  b_charge = tree->GetBranch("charge");
  b_charge->SetAddress(d_charge);
  b_phi = tree->GetBranch("phi");
  b_phi->SetAddress(d_phi);
  b_pt = tree->GetBranch("pt");
  b_pt->SetAddress(d_pt);


  int nevt = (int)tree->GetEntries(); // number of events in tree
  for(int ievt=0; ievt<nevt; ievt++) // loop over events
    {

      if(ievt%1000 == 0) cout << "processing event number " << ievt << endl;
      //if(ievt>1000) break;

      b_mult->GetEntry(ievt);
      b_psi2->GetEntry(ievt);
      b_charge->GetEntry(ievt);
      b_phi->GetEntry(ievt);
      b_pt->GetEntry(ievt);

      int mult = d_mult;
      th1d_mult->Fill(mult);

      float psi2true = d_psi2;
      th1d_psi2->Fill(psi2true);

      // --- first track loop, q-vectors
      float Q2x = 0;
      float Q2y = 0;
      for(int itrk=0; itrk<mult; itrk++)
	{
	  float phi = d_phi[itrk];
	  Q2x += cos(2*phi);
	  Q2y += sin(2*phi);
	} // End of track loop
      float psi2reco = atan2(Q2y,Q2x);
      th1d_psi2_reco->Fill(psi2reco);
      th1d_psi2_tmr->Fill(psi2true-psi2reco);
      float two = ( Q2x*Q2x + Q2y*Q2y ) / (mult*mult - mult);

      tp1d_c2->Fill(1,two);
      tp1d_c2mult->Fill(mult,two);

      for(int itrk=0; itrk<mult; itrk++)
	{
	  int charge = d_charge[itrk];
	  float phi = d_phi[itrk];
	  float pt = d_pt[itrk];
	  th1d_phi->Fill(phi);
	  th1d_pt->Fill(pt);
	  float v2_track_true = cos(2*phi-psi2true);
	  tp1d_v2pT_true->Fill(pt,v2_track_true);
	  float v2_track_reco = cos(2*phi-psi2reco);
	  tp1d_v2pT_reco->Fill(pt,v2_track_reco);

	  float u2x = cos(2*phi);
	  float u2y = sin(2*phi);

	  float twoprime = ( u2x*Q2x + u2y*Q2y ) / (mult - 1);
	  tp1d_d2pT->Fill(pt,twoprime);

	  ntracks++; // count total number of tracks
	} // End of track loop

      nevents++; // count total number of events

    } // End of event loop

  delete tree;
  file->Close();
  delete file;



  mData->Write();
  mData->Close();

  cout<<"Number of events: "<<nevents<<endl;
  cout<<"Number of tracks: "<<ntracks<<endl;

  gettimeofday(&Time,0);
  int endtime = Time.tv_sec;
  //cout<<"endtime is "<<endtime<<endl;

  int tdiff = endtime-begintime;

  cout<<"End of program."<<endl;
  cout<<"Execution time: "<<tdiff<<" seconds"<<endl;

  return 0;

}
