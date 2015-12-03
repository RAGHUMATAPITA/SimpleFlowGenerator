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

float calc4_event(float, float, float, float, float);
float calc4_track(float, float, float, float, float, float, float, float, float);

float calc4_event_YZ(float, float, float, float, float);
float calc4_track_YZ(float, float, float, float, float, float, float, float, float);



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

  TProfile *tp1d_v2pT_true = new TProfile("tp1d_v2pT_true","",20,0,2,-1e10,1e10,"");
  TProfile *tp1d_v2pT_reco = new TProfile("tp1d_v2pT_reco","",20,0,2,-1e10,1e10,"");

  TProfile *tp1d_d2pT = new TProfile("tp1d_d2pT","",20,0,2,-1e10,1e10,"");
  TProfile *tp1d_c2 = new TProfile("tp1d_c2","",1,0,1000,-1e10,1e10,"");
  TProfile *tp1d_c2mult = new TProfile("tp1d_c2mult","",1000,0,1000,-1e10,1e10,"");

  TProfile *tp1d_p4pT = new TProfile("tp1d_p4pT","",20,0,2,-1e10,1e10,"");
  TProfile *tp1d_p4 = new TProfile("tp1d_p4","",1,0,1000,-1e10,1e10,"");
  TProfile *tp1d_p4mult = new TProfile("tp1d_p4mult","",1000,0,1000,-1e10,1e10,"");

  TProfile *tp1d_YZ_p4pT = new TProfile("tp1d_YZ_p4pT","",20,0,2,-1e10,1e10,"");
  TProfile *tp1d_YZ_p4 = new TProfile("tp1d_YZ_p4","",1,0,1000,-1e10,1e10,"");
  TProfile *tp1d_YZ_p4mult = new TProfile("tp1d_YZ_p4mult","",1000,0,1000,-1e10,1e10,"");


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
      float Q4x = 0;
      float Q4y = 0;
      for(int itrk=0; itrk<mult; itrk++)
	{
	  float phi = d_phi[itrk];
	  Q2x += cos(2*phi);
	  Q2y += sin(2*phi);
	  Q4x += cos(4*phi);
	  Q4y += sin(4*phi);
	} // End of track loop
      float psi2reco = atan2(Q2y,Q2x);
      th1d_psi2_reco->Fill(psi2reco);
      th1d_psi2_tmr->Fill(psi2true-psi2reco);
      float two = ( Q2x*Q2x + Q2y*Q2y ) / (mult*mult - mult);
      float four = calc4_event(Q2x,Q2y,Q4x,Q4y,mult);
      float fourYZ = calc4_event_YZ(Q2x,Q2y,Q4x,Q4y,mult);

      tp1d_c2->Fill(1,two);
      tp1d_c2mult->Fill(mult,two);
      tp1d_p4->Fill(1,four);
      tp1d_p4mult->Fill(mult,four);
      tp1d_YZ_p4->Fill(1,fourYZ);
      tp1d_YZ_p4mult->Fill(mult,fourYZ);

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

	  float u4x = cos(4*phi);
	  float u4y = sin(4*phi);

	  float fourprime = calc4_track(u2x,u2y,u4x,u4y,Q2x,Q2y,Q4x,Q4y,mult);
	  tp1d_p4pT->Fill(pt,fourprime);

	  float fourprimeYZ = calc4_track_YZ(u2x,u2y,u4x,u4y,Q2x,Q2y,Q4x,Q4y,mult);
	  tp1d_YZ_p4pT->Fill(pt,fourprimeYZ);

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



float calc4_event(float Xn, float Yn, float X2n, float Y2n, float M)
{

  float Qn2 = Xn*Xn+Yn*Yn;
  float Q2n2 = X2n*X2n+Y2n*Y2n;
  float Qn4 = Qn2*Qn2;
  float Qn2d = Xn*Xn-Yn*Yn;

  float first = Qn4+Q2n2-(2*(X2n*Qn2d));
  float second = 2*(2*(M-2)*Qn2)-(M*(M-3));

  float W_4 = M*(M-1)*(M-2)*(M-3);

  return (first-second)/W_4;

}



float calc4_track(float xn, float yn, float x2n, float y2n, float Xn, float Yn, float X2n, float Y2n, float M)
{

  // --- this code based on a simplified version of the analytical expression
  // --- this code obviously has enormous room for improvement and cleanup, which is welcomed
  // --- also it is not clear if this code works correctly, but I think it doesn't
  float one   = (xn*Xn + yn*Yn)*(Xn*Xn + Yn*Yn);
  float two   = x2n*Xn*Xn - x2n*Yn*Yn + 2*y2n*Xn*Yn;
  float three = xn*Xn*X2n + xn*Yn*Y2n - yn*(X2n*Yn - Xn*Y2n);
  float four  = 2*M*(xn*Xn + yn*Yn);
  float five  = 2*(Xn*Xn + Yn*Yn);
  float six   = 7*(xn*Xn + yn*Yn);
  float seven = xn*Xn + yn*Yn;
  float eight = x2n*X2n + y2n*Y2n;
  float nine = 2*(xn*Xn + yn*Yn);

  float numerator = one - two - three - four - five + six - seven + eight + nine +2*M - 6;
  float denominator = (M-1)*(M-2)*(M-3);

  return numerator/denominator;

}


float calc4_event_YZ(float QTx, float QTy, float QT2x, float QT2y, float MQT)
{

  float cn4 = ((QTx*QTx + QTy*QTy)*(QTx*QTx + QTy*QTy) + (QT2x*QT2x + QT2y*QT2y) -2*(QT2x*QTx*QTx + 2*QT2y*QTx*QTy - QT2x*QTy*QTy) -2*(2*(MQT-2)*(QTx*QTx + QTy*QTy) - MQT*(MQT-3)) )/ (MQT*(MQT-1)*(MQT-2)*(MQT-3));

  return cn4;

}


float calc4_track_YZ(float pnx, float pny, float p2nx, float p2ny, float QTx, float QTy, float QT2x, float QT2y, float MQT)
{

  int mp = 1;

  float dn4 = (((pnx*QTx + pny*QTy)*(QTx*QTx + QTy*QTy)) -
	       (p2nx*QTx*QTx - p2nx*QTy*QTy + 2*p2ny*QTx*QTy)  -
	       (pnx*QTx*QT2x - pny*QTy*QT2x + pnx*QTy*QT2y + pny*QTx*QT2y) -
	       2*MQT*(pnx*QTx + pny*QTy) -
	       2*mp*(QTx*QTx + QTy*QTy) +
	       7*(pnx*QTx + pny*QTy) -
	       (QTx*pnx + QTy*pny) +
	       (p2nx*QT2x + p2ny*QT2y) +
	       2*(pnx*QTx + pny*QTy) + 2*mp*MQT - 6*mp ) / (mp*(MQT-1)*(MQT-2)*(MQT-3));

  return dn4;

}

