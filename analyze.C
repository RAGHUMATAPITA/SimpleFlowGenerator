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
#include "TMath.h"
#include "TComplex.h"

using namespace std;



Long64_t nevents = 0;
Long64_t ntracks = 0;


const float twopi = 2*TMath::Pi();


const int maxn = 5000;

int d_mult;
float d_psi2;
int d_charge[maxn];
float d_pt[maxn];
float d_v2[maxn];
float d_phi[maxn];

TBranch *b_mult;
TBranch *b_psi2;
TBranch *b_charge;
TBranch *b_v2;
TBranch *b_pt;
TBranch *b_phi;

double calc4_event(double, double, double, double, double);
double calc4_track(double, double, double, double, double, double, double, double, double);

double calc4_event_YZ(double, double, double, double, double);
double calc4_track_YZ(double, double, double, double, double, double, double, double, double);

double calc6_event(TComplex&, TComplex&, TComplex&, float);

// --- from generic forumulas ----------------------------------------------------
TComplex Q(int, int); // forward declaration...
TComplex Recursion(int, int*); // forward declaration...
TComplex Recursion(int, int*, int, int); // forward declaration...
//const int h1=1, h2=3, h3=5, h4=0, h5=-2, h6=-4, h7=-1, h8=-6; // from the code
const int h1=2, h2=-2, h3=2, h4=-2, h5=2, h6=-2, h7=2, h8=-2; // simplifed to v2
const int sum = (h1<0?-1*h1:h1)+(h2<0?-1*h2:h2)+(h3<0?-1*h3:h3)+(h4<0?-1*h4:h4)
                + (h5<0?-1*h5:h5)+(h6<0?-1*h6:h6)+(h7<0?-1*h7:h7)+(h8<0?-1*h8:h8);
const int maxCorrelator = 8; // We will not go beyond 8-p correlations
const int maxHarmonic = sum+1;
const int maxPower = maxCorrelator+1;
TComplex Qvector[maxHarmonic][maxPower]; // All needed Q-vector components
// -------------------------------------------------------------------------------


// Main part of program
int main(int argc, char *argv[])
{

  //Char_t inFile[100];
  //Char_t outFile[100];
  char inFile[100];
  char outFile[100];

  char c_offset[100];
  double offset = 0;

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
  else if(argc==4)
    {
      strcpy(inFile,argv[1]);
      strcpy(outFile,argv[2]);
      strcpy(c_offset,argv[3]);
      cout<<"Now beginning program"<<endl;
      cout<<"Program name is "<<argv[0]<<endl;
      cout<<"Input file is "<<inFile<<endl;
      cout<<"Output file is "<<outFile<<endl;
      string s_offset(c_offset);
      offset = stof(s_offset);
      cout << "Numerical offset is " << offset << endl;
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

  TProfile *tp1d_v2pT_pure = new TProfile("tp1d_v2pT_pure","",20,0,2,-1e10,1e10,"");
  TProfile *tp1d_v2pT_true = new TProfile("tp1d_v2pT_true","",20,0,2,-1e10,1e10,"");
  TProfile *tp1d_v2pT_reco = new TProfile("tp1d_v2pT_reco","",20,0,2,-1e10,1e10,"");
  TProfile *tp1d_v2pT_prim = new TProfile("tp1d_v2pT_prim","",20,0,2,-1e10,1e10,"");
  TProfile *tp1d_v2pT_primereco = new TProfile("tp1d_v2pT_primereco","",20,0,2,-1e10,1e10,"");
  TProfile *tp1d_v2pT_doubleprimereco = new TProfile("tp1d_v2pT_doubleprimereco","",20,0,2,-1e10,1e10,"");

  TProfile *tp1d_d2pT = new TProfile("tp1d_d2pT","",20,0,2,-1e10,1e10,"");
  TProfile *tp1d_c2 = new TProfile("tp1d_c2","",1,0,1000,-1e10,1e10,"");
  TProfile *tp1d_c2mult = new TProfile("tp1d_c2mult","",1000,0,1000,-1e10,1e10,"");

  TProfile *tp1d_p4pT = new TProfile("tp1d_p4pT","",20,0,2,-1e10,1e10,"");
  TProfile *tp1d_p4 = new TProfile("tp1d_p4","",1,0,1000,-1e10,1e10,"");
  TProfile *tp1d_p4mult = new TProfile("tp1d_p4mult","",1000,0,1000,-1e10,1e10,"");

  TProfile *tp1d_YZ_p4pT = new TProfile("tp1d_YZ_p4pT","",20,0,2,-1e10,1e10,"");
  TProfile *tp1d_YZ_W = new TProfile("tp1d_YZ_W","",20,0,2,-1e10,1e10,"");
  TProfile *tp1d_YZ_p4 = new TProfile("tp1d_YZ_p4","",1,0,1000,-1e10,1e10,"");
  TProfile *tp1d_YZ_p4mult = new TProfile("tp1d_YZ_p4mult","",1000,0,1000,-1e10,1e10,"");
  TProfile *tp1d_YZ_d2pT = new TProfile("tp1d_YZ_d2pT","",20,0,2,-1e10,1e10,"");
  TProfile *tp1d_YZ_v2pT_true = new TProfile("tp1d_YZ_v2pT_true","",20,0,2,-1e10,1e10,"");
  TProfile *tp1d_YZ_v2pT_reco = new TProfile("tp1d_YZ_v2pT_reco","",20,0,2,-1e10,1e10,"");


  TH1D *th1d_psi2_reco = new TH1D("th1d_psi2_reco","",100,-twopi,twopi);
  TH1D *th1d_psi2_tmr = new TH1D("th1d_psi2_tmr","",100,-twopi,twopi);
  TProfile *tp1d_psi2_tmr = new TProfile("tp1d_psi2_tmr","",1,-twopi,twopi);

  TH1D *th1d_psi2prime_reco = new TH1D("th1d_psi2prime_reco","",100,-twopi,twopi);
  TH1D *th1d_psi2prime_tmr = new TH1D("th1d_psi2prime_tmr","",100,-twopi,twopi);
  TProfile *tp1d_psi2prime_tmr = new TProfile("tp1d_psi2prime_tmr","",1,-twopi,twopi);

  // --- from generic formulas
  TProfile *recursion[2][maxCorrelator] = {{NULL}}; // Correlations calculated from Q-vector components using recursive algorithm
  for(int cs=0;cs<2;cs++)
    {
      for(int c=0;c<maxCorrelator;c++)
        {
          // correlations[cs][c] = new TProfile("","",1,0.,1.);
          // correlations[cs][c]->Sumw2();
          recursion[cs][c] = new TProfile(Form("tp1d_rescursion_%d_%d",cs,c),"",1,0.,1.);
          recursion[cs][c]->Sumw2();
          // nestedLoops[cs][c] = new TProfile("","",1,0.,1.);
          // nestedLoops[cs][c]->Sumw2();
        } // end of for(int c=0;c<maxCorrelator;c++)
    } // end of for(int cs=0;cs<2;cs++)

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
  b_v2 = tree->GetBranch("v2");
  b_v2->SetAddress(d_v2);


  int nevt = (int)tree->GetEntries(); // number of events in tree
  for(int ievt=0; ievt<nevt; ievt++) // loop over events
    {

      bool say_event = ( ievt%1000 == 0) ;
      if ( say_event ) cout << "processing event number " << ievt << endl;
      if ( ievt > 500000 ) break;

      b_mult->GetEntry(ievt);
      b_psi2->GetEntry(ievt);
      b_charge->GetEntry(ievt);
      b_phi->GetEntry(ievt);
      b_pt->GetEntry(ievt);

      int mult = d_mult;
      th1d_mult->Fill(mult);

      double psi2true = d_psi2;
      th1d_psi2->Fill(psi2true);

      // --- for the generic formulas ---------
      for(int h=0;h<maxHarmonic;h++)
        {
          for(int p=0;p<maxPower;p++)
            {
              Qvector[h][p] = TComplex(0.,0.);
            } //  for(int p=0;p<maxPower;p++)
        } // for(int h=0;h<maxHarmonic;h++)
      // --------------------------------------
      // --- first track loop, q-vectors
      double Q2x = 0;
      double Q2y = 0;
      double Q4x = 0;
      double Q4y = 0;
      double Q6x = 0;
      double Q6y = 0;
      double Q2x_sub = 0;
      double Q2y_sub = 0;
      double Q4x_sub = 0;
      double Q4y_sub = 0;
      double Q2xprime = 0;
      double Q2yprime = 0;
      double Q4xprime = 0;
      double Q4yprime = 0;
      double Q2x_subprime = 0;
      double Q2y_subprime = 0;
      double Q4x_subprime = 0;
      double Q4y_subprime = 0;
      for(int itrk=0; itrk<mult; itrk++)
	{
	  double phi = d_phi[itrk];
	  double x  = cos(phi);
	  double y  = sin(phi);
	  double xoff = offset;
	  double yoff = offset;
	  double xprime = x - xoff;
	  double yprime = y - yoff;
	  double phiprime = atan2(yprime,xprime);
          // --- from generic formulas ----------------------------------------------------------------------
          double dPhi = 0.0; // particle angle
          double wPhi = 1.0; // particle weight
          double wPhiToPowerP = 1.0; // particle weight raised to power p
          dPhi = phi; // minimal change from me to match the generic forumlas code
          for(int h=0;h<maxHarmonic;h++)
            {
              for(int p=0;p<maxPower;p++)
                {
                  //if(bUseWeights){wPhiToPowerP = pow(wPhi,p);} // no weights for us...
                  Qvector[h][p] += TComplex(wPhiToPowerP*TMath::Cos(h*dPhi),wPhiToPowerP*TMath::Sin(h*dPhi));
                } //  for(int p=0;p<maxPower;p++)
            } // for(int h=0;h<maxHarmonic;h++)
          // ------------------------------------------------------------------------------------------------
	  Q2x += cos(2*phi);
	  Q2y += sin(2*phi);
	  Q4x += cos(4*phi);
	  Q4y += sin(4*phi);
	  Q6x += cos(6*phi);
	  Q6y += sin(6*phi);
	  Q2xprime += cos(2*phiprime);
	  Q2yprime += sin(2*phiprime);
	  Q4xprime += cos(4*phiprime);
	  Q4yprime += sin(4*phiprime);
	  if(itrk>mult/2) continue;
	  Q2x_sub += cos(2*phi);
	  Q2y_sub += sin(2*phi);
	  Q4x_sub += cos(4*phi);
	  Q4y_sub += sin(4*phi);
	  Q2x_subprime += cos(2*phiprime);
	  Q2y_subprime += sin(2*phiprime);
	  Q4x_subprime += cos(4*phiprime);
	  Q4y_subprime += sin(4*phiprime);
	} // End of track loop
      //double psi2reco = atan2(Q2y,Q2x);
      double psi2reco = atan2(Q2y_sub,Q2x_sub)/2;
      th1d_psi2_reco->Fill(psi2reco);
      th1d_psi2_tmr->Fill(cos(2*(psi2true-psi2reco)));
      tp1d_psi2_tmr->Fill(1,cos(2*(psi2true-psi2reco)));
      double psi2primereco = atan2(Q2y_subprime,Q2x_subprime)/2;
      th1d_psi2prime_reco->Fill(psi2primereco);
      th1d_psi2prime_tmr->Fill(cos(2*(psi2true-psi2primereco)));
      tp1d_psi2prime_tmr->Fill(1,cos(2*(psi2true-psi2primereco)));
      double two = ( Q2x*Q2x + Q2y*Q2y - mult) / (mult*mult - mult);
      double four = calc4_event(Q2x,Q2y,Q4x,Q4y,mult);
      double fourYZ = calc4_event_YZ(Q2x,Q2y,Q4x,Q4y,mult);

      tp1d_c2->Fill(1,two);
      tp1d_c2mult->Fill(mult,two);
      tp1d_p4->Fill(1,four);
      tp1d_p4mult->Fill(mult,four);
      tp1d_YZ_p4->Fill(1,fourYZ);
      tp1d_YZ_p4mult->Fill(mult,fourYZ);

      TComplex tc_Q2(Q2x,Q2y);
      TComplex tc_Q4(Q4x,Q4y);
      TComplex tc_Q6(Q6x,Q6y);
      double six = calc6_event(tc_Q2,tc_Q4,tc_Q6,mult);

      // --- from generic formulas ----------------------------------------------------------------------------
      //  2-p correlations:
      //cout<<" => Calculating 2-p correlations (using recursion)...       \r"<<flush;
      int harmonics_Two_Num[2] = {h1,h2}; // 2, -2
      int harmonics_Two_Den[2] = {0,0}; // recursion gives right combinatorics
      TComplex twoRecursion = Recursion(2,harmonics_Two_Num)/Recursion(2,harmonics_Two_Den).Re();
      double wTwoRecursion = Recursion(2,harmonics_Two_Den).Re();
      recursion[0][0]->Fill(0.5,twoRecursion.Re(),wTwoRecursion); // <<cos(h1*phi1+h2*phi2)>>
      recursion[1][0]->Fill(0.5,twoRecursion.Im(),wTwoRecursion); // <<sin(h1*phi1+h2*phi2)>>
      //  4-p correlations:
      //cout<<" => Calculating 4-p correlations (using recursion)...       \r"<<flush;
      int harmonics_Four_Num[4] = {h1,h2,h3,h4}; // 2, -2, 2, -2 // what about 2, 2, -2, -2?
      int harmonics_Four_Den[4] = {0,0,0,0}; // recursion gives right combinatorics
      TComplex fourRecursion = Recursion(4,harmonics_Four_Num)/Recursion(4,harmonics_Four_Den).Re();
      double wFourRecursion = Recursion(4,harmonics_Four_Den).Re();
      recursion[0][2]->Fill(0.5,fourRecursion.Re(),wFourRecursion); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4)>>
      recursion[1][2]->Fill(0.5,fourRecursion.Im(),wFourRecursion); // <<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4)>>
      //  6-p correlations:
      //cout<<" => Calculating 6-p correlations (using recursion)...       \r"<<flush;
      int harmonics_Six_Num[6] = {h1,h2,h3,h4,h5,h6};
      int harmonics_Six_Den[6] = {0,0,0,0,0,0};
      TComplex sixRecursion = Recursion(6,harmonics_Six_Num)/Recursion(6,harmonics_Six_Den).Re();
      double wSixRecursion = Recursion(6,harmonics_Six_Den).Re();
      recursion[0][4]->Fill(0.5,sixRecursion.Re(),wSixRecursion); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>>
      recursion[1][4]->Fill(0.5,sixRecursion.Im(),wSixRecursion); // <<<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>>
      //  8-p correlations:
      //cout<<" => Calculating 6-p correlations (using recursion)...       \r"<<flush;
      int harmonics_Eight_Num[8] = {h1,h2,h3,h4,h5,h6,h7,h8};
      int harmonics_Eight_Den[8] = {0,0,0,0,0,0,0,0};
      TComplex eightRecursion = Recursion(8,harmonics_Eight_Num)/Recursion(8,harmonics_Eight_Den).Re();
      double wEightRecursion = Recursion(8,harmonics_Eight_Den).Re();
      recursion[0][6]->Fill(0.5,eightRecursion.Re(),wEightRecursion); // <<cos(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>>
      recursion[1][6]->Fill(0.5,eightRecursion.Im(),wEightRecursion); // <<<sin(h1*phi1+h2*phi2+h3*phi3+h4*phi4+h5*phi5+h6*phi6)>>
      // ------------------------------------------------------------------------------------------------------

      if ( say_event )
        {
          cout << "conventional two particle calculation: " << two << endl;
          cout << "recursive two particle calculation: " << twoRecursion.Re() << endl;
          cout << "conventional four particle calculation: " << four << endl;
          cout << "recursive four particle calculation: " << fourRecursion.Re() << endl;
          cout << "conventional six particle calculation: " << six << endl;
          cout << "recursive six particle calculation: " << sixRecursion.Re() << endl;
          cout << "recursive eight particle calculation: " << eightRecursion.Re() << endl;
        }

      for(int itrk=0; itrk<mult; itrk++)
	{
	  int charge = d_charge[itrk];
	  double phi = d_phi[itrk];
	  double pt = d_pt[itrk];
	  double v2 = d_v2[itrk];
	  th1d_phi->Fill(phi);
	  th1d_pt->Fill(pt);
	  double v2_track_true = cos(2*phi-2*psi2true);
	  double v2_track_reco = cos(2*phi-2*psi2reco);
	  tp1d_v2pT_true->Fill(pt,v2_track_true);
	  if(itrk>mult/2) tp1d_v2pT_reco->Fill(pt,v2_track_reco);
	  tp1d_v2pT_pure->Fill(pt,v2);
	  tp1d_YZ_v2pT_true->Fill(pt,v2_track_true*pt);
	  tp1d_YZ_v2pT_reco->Fill(pt,v2_track_reco*pt);

	  double x  = cos(phi);
	  double y  = sin(phi);
	  double xoff = offset;
	  double yoff = offset;
	  double xprime = x - xoff;
	  double yprime = y - yoff;
	  double phiprime = atan2(yprime,xprime);

	  double v2_primetrack_reco = cos(2*phiprime-2*psi2reco);
	  if(itrk>mult/2) tp1d_v2pT_prim->Fill(pt,v2_primetrack_reco);

	  double v2_track_primereco = cos(2*phi-2*psi2primereco);
	  if(itrk>mult/2) tp1d_v2pT_primereco->Fill(pt,v2_track_primereco);

	  double v2_track_doubleprimereco = cos(2*phiprime-2*psi2primereco);
	  if(itrk>mult/2) tp1d_v2pT_doubleprimereco->Fill(pt,v2_track_doubleprimereco);

	  double u2x = cos(2*phi);
	  double u2y = sin(2*phi);

	  double twoprime = ( u2x*Q2x + u2y*Q2y - 1) / (mult - 1);
	  tp1d_d2pT->Fill(pt,twoprime);
	  tp1d_YZ_d2pT->Fill(pt,twoprime*pt);
	  tp1d_YZ_W->Fill(pt,pt);

	  double u4x = cos(4*phi);
	  double u4y = sin(4*phi);

	  double fourprime = calc4_track(u2x,u2y,u4x,u4y,Q2x,Q2y,Q4x,Q4y,mult);
	  tp1d_p4pT->Fill(pt,fourprime);

	  double fourprimeYZ = calc4_track_YZ(u2x,u2y,u4x,u4y,Q2x,Q2y,Q4x,Q4y,mult);
	  tp1d_YZ_p4pT->Fill(pt,fourprimeYZ*pt);

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



double calc4_event(double Xn, double Yn, double X2n, double Y2n, double M)
{

  double Qn2 = Xn*Xn+Yn*Yn;
  double Qn2d = Xn*Xn-Yn*Yn;

  double one   = Qn2*Qn2;
  double two   = X2n*X2n+Y2n*Y2n;
  double three = (2*(X2n*Qn2d + 2*Y2n*Xn*Yn));
  double four  = 2*(2*(M-2)*Qn2);
  double five  = 2*M*(M-3);

  double numerator = one + two - three - four + five;
  double denominator = M*(M-1)*(M-2)*(M-3);

  return numerator/denominator;

}



double calc4_track(double xn, double yn, double x2n, double y2n, double Xn, double Yn, double X2n, double Y2n, double M)
{

  double one   = (xn*Xn + yn*Yn)*(Xn*Xn + Yn*Yn);
  double two   = x2n*Xn*Xn - x2n*Yn*Yn + 2*y2n*Xn*Yn;
  double three = xn*Xn*X2n + xn*Yn*Y2n - yn*(X2n*Yn - Xn*Y2n);
  double four  = 2*M*(xn*Xn + yn*Yn);
  double five  = 2*(Xn*Xn + Yn*Yn);
  double six   = 7*(xn*Xn + yn*Yn);
  double seven = xn*Xn + yn*Yn;
  double eight = x2n*X2n + y2n*Y2n;
  double nine = 2*(xn*Xn + yn*Yn);

  double numerator = one - two - three - four - five + six - seven + eight + nine +2*M - 6;
  double denominator = (M-1)*(M-2)*(M-3);

  return numerator/denominator;

}



double calc4_event_YZ(double QTx, double QTy, double QT2x, double QT2y, double MQT)
{

  double one   = (QTx*QTx + QTy*QTy)*(QTx*QTx + QTy*QTy);
  double two   = (QT2x*QT2x + QT2y*QT2y);
  double three = 2*(QT2x*QTx*QTx + 2*QT2y*QTx*QTy - QT2x*QTy*QTy);
  double four  = 2*(2*(MQT-2)*(QTx*QTx + QTy*QTy));
  double five  = 2*MQT*(MQT-3);

  double cn4 = (one + two - three - four + five) / (MQT*(MQT-1)*(MQT-2)*(MQT-3));

  return cn4;

}


double calc4_track_YZ(double pnx, double pny, double p2nx, double p2ny, double QTx, double QTy, double QT2x, double QT2y, double MQT)
{

  int mp = 1;

  double one   = ((pnx*QTx + pny*QTy)*(QTx*QTx + QTy*QTy));
  double two   = (p2nx*QTx*QTx - p2nx*QTy*QTy + 2*p2ny*QTx*QTy);
  double three = (pnx*QTx*QT2x - pny*QTy*QT2x + pnx*QTy*QT2y + pny*QTx*QT2y);
  double four  = 2*MQT*(pnx*QTx + pny*QTy);
  double five  = 2*mp*(QTx*QTx + QTy*QTy);
  double six   = 7*(pnx*QTx + pny*QTy);
  double seven = (QTx*pnx + QTy*pny);
  double eight = (p2nx*QT2x + p2ny*QT2y);
  double nine  = 2*(pnx*QTx + pny*QTy);

  double dn4 = (one - two - three - four - five + six - seven + eight + nine + 2*mp*MQT - 6*mp) / (mp*(MQT-1)*(MQT-2)*(MQT-3));

  return dn4;

}

double calc6_event(TComplex& qn, TComplex& q2n, TComplex& q3n, float M)
{

  if ( M < 6 ) return -9999;

  // TComplex qn, q2n, q3n;
  // qn = TComplex(Q2x,Q2y);
  // q2n = TComplex(Q4x,Q4y);
  // q3n = TComplex(Q6x,Q6y);

  TComplex temp1;

  // first term
  // |Qn|^6 + 9*|Q2n|^2|Qn|^2 - 6 x Re[Q2n x Qn x Qn* x Qn* x Qn*] / (Mx(M-1)x(M-2)x(M-3)x(M-4)x(M-5)
  double term1a = TMath::Power((qn*TComplex::Conjugate(qn)),3);
  double term1b = 9.0 * q2n*TComplex::Conjugate(q2n) * qn*TComplex::Conjugate(qn);
  temp1 = q2n * qn * TComplex::Conjugate(qn) * TComplex::Conjugate(qn) * TComplex::Conjugate(qn);
  double term1c = -6.0 * temp1.Re();
  double term1 = (term1a+term1b+term1c)/(M*(M-1)*(M-2)*(M-3)*(M-4)*(M-5));

  // second term
  // 4 * [Re[Q3nQn*Qn*Qn*] - 3 Re[Q3nQ2n*Qn*]] / (M(M-1)(M-2)(M-3)(M-4)(M-5)
  temp1 = q3n * TComplex::Conjugate(qn) * TComplex::Conjugate(qn) * TComplex::Conjugate(qn);
  double term2a = temp1.Re();
  temp1 = q3n * TComplex::Conjugate(q2n) * TComplex::Conjugate(qn);
  double term2b = -3.0 * temp1.Re();
  double term2 = 4.0 * (term2a+term2b)/(M*(M-1)*(M-2)*(M-3)*(M-4)*(M-5));

  // third term
  // +2 * (9*(M-4)*Re[Q2nQn*qn*] + 2 |Q3n|^2) / ((M(M-1)(M-2)(M-3)(M-4)(M-5))
  temp1 = q2n*TComplex::Conjugate(qn)*TComplex::Conjugate(qn);
  double term3a = 9.0*(M-4)*temp1.Re();
  double term3b = 2.0*q3n*TComplex::Conjugate(q3n);
  double term3 = 2.0 * (term3a + term3b) / (M*(M-1)*(M-2)*(M-3)*(M-4)*(M-5));

  // fourth term
  //double term4 = -9.0 * (TMath::Power(qn*TComplex::Conjugate(qn),2)+q2n*TComplex::Conjugate(q2n)) / (M*(M-1)*(M-2)*(M-3)*(M-5));
  double term4 = -9.0 * (TMath::Power(qn*TComplex::Conjugate(qn),2)+q2n*TComplex::Conjugate(q2n)) ;
  term4 /= (M*(M-1)*(M-2)*(M-3)*(M-5));

  // fifth term
  //double term5 = 18.0 * qn*TComplex::Conjugate(qn) / (M*(M-1)*(M-3)*(M-4));
  double term5 = 18.0 * qn*TComplex::Conjugate(qn) ;
  term5 /=  (M*(M-1)*(M-3)*(M-4));

  // sixth term
  double term6 = -6.0/((M-1)*(M-2)*(M-3));

  // cos(n(phi1+phi2+phi3-phi4-phi5-phi6))
  double six = term1 + term2 + term3 + term4 + term5 + term6;

  return six;

}

// --- from generic forumulas ----------------------------------------------------
//TComplex Recursion(int n, int* harmonic, int mult = 1, int skip = 0) // ROOT allows you to do dumb shit that C++ doesn't
// need to implement separate two argument version
TComplex Recursion(int n, int* harmonic)
{
  return Recursion(n,harmonic,1,0); // 1 and 0 are defaults from above
}
// remove default arguments in list
TComplex Recursion(int n, int* harmonic, int mult, int skip)
{
 // Calculate multi-particle correlators by using recursion (an improved faster version) originally developed by
 // Kristjan Gulbrandsen (gulbrand@nbi.dk).

  int nm1 = n-1;
  TComplex c(Q(harmonic[nm1], mult));
  if (nm1 == 0) return c;
  c *= Recursion(nm1, harmonic);
  if (nm1 == skip) return c;

  int multp1 = mult+1;
  int nm2 = n-2;
  int counter1 = 0;
  int hhold = harmonic[counter1];
  harmonic[counter1] = harmonic[nm2];
  harmonic[nm2] = hhold + harmonic[nm1];
  TComplex c2(Recursion(nm1, harmonic, multp1, nm2));
  int counter2 = n-3;
  while (counter2 >= skip) {
    harmonic[nm2] = harmonic[counter1];
    harmonic[counter1] = hhold;
    ++counter1;
    hhold = harmonic[counter1];
    harmonic[counter1] = harmonic[nm2];
    harmonic[nm2] = hhold + harmonic[nm1];
    c2 += Recursion(nm1, harmonic, multp1, counter2);
    --counter2;
  }
  harmonic[nm2] = harmonic[counter1];
  harmonic[counter1] = hhold;

  if (mult == 1) return c-c2;
  return c-double(mult)*c2;

}

TComplex Q(int n, int p)
{
 // Using the fact that Q{-n,p} = Q{n,p}^*.

 if(n>=0){return Qvector[n][p];}
 return TComplex::Conjugate(Qvector[-n][p]);

} // TComplex Q(int n, int p)
// -------------------------------------------------------------------------------
