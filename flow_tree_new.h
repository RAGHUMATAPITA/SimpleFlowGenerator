#include "TROOT.h"
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "math.h"
#include <vector>
#include <TVector3.h>

const double pi = TMath::Pi();
const double v2prefactor = 0.1;
char *FileOutput;
int d_mult;
float d_psi2;
int d_charge[50000];
float d_v2[50000];
float d_pt[50000];
float d_phi[50000];
