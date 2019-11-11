#include "TROOT.h"  
#include "TObject.h" 
#include "TFile.h"
#include "TTree.h"
#include "math.h"
#include <vector>
#include <TVector3.h>
#include "flow_tree_new.h"
#include "TDatabasePDG.h"
struct lead_Struct
{
    Int_t        nevent;
    Int_t        nrun;
    Int_t        multi;
    Float_t      impactpar;
    Int_t        NpartP;
    Int_t        NpartT;
    Int_t        NELP;
    Int_t        NINP;
    Int_t        NELT;
    Int_t        NINT;
    Int_t        ev;
    Int_t        id[50000];
    Float_t      px[50000],py[50000],pz[50000],mass[50000],x[50000],y[50000],z[50000],t[50000];
};


void flow_tree_new()
{

  const Int_t kFile = 51;
  TString FileDir[kFile] = {"Pb_ampt0000001","Pb_ampt0000002","Pb_ampt0000003","Pb_ampt0000004","Pb_ampt0000005","Pb_ampt0000006","Pb_ampt0000007","Pb_ampt0000008","Pb_ampt0000009","Pb_ampt0000010","Pb_ampt0000011","Pb_ampt0000012","Pb_ampt0000013","Pb_ampt0000014","Pb_ampt0000015","Pb_ampt0000016","Pb_ampt0000017","Pb_ampt0000018","Pb_ampt0000019","Pb_ampt0000020","Pb_ampt0000021","Pb_ampt0000022","Pb_ampt0000023","Pb_ampt0000024","Pb_ampt0000025","Pb_ampt0000026","Pb_ampt0000027","Pb_ampt0000028","Pb_ampt0000029","Pb_ampt0000030","Pb_ampt0000031","Pb_ampt0000032","Pb_ampt0000033","Pb_ampt0000034","Pb_ampt0000035","Pb_ampt0000036","Pb_ampt0000037","Pb_ampt0000038","Pb_ampt0000040","Pb_ampt0000041","Pb_ampt0000042","Pb_ampt0000043","Pb_ampt0000044","Pb_ampt0000045","Pb_ampt0000046","Pb_ampt0000047","Pb_ampt0000048","Pb_ampt0000049","Pb_ampt0000050"};
  for(Int_t ifile = 0; ifile < kFile; ifile++)
    {
      
      char outfile[256];
      sprintf(outfile,"cumu_flow_pre%d.root",ifile);
      cout<<"outfile is ="<<outfile<<endl;      
      TString filePath = "AMPT_PbPb_Data/";
      filePath += FileDir[ifile];
      filePath += "/ampt_afterART";
      filePath += ".root";
      cout<<"file name is ="<<filePath<<endl;
      if(ifile == 2 )break;
      
      TFile *file = new TFile(filePath); 
      //TFile *file = new TFile("ampt_PPb_2_76TeV.root");
      //TFile *file = new TFile("ampt_afterART1.root");
      TTree *t =(TTree*)file->Get("AMPT");
      lead_Struct mystruct;
      
      TBranch *eve = t->GetBranch("Event");
      eve->SetAddress(&mystruct);
      
      TBranch* id = t->GetBranch("ID");
      id->SetAddress(&mystruct.id);
      
      TBranch* px = t->GetBranch("Px");
      px->SetAddress(&mystruct.px);
      
      TBranch* py = t->GetBranch("Py");
      py->SetAddress(&mystruct.py);
      
      TBranch* pz = t->GetBranch("Pz");
      pz->SetAddress(&mystruct.pz);
      
      TBranch* ma = t->GetBranch("Mass");
      ma->SetAddress(&mystruct.mass);
       
      TFile *file1 = new TFile(outfile,"recreate");
      
      TH1D *hmult = new TH1D("hmult","",30000,0,30000);
      
      TTree *tree = new TTree("simpletree","simpletree");
      tree->Branch("mult",&d_mult,"mult/I");
      tree->Branch("psi2",&d_psi2,"psi2/F");
      tree->Branch("charge",d_charge,"charge[mult]/I");
      tree->Branch("v2",d_v2,"v2[mult]/F");
      tree->Branch("pt",d_pt,"pt[mult]/F");
      tree->Branch("phi",d_phi,"phi[mult]/F");
      
      int entries = t->GetEntries();
      for(int nevent_sig=0; nevent_sig < entries; nevent_sig++)   //~~~~~~~~~~~~~~~~~~~ start event loop
	{
	  t->GetEntry(nevent_sig);
	  
	  Int_t Ntrack_sig = mystruct.multi;
	  
	  float psi2 = gRandom->Uniform(-pi,pi);
	  
	  for(Int_t multiplicity=0; multiplicity<Ntrack_sig ;multiplicity++)     //~~~~~~~~~start particle loop
	    {
	      float px_sig = mystruct.px[multiplicity];
	      float py_sig = mystruct.py[multiplicity];
	      float pz_sig = mystruct.pz[multiplicity];
	      int id = mystruct.id[multiplicity];
	      TVector3 particle_sig (px_sig,py_sig,pz_sig);
	      float pt_sig = particle_sig.Pt();
	      if( pt_sig > 2.0 ) continue;
	      float phi_sig = particle_sig.Phi();
	      float v2 = v2prefactor*pt_sig;
	      TDatabasePDG *db = TDatabasePDG::Instance();
	      TParticlePDG *part = 0x0;
	      if(!db) continue;
	      part = db->GetParticle(id);
	      if( !part) continue;
	      int charge = (part->Charge())/3;
	      if (charge == 0 ) continue;
	      phi_sig += psi2;
	      if(phi_sig>pi) phi_sig -= 2*pi;
	      if(phi_sig<-pi) phi_sig += 2*pi;
	      d_charge[multiplicity] = charge;
	      d_v2[multiplicity] = v2;
	      d_pt[multiplicity] = pt_sig;
	      d_phi[multiplicity] = phi_sig;
	    }
	  d_psi2 = psi2;
	  d_mult = Ntrack_sig;
	  tree->Fill();
	  hmult->Fill(d_mult);
	}
      
      file1->Write();
      file1->Close();
      
    }
  
}
