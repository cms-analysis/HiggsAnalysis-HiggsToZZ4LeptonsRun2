#if !defined(__CINT__) || defined(__MAKECINT__)

#include "HZZ4LeptonsAnalysis_2e2mu.h"
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TROOT.h>
//#include <string>
#include <iostream>
#include <TSystem.h>
#include <TH2.h>
#include "TChain.h"
#include <stdlib.h>


#endif

using namespace std;


int main (int argc, char ** argv){

  string site=argv[7];
  cout << site << " configuration" <<endl;
  
  string dataconf=argv[8];
  cout << dataconf << " data configuration" <<endl;

  string mcconf=argv[9];
  cout << mcconf << " mc configuration" <<endl;


  // bkg
  ifstream fbkg;
  fbkg.open(argv[3]);

  const int nlines = atoi(argv[4]);

  string bkgsamples[nlines];
  float bkgninput[nlines];
  float bkgnhlt[nlines];
  float bkgnskim[nlines];
  float bkgxsection[nlines];
  
  cout << "Bkg samples" << endl;

  for(int i=0;i<nlines;i++){
    fbkg >> bkgsamples[i] >> bkgninput[i] >> bkgnhlt[i] >> bkgnskim[i] >> bkgxsection[i];
    cout << "Sample=" << bkgsamples[i] << " Ninput=" << bkgninput[i] << " NHLT=" << bkgnhlt[i] << " NSkim=" << bkgnskim[i] << " Xsection(pb)=" << bkgxsection[i] << endl;
  }


  // 
  float lumifb=0.;

  if (mcconf.find("Fall11")<5) lumifb=5051./1000.;
  if (mcconf.find("Summer12")<5){
    lumifb=19712./1000.;
    if (dataconf.find("2012_postICHEP")<15) lumifb=12200./1000.;
  }
  if (mcconf.find("Phys14")<5) lumifb=300000./1000.;


  Double_t mH=150.;
  cout << "mH= " << mH << endl;
    
  // Run on bkg
  bool runonbkg=true;
   
  if (runonbkg==true){

      bool itera=false;
      
      //if (mH<=160.){
      //	if (i==0) itera=true;  
      //}
      //else {
      itera=true;
      //}
      
      for(int i=0;i<nlines && itera==true;i++){
	
	string name= "roottree_leptons_"+bkgsamples[i]+".root";
	
	TString dirInput;
	if (site.find("CERN")<5){
          if (mcconf.find("Fall11")<5) dirInput="/castor/cern.ch/user/n/ndefilip/Paper/MCFall11";    // to run at CERN
	  else dirInput="/castor/cern.ch/user/n/ndefilip/Paper/MCSummer12";    // to run at CERN
	}
	else if (site.find("DESY")<5){
	  dirInput="/nfs/dust/test/cmsdas/school16/HZZ4lExercise/bkg"; //to run at DESY
	}
        else if (site.find("FNAL")<5 && mcconf.find("Fall11")<5){
          dirInput="dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/cmsdas/2014/HZZ4lExercise/bkg/Fall11";
        }
        else if (site.find("FNAL")<5 && mcconf.find("Summer12")<5 ){
          dirInput="dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/cmsdas/2014/HZZ4lExercise/bkg/Summer12";
        }
	else if (mcconf.find("Fall11")<5){
          dirInput="/lustre/cms/store/user/defilip/Fall11_445_paper_step_analysis_merged";
	}
	else if (mcconf.find("Summer12")<5){
          dirInput="/lustre/cms/store/user/defilip/Summer12_53X_paper_step_analysis_merged";
	}
        else if (mcconf.find("Phys14")<5){
	  dirInput="/lustre/cms/store/user/defilip/MonoHiggs/Phys14_720_merged";  
	}
	

	TString bkgFile=name;
	
	Char_t nome[300];
	sprintf(nome,"%s/%s",dirInput.Data(),bkgFile.Data());
	
	float weight= lumifb*(bkgxsection[i]*1000.*bkgnskim[i]/bkgninput[i])/bkgnskim[i];
        cout << "weight is " << weight << endl;	

	TFile *file3;
	TTree *tree3;

        if (bkgFile.Contains("TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola")){
          TChain* chain = new TChain("HZZ4LeptonsAnalysis","");
          chain->Add("/lustre/cms/store/user/defilip/MonoHiggs/Phys14_720_merged/roottree_leptons_TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola.root");
	  chain->Add("/lustre/cms/store/user/defilip/MonoHiggs/Phys14_720_merged/roottree_leptons_TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola_1.root");
	  chain->Add("/lustre/cms/store/user/defilip/MonoHiggs/Phys14_720_merged/roottree_leptons_TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola_2.root");
	  chain->Add("/lustre/cms/store/user/defilip/MonoHiggs/Phys14_720_merged/roottree_leptons_TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola_3.root");
          tree3 = chain;
        }
	//else if (bkgFile.Contains("TTTo2L2Nu2B_8TeV-powheg-pythia6")){
	//  TChain* chain = new TChain("HZZ4LeptonsAnalysis","");
	//  chain->Add("/lustre/cms/store/user/ndefilip/Summer12_52X_merged/roottree_leptons_TTTo2L2Nu2B_8TeV-powheg-pythia6.root");
	//  chain->Add("/lustre/cms/store/user/ndefilip/Summer12_52X_merged/roottree_leptons_TTTo2L2Nu2B_8TeV-powheg-pythia6_1.root");
	//  tree3 = chain;
	//}
	//else if (bkgFile.Contains("DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola")){
	//  TChain* chain = new TChain("HZZ4LeptonsAnalysis","");
	//  chain->Add("/lustre/cms/store/user/defilip/Fall11_445_paper_step_analysis_merged/roottree_leptons_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola.root");
	//  chain->Add("/lustre/cms/store/user/defilip/Fall11_445_paper_step_analysis_merged/roottree_leptons_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_1.root"); 
  	//  tree3 = chain;
	//}
	else {	          
	  file3 = TFile::Open(nome);
	  cout << "Read file with name: " << nome << endl;
	  tree3 = (TTree*)file3->Get("HZZ4LeptonsAnalysis");	  
	}

	HZZ4LeptonsAnalysis make3(tree3,weight,dataconf,mcconf);

	char *path=NULL;
	size_t size=300;
	path=getcwd(path,size);

	// cout << "Path is " << path << endl;
	sprintf(nome,"%s/output_%s.root",path,bkgsamples[i].c_str());

	// cout << "Output file name is " << nome << endl;
	make3.Loop(nome);
		
	cout << "Create file with name: " << nome << endl;
	delete tree3;
	file3 -> Close();
	
      }

  }

  return 0; 

}

