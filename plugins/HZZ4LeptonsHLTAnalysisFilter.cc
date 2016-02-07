// system include files
#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsHLTAnalysisFilter.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <memory>
#include <iostream>
#include <fstream>
#include <vector>
#include "TString.h"
#include "FWCore/Framework/interface/FileBlock.h"
#include "DataFormats/MuonReco/interface/Muon.h"

// namespaces
using namespace edm;
using namespace std;

// Constructor
HZZ4LeptonsHLTAnalysisFilter::HZZ4LeptonsHLTAnalysisFilter(const edm::ParameterSet& pset) {

  HLTInfoFired              = pset.getParameter<edm::InputTag>("HLTInfoFired");
  
}


// Destructor
HZZ4LeptonsHLTAnalysisFilter::~HZZ4LeptonsHLTAnalysisFilter() {

}


// Filter Run Event
bool HZZ4LeptonsHLTAnalysisFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<vector<std::string> > HLTfired_;
  iEvent.getByLabel(HLTInfoFired,HLTfired_);
  
  vector<string> HLTimported;
  string tmpstring="";
  
  for (vector<string>::const_iterator cand=HLTfired_->begin(); cand!=HLTfired_->end(); ++cand){
    unsigned int i=cand-HLTfired_->begin();
    HLTimported.push_back(cand->c_str());
    string newstr=HLTimported.at(i) + ":" + tmpstring;
    tmpstring=newstr;
  }
  
  cout << "HLTFiredString= " << tmpstring.c_str() << endl;

  char HLTPathsFired[20000];
  sprintf(HLTPathsFired,tmpstring.c_str());
  
  stringstream ss (stringstream::in | stringstream::out);
  ss << HLTPathsFired;
  TString hlt(ss.str());

  TString out = inputfileName;

  bool debug=true;

  if( out.Contains("2013")){
      
    if( out.Contains("DoubleElectron")){
      
      if( debug ){ cout << "\n ** Step 2 (Trigger): "<< endl ;
	
	cout << "This is HLT in data" << endl;
	cout<<" HLTPathsFired... "<<hlt<<endl;
      }
      
      if(
	 !hlt.Contains("HLT_Ele23_Ele12_CaloId_TrackId_Iso_v") &&    // di-electron trigger
	 !hlt.Contains("HLT_Ele17_Ele12_Ele10_CaloId_TrackId_v") && // Triele
	 !hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v") &&     // di-muon trigger
	 !hlt.Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v") &&   // di-muon trigger
	 !hlt.Contains("HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v")   && // MuEle
	 !hlt.Contains("HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v")     // MuEle
	 ) {
	if( debug )cout << "Event not passing the HLT trigger paths" << endl;
	return false;
      }
      
    }
    else if( out.Contains("DoubleMu") ){
      
      if( debug ){ cout << "\n ** Step 2 (Trigger): "<< endl ;
	
	cout << "This is HLT in data" << endl;
	cout<<" HLTPathsFired... "<<hlt<<endl;
      }
      
      if(	 
	 !hlt.Contains("HLT_Ele23_Ele12_CaloId_TrackId_Iso_v") &&    // di-electron trigger
	 !hlt.Contains("HLT_Ele17_Ele12_Ele10_CaloId_TrackId_v") && // Triele
	 !hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v") &&     // di-muon trigger
	 !hlt.Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v") &&   // di-muon trigger
	 !hlt.Contains("HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v")   && // MuEle
	 !hlt.Contains("HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v")     // MuEle
		 ) {
	if( debug )cout << "Event not passing the HLT trigger paths" << endl;
	return false;
      }

      if(
	 hlt.Contains("HLT_Ele23_Ele12_CaloId_TrackId_Iso_v") || 
	 hlt.Contains("HLT_Ele17_Ele12_Ele10_CaloId_TrackId_v")  // di-ele and tri-ele excluded
	 ) {
	if( debug )cout << "Event not passing the HLT trigger vetos" << endl;
	return false;
      }
    }    
    else if( out.Contains("MuEG") ){
      if( debug ){ cout << "\n ** Step 2 (Trigger): 2011 DoubleElectron"<< endl ;
	
	cout << "This is HLT in data" << endl;
	cout<<" HLTPathsFired... "<<hlt<<endl;
      }
      
      if(
	 !hlt.Contains("HLT_Ele23_Ele12_CaloId_TrackId_Iso_v") &&   // di-electron trigger
	 !hlt.Contains("HLT_Ele17_Ele12_Ele10_CaloId_TrackId_v") && // Triele
	 !hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v") &&    // di-muon trigger
	 !hlt.Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v") &&  // di-muon trigger
	 !hlt.Contains("HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v")   && // MuEle
	 !hlt.Contains("HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v")    // MuEle
	 ) {
	if( debug )cout << "Event not passing the HLT trigger paths" << endl;
	return false;	      
      }
      
      if(
	 hlt.Contains("HLT_Ele23_Ele12_CaloId_TrackId_Iso_v") || 
	 hlt.Contains("HLT_Ele17_Ele12_Ele10_CaloId_TrackId_v") || 
	 hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v")   ||
	 hlt.Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v") // di-ele and tri-ele, di-muon excluded
	 ) {
	if( debug )cout << "Event not passing the HLT trigger vetos" << endl;
	return false;
      }
    }          
  }
  else if( out.Contains("Phys14")){
    if( debug ){ cout << "\n ** Step 2 (Trigger): "<< endl ;

      cout << "This is HLT in MC" << endl;
      cout<<" HLTPathsFired... "<<hlt<<endl;
    }
    
    if(
       !hlt.Contains("HLT_Ele23_Ele12_CaloId_TrackId_Iso_v") &&   // di-electron trigger
       !hlt.Contains("HLT_Ele17_Ele12_Ele10_CaloId_TrackId_v") && // Triele
       !hlt.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v") &&    // di-muon trigger
       !hlt.Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v") &&  // di-muon trigger
       !hlt.Contains("HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v")   && // MuEle
       !hlt.Contains("HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v")    // MuEle
       ) {
      if( debug )cout << "Event not passing the HLT trigger paths" << endl;
      return false;
    }
    
  }
  
  return true;

}

void HZZ4LeptonsHLTAnalysisFilter::respondToOpenInputFile(edm::FileBlock const& fb) {
  inputfileName = fb.fileName();
  cout << "Input Filename is=" << inputfileName.c_str() << endl;
  
}
