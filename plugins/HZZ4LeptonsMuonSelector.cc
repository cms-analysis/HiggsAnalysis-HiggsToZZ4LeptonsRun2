/**\class HZZ4LeptonsMuonSelector.cc
 *
 * Original Author:  Nicola De Filippis 
 *
 */

#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsMuonSelector.h"

#include "FWCore/Framework/interface/ESHandle.h"

// Muons:
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

// Candidate handling
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include "DataFormats/Common/interface/AssociationVector.h"

// system include files
#include <Math/VectorUtil.h>
#include <memory>
#include <vector>

using namespace edm;
using namespace std;
using namespace reco;

// constructor
HZZ4LeptonsMuonSelector::HZZ4LeptonsMuonSelector(const edm::ParameterSet& pset) {
  isGlobalMuon     = pset.getParameter<bool>("isGlobalMuon");
  isTrackerMuon    = pset.getParameter<bool>("isTrackerMuon");
  muonLabel        = consumes<edm::View<reco::Muon> >(pset.getParameter<edm::InputTag>("muonCollection"));
  muonPtMin        = pset.getParameter<double>("muonPtMin");
  muonEtaMax       = pset.getParameter<double>("muonEtaMax");

  string alias;
  produces<reco::MuonCollection>(); 

}


HZZ4LeptonsMuonSelector::~HZZ4LeptonsMuonSelector() {
 
}


//
// member functions
//
void HZZ4LeptonsMuonSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {



  // muons
  auto_ptr<reco::MuonCollection> Gmuon( new reco::MuonCollection );
  edm::Handle<edm::View<Muon> > muons;
  edm::View<reco::Muon>::const_iterator mIter;
    
  iEvent.getByToken(muonLabel, muons);

  // Loop over muons
  for (mIter = muons->begin(); mIter != muons->end(); ++mIter ) {

    bool matchglb=false, matchtrk=false;

    if(isGlobalMuon==true && mIter->isGlobalMuon()){
      matchglb=true;
    }
    else if (isGlobalMuon==false) {
      matchglb=true;
    }

    if(isTrackerMuon==true && mIter->isTrackerMuon()){
      matchtrk=true;
    }
    else if (isTrackerMuon==false) {
      matchtrk=true;
    }


    if (matchglb==true || matchtrk==true){
      
      if ( fabs( mIter->eta() ) > muonEtaMax ) continue;
      if ( mIter->pt() > muonPtMin) Gmuon->push_back( *mIter );   
    }

    cout << "Selected a muon with pT= " << mIter->pt() << endl;
  }

  
  const string iName = "";
  iEvent.put( Gmuon, iName );

}

