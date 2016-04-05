/**\class HZZ4LeptonsMuonCalibrator.cc
 *
 * Original Author:  Nicola De Filippis 
 *
 */

#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsMuonCalibrator.h"

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

#include "KaMuCa/Calibration/interface/KalmanMuonCalibrator.h"

// system include files
#include <Math/VectorUtil.h>
#include <memory>
#include <vector>


using namespace edm;
using namespace std;
using namespace reco;
using namespace math;

// constructor
HZZ4LeptonsMuonCalibrator::HZZ4LeptonsMuonCalibrator(const edm::ParameterSet& pset) {
  isData           = pset.getParameter<bool>("isData");
  muonLabel        = consumes<edm::View<reco::Muon> >(pset.getParameter<edm::InputTag>("muonCollection"));

  string alias;
  produces<reco::MuonCollection>(); 

}


HZZ4LeptonsMuonCalibrator::~HZZ4LeptonsMuonCalibrator() {
 
}


//
// member functions
//
void HZZ4LeptonsMuonCalibrator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {



  // muons
  auto_ptr<reco::MuonCollection> Gmuon( new reco::MuonCollection );
  edm::Handle<edm::View<Muon> > muons;
  edm::View<reco::Muon>::const_iterator mIter;
    
  iEvent.getByToken(muonLabel, muons);

  KalmanMuonCalibrator calibrator;
  double corrPt=0.,corrPtError=0.;
  double smearedPt=0., smearedPtError=0.;
    
  // Loop over muons
  for (mIter = muons->begin(); mIter != muons->end(); ++mIter ) {
    
    reco::Muon* calibmu = mIter->clone();
    
    if (isData){    
      corrPt = calibrator.getCorrectedPt(mIter->pt(), mIter->eta(), mIter->phi(), mIter->charge());
      corrPtError = corrPt * calibrator.getCorrectedError(corrPt, mIter->eta(), mIter->bestTrack()->ptError()/corrPt );
      smearedPt=corrPt; // no smearing on data
      smearedPtError=corrPtError; // no smearing on data
    }
    else { // isMC - calibration from data + smearing 
      corrPt = calibrator.getCorrectedPt(mIter->pt(), mIter->eta(), mIter->phi(), mIter->charge());
      //corrPtError = corrPt * calibrator.getCorrectedError(corrPt, mIter->eta(), mIter->bestTrack()->ptError()/corrPt );
      smearedPt = calibrator.smear(corrPt, mIter->eta());
      smearedPtError = smearedPt * calibrator.getCorrectedErrorAfterSmearing(smearedPt, mIter->eta(), corrPtError /smearedPt );
    }

    cout << "Corrected Muon pT= " << smearedPt << " and pT error= " << smearedPtError << endl;
    reco::Candidate::PolarLorentzVector p4Polar_;
    p4Polar_ = reco::Candidate::PolarLorentzVector(smearedPt, mIter->eta(), mIter->phi(), mIter->mass());
    
    calibmu->setP4(p4Polar_);             
    Gmuon->push_back( *calibmu );
  }

  
  const string iName = "";
  iEvent.put( Gmuon, iName );

}

