/**\class HZZ4LeptonsPFJetSelector.cc
 *
 * Original Author:  Nicola De Filippis
 *
 */

#include "HiggsAnalysis/HiggsToZZ4Leptons/plugins/HZZ4LeptonsPFJetSelector.h"

#include "FWCore/Framework/interface/ESHandle.h"

// PFJets:
#include <DataFormats/JetReco/interface/PFJet.h>

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
HZZ4LeptonsPFJetSelector::HZZ4LeptonsPFJetSelector(const edm::ParameterSet& pset) {
  isLoosePFJetID     = pset.getParameter<bool>("isLoosePFJetID");
  //isMediumPFJetID     = pset.getParameter<bool>("isMediumPFJetID");
  //isTightPFJetID     = pset.getParameter<bool>("isTightPFJetID");
  pfjetsLabel        = pset.getParameter<edm::InputTag>("PFJetCollection");

  string alias;
  produces<reco::PFJetCollection>(); 

}


HZZ4LeptonsPFJetSelector::~HZZ4LeptonsPFJetSelector() {
 
}


//
// member functions
//
void HZZ4LeptonsPFJetSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {



  // PFJets
  auto_ptr<reco::PFJetCollection> GPFJet( new reco::PFJetCollection );
  edm::Handle<edm::View<PFJet> > pfjets;
  edm::View<reco::PFJet>::const_iterator mIter;
    
  iEvent.getByLabel(pfjetsLabel.label(), pfjets);

  // Loop over PFJets
  for (mIter = pfjets->begin(); mIter != pfjets->end(); ++mIter ) {

    //cout << "Jet with pT=" << mIter->pt() << " and eta=" << mIter->eta() << endl;
    if(isLoosePFJetID ==true){
      if (fabs(mIter->eta()) <= 3.0  ){
	if ( mIter->neutralHadronEnergyFraction() < 0.99 && 
	     mIter->neutralEmEnergyFraction() < 0.99 && 
	     (mIter->chargedMultiplicity()+mIter->neutralMultiplicity()) > 1) {
	  if (fabs(mIter->eta()) <= 2.4  ){
	    if (mIter->chargedHadronEnergyFraction() > 0. && 
		mIter->chargedMultiplicity() > 0. && 
		mIter->chargedEmEnergyFraction() < 0.99 ){
	      cout << "Jet passing the loose ID with pT=" << mIter->pt() << " and eta=" << mIter->eta() << endl; 
	      GPFJet->push_back( *mIter );
	    }
	  }
	  else {
	    cout << "Jet passing the loose ID with pT=" << mIter->pt() << " and eta=" << mIter->eta() << endl;
	    GPFJet->push_back( *mIter );
	  }
	}
      }
      else if (mIter->neutralEmEnergyFraction() < 0.90 && 
	       mIter->neutralMultiplicity() >10 ){
	cout << "Jet passing the loose ID with pT=" << mIter->pt() << " and eta=" << mIter->eta() << endl; 
	GPFJet->push_back( *mIter );
      }
    }
    else{
      GPFJet->push_back( *mIter );
    }
  }

    const string iName = "";
    iEvent.put( GPFJet, iName );

}

