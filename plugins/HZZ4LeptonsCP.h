#ifndef HZZ4LeptonsCP_h
#define HZZ4LeptonsCP_h

/** \class HZZ4LeptonsCP
 *
 * Original Author:  Nicola De Filippis
 *
 */

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "TLorentzVector.h"
#include <TVector3.h>

// #include "RooRealVar.h"
// #include "JHU/MELA/src/src/AngularPdfFactory.cc"
// #include "JHU/MELA/src/PDFs/RooqqZZ_JHU.h"

// #include "JHU/MELA/src/scripts/MELA.h"

class HZZ4LeptonsCP : public edm::EDProducer {
public:
	//typedef math::XYZTLorentzVector TLorentzVector;
	
	explicit HZZ4LeptonsCP(const edm::ParameterSet&);
	~HZZ4LeptonsCP();
	
private:
	virtual void beginJob() ;
	virtual void produce(edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;
	
	double Distance( const reco::Candidate & c1, const reco::Candidate & c2 );
	double DistancePhi( const reco::Candidate & c1, const reco::Candidate & c2 );
	
	void calculateAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4Lep11, TLorentzVector thep4Lep12, TLorentzVector thep4Z2, TLorentzVector thep4Lep21, TLorentzVector thep4Lep22, double& costheta1, double& costheta2, double& phi, double& costhetastar, double& phistar1, double& phistar2, double& phistar12, double& phi1, double& phi2);


	//std::pair<double,double> likelihoodDiscriminant (double mZZ, double m1, double m2, double costhetastar, double costheta1, double costheta2, double phi, double phi1,double scaleFactor);

	//template <typename U>
	//void checkZorder(U& z1mass, U& z2mass,
	//		 U& costhetastar, U& costheta1,
	//		 U& costheta2, U& phi,
	//		 U& phistar1);

	bool debug;
	
	/***NVT****/
	edm::InputTag RECOcollName;  
	std::string decayChain_;
		
};
	
#endif
