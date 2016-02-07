#ifndef HZZ4LeptonsIpToVtxProducer_h
#define HZZ4LeptonsIpToVtxProducer_h

/**\class HZZ4LeptonsIpToVtxProducer (from HZZ2e2muIpToVtxProducer)
 *
 * Original Author:  Alexis Pompili   - Bari
 *
 * Computes for each lepton the IP3D significance
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

// using namespace edm;

class HZZ4LeptonsIpToVtxProducer : public edm::EDProducer {

 public:

  explicit HZZ4LeptonsIpToVtxProducer(const edm::ParameterSet&);
  ~HZZ4LeptonsIpToVtxProducer();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  std::string decaychannel;
  edm::InputTag muonTag_, electronTag_,vertexTag_;
  bool useBeamSpot_;

  // PG and FRC 06-07-11 try to reduce printout!
	bool debug;
  //
};

#endif
