#ifndef  HZZ4LeptonsElectronSelector_h
#define  HZZ4LeptonsElectronSelector_h

/**\class HZZ4LeptonsElectronSelector
 *
 *
 * Original Author: Nicola De Filippis
 * Refine electron collection to begin with
 *
 */

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include <DataFormats/EgammaCandidates/interface/GsfElectron.h>
#include <DataFormats/EgammaCandidates/interface/GsfElectronFwd.h>

class HZZ4LeptonsElectronSelector : public edm::EDProducer {
 public:
  explicit HZZ4LeptonsElectronSelector(const edm::ParameterSet& );
  ~HZZ4LeptonsElectronSelector();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  edm::EDGetTokenT<edm::View<reco::GsfElectron> >  elecLabel;

  float elecPtMin;
  float elecEtaMax;
  
  int counterelectron;
 
};
  
#endif

