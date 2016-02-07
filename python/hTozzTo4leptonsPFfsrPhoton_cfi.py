import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsPFfsrPhoton = cms.EDProducer("HZZ4LeptonsPFfsrPhoton",
    pfCollection = cms.InputTag("particleFlow"),                               
)


