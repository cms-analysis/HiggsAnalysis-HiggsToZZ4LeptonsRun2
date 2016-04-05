import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsMuonCalibrator = cms.EDProducer("HZZ4LeptonsMuonCalibrator",
    muonCollection = cms.InputTag("muons"),
    isData         = cms.bool(False)                               
)


