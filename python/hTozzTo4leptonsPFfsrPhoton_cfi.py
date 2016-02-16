import FWCore.ParameterSet.Config as cms

# hTozzTo4leptonsPFfsrPhoton = cms.EDProducer("HZZ4LeptonsPFfsrPhoton",
#     pfCollection = cms.InputTag("particleFlow"),                               
# )

hTozzTo4leptonsPFfsrPhoton = cms.EDFilter(
    "GenericPFCandidateSelector",
    src = cms.InputTag("particleFlow","","RECO"),
    cut = cms.string("pdgId=22 && pt>2. && abs(eta)<2.4")
)

