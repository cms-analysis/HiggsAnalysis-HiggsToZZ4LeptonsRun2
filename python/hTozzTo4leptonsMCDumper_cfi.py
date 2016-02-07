import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsMCDumper = cms.EDFilter("PdgIdAndStatusCandViewSelector",
                           src = cms.InputTag("genParticles"),
                           pdgId = cms.vint32( 25 ),
                           status = cms.vint32(3),
                           filter = cms.bool(False)                         
                           )
