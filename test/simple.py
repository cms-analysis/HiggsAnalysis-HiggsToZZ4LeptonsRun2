import FWCore.ParameterSet.Config as cms

process = cms.Process('TEST')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
 	'file:3295EF7C-2070-E411-89C4-7845C4FC35DB.root'
    )
)


#process.output = cms.OutputModule("PoolOutputModule",
#    splitLevel = cms.untracked.int32(0),
#    fileName = cms.untracked.string('test.root'),
#    dataset = cms.untracked.PSet(
#          dataTier = cms.untracked.string('GEN-SIM-RAW'),
#          filterName = cms.untracked.string('Nicola')
#    ),
#    SelectEvents = cms.untracked.PSet(
#        SelectEvents = cms.vstring('')
#    )    
#)


# process.out_step = cms.EndPath(process.output)


