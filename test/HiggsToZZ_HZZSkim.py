import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")

## Source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    '/store/mc/Spring11/GluGluToHToZZTo4L_M-300_7TeV-powheg-pythia6/AODSIM/PU_S1_START311_V1G1-v1/0006/24F25368-D94E-E011-BA29-003048D45F76.root'
#    '/store/relval/CMSSW_4_2_1/Electron/RECO/GR_R_42_V10_RelVal_electron2010B-v1/0000/6E32B6B8-BB66-E011-A269-0030486790BE.root'
    )
                            )

# process
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_R_42_V10::All'
process.maxEvents  = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000


#load the skim
process.load('Configuration.Skimming.PDWG_HZZSkim_cff')
process.eePath   = cms.Path(process.zzdiElectronSequence)
process.mumuPath = cms.Path(process.zzdiMuonSequence)
process.emuPath  = cms.Path(process.zzeleMuSequence)

#output
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('HZZSkim.root'),
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('eePath', 'mumuPath', 'emuPath') ),
                               outputCommands = cms.untracked.vstring('keep *')
                               )
process.outpath = cms.EndPath(process.out)


#schedule
process.schedule = cms.Schedule( process.eePath, process.mumuPath, process.emuPath, process.outpath )
