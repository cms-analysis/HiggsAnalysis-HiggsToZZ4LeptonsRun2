import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")

## Source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#    '/store/mc/Summer11/GluGluToHToZZTo4L_M-115_7TeV-powheg-pythia6/AODSIM/PU_S4_START42_V11-v1/0000/1210EDE5-4E98-E011-9BFE-00215E2223D0.root'
'file:/lustre/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V5-v2/0000/008CA5C0-CD7A-E111-987D-001A64789D14.root'
    )
                            )

# process
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
# import of standard configurations
process.load("Configuration/StandardSequences/Reconstruction_cff")
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

# process.GlobalTag.globaltag = "START42_V13::All"
# process.GlobalTag.globaltag = "START44_V10::All"
process.GlobalTag.globaltag = "START52_V5::All"
process.maxEvents  = cms.untracked.PSet( input = cms.untracked.int32(100) )

# message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000


#load the skim
process.load('Configuration.Skimming.PDWG_HZZSkim_cff')

# load HZZ analysis
usePAT='false'

# Preselection analysis sequence
if usePAT == 'true':
    process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsPreselectionPAT_2e2mu_cff')
else:
    process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsPreselection_data_hzzskim_cff') 
    process.hTozzTo4leptonsMuonSelector.isTrackerMuon=cms.bool(True)
    process.hTozzTo4leptonsMuonSelector.isGlobalMuon=cms.bool(False) 
    process.higgsToZZ4LeptonsHLTAnalysisData.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
    process.hTozzTo4leptonsHLTAnalysisData.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
    process.hTozzTo4leptonsHLTInfo.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
    process.higgsToZZ4LeptonsSkimFilterData.useHLT  = cms.untracked.bool(False)
    process.patTrigger.processName=cms.string("HLT")
    process.hTozzTo4leptonsHLTInfo.TriggerResultsTag=cms.InputTag("TriggerResults","","HLT")
    process.hTozzTo4leptonsCommonRootTreePresel.triggerEvent  = cms.InputTag("hltTriggerSummaryAOD","","HLT")
    process.hTozzTo4leptonsCommonRootTreePresel.fillPUinfo = True
    process.hTozzTo4leptonsCommonRootTreePresel.fillHLTinfo = cms.untracked.bool(False)                                           
    process.hTozzTo4leptonsCommonRootTreePresel.triggerFilter = cms.string('hltDiMuonL3PreFiltered7')
    process.hTozzTo4leptonsCommonRootTreePresel.triggerFilterAsym = cms.vstring('hltDiMuonL3PreFiltered8','hltDiMuonL3p5PreFiltered8')
    process.hTozzTo4leptonsCommonRootTreePresel.fillMCTruth  = cms.untracked.bool(True)  
    
    
process.genanalysis= cms.Sequence(
        process.hTozzTo4leptonsGenSequence                                 *
#       process.hTozzTo4leptonsMCGenFilter2e2mu             *
#       process.hTozzTo4leptonsMCGenParticleListDrawer2e2mu *
        process.hTozzTo4leptonsMCDumper                     *
        process.hTozzTo4leptonsMCCP                         )

process.hTozzTo4leptonsSelectionPath = cms.Path(
  process.genanalysis *
  process.hTozzTo4leptonsSelectionSequenceData *
  process.hTozzTo4leptonsMatchingSequence *
  process.hTozzTo4leptonsCommonRootTreePresel
  )


process.eePath   = cms.Path(process.zzdiElectronSequence * process.hTozzTo4leptonsSelectionSequenceData)
process.mumuPath = cms.Path(process.zzdiMuonSequence     * process.hTozzTo4leptonsSelectionSequenceData )
process.emuPath  = cms.Path(process.zzeleMuSequence      * process.hTozzTo4leptonsSelectionSequenceData )



#output
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('HZZSkim.root'),
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('eePath', 'mumuPath', 'emuPath') ),
                               outputCommands = cms.untracked.vstring('keep *')
                               )
process.outpath = cms.EndPath(process.out)


# Schedule
# process.schedule = cms.Schedule( process.eePath, process.mumuPath, process.emuPath, process.outpath )
