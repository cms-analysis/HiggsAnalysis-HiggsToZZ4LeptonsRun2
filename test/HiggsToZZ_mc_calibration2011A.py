import FWCore.ParameterSet.Config as cms

process = cms.Process('PreselHZZdata')

# Complete Preselection Sequence for 2e2mu analysis

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')

process.load("Configuration/StandardSequences/Reconstruction_cff")

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

# process.GlobalTag.globaltag = 'START311_V1G1::All'
process.GlobalTag.globaltag = "START42_V13::All"

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedGsfElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(1),
        engineName = cms.untracked.string('TRandom3')
    )
)


usePAT='false'

# Preselection analysis sequence
if usePAT == 'true':
  process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsPreselectionPAT_2e2mu_cff')
else:
  process.load("EgammaCalibratedGsfElectrons.CalibratedElectronProducers.calibratedGsfElectrons_cfi")
  process.calibratedGsfElectrons.isMC = cms.bool(True)
  process.calibratedGsfElectrons.inputDataset = cms.string("Summer11")
  process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsPreselection_data_cff')
  process.hTozzTo4leptonsElectronSelector.electronCollection = cms.InputTag("calibratedGsfElectrons")
  process.vetoElectrons.src = cms.InputTag("calibratedGsfElectrons")

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
        process.hTozzTo4leptonsGenSequence                  *
#       process.hTozzTo4leptonsMCGenFilter2e2mu             *
#       process.hTozzTo4leptonsMCGenParticleListDrawer2e2mu *
        process.hTozzTo4leptonsMCDumper                     *                
        process.hTozzTo4leptonsMCCP                         )
        
process.hTozzTo4leptonsSelectionPath = cms.Path(
  process.genanalysis *
  process.calibratedGsfElectrons *
  process.hTozzTo4leptonsSelectionSequenceData *
  process.hTozzTo4leptonsMatchingSequence *
  process.hTozzTo4leptonsCommonRootTreePresel
  )


# process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsOutputModule_cff')
# from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsOutputModule_cff import *
# process.hTozzTo4leptonsSelectionOutputModuleNew = hTozzTo4leptonsSelectionOutputModule.clone()
# process.hTozzTo4leptonsSelectionOutputModuleNew.fileName = "hTozzToLeptons.root"

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(500))

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#'/store/mc/Summer11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v1/0000/3298E7AE-B19C-E011-86B7-90E6BA442F24.root'
#'file:/lustre/cms/store/mc/Summer11/GluGluToHToZZTo4L_M-190_7TeV-powheg-pythia6/AODSIM/PU_S4_START42_V11-v1/0000/709585F6-4498-E011-B509-00215E21EBCC.root'
#'/store/mc/Summer11/GluGluToHToZZTo4L_M-115_7TeV-powheg-pythia6/AODSIM/PU_S4_START42_V11-v1/0000/1210EDE5-4E98-E011-9BFE-00215E2223D0.root'
#'file:hTozzTo4leptons_run.root'
'file:GluGluToHToZZTo4L_Fall11_M-120.root'
                             )
                           )


# Endpath
# process.o = cms.EndPath ( process.hTozzTo4leptonsSelectionOutputModuleNew )

