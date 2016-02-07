import FWCore.ParameterSet.Config as cms

process = cms.Process('PreselHZZdata')


process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.debugModules=cms.untracked.vstring('muIsoDepositTkNew')

process.load("Configuration/StandardSequences/Reconstruction_cff")

# import of standard configurations
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

# process.GlobalTag.globaltag = "GR_R_42_V18::All"  
process.GlobalTag.globaltag = "GR_P_V32::All"

process.goodOfflinePrimaryVertices = cms.EDFilter("VertexSelector",
                                            src = cms.InputTag('offlinePrimaryVertices'),
                                            cut = cms.string('!isFake && isValid && ndof >= 4.0 && position.Rho < 2.0 && abs(z) < 24'),
                                            filter = cms.bool(True)
                                        )

usePAT='false'

# Preselection analysis sequence
if usePAT == 'true':
  process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsPreselectionPAT_2e2mu_cff')
else:
  process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsPreselection_data_noskim_cff')
  process.hTozzTo4leptonsCommonRootTreePresel.use2011EA = cms.untracked.bool(False)
  process.hTozzTo4leptonsCommonRootTreePresel.fillPUinfo = False
  process.patTrigger.processName=cms.string("HLT")
  process.hTozzTo4leptonsHLTInfo.TriggerResultsTag=cms.InputTag("TriggerResults","","HLT")
  process.hTozzTo4leptonsCommonRootTreePresel.triggerEvent = cms.InputTag("hltTriggerSummaryAOD","","HLT")
  process.hTozzTo4leptonsCommonRootTreePresel.flagHLTnames=cms.VInputTag(cms.InputTag("flagHLTDoubleMu7v1"),cms.InputTag("flagHLTL1DoubleMu0v1"),cms.InputTag("flagHLTDoubleMu4Acoplanarity03v1"),cms.InputTag("flagHLTL2DoubleMu23NoVertexv1"),cms.InputTag("flagHLTL2DoubleMu0v2"),cms.InputTag("flagHLTDoubleMu3v3"),cms.InputTag("flagHLTDoubleMu6v1"),cms.InputTag("flagHLTMu8Jet40v3"),cms.InputTag("flagHLTTripleMu5v2"),cms.InputTag("flagHLTEle17CaloIdLCaloIsoVLEle8CaloIdLCaloIsoVLv2"),cms.InputTag("flagHLTEle17CaloIdLCaloIsoVLv2"),cms.InputTag("flagHLTEle8CaloIdLCaloIsoVLJet40v2"),cms.InputTag("flagHLTEle17CaloIdTTrkIdVLCaloIsoVLTrkIsoVLEle8CaloIdTTrkIdVLCaloIsoVLTrkIsoVLv2"),cms.InputTag("flagHLTPhoton20CaloIdVTIsoTEle8CaloIdLCaloIsoVLv2"),cms.InputTag("flagHLTEle32CaloIdLCaloIsoVLSC17v2"),cms.InputTag("flagHLTTripleEle10CaloIdLTrkIdVLv2"),cms.InputTag("flagHLTEle17CaloIdLCaloIsoVLEle15HFLv2"),cms.InputTag("flagHLTEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8Mass30v2"),cms.InputTag("flagHLTEle8CaloIdLCaloIsoVLv2"),cms.InputTag("flagHLTDoubleEle10CaloIdLTrkIdVLEle10v2"),cms.InputTag("flagHLTEle8CaloIdLTrkIdVLv2"),cms.InputTag("flagHLTEle8v2"),cms.InputTag("flagHLTaccept"))
  process.hTozzTo4leptonsCommonRootTreePresel.triggerFilter = cms.string('hltDiMuonL3PreFiltered7')
  process.hTozzTo4leptonsCommonRootTreePresel.triggerFilterAsym = cms.vstring('hltDiMuonL3PreFiltered8','hltDiMuonL3p5PreFiltered8')


process.hTozzTo4leptonsSelectionPath = cms.Path(
  process.goodOfflinePrimaryVertices *
  process.hTozzTo4leptonsSelectionSequenceData *
  process.hTozzTo4leptonsCommonRootTreePresel
  )



process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsOutputModule_cff')
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsOutputModule_cff import *
process.hTozzTo4leptonsSelectionOutputModuleNew = hTozzTo4leptonsSelectionOutputModule.clone()
process.hTozzTo4leptonsSelectionOutputModuleNew.fileName = "hTozzToLeptons.root"
## process.hTozzTo4leptonsSelectionOutputModuleNew.SelectEvents.SelectEvents = cms.vstring('goodvertex','l1tcollpath'')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
             #               lumisToProcess = cms.untracked.VLuminosityBlockRange('138572:161-138572:168'),
#                             lumisToProcess = cms.untracked.VLuminosityBlockRange('191830:222'),
                            fileNames = cms.untracked.vstring(
#'/store/data/Run2012A/DoubleMu/AOD/PromptReco-v1/000/190/645/AA889834-8B82-E111-8815-002481E94C7E.root',
#'/store/data/Run2012A/DoubleMu/AOD/PromptReco-v1/000/190/645/0000D9F3-9082-E111-BC55-003048F117B6.root'
#'/store/data/Run2012A/DoubleElectron/AOD/PromptReco-v1/000/191/830/0ABFBFF7-B98C-E111-B83E-002481E0D7EC.root'
#'/store/data/Run2012A/DoubleElectron/AOD/PromptReco-v1/000/191/849/D84FFA26-E28C-E111-966E-001D09F2AD84.root'
#'file:/tmp/khurana/pickevents2012Extra.root'
#'file:/lustre/cms//store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/194/210/566AB49C-8DA0-E111-A649-003048F11942.root'
'file:hTozzTo4leptons_run.root'
#'file:/lustre/cms//store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/195/016/8C077D1C-50A9-E111-937A-001D09F29321.root'
)
                           )


# Endpath
process.o = cms.EndPath(process.hTozzTo4leptonsSelectionOutputModuleNew)

