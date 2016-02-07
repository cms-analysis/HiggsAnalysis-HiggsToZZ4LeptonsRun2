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

process.GlobalTag.globaltag = "GR_R_42_V18::All"  


process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(False),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
                                  )


usePAT='false'

# Preselection analysis sequence
if usePAT == 'true':
  process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsPreselectionPAT_2e2mu_cff')
else:
  process.load("EgammaCalibratedGsfElectrons.CalibratedElectronProducers.calibratedGsfElectrons_cfi")
  process.calibratedGsfElectrons.isMC = cms.bool(False)
  process.calibratedGsfElectrons.inputDataset = cms.string("ReReco") 
  process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsPreselection_data_cff')
  process.hTozzTo4leptonsElectronSelector.electronCollection = cms.InputTag("calibratedGsfElectrons")
  process.vetoElectrons.src = cms.InputTag("calibratedGsfElectrons")

  process.hTozzTo4leptonsCommonRootTreePresel.fillPUinfo = False
  process.higgsToZZ4LeptonsHLTAnalysisData.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
  process.higgsToZZ4LeptonsHLTAnalysisData.HLTPaths= cms.vstring('HLT_DoubleMu7_v1','HLT_L1DoubleMu0_v1','HLT_DoubleMu4_Acoplanarity03_v1','HLT_L2DoubleMu23_NoVertex_v1','HLT_L2DoubleMu0_v2','HLT_DoubleMu3_v3','HLT_DoubleMu6_v1','HLT_Mu8_Jet40_v3','HLT_TripleMu5_v2','HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2','HLT_Ele17_CaloIdL_CaloIsoVL_v2','HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2','HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2','HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2','HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2','HLT_TripleEle10_CaloIdL_TrkIdVL_v2','HLT_Ele17_CaloIdL_CaloIsoVL_Ele15_HFL_v2','HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v2','HLT_Ele8_CaloIdL_CaloIsoVL_v2','HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_v2','HLT_Ele8_CaloIdL_TrkIdVL_v2','HLT_Ele8_v2');
  process.hTozzTo4leptonsHLTAnalysisData.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
  process.hTozzTo4leptonsHLTAnalysisData.HLTPaths= cms.vstring('HLT_DoubleMu7_v1','HLT_L1DoubleMu0_v1','HLT_DoubleMu4_Acoplanarity03_v1','HLT_L2DoubleMu23_NoVertex_v1','HLT_L2DoubleMu0_v2','HLT_DoubleMu3_v3','HLT_DoubleMu6_v1','HLT_Mu8_Jet40_v3','HLT_TripleMu5_v2','HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2','HLT_Ele17_CaloIdL_CaloIsoVL_v2','HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2','HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2','HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2','HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2','HLT_TripleEle10_CaloIdL_TrkIdVL_v2','HLT_Ele17_CaloIdL_CaloIsoVL_Ele15_HFL_v2','HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v2','HLT_Ele8_CaloIdL_CaloIsoVL_v2','HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_v2','HLT_Ele8_CaloIdL_TrkIdVL_v2','HLT_Ele8_v2');
  process.higgsToZZ4LeptonsSkimFilterData.useHLT  = cms.untracked.bool(False)
  process.patTrigger.processName=cms.string("HLT")
  process.hTozzTo4leptonsHLTInfo.TriggerResultsTag=cms.InputTag("TriggerResults","","HLT")
  process.hTozzTo4leptonsCommonRootTreePresel.triggerEvent = cms.InputTag("hltTriggerSummaryAOD","","HLT")
  process.hTozzTo4leptonsCommonRootTreePresel.flagHLTnames=cms.VInputTag(cms.InputTag("flagHLTDoubleMu7v1"),cms.InputTag("flagHLTL1DoubleMu0v1"),cms.InputTag("flagHLTDoubleMu4Acoplanarity03v1"),cms.InputTag("flagHLTL2DoubleMu23NoVertexv1"),cms.InputTag("flagHLTL2DoubleMu0v2"),cms.InputTag("flagHLTDoubleMu3v3"),cms.InputTag("flagHLTDoubleMu6v1"),cms.InputTag("flagHLTMu8Jet40v3"),cms.InputTag("flagHLTTripleMu5v2"),cms.InputTag("flagHLTEle17CaloIdLCaloIsoVLEle8CaloIdLCaloIsoVLv2"),cms.InputTag("flagHLTEle17CaloIdLCaloIsoVLv2"),cms.InputTag("flagHLTEle8CaloIdLCaloIsoVLJet40v2"),cms.InputTag("flagHLTEle17CaloIdTTrkIdVLCaloIsoVLTrkIsoVLEle8CaloIdTTrkIdVLCaloIsoVLTrkIsoVLv2"),cms.InputTag("flagHLTPhoton20CaloIdVTIsoTEle8CaloIdLCaloIsoVLv2"),cms.InputTag("flagHLTEle32CaloIdLCaloIsoVLSC17v2"),cms.InputTag("flagHLTTripleEle10CaloIdLTrkIdVLv2"),cms.InputTag("flagHLTEle17CaloIdLCaloIsoVLEle15HFLv2"),cms.InputTag("flagHLTEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8Mass30v2"),cms.InputTag("flagHLTEle8CaloIdLCaloIsoVLv2"),cms.InputTag("flagHLTDoubleEle10CaloIdLTrkIdVLEle10v2"),cms.InputTag("flagHLTEle8CaloIdLTrkIdVLv2"),cms.InputTag("flagHLTEle8v2"),cms.InputTag("flagHLTaccept"))
  process.hTozzTo4leptonsCommonRootTreePresel.triggerFilter = cms.string('hltDiMuonL3PreFiltered7')
  process.hTozzTo4leptonsCommonRootTreePresel.triggerFilterAsym = cms.vstring('hltDiMuonL3PreFiltered8','hltDiMuonL3p5PreFiltered8')


# process.hTozzTo4leptonsSelectionPath = cms.Path(process.L1T1coll * process.primaryVertexFilter * process.noscraping * process.hTozzTo4leptonsSelectionSequenceData * process.edmLumi)

process.hTozzTo4leptonsSelectionPath = cms.Path(process.calibratedGsfElectrons * process.hTozzTo4leptonsSelectionSequenceData * process.hTozzTo4leptonsCommonRootTreePresel)
# process.hTozzTo4leptonsSelectionPath = cms.Path(process.primaryVertexFilter * process.noscraping * process.hTozzTo4leptonsSelectionSequenceData)


#process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsOutputModule_cff')
#from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsOutputModule_cff import *
#process.hTozzTo4leptonsSelectionOutputModuleNew = hTozzTo4leptonsSelectionOutputModule.clone()
#process.hTozzTo4leptonsSelectionOutputModuleNew.fileName = "hTozzToLeptons.root"
## process.hTozzTo4leptonsSelectionOutputModuleNew.SelectEvents.SelectEvents = cms.vstring('goodvertex','l1tcollpath'')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            #lumisToProcess = cms.untracked.VLuminosityBlockRange('138572:113-138572:120'),
                            #lumisToProcess = cms.untracked.VLuminosityBlockRange('138572:137-138572:144'),
             #               lumisToProcess = cms.untracked.VLuminosityBlockRange('138572:161-138572:168'),
                            fileNames = cms.untracked.vstring(
#'/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/176/206/2A916503-D7E0-E011-9657-BCAEC5329707.root'
'/store/data/Run2011B/DoubleElectron/AOD/PromptReco-v1/000/180/250/F0BDD650-1F05-E111-BF4F-003048F118D2.root'
#'/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/176/163/6AE792CC-41DF-E011-A736-003048678110.root'
)
                           )


# Endpath
# process.o = cms.EndPath(process.hTozzTo4leptonsSelectionOutputModuleNew)

