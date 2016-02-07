import FWCore.ParameterSet.Config as cms

process = cms.Process('PreselHZZdata')


process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')

process.load("Configuration/StandardSequences/Reconstruction_cff")

# import of standard configurations
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

# process.GlobalTag.globaltag = "GR_R_42_V18::All"  
process.GlobalTag.globaltag = "GR_P_V32::All"  
# process.GlobalTag.globaltag = "GR_H_V29::All"

process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(False),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
                                  )

process.goodOfflinePrimaryVertices = cms.EDFilter("VertexSelector",
                                            src = cms.InputTag('offlinePrimaryVertices'),
                                            cut = cms.string('!isFake && isValid && ndof >= 4.0 && position.Rho < 2.0 && abs(z) < 24'),
                                            filter = cms.bool(True)
                                        )
        

#load the skim
process.load('Configuration.Skimming.PDWG_HZZSkim_cff')

usePAT='false'

# Preselection analysis sequence
if usePAT == 'true':
  process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsPreselectionPAT_2e2mu_cff')
else:
  process.load("EgammaCalibratedGsfElectrons.CalibratedElectronProducers.calibratedGsfElectrons_cfi")
  process.calibratedGsfElectrons.isMC = cms.bool(False)
  process.calibratedGsfElectrons.isAOD = cms.bool(True)
  process.calibratedGsfElectrons.updateEnergyError = cms.bool(True)
  process.calibratedGsfElectrons.inputDataset = cms.string("ICHEP2012") 
  process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsPreselection_data_hzzskim_cff')
  process.hTozzTo4leptonsElectronSelector.electronCollection = cms.InputTag("calibratedGsfElectrons")
#  process.vetoElectrons.src = cms.InputTag("calibratedGsfElectrons")
  process.hTozzTo4leptonsMuonSelector.isTrackerMuon=cms.bool(True)
  process.hTozzTo4leptonsMuonSelector.isGlobalMuon=cms.bool(False) 
  process.hTozzTo4leptonsCommonRootTreePresel.use2011EA = cms.untracked.bool(False)
  process.hTozzTo4leptonsCommonRootTreePresel.fillPUinfo = False
  process.patTrigger.processName=cms.string("HLT")
  process.hTozzTo4leptonsHLTInfo.TriggerResultsTag=cms.InputTag("TriggerResults","","HLT")
  process.hTozzTo4leptonsCommonRootTreePresel.triggerEvent = cms.InputTag("hltTriggerSummaryAOD","","HLT")
  process.hTozzTo4leptonsCommonRootTreePresel.flagHLTnames=cms.VInputTag(cms.InputTag("flagHLTDoubleMu7v1"),cms.InputTag("flagHLTL1DoubleMu0v1"),cms.InputTag("flagHLTDoubleMu4Acoplanarity03v1"),cms.InputTag("flagHLTL2DoubleMu23NoVertexv1"),cms.InputTag("flagHLTL2DoubleMu0v2"),cms.InputTag("flagHLTDoubleMu3v3"),cms.InputTag("flagHLTDoubleMu6v1"),cms.InputTag("flagHLTMu8Jet40v3"),cms.InputTag("flagHLTTripleMu5v2"),cms.InputTag("flagHLTEle17CaloIdLCaloIsoVLEle8CaloIdLCaloIsoVLv2"),cms.InputTag("flagHLTEle17CaloIdLCaloIsoVLv2"),cms.InputTag("flagHLTEle8CaloIdLCaloIsoVLJet40v2"),cms.InputTag("flagHLTEle17CaloIdTTrkIdVLCaloIsoVLTrkIsoVLEle8CaloIdTTrkIdVLCaloIsoVLTrkIsoVLv2"),cms.InputTag("flagHLTPhoton20CaloIdVTIsoTEle8CaloIdLCaloIsoVLv2"),cms.InputTag("flagHLTEle32CaloIdLCaloIsoVLSC17v2"),cms.InputTag("flagHLTTripleEle10CaloIdLTrkIdVLv2"),cms.InputTag("flagHLTEle17CaloIdLCaloIsoVLEle15HFLv2"),cms.InputTag("flagHLTEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8Mass30v2"),cms.InputTag("flagHLTEle8CaloIdLCaloIsoVLv2"),cms.InputTag("flagHLTDoubleEle10CaloIdLTrkIdVLEle10v2"),cms.InputTag("flagHLTEle8CaloIdLTrkIdVLv2"),cms.InputTag("flagHLTEle8v2"),cms.InputTag("flagHLTaccept"))
  process.hTozzTo4leptonsCommonRootTreePresel.triggerFilter = cms.string('hltDiMuonL3PreFiltered7')
  process.hTozzTo4leptonsCommonRootTreePresel.triggerFilterAsym = cms.vstring('hltDiMuonL3PreFiltered8','hltDiMuonL3p5PreFiltered8')


process.hTozzTo4leptonsSelectionPath = cms.Path(
  process.noscraping                                  *
  process.goodOfflinePrimaryVertices                  *
  process.calibratedGsfElectrons *
  process.hTozzTo4leptonsSelectionSequenceData *
  process.hTozzTo4leptonsCommonRootTreePresel 
  )

process.HZZ4ePath      = cms.Path( process.zz4eSequence   * process.hTozzTo4leptonsSelectionSequenceData)
process.HZZ2e2mPath    = cms.Path( process.zz2e2mSequence * process.hTozzTo4leptonsSelectionSequenceData) 
process.HZZ2m2ePath    = cms.Path( process.zz2m2eSequence * process.hTozzTo4leptonsSelectionSequenceData) 
process.HZZ4mPath      = cms.Path( process.zz4mSequence   * process.hTozzTo4leptonsSelectionSequenceData)
process.HZZem2ePath    = cms.Path( process.zzem2eSequence * process.hTozzTo4leptonsSelectionSequenceData)
process.HZZem2mPath    = cms.Path( process.zzem2mSequence * process.hTozzTo4leptonsSelectionSequenceData)

process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsOutputModule_cff')
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsOutputModule_cff import *
process.hTozzTo4leptonsSelectionOutputModuleNew = hTozzTo4leptonsSelectionOutputModule.clone()
process.hTozzTo4leptonsSelectionOutputModuleNew.fileName = "hTozzToLeptons.root"
process.hTozzTo4leptonsSelectionOutputModuleNew.SelectEvents.SelectEvents = cms.vstring('HZZ4ePath', 'HZZ2e2mPath', 'HZZ2m2ePath','HZZ4mPath', 'HZZem2ePath', 'HZZem2mPath')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
                            #lumisToProcess = cms.untracked.VLuminosityBlockRange('138572:113-138572:120'),
                            #lumisToProcess = cms.untracked.VLuminosityBlockRange('138572:137-138572:144'),
                            # lumisToProcess = cms.untracked.VLuminosityBlockRange('191090:1-191090:13,191090:32-191090:33,191090:70'),
                            # lumisToProcess = cms.untracked.VLuminosityBlockRange('191830:222'),
                            fileNames = cms.untracked.vstring(
#'/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/176/206/2A916503-D7E0-E011-9657-BCAEC5329707.root'
#'/store/data/Run2010A/Mu/AOD/Apr21ReReco-v1/0000/764547F6-5470-E011-BDED-1CC1DE1D1FE6.root'
#'/store/data/Run2012A/DoubleMu/AOD/PromptReco-v1/000/190/645/AA889834-8B82-E111-8815-002481E94C7E.root',
#'/store/data/Run2012A/DoubleMu/AOD/PromptReco-v1/000/190/645/0000D9F3-9082-E111-BC55-003048F117B6.root'
#'/store/data/Run2012A/DoubleElectron/AOD/PromptReco-v1/000/191/830/0ABFBFF7-B98C-E111-B83E-002481E0D7EC.root'
#'/store/data/Run2012A/DoubleMu/AOD/PromptReco-v1/000/193/557/E2898458-C199-E111-AF72-003048D3C944.root'
#'file:/lustre/cms//store/data/Run2012A/DoubleMu/AOD/PromptReco-v1/000/191/247/F28986C0-5888-E111-956D-0025901D626C.root'
#'/store/data/Run2012A/DoubleElectron/AOD/PromptReco-v1/000/191/090/C6A670B4-0887-E111-BB7D-0025B32034EA.root'
#'file:0265BEE4-3E9C-E111-8A43-BCAEC5329702.root'
#'file:hTozzTo4leptons_run.root'
'file:/lustre/cms//store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/194/210/566AB49C-8DA0-E111-A649-003048F11942.root'
)
                           )

# Schedule
# process.schedule = cms.Schedule( process.HZZ4ePath, process.HZZ2e2mPath, process.HZZ2m2ePath, process.HZZ4mPath, process.HZZem2ePath, process.HZZem2mPath )

# Endpath
# process.o = cms.EndPath(process.hTozzTo4leptonsSelectionOutputModuleNew)

