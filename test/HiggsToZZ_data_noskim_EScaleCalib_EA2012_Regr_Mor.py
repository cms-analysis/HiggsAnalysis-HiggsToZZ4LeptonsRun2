import FWCore.ParameterSet.Config as cms

process = cms.Process('PreselHZZdata')

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')

process.load("Configuration/StandardSequences/Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load('Configuration/Geometry/GeometryDB_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

# process.GlobalTag.globaltag = "GR_R_42_V18::All"  
# process.GlobalTag.globaltag = "GR_P_V32::All"
# process.GlobalTag.globaltag = "GR_P_V40::All"
# process.GlobalTag.globaltag = "GR_P_V41_AN2::All"
# process.GlobalTag.globaltag = "FT_P_V43E::All"
# process.GlobalTag.globaltag = "GR_P_V41_AN1::All"
# process.GlobalTag.globaltag = "GR_P_V42_AN2::All"
# process.GlobalTag.globaltag = "FT_53_V6C_AN1::All"
# process.GlobalTag.globaltag = "FT_53_V6_AN2::All"
# process.GlobalTag.globaltag = "FT_53_V10_AN2::All"
process.GlobalTag.globaltag = "GR_P_V42_AN2::All"
# process.GlobalTag.globaltag = "GR_P_V42C::All"

process.es_prefer_calotower       = cms.ESPrefer("CaloTowerGeometryFromDBEP","")
process.es_prefer_calocastor      = cms.ESPrefer("CastorGeometryFromDBEP","")
process.es_prefer_caloecalbarrel  = cms.ESPrefer("EcalBarrelGeometryFromDBEP","")
process.es_prefer_caloecalendcap  = cms.ESPrefer("EcalEndcapGeometryFromDBEP","")
process.es_prefer_caloecalpreshow = cms.ESPrefer("EcalPreshowerGeometryFromDBEP","")
process.es_prefer_calohcal        = cms.ESPrefer("HcalGeometryFromDBEP","")
process.es_prefer_calozdc         = cms.ESPrefer("ZdcGeometryFromDBEP","")


# Random generator
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(1),
        engineName = cms.untracked.string('TRandom3')
    )
)

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

  # Electron ordering in pT
  process.hTozzTo4leptonsElectronOrdering = cms.EDProducer("HZZ4LeptonsElectronOrdering",
                                                           electronCollection = cms.InputTag("gsfElectrons")
                                                           )
  
  # Electron Regression
  process.load('EgammaAnalysis.ElectronTools.electronRegressionEnergyProducer_cfi')
  process.eleRegressionEnergy.inputElectronsTag = cms.InputTag('hTozzTo4leptonsElectronOrdering')
  process.eleRegressionEnergy.inputCollectionType = cms.uint32(0)
  process.eleRegressionEnergy.useRecHitCollections = cms.bool(True)
  process.eleRegressionEnergy.produceValueMaps = cms.bool(True)
  # process.eleRegressionEnergy.debug = cms.untracked.bool(True)
  
  process.load("EgammaAnalysis/ElectronTools/calibratedElectrons_cfi")  
  process.calibratedElectrons.inputElectronsTag = cms.InputTag('hTozzTo4leptonsElectronOrdering')
  process.calibratedElectrons.isMC = cms.bool(False)
  process.calibratedElectrons.inputDataset = cms.string("Moriond2013")
  process.calibratedElectrons.isAOD = cms.bool(True)
  process.calibratedElectrons.updateEnergyError = cms.bool(True)
  process.calibratedElectrons.applyCorrections = cms.int32(1)
  process.calibratedElectrons.verbose = cms.bool(True)
  process.calibratedElectrons.synchronization = cms.bool(False)
  process.calibratedElectrons.smearingRatio = cms.double(0.607)
  # process.calibratedElectrons.synchronization = cms.bool(True) for sync
  # process.calibratedElectrons.smearingRatio = cms.double(0.0) for sync
  
  process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsPreselection_data_noskim_cff')
  process.hTozzTo4leptonsElectronSelector.electronCollection = cms.InputTag("calibratedElectrons","calibratedGsfElectrons")
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
  process.hTozzTo4leptonsElectronOrdering *
  process.eleRegressionEnergy *
  process.calibratedElectrons *
  process.hTozzTo4leptonsSelectionSequenceData *
  process.hTozzTo4leptonsCommonRootTreePresel
  )


process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsOutputModuleReduced_cff')
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsOutputModuleReduced_cff import *
process.hTozzTo4leptonsSelectionOutputModuleNew = hTozzTo4leptonsSelectionOutputModuleReduced.clone()
process.hTozzTo4leptonsSelectionOutputModuleNew.fileName = "hTozzTo4leptons.root"
## process.hTozzTo4leptonsSelectionOutputModuleNew.SelectEvents.SelectEvents = cms.vstring('goodvertex','l1tcollpath'')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            #lumisToProcess = cms.untracked.VLuminosityBlockRange('207922:85'),
                            #lumisToProcess = cms.untracked.VLuminosityBlockRange('207515:490'),
                            # lumisToProcess = cms.untracked.VLuminosityBlockRange('207905:422'),
                            # lumisToProcess = cms.untracked.VLuminosityBlockRange('208307:769'),
                            # lumisToProcess = cms.untracked.VLuminosityBlockRange('206446:368'),
# 	                    lumisToProcess = cms.untracked.VLuminosityBlockRange('207273:5,207273:7-207273:8,207273:17,207273:24-207273:26,207273:28,207273:30,207273:37,207273:98,207273:131,207273:140-207273:143,207273:148,207273:151,207273:153-207273:154,207273:157,207273:159-207273:160,207273:241,207273:256-207273:257,207273:259,207273:262-207273:263,207273:272'),
                            fileNames = cms.untracked.vstring(
#'/store/data/Run2012D/DoubleMu/AOD/PromptReco-v1/000/206/446/125284E8-CF25-E211-9830-0025901D626C.root'
#'file:hTozzTo4leptons_run.root'
#'file:hTozzTo4leptons.root'
'file:alpha',
#'/store/user/defilip/DoubleMu/Data2012_step1_Moriond/d7a85660e9c467da3e0de206f87fb5ab/hTozzTo4leptons_1258_2_OER.root'
#'file:///lustre/cms//store/user/defilip/DoubleElectron/Data2012_step1_Moriond/d7a85660e9c467da3e0de206f87fb5ab/hTozzTo4leptons_2433_1_RQe.root'
#'/store/data/Run2012D/DoubleMu/AOD/PromptReco-v1/000/208/353/82ADD753-883D-E211-A47E-003048D37694.root'
#'/store/data/Run2012D/DoubleMu/AOD/PromptReco-v1/000/207/273/465A4B08-F230-E211-943C-001D09F253C0.root'
#'file:/lustre/cms//store/data/Run2012B/DoubleMu/AOD/PromptReco-v1/000/195/016/8C077D1C-50A9-E111-937A-001D09F29321.root'
#'file:/tmp/simranjit/DoubleElectron.root'
#'file:/lustre/cms/store/data/Run2012D/DoubleMu/AOD/16Jan2013-v1/10000/72FDB725-B960-E211-9B07-E0CB4E5536F2.root'
# '/store/data/Run2012D/DoubleMu/AOD/16Jan2013-v1/10000/4EE889C2-A160-E211-B470-90E6BA442EEC.root'
#'/store/data/Run2012D/DoubleElectron/AOD/16Jan2013-v1/10000/5CD5E241-AB60-E211-BF70-003048679000.root'
#'/store/data/Run2012D/DoubleMu/AOD/16Jan2013-v1/10000/CA7A5EA3-9C60-E211-A429-E0CB4E1A117F.root'
#'/store/data/Run2012D/DoubleMu/AOD/PromptReco-v1/000/208/307/EE9D4877-153D-E211-8464-5404A63886C5.root'
#'/store/data/Run2012D/DoubleMu/AOD/PromptReco-v1/000/207/905/9EBA9FD4-3839-E211-8806-5404A63886B6.root'
#'/store/data/Run2012D/DoubleMu/AOD/PromptReco-v1/000/207/922/5440219B-9639-E211-85E1-003048D2C0F2.root'
#'/store/data/Run2012D/DoubleMu/AOD/16Jan2013-v1/10000/BAD81938-9A60-E211-B5BF-E0CB4E19F99E.root'
#'/store/data/Run2012D/DoubleElectron/AOD/PromptReco-v1/000/207/515/78FFD0C4-6134-E211-B2B4-5404A6388692.root'
#'/store/data/Run2012D/DoubleElectron/AOD/PromptReco-v1/000/207/492/0A1E0B39-BD33-E211-BD79-0030486780EC.root'
#'/store/data/Run2012D/DoubleElectron/AOD/PromptReco-v1/000/207/487/AEF24D0C-AC33-E211-9D27-001D09F24FEC.root'
#'/store/data/Run2012D/DoubleElectron/AOD/PromptReco-v1/000/207/492/A0FE63C1-CE33-E211-B195-003048D2BE08.root'
#'file:/afs/cern.ch/user/i/iross/public/ee187920671.root'
#'file:/lustre/cms//store/data/Run2012D/DoubleElectron/AOD/16Jan2013-v1/10000/66E6A449-9960-E211-9513-003048678FA0.root'
)
                           )


# Endpath
# process.o = cms.EndPath(process.hTozzTo4leptonsSelectionOutputModuleNew)

