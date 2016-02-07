import FWCore.ParameterSet.Config as cms

process = cms.Process('HZZdata')

# Complete Preselection Sequence for 2e2mu analysis

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')

# import of standard configurations
process.load("Configuration/StandardSequences/Reconstruction_cff")
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/Geometry/GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

# Random generator
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(1),
        engineName = cms.untracked.string('TRandom3')
    )
)

# process.GlobalTag.globaltag = "GR_R_42_V18::All"
# process.GlobalTag.globaltag = "GR_P_V32::All"
# process.GlobalTag.globaltag = "GR_P_V40::All"
# process.GlobalTag.globaltag = "GR_P_V41_AN2::All"
# process.GlobalTag.globaltag = "FT_P_V43E::All"
# process.GlobalTag.globaltag = "GR_P_V41_AN1::All"
# process.GlobalTag.globaltag = "GR_P_V42_AN2::All"
# process.GlobalTag.globaltag = "FT_53_V6C_AN1::All"
# process.GlobalTag.globaltag = "FT_53_V6_AN2::All"
# process.GlobalTag.globaltag = "FT_53_V21_AN3::All"
# process.GlobalTag.globaltag = "GR_P_V42_AN2::All"
# process.GlobalTag.globaltag = "GR_P_V42C::All"
process.GlobalTag.globaltag = "FT_53_V21_AN4::All"

process.es_prefer_calotower       = cms.ESPrefer("CaloTowerGeometryFromDBEP","")
process.es_prefer_calocastor      = cms.ESPrefer("CastorGeometryFromDBEP","")
process.es_prefer_caloecalbarrel  = cms.ESPrefer("EcalBarrelGeometryFromDBEP","")
process.es_prefer_caloecalendcap  = cms.ESPrefer("EcalEndcapGeometryFromDBEP","")
process.es_prefer_caloecalpreshow = cms.ESPrefer("EcalPreshowerGeometryFromDBEP","")
process.es_prefer_calohcal        = cms.ESPrefer("HcalGeometryFromDBEP","")
process.es_prefer_calozdc         = cms.ESPrefer("ZdcGeometryFromDBEP","")


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
  process.load('EgammaAnalysis/ElectronTools/electronRegressionEnergyProducer_cfi')
  process.eleRegressionEnergy.inputElectronsTag = cms.InputTag('hTozzTo4leptonsElectronOrdering')
  # process.eleRegressionEnergy.inputElectronsTag = cms.InputTag('gsfElectrons')
  process.eleRegressionEnergy.energyRegressionType = cms.uint32(2)
  process.eleRegressionEnergy.inputCollectionType = cms.uint32(0)
  process.eleRegressionEnergy.useRecHitCollections = cms.bool(True)
  process.eleRegressionEnergy.produceValueMaps = cms.bool(True)
  process.eleRegressionEnergy.debug = cms.untracked.bool(True)
  
  process.load("EgammaAnalysis/ElectronTools/calibratedElectrons_cfi")
  process.calibratedElectrons.inputElectronsTag = cms.InputTag('hTozzTo4leptonsElectronOrdering')
  # process.calibratedElectrons.inputElectronsTag = cms.InputTag('gsfElectrons')
  process.calibratedElectrons.isMC = cms.bool(False)
  process.calibratedElectrons.inputDataset = cms.string("22Jan2013ReReco")
  process.calibratedElectrons.isAOD = cms.bool(True)
  process.calibratedElectrons.updateEnergyError = cms.bool(True)
  process.calibratedElectrons.applyCorrections = cms.int32(1)
  process.calibratedElectrons.verbose = cms.bool(True)
  process.calibratedElectrons.correctionsType = cms.int32(2)
  process.calibratedElectrons.combinationType = cms.int32(3)
  process.calibratedElectrons.lumiRatio = cms.double(1.0)
  process.calibratedElectrons.synchronization = cms.bool(False)
  # process.calibratedElectrons.smearingRatio = cms.double(0.0)
  process.calibratedElectrons.applyLinearityCorrection = cms.bool(True)

  process.MuScleFit = cms.EDProducer("MuScleFitMuonCorrector",
                                     src = cms.InputTag("muons"),
                                     debug = cms.bool(True),
                                     identifier = cms.string("Data2012_53X_ReReco"),
                                     applySmearing = cms.bool(False),
                                     fakeSmearing = cms.bool(False),
                                     )
  
  process.load('MuonAnalysis/MuonAssociators/muonCleanerBySegments_cfi')
  process.cleanMuonsBySegments.src = cms.InputTag("MuScleFit")
  
  process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsPreselection_data_noskim_cff') 
  process.hTozzTo4leptonsElectronSelector.electronCollection = cms.InputTag("calibratedElectrons","calibratedGsfElectrons")
  # process.vetoElectrons.src = cms.InputTag("calibratedElectrons")  
  process.hTozzTo4leptonsHLTInfo.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
  process.patTrigger.processName=cms.string("HLT")
  process.hTozzTo4leptonsHLTInfo.TriggerResultsTag=cms.InputTag("TriggerResults","","HLT")
  process.hTozzTo4leptonsCommonRootTreePresel.use2011EA = cms.untracked.bool(False)
  process.hTozzTo4leptonsCommonRootTreePresel.triggerEvent  = cms.InputTag("hltTriggerSummaryAOD","","HLT")
  process.hTozzTo4leptonsCommonRootTreePresel.fillPUinfo = False
  process.hTozzTo4leptonsCommonRootTreePresel.fillHLTinfo = cms.untracked.bool(False)                                           
  process.hTozzTo4leptonsCommonRootTreePresel.triggerFilter = cms.string('hltDiMuonL3PreFiltered7')
  process.hTozzTo4leptonsCommonRootTreePresel.triggerFilterAsym = cms.vstring('hltDiMuonL3PreFiltered8','hltDiMuonL3p5PreFiltered8')
  process.hTozzTo4leptonsCommonRootTreePresel.fillMCTruth  = cms.untracked.bool(False)    

        
process.hTozzTo4leptonsSelectionPath = cms.Path(
 process.goodOfflinePrimaryVertices      *
 process.hTozzTo4leptonsElectronOrdering *
 process.eleRegressionEnergy *
 process.calibratedElectrons *
 process.MuScleFit *
 process.hTozzTo4leptonsSelectionSequenceData *
 process.hTozzTo4leptonsCommonRootTreePresel
  )


process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsOutputModule_cff')
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsOutputModule_cff import *
process.hTozzTo4leptonsSelectionOutputModuleNew = hTozzTo4leptonsSelectionOutputModule.clone()
process.hTozzTo4leptonsSelectionOutputModuleNew.fileName = "hTozzToLeptons.root"

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#'file:/lustre/cms/store/mc/Summer12/VBF_HToZZTo4L_M-200_8TeV-powheg-pythia6/AODSIM/PU_S7_START52_V9-v1/0000/C45663D6-3899-E111-9D12-0018F3D096A2.root'
#'file:/lustre/cms/store/user/ndefilip/GluGluToHToZZTo4L_M-125_8TeV-powheg-pythia6_PU_S10_START53_V7A_syncr.root'
#'file:test.root'
#'file:/lustre/cms/store/user/defilip/FEEEEFFF-7FFB-E111-8FE2-002618943810.root'
#'file:/lustre/cms/store/user/defilip/4228BABE-70FA-E111-941B-001A92971B26.root'
#'file:/lustre/cms/store/user/defilip/DoubleMu/Data2012_step1_paper_DoubleMu_Run2012A-22Jan2013-v1/d7a85660e9c467da3e0de206f87fb5ab/hTozzTo4leptons_172_1_3N6.root'  
'file:alpha'
)
                           )


## # Endpath
#process.o = cms.EndPath ( process.hTozzTo4leptonsSelectionOutputModuleNew )

