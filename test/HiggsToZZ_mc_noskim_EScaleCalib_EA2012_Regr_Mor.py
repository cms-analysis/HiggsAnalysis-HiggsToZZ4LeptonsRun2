import FWCore.ParameterSet.Config as cms

process = cms.Process('HZZ4lanalysis')

# Complete Preselection Sequence for 4l analysis

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')

# import of standard configurations
process.load("Configuration/StandardSequences/Reconstruction_cff")
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryDB_cff')
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


# process.GlobalTag.globaltag = "START42_V13::All"
# process.GlobalTag.globaltag = "START44_V10::All"
# process.GlobalTag.globaltag = "START52_V5::All"
process.GlobalTag.globaltag = "START53_V7A::All"
# process.GlobalTag.globaltag = "START53_V10::All"


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
 
  process.load('EgammaAnalysis.ElectronTools.electronRegressionEnergyProducer_cfi')
  process.eleRegressionEnergy.inputElectronsTag = cms.InputTag('hTozzTo4leptonsElectronOrdering')
  process.eleRegressionEnergy.inputCollectionType = cms.uint32(0)
  process.eleRegressionEnergy.useRecHitCollections = cms.bool(True)
  process.eleRegressionEnergy.produceValueMaps = cms.bool(True)
  # process.eleRegressionEnergy.debug = cms.untracked.bool(True)
  
  process.load("EgammaAnalysis/ElectronTools/calibratedElectrons_cfi")
  process.calibratedElectrons.inputElectronsTag = cms.InputTag('hTozzTo4leptonsElectronOrdering')
  process.calibratedElectrons.isMC = cms.bool(True)
  process.calibratedElectrons.inputDataset = cms.string("Summer12_DR53X_HCP2012")
  process.calibratedElectrons.isAOD = cms.bool(True)
  process.calibratedElectrons.updateEnergyError = cms.bool(True)
  process.calibratedElectrons.applyCorrections = cms.int32(1)
  process.calibratedElectrons.verbose = cms.bool(True)
  process.calibratedElectrons.synchronization = cms.bool(False)
  process.calibratedElectrons.smearingRatio = cms.double(0.607)
  # process.calibratedElectrons.synchronization = cms.bool(True) for sync
  # process.calibratedElectrons.smearingRatio = cms.double(0.0) for sync

  process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsPreselection_data_noskim_cff') 
  # process.hTozzTo4leptonsElectronSelector.electronCollection = cms.InputTag("calibratedElectrons")
  process.hTozzTo4leptonsElectronSelector.electronCollection = cms.InputTag("calibratedElectrons","calibratedGsfElectrons")
#  process.vetoElectrons.src = cms.InputTag("calibratedElectrons")  
  process.hTozzTo4leptonsHLTInfo.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
  process.patTrigger.processName=cms.string("HLT")
  process.hTozzTo4leptonsHLTInfo.TriggerResultsTag=cms.InputTag("TriggerResults","","HLT")
  process.hTozzTo4leptonsCommonRootTreePresel.use2011EA = cms.untracked.bool(False)
  process.hTozzTo4leptonsCommonRootTreePresel.triggerEvent  = cms.InputTag("hltTriggerSummaryAOD","","HLT")
  process.hTozzTo4leptonsCommonRootTreePresel.fillPUinfo = True
  process.hTozzTo4leptonsCommonRootTreePresel.fillHLTinfo = cms.untracked.bool(False)                                           
  process.hTozzTo4leptonsCommonRootTreePresel.triggerFilter = cms.string('hltDiMuonL3PreFiltered7')
  process.hTozzTo4leptonsCommonRootTreePresel.triggerFilterAsym = cms.vstring('hltDiMuonL3PreFiltered8','hltDiMuonL3p5PreFiltered8')
  process.hTozzTo4leptonsCommonRootTreePresel.fillMCTruth  = cms.untracked.bool(True)    


process.genanalysis= cms.Sequence(
        process.goodOfflinePrimaryVertices                  *
        process.hTozzTo4leptonsGenSequence                  *
#       process.hTozzTo4leptonsMCGenFilter2e2mu             *
#       process.hTozzTo4leptonsMCGenParticleListDrawer2e2mu *
        process.hTozzTo4leptonsMCDumper                     *                
        process.hTozzTo4leptonsMCCP                         )
        
process.hTozzTo4leptonsSelectionPath = cms.Path(
  process.genanalysis *
  process.hTozzTo4leptonsElectronOrdering *
  process.eleRegressionEnergy *
  process.calibratedElectrons *
  process.hTozzTo4leptonsSelectionSequenceData *
  process.hTozzTo4leptonsMatchingSequence *
  process.hTozzTo4leptonsCommonRootTreePresel
  )


process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsOutputModuleReduced_cff')
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsOutputModuleReduced_cff import *
process.hTozzTo4leptonsSelectionOutputModuleNew = hTozzTo4leptonsSelectionOutputModuleReduced.clone()
# process.hTozzTo4leptonsSelectionOutputModuleNew.fileName = "hTozzTo4leptons.root"

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#'file:/lustre/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V5-v2/0000/008CA5C0-CD7A-E111-987D-001A64789D14.root'
#'file:/lustre/cms/store/user/ndefilip/0CAA68E2-3491-E111-9F03-003048FFD760.root'
#'file:hTozzTo4leptons.root'
'file:alpha'
#'file:/lustre/cms/store/mc/Summer12_DR53X/GluGluToHToZZTo4L_M-1000_8TeV-powheg-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/BA5C279E-79FD-E111-89B5-00261894390C.root'
#'file:/lustre/cms/store/mc/Summer12/VBF_HToZZTo4L_M-200_8TeV-powheg-pythia6/AODSIM/PU_S7_START52_V9-v1/0000/C45663D6-3899-E111-9D12-0018F3D096A2.root'
#'file:/lustre/cms/store/user/ndefilip/GluGluToHToZZTo4L_M-125_8TeV-powheg-pythia6_PU_S10_START53_V7A_syncr.root'
#'file:/lustre/cms/store/mc/Summer12_DR53X/GluGluToHToZZTo4L_M-210_8TeV-powheg-pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/C0CF249B-9711-E211-9188-002618943904.root'
#'file:twoevents.root'
#'file:threeevents.root'  
                             )
                           )


# Endpath
# process.o = cms.EndPath ( process.hTozzTo4leptonsSelectionOutputModuleNew )

