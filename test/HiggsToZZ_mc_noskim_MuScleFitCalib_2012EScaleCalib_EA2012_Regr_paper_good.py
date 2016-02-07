
import FWCore.ParameterSet.Config as cms

process = cms.Process('MonoHiggs')

# Complete Preselection Sequence for 4l analysis

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/Geometry/GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration/EventContent/EventContent_cff')

# Random generator
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(1),
        engineName = cms.untracked.string('TRandom3')
    )
)


from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PHYS14_25_V1', '')


#process.es_prefer_calotower       = cms.ESPrefer("CaloTowerGeometryFromDBEP","")
#process.es_prefer_calocastor      = cms.ESPrefer("CastorGeometryFromDBEP","")
#process.es_prefer_caloecalbarrel  = cms.ESPrefer("EcalBarrelGeometryFromDBEP","")
#process.es_prefer_caloecalendcap  = cms.ESPrefer("EcalEndcapGeometryFromDBEP","")
#process.es_prefer_caloecalpreshow = cms.ESPrefer("EcalPreshowerGeometryFromDBEP","")
#process.es_prefer_calohcal        = cms.ESPrefer("HcalGeometryFromDBEP","")
#process.es_prefer_calozdc         = cms.ESPrefer("ZdcGeometryFromDBEP","")


process.goodOfflinePrimaryVerticesNew = cms.EDFilter("VertexSelector",
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

  process.load('RecoJets.JetProducers.kt4PFJets_cfi')
  process.kt6PFJets=process.kt4PFJets.clone()
  process.kt6PFJets.rParam = cms.double(0.6)
  process.kt6PFJets.doRhoFastjet = cms.bool(True) 
  process.kt6PFJets.Rho_EtaMax = cms.double(2.5)
  process.kt6PFJets.Ghost_EtaMax = cms.double(2.5)


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
  process.calibratedElectrons.isMC = cms.bool(True)
  process.calibratedElectrons.inputDataset = cms.string("Summer12_LegacyPaper")
  process.calibratedElectrons.isAOD = cms.bool(True)
  process.calibratedElectrons.updateEnergyError = cms.bool(True)
  process.calibratedElectrons.applyCorrections = cms.int32(1)
  process.calibratedElectrons.verbose = cms.bool(True)
  process.calibratedElectrons.correctionsType = cms.int32(2)
  process.calibratedElectrons.combinationType = cms.int32(3)
  process.calibratedElectrons.lumiRatio = cms.double(1.0)
  process.calibratedElectrons.synchronization = cms.bool(False)
  #process.calibratedElectrons.smearingRatio = cms.double(0.0)
  process.calibratedElectrons.applyLinearityCorrection = cms.bool(True)
    
  process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsPreselection_data_noskim_cff') 
  process.hTozzTo4leptonsElectronSelector.electronCollection = cms.InputTag("calibratedElectrons","calibratedGsfElectrons")
  # process.vetoElectrons.src = cms.InputTag("calibratedElectrons")  
  process.hTozzTo4leptonsHLTInfo.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
  process.patTrigger.processName=cms.string("HLT")
  process.hTozzTo4leptonsHLTInfo.TriggerResultsTag=cms.InputTag("TriggerResults","","HLT")
  process.hTozzTo4leptonsCommonRootTreePresel.use2011EA = cms.untracked.bool(False)
  process.hTozzTo4leptonsCommonRootTreePresel.triggerEvent  = cms.InputTag("hltTriggerSummaryAOD","","HLT")
  process.hTozzTo4leptonsCommonRootTreePresel.fillPUinfo = True
  process.hTozzTo4leptonsCommonRootTreePresel.fillHLTinfo = cms.untracked.bool(False)                                           
  process.hTozzTo4leptonsCommonRootTreePresel.triggerFilter = cms.string('hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q')
  #process.hTozzTo4leptonsCommonRootTreePresel.triggerFilterAsym = cms.vstring('hltDiMuonL3PreFiltered8','hltDiMuonL3p5PreFiltered8')
  process.hTozzTo4leptonsCommonRootTreePresel.fillMCTruth  = cms.untracked.bool(True)    
  process.hTozzTo4leptonsCommonRootTreePresel.isVBF  = cms.bool(False)

process.genanalysis= cms.Sequence(
        process.hTozzTo4leptonsGenSequence                  *
#       process.hTozzTo4leptonsMCGenFilter2e2mu             *
#       process.hTozzTo4leptonsMCGenParticleListDrawer2e2mu *
        process.hTozzTo4leptonsMCDumper                     *                
        process.hTozzTo4leptonsMCCP                         )
        
process.hTozzTo4leptonsSelectionPath = cms.Path(
#  process.goodOfflinePrimaryVerticesNew      *
  process.genanalysis *
  process.hTozzTo4leptonsElectronOrdering *
  process.kt6PFJets *
  process.eleRegressionEnergy *
 process.calibratedElectrons *
 process.hTozzTo4leptonsSelectionSequenceData *
 process.hTozzTo4leptonsMatchingSequence *
 process.hTozzTo4leptonsCommonRootTreePresel
  )


process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsOutputModule_cff')
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsOutputModule_cff import *
process.hTozzTo4leptonsSelectionOutputModuleNew = hTozzTo4leptonsSelectionOutputModule.clone()
process.hTozzTo4leptonsSelectionOutputModuleNew.fileName = "hTozzToLeptons.root"

#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
#readFiles = cms.untracked.vstring(
#'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_AODSIM/pickevents_1_1_Re9.root',
#'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_AODSIM/pickevents_2_3_iv3.root',
#'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_AODSIM/pickevents_3_1_WDQ.root',
#'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_AODSIM/pickevents_4_1_yES.root',
#'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_AODSIM/pickevents_5_3_XQx.root',
#'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_AODSIM/pickevents_6_4_AVW.root',
#'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_AODSIM/pickevents_7_1_CS9.root',
#'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_AODSIM/pickevents_8_1_deB.root'
#)

readFiles = cms.untracked.vstring(
'file:pickevents_1_1_69.root'
)

process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)


## # Endpath
process.o = cms.EndPath ( process.hTozzTo4leptonsSelectionOutputModuleNew )
