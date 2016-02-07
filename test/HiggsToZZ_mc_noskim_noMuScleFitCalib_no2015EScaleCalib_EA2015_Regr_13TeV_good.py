
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


from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PHYS14_25_V1', '')


#process.es_prefer_calotower       = cms.ESPrefer("CaloTowerGeometryFromDBEP","")
#process.es_prefer_calocastor      = cms.ESPrefer("CastorGeometryFromDBEP","")
#process.es_prefer_caloecalbarrel  = cms.ESPrefer("EcalBarrelGeometryFromDBEP","")
#process.es_prefer_caloecalendcap  = cms.ESPrefer("EcalEndcapGeometryFromDBEP","")
#process.es_prefer_caloecalpreshow = cms.ESPrefer("EcalPreshowerGeometryFromDBEP","")
#process.es_prefer_calohcal        = cms.ESPrefer("HcalGeometryFromDBEP","")
#process.es_prefer_calozdc         = cms.ESPrefer("ZdcGeometryFromDBEP","")


process.goodOfflinePrimaryVertices = cms.EDFilter("VertexSelector",
                                            src = cms.InputTag('offlinePrimaryVertices'),
                                            cut = cms.string('!isFake && ndof > 4.0 && position.Rho <= 2.0 && abs(z) <= 24'),
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

  
  process.load('MuonAnalysis/MuonAssociators/muonCleanerBySegments_cfi')
  process.cleanMuonsBySegments.src = cms.InputTag("muons")

    
  process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsPreselection_data_noskim_cff') 
  process.hTozzTo4leptonsElectronSelector.electronCollection = cms.InputTag("gedGsfElectrons")
  # process.vetoElectrons.src = cms.InputTag("calibratedElectrons")  
  process.hTozzTo4leptonsHLTInfo.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
  process.patTrigger.processName=cms.string("HLT")
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
  process.goodOfflinePrimaryVertices     *
  process.genanalysis *
  process.hTozzTo4leptonsElectronOrdering *
  process.hTozzTo4leptonsSelectionSequenceData *
  process.hTozzTo4leptonsMatchingSequence *
  process.hTozzTo4leptonsCommonRootTreePresel
  )


#process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsOutputModule_cff')
#from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsOutputModule_cff import *
#process.hTozzTo4leptonsSelectionOutputModuleNew = hTozzTo4leptonsSelectionOutputModule.clone()
#process.hTozzTo4leptonsSelectionOutputModuleNew.fileName = "hTozzToLeptons.root"

#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
readFiles = cms.untracked.vstring(
#'file:run.root'
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_AODSIM/pickevents_1_1_Re9.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_AODSIM/pickevents_2_3_iv3.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_AODSIM/pickevents_3_1_WDQ.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_AODSIM/pickevents_4_1_yES.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_AODSIM/pickevents_5_3_XQx.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_AODSIM/pickevents_6_4_AVW.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_AODSIM/pickevents_7_1_CS9.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_AODSIM/pickevents_8_1_deB.root'
)


process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)


## # Endpath
#process.o = cms.EndPath ( process.hTozzTo4leptonsSelectionOutputModuleNew )
