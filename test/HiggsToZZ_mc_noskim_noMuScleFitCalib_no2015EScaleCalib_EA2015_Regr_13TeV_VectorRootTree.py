
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


# Global tag
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_74_V9', '')
# process.GlobalTag = GlobalTag(process.GlobalTag, '74X_dataRun2_Prompt_v0', '')



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

  
  #process.load('MuonAnalysis/MuonAssociators/muonCleanerBySegments_cfi')
  #process.cleanMuonsBySegments.src = cms.InputTag("muons")

    
  process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsPreselection_data_noskim_vector_cff') 
  process.hTozzTo4leptonsElectronSelector.electronCollection = cms.InputTag("gedGsfElectrons")
  # process.vetoElectrons.src = cms.InputTag("calibratedElectrons")  
  process.hTozzTo4leptonsHLTInfo.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
  process.patTrigger.processName=cms.string("HLT")
  process.hTozzTo4leptonsRootTreePresel.use2011EA = cms.untracked.bool(False)
  process.hTozzTo4leptonsRootTreePresel.triggerEvent  = cms.InputTag("hltTriggerSummaryAOD","","HLT")
  process.hTozzTo4leptonsRootTreePresel.fillPUinfo = True
  process.hTozzTo4leptonsRootTreePresel.fillHLTinfo = cms.untracked.bool(False)                                           
  process.hTozzTo4leptonsRootTreePresel.triggerFilter = cms.string('hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q')
  #process.hTozzTo4leptonsRootTreePresel.triggerFilterAsym = cms.vstring('hltDiMuonL3PreFiltered8','hltDiMuonL3p5PreFiltered8')
  process.hTozzTo4leptonsRootTreePresel.fillMCTruth  = cms.untracked.bool(True)    
  process.hTozzTo4leptonsRootTreePresel.isVBF  = cms.bool(False)

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
  process.hTozzTo4leptonsRootTreePresel
  )


#process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsOutputModule_cff')
#from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsOutputModule_cff import *
#process.hTozzTo4leptonsSelectionOutputModuleNew = hTozzTo4leptonsSelectionOutputModule.clone()
#process.hTozzTo4leptonsSelectionOutputModuleNew.fileName = "hTozzToLeptons.root"

#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
readFiles = cms.untracked.vstring(
#'file:zz.root'
#'file:pickevents_251562_143_130298882.root'
#'file:ZZTo4L_13TeV_powheg_pythia8_003264D6-C815-E511-8D7A-00A0D1EE26D0.root'
#'file:reco_RAW2DIGI_L1Reco_RECO.root'
#'/store/mc/RunIISpring15DR74/ZZ_TuneCUETP8M1_13TeV-pythia8/AODSIM/Startup25ns_EXOReReco_74X_Spring15_mcRun2_startup25ns_v0-v1/10000/04F5B26D-107B-E511-9BF2-20CF3027A5FD.root'
#'root://xrootd.ba.infn.it//store/mc/RunIISpring15DR74/ZZ_TuneCUETP8M1_13TeV-pythia8/AODSIM/Startup25ns_EXOReReco_74X_Spring15_mcRun2_startup25ns_v0-v1/10000/04F5B26D-107B-E511-9BF2-20CF3027A5FD.root'
'file:/lustre/cms/store/user/defilip/MonoHiggs/jobsreco_MZP1200_MA0300/reco_MZP1200_MA0300_2501.root'
)


process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)


## # Endpath
#process.o = cms.EndPath ( process.hTozzTo4leptonsSelectionOutputModuleNew )
