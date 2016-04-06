
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
process.GlobalTag = GlobalTag(process.GlobalTag, '76X_mcRun2_asymptotic_v12', '')

# Random generator
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(1),
        engineName = cms.untracked.string('TRandom3')
    )
)

process.goodOfflinePrimaryVertices = cms.EDFilter("VertexSelector",
                                            src = cms.InputTag('offlinePrimaryVertices'),
					    cut = cms.string('!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2'),
                                            filter = cms.bool(True)
                                        )
        

process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsMuonCalibrator_cfi')
process.hTozzTo4leptonsMuonCalibrator.isData = cms.bool(False) 

process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')
process.calibratedElectrons.isMC = cms.bool(True)

process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsPreselection_data_noskim_cff') 
process.hTozzTo4leptonsHLTInfo.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
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


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
readFiles = cms.untracked.vstring(
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_10.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_11.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_12.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_13.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_14.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_15.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_16.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_17.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_18.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_19.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_1.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_20.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_21.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_22.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_23.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_24.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_25.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_2.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_3.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_4.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_5.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_6.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_7.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_8.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_201428/0000/pickevents_9.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_222700/0000/pickevents_1.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_222700/0000/pickevents_2.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_222700/0000/pickevents_3.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_76x_v1/crab_pickEvents/160309_222700/0000/pickevents_4.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_76x/crab_pickEvents/160309_222348/0000/pickevents_10.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_76x/crab_pickEvents/160309_222348/0000/pickevents_11.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_76x/crab_pickEvents/160309_222348/0000/pickevents_12.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_76x/crab_pickEvents/160309_222348/0000/pickevents_13.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_76x/crab_pickEvents/160309_222348/0000/pickevents_14.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_76x/crab_pickEvents/160309_222348/0000/pickevents_15.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_76x/crab_pickEvents/160309_222348/0000/pickevents_16.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_76x/crab_pickEvents/160309_222348/0000/pickevents_17.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_76x/crab_pickEvents/160309_222348/0000/pickevents_18.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_76x/crab_pickEvents/160309_222348/0000/pickevents_1.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_76x/crab_pickEvents/160309_222348/0000/pickevents_2.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_76x/crab_pickEvents/160309_222348/0000/pickevents_3.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_76x/crab_pickEvents/160309_222348/0000/pickevents_4.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_76x/crab_pickEvents/160309_222348/0000/pickevents_5.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_76x/crab_pickEvents/160309_222348/0000/pickevents_6.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_76x/crab_pickEvents/160309_222348/0000/pickevents_7.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_76x/crab_pickEvents/160309_222348/0000/pickevents_8.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_76x/crab_pickEvents/160309_222348/0000/pickevents_9.root'
  )


process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)


## # Endpath
# process.o = cms.EndPath ( process.hTozzTo4leptonsSelectionOutputModuleNew )
