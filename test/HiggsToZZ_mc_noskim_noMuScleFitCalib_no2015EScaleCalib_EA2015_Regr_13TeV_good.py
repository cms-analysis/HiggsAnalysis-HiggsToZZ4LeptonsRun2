
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
process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_74_V9', '')


process.goodOfflinePrimaryVertices = cms.EDFilter("VertexSelector",
                                            src = cms.InputTag('offlinePrimaryVertices'),
					    cut = cms.string('!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2'),
                                            filter = cms.bool(True)
                                        )
        

# Electron ordering in pT
process.hTozzTo4leptonsElectronOrdering = cms.EDProducer("HZZ4LeptonsElectronOrdering",
                                                         electronCollection = cms.InputTag("gsfElectrons")
                                                         )


process.load('HiggsAnalysis/HiggsToZZ4Leptons/muonCleanerBySegments_cfi')
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
#  'file:run.root'
#  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/pickevents_1_3_d5z.root',
#  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/pickevents_2_1_yGP.root',
#  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/pickevents_1_2_Ndh.root',
#  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/pickevents_2_1_d3s.root',
#  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/pickevents_3_1_h1h.root',
#  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/pickevents_4_1_TwB.root',
#  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/pickevents_5_2_g16.root',
#  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/pickevents_6_1_JtO.root',
#  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/pickevents_7_1_0kb.root',
#  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WplusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/pickevents_1_2_6qB.root',
#  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WplusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/pickevents_2_3_zsk.root',
#  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WplusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/pickevents_3_1_blk.root',
#  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WplusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/pickevents_4_1_0nb.root',
#  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WplusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/pickevents_5_1_flo.root',
#  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WplusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/pickevents_6_1_S6y.root',
#  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/ZH_HToZZ_4LFilter_M125_13TeV_powheg-minlo-HZJ_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/pickevents_1_1_R2j.root' 
#
  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/ttH_HToZZ_4LFilter_M125_13TeV_powheg_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_VERSION2/pickevents_1_1_FUg.root',
  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/ttH_HToZZ_4LFilter_M125_13TeV_powheg_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_VERSION2/pickevents_2_1_vR1.root',
#  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_VERSION2/pickevents_1_5_b5b.root',
#  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_VERSION2/pickevents_2_1_igP.root',
#  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/VBF_HToZZTo4L_M125_13TeV_powheg_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_VERSION2/pickevents_3_3_cdD.root',
#  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2_VERSION2/pickevents_1_1_D7J.root',
#  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2_VERSION2/pickevents_2_1_DpJ.root',
#  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WminusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2_VERSION2/pickevents_3_1_DVg.root',
#  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WplusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_VERSION2/pickevents_1_1_hp9.root',
#  'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/WplusH_HToZZTo4L_M125_13TeV_powheg-minlo-HWJ_JHUgen_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1_VERSION2/pickevents_2_1_2Rt.root'
  )


process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)


## # Endpath
#process.o = cms.EndPath ( process.hTozzTo4leptonsSelectionOutputModuleNew )
