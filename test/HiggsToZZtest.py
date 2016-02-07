
import FWCore.ParameterSet.Config as cms

process = cms.Process('TestTest')


process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')

process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsRunEventFilter_cfi')
process.hTozzTo4leptonsPath = cms.Path(process.hTozzTo4leptonsRunEventFilter)

#process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsMCGenParticleTreeDrawer_cfi')
#process.hTozzTo4leptonsPath = cms.Path(process.hTozzTo4leptonsMCGenParticleTreeDrawer)

#process.printTree1 = cms.EDAnalyzer("ParticleListDrawer",
#    src = cms.InputTag("genParticles"),
#    maxEventsToPrint  = cms.untracked.int32(33)
#)

#process.printTree2 = cms.EDAnalyzer("ParticleTreeDrawer",
#    src = cms.InputTag("genParticles"),
#    printP4 = cms.untracked.bool(False),
#    printPtEtaPhi = cms.untracked.bool(False),
#    printVertex = cms.untracked.bool(False),
#    printStatus = cms.untracked.bool(False),
#    printIndex  = cms.untracked.bool(False)
#)

#process.hTozzTo4leptonsPath = cms.Path(process.printTree1)


# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('run.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string('higgsToZZtest')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('hTozzTo4leptonsPath')
    )                               
 )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_AODSIM/pickevents_1_1_Re9.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_AODSIM/pickevents_2_3_iv3.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_AODSIM/pickevents_3_1_WDQ.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_AODSIM/pickevents_4_1_yES.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_AODSIM/pickevents_5_3_XQx.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_AODSIM/pickevents_6_4_AVW.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_AODSIM/pickevents_7_1_CS9.root',
'file:/lustre/cms/store/user/defilip/MonoHiggs/Syncr13TeV/GluGluToHToZZTo4L_M-125_13TeV-powheg-pythia6_Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1_AODSIM/pickevents_8_1_deB.root'
                             )
                           )


# Endpath
process.o = cms.EndPath ( process.output )
