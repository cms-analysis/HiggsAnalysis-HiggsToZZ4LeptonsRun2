import FWCore.ParameterSet.Config as cms

process = cms.Process('step2')

# Complete Preselection Sequence for 4l analysis

# import of standard configurations
process.load("Configuration/StandardSequences/Reconstruction_cff")
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)


# process.GlobalTag.globaltag = "START42_V13::All"
# process.GlobalTag.globaltag = "START44_V10::All"
# process.GlobalTag.globaltag = "START52_V5::All"
process.GlobalTag.globaltag = "START53_V7G::All"

process.es_prefer_calotower       = cms.ESPrefer("CaloTowerGeometryFromDBEP","")
process.es_prefer_calocastor      = cms.ESPrefer("CastorGeometryFromDBEP","")
process.es_prefer_caloecalbarrel  = cms.ESPrefer("EcalBarrelGeometryFromDBEP","")
process.es_prefer_caloecalendcap  = cms.ESPrefer("EcalEndcapGeometryFromDBEP","")
process.es_prefer_caloecalpreshow = cms.ESPrefer("EcalPreshowerGeometryFromDBEP","")
process.es_prefer_calohcal        = cms.ESPrefer("HcalGeometryFromDBEP","")
process.es_prefer_calozdc         = cms.ESPrefer("ZdcGeometryFromDBEP","")



#PFJet ID
process.load('HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsPFJetSelector_cfi')
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsPFJetSelector_cfi import *

#PFJet Energy Corrections
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
from JetMETCorrections.Configuration.DefaultJEC_cff import *
process.load('JetMETCorrections.Configuration.JetCorrectionServices_cff')
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *
process.load('JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff')
from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *

process.ak5PFJetsCorrection = ak5PFJetsL1FastL2L3.clone()
process.ak5PFJetsCorrection.src = cms.InputTag('hTozzTo4leptonsPFJetSelector')
process.ak5PFJetsCorrectionData = ak5PFJetsL1FastL2L3Residual.clone()
process.ak5PFJetsCorrectionData.src = cms.InputTag('hTozzTo4leptonsPFJetSelector')

process.load('CMGTools.External.pujetidsequence_cff')
from CMGTools.External.pujetidsequence_cff import puJetId,puJetMva

process.recoPuJetIdMC = puJetId.clone(
  jets = cms.InputTag("ak5PFJetsCorrection"),
)
process.recoPuJetMvaMC = puJetMva.clone(
  jets = cms.InputTag("ak5PFJetsCorrection"),
  jetids = cms.InputTag("recoPuJetIdMC"),
)
process.recoPuJetIdMCsequence=cms.Sequence(process.recoPuJetIdMC * process.recoPuJetMvaMC) 

process.recoPuJetIdData = puJetId.clone(
  jets = cms.InputTag("ak5PFJetsCorrectionData"),
)
process.recoPuJetMvaData = puJetMva.clone(
  jets = cms.InputTag("ak5PFJetsCorrectionData"),
  jetids = cms.InputTag("recoPuJetIdData"),
)
process.recoPuJetIdDatasequence=cms.Sequence(process.recoPuJetIdData * process.recoPuJetMvaData)
        
process.hTozzTo4leptonsSelectionPath = cms.Path(
    process.hTozzTo4leptonsPFJetSelector +
    process.ak5PFJetsCorrection          +
    process.ak5PFJetsCorrectionData      +    
    process.recoPuJetIdMCsequence        +
    process.recoPuJetIdDatasequence
    )

process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsOutputModuleReduced_cff')
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsOutputModuleReduced_cff import *
process.hTozzTo4leptonsSelectionOutputModuleNew = hTozzTo4leptonsSelectionOutputModuleReduced.clone()
process.hTozzTo4leptonsSelectionOutputModuleNew.fileName = "hTozzTo4leptons.root"

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
'file:alpha'
#'file:/lustre/cms/store/user/defilip/DoubleMu/Data2012_step1_paper_DoubleMu_Run2012A-22Jan2013-v1/d7a85660e9c467da3e0de206f87fb5ab/hTozzTo4leptons_1_1_1kr.root'
                             )
                           )


# Endpath
process.o = cms.EndPath ( process.hTozzTo4leptonsSelectionOutputModuleNew )

