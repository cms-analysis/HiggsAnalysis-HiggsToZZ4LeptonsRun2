import FWCore.ParameterSet.Config as cms

process = cms.Process('step1')

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

# Random generator
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(1),
        engineName = cms.untracked.string('TRandom3')
    )
)


# process.GlobalTag.globaltag = "START42_V13::All"
# process.GlobalTag.globaltag = "START44_V10::All"
process.GlobalTag.globaltag = "START52_V5::All"

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


# Muon ghost cleaning
process.load('MuonAnalysis/MuonAssociators/muonCleanerBySegments_cfi')
from MuonAnalysis.MuonAssociators.muonCleanerBySegments_cfi import *
    
# Muon relaxed selection
process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsMuonSelector_cfi')
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMuonSelector_cfi import *
process.hTozzTo4leptonsMuonSelector=hTozzTo4leptonsMuonSelector.clone()
process.hTozzTo4leptonsMuonSelector.muonCollection = cms.InputTag("cleanMuonsBySegments")
process.hTozzTo4leptonsMuonSelector.isGlobalMuon=cms.bool(False)
process.hTozzTo4leptonsMuonSelector.isTrackerMuon=cms.bool(True)
process.hTozzTo4leptonsMuonSelector.muonPtMin=cms.double(3.)
process.hTozzTo4leptonsMuonSelector.muonEtaMax=cms.double(2.5)

        
process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsHLTInfo_cfi')
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsHLTInfo_cfi import *

process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsHLTAnalysisFilter_cfi')
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsHLTAnalysisFilter_cfi import *

process.hTozzTo4leptonsSelectionPath = cms.Path(
    process.goodOfflinePrimaryVertices +
    process.cleanMuonsBySegments +
    process.hTozzTo4leptonsMuonSelector +
    process.hTozzTo4leptonsHLTInfo +
    process.hTozzTo4leptonsHLTAnalysisFilter
    )

process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsOutputModuleReduced_cff')
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsOutputModuleReduced_cff import *
process.hTozzTo4leptonsSelectionOutputModuleNew = hTozzTo4leptonsSelectionOutputModuleReduced.clone()
process.hTozzTo4leptonsSelectionOutputModuleNew.fileName = "hTozzTo4leptons.root"

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(50))

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#'file:/lustre/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V5-v2/0000/008CA5C0-CD7A-E111-987D-001A64789D14.root'
#'file:/lustre/cms/store/user/ndefilip/0CAA68E2-3491-E111-9F03-003048FFD760.root'
#'file:hTozzTo4leptons.root'
#'file:/lustre/cms/store/mc/Summer12/VBF_HToZZTo4L_M-200_8TeV-powheg-pythia6/AODSIM/PU_S7_START52_V9-v1/0000/C45663D6-3899-E111-9D12-0018F3D096A2.root'
#'file:/lustre/cms/store/user/ndefilip/GluGluToHToZZTo4L_M-125_8TeV-powheg-pythia6_PU_S10_START53_V7A_syncr.root'
'file:/lustre/cms/store/mc/Summer12_DR53X/GluGluToHToZZTo4L_M-210_8TeV-powheg-pythia6/AODSIM/PU_S10_START53_V7A-v1/00000/C0CF249B-9711-E211-9188-002618943904.root'
#'file:twoevents.root'
#'file:threeevents.root'  
                             )
                           )


# Endpath
process.o = cms.EndPath ( process.hTozzTo4leptonsSelectionOutputModuleNew )

