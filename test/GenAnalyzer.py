import FWCore.ParameterSet.Config as cms

process = cms.Process("GenAnalysis")

process.load('Configuration/StandardSequences/Services_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.threshold = 'INFO'

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

#process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsMCGenFilter_cfi')
#process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsMCGenParticleListDrawer_cfi')
#process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsMCGenParticleTreeDrawer_cfi')

#process.printList = cms.EDFilter("ParticleListDrawer",
#                                src = cms.InputTag("genParticles"),
#                               maxEventsToPrint = cms.untracked.int32(10)
#                              )

#process.printTree = cms.EDFilter("ParticleTreeDrawer",
#                                src = cms.InputTag("genParticles"),
#                               printP4 = cms.untracked.bool(False),
#                              printPtEtaPhi = cms.untracked.bool(False),
#                             printVertex = cms.untracked.bool(True),
#                            printStatus = cms.untracked.bool(True),
#                           printIndex = cms.untracked.bool(False),
#                          status = cms.untracked.vint32(1,2,3)
#                         )

# select only Z, and save clones
process.Z = cms.EDFilter("PdgIdAndStatusCandViewSelector",
                          src = cms.InputTag("genParticles"),
                           pdgId = cms.vint32( 23 ),
                           status = cms.vint32(3)                       
                           )


# di-Z
process.diZ = cms.EDProducer("CandViewShallowCloneCombiner",
                                   decay = cms.string("Z Z"),
				   checkCharge= cms.bool(False),
                                   cut = cms.string("0< mass")
                                   )


# leptons
process.leptons = cms.EDFilter("PdgIdAndStatusCandViewSelector",
                                 src = cms.InputTag("genParticles"),
                                 pdgId = cms.vint32(11,13,15),
                                 status = cms.vint32(1)                     
                                 )

# di-leptons
process.dileptons = cms.EDProducer("CandViewShallowCloneCombiner",
                                   decay = cms.string("leptons@+ leptons@-"),
                                   cut = cms.string("0< mass && (daughter(0).charge>0 && daughter(1).charge<0)")
                                   )
## trileptons
process.trileptons = cms.EDProducer("CandViewShallowCloneCombiner",
                                    decay = cms.string("leptons@+ leptons@- leptons@+"),
                                    cut = cms.string("0 < mass && (daughter(0).charge>0 && daughter(1).charge<0 && daughter(2).charge>0)")
                                    )

#fourleptons
process.fourleptons = cms.EDProducer("CandViewShallowCloneCombiner",
                                     decay = cms.string("trileptons@+ leptons@-"),
                                     cut = cms.string("0 < mass && (daughter(0).charge>0 && daughter(1).charge<0)")
                                     )

process.GenAnalyzer = cms.EDAnalyzer("GenAnalyzer",
                                        src = cms.InputTag('genParticles')
                                        )


process.p = cms.Path(process.Z * process.diZ * process.leptons * process.dileptons * process.trileptons * process.fourleptons * process.GenAnalyzer)



#process.output = cms.OutputModule("PoolOutputModule",
#                                  fileName = cms.untracked.string('GenAnal.root'),
#                                  dataset = cms.untracked.PSet(
#                                    filterName = cms.untracked.string('higgsToZZtest')
# ),
#                                  SelectEvents = cms.untracked.PSet(
#    SelectEvents = cms.vstring('p')
#    )                               
#                                  )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#	'/store/mc/Spring11/GluGluToHToZZTo4L_M-200_7TeV-powheg-pythia6/AODSIM/PU_S1_START311_V1G1-v1/0004/243B12AC-8E4E-E011-A717-0025B3E05DCA.root'
#'file:hTozzTo4leptons_test.root'
#'file:/cmshome/nicola/tmp/CMSSW_4_1_4/src/ZZTo2e2mu_TuneZ2_7TeV_powheg_pythia6_cff_py_GEN.root'
'file:/cmshome/nicola/slc6/MonoHiggs/CMSSW_7_2_0/src/Hadronizer_MgmMatchTune4C_13TeV_madgraph_pythia8_Tauola_cff_py_GEN.root'
#'file:ZZ_madgraph_FC32337A-4F5B-E011-951A-0024E8768867.root'
    								)
                            )
