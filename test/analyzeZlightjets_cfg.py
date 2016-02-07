import FWCore.ParameterSet.Config as cms

process = cms.Process("MatchingAnalysis")

process.load('Configuration/StandardSequences/Services_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.threshold = 'INFO'
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('Gen_775_1.root'),
                            skipEvents = cms.untracked.vint32(0)
                            )

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load('HiggsAnalysis/HiggsToZZ4Leptons/hTozzTo4leptonsMCGenFilter_cfi')
from HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenFilter_cfi import *
import HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenFilter_cfi
process.hTozzTo4leptonsMCGenFilterNEW = HiggsAnalysis.HiggsToZZ4Leptons.hTozzTo4leptonsMCGenFilter_cfi.hTozzTo4leptonsMCGenFilter.clone()
process.hTozzTo4leptonsMCGenFilterNEW.HZZ4LeptonsMCFilterLeptonFlavour = cms.int32(7)


# process.printList = cms.EDFilter("ParticleListDrawer",
#                                src = cms.InputTag("genParticles"),
#                               maxEventsToPrint = cms.untracked.int32(10)
#                              )
# process.printTree = cms.EDFilter("ParticleTreeDrawer",
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

process.c = cms.EDFilter("PdgIdAndStatusCandViewSelector",
                         src = cms.InputTag("genParticles"),
                         pdgId = cms.vint32( 4 ),
                         status = cms.vint32(2)                               
                         )
process.b = cms.EDFilter("PdgIdAndStatusCandViewSelector",
                         src = cms.InputTag("genParticles"),
                         pdgId = cms.vint32( 5 ),
                         status = cms.vint32(2)                               
                         )
## light jets
process.lj = cms.EDFilter("PdgIdAndStatusCandViewSelector",
                          src = cms.InputTag("genParticles"),
                          pdgId = cms.vint32(1,2,3),
                          status = cms.vint32(2)                               
                          )
## pair of same flavor opoosite charge  light jets (lj) is written as 'jj' throughout whole file

process.jjCands = cms.EDProducer("CandViewShallowCloneCombiner",
                                 decay = cms.string('lj lj'),                                 
                                 cut = cms.string("mass > 0 && ((daughter(0).pdgId())!= (daughter(1).pdgId()))"),
                                 checkOverlap = cms.bool(False)
                                 )
## bBar cleaned
process.jjCandscleaned = cms.EDProducer("HZZ4LeptonsCandViewCleaner",
                                        srcObject = cms.InputTag("jjCands"),
                                        srcObjectsToRemove = cms.InputTag("jjCands"),
                                        module_label=cms.string("jjbar")                                              
                                        )
# make z+jets
process.lljjCands = cms.EDProducer("CandViewShallowCloneCombiner",
                                   decay = cms.string("Z jjCandscleaned"),
                                   cut = cms.string("0 < mass && (((daughter(0).pdgId!=1) || (daughter(0).pdgId!=-1)) && ((daughter(0).pdgId!=2) || (daughter(0).pdgId!=-2)) && ((daughter(0).pdgId!=3) || (daughter(0).pdgId!=-3))) "),
                                   )

# make bbBar
process.bBCands = cms.EDProducer("CandViewShallowCloneCombiner",
                                 decay = cms.string("b b"),
                                 cut = cms.string("mass > 0 && ((daughter(0).pdgId())!= (daughter(1).pdgId()))"),
                                 checkOverlap = cms.bool(False)
                                 )

# bBar cleaned
process.bBCandscleaned = cms.EDProducer("HZZ4LeptonsCandViewCleaner",
                                        srcObject = cms.InputTag("bBCands"),
                                        srcObjectsToRemove = cms.InputTag("bBCands"),
                                        module_label=cms.string("bbbar")                                              
                                        )

# make ZbBar
process.llbBCands = cms.EDProducer("CandViewShallowCloneCombiner",
                                   decay = cms.string("Z bBCandscleaned"),
                                   cut = cms.string("0 < mass && ((daughter(0).pdgId!=5) || (daughter(0).pdgId!=-5))"),
                                   )

# make ccBar
process.ccCands = cms.EDProducer("CandViewShallowCloneCombiner",
                                 decay = cms.string("c c"),
                                 cut = cms.string("mass > 0 && ((daughter(0).pdgId())!= (daughter(1).pdgId()))"),
                                 checkOverlap = cms.bool(False)
                                 )
# ccbar cleaned
process.ccCandscleaned = cms.EDProducer("HZZ4LeptonsCandViewCleaner",
                                        srcObject = cms.InputTag("ccCands"),
                                        srcObjectsToRemove = cms.InputTag("ccCands"),
                                        module_label=cms.string("ccbar")                                              
                                        )
# make Zccbar
process.llccCands = cms.EDProducer("CandViewShallowCloneCombiner",
                                   decay = cms.string("Z ccCandscleaned"),
                                   cut = cms.string("0 < mass && ((daughter(0).pdgId!=4) || (daughter(0).pdgId!=-4)) "),
                                   )

# leptons
process.leptons = cms.EDFilter("PdgIdAndStatusCandViewSelector",
                               src = cms.InputTag("genParticles"),
                               pdgId = cms.vint32(11,13),
                               status = cms.vint32(1)                     
                               )
# di-leptons
process.dileptons = cms.EDProducer("CandViewShallowCloneCombiner",
                                   decay = cms.string("leptons@+ leptons@-"),
                                   cut = cms.string("0< mass")
                                   )

## trileptons
process.trileptons = cms.EDProducer("CandViewShallowCloneCombiner",
                                    decay = cms.string("leptons@+ leptons@- leptons@+"),
                                    cut = cms.string("0 < mass")
                                    )
# trileptons cleaned
process.trileptonscleaned = cms.EDProducer("HZZ4LeptonsCandViewCleaner",
                                           srcObject = cms.InputTag("trileptons"),
                                           srcObjectsToRemove = cms.InputTag("trileptons"),
                                           module_label=cms.string("tri-leptons")                                              
                                           )

# fourleptons
process.fourleptons = cms.EDProducer("CandViewShallowCloneCombiner",
                                     decay = cms.string("trileptonscleaned@+ leptons@-"),
                                     cut = cms.string("0 < mass")
                                     )

# fourleptons cleaned
process.fourleptonscleaned = cms.EDProducer("HZZ4LeptonsCandViewCleaner",
                                            srcObject = cms.InputTag("fourleptons"),
                                            srcObjectsToRemove = cms.InputTag("fourleptons"),
                                            module_label=cms.string("four-leptons")                                              
                                            )

process.ZjetsGenAnalyzer = cms.EDAnalyzer("ZjetsGenAnalyzer",
                                          src = cms.InputTag('genParticles')
                                          )

process.p = cms.Path(process.Z * process.c * process.b * process.lj * process.jjCands * process.jjCandscleaned * process.lljjCands * process.bBCands * process.bBCandscleaned * process.llbBCands * process.ccCands * process.ccCandscleaned * process.llccCands * process.leptons * process.dileptons * process.trileptons * process.trileptonscleaned * process.fourleptons * process.fourleptonscleaned * process.hTozzTo4leptonsMCGenFilterNEW * process.ZjetsGenAnalyzer)

# Output definition
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('ZjetsGenAnalysisDY50.root')
                                   )

process.output = cms.OutputModule("PoolOutputModule",
                                  fileName = cms.untracked.string('ZjetsGenDY50.root'),
                                  dataset = cms.untracked.PSet(
    filterName = cms.untracked.string('Zjetsfilter')
    ),
                                  SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    )                               
                                  )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'file:F4BF2B2E-BE9C-E011-976E-001A4BA81FB8.root',
    'file:FE0F6638-C09C-E011-BCCA-001A4BD25578.root',
    'file:FE3430A0-BE9C-E011-AF2D-E0CB4E19F972.root'
    ),
                            )
process.out_step = cms.EndPath(process.output) # edm generation 
