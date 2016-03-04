import FWCore.ParameterSet.Config as cms

hTozzTo4leptonsRootTree = cms.EDAnalyzer("HZZ4LeptonsRootTree",
    decaychannel = cms.string('2e2mu'),
    rootFileName = cms.untracked.string('roottree_2e2mu.root'),
    useRECOformat = cms.untracked.bool(False),                                         

    module_to_search =  cms.untracked.vstring('a','b'),                                     
    par_to_search = cms.untracked.string('filename'),
    # PU
    fillPUinfo = cms.untracked.bool(True),
    PileupSrc  = cms.InputTag("addPileupInfo"),
                                               
    # HLT
    fillHLTinfo  = cms.untracked.bool(False),
    HLTInfoFired = cms.InputTag("hTozzTo4leptonsHLTInfo"),                                           
    HLTAnalysisinst = cms.string('hTozzTo4leptonsHLTAnalysis'),
    flagHLTnames=cms.VInputTag(cms.InputTag("flagHLTIsoMu11"), cms.InputTag("flagHLTMu15"),cms.InputTag("flagHLTDoubleMu3"),cms.InputTag("flagHLTIsoEle15L1I"),cms.InputTag("flagHLTIsoEle18L1R"),cms.InputTag("flagHLTDoubleIsoEle10L1I"), cms.InputTag("flagHLTDoubleIsoEle12L1R"), cms.InputTag("flagHLTaccept")),

    # Trigger matching                                           
    triggerEvent  = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    triggerFilter = cms.string('hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q'),
    triggerMatchObject   =  cms.InputTag("muonTriggerMatchHLT"),
    # triggerFilterAsym   =  cms.vstring('hltDiMuonL3PreFiltered8','hltDiMuonL3p5PreFiltered8'),
    # triggerFilterAsym   =  cms.vstring('hltMu40eta2p1CentralPFJet200MuCleaned','hltMu40eta2p1DiCentralPFJet50MuCleaned'
    triggerFilterAsym   =  cms.vstring('hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q'),
    triggerMatchObject_asym =  cms.InputTag("muonTriggerMatchHLTasym"),
                                               
    triggerHLTcollection =  cms.string("hltL3MuonCandidates"),                                              
    triggerMatchObjectEle =  cms.InputTag("electronTriggerMatchHLT"),                                           

    # Skim Early Data Higgs WG
    useSkimEarlyData=cms.untracked.bool(False),
    SkimEarlyDataAnalysisinst = cms.string('hTozzTo4leptonsSkimEarlyDataAnalysis'),
    flagSkimEarlyDatanames=cms.VInputTag(cms.InputTag("Flag_spikes"), cms.InputTag("Skim_highEnergyMuons"),cms.InputTag("Skim_highEnergyElectrons"),cms.InputTag("Skim_recoWMNfromPf"),cms.InputTag("Skim_recoWMNfromTc"),cms.InputTag("Skim_recoWENfromPf"),cms.InputTag("Skim_recoWENfromTc"),cms.InputTag("Skim_diMuonsJPsi"),cms.InputTag("Skim_diMuonsZ"),cms.InputTag("Skim_diElectronsZ"),cms.InputTag("Skim_triLeptonsMuMuMu"),cms.InputTag("Skim_triLeptonsMuMuEl"),cms.InputTag("Skim_triLeptonsMuElEl"),cms.InputTag("Skim_triLeptonsElElEl"),cms.InputTag("Skim_quadLeptons4Mu"),cms.InputTag("Skim_quadLeptons2Mu2El"),cms.InputTag("Skim_quadLeptons4El")),                                               

    # MC truth
    fillMCTruth = cms.untracked.bool(False),
    MCcollName = cms.InputTag("hTozzTo4leptonsMCDumper"),

    # RECO                                               
    RECOcollNameBest2e2mu = cms.VInputTag(cms.InputTag("hToZZTo4LeptonsBestCandidateOffselMother"), cms.InputTag("hToZZTo4LeptonsBestCandidateOffselBoson0"), cms.InputTag("hToZZTo4LeptonsBestCandidateOffselBoson1")),
    RECOcollNameBest4mu = cms.VInputTag(cms.InputTag("hToZZTo4LeptonsBestCandidateOffselMother"), cms.InputTag("hToZZTo4LeptonsBestCandidateOffselBoson0"), cms.InputTag("hToZZTo4LeptonsBestCandidateOffselBoson1")),
    RECOcollNameBest4e = cms.VInputTag(cms.InputTag("hToZZTo4LeptonsBestCandidateOffselMother"), cms.InputTag("hToZZTo4LeptonsBestCandidateOffselBoson0"), cms.InputTag("hToZZTo4LeptonsBestCandidateOffselBoson1")),

    # Additional RECO    
    useAdditionalRECO  = cms.untracked.bool(False),
    use2011EA          = cms.untracked.bool(True),                                     
    RECOcollNameZ      = cms.VInputTag(cms.InputTag("zToMuMu"), cms.InputTag("zToEE")),
    RECOcollNameZss    = cms.VInputTag(cms.InputTag("zToMuMussmerge"),cms.InputTag("zToEEssmerge"),cms.InputTag("zToCrossLeptons")),
    RECOcollNameDiLep  = cms.InputTag("dileptons"),                                                                                     
    RECOcollNameEEMM   = cms.VInputTag(cms.InputTag("hTozzTo4leptonsLooseIsol")),
    RECOcollNameMMMM   = cms.VInputTag(cms.InputTag("hTozzTo4leptonsMMMMLooseIsol")),
    RECOcollNameEEEE   = cms.VInputTag(cms.InputTag("hTozzTo4leptonsEEEELooseIsol")),
    RECOcollNameLLLLss = cms.VInputTag(cms.InputTag("quadLeptons4Mu"),cms.InputTag("quadLeptons2Mu2E"),cms.InputTag("quadLeptons4E")),
    RECOcollNameLLL    = cms.VInputTag(cms.InputTag("triLeptonsMuMuMu"),cms.InputTag("triLeptonsMuMuE"),cms.InputTag("triLeptonsMuEE"),cms.InputTag("triLeptonsEEE")),
    RECOcollNameLLLl   = cms.VInputTag(cms.InputTag("quadLeptons3Mu1E"),cms.InputTag("quadLeptons3E1Mu")),
    RECOcollNameLLLLssos = cms.VInputTag(cms.InputTag("quadLeptonsSSOSele"),cms.InputTag("quadLeptonsSSOSmu"),cms.InputTag("quadLeptonsSSOSelemu"),cms.InputTag("quadLeptonsSSOSmuele")),
    RECOcollNameLLLL   = cms.InputTag("allLLLL"),
                                               
    # isolation Tk, Ecal and Hcal
    SuperClustersLabel       = cms.InputTag("hTozzTo4leptonsMergedSuperClusters"),
    GsfTracksElectronsLabel  = cms.InputTag("electronGsfTracks"),
    ElectronsLabel           = cms.InputTag("hTozzTo4leptonsElectronIsolationProducer"),
    ElectronsMapLabel        = cms.InputTag("hTozzTo4leptonsHadIsolationProducer"),
    ElectronsTkMapLabel      = cms.InputTag("hTozzTo4leptonsHadIsolationProducer:Tk"),
    ElectronsEcalMapLabel    = cms.InputTag("hTozzTo4leptonsHadIsolationProducer:Ecal"),                                              
    ElectronsHcalMapLabel    = cms.InputTag("hTozzTo4leptonsHadIsolationProducer:Hcal"),
    ElectronsEgmLabel        = cms.InputTag("hTozzTo4leptonsElectronIsolationProducerEgamma"),
    ElectronsEgmTkMapLabel   = cms.InputTag("eleIsoFromDepsTkOptimized"),
    ElectronsEgmEcalMapLabel = cms.InputTag("eleIsoFromDepsEcalFromHitsByCrystalOptimized"),
    ElectronsEgmHcalMapLabel = cms.InputTag("eleIsoFromDepsHcalFromTowersOptimized"),

    MuonsLabel               = cms.InputTag("hTozzTo4leptonsMuonIsolationProducer"),
    MuonsMapLabel            = cms.InputTag("hTozzTo4leptonsMuonIsolationProducerOffsel"),
    MuonsTkMapLabel          = cms.InputTag("hTozzTo4leptonsMuonIsolationProducerOffsel:Tk"),                         
    MuonsEcalMapLabel        = cms.InputTag("hTozzTo4leptonsMuonIsolationProducerOffsel:Ecal"),
    MuonsHcalMapLabel        = cms.InputTag("hTozzTo4leptonsMuonIsolationProducerOffsel:Hcal"),

    # PF muons
    PFMuonsLabel             = cms.InputTag("hTozzTo4leptonsPFtoRECOMuon"),                                           

    # Particle Flow Isolation
    MuonPFIsoValueChargedAll    = cms.InputTag("muPFIsoValueChargedAll04PFBRECO"),
    MuonPFIsoValueCharged       = cms.InputTag("muPFIsoValueCharged04PFBRECO"),
    MuonPFIsoValueNeutral       = cms.InputTag("muPFIsoValueNeutral04PFBRECO"),
    MuonPFIsoValueGamma         = cms.InputTag("muPFIsoValueGamma04PFBRECO"),
    MuonPFIsoValuePU            = cms.InputTag("muPFIsoValuePU04PFBRECO"),


    ElectronPFIsoValueChargedAll = cms.InputTag("elPFIsoValueChargedAll04PFIdPFBRECO"),
    ElectronPFIsoValueCharged    = cms.InputTag("elPFIsoValueCharged04PFIdPFBRECO"),
    ElectronPFIsoValueNeutral    = cms.InputTag("elPFIsoValueNeutral04PFIdPFBRECO"),
    ElectronPFIsoValueGamma      = cms.InputTag("elPFIsoValueGamma04PFIdPFBRECO"),
    ElectronPFIsoValuePU         = cms.InputTag("elPFIsoValuePU04PFIdPFBRECO"),

    PFPhotonsLabel             = cms.InputTag("hTozzTo4leptonsPFfsrPhoton"),                                           
    PFpterrorLabel             = cms.InputTag("hTozzTo4leptonsPFfsrPhoton:ErrorMap"),

    PhotonPFIsoValueChargedAll = cms.InputTag("phPFIsoValueChargedAll03PFId"),
    PhotonPFIsoValueCharged    = cms.InputTag("phPFIsoValueCharged03PFId"),
    PhotonPFIsoValueNeutral    = cms.InputTag("phPFIsoValueNeutral03PFId"),
    PhotonPFIsoValueGamma      = cms.InputTag("phPFIsoValueGamma03PFId"),
    PhotonPFIsoValuePU         = cms.InputTag("phPFIsoValuePU03PFId"),
                                           

    
    # vertexing w.r.t primary vertex DA
    MuonsLabelVert           = cms.InputTag("hTozzTo4leptonsMuonSelector"),    
    MuonsMapLabelVert        = cms.InputTag("hTozzTo4leptonsIpToVtxProducer:VertexMuMap"),
    MuonsMapLabelVertValue   = cms.InputTag("hTozzTo4leptonsIpToVtxProducer:VertexValueMuMap"),
    MuonsMapLabelVertError   = cms.InputTag("hTozzTo4leptonsIpToVtxProducer:VertexErrorMuMap"),

    # vertexing w.r.t primary vertex KF
    MuonsMapLabelVertKF       = cms.InputTag("hTozzTo4leptonsIpToVtxProducerKF:VertexMuMap"),
    MuonsMapLabelVertValueKF  = cms.InputTag("hTozzTo4leptonsIpToVtxProducerKF:VertexValueMuMap"),
    MuonsMapLabelVertErrorKF  = cms.InputTag("hTozzTo4leptonsIpToVtxProducerKF:VertexErrorMuMap"),                                           

    # vertexing w.r.t GD, standard kalman and kinematic fit
    MuonsMapLabelVertGD       = cms.InputTag("hTozzTo4leptonsIpToVtxProducerGD:VertexMuMap"),
    MuonsMapLabelVertGDMMMM   = cms.InputTag("hTozzTo4leptonsIpToVtxProducerGDMMMM:VertexMuMap"),
    MuonsMapLabelVertStd      = cms.InputTag("hTozzTo4leptonsIpToVtxProducerStd:VertexMuMap"),
    MuonsMapLabelVertStdMMMM  = cms.InputTag("hTozzTo4leptonsIpToVtxProducerStdMMMM:VertexMuMap"),
    MuonsMapLabelVertKin      = cms.InputTag("hTozzTo4leptonsIpToVtxProducerKin:VertexMuMap"),
    MuonsMapLabelVertKinMMMM  = cms.InputTag("hTozzTo4leptonsIpToVtxProducerKinMMMM:VertexMuMap"),

    # vertexing w.r.t primary vertex DA                                     
    ElectronsLabelVert         = cms.InputTag("hTozzTo4leptonsElectronIsolationProducerEgamma"),
    ElectronsMapLabelVert      = cms.InputTag("hTozzTo4leptonsIpToVtxProducer:VertexEleMap"),
    ElectronsMapLabelVertValue = cms.InputTag("hTozzTo4leptonsIpToVtxProducer:VertexValueEleMap"), 
    ElectronsMapLabelVertError = cms.InputTag("hTozzTo4leptonsIpToVtxProducer:VertexErrorEleMap"),

    # vertexing w.r.t primary vertex KF
    ElectronsMapLabelVertKF      = cms.InputTag("hTozzTo4leptonsIpToVtxProducerKF:VertexEleMap"),
    ElectronsMapLabelVertValueKF = cms.InputTag("hTozzTo4leptonsIpToVtxProducerKF:VertexValueEleMap"), 
    ElectronsMapLabelVertErrorKF = cms.InputTag("hTozzTo4leptonsIpToVtxProducerKF:VertexErrorEleMap"),                                          
                                               

    # vertexing w.r.t GD, standard kalman and kinematic fit
    ElectronsMapLabelVertGD      = cms.InputTag("hTozzTo4leptonsIpToVtxProducerGD:VertexEleMap"),
    ElectronsMapLabelVertGDEEEE  = cms.InputTag("hTozzTo4leptonsIpToVtxProducerGDEEEE:VertexEleMap"),
    ElectronsMapLabelVertStd     = cms.InputTag("hTozzTo4leptonsIpToVtxProducerStd:VertexEleMap"),
    ElectronsMapLabelVertStdEEEE = cms.InputTag("hTozzTo4leptonsIpToVtxProducerStdEEEE:VertexEleMap"),
    ElectronsMapLabelVertKin     = cms.InputTag("hTozzTo4leptonsIpToVtxProducerKin:VertexEleMap"),
    ElectronsMapLabelVertKinEEEE = cms.InputTag("hTozzTo4leptonsIpToVtxProducerKinEEEE:VertexEleMap"),
                                               

    # vertexing with respect to primary vertex                                           
    MuonsSTIPMapLabelVert          = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:TipMuMap"),
    MuonsSLIPMapLabelVert          = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:LipMuMap"),
    MuonsSTIPMapLabelVertValue     = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:TipValueMuMap"),
    MuonsSLIPMapLabelVertValue     = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:LipValueMuMap"),
    MuonsSTIPMapLabelVertError     = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:TipErrorMuMap"),
    MuonsSLIPMapLabelVertError     = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:LipErrorMuMap"),
    	
    ElectronsSTIPMapLabelVert  = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:TipEleMap"),
    ElectronsSLIPMapLabelVert  = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:LipEleMap"),   	
    ElectronsSTIPMapLabelVertValue     = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:TipValueEleMap"),
    ElectronsSLIPMapLabelVertValue     = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:LipValueEleMap"),
    ElectronsSTIPMapLabelVertError     = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:TipErrorEleMap"),
    ElectronsSLIPMapLabelVertError     = cms.InputTag("hTozzTo4leptonsTipLipToVtxProducer:LipErrorEleMap"),

    # electron ID
    eleIDLabel                = cms.VInputTag(cms.InputTag("IsoleidClassLoose"),cms.InputTag("IsoleidClassMedium")),

    # electron regression
    eleRegressionEnergyErrorLabel  = cms.InputTag("eleRegressionEnergy:eneErrorRegForGsfEle"),
    eleRegressionEnergyLabel       = cms.InputTag("eleRegressionEnergy:eneRegForGsfEle"),          

    # CiC                                           
    eleID_VeryLooseTag        = cms.InputTag("eidVeryLoose"),
    eleID_LooseTag            = cms.InputTag("eidLoose"),
    eleID_MediumTag           = cms.InputTag("eidMedium"),
    eleID_TightTag            = cms.InputTag("eidTight"),                                           

    # Cic HZZ                                           
    eleID_HZZVeryLooseTag     = cms.InputTag("eidHZZVeryLoose"),
    eleID_HZZLooseTag         = cms.InputTag("eidHZZLoose"),
    eleID_HZZMediumTag        = cms.InputTag("eidHZZMedium"),
    eleID_HZZTightTag         = cms.InputTag("eidHZZHyperTight1"),

    # MVA ele ID BDT
    mvaTrigV0MapTag           = cms.InputTag("mvaTrigV025nsPHYS14"),
    mvaNonTrigV0MapTag        = cms.InputTag("mvaNonTrigV025nsPHYS14"),         

    # GD                                          
    ftsigmaVert               = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducer:ftsigma"),
    ftsigmalagVert            = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducer:ftsigmalag"),
    gdX_Vert                  = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducer:gdX"),
    gdY_Vert                  = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducer:gdY"),                                           
    gdZ_Vert                  = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducer:gdZ"),
    gdlagX_Vert               = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducer:gdlagX"),
    gdlagY_Vert               = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducer:gdlagY"),                                           
    gdlagZ_Vert               = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducer:gdlagZ"),
    gdlagProb_Vert            = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducer:gdlagProb"),
    gdlagNdof_Vert            = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducer:gdlagNdof"),                                           
    ftsigmaVertMMMM           = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerMMMM:ftsigma"),
    ftsigmalagVertMMMM        = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerMMMM:ftsigmalag"),
    gdX_VertMMMM              = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerMMMM:gdX"),
    gdY_VertMMMM              = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerMMMM:gdY"),                                           
    gdZ_VertMMMM              = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerMMMM:gdZ"),
    gdlagX_VertMMMM           = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerMMMM:gdlagX"),
    gdlagY_VertMMMM           = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerMMMM:gdlagY"),                                           
    gdlagZ_VertMMMM           = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerMMMM:gdlagZ"),
    gdlagProb_VertMMMM        = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerMMMM:gdlagProb"),
    gdlagNdof_VertMMMM        = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerMMMM:gdlagNdof"),                                             
    ftsigmaVertEEEE           = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerEEEE:ftsigma"),
    ftsigmalagVertEEEE        = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerEEEE:ftsigmalag"),
    gdX_VertEEEE              = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerEEEE:gdX"),
    gdY_VertEEEE              = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerEEEE:gdY"),                                           
    gdZ_VertEEEE              = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerEEEE:gdZ"),
    gdlagX_VertEEEE           = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerEEEE:gdlagX"),
    gdlagY_VertEEEE           = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerEEEE:gdlagY"),                                           
    gdlagZ_VertEEEE           = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerEEEE:gdlagZ"),
    gdlagProb_VertEEEE        = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerEEEE:gdlagProb"),
    gdlagNdof_VertEEEE        = cms.InputTag("hTozzTo4leptonsGeomDiscrimProducerEEEE:gdlagNdof"),                                             

    # ConstraintFit for 4l vertex                                        
    StandardFitVertex          = cms.InputTag("hTozzTo4leptonsConstraintFitProducer:StandardFitVertex"),
    StandardFitVertexMMMM      = cms.InputTag("hTozzTo4leptonsConstraintFitProducerMMMM:StandardFitVertex"),                                     
    StandardFitVertexEEEE      = cms.InputTag("hTozzTo4leptonsConstraintFitProducerEEEE:StandardFitVertex"),
    KinematicFitVertex         = cms.InputTag("hTozzTo4leptonsConstraintFitProducer:KinematicFitVertex"),
    KinematicFitVertexMMMM     = cms.InputTag("hTozzTo4leptonsConstraintFitProducerMMMM:KinematicFitVertex"),                                     
    KinematicFitVertexEEEE     = cms.InputTag("hTozzTo4leptonsConstraintFitProducerEEEE:KinematicFitVertex"),                                           
    RefittedMass               = cms.InputTag("hTozzTo4leptonsConstraintFitProducer:RefittedMass"),
    RefittedMassMMMM           = cms.InputTag("hTozzTo4leptonsConstraintFitProducerMMMM:RefittedMass"),
    RefittedMassEEEE           = cms.InputTag("hTozzTo4leptonsConstraintFitProducerEEEE:RefittedMass"),                                           

    # ConstraintFit for dilepton vertex
    StandardFitVertexDiLep     = cms.InputTag("hTozzTo4leptonsConstraintFitProducerDiLeptons:StandardFitVertex"),

    # ConstraintFit for 3l vertex                                               
    StandardFitVertexMMM       = cms.InputTag("hTozzTo4leptonsConstraintFitProducerTriLeptonsMMM:StandardFitVertex"),
    StandardFitVertexMME       = cms.InputTag("hTozzTo4leptonsConstraintFitProducerTriLeptonsMME:StandardFitVertex"),                                     
    StandardFitVertexEEE       = cms.InputTag("hTozzTo4leptonsConstraintFitProducerTriLeptonsEEE:StandardFitVertex"),
    StandardFitVertexMEE       = cms.InputTag("hTozzTo4leptonsConstraintFitProducerTriLeptonsMEE:StandardFitVertex"),                                              
                                               
                                               
    # CP
    MCCP_PhiLabel             = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsMCBestCandidateCPPhi"),
    MCCP_Phi1Label            = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsMCBestCandidateCPPhi1"),
    MCCP_Phi2Label            = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsMCBestCandidateCPPhi2"),
    MCCP_phi1RFLabel          = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsMCBestCandidateCPphi1RF"),
    MCCP_phi2RFLabel          = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsMCBestCandidateCPphi2RF"),
    MCCP_cosThetaStarLabel    = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsMCBestCandidateCPcosThetaStar"),
    MCCP_cosTheta1Label       = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsMCBestCandidateCPcosTheta1"),
    MCCP_cosTheta2Label       = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsMCBestCandidateCPcosTheta2"),
    MCCP_MELALabel            = cms.InputTag("hTozzTo4leptonsMCCP:hToZZTo4LeptonsMCBestCandidateCPMELA"),
 
    CP2e2mu_PhiLabel          = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsBestCandidateCPPhi"),
    CP2e2mu_Phi1Label         = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsBestCandidateCPPhi1"),
    CP2e2mu_Phi2Label         = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsBestCandidateCPPhi2"),
    CP2e2mu_phi1RFLabel       = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsBestCandidateCPphi1RF"),
    CP2e2mu_phi2RFLabel       = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsBestCandidateCPphi2RF"),
    CP2e2mu_cosThetaStarLabel = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsBestCandidateCPcosThetaStar"),
    CP2e2mu_cosTheta1Label    = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsBestCandidateCPcosTheta1"),
    CP2e2mu_cosTheta2Label    = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsBestCandidateCPcosTheta2"),
    CP2e2mu_MELALabel         = cms.InputTag("hTozzTo4leptonsCP:hToZZTo4LeptonsBestCandidateCPMELA"),
                                               
    CP4mu_PhiLabel            = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsBestCandidateCPPhi"),
    CP4mu_Phi1Label           = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsBestCandidateCPPhi1"),
    CP4mu_Phi2Label           = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsBestCandidateCPPhi2"),
    CP4mu_phi1RFLabel         = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsBestCandidateCPphi1RF"),
    CP4mu_phi2RFLabel         = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsBestCandidateCPphi2RF"),
    CP4mu_cosThetaStarLabel   = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsBestCandidateCPcosThetaStar"),
    CP4mu_cosTheta1Label      = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsBestCandidateCPcosTheta1"),
    CP4mu_cosTheta2Label      = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsBestCandidateCPcosTheta2"),
    CP4mu_MELALabel           = cms.InputTag("hTozzTo4leptonsCPMMMM:hToZZTo4LeptonsBestCandidateCPMELA"),
                                               
    CP4e_PhiLabel             = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsBestCandidateCPPhi"),
    CP4e_Phi1Label            = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsBestCandidateCPPhi1"),
    CP4e_Phi2Label            = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsBestCandidateCPPhi2"),
    CP4e_phi1RFLabel          = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsBestCandidateCPphi1RF"),
    CP4e_phi2RFLabel          = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsBestCandidateCPphi2RF"),
    CP4e_cosThetaStarLabel    = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsBestCandidateCPcosThetaStar"),
    CP4e_cosTheta1Label       = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsBestCandidateCPcosTheta1"),
    CP4e_cosTheta2Label       = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsBestCandidateCPcosTheta2"), 
    CP4e_MELALabel            = cms.InputTag("hTozzTo4leptonsCPEEEE:hToZZTo4LeptonsBestCandidateCPMELA"),


    # Other Objetcs
    PhotonsLabel       = cms.InputTag("photons"),
    TracksLabel        = cms.InputTag("generalTracks"),
    JetsLabel          = cms.InputTag("ak5PFJetsCorrection"),
    PuJetMvaMCfullDiscrLabel   = cms.InputTag("recoPuJetMvaMC:full53xDiscriminant"),
    PuJetMvaMCfullIdLabel      = cms.InputTag("recoPuJetMvaMC:full53xId"),                                              
    PuJetMvaDatafullDiscrLabel = cms.InputTag("recoPuJetMvaData:full53xDiscriminant"),
    PuJetMvaDatafullIdLabel    = cms.InputTag("recoPuJetMvaData:full53xId"),                                                             
    # RhoJetsLabel       = cms.InputTag("kt6corPFJets:rho"),
    RhoJetsLabel       = cms.InputTag("kt6PFJetsCentral:rho"),                                           
    VerticesLabel      = cms.InputTag("offlinePrimaryVertices"),

    # Gen MET
    GenMETLabel        = cms.InputTag("genMetTrue"),
    # Tracker MET                                           
    TrackerMETLabel    = cms.InputTag("tcMet"),                                           
    # Calo MET
    CaloMETLabel           = cms.InputTag("met"),                                           
    CaloMET_NoHFLabel      = cms.InputTag("metNoHF"),	

    useAdditionalMET       = cms.untracked.bool(False),
    CaloMET_HOLabel        = cms.InputTag("metHO"),
    CaloMET_OptLabel       = cms.InputTag("metOpt"), 
    CaloMET_OptNoHFLabel   = cms.InputTag("metOptNoHF"), 
    CaloMET_OptNoHFHOLabel = cms.InputTag("metOptNoHFHO"), 
    CaloMET_OptHOLabel     = cms.InputTag("metOptHO"), 
    CaloMET_NoHFHOLabel    = cms.InputTag("metNoHFHO"),
 
    # PF MET
    PfMETLabel             = cms.InputTag("pfMet"), 
    # HT MET                                          
    HtMET_IC5Label         = cms.InputTag("htMetIC5"), 
    HtMET_KT4Label         = cms.InputTag("htMetKT4"),
    HtMET_KT6Label         = cms.InputTag("htMetKT6"),
    #HtMET_SC5Label         = cms.InputTag("htMetSC5"),
    #HtMET_SC7Label         = cms.InputTag("htMetSC7"),
    HtMET_SC5Label         = cms.InputTag("htMetAK5"),
    HtMET_SC7Label         = cms.InputTag("htMetAK7"),                                           
    # JES correction                                           
    MET_JESCorIC5CaloJetLabel = cms.InputTag("metJESCorIC5CaloJet"),
    MET_JESCorKT4CaloJetLabel = cms.InputTag("metJESCorKT4CaloJet"),
    MET_JESCorKT6CaloJetLabel = cms.InputTag("metJESCorKT6CaloJet"),
    MET_JESCorSC5CaloJetLabel = cms.InputTag("metJESCorSC5CaloJet"),
    MET_JESCorSC7CaloJetLabel = cms.InputTag("metJESCorSC7CaloJet"),
    # MET correction for muons
    CorMETGlobalMuLabel       = cms.InputTag("corMetGlobalMuons"),

    # btagging                                           
    tCHighEff_bTagLabel  = cms.InputTag("trackCountingHighEffBJetTags"),
    tCHighPur_bTagLabel  = cms.InputTag("trackCountingHighPurBJetTags"),
    jPHighEff_bTagLabel  = cms.InputTag("jetProbabilityBJetTags"),
    jBP_bTagLabel        = cms.InputTag("jetBProbabilityBJetTags"),
    sSVHighEff_bTagLabel = cms.InputTag("simpleSecondaryVertexHighEffBJetTags"),
    sSVHighPur_bTagLabel = cms.InputTag("simpleSecondaryVertexHighPurBJetTags"),
    cSV_bTagLabel        = cms.InputTag("combinedSecondaryVertexBJetTags"),
    cSVMVA_bTagLabel     = cms.InputTag("combinedSecondaryVertexMVABJetTags"),
    sEByIP3d_bTagLabel   = cms.InputTag("softElectronByIP3dBJetTags"),
    sEByPt_bTagLabel     = cms.InputTag("softElectronByPtBJetTags"),
    sM_bTagLabel         = cms.InputTag("softMuonBJetTags"),
    sMByIP3d_bTagLabel   = cms.InputTag("softMuonByIP3dBJetTags"),
    sMByPt_bTagLabel     = cms.InputTag("softMuonByPtBJetTags"),
                                                

    # Conversion finder
    ConvMapDist          = cms.InputTag("ConvValueMapProd:dist"),                          
    ConvMapDcot          = cms.InputTag("ConvValueMapProd:dcot"),

    # VBF
    isVBF                = cms.bool(False)                               
)
