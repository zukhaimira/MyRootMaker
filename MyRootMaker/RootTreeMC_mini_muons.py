import FWCore.ParameterSet.Config as cms
#from MyRootMaker.MyRootMaker.RootMakerTemplateMC_mini_cfg import *
#from L1Trigger.GlobalTrigger.gtDigis_cfi import *

process = cms.Process("ROOTMAKER")

process.load('FastSimulation.HighLevelTrigger.HLTFastReco_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(3000) )

process.MessageLogger.cerr.FwkReport.reportEvery = 10
#process.MessageLogger.categories.append('PATSummaryTables')
#process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
#process.MessageLogger.categories.extend(["GetManyWithoutRegistration","GetByLabelWithoutRegistration"])
#_messageSettings = cms.untracked.PSet(
#                reportEvery = cms.untracked.int32(1),
#                            optionalPSet = cms.untracked.bool(True),
#                            limit = cms.untracked.int32(10000000)
#                        )
#process.MessageLogger.cerr.GetManyWithoutRegistration = _messageSettings
#process.MessageLogger.cerr.GetByLabelWithoutRegistration = _messageSettings

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/ZZTo4L_Tune4C_13TeV-powheg-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/241933C6-E269-E411-B532-7845C4FC3A19.root'
    )
)
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions
process.GlobalTag.globaltag = cms.string('POSTLS172_V3::All')

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('AC1B_test.root')
)

process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))

# ROOTMAKER #########################################################################################
process.makeroottree = cms.EDAnalyzer("RootMakerMini",
    
    packedPfCands = cms.InputTag("packedPFCandidates"),
    hcalNoiseInfo = cms.InputTag("hcalnoise", "", "RECO"),

# TRIGGER #####################################################
    triggerBits = cms.InputTag("TriggerResults","","HLT"),
    triggerPrescales = cms.InputTag("patTrigger"),
    triggerObjects = cms.InputTag("selectedPatTrigger"),
    puInfo = cms.InputTag("addPileupInfo", "", "HLT"),

    metFilterBits = cms.InputTag("TriggerResults", "", "PAT"),

    HLTriggerSelection = cms.untracked.vstring(),
    Trigger = cms.untracked.bool(False),

# MUONS #######################################################
    muons = cms.InputTag("slimmedMuons"),

    RecMuonHLTriggerMatching = cms.untracked.vstring(
        'HLT_Mu17_Mu8_v.*:FilterTrue',
        'HLT_Mu17_TkMu8_v.*:FilterTrue',
        'HLT_Mu22_TkMu22_v.*:FilterTrue',
        'HLT_Mu(8|17)_v.*',
        'HLT_IsoMu(24|30)_v.*',
        'HLT_IsoMu(24|30)_eta2p1_v.*',
        'HLT_Mu40_v.*',
        'HLT_Mu50_eta2p1_v.*',
        'HLT_Mu40_eta2p1_v.*'
    ),

    RecMuon = cms.untracked.bool(True),
    RecMuonPtMin = cms.untracked.double(20),
#    RecMuonTrackIso = cms.untracked.double(1000000),
    RecMuonEtaMax = cms.untracked.double(2.5),
    RecMuonNum = cms.untracked.int32(1000),
    
# ELECTRONS and PHOTONS #######################################

    electrons = cms.InputTag("slimmedElectrons"),
    photons = cms.InputTag("slimmedPhotons"),

    ebRecHits = cms.InputTag("reducedEgamma", "reducedEBRecHits", "PAT"),
    eeRecHits = cms.InputTag("reducedEgamma", "reducedEERecHits", "PAT"),
    esRecHits = cms.InputTag("reducedEgamma", "reducedESRecHits", "PAT"),
    ebeeClusters = cms.InputTag("reducedEgamma", "reducedEBEEClusters", "PAT"),
    esClusters = cms.InputTag("reducedEgamma", "reducedESClusters", "PAT"),
    conversions = cms.InputTag("reducedEgamma", "reducedConversions", "PAT"),
    singleLegConversions = cms.InputTag("reducedEgamma", "reducedSingleLegConversions", "PAT"),
    gedGsfElectronCores = cms.InputTag("reducedEgamma", "reducedGedGsfElectronCores", "PAT"),
    gedPhotonCores = cms.InputTag("reducedEgamma", "reducedGedPhotonCores", "PAT"),
    superClusters = cms.InputTag("reducedEgamma", "reducedSuperClusters", "PAT"),

    RecElectronHLTriggerMatching = cms.untracked.vstring(
        'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v.*:FilterTrue',
        'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v.*:FilterTrue',
        'HLT_Ele27_WP80_v.*',
        'HLT_Ele80_CaloIdVT_TrkIdT_v.*',
        'HLT_Ele80_CaloIdVT_GsfTrkIdT_v.*'
    ),

    RecElectron = cms.untracked.bool(True),
    RecElectronPtMin = cms.untracked.double(20.),
#    RecElectronTrackIso = cms.untracked.double(1000000.),
    RecElectronEta = cms.untracked.double(2.5),
    RecElectronNum = cms.untracked.int32(100000),
    RecElectronFilterPtMin = cms.untracked.double(20.),

    RecPhotonHLTriggerMatching = cms.untracked.vstring(
    ),
    RecPhoton = cms.untracked.bool(True),
    RecPhotonPtMin = cms.untracked.double(10.),
    RecPhotonEtaMax = cms.untracked.double(2.5),
    RecPhotonNum = cms.untracked.int32(100000),
    RecPhotonFilterPtMin = cms.untracked.double(10),
    
# TAUS ########################################################
    taus = cms.InputTag("slimmedTaus"),

    RecTauHLTriggerMatching = cms.untracked.vstring(
    ),

    RecTau = cms.untracked.bool(True),
    RecTauDiscriminators = cms.untracked.vstring(
        'hpsPFTauDiscriminationByDecayModeFinding',
        'hpsPFTauDiscriminationByLooseElectronRejection',
        'hpsPFTauDiscriminationByLooseIsolation',
        'hpsPFTauDiscriminationByLooseMuonRejection',
        'hpsPFTauDiscriminationByMediumElectronRejection',
        'hpsPFTauDiscriminationByMediumIsolation',
        'hpsPFTauDiscriminationByTightElectronRejection',
        'hpsPFTauDiscriminationByTightIsolation',
        'hpsPFTauDiscriminationByTightMuonRejection',
        'hpsPFTauDiscriminationByVLooseIsolation'
    ),
    RecTauPtMin = cms.untracked.double(0.),
    RecTauEta = cms.untracked.double(0.),
    RecTauNum = cms.untracked.int32(100000),

# JETS ########################################################
    jets = cms.InputTag("slimmedJets"),
    jetsAK8 = cms.InputTag("slimmedJetsAK8"),

    rhoAll = cms.InputTag("fixedGridRhoAll", "", "RECO"),
    rhoFastjetAll = cms.InputTag("fixedGridRhoFastjetAll", "", "RECO"),
    rhoFastjetAllCalo = cms.InputTag("fixedGridRhoFastjetAllCalo", "", "RECO"),
    rhoFastjetCentralCalo = cms.InputTag("fixedGridRhoFastjetCentralCalo", "", "RECO"),
    rhoFastjetCentralChargedPileUp = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp", "", "RECO"),
    rhoFastjetCentralNeutral = cms.InputTag("fixedGridRhoFastjetCentralNeutral", "", "RECO"),

    RecJetHLTriggerMatching = cms.untracked.vstring(
    ),

    JetCorrection = cms.untracked.string('L1FastL2L3'),#MC
    
    RecAK4CaloJet = cms.untracked.bool(True),
    RecAK4CaloPtMin = cms.untracked.double(20.),
    RecAK4CaloEtaMax = cms.untracked.double(2.4),
    RecAK4CaloNum = cms.untracked.int32(100000),
    RecAK4CaloFilterPtMin = cms.untracked.double(20.),
    
    RecAK4JPTJet = cms.untracked.bool(True),
    RecAK4JPTPtMin = cms.untracked.double(20.),
    RecAK4JPTEtaMax = cms.untracked.double(2.4),
    RecAK4JPTNum = cms.untracked.int32(100000),
    RecAK4JPTFilterPtMin = cms.untracked.double(20.),
    
    RecAK4PFJet = cms.untracked.bool(True),
    RecAK4PFPtMin = cms.untracked.double(20.),
    RecAK4PFEtaMax = cms.untracked.double(3.0),
    RecAK4PFNum = cms.untracked.int32(100000),
    RecAK4PFFilterPtMin = cms.untracked.double(20.),
    
    RecAK4PFCHSJet = cms.untracked.bool(False),
    RecAK4PFCHSPtMin = cms.untracked.double(20.),
    RecAK4PFCHSEtaMax = cms.untracked.double(3.0),
    RecAK4PFCHSNum = cms.untracked.int32(100000),
    RecAK4PFCHSFilterPtMin = cms.untracked.double(20.),
    
# GEN PARTICLES ###############################################
    packedGenParticles = cms.InputTag("packedGenParticles"),
    prunedGenParticles = cms.InputTag("prunedGenParticles"),
    genJets = cms.InputTag("slimmedGenJets", "", "PAT"),

    genInfo = cms.InputTag("generator", "", "SIM"),

    GenAllParticles = cms.untracked.bool(False),
    GenSomeParticles = cms.untracked.bool(True),
    GenAK4Jets = cms.untracked.bool(True),

# MET #########################################################
    mets = cms.InputTag("slimmedMETs"),

    RecPFMet = cms.untracked.bool(True),
    
# TRACKS ######################################################
    lostTracks = cms.InputTag("lostTracks", "", "PAT"),
    #generalTracks = cms.InputTag("generalTracks", "", "RECO"),

#    RecTrack = cms.untracked.bool(False),
#    RecTrackPtMin = cms.untracked.double(10.),
#    RecTrackEtaMax = cms.untracked.double(2.5),
#    RecTrackNum = cms.untracked.int32(100000),
#    RecTrackFilterPtMin = cms.untracked.double(18.),
    
# SUPERCLUSTER ################################################
    RecSuperCluster = cms.untracked.bool(True),
    RecSuperClusterFilterPtMin = cms.untracked.double(8.),
    RecSuperClusterBasicCluster = cms.untracked.bool(False),
    RecSuperClusterHit = cms.untracked.bool(False),

# VERTICES ####################################################
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    secondaryVertices = cms.InputTag("slimmedSecondaryVertices", "", "PAT"),

    beamSpot = cms.InputTag("offlineBeamSpot", "", "RECO"),

    RecBeamSpot = cms.untracked.bool(True),
    RecPrimVertex = cms.untracked.bool(True),
    RecSecVertices = cms.untracked.bool(True),
    RecVertexTRKChi2 = cms.untracked.double(5),
    RecVertexTRKHitsMin = cms.untracked.int32(6),
    RecVertexChi2 = cms.untracked.double(3),
    RecVertexSig2D = cms.untracked.double(15),
    RecKaonMasswin = cms.untracked.double(0.05),
    RecLambdaMasswin = cms.untracked.double(0.02)

)

# examples of things that go in the short version
#process.makeroottree.RecMuonNum = cms.untracked.int32(0)
#process.makeroottree.HLTriggerSelection = cms.untracked.vstring()
#process.makeroottree.GenAllParticles = cms.untracked.bool(False)
#process.makeroottree.GenSomeParticles = cms.untracked.bool(False)
#process.makeroottree.GenAK4Jets = cms.untracked.bool(False)


#######Electron ID
#process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_cfi")
#process.electron_step = cms.Path(process.eidHyperTight1MC*process.eidLooseMC)
process.p = cms.Path(
    #process.eidHyperTight1MC*
    #process.eidLooseMC*
    process.makeroottree
)

##################################################################################################
#from MyRootMaker.MyRootMaker.RootMakerTemplateMC_cfg import *
#
##process.GlobalTag.toGet = cms.VPSet(
##		cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
##			tag = cms.string("TrackProbabilityCalibration_2D_MC53X_v2"),
##			connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
##		cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
##			tag = cms.string("TrackProbabilityCalibration_3D_MC53X_v2"),
##			connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
##		)
#
#process.load('RecoBTag/Configuration/RecoBTag_cff')
#process.btag = cms.Path(process.btagging)
#
#process.makeroottree.RecMuonNum = cms.untracked.int32(0)
#process.makeroottree.HLTriggerSelection = cms.untracked.vstring()
#process.patJetCorrFactors.levels=cms.vstring('L1FastJet','L2Relative', 'L3Absolute')
#process.makeroottree.GenAllParticles = cms.untracked.bool(True)
#process.makeroottree.GenSomeParticles = cms.untracked.bool(False)
#process.makeroottree.GenAK5Jets = cms.untracked.bool(True)
#
################################################################################
################################################################################
################################################################################

#process.MessageLogger.cerr.threshold = 'INFO'
#process.GlobalTag.globaltag = cms.string('START70_V7::All')
#
## The Good vertices collection _____________________________________________||
#process.goodVertices = cms.EDFilter(
#                "VertexSelector",
#                filter = cms.bool(False),
#                src = cms.InputTag("offlinePrimaryVertices"),
#                cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
#                )
#process.vertex_step=cms.Path(process.goodVertices)
######Filter
### The good primary vertex filter ____________________________________________||
#process.primaryVertexFilter = cms.EDFilter(
#                "VertexSelector",
#                src = cms.InputTag("offlinePrimaryVertices"),
#                cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
#                filter = cms.bool(True)
#                )
### The beam scraping filter __________________________________________________||
#process.noscraping = cms.EDFilter(
#                "FilterOutScraping",
#                applyfilter = cms.untracked.bool(True),
#                debugOn = cms.untracked.bool(False),
#                numtrack = cms.untracked.uint32(10),
#                thresh = cms.untracked.double(0.25)
#                )
### The iso-based HBHE noise filter ___________________________________________||
#process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')
### The CSC beam halo tight filter ____________________________________________||
#process.load('RecoMET.METFilters.CSCTightHaloFilter_cfi')
### The HCAL laser filter _____________________________________________________||
#process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")
### The ECAL dead cell trigger primitive filter _______________________________||
#process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
### The EE bad SuperCrystal filter ____________________________________________||
#process.load('RecoMET.METFilters.eeBadScFilter_cfi')
### The ECAL laser correction filter
#process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')
### The tracking failure filter _______________________________________________||
#process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
#
#process.filters_step = cms.Path(
#                process.primaryVertexFilter *
#                process.noscraping *
#                process.HBHENoiseFilter *
#                process.CSCTightHaloFilter *
#                process.hcalLaserEventFilter *
#                process.EcalDeadCellTriggerPrimitiveFilter *
#                process.trackingFailureFilter *
#                process.eeBadScFilter *
#                process.ecalLaserCorrFilter
#                )
#
######## Jet MET corrections
#process.load('RecoJets.Configuration.RecoPFJets_cff')
#process.kt6PFJets.doRhoFastjet = True
#process.kt6PFJets.doAreaFastjet = True
#process.kt6PFJets.voronoiRfact = 0.9
#
#process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#process.load('JetMETCorrections.Configuration.JetCorrectionServices_cff')
#process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")
##process.load('JetMETCorrections.Type1MET.pfMETCorrections_cff')
##process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3") #MC
##process.pfJetMETcorr.jetCorrLabel = cms.string("ak4PFL1FastL2L3") #MC
#
##from JetMETCorrections.Type1MET.pfMETCorrections_cff import pfType1CorrectedMet
##process.pfType0Type1CorrectedMet = pfType1CorrectedMet.clone(
##applyType0Corrections = cms.bool(False),
##srcType1Corrections = cms.VInputTag(
##    cms.InputTag('pfMETcorrType0'),
##    cms.InputTag('pfJetMETcorr', 'type1')
##)
##)
##process.metAnalysisSequence=cms.Sequence(process.type0PFMEtCorrection*process.producePFMETCorrections*process.pfType0Type1CorrectedMet)
##process.jet_step = cms.Path(process.kt6PFJets*process.metAnalysisSequence)
#
#######PF ISO calculation for Electrons
#from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFPhotonIso
##from CommonTools.ParticleFlow.Isolation.pfElectronIsolation_cff import *
#
## NEW STUFF ##################
#process.stdElectronSequencePFIso = setupPFElectronIso(process, 'gedGsfElectrons')#'gsfElectrons')
#process.stdPhotonSequencePFIso = setupPFPhotonIso(process, 'photons')
#process.pfiso_step = cms.Path( process.pfParticleSelectionSequence +
#                               process.stdElectronSequencePFIso +
#                               process.stdPhotonSequencePFIso)
###############################
#
## OLD STUFF #########################
##process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
##process.phoIsoSequence = setupPFPhotonIso(process, 'photons')
##process.pfiso_step = cms.Path( process.pfParticleSelectionSequence + process.eleIsoSequence + process.phoIsoSequence)
######################################
#
#
#######Matching Partons to Jets
#process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")
##process.AK5byRef.jets = cms.InputTag("ak5PFJets")
##process.jetflavour_step = cms.Path(process.myPartons * process.AK5Flavour)
#process.AK4byRef.jets = cms.InputTag("slimmedJets")
#process.jetflavour_step = cms.Path(process.myPartons * process.AK4Flavour)
#
#######PAT
#process.load("PhysicsTools.PatAlgos.patSequences_cff")
#from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
#process.out = cms.OutputModule("PoolOutputModule",
#		fileName = cms.untracked.string('patTuple.root'),
## save only events passing the full path
#		SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
## save PAT Layer 1 output; you need a '*' to
## unpack the list of commands 'patEventContent'
#		outputCommands = cms.untracked.vstring('drop *', *patEventContent )
#		)
##from PhysicsTools.PatAlgos.tools.coreTools import *
##removeAllPATObjectsBut(process, ['Jets'])
#from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT
#postfix = "PFlow"
#usePF2PAT(process,runPF2PAT=True,
#		jetAlgo='AK5', runOnMC=False, postfix=postfix,
#                jetCorrections=('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'],'Type-1'),
#		pvCollection=cms.InputTag('offlinePrimaryVertices')
#	 )
#process.pfPileUpPFlow.checkClosestZVertex = False
#getattr(process, "patPF2PATSequence"+postfix).remove(getattr(process, "cleanPatCandidates"+postfix))
#getattr(process, "patPF2PATSequence"+postfix).remove(getattr(process, "countPatCandidates"+postfix))
#
#process.patseq = cms.Sequence(getattr(process,"patPF2PATSequence"+postfix))
#process.pat_step = cms.Path(process.patseq)
#
#######Pileup Jet ID (Hendrik)
##https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetID#Running_the_algorithm_on_reco_je
#from CMGTools.External.pujetidsequence_cff import puJetId, puJetMva
#
#process.recoPuJetId = puJetId.clone(
#		jets = cms.InputTag("ak5PFJets"),
#		applyJec = cms.bool(True),
#		inputIsCorrected = cms.bool(False),                
#		)
#
#process.recoPuJetMva = puJetMva.clone(
#		jets = cms.InputTag("ak5PFJets"),
#		jetids = cms.InputTag("recoPuJetId"),
#		applyJec = cms.bool(True),
#		inputIsCorrected = cms.bool(False),                
#		)
#process.recoPuJetIdSqeuence = cms.Sequence(process.recoPuJetId * process.recoPuJetMva)
#process.jetpuid_step = cms.Path(process.recoPuJetIdSqeuence)
#
##process.genPlusSimParticles = cms.EDProducer("GenPlusSimParticleProducer",
##src = cms.InputTag("g4SimHits"),
##setStatus = cms.int32(8),
##filter = cms.vstring("pt > 0.0"),
##genParticles = cms.InputTag("genParticles")
##)
#
#process.roottree_step = cms.EndPath(process.makeroottree)
#
