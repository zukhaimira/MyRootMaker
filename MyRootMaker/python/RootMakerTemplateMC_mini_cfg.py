import FWCore.ParameterSet.Config as cms
from L1Trigger.GlobalTrigger.gtDigis_cfi import *

process = cms.Process("ROOTMAKER")

# initialize MessageLogger and output report
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.categories.append('PATSummaryTables')

process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load('FastSimulation.HighLevelTrigger.HLTFastReco_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
    )
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string(
        'AC1B_out.root'
    )
)

#
#
#
#
#
## Sequences for AOD #####################################################################################
#if (!isMini):
#    # The Good Vertices collection ######################################################################
#    process.goodVertices = cms.EDFilter(
#        "VertexSelector",
#        filter = cms.bool(False),
#        src = cms.InputTag("offlinePrimaryVertices"),
#        cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
#    )
#    process.vertex_step = cms.Path(process.goodVertices)
#    
#    # The good primary vertex filter ####################################################################
#    process.primaryVertexFilter = cms.EDFilter(
#        "VertexSelector",
#        src = cms.InputTag("offlinePrimaryVertices"),
#        cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
#        filter = cms.bool(True)
#    )
#    # The beam scraping filter ##########################################################################
#    process.noscraping = cms.EDFilter(
#        "FilterOutScraping",
#        applyfilter = cms.untracked.bool(True),
#        debugOn = cms.untracked.bool(False),
#        numtrack = cms.untracked.uint32(10),
#        thresh = cms.untracked.double(0.25)
#    )
#    # The iso-based HBHE noise filter ###################################################################
#    process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')
#    # The CSC beam halo tight filter ####################################################################
#    process.load('RecoMET.METFilters.CSCTightHaloFilter_cfi')
#    # The HCAL laser filter #############################################################################
#    process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")
#    # The ECAL dead cell trigger primitive filter #######################################################
#    process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
#    # The EE bad SuperCrystal filter ####################################################################
#    process.load('RecoMET.METFilters.eeBadScFilter_cfi')
#    # The ECAL laser correction filter ##################################################################
#    process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')
#    # The tracking failure filter #######################################################################
#    process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
#    
#    process.filters_step = cms.Path(
#        process.primaryVertexFilter *
#        process.noscraping *
#        process.HBHENoiseFilter *
#        process.CSCTightHaloFilter *
#        process.hcalLaserEventFilter *
#        process.EcalDeadCellTriggerPrimitiveFilter *
#        process.trackingFailureFilter *
#        process.eeBadScFilter *
#        process.ecalLaserCorrFilter
#    )
#    
#    # Jet MET corrections ###############################################################################
#    process.load('RecoJets.Configuration.RecoPFJets_cff')
#    process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#    process.load('JetMETCorrections.Configuration.JetCorrectionServices_cff')
#    process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")
#    process.load("JetMETCorrections.Type1MET.correctedMet_cff")
#    process.load("PhysicsTools.PatAlgos.recoLayer0.metCorrections_cff")
#    from JetMETCorrections.Type1MET.correctedMet_cff import pfMetT1
#    
#    process.kt6PFJets.doRhoFastjet = True
#    process.kt6PFJets.doAreaFastjet = True
#    process.kt6PFJets.voronoiRfact = 0.9
#    
#    #process.pfJetMETcorr.jetCorrLabel = cms.string("ak4PFL1FastL2L3") #MC
#    
#    process.pfType0Type1CorrectedMet = pfMetT1.clone(
#        applyType0Corrections = cms.bool(False),
#        srcType1Corrections = cms.VInputTag(
#            cms.InputTag('pfMETcorrType0'),
#            cms.InputTag('pfJetMETcorr', 'type1')
#        )
#    )
#    process.metAnalysisSequence = cms.Sequence(
#        process.type0PFMEtCorrection*
#        process.patMETCorrections*
#        process.pfType0Type1CorrectedMet
#    )
#    process.jet_step = cms.Path(process.kt6PFJets*process.metAnalysisSequence)
#    
#    # PF Iso calculation for electrons ##################################################################
#    process.load('CommonTools.ParticleFlow.Isolation.pfElectronIsolation_cff')
#    from CommonTools.ParticleFlow.Isolation.pfElectronIsolation_cff import *
#    from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFPhotonIso
#    
#    process.stdElectronSequencePFIso = setupPFElectronIso(process, 'gedGsfElectrons')
#    process.stdPhotonSequencePFIso = setupPFPhotonIso(process, 'gedPhotons')
#    process.pfiso_step = cms.Path(
#        process.pfParticleSelectionSequence+ 
#        process.stdElectronSequencePFIso+ 
#        process.stdPhotonSequencePFIso
#    )
#    
#    # Electron ID #######################################################################################
#    process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_cfi")
#    process.electron_step = cms.Path(
#        process.eidHyperTight1MC*
#        process.eidLooseMC
#    )
#    
#    # Matching partons to jets ##########################################################################
#    process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")
#    process.AK4byRef.jets = cms.InputTag("ak4PFJets")
#    process.jetflavour_step = cms.Path(
#        process.myPartons* 
#        process.AK4Flavour
#    )
#    
#    # PAT ###############################################################################################
#    process.load("PhysicsTools.PatAlgos.patSequences_cff")
#    from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
#    process.out = cms.OutputModule("PoolOutputModule",
#        fileName = cms.untracked.string('patTuple.root'),
#        # save only events passing the full path
#        SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
#        # save PAT Layer 1 output; you need a '*' to
#        # unpack the list of commands 'patEventContent'
#        outputCommands = cms.untracked.vstring('drop *', *patEventContent )
#    )
#    from PhysicsTools.PatAlgos.tools.coreTools import *
#    from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT
#    postfix = "PFlow"
#    usePF2PAT(process,runPF2PAT=True,
#        jetAlgo='AK4', runOnMC=True, postfix=postfix,
#        jetCorrections=('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'],'Type-1'),
#        pvCollection=cms.InputTag('offlinePrimaryVertices')
#    )
#    process.pfPileUp.checkClosestZVertex = False
#    
#    ####getattr(process, "patPF2PATSequence"+postfix).remove(getattr(process, "cleanPatCandidates"+postfix))
#    ####getattr(process, "patPF2PATSequence"+postfix).remove(getattr(process, "countPatCandidates"+postfix))
#    ####process.patseq = cms.Sequence(getattr(process,"patPF2PATSequence"+postfix))
#    
#    #process.pat_step = cms.Path(process.patseq)
#    
#    # Pileup Jet ID (Hendrik) ###########################################################################
#    #https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetID#Running_the_algorithm_on_reco_je
#    #from CMGTools.External.pujetidsequence_cff import puJetId, puJetMva
#    #
#    #from RecoJets.JetProducers.puJetIDAlgo_cff import * 
#    #from RecoJets.JetProducers.puJetIDParams_cfi import * 
#    from RecoJets.JetProducers.PileupJetID_cfi import pileupJetIdProducer, _stdalgos_4x, _stdalgos_5x, _stdalgos, cutbased, _chsalgos_4x, _chsalgos_5x, _chsalgos 
#    from RecoJets.JetProducers.PileupJetID_cfi import * 
#    #process.load("RecoJets.JetProducers.PileupJetID_cfi")
#    #process.load("RecoJets.JetProducers.puJetIDAlgo_cff")
#    
#    process.recoPuJetId = pileupJetIdProducer.clone(
#        jets = cms.InputTag("ak4PFJets"),
#        applyJec = cms.bool(True),
#        inputIsCorrected = cms.bool(False), 
#        produceJetIds = cms.bool(True),
#        jetids = cms.InputTag(""),
#        runMvas = cms.bool(False),
#        vertexes = cms.InputTag("offlinePrimaryVertices"),
#        algos = cms.VPSet(cutbased)
#    )
#    
#    process.recoPuJetMva = pileupJetIdProducer.clone(
#        jets = cms.InputTag("ak4PFJets"),
#        jetids = cms.InputTag("recoPuJetId"),
#        applyJec = cms.bool(True),
#        produceJetIds = cms.bool(False),
#        runMvas = cms.bool(True),
#        vertexes = cms.InputTag("offlinePrimaryVertices"),
#        algos = cms.VPSet(_stdalgos),
#        inputIsCorrected = cms.bool(True),                                     
#    )
#    
#    
#    #from CMGTools.External.pujetidproducer_cfi import stdalgos_4x, stdalgos_5x, stdalgos, cutbased, chsalgos_4x, chsalgos_5x, chsalgos
#    #pileupJetIdProducer = cms.EDProducer('PileupJetIdProducer',
#    #                         produceJetIds = cms.bool(True),
#    #                         jetids = cms.InputTag(""),
#    #                         runMvas = cms.bool(True),
#    #                         jets = cms.InputTag("selectedPatJetsPFlow"),
#    #                         vertexes = cms.InputTag("offlinePrimaryVertices"),
#    #                         algos = cms.VPSet(stdalgos),
#    #                                     
#    #                         rho     = cms.InputTag("fixedGridRhoFastjetAll"),
#    #                         jec     = cms.string("AK4PF"),
#    #                         applyJec = cms.bool(False),
#    #                         inputIsCorrected = cms.bool(True),                                     
#    #                         residualsFromTxt = cms.bool(False),
#    #                       #  residualsTxt     = cms.FileInPath("RecoJets/JetProducers/data/dummy.txt"),
#    #)
#    #
#    #
#    #
#    #puJetId = pileupJetIdProducer.clone(
#    #    produceJetIds = cms.bool(True),
#    #    jetids = cms.InputTag(""),
#    #    runMvas = cms.bool(False),
#    #    jets = cms.InputTag("selectedPatJets"),
#    #    vertexes = cms.InputTag("offlinePrimaryVertices"),
#    #    algos = cms.VPSet(cutbased)
#    #)
#    #puJetMva = pileupJetIdProducer.clone(
#    #    produceJetIds = cms.bool(False),
#    #    jetids = cms.InputTag("puJetId"),
#    #    runMvas = cms.bool(True),
#    #    jets = cms.InputTag("selectedPatJets"),
#    #    vertexes = cms.InputTag("offlinePrimaryVertices"),
#    #    algos = stdalgos
#    #)
#    #
#    #process.recoPuJetId = puJetId.clone(
#    #    jets = cms.InputTag("ak4PFJets"),
#    #    applyJec = cms.bool(True),
#    #    inputIsCorrected = cms.bool(False),
#    #)
#    #process.recoPuJetMva = puJetMva.clone(
#    #    jets = cms.InputTag("ak4PFJets"),
#    #    jetids = cms.InputTag("recoPuJetId"),
#    #    applyJec = cms.bool(True),
#    #    inputIsCorrected = cms.bool(False),
#    #)
#    #
#    
#    
#    process.recoPuJetIdSequence = cms.Sequence(process.recoPuJetId * process.recoPuJetMva)
#    process.jetpuid_step = cms.Path(process.recoPuJetIdSequence)
#    
#












#process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))

# ROOTMAKER #########################################################################################
    process.makeroottree = cms.EDAnalyzer("RootMaker", 
    isMiniAOD = cms.untracked.bool(True),
    debug = cms.untracked.bool(False),

    # TRIGGER #####################################################
    #triggerPrescales = cms.InputTag("patTrigger"),
    #triggerObjects = cms.InputTag("selectedPatTrigger"),
    HLTriggerSelection = cms.untracked.vstring(),
    Trigger = cms.untracked.bool(True),

    # MUONS #######################################################
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
    RecMuonTrackIso = cms.untracked.double(1000000),
    RecMuonEtaMax = cms.untracked.double(2.5),
    RecMuonNum = cms.untracked.int32(1000),
    
    # ELECTRONS #######################################
    RecElectronHLTriggerMatching = cms.untracked.vstring(
        'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v.*:FilterTrue',
        'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v.*:FilterTrue',
        'HLT_Ele27_WP80_v.*',
        'HLT_Ele80_CaloIdVT_TrkIdT_v.*',
        'HLT_Ele80_CaloIdVT_GsfTrkIdT_v.*'
    ),
    RecElectron = cms.untracked.bool(True),
    RecElectronPtMin = cms.untracked.double(20.),
    RecElectronTrackIso = cms.untracked.double(1000000.),
    RecElectronEta = cms.untracked.double(2.5),
    RecElectronNum = cms.untracked.int32(100000),
    RecElectronFilterPtMin = cms.untracked.double(20.),

    # PHOTONS #######################################
    RecPhotonHLTriggerMatching = cms.untracked.vstring(
    ),
    RecPhoton = cms.untracked.bool(True),
    RecPhotonPtMin = cms.untracked.double(10.),
    RecPhotonEtaMax = cms.untracked.double(2.5),
    RecPhotonNum = cms.untracked.int32(100000),
    RecPhotonFilterPtMin = cms.untracked.double(10),
    
    # TAUS ########################################################
    RecTauHLTriggerMatching = cms.untracked.vstring(
    ),

    RecTau = cms.untracked.bool(False),
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
    RecJetHLTriggerMatching = cms.untracked.vstring(
    ),
    JetCorrection = cms.untracked.string('L1FastL2L3'),#MC

    RecAK4CALOJet = cms.untracked.bool(False),
    RecAK4JPTJet = cms.untracked.bool(False),
    RecAK4PFJet = cms.untracked.bool(False),
    RecAK4PFCHSJet = cms.untracked.bool(True),

    RecAK4PFCHSPtMin = cms.untracked.double(20.),
    RecAK4PFCHSEtaMax = cms.untracked.double(3.0),
    RecAK4PFCHSNum = cms.untracked.int32(100000),
    RecAK4PFCHSFilterPtMin = cms.untracked.double(20.),
    
    # GEN PARTICLES ###############################################
    GenAllParticles = cms.untracked.bool(False),
    GenSomeParticles = cms.untracked.bool(True),
    GenAK4Jets = cms.untracked.bool(True),

    # MET #########################################################
    RecPFMet = cms.untracked.bool(True),
    
    # TRACKS ######################################################
    RecTrack = cms.untracked.bool(False),
    RecTrackPtMin = cms.untracked.double(10.),
    RecTrackEtaMax = cms.untracked.double(2.5),
    RecTrackNum = cms.untracked.int32(100000),
    RecTrackFilterPtMin = cms.untracked.double(18.),
    
    # SUPERCLUSTER ################################################
    RecSuperCluster = cms.untracked.bool(True),
    RecSuperClusterFilterPtMin = cms.untracked.double(8.),
    RecSuperClusterBasicCluster = cms.untracked.bool(False),
    RecSuperClusterHit = cms.untracked.bool(False),

    # VERTICES ####################################################
    RecBeamSpot = cms.untracked.bool(True),
    RecPrimVertex = cms.untracked.bool(True),
    RecSecVertices = cms.untracked.bool(False),

    # INPUT TAGS ##################################################
    triggerBits = cms.InputTag("TriggerResults","","HLT"),
    metFilterBits = cms.InputTag("TriggerResults", "", "PAT"),
    lostTracks = cms.InputTag("lostTracks", "", "PAT"), # mini only
    generalTracks = cms.InputTag("generalTracks"), # AOD only
    unpackedTracks = cms.InputTag("unpackedTracksAndVertices"), # mini only
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    secondaryVertices = cms.InputTag("slimmedSecondaryVertices", "", "PAT"),
    beamSpot = cms.InputTag("offlineBeamSpot", "", "RECO"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    photons = cms.InputTag("slimmedPhotons"),
    taus = cms.InputTag("slimmedTaus"),
    ak4calojets = cms.InputTag("patJetsAK4Calo"), # AOD only
    ak4jptjets = cms.InputTag("patJetsAK4JPT"), # AOD only
    ak4pfjets = cms.InputTag("ak4PFJets"), # AOD only
    ak4pfchsjets = cms.InputTag("slimmedJets"),
    jetsAK8 = cms.InputTag("slimmedJetsAK8"), # mini only
    mets = cms.InputTag("slimmedMETs"), # mini only
    pfmet = cms.InputTag("pfMet"), # AOD only
    pfmett1 = cms.InputTag("pfMetT1"), # AOD only
    pfmett1t0 = cms.InputTag("pfType0Type1CorrectedMet"), # AOD only
    packedGenParticles = cms.InputTag("packedGenParticles"), # mini only
    genParticles = cms.InputTag("prunedGenParticles"), 
    genJets = cms.InputTag("slimmedGenJets", "", "PAT"),
    genInfo = cms.InputTag("generator", "", "SIM"),
    puInfo = cms.InputTag("addPileupInfo", "", "HLT"),
    PfCands = cms.InputTag("PFCandidates"), # AOD only
    packedPfCands = cms.InputTag("packedPFCandidates"), # mini only
    barrelHits = cms.InputTag("reducedEgamma", "reducedEBRecHits", "PAT"),
    endcapHits = cms.InputTag("reducedEgamma", "reducedEERecHits", "PAT"),
    esHits = cms.InputTag("reducedEgamma", "reducedESRecHits", "PAT"),
    superClusters = cms.InputTag("reducedEgamma", "reducedSuperClusters", "PAT"),
    ebeeClusters = cms.InputTag("reducedEgamma", "reducedEBEEClusters", "PAT"),
    esClusters = cms.InputTag("reducedEgamma", "reducedESClusters", "PAT"),
    conversions = cms.InputTag("reducedEgamma", "reducedConversions", "PAT"),
    gedGsfElectronCores = cms.InputTag("reducedEgamma", "reducedGedGsfElectronCores", "PAT"), # mini only
    gedPhotonCores = cms.InputTag("reducedEgamma", "reducedGedPhotonCores", "PAT"), # mini only

)

#######Electron ID
process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_cfi")
#process.electron_step = cms.Path(process.eidHyperTight1MC*process.eidLooseMC)

#######Tracks
process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')
process.load('RecoJets.JetAssociationProducers.ak4JTA_cff')
process.ak4JetTracksAssociatorAtVertexPF.jets = cms.InputTag("ak4PFJetsCHS")
process.ak4JetTracksAssociatorAtVertexPF.tracks = cms.InputTag("unpackedTracksAndVertices")

#process.makeroottree.isMiniAOD = cms.untracked.bool(isMini)
process.p = cms.Path(
    process.unpackedTracksAndVertices *
#    process.eidHyperTight1MC*
#    process.eidLooseMC*
    process.makeroottree
)

