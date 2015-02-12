import FWCore.ParameterSet.Config as cms
#from MyRootMaker.MyRootMaker.RootMakerTemplateMC_mini_cfg import *
from L1Trigger.GlobalTrigger.gtDigis_cfi import *

process = cms.Process("ROOTMAKER")

process.load('FastSimulation.HighLevelTrigger.HLTFastReco_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) ) # 13000

process.MessageLogger.cerr.FwkReport.reportEvery = 100
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
        'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/MINIAODSIM/PU40bx25_PHYS14_25_V1-v1/00000/C20B68E7-0277-E411-85E5-001E67396A22.root' # 35400 events
    )
)
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions
process.GlobalTag.globaltag = cms.string('PHYS14_25_V1::All')

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('AC1B_test.root')
)

#process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))

# ROOTMAKER #########################################################################################
process.makeroottree = cms.EDAnalyzer("RootMaker",
    
    isMiniAOD = cms.untracked.bool(True),
    debug = cms.untracked.bool(False),

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
#    RecMuonNum = cms.untracked.int32(1000),
    
# ELECTRONS and PHOTONS #######################################

    electrons = cms.InputTag("slimmedElectrons"),
    photons = cms.InputTag("slimmedPhotons"),

    barrelHits = cms.InputTag("reducedEgamma", "reducedEBRecHits", "PAT"),
    endcapHits = cms.InputTag("reducedEgamma", "reducedEERecHits", "PAT"),
    esHits = cms.InputTag("reducedEgamma", "reducedESRecHits", "PAT"),
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
#    RecElectronNum = cms.untracked.int32(100000),
    RecElectronFilterPtMin = cms.untracked.double(20.),

    RecPhotonHLTriggerMatching = cms.untracked.vstring(
    ),
    RecPhoton = cms.untracked.bool(True),
    RecPhotonPtMin = cms.untracked.double(10.),
    RecPhotonEtaMax = cms.untracked.double(2.5),
#    RecPhotonNum = cms.untracked.int32(100000),
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
#    RecTauNum = cms.untracked.int32(100000),

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
    
#    RecAK4CaloJet = cms.untracked.bool(True),
#    RecAK4CaloPtMin = cms.untracked.double(20.),
#    RecAK4CaloEtaMax = cms.untracked.double(2.4),
#    RecAK4CaloNum = cms.untracked.int32(100000),
#    RecAK4CaloFilterPtMin = cms.untracked.double(20.),
    
#    RecAK4JPTJet = cms.untracked.bool(True),
#    RecAK4JPTPtMin = cms.untracked.double(20.),
#    RecAK4JPTEtaMax = cms.untracked.double(2.4),
#    RecAK4JPTNum = cms.untracked.int32(100000),
#    RecAK4JPTFilterPtMin = cms.untracked.double(20.),
    
#    RecAK4PFJet = cms.untracked.bool(True),
#    RecAK4PFPtMin = cms.untracked.double(20.),
#    RecAK4PFEtaMax = cms.untracked.double(3.0),
#    RecAK4PFNum = cms.untracked.int32(100000),
#    RecAK4PFFilterPtMin = cms.untracked.double(20.),
    
    RecAK4PFCHSJet = cms.untracked.bool(True),
    RecAK4PFCHSPtMin = cms.untracked.double(20.),
    RecAK4PFCHSEtaMax = cms.untracked.double(3.0),
#    RecAK4PFCHSNum = cms.untracked.int32(100000),
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

    RecTrack = cms.untracked.bool(False),
    RecTrackPtMin = cms.untracked.double(10.),
    RecTrackEtaMax = cms.untracked.double(2.5),
#    RecTrackNum = cms.untracked.int32(100000),
    RecTrackFilterPtMin = cms.untracked.double(18.),
    
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
process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_cfi")
#process.electron_step = cms.Path(process.eidHyperTight1MC*process.eidLooseMC)
process.p = cms.Path(
#    process.eidHyperTight1MC*
#    process.eidLooseMC*
    process.makeroottree
)

