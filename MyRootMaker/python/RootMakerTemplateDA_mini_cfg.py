import FWCore.ParameterSet.Config as cms
from L1Trigger.GlobalTrigger.gtDigis_cfi import *

process = cms.Process("ROOTMAKER")

process.load('PhysicsTools.PatUtils.patPFMETCorrections_cff')
from PhysicsTools.PatUtils.patPFMETCorrections_cff import *
process.load('FastSimulation.HighLevelTrigger.HLTFastReco_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')

process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.categories.append('PATSummaryTables')
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))
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
        'root://cms-xrd-global.cern.ch//store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/252/00000/9ADEE140-9C27-E511-919A-02163E011D23.root'
    )
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1) 
)

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions
#process.GlobalTag.globaltag = cms.string('PHYS14_25_V1::All')
process.GlobalTag.globaltag = cms.string('74X_dataRun2_Prompt_v0')

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('AC1B_DA_mini.root')
)

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be DataFormat.AOD or DataFormat.MiniAOD, as appropriate
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
# define which IDs we want to produce
my_id_modules = [
    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_PHYS14_PU20bx25_nonTrig_V1_cff',
    'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff',
    'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff'
]
#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


# ROOTMAKER #########################################################################################
process.makeroottree = cms.EDAnalyzer("RootMaker",
    
    isMiniAOD = cms.untracked.bool(True),
    debug = cms.untracked.bool(False),

    # TRIGGER #####################################################
    HLTriggerSelection = cms.untracked.vstring(),
    TriggerProcess = cms.untracked.string('HLT'), #REDIGI311X'),
    Trigger = cms.untracked.bool(True),

    # MUONS #######################################################
    RecMuonHLTriggerMatching = cms.untracked.vstring(
         'HLT_IsoMu20_v.*:FilterTrue', 
         'HLT_IsoTkMu20_v.*:FilterTrue'
#        'HLT_Mu17_Mu8_v.*:FilterTrue',
#        'HLT_Mu17_TkMu8_v.*:FilterTrue',
#        'HLT_Mu22_TkMu22_v.*:FilterTrue',
#        'HLT_Mu(8|17)_v.*',
#        'HLT_IsoMu(24|30)_v.*',
#        'HLT_IsoMu(24|30)_eta2p1_v.*',
#        'HLT_Mu40_v.*',
#        'HLT_Mu50_eta2p1_v.*',
#        'HLT_Mu17_Mu8_v.*',
#        'HLT_Mu17_TkMu8_v.*',
#        'HLT_Mu30_TkMu11_v.*',
#        'HLT_IsoTkMu24_IterTrk02_v.*',
#        'HLT_IsoTkMu24_eta2p1_IterTrk02_v.*',
#        'HLT_DoubleMu33NoFiltersNoVtx_v.*',
#        'HLT_Mu40_eta2p1_v.*'
    ),
    RecMuon = cms.untracked.bool(True),
    RecMuonPtMin = cms.untracked.double(20),
    RecMuonTrackIso = cms.untracked.double(1000000),
    RecMuonEtaMax = cms.untracked.double(2.5),
    RecMuonNum = cms.untracked.int32(1000),
    
    # ELECTRONS and PHOTONS #######################################
    RecElectronHLTriggerMatching = cms.untracked.vstring(
        'HLT_Ele27_eta2p1_WPLoose_Gsf_v.*:FilterTrue'
#        'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v.*:FilterTrue',
#        'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v.*:FilterTrue',
#        'HLT_Ele27_WP80_v.*',
#        'HLT_Ele80_CaloIdVT_TrkIdT_v.*',
#        'HLT_Ele80_CaloIdVT_GsfTrkIdT_v.*'
    ),

    RecElectron = cms.untracked.bool(True),
    RecElectronPtMin = cms.untracked.double(20.),
    RecElectronTrackIso = cms.untracked.double(1000000.),
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
    eleVetoIdMap    = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto"),
    eleLooseIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose"),
    eleMediumIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium"),
    eleTightIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight"),
    eleHeepV60IdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
    eleMVAIdMap_wp80= cms.InputTag("egmGsfElectronIDs:mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp80"),
    eleMVAIdMap_wp90= cms.InputTag("egmGsfElectronIDs:mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp90"),

    # TAUS ########################################################
    RecTauHLTriggerMatching = cms.untracked.vstring(
    ),

    RecTau = cms.untracked.bool(True),
    RecTauDiscriminators = cms.untracked.vstring(
        #'decayModeFindingOldDMs',
        'byLooseCombinedIsolationDeltaBetaCorr3Hits',
        'byMediumCombinedIsolationDeltaBetaCorr3Hits',
        'byTightCombinedIsolationDeltaBetaCorr3Hits',
        'byCombinedIsolationDeltaBetaCorrRaw3Hits',
        'againstMuonLoose3',
        'againstMuonTight3',
        'againstMuonLooseMVA',
        'againstMuonMediumMVA',
        'againstMuonTightMVA',
        'againstElectronVLooseMVA5',
        'againstElectronLooseMVA5',
        'againstElectronMediumMVA5',
        'againstElectronTightMVA5',
        'againstElectronVTightMVA5'
        #'hpsPFTauDiscriminationByDecayModeFinding',
        #'hpsPFTauDiscriminationByLooseElectronRejection',
        #'hpsPFTauDiscriminationByLooseIsolation',
        #'hpsPFTauDiscriminationByLooseMuonRejection',
        #'hpsPFTauDiscriminationByMediumElectronRejection',
        #'hpsPFTauDiscriminationByMediumIsolation',
        #'hpsPFTauDiscriminationByTightElectronRejection',
        #'hpsPFTauDiscriminationByTightIsolation',
        #'hpsPFTauDiscriminationByTightMuonRejection',
        #'hpsPFTauDiscriminationByVLooseIsolation'
    ),
    RecTauPtMin = cms.untracked.double(0.),
    RecTauEta = cms.untracked.double(0.),
    RecTauNum = cms.untracked.int32(100000),

    # JETS ########################################################
    rhoAll = cms.InputTag("fixedGridRhoAll", "", "RECO"),

    RecJetHLTriggerMatching = cms.untracked.vstring(
    ),

    JetCorrection = cms.untracked.string('L1FastL2L3'),#MC

    RecAK4PFCHSJet = cms.untracked.bool(True),
    RecAK4PFCHSPtMin = cms.untracked.double(20.),
    RecAK4PFCHSEtaMax = cms.untracked.double(4.7),
    RecAK4PFCHSNum = cms.untracked.int32(100000),
    RecAK4PFCHSFilterPtMin = cms.untracked.double(20.),

    RecAK4PFCHSPuppiJet = cms.untracked.bool(True),
    RecAK4PFCHSPuppiPtMin = cms.untracked.double(20.),
    RecAK4PFCHSPuppiEtaMax = cms.untracked.double(4.7),
    RecAK4PFCHSPuppiNum = cms.untracked.int32(100000),
    RecAK4PFCHSPuppiFilterPtMin = cms.untracked.double(20.),
    
    # GEN PARTICLES ###############################################
    GenAllParticles = cms.untracked.bool(False),
    GenSomeParticles = cms.untracked.bool(True),
    GenAK4Jets = cms.untracked.bool(True),

    # MET #########################################################
    RecPFMet = cms.untracked.bool(True),
    RecPFMetPuppi = cms.untracked.bool(True),
    
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
    RecSecVertices = cms.untracked.bool(True),
    RecVertexTRKChi2 = cms.untracked.double(5),
    RecVertexTRKHitsMin = cms.untracked.int32(6),
    RecVertexChi2 = cms.untracked.double(3),
    RecVertexSig2D = cms.untracked.double(15),
    RecKaonMasswin = cms.untracked.double(0.05),
    RecLambdaMasswin = cms.untracked.double(0.02),

    # INPUT TAGS ##################################################
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("selectedPatTrigger"),

    dEdxharmonic2 = cms.InputTag("dEdxharmonic2"),
    l1trigger = cms.InputTag("gtDigis"),
#    metFilterBits = cms.InputTag("TriggerResults", "", "PAT"),
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
    taujets = cms.InputTag("slimmedJets"),
    ak4calojets = cms.InputTag("patJetsAK4Calo"), # AOD only
    ak4jptjets = cms.InputTag("patJetsAK4JPT"), # AOD only
    ak4pfjets = cms.InputTag("ak4PFJets"), # AOD only
    ak4pfchsjets = cms.InputTag("slimmedJets"),# mini
    ak4pfchsjetspuppi = cms.InputTag("slimmedJetsPuppi"),# mini, jets corrected for the pileup response
    jetsAK8 = cms.InputTag("slimmedJetsAK8"), # mini only
    pfMetMVAEMT = cms.InputTag("pfMetMVAEMT"), # mini only
    pfMetMVAEM  = cms.InputTag("pfMetMVAEM"), # mini only
    pfMetMVAET  = cms.InputTag("pfMetMVAET"), # mini only
    pfMetMVAMT  = cms.InputTag("pfMetMVAMT"), # mini only
    pfMetMVATT  = cms.InputTag("pfMetMVATT"), # mini only
    #newCorrectedSlimmedMetLabel = cms.InputTag("slimmedMETs","","RERUN"),
    patMETs = cms.InputTag("slimmedMETs"),
    #cms.InputTag("patMETs"),
    mets = cms.InputTag("slimmedMETs"), # mini only
    metspuppi = cms.InputTag("slimmedMETsPuppi"), # mini only
    pfmet = cms.InputTag("pfMet"), # AOD only
    pfmett1 = cms.InputTag("pfMetT1"), # AOD only
    pfmett1t0 = cms.InputTag("pfType0Type1CorrectedMet"), # AOD only
    packedGenParticles = cms.InputTag("packedGenParticles"), # mini only
    genSimParticles = cms.InputTag("prunedGenParticles"), 
    genParticles = cms.InputTag("prunedGenParticles"), 
    genJets = cms.InputTag("slimmedGenJets", "", "PAT"),
    genInfo = cms.InputTag("generator", "", "SIM"),
    puInfo = cms.InputTag("addPileupInfo", "", "HLT"),
    PfCands = cms.InputTag("PFCandidates"), # AOD only
    packedPfCands = cms.InputTag("packedPFCandidates"), # mini only
    #barrelHits = cms.InputTag("reducedEgamma", "reducedEBRecHits", "PAT"),
    #endcapHits = cms.InputTag("reducedEgamma", "reducedEERecHits", "PAT"),
    #esHits = cms.InputTag("reducedEgamma", "reducedESRecHits", "PAT"),
    barrelHits = cms.InputTag("reducedEgamma", "reducedEBRecHits", "RECO"),
    endcapHits = cms.InputTag("reducedEgamma", "reducedEERecHits", "RECO"),
    esHits = cms.InputTag("reducedEgamma", "reducedESRecHits", "RECO"),
    superClusters = cms.InputTag("reducedEgamma", "reducedSuperClusters", "PAT"),
    ebeeClusters = cms.InputTag("reducedEgamma", "reducedEBEEClusters", "PAT"),
    esClusters = cms.InputTag("reducedEgamma", "reducedESClusters", "PAT"),
    #conversions = cms.InputTag("reducedEgamma", "reducedConversions", "PAT"),
    conversions = cms.InputTag("reducedEgamma", "reducedConversions", "RECO"),
    gedGsfElectronCores = cms.InputTag("reducedEgamma", "reducedGedGsfElectronCores", "PAT"), # mini only
    gedPhotonCores = cms.InputTag("reducedEgamma", "reducedGedPhotonCores", "PAT"), # mini only

)

process.p = cms.Path(
    process.egmGsfElectronIDSequence *
    process.makeroottree
)
