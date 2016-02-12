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
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch//store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/252/00000/9ADEE140-9C27-E511-919A-02163E011D23.root'
    )
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1) 
)

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions
process.GlobalTag.globaltag = cms.string('74X_dataRun2_Prompt_v0')

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('AC1B_76DATA.root')
)


# ELECTRON ID #######################################################################################
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


# MVA MET ###########################################################################################
from RecoMET.METPUSubtraction.objectSelection_cff import *
from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *
from JetMETCorrections.Configuration.DefaultJEC_cff import *
from RecoJets.JetProducers.PileupJetIDParams_cfi import JetIdParams

##initial miniaod ak4pfjets have to be reprocessed
process.load("RecoJets.JetProducers.ak4PFJets_cfi")
process.calibratedAK4PFJetsForPFMVAMEt = cms.EDProducer('PFJetCorrectionProducer',
    src = cms.InputTag('packedPFCandidates'),
    correctors = cms.vstring("ak4PFL1FastL2L3Residual") # Data
)

#from JetMETCorrections.Configuration.DefaultJEC_cff import ak4PFJetsL1FastL2L3 # MC
from JetMETCorrections.Configuration.DefaultJEC_cff import ak4PFL1FastL2L3Residual # Data
from RecoJets.JetProducers.pileupjetidproducer_cfi import pileupJetIdEvaluator
from RecoJets.JetProducers.PileupJetIDParams_cfi import JetIdParams
process.puJetIdForPFMVAMEt = pileupJetIdEvaluator.clone(
    algos = cms.VPSet(
        cms.PSet(
        tmvaVariables = cms.vstring(
            "nvtx",
            "jetPt",
            "jetEta",
            "jetPhi",
            "dZ",
            "beta",
            "betaStar",
            "nCharged",
            "nNeutrals",
            "dR2Mean",
            "ptD",
            "frac01",
            "frac02",
            "frac03",
            "frac04",
            "frac05"
            ),

        tmvaWeights = cms.string("RecoJets/JetProducers/data/TMVAClassificationCategory_JetID_MET_53X_Dec2012.weights.xml.gz"),
        tmvaMethod = cms.string("JetID"),
        tmvaSpectators = cms.vstring(),
        JetIdParams = JetIdParams,
        impactParTkThreshold = cms.double(0.),
        version = cms.int32(-1),
        cutBased = cms.bool(False),
        label = cms.string("full")
        )
        ),
    produceJetIds = cms.bool(True),
    runMvas = cms.bool(True),
    rho     = cms.InputTag("fixedGridRhoFastjetAll"),
    vertexes = cms.InputTag("offlinePrimaryVertices"),#offlineSlimmedPrimaryVertices
    jets = cms.InputTag("calibratedAK4PFJetsForPFMVAMEt"),
    applyJec = cms.bool(True),
    inputIsCorrected = cms.bool(True),
    jec     = cms.string("AK4PF"),
)

## mvamet module
process.load("RecoMET.METPUSubtraction.mvaPFMET_cff")
process.pfMVAMEt.srcPFCandidates = cms.InputTag("packedPFCandidates")
process.pfMVAMEt.srcVertices = cms.InputTag("offlineSlimmedPrimaryVertices")
process.puJetIdForPFMVAMEt.vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")

## miniAOD settings
process.ak4PFJets.src = cms.InputTag("packedPFCandidates")
process.pfMVAMEt.srcPFCandidates = cms.InputTag('packedPFCandidates')
process.pfMVAMEt.srcVertices = cms.InputTag('offlineSlimmedPrimaryVertices')

## muons as input for mvaMET producer
process.mvaMETMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag('slimmedMuons'),#patMuons
    cut = cms.string(
       "abs(eta)<2.6 & pt>9.5"                                      +
        ## muon ID
        "& isTrackerMuon"                                           +
        "& isPFMuon"                                                +
        "& globalTrack.isNonnull"                                   +
        "& innerTrack.hitPattern.numberOfValidPixelHits    >  0"    +
        "& innerTrack.normalizedChi2                       < 10"    +
        "& numberOfMatches                                 >  0"    +
        "& innerTrack.hitPattern.numberOfValidTrackerHits  >  5"    +
        "& globalTrack.hitPattern.numberOfValidHits        >  0"    +
        "& abs(innerTrack().dxy)                           <2.0"    +
        ## muon isolation (w/o deltaBeta, therefore weaker selection criteria)
        "& (pfIsolationR03.sumChargedHadronPt+pfIsolationR03.sumNeutralHadronEt+pfIsolationR03.sumPhotonEt)/pt < 0.3"
        ),
    filter = cms.bool(False)
)

## electrons as input for mvaMET producer
process.mvaMETElectrons = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag('slimmedElectrons'),#patElectrons
    cut = cms.string(
        "abs(eta) < 2.6 && pt > 9.5"                                 +
        #"&& gsfTrack.trackerExpectedHitsInner.numberOfHits == 0"    + #to be adjusted for 72X releses
        ## electron ID for barrel electrons
        "&& ((abs(eta) < 1.4442  "                                  +
        "&& abs(deltaEtaSuperClusterTrackAtVtx)            < 0.007" +
        "&& abs(deltaPhiSuperClusterTrackAtVtx)            < 0.8"   +
        "&& sigmaIetaIeta                                  < 0.01"  +
        "&& hcalOverEcal                                   < 0.15"  +
        "&& abs(1./superCluster.energy - 1./p)             < 0.05)" +
        ## electron ID for endcap electrons
        "|| (abs(eta)  > 1.566 "                                    +
        "&& abs(deltaEtaSuperClusterTrackAtVtx)            < 0.009 "+
        "&& abs(deltaPhiSuperClusterTrackAtVtx)            < 0.10"  +
        "&& sigmaIetaIeta                                  < 0.03"  +
        "&& hcalOverEcal                                   < 0.10"  +
        "&& abs(1./superCluster.energy - 1./p)             < 0.05))"#+
        ## electron isolation (w/o deltaBeta, therefore weaker selection criteria)
        #"&& (pfIsolationVariables.chargedHadronIso+pfIsolationVariables.neutralHadronIso)/et < 0.3"
        ),
   filter = cms.bool(False)
)

## taus as input for mvaMET producer
##  - NOTE that the selection for taus depends on the final state!!! to be checked!
process.mvaMETTausET = cms.EDFilter("PATTauRefSelector",
    src = cms.InputTag('slimmedTaus'),# patTaus
        cut = cms.string(
        'abs(eta) < 2.6 && pt > 20'                                       +
        ' & tauID("decayModeFinding") > 0.5'                              +
        ' & tauID("againstMuonLoose3") > 0.5'                             +
        ' & tauID("againstElectronMediumMVA5") > 0.5'                     +
        ' & tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5'
        ),
    filter = cms.bool(False)
)

process.mvaMETTausMT = cms.EDFilter("PATTauRefSelector",
    src = cms.InputTag('slimmedTaus'),#patTaus
        cut = cms.string(
        'abs(eta) < 2.6 && pt > 20'                                       +
        ' & tauID("decayModeFinding") > 0.5'                              +
        ' & tauID("againstMuonLoose3") > 0.5'                             +
        ' & tauID("againstElectronMediumMVA5") > 0.5'                     +
        ' & tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5'
        ),
    filter = cms.bool(False)
)

process.mvaMETTausTT = cms.EDFilter("PATTauRefSelector",
    src = cms.InputTag('slimmedTaus'),#patTaus
        cut = cms.string(
        'abs(eta) < 2.6 && pt > 20'                                       +
        ' & tauID("decayModeFinding") > 0.5'                              +
        ' & tauID("againstMuonLoose3") > 0.5'                             +
        ' & tauID("againstElectronMediumMVA5") > 0.5'                     +
        ' & tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5'
        ),
    filter = cms.bool(False)
)

## specify the leptons similar to those used in the analysis (channel specific)
process.pfMetMVAEMT = process.pfMVAMEt.clone(srcLeptons = cms.VInputTag("mvaMETElectrons", "mvaMETMuons", "mvaMETTausTT" ))
## electron selection should be specified
process.pfMetMVAEM = process.pfMVAMEt.clone(srcLeptons = cms.VInputTag("mvaMETElectrons", "mvaMETMuons" ))
process.pfMetMVAET = process.pfMVAMEt.clone(srcLeptons = cms.VInputTag("mvaMETElectrons", "mvaMETTausET"))
process.pfMetMVAMT = process.pfMVAMEt.clone(srcLeptons = cms.VInputTag("mvaMETMuons"    , "mvaMETTausMT"))
process.pfMetMVATT = process.pfMVAMEt.clone(srcLeptons = cms.VInputTag("mvaMETTausTT"))

process.pfMetMVAEMT.srcPFCandidates = cms.InputTag('packedPFCandidates')
process.pfMetMVAEMT.srcVertices = cms.InputTag('offlineSlimmedPrimaryVertices')
process.pfMetMVAEM.srcPFCandidates = cms.InputTag('packedPFCandidates')
process.pfMetMVAEM.srcVertices = cms.InputTag('offlineSlimmedPrimaryVertices')
process.pfMetMVAET.srcPFCandidates = cms.InputTag('packedPFCandidates')
process.pfMetMVAET.srcVertices = cms.InputTag('offlineSlimmedPrimaryVertices')
process.pfMetMVAMT.srcPFCandidates = cms.InputTag('packedPFCandidates')
process.pfMetMVAMT.srcVertices = cms.InputTag('offlineSlimmedPrimaryVertices')
process.pfMetMVATT.srcPFCandidates = cms.InputTag('packedPFCandidates')
process.pfMetMVATT.srcVertices = cms.InputTag('offlineSlimmedPrimaryVertices')
process.mvaMETMuons.src = cms.InputTag('slimmedMuons')
process.mvaMETElectrons.src = cms.InputTag('slimmedElectrons')
process.mvaMETTausET.src = cms.InputTag('slimmedTaus')
process.mvaMETTausMT.src = cms.InputTag('slimmedTaus')
process.mvaMETTausTT.src = cms.InputTag('slimmedTaus')

# define mva sequence
process.mvaMET = cms.Sequence(
    process.ak4PFJets*
    process.calibratedAK4PFJetsForPFMVAMEt*
    process.puJetIdForPFMVAMEt*
    process.mvaMETMuons *
    process.mvaMETElectrons *
    process.mvaMETTausET *
    process.mvaMETTausMT *
    process.mvaMETTausTT *
    process.pfMetMVAEMT *
    process.pfMetMVAEM *
    process.pfMetMVAET *
    process.pfMetMVAMT *
    process.pfMetMVATT
)

# MET CORRECTIONS/UNCERTAINTIES #####################################################################
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
#default configuration for miniAOD reprocessing, change the isData flag to run on data
# runMetCorAndUncFromMiniAOD(process, isData=False)#MC
runMetCorAndUncFromMiniAOD(process, isData=True)#Data

### the lines below remove the L2L3 residual uncertainties when processing data
process.patPFMetT1T2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
process.patPFMetT1T2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
process.patPFMetT2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
process.patPFMetT2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
process.shiftedPatJetEnDown.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
process.shiftedPatJetEnUp.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")

from PhysicsTools.PatUtils.tools.runType1PFMEtUncertainties import runType1PFMEtUncertainties
runType1PFMEtUncertainties(process,addToPatDefaultSequence=False,
                           photonCollection="slimmedPhotons",
                           jetCollection="slimmedJets",
                           electronCollection="slimmedElectrons",
                           muonCollection="slimmedMuons",
                           tauCollection="slimmedTaus")


# ROOTMAKER #########################################################################################
process.makeroottree = cms.EDAnalyzer("RootMaker",
    
    isMC = cms.untracked.bool(False),
    debug = cms.untracked.bool(False),

    # TRIGGER #####################################################
    HLTriggerSelection = cms.untracked.vstring(),
    TriggerProcess = cms.untracked.string('HLT'), #REDIGI311X'),
    Trigger = cms.untracked.bool(True),

    # MUONS #######################################################
    RecMuonHLTriggerMatching = cms.untracked.vstring(
         'HLT_IsoMu20_v.*:FilterTrue', 
         'HLT_IsoTkMu20_v.*:FilterTrue'
    ),
    RecMuon = cms.untracked.bool(True),
    RecMuonPtMin = cms.untracked.double(20),
    RecMuonTrackIso = cms.untracked.double(1000000),
    RecMuonEtaMax = cms.untracked.double(2.5),
    RecMuonNum = cms.untracked.int32(1000),
    
    # ELECTRONS and PHOTONS #######################################
    eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
    eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
    eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
    #eleHeepV60IdMap  = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
    #eleMVAIdMap_wp80 = cms.InputTag("egmGsfElectronIDs:mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp80"),
    #eleMVAIdMap_wp90 = cms.InputTag("egmGsfElectronIDs:mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp90"),

    phoLooseIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-loose"),
    phoMediumIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-medium"),
    phoTightIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-tight"),

    RecElectronHLTriggerMatching = cms.untracked.vstring(
        'HLT_Ele27_eta2p1_WPLoose_Gsf_v.*:FilterTrue'
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

    # TAUS ########################################################
    RecTauHLTriggerMatching = cms.untracked.vstring(
    ),

    RecTau = cms.untracked.bool(True),
    RecTauDiscriminators = cms.untracked.vstring(
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
    lostTracks = cms.InputTag("lostTracks", "", "PAT"),
    unpackedTracks = cms.InputTag("unpackedTracksAndVertices"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    secondaryVertices = cms.InputTag("slimmedSecondaryVertices", "", "PAT"),
    beamSpot = cms.InputTag("offlineBeamSpot", "", "RECO"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    photons = cms.InputTag("slimmedPhotons"),
    taus = cms.InputTag("slimmedTaus"),
    taujets = cms.InputTag("slimmedJets"),
    ak4pfchsjets = cms.InputTag("slimmedJets"),
    ak4pfchsjetspuppi = cms.InputTag("slimmedJetsPuppi"),# mini, jets corrected for the pileup response
    jetsAK8 = cms.InputTag("slimmedJetsAK8"),
    pfMetMVAEMT = cms.InputTag("pfMetMVAEMT"),
    pfMetMVAEM  = cms.InputTag("pfMetMVAEM"),
    pfMetMVAET  = cms.InputTag("pfMetMVAET"),
    pfMetMVAMT  = cms.InputTag("pfMetMVAMT"),
    pfMetMVATT  = cms.InputTag("pfMetMVATT"),
    #newCorrectedSlimmedMetLabel = cms.InputTag("slimmedMETs","","RERUN"),
    patMETs = cms.InputTag("slimmedMETs"), #cms.InputTag("patMETs"),
    mets = cms.InputTag("slimmedMETs"),
    metspuppi = cms.InputTag("slimmedMETsPuppi"),
    packedGenParticles = cms.InputTag("packedGenParticles"),
    genSimParticles = cms.InputTag("prunedGenParticles"), 
    genParticles = cms.InputTag("prunedGenParticles"), 
    genJets = cms.InputTag("slimmedGenJets", "", "PAT"),
    genInfo = cms.InputTag("generator", "", "SIM"),
    puInfo = cms.InputTag("addPileupInfo", "", "HLT"),
    packedPfCands = cms.InputTag("packedPFCandidates"),
    barrelHits = cms.InputTag("reducedEgamma", "reducedEBRecHits", "RECO"),
    endcapHits = cms.InputTag("reducedEgamma", "reducedEERecHits", "RECO"),
    esHits = cms.InputTag("reducedEgamma", "reducedESRecHits", "RECO"),
    superClusters = cms.InputTag("reducedEgamma", "reducedSuperClusters", "PAT"),
    ebeeClusters = cms.InputTag("reducedEgamma", "reducedEBEEClusters", "PAT"),
    esClusters = cms.InputTag("reducedEgamma", "reducedESClusters", "PAT"),
    conversions = cms.InputTag("reducedEgamma", "reducedConversions", "RECO"),
    gedGsfElectronCores = cms.InputTag("reducedEgamma", "reducedGedGsfElectronCores", "PAT"),
    gedPhotonCores = cms.InputTag("reducedEgamma", "reducedGedPhotonCores", "PAT"),
)

process.p = cms.Path(
    process.mvaMET *
    process.egmGsfElectronIDSequence *
    process.makeroottree
)
