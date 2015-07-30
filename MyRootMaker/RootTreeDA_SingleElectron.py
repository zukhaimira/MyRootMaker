from MyRootMaker.MyRootMaker.RootMakerTemplateDA_mini_cfg import *

process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load('RecoJets.Configuration.RecoJetAssociations_cff')
process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')
process.load('RecoBTag.Configuration.RecoBTag_cff')
process.load('RecoJets.Configuration.RecoJetAssociations_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')

process.ak4JetTracksAssociatorAtVertexPF.jets = cms.InputTag("ak4PFJetsCHS")
process.ak4JetTracksAssociatorAtVertexPF.tracks = cms.InputTag("unpackedTracksAndVertices")
process.impactParameterTagInfos.primaryVertex = cms.InputTag("unpackedTracksAndVertices")
process.inclusiveSecondaryVertexFinderTagInfos.extSVCollection = cms.InputTag("unpackedTracksAndVertices","secondary","")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch//store/data/Run2015B/SingleElectron/MINIAOD/PromptReco-v1/000/251/096/00000/22D22D7F-5626-E511-BDE3-02163E011FAB.root'
    )
)
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions
#process.GlobalTag.globaltag = cms.string('MCRUN2_74_V9')
process.GlobalTag.globaltag = cms.string('74X_dataRun2_Prompt_v0')
process.makeroottree.isMiniAOD = cms.untracked.bool(True)

process.maxEvents = cms.untracked.PSet(
    #input = cms.untracked.int32(1000) 
    input = cms.untracked.int32(-1) 
)
#process.makeroottree.debug = cms.untracked.bool(True)

process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.TFileService = cms.Service("TFileService",
	fileName = cms.string('AC1B_DA_SingleElectron_mini_test.root')
)
#process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))
####
#process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
#process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")
#process.load("PhysicsTools.PatUtils.patPFMETCorrections_cff")
#from PhysicsTools.PatAlgos.tools.jetTools import switchJetCollection
#switchJetCollection(process,
#                    jetSource = cms.InputTag('slimmedJets'),
#                    jetCorrections = ('AK4PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], '')
#                    )

# apply type I/type I + II PFMEt corrections to pat::MET object
# and estimate systematic uncertainties on MET
#from PhysicsTools.PatUtils.tools.runType1PFMEtUncertainties import runType1PFMEtUncertainties
#runType1PFMEtUncertainties(process,addToPatDefaultSequence=False,
#                           photonCollection="slimmedPhotons",
#                           jetCollection="slimmedJets",
#                           electronCollection="slimmedElectrons",
#                           muonCollection="slimmedMuons",
#                           tauCollection="slimmedTaus")
####

from RecoMET.METPUSubtraction.objectSelection_cff import *
from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *
from JetMETCorrections.Configuration.DefaultJEC_cff import *
from RecoJets.JetProducers.PileupJetIDParams_cfi import JetIdParams

##initial miniaod ak4pfjets have to be reprocessed
process.load("RecoJets.JetProducers.ak4PFJets_cfi")
process.calibratedAK4PFJetsForPFMVAMEt = cms.EDProducer('PFJetCorrectionProducer',
    src = cms.InputTag('packedPFCandidates'),#ak4PFJets 
    #correctors = cms.vstring("ak4PFL1FastL2L3") # NOTE: use "ak5PFL1FastL2L3" for MC / "ak5PFL1FastL2L3Residual" for Data
    correctors = cms.vstring("ak4PFL1FastL2L3Residual") # NOTE: use "ak5PFL1FastL2L3Residual" for Data
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

###mvamet module
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
###electron selection should be specified
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

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
#default configuration for miniAOD reprocessing, change the isData flag to run on data
# runMetCorAndUncFromMiniAOD(process, isData=False)#MC
runMetCorAndUncFromMiniAOD(process, isData=True)#Data

### -------------------------------------------------------------------
### the lines below remove the L2L3 residual uncertainties when processing data
### -------------------------------------------------------------------
process.patPFMetT1T2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
process.patPFMetT1T2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
process.patPFMetT2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
process.patPFMetT2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
process.shiftedPatJetEnDown.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
process.shiftedPatJetEnUp.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
### ------------------------------------------------------------------
from PhysicsTools.PatUtils.tools.runType1PFMEtUncertainties import runType1PFMEtUncertainties
runType1PFMEtUncertainties(process,addToPatDefaultSequence=False,
                           photonCollection="slimmedPhotons",
                           jetCollection="slimmedJets",
                           electronCollection="slimmedElectrons",
                           muonCollection="slimmedMuons",
                           tauCollection="slimmedTaus")

process.makeroottree.RecMuonNum = cms.untracked.int32(0)
process.makeroottree.RecElectronNum = cms.untracked.int32(0)
process.makeroottree.RecAK4PFCHSNum = cms.untracked.int32(0)
process.makeroottree.RecAK4PFCHSPuppiNum = cms.untracked.int32(0)
process.makeroottree.HLTriggerSelection = cms.untracked.vstring()
process.makeroottree.GenAllParticles = cms.untracked.bool(True)
process.makeroottree.GenSomeParticles = cms.untracked.bool(True)
process.makeroottree.GenAK4Jets = cms.untracked.bool(True)

process.p = cms.Path(
    process.mvaMET *
    process.egmGsfElectronIDSequence *
    process.makeroottree
)

