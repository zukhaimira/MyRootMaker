import FWCore.ParameterSet.Config as cms

process = cms.Process("ROOTMAKER")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")

	
process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('PATSummaryTables')
process.GlobalTag.globaltag = cms.string('START70_V7::All')

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

process.source = cms.Source("PoolSource", 
     fileNames = cms.untracked.vstring(
#'/store/mc/Summer12_DR53X/T_t-channel_TuneZ2star_8TeV-powheg-tauola/AODSIM/PU_S10_START53_V7A-v3/00000/90B3B34B-B00F-E211-B13A-003048678E24.root'
'/store/mc/Spring14dr/GluGluToHToMuMu_M-125_13TeV-powheg-pythia6/AODSIM/PU_S14_POSTLS170_V6-v1/00000/1E93B8DB-7CFD-E311-BFC7-7845C4FC346A.root'
)
#     noEventSort = cms.untracked.bool(True),
#     duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1) 
)
# The Good vertices collection _____________________________________________||
process.goodVertices = cms.EDFilter(
                "VertexSelector",
                filter = cms.bool(False),
                src = cms.InputTag("offlinePrimaryVertices"),
                cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
                )
process.vertex_step=cms.Path(process.goodVertices)
#####Filter
## The good primary vertex filter ____________________________________________||
process.primaryVertexFilter = cms.EDFilter(
                "VertexSelector",
                src = cms.InputTag("offlinePrimaryVertices"),
                cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
                filter = cms.bool(True)
                )
## The beam scraping filter __________________________________________________||
process.noscraping = cms.EDFilter(
                "FilterOutScraping",
                applyfilter = cms.untracked.bool(True),
                debugOn = cms.untracked.bool(False),
                numtrack = cms.untracked.uint32(10),
                thresh = cms.untracked.double(0.25)
                )
## The iso-based HBHE noise filter ___________________________________________||
process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')
## The CSC beam halo tight filter ____________________________________________||
process.load('RecoMET.METFilters.CSCTightHaloFilter_cfi')
## The HCAL laser filter _____________________________________________________||
process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")
## The ECAL dead cell trigger primitive filter _______________________________||
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
## The EE bad SuperCrystal filter ____________________________________________||
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
## The ECAL laser correction filter
process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')
## The tracking failure filter _______________________________________________||
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')

process.filters_step = cms.Path(
                process.primaryVertexFilter *
                process.noscraping *
                process.HBHENoiseFilter *
                process.CSCTightHaloFilter *
                process.hcalLaserEventFilter *
                process.EcalDeadCellTriggerPrimitiveFilter *
                process.trackingFailureFilter *
                process.eeBadScFilter *
                process.ecalLaserCorrFilter
                )

####### Jet MET corrections
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.kt6PFJets.doRhoFastjet = True
process.kt6PFJets.doAreaFastjet = True
process.kt6PFJets.voronoiRfact = 0.9

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('JetMETCorrections.Configuration.JetCorrectionServices_cff')
process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")

#process.pfJetMETcorr.jetCorrLabel = cms.string("ak4PFL1FastL2L3") #MC
#from JetMETCorrections.Type1MET.pfMETCorrections_cff import pfType1CorrectedMet
process.load("JetMETCorrections.Type1MET.correctedMet_cff")
#process.pfType0Type1CorrectedMet = pfType1CorrectedMet.clone(
#    applyType0Corrections = cms.bool(False),
#    srcType1Corrections = cms.VInputTag(
#        cms.InputTag('pfMETcorrType0'),
#        cms.InputTag('pfJetMETcorr', 'type1')
#    )
#)
#process.metAnalysisSequence=cms.Sequence(process.type0PFMEtCorrection*process.producePFMETCorrections*process.pfType0Type1CorrectedMet)
from JetMETCorrections.Type1MET.correctedMet_cff import pfMetT1
pfMetT1seq = pfMetT1.clone (
    src = cms.InputTag('ak4PFJets'),
    srcType1Corrections = cms.VInputTag(
        cms.InputTag('pfMETcorrType0'),
        cms.InputTag('pfJetMETcorr', 'type1')
    ),
    offsetCorrLabel = cms.string("ak4PFL1Fastjet"),
    jetCorrLabel = cms.string("ak4PFL1FastL2L3"), # NOTE: use "ak4PFL1FastL2L3" for MC / "ak4PFL1FastL2L3Residual" for Data
    jetCorrEtaMax = cms.double(9.9),
    type1JetPtThreshold = cms.double(10.0),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.90),
    skipMuons = cms.bool(True),
)
process.metAnalysisSequence=cms.Sequence(process.pfMetT1)
process.jet_step = cms.Path(process.kt6PFJets*process.metAnalysisSequence)

######PF ISO calculation for Electrons
process.load('CommonTools.ParticleFlow.Isolation.pfElectronIsolation_cff')
from CommonTools.ParticleFlow.Isolation.pfElectronIsolation_cff import *
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFPhotonIso

process.stdElectronSequencePFIso = setupPFElectronIso(process, 'gedGsfElectrons')
process.stdPhotonSequencePFIso = setupPFPhotonIso(process, 'gedPhotons')
process.pfiso_step = cms.Path( process.pfParticleSelectionSequence + process.stdElectronSequencePFIso + process.stdPhotonSequencePFIso)


######Electron ID
process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_cfi")
process.electron_step = cms.Path(process.eidHyperTight1MC*process.eidLooseMC)

######Matching Partons to Jets
process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")
process.AK4byRef.jets = cms.InputTag("ak4PFJets")
process.jetflavour_step = cms.Path(process.myPartons * process.AK4Flavour)

######PAT
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
		fileName = cms.untracked.string('patTuple.root'),
# save only events passing the full path
		SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
# save PAT Layer 1 output; you need a '*' to
# unpack the list of commands 'patEventContent'
		outputCommands = cms.untracked.vstring('drop *', *patEventContent )
		)
#from PhysicsTools.PatAlgos.tools.coreTools import *
#removeAllPATObjectsBut(process, ['Jets'])
from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT
postfix = "PFlow"
usePF2PAT(process,runPF2PAT=True,
		jetAlgo='AK4', runOnMC=False, postfix=postfix,
		#jetCorrections=('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute']),
                jetCorrections=('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'],'Type-1'),

		pvCollection=cms.InputTag('offlinePrimaryVertices')
	 )
process.pfPileUpPFlow.checkClosestZVertex = False
#getattr(process, "patPF2PATSequence"+postfix).remove(getattr(process, "cleanPatCandidates"+postfix))
#getattr(process, "patPF2PATSequence"+postfix).remove(getattr(process, "countPatCandidates"+postfix))
#process.patseq = cms.Sequence(getattr(process,"patPF2PATSequence"+postfix))
#process.pat_step = cms.Path(process.patseq)

######Pileup Jet ID (Hendrik)
#https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetID#Running_the_algorithm_on_reco_je
#from CMGTools.External.pujetidsequence_cff import puJetId, puJetMva
#
#from RecoJets.JetProducers.puJetIDAlgo_cff import * 
#from RecoJets.JetProducers.puJetIDParams_cfi import * 
from RecoJets.JetProducers.PileupJetID_cfi import pileupJetIdProducer, _stdalgos_4x, _stdalgos_5x, _stdalgos, cutbased, _chsalgos_4x, _chsalgos_5x, _chsalgos 
from RecoJets.JetProducers.PileupJetID_cfi import * 
#process.load("RecoJets.JetProducers.PileupJetID_cfi")
#process.load("RecoJets.JetProducers.puJetIDAlgo_cff")

process.recoPuJetId = pileupJetIdProducer.clone(
    jets = cms.InputTag("ak4PFJets"),
    applyJec = cms.bool(True),
    inputIsCorrected = cms.bool(False), 
    produceJetIds = cms.bool(True),
    jetids = cms.InputTag(""),
    runMvas = cms.bool(False),
    vertexes = cms.InputTag("offlinePrimaryVertices"),
    algos = cms.VPSet(cutbased)
)

process.recoPuJetMva = pileupJetIdProducer.clone(
    jets = cms.InputTag("ak4PFJets"),
    jetids = cms.InputTag("recoPuJetId"),
    applyJec = cms.bool(True),
    produceJetIds = cms.bool(False),
    runMvas = cms.bool(True),
    vertexes = cms.InputTag("offlinePrimaryVertices"),
    algos = cms.VPSet(_stdalgos),
    inputIsCorrected = cms.bool(True),                                     
)


#process.recoPuJetIdSequence = cms.Sequence(process.recoPuJetId * process.recoPuJetMva)
process.recoPuJetIdSequence = cms.Sequence(process.recoPuJetId)
process.jetpuid_step = cms.Path(process.recoPuJetIdSequence)

#from CMGTools.External.pujetidproducer_cfi import stdalgos_4x, stdalgos_5x, stdalgos, cutbased, chsalgos_4x, chsalgos_5x, chsalgos
#pileupJetIdProducer = cms.EDProducer('PileupJetIdProducer',
#                         produceJetIds = cms.bool(True),
#                         jetids = cms.InputTag(""),
#                         runMvas = cms.bool(True),
#                         jets = cms.InputTag("selectedPatJetsPFlow"),
#                         vertexes = cms.InputTag("offlinePrimaryVertices"),
#                         algos = cms.VPSet(stdalgos),
#                                     
#                         rho     = cms.InputTag("fixedGridRhoFastjetAll"),
#                         jec     = cms.string("AK4PF"),
#                         applyJec = cms.bool(False),
#                         inputIsCorrected = cms.bool(True),                                     
#                         residualsFromTxt = cms.bool(False),
#                       #  residualsTxt     = cms.FileInPath("RecoJets/JetProducers/data/dummy.txt"),
#)
#
#
#
#puJetId = pileupJetIdProducer.clone(
#    produceJetIds = cms.bool(True),
#    jetids = cms.InputTag(""),
#    runMvas = cms.bool(False),
#    jets = cms.InputTag("selectedPatJets"),
#    vertexes = cms.InputTag("offlinePrimaryVertices"),
#    algos = cms.VPSet(cutbased)
#)
#puJetMva = pileupJetIdProducer.clone(
#    produceJetIds = cms.bool(False),
#    jetids = cms.InputTag("puJetId"),
#    runMvas = cms.bool(True),
#    jets = cms.InputTag("selectedPatJets"),
#    vertexes = cms.InputTag("offlinePrimaryVertices"),
#    algos = stdalgos
#)
#
#process.recoPuJetId = puJetId.clone(
#    jets = cms.InputTag("ak4PFJets"),
#    applyJec = cms.bool(True),
#    inputIsCorrected = cms.bool(False),
#)
#process.recoPuJetMva = puJetMva.clone(
#    jets = cms.InputTag("ak4PFJets"),
#    jetids = cms.InputTag("recoPuJetId"),
#    applyJec = cms.bool(True),
#    inputIsCorrected = cms.bool(False),
#)
#



process.recoPuJetIdSequence = cms.Sequence(process.recoPuJetId)
process.jetpuid_step = cms.Path(process.recoPuJetIdSequence)



######ROOTMAKER 
process.makeroottree = cms.EDAnalyzer("RootMaker",

    ebRecHits = cms.InputTag("reducedEcalRecHitsEB"),
    eeRecHits = cms.InputTag("reducedEcalRecHitsEB"),
    esRecHits = cms.InputTag("reducedEcalRecHitsES"),

    GenAllParticles = cms.untracked.bool(False),
    GenSomeParticles = cms.untracked.bool(True),
    GenAK4Jets = cms.untracked.bool(True),
    
    Trigger = cms.untracked.bool(True),
    HLTriggerSelection = cms.untracked.vstring(),
    #TriggerProcess = cms.untracked.string('REDIGI311X'),
    RecPrimVertex = cms.untracked.bool(True),
    RecBeamSpot = cms.untracked.bool(True),
    
    RecPFMet = cms.untracked.bool(True),
    
    RecTrack = cms.untracked.bool(True),
    RecTrackPtMin = cms.untracked.double(10.),
    RecTrackEtaMax = cms.untracked.double(2.5),
    RecTrackNum = cms.untracked.int32(100000),
    RecTrackFilterPtMin = cms.untracked.double(18.),
    
    RecSuperCluster = cms.untracked.bool(True),
    RecSuperClusterFilterPtMin = cms.untracked.double(8.),
    RecSuperClusterBasicCluster = cms.untracked.bool(False),
    RecSuperClusterHit = cms.untracked.bool(False),
    
    RecMuon = cms.untracked.bool(True),
    RecMuonPtMin = cms.untracked.double(20),
    RecMuonTrackIso = cms.untracked.double(1000000),
    RecMuonEtaMax = cms.untracked.double(2.5),
    RecMuonNum = cms.untracked.int32(1000),
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
    #RecMassMuMuMin = cms.untracked.double(2.6),
    #RecMassMuMuMax = cms.untracked.double(3.5),
    
    RecPhoton = cms.untracked.bool(True),
    RecPhotonHLTriggerMatching = cms.untracked.vstring(),
    RecPhotonPtMin = cms.untracked.double(10.),
    RecPhotonEtaMax = cms.untracked.double(2.5),
    RecPhotonNum = cms.untracked.int32(100000),
    RecPhotonFilterPtMin = cms.untracked.double(10),
    
    RecElectron = cms.untracked.bool(True),
    RecElectronHLTriggerMatching = cms.untracked.vstring(
    'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v.*:FilterTrue',
    'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v.*:FilterTrue',
    'HLT_Ele27_WP80_v.*',
    'HLT_Ele80_CaloIdVT_TrkIdT_v.*',
    'HLT_Ele80_CaloIdVT_GsfTrkIdT_v.*'
    ),
    RecElectronPtMin = cms.untracked.double(20.),
    RecElectronTrackIso = cms.untracked.double(1000000.),
    RecElectronEta = cms.untracked.double(2.5),
    RecElectronNum = cms.untracked.int32(100000),
    RecElectronFilterPtMin = cms.untracked.double(20.),
    
    RecTau = cms.untracked.bool(False),
    RecTauHLTriggerMatching = cms.untracked.vstring(),
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
    
    RecAK4CaloJet = cms.untracked.bool(False),
    RecAK4CaloPtMin = cms.untracked.double(20.),
    RecAK4CaloEtaMax = cms.untracked.double(2.4),
    RecAK4CaloNum = cms.untracked.int32(100000),
    RecAK4CaloFilterPtMin = cms.untracked.double(20.),
    
    RecAK4JPTJet = cms.untracked.bool(False),
    RecAK4JPTPtMin = cms.untracked.double(20.),
    RecAK4JPTEtaMax = cms.untracked.double(2.4),
    RecAK4JPTNum = cms.untracked.int32(100000),
    RecAK4JPTFilterPtMin = cms.untracked.double(20.),
    
    #JetCorrection = cms.untracked.string('L1FastL2L3Residual'),#Data
    JetCorrection = cms.untracked.string('L1FastL2L3'),#MC
    RecJetHLTriggerMatching = cms.untracked.vstring(),
    
    RecAK4PFJet = cms.untracked.bool(True),
    RecAK4PFPtMin = cms.untracked.double(20.),
    RecAK4PFEtaMax = cms.untracked.double(3.0),
    RecAK4PFNum = cms.untracked.int32(100000),
    RecAK4PFFilterPtMin = cms.untracked.double(20.),
    
    RecAK4PFCHSJet = cms.untracked.bool(True),
    RecAK4PFCHSPtMin = cms.untracked.double(20.),
    RecAK4PFCHSEtaMax = cms.untracked.double(3.0),
    RecAK4PFCHSNum = cms.untracked.int32(100000),
    RecAK4PFCHSFilterPtMin = cms.untracked.double(20.),
    
    RecSecVertices = cms.untracked.bool(False),
    RecVertexTRKChi2 = cms.untracked.double(5),
    RecVertexTRKHitsMin = cms.untracked.int32(6),
    RecVertexChi2 = cms.untracked.double(3),
    RecVertexSig2D = cms.untracked.double(15),
    RecKaonMasswin = cms.untracked.double(0.05),
    RecLambdaMasswin = cms.untracked.double(0.02)
)

#process.genPlusSimParticles = cms.EDProducer("GenPlusSimParticleProducer",
#src = cms.InputTag("g4SimHits"),
#setStatus = cms.int32(8),
#filter = cms.vstring("pt > 0.0"),
#genParticles = cms.InputTag("genParticles")
#)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('AC1B_test.root')
)

process.roottree_step = cms.EndPath(process.makeroottree)

