import FWCore.ParameterSet.Config as cms
#from MyRootMaker.MyRootMaker.RootMakerTemplateMC_mini_cfg import *
#from L1Trigger.GlobalTrigger.gtDigis_cfi import *

process = cms.Process("ROOTMAKER")

#process.load('FastSimulation.HighLevelTrigger.HLTFastReco_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")


###################################################
#              _   _                 
#   ___  _ __ | |_(_) ___  _ __  ___ 
#  / _ \| '_ \| __| |/ _ \| '_ \/ __|
# | (_) | |_) | |_| | (_) | | | \__ \
#  \___/| .__/ \__|_|\___/|_| |_|___/
#       |_|                          
###################################################
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(3000) )

#process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))

#process.options.allowUnscheduled = cms.untracked.bool(True)
process.options = cms.untracked.PSet(allowUnscheduled = cms.untracked.bool(True))

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

# examples of things that go in the short version
#process.makeroottree.RecMuonNum = cms.untracked.int32(0)
#process.makeroottree.HLTriggerSelection = cms.untracked.vstring()
#process.makeroottree.GenAllParticles = cms.untracked.bool(False)
#process.makeroottree.GenSomeParticles = cms.untracked.bool(False)
#process.makeroottree.GenAK4Jets = cms.untracked.bool(False)


###################################################
#  _                   _          __          _               _   
# (_)_ __  _ __  _   _| |_ ___   / /__  _   _| |_ _ __  _   _| |_ 
# | | '_ \| '_ \| | | | __/ __| / / _ \| | | | __| '_ \| | | | __|
# | | | | | |_) | |_| | |_\__ \/ / (_) | |_| | |_| |_) | |_| | |_ 
# |_|_| |_| .__/ \__,_|\__|___/_/ \___/ \__,_|\__| .__/ \__,_|\__|
#         |_|                                    |_|              
###################################################
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/ZZTo4L_Tune4C_13TeV-powheg-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/08E67F10-0069-E411-B1B1-00266CF9BEE4.root'
    )
)
# is the input file a miniAOD?
process.options = cms.untracked.PSet(isMiniAOD = cms.untracked.bool(False))
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions
process.GlobalTag.globaltag = cms.string('PHYS14_25_V1::All')

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('AC1B_test.root')
)

                                               
###################################################
#  ___  ___  __ _ _   _  ___ _ __   ___ ___  ___ 
# / __|/ _ \/ _` | | | |/ _ \ '_ \ / __/ _ \/ __|
# \__ \  __/ (_| | |_| |  __/ | | | (_|  __/\__ \
# |___/\___|\__, |\__,_|\___|_| |_|\___\___||___/
#              |_|                               
###################################################

# The Good vertices collection _____________________________________________||
process.goodVertices = cms.EDFilter(
    "VertexSelector",
    filter = cms.bool(False),
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
)
process.vertex_step=cms.Path(process.goodVertices)
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

######## Jet MET corrections
process.load('RecoJets.Configuration.RecoPFJets_cff')
#from RecoJets.Configuration.RecoPFJets_cff import kt6PFJets
#process.kt6PFJets = kt4PFJets.clone( rParam = 0.6 )
process.kt6PFJets.doRhoFastjet = True
process.kt6PFJets.doAreaFastjet = True
process.kt6PFJets.voronoiRfact = 0.9

#from JetMETCorrections.Configuration.JetCorrectionServices_cff import *
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('JetMETCorrections.Configuration.JetCorrectionServices_cff')
process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")
#process.load('JetMETCorrections.Type1MET.pfMETCorrections_cff')


process.pfJetMETcorr = cms.EDProducer("PFJetMETcorrInputProducer",
    src = cms.InputTag('ak4PFJets'),
    offsetCorrLabel = cms.string("ak4PFL1Fastjet"),
    jetCorrLabel = cms.string("ak4PFL1FastL2L3"), # NOTE: use "ak4PFL1FastL2L3" for MC / "ak4PFL1FastL2L3Residual" for Data
    jetCorrEtaMax = cms.double(9.9),
    type1JetPtThreshold = cms.double(10.0),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.90),
    skipMuons = cms.bool(True),
    skipMuonSelection = cms.string("isGlobalMuon | isStandAloneMuon")
)
process.pfType1CorrectedMet = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag('pfMet'),
    applyType0Corrections = cms.bool(False),
    srcCHSSums = cms.VInputTag(
        cms.InputTag('pfchsMETcorr', 'type0')
    ),
    type0Rsoft = cms.double(0.6),
    applyType1Corrections = cms.bool(True),
    srcType1Corrections = cms.VInputTag(
        cms.InputTag('pfJetMETcorr', 'type1')
    ),
    applyType2Corrections = cms.bool(False)
)  

#from JetMETCorrections.Type1MET.correctedMet_cff import *
process.pfJetMETcorr.jetCorrLabel = cms.string("ak4PFL1FastL2L3") # ? ak5PFL1FastL2L3") #MC

#from JetMETCorrections.Type1MET.pfMETCorrections_cff import pfType1CorrectedMet
process.pfType0Type1CorrectedMet = process.pfType1CorrectedMet.clone(
    applyType0Corrections = cms.bool(False),
    srcType1Corrections = cms.VInputTag(
        cms.InputTag('pfMETcorrType0'),
        cms.InputTag('pfJetMETcorr', 'type1')
    )
)
process.producePFMETCorrections=cms.Sequence(process.pfJetMETcorr*process.pfType1CorrectedMet)
process.metAnalysisSequence=cms.Sequence(process.type0PFMEtCorrection*process.producePFMETCorrections*process.pfType0Type1CorrectedMet)
#process.jet_step = cms.Path(process.kt6PFJets*process.metAnalysisSequence)

######PF ISO calculation for Electrons
from CommonTools.ParticleFlow.Isolation.pfElectronIsolation_cff import *
process.load('CommonTools.ParticleFlow.Isolation.pfElectronIsolation_cff')
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFPhotonIso
process.stdElectronSequencePFIso = setupPFElectronIso(process, 'gedGsfElectrons')#'gsfElectrons')
process.stdPhotonSequencePFIso = setupPFPhotonIso(process, 'photons')
process.pfiso_step = cms.Path( process.pfParticleSelectionSequence + process.stdElectronSequencePFIso + process.stdPhotonSequencePFIso)
##############################

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
    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
# save PAT Layer 1 output; you need a '*' to unpack the list of commands 'patEventContent'
    outputCommands = cms.untracked.vstring('drop *', *patEventContent )
)


from PhysicsTools.PatAlgos.tools.coreTools import *
#process.load("PhysicsTools.PatAlgos.patTemplate_cfg") #from PhysicsTools.PatAlgos.patTemplate_cfg import *
#removeAllPATObjectsBut(process, ['Jets'])
#from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT
from PhysicsTools.PatAlgos.tools.pfTools import *

postfix = "PFlow"

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")

usePF2PAT(process,runPF2PAT=True, jetAlgo='AK4', runOnMC=False, postfix=postfix,
    #jetCorrections=('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute']),
    jetCorrections=('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'],'Type-1'),
    pvCollection=cms.InputTag('offlinePrimaryVertices')
)
process.pfPileUpPFlow.checkClosestZVertex = False
getattr(process, "patPF2PATSequence"+postfix).remove(getattr(process, "cleanPatCandidates"+postfix))
getattr(process, "patPF2PATSequence"+postfix).remove(getattr(process, "countPatCandidates"+postfix))
process.patseq = cms.Sequence(getattr(process,"patPF2PATSequence"+postfix))
process.pat_step = cms.Path(process.patseq)

#######Pileup Jet ID (Hendrik)
##https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetID#Running_the_algorithm_on_reco_je
#
#from CMGTools.External.pujetidsequence_cff import puJetId, puJetMva
#
#process.recoPuJetId = puJetId.clone(
#		jets = cms.InputTag("ak4PFJets"),
#		applyJec = cms.bool(True),
#		inputIsCorrected = cms.bool(False),                
#		)
#
#process.recoPuJetMva = puJetMva.clone(
#		jets = cms.InputTag("ak4PFJets"),
#		jetids = cms.InputTag("recoPuJetId"),
#		applyJec = cms.bool(True),
#		inputIsCorrected = cms.bool(False),                
#		)
#process.recoPuJetIdSqeuence = cms.Sequence(process.recoPuJetId * process.recoPuJetMva)
#process.jetpuid_step = cms.Path(process.recoPuJetIdSqeuence)

#process.load("RecoJets.JetProducers.PileupJetID_cfi")
pileupJetIdProducer = cms.EDProducer('PileupJetIdProducer')
process.recoPuJetId = pileupJetIdProducer.clone(
    jets = cms.InputTag("ak4PFJets"),
    jetids = cms.InputTag("recoPuJetId"),
    applyJec = cms.bool(True),
    inputIsCorrected = cms.bool(False)
)

###################################################
#                  _                   _             
#  _ __ ___   ___ | |_ _ __ ___   __ _| | _____ _ __ 
# | '__/ _ \ / _ \| __| '_ ` _ \ / _` | |/ / _ \ '__|
# | | | (_) | (_) | |_| | | | | | (_| |   <  __/ |   
# |_|  \___/ \___/ \__|_| |_| |_|\__,_|_|\_\___|_|   
#                                                                        
###################################################
process.makeroottree = cms.EDAnalyzer("RootMaker",
    
    packedPfCands = cms.InputTag("packedPFCandidates"),
    # not used in AOD

    hcalNoiseInfo = cms.InputTag("hcalnoise", "", "RECO"),
    # same in MINI and AOD

# TRIGGER #####################################################
    triggerBits = cms.InputTag("TriggerResults","","HLT"),
    # same in mini and aod

    triggerPrescales = cms.InputTag("patTrigger"),
    # 

    triggerObjects = cms.InputTag("selectedPatTrigger"),
    #

    puInfo = cms.InputTag("addPileupInfo", "", "HLT"),
    # same in mini and aod

    metFilterBits = cms.InputTag("TriggerResults", "", "PAT"),

    HLTriggerSelection = cms.untracked.vstring(),
    Trigger = cms.untracked.bool(False),

# MUONS #######################################################
    #muons = cms.InputTag("slimmedMuons"),
    muons = cms.InputTag("muons"),

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

    #electrons = cms.InputTag("slimmedElectrons"),
    #photons = cms.InputTag("slimmedPhotons"),
    electrons = cms.InputTag("gedGsfElectrons"),
    photons = cms.InputTag("gedPhotons"),

    #ebRecHits = cms.InputTag("reducedEgamma", "reducedEBRecHits", "PAT"),
    #eeRecHits = cms.InputTag("reducedEgamma", "reducedEERecHits", "PAT"),
    #esRecHits = cms.InputTag("reducedEgamma", "reducedESRecHits", "PAT"),
    ebRecHits = cms.InputTag("reducedEcalRecHitsEB"),
    eeRecHits = cms.InputTag("reducedEcalRecHitsEB"),
    esRecHits = cms.InputTag("reducedEcalRecHitsES"),

    #ebeeClusters = cms.InputTag("reducedEgamma", "reducedEBEEClusters", "PAT"),
    #esClusters = cms.InputTag("reducedEgamma", "reducedESClusters", "PAT"),
    #conversions = cms.InputTag("reducedEgamma", "reducedConversions", "PAT"),
    ebeeClusters = cms.InputTag("EBEEClusters"),
    esClusters = cms.InputTag("ESClusters"),
    conversions = cms.InputTag("conversions"),

    singleLegConversions = cms.InputTag("reducedEgamma", "reducedSingleLegConversions", "PAT"),
    #gedGsfElectronCores = cms.InputTag("reducedEgamma", "reducedGedGsfElectronCores", "PAT"),
    #gedPhotonCores = cms.InputTag("reducedEgamma", "reducedGedPhotonCores", "PAT"),
    #singleLegConversions = cms.InputTag(""),
    gedGsfElectronCores = cms.InputTag("gedGsfElectronCores"),
    gedPhotonCores = cms.InputTag("gedPhotonCore"),

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
    #taus = cms.InputTag("slimmedTaus"),
    taus = cms.InputTag("taus"),

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
    #jets = cms.InputTag("slimmedJets"),
    jets = cms.InputTag("ak4PFJets"),
    #jetsAK8 = cms.InputTag("slimmedJetsAK8"),
    jetsAK8 = cms.InputTag("ak8PFJets"),

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
    genParticles = cms.InputTag("genParticles"),
    #genJets = cms.InputTag("slimmedGenJets", "", "PAT"),
    genJets = cms.InputTag("ak4GenJets"),

    genInfo = cms.InputTag("generator", "", "SIM"),

    GenAllParticles = cms.untracked.bool(False),
    GenSomeParticles = cms.untracked.bool(True),
    GenAK4Jets = cms.untracked.bool(True),

# MET #########################################################
    #mets = cms.InputTag("slimmedMETs"),
    mets = cms.InputTag("pfMet"),

    RecPFMet = cms.untracked.bool(True),
    
# TRACKS ######################################################
    lostTracks = cms.InputTag("lostTracks", "", "PAT"),
    generalTracks = cms.InputTag("generalTracks", "", "RECO"),
    # same in aod and mini

#    RecTrack = cms.untracked.bool(False),
#    RecTrackPtMin = cms.untracked.double(10.),
#    RecTrackEtaMax = cms.untracked.double(2.5),
#    RecTrackNum = cms.untracked.int32(100000),
#    RecTrackFilterPtMin = cms.untracked.double(18.),
    
# SUPERCLUSTER ################################################
    #superClusters = cms.InputTag("reducedEgamma", "reducedSuperClusters", "PAT"),
    superClusters = cms.InputTag(""),

    RecSuperCluster = cms.untracked.bool(True),
    RecSuperClusterFilterPtMin = cms.untracked.double(8.),
    RecSuperClusterBasicCluster = cms.untracked.bool(False),
    RecSuperClusterHit = cms.untracked.bool(False),

# VERTICES ####################################################
    #vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    vertices = cms.InputTag("offlinePrimaryVertices"),
    secondaryVertices = cms.InputTag("slimmedSecondaryVertices", "", "PAT"),

    beamSpot = cms.InputTag("offlineBeamSpot", "", "RECO"),
    # same in mini and aod

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


###################################################
#  _ __  _ __ ___   ___ ___  ___ ___ 
# | '_ \| '__/ _ \ / __/ _ \/ __/ __|
# | |_) | | | (_) | (_|  __/\__ \__ \
# | .__/|_|  \___/ \___\___||___/___/
# |_|                                
###################################################
process.p = cms.Path(
    process.goodVertices*
    process.eidHyperTight1MC*
    process.eidLooseMC*
    process.pfParticleSelectionSequence*
    process.stdElectronSequencePFIso*
    process.stdPhotonSequencePFIso*
    process.makeroottree
)
