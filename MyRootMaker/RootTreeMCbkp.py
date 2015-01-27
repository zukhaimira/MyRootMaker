from MyRootMaker.MyRootMaker.RootMakerTemplateMC_cfg import *

# if using crab, doesn't matter what this is, input is set in crab.cfg
# if using cmsRun, this is used for input
process.source.fileNames = cms.untracked.vstring(
    #'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/ZZTo4L_Tune4C_13TeV-powheg-pythia8/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/08E67F10-0069-E411-B1B1-00266CF9BEE4.root'
    #'root://cms-xrd-global.cern.ch//store/mc/Spring14dr/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/AODSIM/PU_S14_POSTLS170_V6-v1/00000/003B7873-00F3-E311-8F81-0025905A48F2.root'
    #'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/DYJetsToEEMuMu_M-120To200_13TeV-madgraph/AODSIM/PU20bx25_PHYS14_25_V1-v2/10000/6E4F7A3D-9B7C-E411-9A04-00259073E34C.root'
    #'root://cms-xrd-global.cern.ch//store/mc/Spring14dr/WH_ZH_HToZZ_4LFilter_M-125_13TeV_pythia6/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/1CCAA6A7-D1C8-E311-BA70-002590A371CA.root'
    #'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/WH_ZH_HToGG_M-125_13TeV_pythia6/AODSIM/PU20bx25_PHYS14_25_V1-v1/00000/2A674041-3E69-E411-A73B-001E67397021.root'
    'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/00000/20329E54-327F-E411-B47B-001E673972E7.root'
)

# also overridden by crab.cfg
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

#process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))

#process.GlobalTag.globaltag = cms.string('POSTLS162_V2::All')
#process.GlobalTag.globaltag = cms.string('POSTLS170_V5::All')
process.GlobalTag.globaltag = cms.string('PHYS14_25_V1::All')

#from Configuration.PyReleaseValidation.autoCond import autoCond
#process.GlobalTag.globaltag = autoCond['startup']

process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
        tag = cms.string("TrackProbabilityCalibration_2D_MC53X_v2"),
        connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
    cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
        tag = cms.string("TrackProbabilityCalibration_3D_MC53X_v2"),
        connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
)

# No data of type "JetTagComputer" with label "negativeTrackCounting3D2nd" in record "JetTagComputerRecord"
# process.load("RecoBTag.ImpactParameter.negativeTrackCounting3D2ndComputer_cfi")
# No data of type "JetCorrectorParametersCollection" with label "AK4PF" in record "JetCorrectionsRecord"
#process.load("RecoJets.JetProducers.ak4PFJets_cfi")

#from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
#process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff")



process.load('RecoBTag/Configuration/RecoBTag_cff')
process.btag = cms.Path(process.btagging)

process.makeroottree.RecMuonNum = cms.untracked.int32(0)
process.makeroottree.HLTriggerSelection = cms.untracked.vstring()
#process.patJetCorrFactors.levels=cms.vstring('L1FastJet','L2Relative', 'L3Absolute')
process.makeroottree.GenAllParticles = cms.untracked.bool(True)
process.makeroottree.GenSomeParticles = cms.untracked.bool(False)
process.makeroottree.GenAK4Jets = cms.untracked.bool(True)

process.schedule = cms.Schedule(
process.vertex_step,
process.filters_step,
process.jet_step,
process.jetpuid_step,
process.btag,
#process.pat_step,
process.electron_step,
process.jetflavour_step,
process.pfiso_step,
process.roottree_step)

process.options.allowUnscheduled = cms.untracked.bool( True )


