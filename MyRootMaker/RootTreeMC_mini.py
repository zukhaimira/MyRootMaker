from MyRootMaker.MyRootMaker.RootMakerTemplateMC_cfg import *

# if using crab, doesn't matter what this is, input is set in crab.cfg
# if using cmsRun, this is used for input
process.source.fileNames = cms.untracked.vstring(
    'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/00000/20329E54-327F-E411-B47B-001E673972E7.root'
)

# also overridden by crab.cfg
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

#process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))

#process.GlobalTag.globaltag = cms.string('POSTLS170_V5::All')
process.GlobalTag.globaltag = cms.string('PHYS14_25_V1::All')

process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
        tag = cms.string("TrackProbabilityCalibration_2D_MC53X_v2"),
        connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
    cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
        tag = cms.string("TrackProbabilityCalibration_3D_MC53X_v2"),
        connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
)

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


