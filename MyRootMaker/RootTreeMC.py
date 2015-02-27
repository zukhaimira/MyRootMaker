from MyRootMaker.MyRootMaker.RootMakerTemplateMC_cfg import *

# if using crab, doesn't matter what this is, input is set in crab.cfg
# if using cmsRun, this is used for input
process.source.fileNames = cms.untracked.vstring(
    'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/0456F76A-EA77-E411-81F8-0025B31E3D3C.root', # 2600 events
    #'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/10C19613-E176-E411-92F7-F04DA23BBCCA.root', # 4200 events
    #'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/1CA79CF9-F876-E411-BD55-001E67396E05.root' # 6400 events
    #'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/DYJetsToEEMuMu_M-120To200_13TeV-madgraph/AODSIM/PU20bx25_PHYS14_25_V1-v2/10000/062D6736-9F7C-E411-B200-E0CB4E19F9BD.root'
)

# also overridden by crab.cfg
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)
debug = cms.untracked.bool(True)

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


