from MyRootMaker.MyRootMaker.RootMakerTemplateMC_cfg import *

# if using crab, doesn't matter what this is, input is set in crab.cfg
# if using cmsRun, this is used for input
process.source.fileNames = cms.untracked.vstring(
    #'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/0456F76A-EA77-E411-81F8-0025B31E3D3C.root', # 2600 events
    #'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/10C19613-E176-E411-92F7-F04DA23BBCCA.root', # 4200 events
    #'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/1CA79CF9-F876-E411-BD55-001E67396E05.root', # 6400 events
    #'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/4059ABE0-E676-E411-81E5-002590A4C6BE.root', # 3900 events
    #'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/422588BA-F376-E411-82B1-001E67397AE4.root', # 6100 events
    #'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/46324C48-F076-E411-8143-001E67396892.root', # 6500 events
    #'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/507A8AF3-D276-E411-87E5-002590A370B2.root', # 4000 events
    #'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/68F5DFAD-FD76-E411-AAE0-002590200AE0.root', # 4600 events 
    #'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/6C4F596C-E576-E411-B749-002590A4C69A.root', # 4900 events
    'root://cms-xrd-global.cern.ch//store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/AODSIM/PU40bx25_PHYS14_25_V1-v1/00000/6E125622-FF76-E411-A4CC-001E67397747.root' # 7400 events
    #'root://cms-xrd-global.cern.ch/'
)

process.GlobalTag.globaltag = cms.string('PHYS14_25_V1::All')

#process.source.skipEvents=cms.untracked.uint32(290)
process.MessageLogger.cerr.FwkReport.reportEvery = 10

# also overridden by crab.cfg
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
#process.makeroottree.debug = cms.untracked.bool(True)



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

# DEBUGGING
#process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))
#process.makeroottree.RecMuon = cms.untracked.bool(False)
#process.makeroottree.RecElectron = cms.untracked.bool(False)
#process.makeroottree.RecPhoton = cms.untracked.bool(False)
#process.makeroottree.RecTau = cms.untracked.bool(False)
#process.makeroottree.RecAK4PFCHSJet = cms.untracked.bool(False)
process.makeroottree.RecAK4CaloJet = cms.untracked.bool(False)
process.makeroottree.RecAK4JPTJet = cms.untracked.bool(False)
#process.makeroottree.RecAK4PFJet = cms.untracked.bool(False)
#process.makeroottree.RecAK4PFCHSJet = cms.untracked.bool(False)    
#process.makeroottree.RecSecVertices = cms.untracked.bool(False)

process.makeroottree.RecMuonNum = cms.untracked.int32(0)
process.makeroottree.HLTriggerSelection = cms.untracked.vstring()
#process.patJetCorrFactors.levels=cms.vstring('L1FastJet','L2Relative', 'L3Absolute')
process.makeroottree.GenAllParticles = cms.untracked.bool(False)
process.makeroottree.GenSomeParticles = cms.untracked.bool(True)
process.makeroottree.GenAK4Jets = cms.untracked.bool(True)

process.makeroottree.isMiniAOD = cms.untracked.bool(False)
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


