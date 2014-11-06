from MyRootMaker.MyRootMaker.RootMakerTemplateMC_cfg import *

# uses this set for a test run
process.source.fileNames = cms.untracked.vstring(
'root://cms-xrd-global.cern.ch//store/mc/Spring14dr/GluGluToHToMuMu_M-125_13TeV-powheg-pythia6/AODSIM/PU_S14_POSTLS170_V6-v1/00000/1E93B8DB-7CFD-E311-BFC7-7845C4FC346A.root'
)

process.maxEvents = cms.untracked.PSet(
#    input = cms.untracked.int32(-1)
    input = cms.untracked.int32(100)
)

process.GlobalTag.globaltag = cms.string('START70_V7::All')

#process.GlobalTag.toGet = cms.VPSet(
#		cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
#			tag = cms.string("TrackProbabilityCalibration_2D_MC53X_v2"),
#			connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
#		cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
#			tag = cms.string("TrackProbabilityCalibration_3D_MC53X_v2"),
#			connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
#		)

process.load('RecoBTag/Configuration/RecoBTag_cff')
process.btag = cms.Path(process.btagging)

process.makeroottree.RecMuonNum = cms.untracked.int32(0)
process.makeroottree.HLTriggerSelection = cms.untracked.vstring()
process.patJetCorrFactors.levels=cms.vstring('L1FastJet','L2Relative', 'L3Absolute')
process.makeroottree.GenAllParticles = cms.untracked.bool(True)
process.makeroottree.GenSomeParticles = cms.untracked.bool(False)
process.makeroottree.GenAK5Jets = cms.untracked.bool(True)

process.schedule = cms.Schedule(
process.vertex_step,
process.filters_step,
process.jet_step,
#process.jetpuid_step,
process.btag,
#process.pat_step,
process.electron_step,
process.jetflavour_step,
process.pfiso_step,
process.roottree_step)

process.options.allowUnscheduled = cms.untracked.bool( True )
