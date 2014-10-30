from MyRootMaker.MyRootMaker.RootMakerTemplateMC_Muons_cfg import *
import FWCore.ParameterSet.Config as cms
process.source.fileNames = cms.untracked.vstring(
'root://cmsxrootd.fnal.gov//store/mc/Spring14dr/GluGluToHToMuMu_M-125_13TeV-powheg-pythia6/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/14646D22-DDFD-E311-BB64-7845C4FC3C4D.root'
)

process.GlobalTag.globaltag = cms.string('START70_V7::All')

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

#process.makeroottree.GenAllParticles = cms.untracked.bool(True)
#process.makeroottree.GenSomeParticles = cms.untracked.bool(False)
#process.makeroottree.GenAK5Jets = cms.untracked.bool(True)
process.makeroottree.RecMuonNum = cms.untracked.int32(0)
process.makeroottree.HLTriggerSelection = cms.untracked.vstring()

#process.patJetCorrFactors.levels=cms.vstring('L1FastJet','L2Relative', 'L3Absolute')
#process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual") #DATA
#process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3")#MC

process.schedule = cms.Schedule(
#process.vertex_step,
#process.filters_step,
#process.jet_step,
#process.jetpuid_step,
#process.btag,
#process.pat_step,
#process.electron_step,
#process.jetflavour_step,
#process.pfiso_step,
process.roottree_step
)

# added 14 october from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATRecipes#CMSSW_7_0_X_dev2014
process.options.allowUnscheduled = cms.untracked.bool( True )
#process.options.SkipEvent = cms.untracked.vstring('ProductNotFound')
#FWCore.ParameterSet.Utilities.convertToUnscheduled(ROOTMAKER)
