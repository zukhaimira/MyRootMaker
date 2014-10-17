from MyRootMaker.MyRootMaker.RootMakerTemplateMC_cfg import *

# uses this set for a test run
process.source.fileNames = cms.untracked.vstring(
#'/GluGluToHToMuMu_M-125_13TeV-powheg-pythia6/Spring14dr-PU20bx25_POSTLS170_V5-v1/AODSIM'
#'/store/mc/Summer12_DR53X/GluGluToHToZG_M-125_8TeV-powheg-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/DEF04071-6EFA-E111-BA18-00266CFFC4D4.root'
'/store/mc/Summer12_DR53X/ZZTo2e2mu_8TeV_ext-powheg-pythia6/AODSIM/PU_S10_START53_V7C-v1/20000/00902F86-D46A-E211-851D-002618943970.root'
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
#process.btag = cms.Path(process.btagging)

process.makeroottree.RecMuonNum = cms.untracked.int32(0)
process.makeroottree.HLTriggerSelection = cms.untracked.vstring()
process.patJetCorrFactors.levels=cms.vstring('L1FastJet','L2Relative', 'L3Absolute')
#process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual") #DATA
#process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3")#MC
#idk if these are needed ??????
process.makeroottree.GenAllParticles = cms.untracked.bool(True)
process.makeroottree.GenSomeParticles = cms.untracked.bool(False)
process.makeroottree.GenAK5Jets = cms.untracked.bool(True)

# load the standard PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")


# load the coreTools of PAT
#from PhysicsTools.PatAlgos.tools.electronTools import *
#addElectronUserIsolation(process)

from PhysicsTools.PatAlgos.tools.jetTools import *
#(process.btestalt,labels)  =runBTagging(process, cms.InputTag('ak5PFJets'),"AOD","alt")

process.schedule = cms.Schedule(
##process.vertex_step,
##process.filters_step,
#process.jet_step,
#process.jetpuid_step,
#process.btag,
#process.pat_step,
#process.electron_step,
##process.jetflavour_step,
#process.pfiso_step,
process.roottree_step)


# added 14 october from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATRecipes#CMSSW_7_0_X_dev2014
process.options.allowUnscheduled = cms.untracked.bool( True )
#process.convertToUnscheduled(process)

