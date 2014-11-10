from MyRootMaker.MyRootMaker.RootMakerTemplateMC_cfg import *

process.source.fileNames = cms.untracked.vstring(
#'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/GluGluToHToZG_M-125_8TeV-powheg-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/DEF04071-6EFA-E111-BA18-00266CFFC4D4.root'
'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU25bx50_START53_V19D-v1/20000/02C8446A-7FCD-E211-8E9F-0025B3E06394.root'
#'root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU25bx25_START53_V19D-v1/20000/0C8546EA-F0CC-E211-A809-002481E154CE.root'
)
process.maxEvents = cms.untracked.PSet(
   input = cms.untracked.int32(1100)
)

process.GlobalTag.globaltag = cms.string('START53_V7G::All')

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

process.makeroottree.GenAllParticles = cms.untracked.bool(True)
process.makeroottree.GenSomeParticles = cms.untracked.bool(False)
process.makeroottree.GenAK5Jets = cms.untracked.bool(True)
process.makeroottree.RecMuonNum = cms.untracked.int32(0)
process.makeroottree.HLTriggerSelection = cms.untracked.vstring()

process.schedule = cms.Schedule(
process.vertex_step,
process.filters_step,
process.jet_step,
process.jetpuid_step,
process.btag,
process.pat_step,
process.electron_step,
#process.jetflavour_step,
process.pfiso_step,
process.roottree_step)
