from MyRootMaker.MyRootMaker.RootMakerTemplateMC_mini_cfg import *

process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load('RecoJets.Configuration.RecoJetAssociations_cff')
process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')
process.load('RecoBTag.Configuration.RecoBTag_cff')
process.load('RecoJets.Configuration.RecoJetAssociations_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

#process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

process.ak4JetTracksAssociatorAtVertexPF.jets = cms.InputTag("ak4PFJetsCHS")
process.ak4JetTracksAssociatorAtVertexPF.tracks = cms.InputTag("unpackedTracksAndVertices")
process.impactParameterTagInfos.primaryVertex = cms.InputTag("unpackedTracksAndVertices")
process.inclusiveSecondaryVertexFinderTagInfos.extSVCollection = cms.InputTag("unpackedTracksAndVertices","secondary","")
#process.combinedSecondaryVertex.trackMultiplicityMin = 1


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch//store/data/Run2015B/DoubleMuon/MINIAOD/PromptReco-v1/000/251/252/00000/9ADEE140-9C27-E511-919A-02163E011D23.root' # DoubleMuon, MINIAOD run 251252
        #'file:/afs/cern.ch/work/e/ekennedy/work/tuplizer/miniAOD/TTbarH_M-125_13TeV_mini_PU40bx25_PHYS14_25_V1_file2.root'
    )
)
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions
process.GlobalTag.globaltag = cms.string('74X_dataRun2_Prompt_v0')
process.makeroottree.isMiniAOD = cms.untracked.bool(True)




process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100) 
    #input = cms.untracked.int32(-1) 
)
#process.makeroottree.debug = cms.untracked.bool(True)

process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.TFileService = cms.Service("TFileService",
	fileName = cms.string('AC1B_test74_data_mini.root')
)

#process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))

process.makeroottree.RecMuonNum = cms.untracked.int32(0)
process.makeroottree.RecElectronNum = cms.untracked.int32(0)
process.makeroottree.RecAK4PFCHSNum = cms.untracked.int32(0)
process.makeroottree.HLTriggerSelection = cms.untracked.vstring()
process.makeroottree.GenAllParticles = cms.untracked.bool(True)
process.makeroottree.GenSomeParticles = cms.untracked.bool(True)
process.makeroottree.GenAK4Jets = cms.untracked.bool(True)

process.p = cms.Path(
    process.makeroottree
)

