from MyRootMaker.MyRootMaker.RootMakerTemplateMC_mini_cfg import *

process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load('RecoJets.Configuration.RecoJetAssociations_cff')
process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')
process.load('RecoBTag.Configuration.RecoBTag_cff')
process.load('RecoJets.Configuration.RecoJetAssociations_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
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
#        'root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/120BB291-6912-E511-82D9-008CFA064770.root',
#        'root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/982C788D-6912-E511-A1D6-008CFA1980B8.root',
#        'root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/7AB7A58E-6912-E511-AAF3-549F35AD8BE3.root'
 #       'root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/026B146D-8716-E511-B16D-00266CF9BDFC.root',
        'root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/08223373-8716-E511-A47C-00266CF9AFF0.root'
#        'root://cms-xrd-global.cern.ch//store/mc/Spring14miniaod/ZZTo4L_Tune4C_13TeV-powheg-pythia8/MINIAODSIM/141029_PU40bx50_PLS170_V6AN2-v1/00000/0649812F-EB6A-E411-90DE-00259094F2E8.root'
        #'file:/afs/cern.ch/work/e/ekennedy/work/tuplizer/miniAOD/TTbarH_M-125_13TeV_mini_PU40bx25_PHYS14_25_V1_file2.root'
    )
)
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions
#process.GlobalTag.globaltag = cms.string('PHYS14_25_V1::All')
#process.GlobalTag.globaltag = cms.string('PLS170_V6AN2::All')
process.GlobalTag.globaltag = cms.string('MCRUN2_74_V9')
process.makeroottree.isMiniAOD = cms.untracked.bool(True)




process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000) 
    #input = cms.untracked.int32(-1) 
)
#process.makeroottree.debug = cms.untracked.bool(True)

process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.TFileService = cms.Service("TFileService",
	fileName = cms.string('AC1B_test74_mini.root')
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
    process.egmGsfElectronIDSequence *
    process.makeroottree
)

