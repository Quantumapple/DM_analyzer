import os
import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
#mylist = FileUtils.loadListFromFile('mjetfalsev2.txt')

process = cms.Process("lambda")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(15) )

process.option = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #*mylist
       #'/DarkHiggs_MonoHs_LO_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-rp_102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM'
       #'/store/mc/RunIIAutumn18MiniAOD/DarkHiggs_MonoHs_LO_TuneCP5_13TeV-madgraph-pythia8/MINIAODSIM/rp_102X_upgrade2018_realistic_v15_ext1-v1/40000/C9868F55-0253-604E-ADF7-D76128BADAAE.root'
       '/store/mc/RunIIAutumn18MiniAOD/DarkHiggs_MonoHs_LO_TuneCP5_13TeV-madgraph-pythia8/MINIAODSIM/rp_102X_upgrade2018_realistic_v15_ext1-v1/40000/FA8F9AF5-CECD-C04D-B9CC-72FA84A35725.root'
    )
)


process.lamb = cms.EDAnalyzer("LambdaAnalyzer",
                              genSet = cms.PSet(
                                                genProduct = cms.InputTag('generator'),
                                                lheProduct = cms.InputTag('externalLHEProducer'),
                                                genParticles = cms.InputTag('prunedGenParticles'),
                                                pdgId = cms.vint32(1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 21, 23, 24, 25, 52, 55, 54), #52 DM; 55 Zp; 54 hs 
                                                ),
                              electronSet = cms.PSet(
                                                electrons = cms.InputTag('slimmedElectrons'),
                                                vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                                ),
                              muonSet = cms.PSet(
                                                muons = cms.InputTag('slimmedMuons'),
                                                ),
                              jetSet = cms.PSet(
                                                vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                                jets = cms.InputTag('slimmedJets'),
                                                met = cms.InputTag('slimmedMETs'),
                                                genjets = cms.InputTag('slimmedGenJets'),
                                                ),
                              #fatJetSet = cms.PSet(
                              #                  jets = cms.InputTag('slimmedJetsAK8'),
                              #                  ),
                              histFile = cms.string('%s/src/DM_analyzer/LambdaAnalyzer/data/HistList.dat' % os.environ['CMSSW_BASE']),
                              )

process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string("offical_test.root"),
                                    closeFileFast = cms.untracked.bool(True)
                                    )

process.p = cms.Path(process.lamb)
