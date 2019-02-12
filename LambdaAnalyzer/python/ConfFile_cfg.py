import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
#mylist = FileUtils.loadListFromFile('list2.txt')

process = cms.Process("lambda")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.option = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )



process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:myfile.root'
        #*mylist
        '/store/mc/RunIISummer16MiniAODv2/Vector_MonoJ_NLO_Mphi-2000_Mchi-400_gSM-0p25_gDM-1p0_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/E2E1930D-69D5-E611-B2B2-FA163E039A43.root'
    )
)

process.ntuple = cms.EDAnalyzer("LambdaAnayzer",
                                genSet = cms.PSet(
                                                  genProduct = cms.InputTag('generator'),
                                                  lheProduct = cms.InputTag('externalLHEProducer' if not isCustom else 'source'),
                                                  genParticles = cms.InputTag('prunedGenParticles'),
                                                  ),
                                jetSet = cms.PSet(
                                                  vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                                  jets = cms.untracked.InputTag('slimmedJets'),
                                                  met = cms.InputTag('slimmedMETs'),
                                                  ),
                                fatJetSet = cms.PSet(
                                                  jets = cms.InputTag('slimmedJetsAK8'),
                                                  ),
                                histFile = cms.string('%s/src/edmAna/LambdaAnayzer/data/HistList.dat' % os.environ['CMSSW_BASE']),
                                )


process.p = cms.Path(process.ntuple)
