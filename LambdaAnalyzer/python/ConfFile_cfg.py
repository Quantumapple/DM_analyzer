import os
import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
#mylist = FileUtils.loadListFromFile('axial_samplelist.txt')
#mylist = FileUtils.loadListFromFile('vector_samplelist.txt')
#mylist = FileUtils.loadListFromFile('darkHiggsMonoJ_samplelist.txt')
mylist = FileUtils.loadListFromFile('subset.txt')

process = cms.Process("lambda")

process.load("FWCore.MessageService.MessageLogger_cfi")

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.option = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:myfile.root'
        *mylist
        #'/store/mc/RunIISummer16MiniAODv2/Axial_MonoJ_NLO_Mphi-500_Mchi-150_gSM-0p25_gDM-1p0_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/74EB2BE6-D5D4-E611-8A99-1CC1DE1CF44E.root'
        #'root://cmseos.fnal.gov://store/user/jongho/zprime/GenProduction/MonojetMatch1JetsDM_LO_MZprime_500_Mhs_50_Mchi_150_gSM_0p25_gDM_1p0_th_0p01_13TeV_madgraph/MonojetMatch1JetsDM_LO_MZprime-500_Mhs-50_Mchi-150_gSM-0p25_gDM-1p0_th_0p01_13TeV-madgraph_27686463_miniaod.root'
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
                              histFile = cms.string('%s/src/EDMnalyzer/LambdaAnalyzer/data/HistList.dat' % os.environ['CMSSW_BASE']),
                              )

process.TFileService = cms.Service("TFileService",
                                    #fileName = cms.string("histo.root"),
                                    #fileName = cms.string("AxialMonoJ_histo.root"),
                                    fileName = cms.string("SUBSET.root"),
                                    #fileName = cms.string("DarkHiggs_histo.root"),
                                    closeFileFast = cms.untracked.bool(True)
                                    )

process.p = cms.Path(process.lamb)
