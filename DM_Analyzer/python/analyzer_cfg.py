import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
mylist = FileUtils.loadListFromFile('list2.txt')

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.option = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use, run_01
    fileNames = cms.untracked.vstring(
        ### Classical vector monojet root file (Mphi500 - Mchi150) ###
        #'dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/jongho/DarkMatter/F4DA5CCC-57D1-E611-95AE-FA163ECB4BF2.root'
	#'/store/mc/RunIISummer16MiniAODv2/Vector_MonoJ_NLO_Mphi-2000_Mchi-400_gSM-0p25_gDM-1p0_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/FEE42DF9-8CD5-E611-B1A2-6CC2173C3E80.root',
        #'/store/mc/RunIISummer16MiniAODv2/Vector_MonoJ_NLO_Mphi-2000_Mchi-400_gSM-0p25_gDM-1p0_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/FA537B6D-79D5-E611-9CB2-0CC47A78A45A.root',
	#'/store/mc/RunIISummer16MiniAODv2/Vector_MonoJ_NLO_Mphi-2000_Mchi-400_gSM-0p25_gDM-1p0_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/E2E1930D-69D5-E611-B2B2-FA163E039A43.root'
	### Our monojet root file (MZprime500 - Mhs50 - Mchi150) ###
	#'dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/jongho/DarkMatter/temp/MonojetMatch1JetsDM_LO_MZprime-500_Mhs-50_Mchi-150_gSM-0p25_gDM-1p0_th_0p01_13TeV-madgraph_10024331_miniaod.root'
	*mylist
	#'dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/jongho/DarkMatter/temp/MonojetMatch1JetsDM_LO_MZprime-500_Mhs-50_Mchi-150_gSM-0p25_gDM-1p0_th_0p01_13TeV-madgraph_10024331_miniaod.root'
    )
)

process.demo = cms.EDAnalyzer('DM_Analyzer',
                              genParticleTag = cms.untracked.InputTag('prunedGenParticles'),
                              metTag = cms.untracked.InputTag('slimmedMETs')
)

#process.TFileService = cms.Service("TFileService", fileName = cms.string('test.root') )
#process.TFileService = cms.Service("TFileService", fileName = cms.string('MonoJClassical_Mphi_2000_Mchi_400.root') )
process.TFileService = cms.Service("TFileService", fileName = cms.string('MonoJdarkHiggs_Mphi_2000_Mhs_50_Mchi_400.root') )


process.p = cms.Path(process.demo)
