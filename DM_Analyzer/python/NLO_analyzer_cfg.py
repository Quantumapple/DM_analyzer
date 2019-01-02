import FWCore.ParameterSet.Config as cms

process = cms.Process("DEMO")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.option = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use, run_01
    fileNames = cms.untracked.vstring(
        #'file:hadronized_DMScalar_mphi_10_mchi_1.root'
        'file:/afs/cern.ch/work/s/shoh/analysis/bb_MET_Analysis_13TeV/signal_generation/CMSSW_7_1_24/src/DMScalar_run_01_mphi_10_mchi_1_NLO.root' 
        #'file:/afs/cern.ch/work/s/shoh/analysis/bb_MET_Analysis_13TeV/signal_generation/CMSSW_7_1_24/src/DMScalar_run_02_mphi_100_mchi_1_NLO.root'
        #'file:/afs/cern.ch/work/s/shoh/analysis/bb_MET_Analysis_13TeV/signal_generation/CMSSW_7_1_24/src/DMScalar_run_03_mphi_500_mchi_1_NLO.root'
    )
)

process.demo = cms.EDAnalyzer('DM_Analyzer',
                              genParticleTag = cms.untracked.InputTag('genParticles')
)

process.TFileService = cms.Service("TFileService", fileName = cms.string('v2_iteration/SAVE_DMScalar_run_01_mphi_10_mchi_1_NLO.root') )
#process.TFileService = cms.Service("TFileService", fileName = cms.string('v2_iteration/SAVE_DMScalar_run_02_mphi_100_mchi_1_NLO.root') )
#process.TFileService = cms.Service("TFileService", fileName = cms.string('v2_iteration/SAVE_DMScalar_run_03_mphi_500_mchi_1_NLO.root') )

process.p = cms.Path(process.demo)
