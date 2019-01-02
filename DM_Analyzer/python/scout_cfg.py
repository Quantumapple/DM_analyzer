import FWCore.ParameterSet.Config as cms

process = cms.Process("SCOUT")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.option = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:matching/Hadronized_roots/matching_phim_100_mchi_1_no_cut_ptllmin_100.root'
        'file:samples/hadronized_files/had_match_mphi_10_mchi_1.root'
    )
)

process.SCOUT = cms.EDAnalyzer('scout',
                                   genParticleTag = cms.untracked.InputTag('genParticles')
                                   )

#process.TFileService = cms.Service("TFileService", fileName = cms.string('scout.root') )

process.p = cms.Path(process.SCOUT)
