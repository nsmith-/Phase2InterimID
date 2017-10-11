import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.StandardSequences.Eras import eras

process = cms.Process("Ntupler", eras.Phase2)
options = VarParsing('analysis')
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load("RecoParticleFlow.PFClusterProducer.particleFlowRecHitHGC_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring(options.inputFiles) )

process.ntupler = cms.EDAnalyzer("Phase2ElectronTupler",
    ecalDrivenElectrons = cms.InputTag("ecalDrivenGsfElectrons"),
    gedGsfElectrons = cms.InputTag("gedGsfElectrons"),
    genParticles = cms.InputTag("genParticles"),
    simClusters = cms.InputTag("mix:MergedCaloTruth"),
    caloParticles = cms.InputTag("mix:MergedCaloTruth"),
    trackingParticles = cms.InputTag("mix:MergedTrackTruth"),
    HGCalIDToolConfig = cms.PSet(
        HGCBHInput = cms.InputTag("HGCalRecHit","HGCHEBRecHits"),
        HGCEEInput = cms.InputTag("HGCalRecHit","HGCEERecHits"),
        HGCFHInput = cms.InputTag("HGCalRecHit","HGCHEFRecHits"),
        HGCPFRecHits = cms.InputTag("particleFlowRecHitHGC::Ntupler"),
        withPileup = cms.bool(True),
        debug = cms.bool(False),
    ),
)
process.ntupler.localRecoMisc = cms.PSet(
    full5x5_sigmaIetaIeta = cms.string("full5x5_sigmaIetaIeta()"),
    hadTowOverEm = cms.string("hadTowOverEm()"),
    chargedHadronIso = cms.string("chargedHadronIso()"),
    neutralHadronIso = cms.string("neutralHadronIso()"),
    photonIso = cms.string("photonIso()"),
    sumPUPt = cms.string("sumPUPt()"),
)
process.ntupler.gedRecoMisc = cms.PSet(
    process.ntupler.localRecoMisc
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
)

process.p = cms.Path(
    #process.particleFlowRecHitHGCSeq+
    process.ntupler
)
