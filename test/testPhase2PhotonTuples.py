import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.StandardSequences.Eras import eras

process = cms.Process("Ntupler", eras.Phase2)
options = VarParsing('analysis')
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load("RecoParticleFlow.PFClusterProducer.particleFlowRecHitHGC_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring(options.inputFiles) )

process.ntupler = cms.EDAnalyzer("Phase2PhotonTupler",
    photons = cms.InputTag("photonsFromMultiCl"),
    gedPhotons = cms.InputTag("gedPhotons"),
    genParticles = cms.InputTag("genParticles"),
    genCut = cms.string("pt>5 && status==1 && (abs(pdgId)==11 || pdgId==22)"),
    simClusters = cms.InputTag("mix:MergedCaloTruth"),
    caloParticles = cms.InputTag("mix:MergedCaloTruth"),
    doPremixContent = cms.bool(False),
    trackingParticles = cms.InputTag("mix:MergedTrackTruth"),
    trackingVertices = cms.InputTag("mix:MergedTrackTruth"),
    caloHits = cms.InputTag("g4SimHits:EcalHitsEB"),
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
    scEta = cms.string("superCluster().eta()"),
    trkSumPtSolidConeDR04 = cms.string("trkSumPtSolidConeDR04()"),
    nTrkSolidConeDR04 = cms.string("nTrkSolidConeDR04()"),
    full5x5_sigmaIetaIeta = cms.string("full5x5_sigmaIetaIeta()"),
    full5x5_sigmaIetaIphi = cms.string("full5x5_showerShapeVariables().sigmaIetaIphi"),
    full5x5_sigmaIphiIphi = cms.string("full5x5_showerShapeVariables().sigmaIphiIphi"),
    hadronicOverEm = cms.string("hadronicOverEm()"),
    hadronicDepth1OverEm = cms.string("hadronicDepth1OverEm()"),
    hadronicDepth2OverEm = cms.string("hadronicDepth2OverEm()"),
    hadTowOverEm = cms.string("hadTowOverEm()"),
    hadTowDepth1OverEm = cms.string("hadTowDepth1OverEm()"),
    hadTowDepth2OverEm = cms.string("hadTowDepth2OverEm()"),
    hasPixelSeed = cms.string("hasPixelSeed()"),
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
