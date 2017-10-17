import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.StandardSequences.Eras import eras

process = cms.Process("Ntupler", eras.Phase2)
options = VarParsing('analysis')
options.parseArguments()

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring(options.inputFiles) )

# phoSrc = cms.InputTag("photons")
phoSrc = cms.InputTag("photonsFromMultiCl")

process.load("EgammaTools.EgammaAnalysis.HGCalPhotonIDValueMap_cfi")
process.HGCalPhotonIDValueMap.photons = phoSrc

process.ntupler = cms.EDAnalyzer("Phase2PhotonTupler",
    photons = phoSrc,
    gedPhotons = cms.InputTag("gedPhotons"),
    genParticles = cms.InputTag("genParticles"),
    genCut = cms.string("pt>5 && status==1 && pdgId==22"),
    simClusters = cms.InputTag("mix:MergedCaloTruth"),
    caloParticles = cms.InputTag("mix:MergedCaloTruth"),
    simTracksSrc = cms.InputTag("g4SimHits"),
    simVerticesSrc = cms.InputTag("g4SimHits"),
    rhoSrc = cms.InputTag("fixedGridRhoFastjetAll"),
    doPremixContent = cms.bool(False),
    trackingParticles = cms.InputTag("mix:MergedTrackTruth"),
    trackingVertices = cms.InputTag("mix:MergedTrackTruth"),
    caloHits = cms.InputTag("g4SimHits:EcalHitsEB"),
)

process.ntupler.gedRecoMisc = cms.PSet(
    scEta = cms.string("superCluster().eta()"),
    trkSumPtSolidConeDR04 = cms.string("trkSumPtSolidConeDR04()"),
    nTrkSolidConeDR04 = cms.string("nTrkSolidConeDR04()"),
    full5x5_sigmaIetaIeta = cms.string("full5x5_sigmaIetaIeta()"),
    full5x5_sigmaIetaIphi = cms.string("full5x5_showerShapeVariables().sigmaIetaIphi"),
    full5x5_sigmaIphiIphi = cms.string("full5x5_showerShapeVariables().sigmaIphiIphi"),
    full5x5_maxEnergyXtal = cms.string("full5x5_maxEnergyXtal()"),
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

process.ntupler.localRecoMisc = cms.PSet(
    process.ntupler.gedRecoMisc,
    sigmaUU = cms.InputTag("HGCalPhotonIDValueMap:sigmaUU"),
    sigmaVV = cms.InputTag("HGCalPhotonIDValueMap:sigmaVV"),
    sigmaEE = cms.InputTag("HGCalPhotonIDValueMap:sigmaEE"),
    sigmaPP = cms.InputTag("HGCalPhotonIDValueMap:sigmaPP"),
    nLayers = cms.InputTag("HGCalPhotonIDValueMap:nLayers"),
    firstLayer = cms.InputTag("HGCalPhotonIDValueMap:firstLayer"),
    lastLayer = cms.InputTag("HGCalPhotonIDValueMap:lastLayer"),
    energyEE = cms.InputTag("HGCalPhotonIDValueMap:energyEE"),
    energyFH = cms.InputTag("HGCalPhotonIDValueMap:energyFH"),
    energyBH = cms.InputTag("HGCalPhotonIDValueMap:energyBH"),
    measuredDepth = cms.InputTag("HGCalPhotonIDValueMap:measuredDepth"),
    expectedDepth = cms.InputTag("HGCalPhotonIDValueMap:expectedDepth"),
    expectedSigma = cms.InputTag("HGCalPhotonIDValueMap:expectedSigma"),
    depthCompatibility = cms.InputTag("HGCalPhotonIDValueMap:depthCompatibility"),
    caloIsoRing0 = cms.InputTag("HGCalPhotonIDValueMap:caloIsoRing0"),
    caloIsoRing1 = cms.InputTag("HGCalPhotonIDValueMap:caloIsoRing1"),
    caloIsoRing2 = cms.InputTag("HGCalPhotonIDValueMap:caloIsoRing2"),
    caloIsoRing3 = cms.InputTag("HGCalPhotonIDValueMap:caloIsoRing3"),
    caloIsoRing4 = cms.InputTag("HGCalPhotonIDValueMap:caloIsoRing4"),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
)

process.p = cms.Path(
    process.HGCalPhotonIDValueMap  +
    process.ntupler
)
