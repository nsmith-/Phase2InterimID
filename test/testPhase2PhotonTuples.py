import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.StandardSequences.Eras import eras

process = cms.Process("Ntupler", eras.Phase2)
options = VarParsing('analysis')
options.parseArguments()

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    numberOfThreads = cms.untracked.uint32(4),
    numberOfStreams = cms.untracked.uint32(0),
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring(options.inputFiles) )

# phoSrc = cms.InputTag("photons")
phoSrc = cms.InputTag("photonsFromMultiCl")

process.load("EgammaTools.EgammaAnalysis.HGCalPhotonIDValueMap_cfi")
process.HGCalPhotonIDValueMap.photons = phoSrc

process.ntupler = cms.EDAnalyzer("Phase2PhotonTupler",
    photons = phoSrc,
    localRecoCut = cms.string("pt>10 && !isEB"),
    gedRecoCut = cms.string("pt>10"),
    gedPhotons = cms.InputTag("gedPhotons"),
    genParticles = cms.InputTag("genParticles"),
    genCut = cms.string("pt>5 && status==1 && (abs(pdgId)==11 || pdgId==22)"),
    simClusters = cms.InputTag("mix:MergedCaloTruth"),
    caloParticles = cms.InputTag("mix:MergedCaloTruth"),
    simTracksSrc = cms.InputTag("g4SimHits"),
    simVerticesSrc = cms.InputTag("g4SimHits"),
    rhoSrc = cms.InputTag("fixedGridRhoFastjetAll"),
    ecalDrivenElectrons = cms.InputTag("ecalDrivenGsfElectrons"),
    gedGsfElectrons = cms.InputTag("gedGsfElectrons"),
    conversions = cms.InputTag("conversions"),
    beamspot = cms.InputTag("offlineBeamSpot"),
    doPremixContent = cms.bool(False),
    trackingParticles = cms.InputTag("mix:MergedTrackTruth"),
    trackingVertices = cms.InputTag("mix:MergedTrackTruth"),
    caloHits = cms.InputTag("g4SimHits:EcalHitsEB"),
)

common = cms.PSet(
    scEta = cms.string("superCluster().eta()"),
    scRawEnergy = cms.string("superCluster().rawEnergy()"),
    hasPixelSeed = cms.string("hasPixelSeed()"),
)

process.ntupler.gedRecoMisc = cms.PSet(
    common,
    full5x5_sigmaIetaIeta = cms.string("full5x5_sigmaIetaIeta()"),
    full5x5_sigmaIetaIphi = cms.string("full5x5_showerShapeVariables().sigmaIetaIphi"),
    full5x5_sigmaIphiIphi = cms.string("full5x5_showerShapeVariables().sigmaIphiIphi"),
    full5x5_maxEnergyXtal = cms.string("full5x5_maxEnergyXtal()"),
    etaWidth = cms.string("superCluster().etaWidth()"),
    phiWidth = cms.string("superCluster().phiWidth()"),
    r9 = cms.string("r9()"),
    s4 = cms.string("showerShapeVariables().e2x2/showerShapeVariables().e5x5"),
    full5x5_r9 = cms.string("full5x5_r9()"),
    full5x5_s4 = cms.string("full5x5_showerShapeVariables().e2x2/full5x5_showerShapeVariables().e5x5"),
    hadronicOverEm = cms.string("hadronicOverEm()"),
    hadronicDepth1OverEm = cms.string("hadronicDepth1OverEm()"),
    hadronicDepth2OverEm = cms.string("hadronicDepth2OverEm()"),
    hadTowOverEm = cms.string("hadTowOverEm()"),
    hadTowDepth1OverEm = cms.string("hadTowDepth1OverEm()"),
    hadTowDepth2OverEm = cms.string("hadTowDepth2OverEm()"),
    chargedHadronIso = cms.string("chargedHadronIso()"),
    neutralHadronIso = cms.string("neutralHadronIso()"),
    photonIso = cms.string("photonIso()"),
    trkSumPtSolidConeDR04 = cms.string("trkSumPtSolidConeDR04()"),
    nTrkSolidConeDR04 = cms.string("nTrkSolidConeDR04()"),
)

process.ntupler.localRecoMisc = cms.PSet(
    common,
    seed_det = cms.string("superCluster().seed().hitsAndFractions().at(0).first.det()"),
    seed_subdet = cms.string("superCluster().seed().hitsAndFractions().at(0).first.subdetId()"),
)
for key in process.HGCalPhotonIDValueMap.variables:
    setattr(process.ntupler.localRecoMisc, key, cms.InputTag("HGCalPhotonIDValueMap", key))

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
)

process.p = cms.Path(
    process.HGCalPhotonIDValueMap  +
    process.ntupler
)

