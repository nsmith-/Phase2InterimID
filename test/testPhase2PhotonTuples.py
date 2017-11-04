import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.StandardSequences.Eras import eras

process = cms.Process("Ntupler", eras.Phase2)
options = VarParsing('analysis')
options.register(
    "phase2",
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "phase 2 sample"
)
options.parseArguments()

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    numberOfThreads = cms.untracked.uint32(4 if options.phase2 else 1),
    numberOfStreams = cms.untracked.uint32(0),
)

process.load("FWCore.MessageService.MessageLogger_cfi")
if options.phase2:
    process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring(options.inputFiles) )

phoSrc = cms.InputTag("gedPhotons")
if options.phase2:
    phoSrc = cms.InputTag("photonsFromMultiCl")
    process.load("RecoEgamma.EgammaTools.hgcalPhotonIDValueMap_cfi")
    process.load("RecoEgamma.Phase2InterimID.hgcalPhotonMVAProducer_cfi")
    process.hgcalPhotonIDValueMap.photons = phoSrc

process.ntupler = cms.EDAnalyzer("Phase2PhotonTupler",
    photons = phoSrc,
    localRecoCut = cms.string("pt>10 && !isEB"),
    gedRecoCut = cms.string("pt>10 && isEB"),
    gedPhotons = cms.InputTag("gedPhotons"),
    genParticles = cms.InputTag("genParticles"),
    genCut = cms.string("pt>5 && status==1 && (abs(pdgId)==11 || pdgId==22)"),
    rhoSrc = cms.InputTag("fixedGridRhoFastjetAll"),
    ecalDrivenElectrons = cms.InputTag("ecalDrivenGsfElectronsFromMultiCl" if options.phase2 else "gedGsfElectrons"),
    gedGsfElectrons = cms.InputTag("gedGsfElectrons"),
    conversions = cms.InputTag("conversions"),
    beamspot = cms.InputTag("offlineBeamSpot"),
    vertices = cms.InputTag("offlinePrimaryVertices"),
    doRecoContent = cms.bool(options.phase2),
    ecalRecHits = cms.InputTag("ecalRecHit:EcalRecHitsEB"),
    simClusters = cms.InputTag("mix:MergedCaloTruth"),
    caloParticles = cms.InputTag("mix:MergedCaloTruth"),
    simTracksSrc = cms.InputTag("g4SimHits"),
    simVerticesSrc = cms.InputTag("g4SimHits"),
    doPremixContent = cms.bool(False),
    trackingParticles = cms.InputTag("mix:MergedTrackTruth"),
    trackingVertices = cms.InputTag("mix:MergedTrackTruth"),
    caloHits = cms.InputTag("g4SimHits:EcalHitsEB"),
)

common = cms.PSet(
    scEta = cms.string("superCluster().eta()"),
    scEnergy = cms.string("superCluster().energy()"),
    scRawEnergy = cms.string("superCluster().rawEnergy()"),
    seedOrigEnergy = cms.string("superCluster().seed().energy()"),
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

if options.phase2:
    process.ntupler.localRecoMisc = cms.PSet(
        common,
        seed_det = cms.string("superCluster().seed().hitsAndFractions().at(0).first.det()"),
        seed_subdet = cms.string("superCluster().seed().hitsAndFractions().at(0).first.subdetId()"),
    )
    for key in process.hgcalPhotonIDValueMap.variables:
        setattr(process.ntupler.localRecoMisc, key, cms.InputTag("hgcalPhotonIDValueMap", key))
    process.ntupler.localRecoMisc.mvaValue = cms.InputTag("hgcalPhotonMVA")
else:
    process.ntupler.localRecoMisc = cms.PSet(process.ntupler.gedRecoMisc)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
)

if options.phase2:
    process.p = cms.Path(
        process.hgcalPhotonIDValueMap  +
        process.hgcalPhotonMVA +
        process.ntupler
    )
else:
    process.p = cms.Path( process.ntupler )

