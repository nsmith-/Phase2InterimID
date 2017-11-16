import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.StandardSequences.Eras import eras

process = cms.Process("Ntupler", eras.Phase2)
options = VarParsing('analysis')
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '93X_upgrade2023_realistic_v2', '')

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
)

process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.source = cms.Source ("PoolSource",
    fileNames  = cms.untracked.vstring("/store/mc/PhaseIITDRFall17DR/DiPhotonJetsBox_MGG-80toInf_14TeV-Sherpa/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/150000/24BB7BB2-A9B4-E711-9DC9-FA163E7FFB3C.root"),
)

process.load("RecoEgamma.Phase2InterimID.phase2EgammaRECO_cff")
process.p = cms.Path( process.phase2Egamma )

process.out = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    fileName = cms.untracked.string('file:reco.root'),
    outputCommands = cms.untracked.vstring(
        "keep *_hgcElectron*_*_*",
        "keep *_hgcPhoton*_*_*",
        "keep *_gedGsfElectrons_*_*",
        "keep *_cleanedEcalDrivenGsfElectronsFromMultiCl_*_*",
        "keep *_gedPhotons_*_*",
        "keep *_photonsFromMultiCl_*_*",
    ),
)
process.outstep = cms.EndPath(process.out)


