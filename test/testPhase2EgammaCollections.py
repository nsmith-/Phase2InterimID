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
process.GlobalTag = GlobalTag(process.GlobalTag, '93X_upgrade2023_realistic_v5', '')

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
)

process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.source = cms.Source ("PoolSource",
    fileNames = cms.untracked.vstring([
        '/store/mc/PhaseIISpr18AODMiniAOD/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v5-v1/100000/744E0163-5B1B-E811-B29F-A0369F7FC540.root',
        '/store/mc/PhaseIISpr18AODMiniAOD/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v5-v1/100000/DA1E877B-631B-E811-89E8-A0369F7FC6EC.root',
        # '/store/mc/PhaseIISpr18AODMiniAOD/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/MINIAODSIM/noPU_93X_upgrade2023_realistic_v5-v1/70000/80052AEC-1D1B-E811-A90D-0025905C53D0.root',
    ]),
)

process.load("RecoEgamma.Phase2InterimID.phase2EgammaPAT_cff")
process.p = cms.Path( process.phase2Egamma )

process.out = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    fileName = cms.untracked.string('file:miniaod.root'),
    outputCommands = cms.untracked.vstring(
        "keep *_phase2Photons_*_*",
        "keep *_phase2Electrons_*_*",
    ),
)
process.outstep = cms.EndPath(process.out)


