import FWCore.ParameterSet.Config as cms

from RecoEgamma.EgammaTools.slimmedEgammaFromMultiCl_cff import *

from RecoEgamma.Phase2InterimID.hgcalElectronMVAProducer_cfi import hgcalElectronMVA
hgcElectronMVAbarrel = hgcalElectronMVA.clone(electrons=cms.InputTag("slimmedElectrons"), usePAT=True)
hgcElectronMVAendcap = hgcalElectronMVA.clone(electrons=cms.InputTag("ecalDrivenGsfElectronsFromMultiCl"), usePAT=False)
from RecoEgamma.Phase2InterimID.hgcalPhotonMVAProducer_cfi import hgcalPhotonMVA
hgcPhotonMVAbarrel = hgcalPhotonMVA.clone(photons=cms.InputTag("slimmedPhotons"), usePAT=True)
hgcPhotonMVAendcap = hgcalPhotonMVA.clone()
# insert at PAT creation time for simplicity
patPhotonsFromMultiCl.userData.userFloats.src.append(cms.InputTag("hgcPhotonMVAendcap"))

phase2Photons = cms.EDProducer("Phase2PhotonMerger",
    barrelPhotons = cms.InputTag("slimmedPhotons"),
    barrelID = cms.InputTag("hgcPhotonMVAbarrel"),
    barrelCut = cms.string(""),
    endcapPhotons = cms.InputTag("slimmedPhotonsFromMultiCl"),
    endcapCut = cms.string(""),
)

phase2EgammaTask = cms.Task(
    slimmedEgammaFromMultiClTask,
    hgcElectronMVAbarrel,
    hgcElectronMVAendcap,
    hgcPhotonMVAbarrel,
    hgcPhotonMVAendcap,
    phase2Photons
)

# Caution: won't work in scheduled mode!
phase2Egamma = cms.Sequence(phase2EgammaTask)
