import FWCore.ParameterSet.Config as cms

from RecoEgamma.Phase2InterimID.hgcalElectronMVAProducer_cfi import hgcalElectronMVA
hgcElectronMVAbarrel = hgcalElectronMVA.clone(electrons=cms.InputTag("slimmedElectrons"), usePAT=True)
hgcElectronMVAendcap = hgcalElectronMVA.clone(electrons=cms.InputTag("slimmedElectronsFromMultiCl"), usePAT=True)
from RecoEgamma.Phase2InterimID.hgcalPhotonMVAProducer_cfi import hgcalPhotonMVA
hgcPhotonMVAbarrel = hgcalPhotonMVA.clone(photons=cms.InputTag("slimmedPhotons"), usePAT=True)
hgcPhotonMVAendcap = hgcalPhotonMVA.clone(photons=cms.InputTag("slimmedPhotonsFromMultiCl"), electronsFromMultiCl=cms.InputTag("slimmedElectronsFromMultiCl"), usePAT=True)

phase2Electrons = cms.EDProducer("Phase2ElectronMerger",
    barrelElectrons = cms.InputTag("slimmedElectrons"),
    barrelID = cms.InputTag("hgcElectronMVAbarrel"),
    barrelCut = cms.string(""),
    endcapElectrons = cms.InputTag("slimmedElectronsFromMultiCl"),
    endcapID = cms.InputTag("hgcElectronMVAendcap"),
    endcapCut = cms.string(""),
)

phase2Photons = cms.EDProducer("Phase2PhotonMerger",
    barrelPhotons = cms.InputTag("slimmedPhotons"),
    barrelID = cms.InputTag("hgcPhotonMVAbarrel"),
    barrelCut = cms.string(""),
    endcapPhotons = cms.InputTag("slimmedPhotonsFromMultiCl"),
    endcapID = cms.InputTag("hgcPhotonMVAendcap"),
    endcapCut = cms.string(""),
)

phase2EgammaTask = cms.Task(
    hgcElectronMVAbarrel,
    hgcElectronMVAendcap,
    hgcPhotonMVAbarrel,
    hgcPhotonMVAendcap,
    phase2Electrons,
    phase2Photons
)

phase2Egamma = cms.Sequence(phase2EgammaTask)
