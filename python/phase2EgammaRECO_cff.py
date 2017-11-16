import FWCore.ParameterSet.Config as cms

from RecoEgamma.EgammaTools.cleanedEcalDrivenGsfElectronsFromMultiCl_cfi import cleanedEcalDrivenGsfElectronsFromMultiCl
from RecoEgamma.EgammaTools.hgcalElectronIDValueMap_cff import hgcalElectronIDValueMap
from RecoEgamma.EgammaTools.hgcalPhotonIDValueMap_cff import hgcalPhotonIDValueMap
from RecoEgamma.Phase2InterimID.hgcalElectronMVAProducer_cfi import hgcalElectronMVA
from RecoEgamma.Phase2InterimID.hgcalPhotonMVAProducer_cfi import hgcalPhotonMVA

# Make sure all of these are in path or task
hgcElectronID = hgcalElectronIDValueMap.clone(
    electrons = cms.InputTag("cleanedEcalDrivenGsfElectronsFromMultiCl"),
)
hgcPhotonID = hgcalPhotonIDValueMap.clone()
hgcElectronMVAbarrel = hgcalElectronMVA.clone(electrons=cms.InputTag("gedGsfElectrons"))
hgcElectronMVAendcap = hgcalElectronMVA.clone(electrons=cms.InputTag("cleanedEcalDrivenGsfElectronsFromMultiCl"))
hgcPhotonMVAbarrel = hgcalPhotonMVA.clone(photons=cms.InputTag("gedPhotons"))
hgcPhotonMVAendcap = hgcalPhotonMVA.clone(photons=cms.InputTag("photonsFromMultiCl"))

phase2EgammaTask = cms.Task(
    cleanedEcalDrivenGsfElectronsFromMultiCl,
    hgcElectronID,
    hgcPhotonID,
    hgcElectronMVAbarrel,
    hgcElectronMVAendcap,
    hgcPhotonMVAbarrel,
    hgcPhotonMVAendcap,
)

phase2Egamma = cms.Sequence(phase2EgammaTask)
