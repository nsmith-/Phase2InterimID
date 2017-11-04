Phase 2 Interim EGamma ID
=========================
This package has been used to develop photon IDs for the phase 2 HGCal TDR

Recipe
------
```bash
cmsrel CMSSW_9_3_2
cd CMSSW_9_3_2/src
cmsenv
git cms-init
git cms-merge-topic -u nsmith-:EgammaFromMultiCl_932
mkdir RecoEgamma && pushd RecoEgamma
git clone -b integrated git@github.com:nsmith-/Phase2InterimID.git
popd
scram b -j 8
```

Running on RECO tier
--------------------
In `RECO`, the electrons and photons are split into two collections for barrel and endcap,

 | Product type                  | Product name                        | Description                                                                                |
 |-------------------------------|-------------------------------------|--------------------------------------------------------------------------------------------|
 | `reco::GsfElectronCollection` | `gedGsfElectrons`                   | Barrel electrons from the particle-flow global event description                           |
 | `reco::GsfElectronCollection` | `ecalDrivenGsfElectronsFromMultiCl` | Endcap electrons using local GSF electron reconstruction seeded by the HGCal multiclusters |
 | `reco::PhotonCollection`      | `gedPhotons`                        | Barrel photons from the particle-flow global event description                             |
 | `reco::PhotonCollection`      | `photonsFromMultiCl`                | Endcap photons using local 'island cluster' reconstruction, seeded by HGCal multiclusters  |

If you want to use RECO objects, you will have to load the two separate collections into your analysis.  If you prefer to use PAT objects, see below.
All electron and photon ID MVA input variables are either part of the `reco::` object or available by running the following ValueMap producers.
```python
from RecoEgamma.EgammaTools.hgcalElectronIDValueMap_cfi import hgcalElectronIDValueMap
from RecoEgamma.EgammaTools.hgcalPhotonIDValueMap_cfi import hgcalPhotonIDValueMap
# TODO: electron
from RecoEgamma.Phase2InterimID.hgcalPhotonMVAProducer_cfi import hgcalPhotonMVA

# Make sure all of these are in path or task
process.hgcElectronID = hgcalElectronIDValueMap.clone()
process.hgcPhotonID = hgcalPhotonIDValueMap.clone()
# TODO: electron
process.hgcPhotonMVAbarrel = hgcalPhotonMVA.clone(photons=cms.InputTag("gedPhotons"))
process.hgcPhotonMVAendcap = hgcalPhotonMVA.clone()

# e.g. 
process.ntupler = cms.EDAnalyzer("MyTuples",
    barrelPhotons = cms.InputTag("gedPhotons"),
    barrelPhoMva  = cms.InputTag("hgcPhotonMVAbarrel"),
    endcapPhotons = cms.InputTag("photonsFromMultiCl"),
    endcapPhoMva  = cms.InputTag("hgcPhotonMVAendcap"),
    ...
)
```

Running on PAT objects
----------------------
No endcap EGamma collections exist in MiniAOD yet.  See [here](https://github.com/cms-sw/cmssw/pull/21037) for status.
However, one can always run with `secondaryInputFiles` to access the RECO collections.  In CRAB, there is a simple `useParent` option.

To aid in this, a `cff` has been made that forms the endcap PAT collections, runs ID on both barrel and endcap photons, and produces a merged collection with embedded MVA values.
It can be imported as follows:
```python
process.load("RecoEgamma.Phase2InterimID.phase2EgammaPAT_cff")
# The phase2Egamma sequence won't work in scheduled mode
process.options.allowUnscheduled = cms.untracked.bool(True)

# e.g.
process.ntupler.patPhotonsSrc = cms.InputTag("phase2Photons")
process.p = cms.Path( process.phase2Egamma + process.ntupler )
```
See `test/testPhase2EgammaCollections.py` for a more complete example.
The ID value is accessible as `userFloat("mvaValue")`

