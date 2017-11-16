Phase 2 Interim EGamma ID
=========================
This package has been used to develop photon IDs for the phase 2 HGCal TDR, and now also serves as a recipe and helper code repository.

Recipe
------
```bash
cmsrel CMSSW_9_3_2
cd CMSSW_9_3_2/src
cmsenv
git cms-init
git cms-merge-topic -u nsmith-:EgammaFromMultiCl_932v2
mkdir -p RecoEgamma && pushd RecoEgamma
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

If you want to use RECO objects, you will have to load the two separate collections into your analysis.  If you prefer to use PAT objects, you can use a combined collection, see below.
All electron and photon ID MVA input variables are either part of the `reco::` object or available by running the following ValueMap producers.
```python
process.load("RecoEgamma.Phase2InterimID.phase2EgammaRECO_cff")

# e.g. 
process.ntupler = cms.EDAnalyzer("MyTuples",
    barrelPhotons = cms.InputTag("gedPhotons"),
    barrelPhoMva  = cms.InputTag("hgcPhotonMVAbarrel"),
    endcapPhotons = cms.InputTag("photonsFromMultiCl"),
    endcapPhoMva  = cms.InputTag("hgcPhotonMVAendcap"),
)
process.p = cms.Path( process.phase2Egamma + process.ntupler )
```
See `test/testPhase2EgammaCollectionsRECO.py` for a more complete example.
It is suggested to save the MVA value and cut later.  See below for working points.

Running on PAT objects
----------------------
No endcap EGamma collections exist in MiniAOD yet.  See [here](https://github.com/cms-sw/cmssw/pull/21037) for status.
However, one can always run with `secondaryInputFiles` to access the RECO collections.  In CRAB, there is a simple `useParent` option.

To aid in this, a `cff` has been made that forms the endcap PAT collections, runs ID on both barrel and endcap photons, and produces a merged collection with embedded MVA values.
It can be imported as follows:
```python
process.load("RecoEgamma.Phase2InterimID.phase2EgammaPAT_cff")

# e.g.
process.ntupler.patPhotonsSrc = cms.InputTag("phase2Photons")
process.ntupler.patElectronsSrc = cms.InputTag("phase2Electrons")
process.p = cms.Path( process.phase2Egamma + process.ntupler )
```
See `test/testPhase2EgammaCollections.py` for a more complete example.
The ID value is accessible as `userFloat("mvaValue")`.  Use `isEB()` to decide whether the photon is from HGCal multiclusters or standard barrel GED.
It is suggested to save the MVA value and cut later.  See below for working points.

Cut working points
------------------
# Electrons
See https://github.com/CMS-HGCAL/EgammaTools/blob/master/ELECTRONBDT.md#recommended-id-cuts

# Photons
Current MVA is `V4`, pass `>=` value.

 | MVA name | Loose WP | Tight WP |
 | -------- | -------- | -------- |
 | barrelV4 |   0.00   |   0.56   |
 | endcapV4 |   0.20   |   0.68   |

