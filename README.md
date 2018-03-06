Phase 2 Interim EGamma ID
=========================
This package has been used to develop photon IDs for the phase 2 HGCal TDR, and now also serves as a recipe and helper code repository for both electrons and photons.

Recipe
------
```bash
cmsrel CMSSW_9_3_5
cd CMSSW_9_3_5/src
cmsenv
git cms-init
mkdir -p RecoEgamma && pushd RecoEgamma
git clone git@github.com:nsmith-/Phase2InterimID.git
popd
scram b -j 8
```

Running on MiniAOD tier
-----------------------
MINIAOD is the simplest and preferred method.
The ID input variables are already computed for you in the `PhaseIISpr18AODMiniAOD` campaign, but the final MVA ID value is not precomputed.
To compute the MVA value, a `cff` has been made which runs ID on barrel and endcap collections, and 
produces merged collections of electrons and photons with embedded MVA values.
It can be imported as follows:
```python
process.load("RecoEgamma.Phase2InterimID.phase2EgammaPAT_cff")

# e.g.
process.ntupler.patPhotonsSrc = cms.InputTag("phase2Photons")
process.ntupler.patElectronsSrc = cms.InputTag("phase2Electrons")
process.p = cms.Path( process.phase2Egamma + process.ntupler )
```
See `test/testPhase2EgammaCollections.py` for a more complete example.
The ID value is accessible as `userFloat("mvaValue")`.  Use `isEB()` to decide whether the object is from HGCal multiclusters or standard barrel GED.
It is suggested to save the MVA value and cut later.  See below for working points.

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
All electron and photon ID MVA input variables are either part of the `reco::` object or available by ValueMap producers.
A `cff` has been made that runs all the relevant producers, and can be imported as follows:
```python
process.load("RecoEgamma.Phase2InterimID.phase2EgammaRECO_cff")

# e.g. 
process.ntupler = cms.EDAnalyzer("MyTuples",
    barrelElectrons = cms.InputTag("gedGsfElectrons"),
    barrelElectronMVA  = cms.InputTag("hgcElectronMVAbarrel"),
    endcapElectrons = cms.InputTag("cleanedEcalDrivenGsfElectronsFromMultiCl"),
    endcapElectronMVA  = cms.InputTag("hgcElectronMVAendcap"),
    barrelPhotons = cms.InputTag("gedPhotons"),
    barrelPhoMVA  = cms.InputTag("hgcPhotonMVAbarrel"),
    endcapPhotons = cms.InputTag("photonsFromMultiCl"),
    endcapPhoMVA  = cms.InputTag("hgcPhotonMVAendcap"),
)
process.p = cms.Path( process.phase2Egamma + process.ntupler )
```
See `test/testPhase2EgammaCollectionsRECO.py` for a more complete example.
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

Recommended Kinematics
----------------------
For barrel photons, the best performing energy is found using an algorithm developed at the time of the barrel TDR, which adds the energy of the top 15 most energetic crystals in the supercluster.  This improves the pileup resistance of the whole sum, while still maintaining a reasonable shower containment.
An example implementation can be found [here](https://github.com/nsmith-/Phase2InterimID/blob/integrated/plugins/Phase2PhotonTupler.cc#L451-L490).
The supercluster position `superCluster()->position()` seems to work in a satisfactory way for barrel photons.

For endcap photons, the best performing energy variable (so far) is the seed multicluster energy, `photon.superCluster()->seed()->energy()`.
In the endcap, the seed cluster position `superCluster()->seed->position()` is probably the best choice for now to compute the momentum vector, after choosing the appropriate vertex.
