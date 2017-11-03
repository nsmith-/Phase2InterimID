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

ValueMap Setup
--------------
All electron and photon ID variables are either part of the `reco::` object or available by running
the following ValueMap producers.
```python
process.load("RecoEgamma.EgammaTools.hgcalElectronIDValueMap_cfi") import hgcalElectronIDValueMap
process.load("RecoEgamma.EgammaTools.hgcalPhotonIDValueMap_cfi") import hgcalPhotonIDValueMap

# example
process.path = cms.Path(process.hgcalElectronIDValueMap + process.hgcalPhotonIDValueMap + process.ntupler)
```

BDT Reader
----------
TODO
