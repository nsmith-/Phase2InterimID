Phase 2 Interim EG ID
=====================
This package pieces together the necessary information from `AODSIM` to make 
sensible but *very preliminary* cuts on forward (`eta>1.4`) electrons (and photons in the future)

Recipe
------
You will need CMSSW 90X, 82X, or higher.
Place this inside RecoEgamma, e.g. with example recipe:
```bash
cmsrel CMSSW_8_2_0_patch1
cd CMSSW_8_2_0_patch1/src
cmsenv
git cms-init
git cms-merge-topic -u nsmith-/phase2_hgcalInterimID
pushd RecoEgamma
git clone git@github.com:nsmith-/Phase2InterimID.git
popd
scram b
```

ValueMap Setup
--------------
TODO

Examples
--------
Sketch example to use the tool:
```c++
#include "RecoEgamma/Phase2InterimID/interface/HGCalIDTool.h"

class AsdfSuperTupler : public edm::EDAnalyzer {
  // ...

  std::unique_ptr<HGCalIDTool> hgcEmId_;
};


AsdfSuperTupler::AsdfSuperTupler(const edm::ParameterSet& iConfig) {
  const edm::ParameterSet& hgcIdCfg = iConfig.getParameterSet("HGCalIDToolConfig");
  auto cc = consumesCollector();
  hgcEmId_.reset( new HGCalIDTool(hgcIdCfg, cc) );
}

void
AsdfSuperTupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  hgcEmId_->getEventSetup(iSetup);
  hgcEmId_->getEvent(iEvent);

  // ...
  
  for (const auto& el : ecalDrivenGsfElectrons) {
    bool isHGCal = hgcEmId_->setElectronPtr(&el);
    if ( isHGCal ) {
      // These are some std::vector<float> variables
      // you can add to your nTuple or whatever...
      hgcId_hadronFraction.push_back( hgcEmId_->getClusterHadronFraction() );
      hgcId_absEndcapShowerZ.push_back( std::abs(hgcEmId_->getClusterStartPosition().z()) );
      hgcId_sigmaEtaEta.push_back( hgcEmId_->getClusterSigmaEtaEta() );
      hgcId_lengthCompatibility.push_back( hgcEmId_->getClusterLengthCompatibility() );
    }
    else {
      // Use run 2 variables
    }
  }
}
```

Example cfg file:
```python
import FWCore.ParameterSet.Config as cms
process = cms.Process("ANA")
# ...

process.load("RecoParticleFlow.PFClusterProducer.particleFlowRecHitHGC_cff")

# Bonus: jurassic track isolation
# Since PF seems to be a bit buggy lately...
process.load("RecoEgamma.EgammaIsolationAlgos.electronTrackIsolationLcone_cfi")
process.electronTrackIsolationLcone.electronProducer = cms.InputTag("ecalDrivenGsfElectrons")
process.electronTrackIsolationLcone.intRadiusBarrel = 0.04
process.electronTrackIsolationLcone.intRadiusEndcap = 0.04

process.ntupler = cms.EDAnalyzer("AsdfSuperTupler",
    HGCalIDToolConfig = cms.PSet(
        HGCBHInput = cms.InputTag("HGCalRecHit","HGCHEBRecHits"),
        HGCEEInput = cms.InputTag("HGCalRecHit","HGCEERecHits"),
        HGCFHInput = cms.InputTag("HGCalRecHit","HGCHEFRecHits"),
        HGCPFRecHits = cms.InputTag("particleFlowRecHitHGC::ANA"),
        withPileup = cms.bool(True),
        debug = cms.bool(False),
    ),
    trackIsoValueMap = cms.InputTag("electronTrackIsolationLcone"),
)

process.p = cms.Path(process.electronTrackIsolationLcone+process.particleFlowRecHitHGCSeq+process.ntupler)
```

Provenance
----------
The main code was ripped from [HGCALShowerBasedEmIdentification](https://github.com/cms-sw/cmssw/blob/CMSSW_6_2_X_SLHC/RecoEcal/EgammaClusterAlgos/src/HGCALShowerBasedEmIdentification.cc) which only exists in `CMSSW_6_2_X_SLHC` at the moment.  Direct porting was not possible to some missing variables in the `PFCandidate` data format.  The `PFRecHitFraction`s had to be rebuilt because the `AODSIM` output dropped all the `PFRecHit`s (there is a "Cleaned" collection but it is empty).
