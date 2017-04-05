Phase 2 Interim EG ID
---------------------
This package pieces together the necessary information from `AODSIM` to make 
sensible but /very preliminary/ cuts on forward (eta>1.4) electrons (and photons in the future)

You will need CMSSW 90X, 82X, or higher.
Place this inside RecoEgamma, e.g. with example recipe:
```bash
cmsrel CMSSW_8_2_0_patch1
cd CMSSW_8_2_0_patch1/src
cmsenv
git cms-init
git cms-addpkg RecoEgamma/ElectronIdentification
pushd RecoEgamma
git clone git@github.com:nsmith-/Phase2InterimID.git
popd
scram b
```

TODO: example setup

sketch:
```c++
#include "RecoEgamma/Phase2InterimID/interface/HGCalIDTool.h"

class Asdf : EDAnalyzer {
  // ...

  std::unique_ptr<HGCalIDTool> hgcEmId_;
};


Asdf::Asdf(const edm::ParameterSet& iConfig) {
  const edm::ParameterSet& hgcIdCfg = iConfig.getParameterSet("HGCalIDToolConfig");
  auto cc = consumesCollector();
  hgcEmId_.reset( new HGCalIDTool(hgcIdCfg, cc) );
}

void
Asdf::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  hgcEmId_->getEventSetup(iSetup);
  hgcEmId_->getEvent(iEvent);

  // ...
  for (const auto& el : ecalDrivenGsfElectrons) {
    if ( el->superCluster().isNonnull() && el->superCluster()->seed().isNonnull() ) {
      const reco::CaloCluster * cluster = el->superCluster()->seed().get();
      bool isHGCal = hgcEmId_->setClusterPtr(cluster);
      if ( isHGCal ) {
        hOverE_hgcalSafe_.push_back( hgcEmId_->getHadronFraction() );
        hgcId_startPosition_.push_back( std::abs(hgcEmId_->getStartPosition().z()) );
        hgcId_sigmaietaieta_.push_back( hgcEmId_->getSigmaEtaEta() );
        hgcId_lengthCompatibility_.push_back( hgcEmId_->getLengthCompatibility() );
      }
      else {
        // Use run 2 variables
      }
    }
  }
}
```
