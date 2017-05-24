// -*- C++ -*-
//
// Package:    Analysis/Phase2ElectronTupler
// Class:      Phase2ElectronTupler
// 
/**\class Phase2ElectronTupler Phase2ElectronTupler.cc Analysis/Phase2ElectronTupler/plugins/Phase2ElectronTupler.cc

 Description: 

 Implementation:
     
*/
//
// Original Author:  Nicholas Charles Smith
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"

class Phase2ElectronTupler : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit Phase2ElectronTupler(const edm::ParameterSet&);
    ~Phase2ElectronTupler();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    edm::EDGetTokenT<reco::GsfElectronCollection> ecalDrivenElectronsToken_;
    edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
    edm::EDGetTokenT<SimClusterCollection> simClustersToken_;
};

Phase2ElectronTupler::Phase2ElectronTupler(const edm::ParameterSet& iConfig):
  ecalDrivenElectronsToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("ecalDrivenElectrons"))),
  genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  simClustersToken_(consumes<SimClusterCollection>(iConfig.getParameter<edm::InputTag>("simClusters")))
{
  usesResource("TFileService");
}


Phase2ElectronTupler::~Phase2ElectronTupler()
{
}


void
Phase2ElectronTupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle<reco::GsfElectronCollection> ecalDrivenElectronsH;
  iEvent.getByToken(ecalDrivenElectronsToken_, ecalDrivenElectronsH);

  Handle<reco::GenParticleCollection> genParticlesH;
  iEvent.getByToken(genParticlesToken_, genParticlesH);

  Handle<SimClusterCollection> simClustersH;
  iEvent.getByToken(simClustersToken_, simClustersH);
}


void 
Phase2ElectronTupler::beginJob()
{
}


void 
Phase2ElectronTupler::endJob() 
{
}


void
Phase2ElectronTupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(Phase2ElectronTupler);
