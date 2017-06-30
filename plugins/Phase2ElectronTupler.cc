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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

namespace {
  struct GenToRecoStruct {
    Long64_t run;
    Long64_t lumi;
    Long64_t event;
    float gen_pt;
    float gen_eta;
    float gen_phi;
    float gen_id;
    float gen_fBrem;
    float gen_nGeantTracks;
    float gen_isPromptFinalState;
    float localReco_deltaR;
    float localReco_pt;
    float localReco_eta;
    float localReco_phi;
    float gedReco_deltaR;
    float gedReco_pt;
    float gedReco_eta;
    float gedReco_phi;
    float gedReco_eSC;
    float gedIsoChargedHadrons;
    float gedIsoNeutralHadrons;
    float gedIsoPhotons;
    float gedIsoChargedFromPU;
  };
}

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
    edm::EDGetTokenT<reco::GsfElectronCollection> gedGsfElectronsToken_;
    edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
    edm::EDGetTokenT<SimClusterCollection> simClustersToken_;
    edm::EDGetTokenT<CaloParticleCollection> caloParticlesToken_;
    edm::EDGetTokenT<TrackingParticleCollection> trackingParticlesToken_;

    TTree * genToRecoTree_;
    GenToRecoStruct genToReco_;
};

Phase2ElectronTupler::Phase2ElectronTupler(const edm::ParameterSet& iConfig):
  ecalDrivenElectronsToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("ecalDrivenElectrons"))),
  gedGsfElectronsToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("gedGsfElectrons"))),
  genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  simClustersToken_(consumes<SimClusterCollection>(iConfig.getParameter<edm::InputTag>("simClusters"))),
  caloParticlesToken_(consumes<CaloParticleCollection>(iConfig.getParameter<edm::InputTag>("caloParticles"))),
  trackingParticlesToken_(consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("trackingParticles")))
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  genToRecoTree_ = fs->make<TTree>("genToReco", "All gen electrons, nearest reco");
  genToRecoTree_->Branch("run", &genToReco_.run);
  genToRecoTree_->Branch("lumi", &genToReco_.lumi);
  genToRecoTree_->Branch("event", &genToReco_.event);
  genToRecoTree_->Branch("gen_pt", &genToReco_.gen_pt);
  genToRecoTree_->Branch("gen_eta", &genToReco_.gen_eta);
  genToRecoTree_->Branch("gen_phi", &genToReco_.gen_phi);
  genToRecoTree_->Branch("gen_id", &genToReco_.gen_id);
  genToRecoTree_->Branch("gen_fBrem", &genToReco_.gen_fBrem);
  genToRecoTree_->Branch("gen_nGeantTracks", &genToReco_.gen_nGeantTracks);
  genToRecoTree_->Branch("gen_isPromptFinalState", &genToReco_.gen_isPromptFinalState);
  genToRecoTree_->Branch("localReco_deltaR", &genToReco_.localReco_deltaR);
  genToRecoTree_->Branch("localReco_pt", &genToReco_.localReco_pt);
  genToRecoTree_->Branch("localReco_eta", &genToReco_.localReco_eta);
  genToRecoTree_->Branch("localReco_phi", &genToReco_.localReco_phi);
  genToRecoTree_->Branch("gedReco_deltaR", &genToReco_.gedReco_deltaR);
  genToRecoTree_->Branch("gedReco_pt", &genToReco_.gedReco_pt);
  genToRecoTree_->Branch("gedReco_eta", &genToReco_.gedReco_eta);
  genToRecoTree_->Branch("gedReco_phi", &genToReco_.gedReco_phi);
  genToRecoTree_->Branch("gedIsoChargedHadrons", &genToReco_.gedIsoChargedHadrons);
  genToRecoTree_->Branch("gedIsoNeutralHadrons", &genToReco_.gedIsoNeutralHadrons);
  genToRecoTree_->Branch("gedIsoPhotons", &genToReco_.gedIsoPhotons);
  genToRecoTree_->Branch("gedIsoChargedFromPU", &genToReco_.gedIsoChargedFromPU);
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

  Handle<reco::GsfElectronCollection> gedGsfElectronsH;
  iEvent.getByToken(gedGsfElectronsToken_, gedGsfElectronsH);

  Handle<reco::GenParticleCollection> genParticlesH;
  iEvent.getByToken(genParticlesToken_, genParticlesH);

  Handle<SimClusterCollection> simClustersH;
  iEvent.getByToken(simClustersToken_, simClustersH);

  Handle<CaloParticleCollection> caloParticlesH;
  iEvent.getByToken(caloParticlesToken_, caloParticlesH);

  Handle<TrackingParticleCollection> trackingParticlesH;
  iEvent.getByToken(trackingParticlesToken_, trackingParticlesH);

  genToReco_.run = iEvent.run();
  genToReco_.lumi = iEvent.luminosityBlock();
  genToReco_.event = iEvent.id().event();

  for (size_t iGp=0; iGp<genParticlesH->size(); ++iGp) {
    const auto& gp = genParticlesH->at(iGp);
    if ( gp.pt() < 5. || gp.status() != 1 || std::abs(gp.pdgId()) != 11 ) continue;
    genToReco_.gen_pt = gp.pt();
    genToReco_.gen_eta = gp.eta();
    genToReco_.gen_phi = gp.phi();
    genToReco_.gen_id = gp.pdgId();
    genToReco_.gen_isPromptFinalState = gp.isPromptFinalState();

    // Find the corresponding TrackingParticle object
    // which contains the GEANT history of the electron
    const TrackingParticle * matchedTp{nullptr};
    for (size_t iTp=0; iTp<trackingParticlesH->size(); ++iTp) {
      const auto& tp = trackingParticlesH->at(iTp);
      if ( tp.genParticles().empty() ) continue;
      if ( tp.genParticles().at(0) == reco::GenParticleRef(genParticlesH, iGp) ) {
        matchedTp = &tp;
        break;
      }
    }
    genToReco_.gen_fBrem = -1.;
    genToReco_.gen_nGeantTracks = -1.;
    if ( matchedTp != nullptr ) {
      genToReco_.gen_fBrem = (matchedTp->g4Tracks().rbegin()->momentum().pt() - matchedTp->g4Tracks().begin()->momentum().pt()) / gp.pt();
      genToReco_.gen_nGeantTracks = matchedTp->g4Tracks().size();
    }


    // Helper function find closest electron to gen
    auto closestElectronPtr = [gp](const reco::GsfElectronCollection& coll) {
      auto ptr = coll.end();
      double minDR = std::numeric_limits<double>::infinity();
      for (auto el=coll.begin(); el!=coll.end(); el++) {
        if ( reco::deltaR(*el, gp) < minDR ) {
          minDR = reco::deltaR(*el, gp);
          ptr = el;
        }
      }
      return ptr;
    };

    
    // Closest ecalDriven electron (i.e. local GSF track reco seeded by superclusters)
    auto ecalDrivenPtr = closestElectronPtr(*ecalDrivenElectronsH);
    genToReco_.localReco_deltaR = 999.;
    genToReco_.localReco_pt = 0.;
    genToReco_.localReco_eta = 0.;
    genToReco_.localReco_phi = 0.;
    if ( ecalDrivenPtr != ecalDrivenElectronsH->end() ) {
      genToReco_.localReco_deltaR = reco::deltaR(*ecalDrivenPtr, gp);
      genToReco_.localReco_pt = ecalDrivenPtr->pt();
      genToReco_.localReco_eta = ecalDrivenPtr->eta();
      genToReco_.localReco_phi = ecalDrivenPtr->phi();
    }


    // Closest GED electron (i.e. particle flow electron)
    auto gedPtr = closestElectronPtr(*gedGsfElectronsH);
    genToReco_.gedReco_deltaR = 999.;
    genToReco_.gedReco_pt = 0.;
    genToReco_.gedReco_eta = 0.;
    genToReco_.gedReco_phi = 0.;
    genToReco_.gedIsoChargedHadrons = 0.;
    genToReco_.gedIsoNeutralHadrons = 0.;
    genToReco_.gedIsoPhotons = 0.;
    genToReco_.gedIsoChargedFromPU = 0.;
    if ( gedPtr != gedGsfElectronsH->end() ) {
      genToReco_.gedReco_deltaR = reco::deltaR(*gedPtr, gp);
      genToReco_.gedReco_pt = gedPtr->pt();
      genToReco_.gedReco_eta = gedPtr->eta();
      genToReco_.gedReco_phi = gedPtr->phi();

      reco::GsfElectron::PflowIsolationVariables pfIso = gedPtr->pfIsolationVariables();
      genToReco_.gedIsoChargedHadrons = pfIso.sumChargedHadronPt;
      genToReco_.gedIsoNeutralHadrons = pfIso.sumNeutralHadronEt;
      genToReco_.gedIsoPhotons = pfIso.sumPhotonEt;
      genToReco_.gedIsoChargedFromPU = pfIso.sumPUPt;
    }


    genToRecoTree_->Fill();
  }

  /*
  int nGeantElectrons{0};
  for (size_t iTp=0; iTp<trackingParticlesH->size(); ++iTp) {
    const auto& tp = trackingParticlesH->at(iTp);
    if ( std::abs(tp.pdgId()) != 11 || tp.status() != 1 || tp.pt() < 5. ) continue;
    // std::cout << tp << std::endl;

    double genIsolation{0.};
    for (size_t jTp=0; jTp<trackingParticlesH->size(); ++jTp) {
      if ( jTp == iTp ) continue;
      const auto& tp2 = trackingParticlesH->at(jTp);
      if ( reco::deltaR(tp, tp2) > 0.3 || tp2.vertex() != tp.vertex() ) continue;
      // std::cout << " -- iso cone" << tp2 << std::endl;
      genIsolation += tp2.pt();
    }
    genIsolation /= tp.pt();
    std::cout << "Gen. isolation = " << genIsolation << std::endl;

    nGeantElectrons++;
  }

  for (size_t iEl=0; iEl<ecalDrivenElectronsH->size(); iEl++) {
    break;
    const auto& el = ecalDrivenElectronsH->at(iEl);
    std::cout << "New ecalDrivenElectron candidate ------------------------------" << std::endl;
    std::cout << "Polar p4: " << el.pt() << ", " << el.eta() << ", " << el.phi() << std::endl;
    std::cout << "Cartesian p4: " << el.p4() << std::endl;
    
    auto matchedGenPtr = genParticlesH->end();
    double minDR = std::numeric_limits<double>::infinity();
    for (auto gp=genParticlesH->begin(); gp!=genParticlesH->end(); gp++) {
      if ( gp->status() != 1 ) continue;
      if ( reco::deltaR(*gp, el) < minDR ) {
        minDR = reco::deltaR(*gp, el);
        matchedGenPtr = gp;
      }
    }
    std::cout << "Nearest genParticle deltaR=" << reco::deltaR(*matchedGenPtr, el) << " info:" << std::endl;
    std::cout << "  genP id, p4, pt, eta: " << matchedGenPtr->pdgId() << ", " << matchedGenPtr->p4() << ", " << matchedGenPtr->pt() << ", " << matchedGenPtr->eta() << std::endl;
    
    std::vector<const SimCluster *> nearbyClusters;
    for (const SimCluster& clu : *simClustersH) {
      if ( reco::deltaR(clu, el) < 0.1 ) {
        nearbyClusters.push_back(&clu);
      }
    }
    std::sort(nearbyClusters.begin(), nearbyClusters.end(), [&el](const SimCluster* a, const SimCluster* b) -> double { return reco::deltaR(*a, el) < reco::deltaR(*b, el); });
    for (auto clu : nearbyClusters) {
      std::cout << "Nearby simCluster dR=" << reco::deltaR(*clu, el) << std::endl;
      std::cout << *clu << std::endl;
    }

    std::vector<const CaloParticle *> nearbyCaloParticles;
    for (const CaloParticle& cp : *caloParticlesH) {
      if ( reco::deltaR(cp, el) < 0.1 ) {
        nearbyCaloParticles.push_back(&cp);
      }
    }
    std::sort(nearbyCaloParticles.begin(), nearbyCaloParticles.end(), [&el](const CaloParticle* a, const CaloParticle* b) -> double { return reco::deltaR(*a, el) < reco::deltaR(*b, el); });
    for (auto cp : nearbyCaloParticles) {
      std::cout << "Nearby caloParticle dR=" << reco::deltaR(*cp, el) << std::endl;
      std::cout << *cp << std::endl;
    }
  }
  */

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
