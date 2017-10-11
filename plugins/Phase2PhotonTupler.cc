// -*- C++ -*-
//
// Package:    Analysis/Phase2PhotonTupler
// Class:      Phase2PhotonTupler
// 
/**\class Phase2PhotonTupler Phase2PhotonTupler.cc Analysis/Phase2PhotonTupler/plugins/Phase2PhotonTupler.cc

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
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "TTree.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

namespace {
  template<typename T>
  class EZBranch {
    public:
      EZBranch(TTree* tree, std::string name, std::string method) :
        ptr_(std::make_unique<std::vector<float>>()),
        branch_(tree->Branch(name.c_str(), ptr_.get())),
        fcn_(method)
      {};

      void push_back(const T& obj) {
        ptr_->push_back(fcn_(obj));
      };

      void clear() { ptr_->clear(); };

    private:
      std::unique_ptr<std::vector<float>> ptr_;
      TBranch * branch_;
      // lazy = false since we are going to be using the derived class
      StringObjectFunction<T, false> fcn_;
  };

  struct EventStruct {
    Long64_t run;
    Long64_t lumi;
    Long64_t event;
    std::vector<float> gen_pt;
    std::vector<float> gen_eta;
    std::vector<float> gen_phi;
    std::vector<int> gen_id;
    std::vector<int> gen_parentId;
    std::vector<float> gen_fBrem;
    std::vector<float> gen_conversionRho;
    std::vector<int> gen_nGeantTracks;
    std::vector<bool> gen_isPromptFinalState;
    std::vector<int> gen_iLocalReco;
    std::vector<float> gen_localRecoDeltaR;
    std::vector<int> gen_iGedReco;
    std::vector<float> gen_gedRecoDeltaR;

    std::vector<float> localReco_pt;
    std::vector<float> localReco_eta;
    std::vector<float> localReco_phi;
    std::vector<EZBranch<reco::Photon>> localReco_misc;

    std::vector<float> gedReco_pt;
    std::vector<float> gedReco_eta;
    std::vector<float> gedReco_phi;
    std::vector<EZBranch<reco::Photon>> gedReco_misc;
  };
}

class Phase2PhotonTupler : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit Phase2PhotonTupler(const edm::ParameterSet&);
    ~Phase2PhotonTupler();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    edm::EDGetTokenT<reco::PhotonCollection> photonsToken_;
    edm::EDGetTokenT<reco::PhotonCollection> gedPhotonsToken_;
    edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
    edm::EDGetTokenT<SimClusterCollection> simClustersToken_;
    edm::EDGetTokenT<CaloParticleCollection> caloParticlesToken_;
    edm::EDGetTokenT<TrackingParticleCollection> trackingParticlesToken_;
    edm::EDGetTokenT<TrackingVertexCollection> trackingVerticesToken_;

    TTree * photonTree_;
    EventStruct event_;

    StringCutObjectSelector<reco::GenParticle, false> genCut_;
};

Phase2PhotonTupler::Phase2PhotonTupler(const edm::ParameterSet& iConfig):
  photonsToken_(consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  gedPhotonsToken_(consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("gedPhotons"))),
  genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  simClustersToken_(consumes<SimClusterCollection>(iConfig.getParameter<edm::InputTag>("simClusters"))),
  caloParticlesToken_(consumes<CaloParticleCollection>(iConfig.getParameter<edm::InputTag>("caloParticles"))),
  trackingParticlesToken_(consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("trackingParticles"))),
  trackingVerticesToken_(consumes<TrackingVertexCollection>(iConfig.getParameter<edm::InputTag>("trackingVertices"))),
  genCut_(iConfig.getParameter<std::string>("genCut"))
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  photonTree_ = fs->make<TTree>("photons", "");
  photonTree_->Branch("run", &event_.run);
  photonTree_->Branch("lumi", &event_.lumi);
  photonTree_->Branch("event", &event_.event);
  photonTree_->Branch("gen_pt", &event_.gen_pt);
  photonTree_->Branch("gen_eta", &event_.gen_eta);
  photonTree_->Branch("gen_phi", &event_.gen_phi);
  photonTree_->Branch("gen_id", &event_.gen_id);
  photonTree_->Branch("gen_parentId", &event_.gen_parentId);
  photonTree_->Branch("gen_fBrem", &event_.gen_fBrem);
  photonTree_->Branch("gen_conversionRho", &event_.gen_conversionRho);
  photonTree_->Branch("gen_nGeantTracks", &event_.gen_nGeantTracks);
  photonTree_->Branch("gen_isPromptFinalState", &event_.gen_isPromptFinalState);
  photonTree_->Branch("gen_iLocalReco", &event_.gen_iLocalReco);
  photonTree_->Branch("gen_localRecoDeltaR", &event_.gen_localRecoDeltaR);
  photonTree_->Branch("gen_iGedReco", &event_.gen_iGedReco);
  photonTree_->Branch("gen_gedRecoDeltaR", &event_.gen_gedRecoDeltaR);

  photonTree_->Branch("localReco_pt", &event_.localReco_pt);
  photonTree_->Branch("localReco_eta", &event_.localReco_eta);
  photonTree_->Branch("localReco_phi", &event_.localReco_phi);
  auto localRecoMisc = iConfig.getParameter<edm::ParameterSet>("localRecoMisc");
  for (auto name : localRecoMisc.getParameterNames()) {
    if ( localRecoMisc.existsAs<std::string>(name) ) {
      event_.localReco_misc.emplace_back(photonTree_, "localReco_"+name, localRecoMisc.getParameter<std::string>(name));
    }
  }

  photonTree_->Branch("gedReco_pt", &event_.gedReco_pt);
  photonTree_->Branch("gedReco_eta", &event_.gedReco_eta);
  photonTree_->Branch("gedReco_phi", &event_.gedReco_phi);
  auto gedRecoMisc = iConfig.getParameter<edm::ParameterSet>("gedRecoMisc");
  for (auto name : gedRecoMisc.getParameterNames()) {
    if ( gedRecoMisc.existsAs<std::string>(name) ) {
      event_.gedReco_misc.emplace_back(photonTree_, "gedReco_"+name, gedRecoMisc.getParameter<std::string>(name));
    }
  }

}


Phase2PhotonTupler::~Phase2PhotonTupler()
{
}


void
Phase2PhotonTupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle<reco::PhotonCollection> photonsH;
  iEvent.getByToken(photonsToken_, photonsH);

  Handle<reco::PhotonCollection> gedPhotonsH;
  iEvent.getByToken(gedPhotonsToken_, gedPhotonsH);

  Handle<reco::GenParticleCollection> genParticlesH;
  iEvent.getByToken(genParticlesToken_, genParticlesH);

  Handle<SimClusterCollection> simClustersH;
  iEvent.getByToken(simClustersToken_, simClustersH);

  Handle<CaloParticleCollection> caloParticlesH;
  iEvent.getByToken(caloParticlesToken_, caloParticlesH);

  Handle<TrackingParticleCollection> trackingParticlesH;
  iEvent.getByToken(trackingParticlesToken_, trackingParticlesH);

  Handle<TrackingVertexCollection> trackingVerticesH;
  iEvent.getByToken(trackingVerticesToken_, trackingVerticesH);

  event_.run = iEvent.run();
  event_.lumi = iEvent.luminosityBlock();
  event_.event = iEvent.id().event();

  event_.localReco_pt.clear();
  event_.localReco_eta.clear();
  event_.localReco_phi.clear();
  for(auto&& b : event_.localReco_misc) b.clear();
  for(size_t iPho=0; iPho<photonsH->size(); ++iPho) {
    const auto& pho = photonsH->at(iPho);

    event_.localReco_pt.push_back(pho.pt());
    event_.localReco_eta.push_back(pho.eta());
    event_.localReco_phi.push_back(pho.phi());
    for(auto&& b : event_.localReco_misc) b.push_back(pho);
  }

  event_.gedReco_pt.clear();
  event_.gedReco_eta.clear();
  event_.gedReco_phi.clear();
  for(auto&& b : event_.gedReco_misc) b.clear();
  for(size_t iPho=0; iPho<gedPhotonsH->size(); ++iPho) {
    const auto& pho = gedPhotonsH->at(iPho);

    event_.gedReco_pt.push_back(pho.pt());
    event_.gedReco_eta.push_back(pho.eta());
    event_.gedReco_phi.push_back(pho.phi());
    for(auto&& b : event_.gedReco_misc) b.push_back(pho);
  }

  auto parentId = [](const reco::GenParticle& p) {
    if ( p.numberOfMothers() == 0 ) return p.pdgId();
    auto mom = p.mother(0);
    while ( mom->pdgId() == p.pdgId() and mom->numberOfMothers() > 0 ) mom = mom->mother(0);
    return mom->pdgId();
  };

  event_.gen_pt.clear();
  event_.gen_eta.clear();
  event_.gen_phi.clear();
  event_.gen_id.clear();
  event_.gen_parentId.clear();
  event_.gen_isPromptFinalState.clear();
  event_.gen_fBrem.clear();
  event_.gen_conversionRho.clear();
  event_.gen_nGeantTracks.clear();
  event_.gen_iLocalReco.clear();
  event_.gen_localRecoDeltaR.clear();
  event_.gen_iGedReco.clear();
  event_.gen_gedRecoDeltaR.clear();
  for (size_t iGp=0; iGp<genParticlesH->size(); ++iGp) {
    const auto& gp = genParticlesH->at(iGp);
    if ( not genCut_(gp) ) continue;
    event_.gen_pt.push_back(gp.pt());
    event_.gen_eta.push_back(gp.eta());
    event_.gen_phi.push_back(gp.phi());
    event_.gen_id.push_back(gp.pdgId());
    event_.gen_parentId.push_back(parentId(gp));
    event_.gen_isPromptFinalState.push_back(gp.isPromptFinalState());

    // Find the corresponding TrackingParticle object
    // which contains the GEANT history of the photons
    const TrackingParticle * matchedTp{nullptr};
    for (size_t iTp=0; iTp<trackingParticlesH->size(); ++iTp) {
      const auto& tp = trackingParticlesH->at(iTp);
      if ( tp.genParticles().empty() ) continue;
      if ( tp.genParticles().at(0) == reco::GenParticleRef(genParticlesH, iGp) ) {
        matchedTp = &tp;
        break;
      }
    }
    float gen_fBrem = -1.;
    float gen_conversionRho = -1.;
    int gen_nGeantTracks = -1.;
    if ( matchedTp != nullptr ) {
      gen_fBrem = (1-(matchedTp->g4Tracks().rbegin()->momentum().pt() - matchedTp->g4Tracks().begin()->momentum().pt())) / gp.pt();
      auto firstHit = matchedTp->decayVertices_begin();
      if ( firstHit != matchedTp->decayVertices_end() ) {
        gen_conversionRho = firstHit->get()->position().rho();
      }
      gen_nGeantTracks = matchedTp->g4Tracks().size();
    }
    event_.gen_fBrem.push_back(gen_fBrem);
    event_.gen_conversionRho.push_back(gen_conversionRho);
    event_.gen_nGeantTracks.push_back(gen_nGeantTracks);

    int iLocalReco = -1;
    float minDrLocalReco = 999.;
    for (size_t iReco=0; iReco<event_.localReco_pt.size(); ++iReco) {
      float drTemp = reco::deltaR(gp.eta(), gp.phi(), event_.localReco_eta[iReco], event_.localReco_phi[iReco]);
      if ( drTemp < minDrLocalReco ) {
        iLocalReco = iReco;
        minDrLocalReco = drTemp;
      }
    }
    event_.gen_iLocalReco.push_back(iLocalReco);
    event_.gen_localRecoDeltaR.push_back(minDrLocalReco);

    int iGedReco = -1;
    float minDrGedReco = 999.;
    for (size_t iReco=0; iReco<event_.gedReco_pt.size(); ++iReco) {
      float drTemp = reco::deltaR(gp.eta(), gp.phi(), event_.gedReco_eta[iReco], event_.gedReco_phi[iReco]);
      if ( drTemp < minDrGedReco ) {
        iGedReco = iReco;
        minDrGedReco = drTemp;
      }
    }
    event_.gen_iGedReco.push_back(iGedReco);
    event_.gen_gedRecoDeltaR.push_back(minDrGedReco);
  }

  photonTree_->Fill();
}


void 
Phase2PhotonTupler::beginJob()
{
}


void 
Phase2PhotonTupler::endJob() 
{
}


void
Phase2PhotonTupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(Phase2PhotonTupler);
