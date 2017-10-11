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
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "TTree.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
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
    float rho;
    int nPV;

    std::vector<float> gen_pt;
    std::vector<float> gen_eta;
    std::vector<float> gen_phi;
    std::vector<int> gen_id;
    std::vector<int> gen_parentId;
    std::vector<float> gen_fBrem;
    std::vector<int> gen_nGeantTracks;
    std::vector<bool> gen_isPromptFinalState;
    std::vector<int> gen_iLocalReco;
    std::vector<float> gen_localRecoDeltaR;
    std::vector<int> gen_iGedReco;
    std::vector<float> gen_gedRecoDeltaR;

    std::vector<float> localReco_pt;
    std::vector<float> localReco_eta;
    std::vector<float> localReco_phi;
    std::vector<EZBranch<reco::GsfElectron>> localReco_misc;

    std::vector<float> gedReco_pt;
    std::vector<float> gedReco_eta;
    std::vector<float> gedReco_phi;
    std::vector<EZBranch<reco::GsfElectron>> gedReco_misc;
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

    edm::EDGetTokenT<double> rhoToken_;
    edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
    edm::EDGetTokenT<reco::GsfElectronCollection> ecalDrivenElectronsToken_;
    edm::EDGetTokenT<reco::GsfElectronCollection> gedGsfElectronsToken_;
    edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
    edm::EDGetTokenT<SimClusterCollection> simClustersToken_;
    edm::EDGetTokenT<CaloParticleCollection> caloParticlesToken_;
    edm::EDGetTokenT<TrackingParticleCollection> trackingParticlesToken_;

    TTree * electronTree_;
    EventStruct event_;
};

Phase2ElectronTupler::Phase2ElectronTupler(const edm::ParameterSet& iConfig):
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexSrc"))),
  ecalDrivenElectronsToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("ecalDrivenElectrons"))),
  gedGsfElectronsToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("gedGsfElectrons"))),
  genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  simClustersToken_(consumes<SimClusterCollection>(iConfig.getParameter<edm::InputTag>("simClusters"))),
  caloParticlesToken_(consumes<CaloParticleCollection>(iConfig.getParameter<edm::InputTag>("caloParticles"))),
  trackingParticlesToken_(consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("trackingParticles")))
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  electronTree_ = fs->make<TTree>("electrons", "");
  electronTree_->Branch("run", &event_.run);
  electronTree_->Branch("lumi", &event_.lumi);
  electronTree_->Branch("event", &event_.event);
  electronTree_->Branch("rho", &event_.rho);
  electronTree_->Branch("nPV", &event_.nPV);

  electronTree_->Branch("gen_pt", &event_.gen_pt);
  electronTree_->Branch("gen_eta", &event_.gen_eta);
  electronTree_->Branch("gen_phi", &event_.gen_phi);
  electronTree_->Branch("gen_id", &event_.gen_id);
  electronTree_->Branch("gen_parentId", &event_.gen_parentId);
  electronTree_->Branch("gen_fBrem", &event_.gen_fBrem);
  electronTree_->Branch("gen_nGeantTracks", &event_.gen_nGeantTracks);
  electronTree_->Branch("gen_isPromptFinalState", &event_.gen_isPromptFinalState);
  electronTree_->Branch("gen_iLocalReco", &event_.gen_iLocalReco);
  electronTree_->Branch("gen_localRecoDeltaR", &event_.gen_localRecoDeltaR);
  electronTree_->Branch("gen_iGedReco", &event_.gen_iGedReco);
  electronTree_->Branch("gen_gedRecoDeltaR", &event_.gen_gedRecoDeltaR);

  electronTree_->Branch("localReco_pt", &event_.localReco_pt);
  electronTree_->Branch("localReco_eta", &event_.localReco_eta);
  electronTree_->Branch("localReco_phi", &event_.localReco_phi);
  auto localRecoMisc = iConfig.getParameter<edm::ParameterSet>("localRecoMisc");
  for (auto name : localRecoMisc.getParameterNames()) {
    if ( localRecoMisc.existsAs<std::string>(name) ) {
      event_.localReco_misc.emplace_back(electronTree_, "localReco_"+name, localRecoMisc.getParameter<std::string>(name));
    }
  }

  electronTree_->Branch("gedReco_pt", &event_.gedReco_pt);
  electronTree_->Branch("gedReco_eta", &event_.gedReco_eta);
  electronTree_->Branch("gedReco_phi", &event_.gedReco_phi);
  auto gedRecoMisc = iConfig.getParameter<edm::ParameterSet>("gedRecoMisc");
  for (auto name : gedRecoMisc.getParameterNames()) {
    if ( gedRecoMisc.existsAs<std::string>(name) ) {
      event_.gedReco_misc.emplace_back(electronTree_, "gedReco_"+name, gedRecoMisc.getParameter<std::string>(name));
    }
  }

}


Phase2ElectronTupler::~Phase2ElectronTupler()
{
}


void
Phase2ElectronTupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle<double> rhoH;
  iEvent.getByToken(rhoToken_, rhoH);

  Handle<reco::VertexCollection> vertexH;
  iEvent.getByToken(vtxToken_, vertexH);

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

  event_.run = iEvent.run();
  event_.lumi = iEvent.luminosityBlock();
  event_.event = iEvent.id().event();
  event_.rho = *rhoH;
  event_.nPV = vertexH->size();

  event_.localReco_pt.clear();
  event_.localReco_eta.clear();
  event_.localReco_phi.clear();
  for(auto&& b : event_.localReco_misc) b.clear();
  for(size_t iPho=0; iPho<ecalDrivenElectronsH->size(); ++iPho) {
    const auto& pho = ecalDrivenElectronsH->at(iPho);

    event_.localReco_pt.push_back(pho.pt());
    event_.localReco_eta.push_back(pho.eta());
    event_.localReco_phi.push_back(pho.phi());
    for(auto&& b : event_.localReco_misc) b.push_back(pho);
  }

  event_.gedReco_pt.clear();
  event_.gedReco_eta.clear();
  event_.gedReco_phi.clear();
  for(auto&& b : event_.gedReco_misc) b.clear();
  for(size_t iPho=0; iPho<gedGsfElectronsH->size(); ++iPho) {
    const auto& pho = gedGsfElectronsH->at(iPho);

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
  event_.gen_nGeantTracks.clear();
  event_.gen_iLocalReco.clear();
  event_.gen_localRecoDeltaR.clear();
  event_.gen_iGedReco.clear();
  event_.gen_gedRecoDeltaR.clear();
  for (size_t iGp=0; iGp<genParticlesH->size(); ++iGp) {
    const auto& gp = genParticlesH->at(iGp);
    if ( gp.pt() < 5. || gp.status() != 1 ) continue;
    if ( std::abs(gp.pdgId()) != 11 ) continue;
    event_.gen_pt.push_back(gp.pt());
    event_.gen_eta.push_back(gp.eta());
    event_.gen_phi.push_back(gp.phi());
    event_.gen_id.push_back(gp.pdgId());
    event_.gen_parentId.push_back(parentId(gp));
    event_.gen_isPromptFinalState.push_back(gp.isPromptFinalState());

    // Find the corresponding TrackingParticle object
    // which contains the GEANT history of the electrons
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
    int gen_nGeantTracks = -1.;
    if ( matchedTp != nullptr ) {
      gen_fBrem = (1-(matchedTp->g4Tracks().rbegin()->momentum().pt() - matchedTp->g4Tracks().begin()->momentum().pt())) / gp.pt();
      gen_nGeantTracks = matchedTp->g4Tracks().size();
    }
    event_.gen_fBrem.push_back(gen_fBrem);
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

  electronTree_->Fill();
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
