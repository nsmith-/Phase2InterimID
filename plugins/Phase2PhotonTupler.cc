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
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "RecoEgamma/EgammaMCTools/interface/ElectronMCTruth.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruth.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruthFinder.h"

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

  template<typename T>
  class ValueMapBranch {
    public:
      ValueMapBranch(TTree* tree, std::string name, edm::EDGetTokenT<edm::ValueMap<float>>&& token) :
        ptr_(std::make_unique<std::vector<float>>()),
        branch_(tree->Branch(name.c_str(), ptr_.get())),
        tok_(token)
      {};

      template<typename CollectionType>
      void clearAndSet(const edm::Event& e, const edm::Handle<CollectionType>& coll) {
        e.getByToken(tok_, map_);
        collId_ = coll.id();
        ptr_->clear();
      };

      void push_back(size_t i) {
        ptr_->push_back(map_->get(collId_, i));
      };

    private:
      std::unique_ptr<std::vector<float>> ptr_;
      TBranch * branch_;
      edm::EDGetTokenT<edm::ValueMap<float>> tok_;
      edm::Handle<edm::ValueMap<float>> map_;
      edm::ProductID collId_;
  };

  struct EventStruct {
    Long64_t run;
    Long64_t lumi;
    Long64_t event;
    float rho;
    std::vector<float> gen_pt;
    std::vector<float> gen_eta;
    std::vector<float> gen_phi;
    std::vector<int> gen_id;
    std::vector<int> gen_parentId;
    std::vector<float> gen_fBrem;
    std::vector<float> gen_conversionRho;
    std::vector<float> gen_conversionZ;
    std::vector<int> gen_nGeantTracks;
    std::vector<bool> gen_isPromptFinalState;
    std::vector<int> gen_iLocalReco;
    std::vector<float> gen_localRecoDeltaR;
    std::vector<int> gen_iGedReco;
    std::vector<float> gen_gedRecoDeltaR;

    std::vector<float> localReco_pt;
    std::vector<float> localReco_eta;
    std::vector<float> localReco_phi;
    std::vector<int> localReco_iGen;
    std::vector<EZBranch<reco::Photon>> localReco_misc;
    std::vector<ValueMapBranch<reco::Photon>> localReco_valuemaps;
    std::vector<float> localReco_matchedGsfChi2;
    std::vector<float> localReco_matchedGsfLostHits;

    std::vector<float> gedReco_pt;
    std::vector<float> gedReco_eta;
    std::vector<float> gedReco_phi;
    std::vector<int> gedReco_iGen;
    std::vector<EZBranch<reco::Photon>> gedReco_misc;
    std::vector<ValueMapBranch<reco::Photon>> gedReco_valuemaps;
    std::vector<float> gedReco_conversionSafeElectronVeto;
    std::vector<float> gedReco_TPmetric;
    std::vector<int> gedReco_TPid;
    std::vector<float> gedReco_TPpt;
    std::vector<float> gedReco_TPeta;
    std::vector<float> gedReco_TPphi;
    std::vector<int> gedReco_TPevent;
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
    edm::EDGetTokenT<std::vector<SimTrack>> simTracksToken_;
    edm::EDGetTokenT<std::vector<SimVertex>> simVerticesToken_;
    edm::EDGetTokenT<double> rhoToken_;
    edm::EDGetTokenT<reco::GsfElectronCollection> ecalDrivenElectronsToken_;
    edm::EDGetTokenT<reco::GsfElectronCollection> gedGsfElectronsToken_;
    edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
    edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;

    bool doPremixContent_;
    edm::EDGetTokenT<TrackingParticleCollection> trackingParticlesToken_;
    edm::EDGetTokenT<TrackingVertexCollection> trackingVerticesToken_;
    edm::EDGetTokenT<edm::PCaloHitContainer> caloHitsToken_;

    TTree * photonTree_;
    EventStruct event_;

    StringCutObjectSelector<reco::Photon, false> localRecoCut_;
    StringCutObjectSelector<reco::Photon, false> gedRecoCut_;
    StringCutObjectSelector<reco::GenParticle, false> genCut_;
};

Phase2PhotonTupler::Phase2PhotonTupler(const edm::ParameterSet& iConfig):
  photonsToken_(consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  gedPhotonsToken_(consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("gedPhotons"))),
  genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  simClustersToken_(consumes<SimClusterCollection>(iConfig.getParameter<edm::InputTag>("simClusters"))),
  caloParticlesToken_(consumes<CaloParticleCollection>(iConfig.getParameter<edm::InputTag>("caloParticles"))),
  simTracksToken_(consumes<std::vector<SimTrack>>(iConfig.getParameter<edm::InputTag>("simTracksSrc"))),
  simVerticesToken_(consumes<std::vector<SimVertex>>(iConfig.getParameter<edm::InputTag>("simVerticesSrc"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoSrc"))),
  ecalDrivenElectronsToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("ecalDrivenElectrons"))),
  gedGsfElectronsToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("gedGsfElectrons"))),
  conversionsToken_(consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conversions"))),
  beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"))),
  doPremixContent_(iConfig.getParameter<bool>("doPremixContent")),
  trackingParticlesToken_(consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("trackingParticles"))),
  trackingVerticesToken_(consumes<TrackingVertexCollection>(iConfig.getParameter<edm::InputTag>("trackingVertices"))),
  caloHitsToken_(consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("caloHits"))),
  localRecoCut_(iConfig.getParameter<std::string>("localRecoCut")),
  gedRecoCut_(iConfig.getParameter<std::string>("gedRecoCut")),
  genCut_(iConfig.getParameter<std::string>("genCut"))
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  photonTree_ = fs->make<TTree>("photons", "");
  photonTree_->Branch("run", &event_.run);
  photonTree_->Branch("lumi", &event_.lumi);
  photonTree_->Branch("event", &event_.event);
  photonTree_->Branch("rho",  &event_.rho, "rho/F");
  photonTree_->Branch("gen_pt", &event_.gen_pt);
  photonTree_->Branch("gen_eta", &event_.gen_eta);
  photonTree_->Branch("gen_phi", &event_.gen_phi);
  photonTree_->Branch("gen_id", &event_.gen_id);
  photonTree_->Branch("gen_parentId", &event_.gen_parentId);
  photonTree_->Branch("gen_fBrem", &event_.gen_fBrem);
  photonTree_->Branch("gen_conversionRho", &event_.gen_conversionRho);
  photonTree_->Branch("gen_conversionZ", &event_.gen_conversionZ);
  photonTree_->Branch("gen_nGeantTracks", &event_.gen_nGeantTracks);
  photonTree_->Branch("gen_isPromptFinalState", &event_.gen_isPromptFinalState);
  photonTree_->Branch("gen_iLocalReco", &event_.gen_iLocalReco);
  photonTree_->Branch("gen_localRecoDeltaR", &event_.gen_localRecoDeltaR);
  photonTree_->Branch("gen_iGedReco", &event_.gen_iGedReco);
  photonTree_->Branch("gen_gedRecoDeltaR", &event_.gen_gedRecoDeltaR);

  photonTree_->Branch("localReco_pt", &event_.localReco_pt);
  photonTree_->Branch("localReco_eta", &event_.localReco_eta);
  photonTree_->Branch("localReco_phi", &event_.localReco_phi);
  photonTree_->Branch("localReco_iGen", &event_.localReco_iGen);
  auto localRecoMisc = iConfig.getParameter<edm::ParameterSet>("localRecoMisc");
  for (auto name : localRecoMisc.getParameterNames()) {
    if ( localRecoMisc.existsAs<std::string>(name) ) {
      event_.localReco_misc.emplace_back(photonTree_, "localReco_"+name, localRecoMisc.getParameter<std::string>(name));
    }
    else if ( localRecoMisc.existsAs<edm::InputTag>(name) ) {
      event_.localReco_valuemaps.emplace_back(photonTree_, "localReco_"+name, consumes<edm::ValueMap<float>>(localRecoMisc.getParameter<edm::InputTag>(name)));
    }
  }
  photonTree_->Branch("localReco_matchedGsfChi2", &event_.localReco_matchedGsfChi2);
  photonTree_->Branch("localReco_matchedGsfLostHits", &event_.localReco_matchedGsfLostHits);

  photonTree_->Branch("gedReco_pt", &event_.gedReco_pt);
  photonTree_->Branch("gedReco_eta", &event_.gedReco_eta);
  photonTree_->Branch("gedReco_phi", &event_.gedReco_phi);
  photonTree_->Branch("gedReco_iGen", &event_.gedReco_iGen);
  auto gedRecoMisc = iConfig.getParameter<edm::ParameterSet>("gedRecoMisc");
  for (auto name : gedRecoMisc.getParameterNames()) {
    if ( gedRecoMisc.existsAs<std::string>(name) ) {
      event_.gedReco_misc.emplace_back(photonTree_, "gedReco_"+name, gedRecoMisc.getParameter<std::string>(name));
    }
    else if ( gedRecoMisc.existsAs<edm::InputTag>(name) ) {
      event_.gedReco_valuemaps.emplace_back(photonTree_, "gedReco_"+name, consumes<edm::ValueMap<float>>(gedRecoMisc.getParameter<edm::InputTag>(name)));
    }
  }
  photonTree_->Branch("gedReco_conversionSafeElectronVeto", &event_.gedReco_conversionSafeElectronVeto);
  if ( doPremixContent_ ) {
    photonTree_->Branch("gedReco_TPmetric", &event_.gedReco_TPmetric);
    photonTree_->Branch("gedReco_TPid", &event_.gedReco_TPid);
    photonTree_->Branch("gedReco_TPpt", &event_.gedReco_TPpt);
    photonTree_->Branch("gedReco_TPeta", &event_.gedReco_TPeta);
    photonTree_->Branch("gedReco_TPphi", &event_.gedReco_TPphi);
    photonTree_->Branch("gedReco_TPevent", &event_.gedReco_TPevent);
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

  Handle<double> rhoH;
  iEvent.getByToken(rhoToken_, rhoH);

  Handle<reco::GsfElectronCollection> ecalDrivenElectronsH;
  iEvent.getByToken(ecalDrivenElectronsToken_, ecalDrivenElectronsH);

  Handle<reco::GsfElectronCollection> gedGsfElectronsH;
  iEvent.getByToken(gedGsfElectronsToken_, gedGsfElectronsH);

  Handle<reco::ConversionCollection> conversionsH;
  iEvent.getByToken(conversionsToken_, conversionsH);

  Handle<reco::BeamSpot> beamspotH;
  iEvent.getByToken(beamspotToken_, beamspotH);

  Handle<TrackingParticleCollection> trackingParticlesH;
  Handle<TrackingVertexCollection> trackingVerticesH;
  Handle<PCaloHitContainer> caloHitsH;
  if ( doPremixContent_ ) {
    iEvent.getByToken(trackingParticlesToken_, trackingParticlesH);
    iEvent.getByToken(trackingVerticesToken_, trackingVerticesH);
    iEvent.getByToken(caloHitsToken_, caloHitsH);
  }

  Handle<std::vector<SimTrack>> simTracks;
  iEvent.getByToken(simTracksToken_ ,simTracks);
  Handle<std::vector<SimVertex>> simVertices;
  iEvent.getByToken(simVerticesToken_ ,simVertices);
  std::unique_ptr<PhotonMCTruthFinder> thePhotonMCTruthFinder_(new PhotonMCTruthFinder());
  std::vector<PhotonMCTruth> mcPhotons;
  mcPhotons = thePhotonMCTruthFinder_->find(*simTracks,  *simVertices);

  event_.run = iEvent.run();
  event_.lumi = iEvent.luminosityBlock();
  event_.event = iEvent.id().event();

  event_.rho = *rhoH;

  event_.localReco_pt.clear();
  event_.localReco_eta.clear();
  event_.localReco_phi.clear();
  for(auto&& b : event_.localReco_misc) b.clear();
  for(auto&& b : event_.localReco_valuemaps) b.clearAndSet(iEvent, photonsH);
  event_.localReco_matchedGsfChi2.clear();
  event_.localReco_matchedGsfLostHits.clear();
  for(size_t iPho=0; iPho<photonsH->size(); ++iPho) {
    const auto& pho = photonsH->at(iPho);
    if ( not localRecoCut_(pho) ) continue;

    event_.localReco_pt.push_back(pho.pt());
    event_.localReco_eta.push_back(pho.eta());
    event_.localReco_phi.push_back(pho.phi());
    for(auto&& b : event_.localReco_misc) b.push_back(pho);
    for(auto&& b : event_.localReco_valuemaps) b.push_back(iPho);
   
    int nel = 0;
    for(const auto& el : *ecalDrivenElectronsH) {
      if ( reco::deltaR(*el.superCluster(), *pho.superCluster()) > 0.01 ) continue;
      nel++;
      if ( nel == 1 ) {
        event_.localReco_matchedGsfChi2.push_back( el.gsfTrack()->chi2() );
        event_.localReco_matchedGsfLostHits.push_back( el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) );
      } else {
        std::cout << "We can find more than one electron from the same SC I guess?" << std::endl;
      }
    }
    if ( nel == 0 ) {
      event_.localReco_matchedGsfChi2.push_back( 1.e5 );
      event_.localReco_matchedGsfLostHits.push_back( -1. );
    }
  }

  event_.gedReco_pt.clear();
  event_.gedReco_eta.clear();
  event_.gedReco_phi.clear();
  for(auto&& b : event_.gedReco_misc) b.clear();
  for(auto&& b : event_.gedReco_valuemaps) b.clearAndSet(iEvent, photonsH);
  event_.gedReco_conversionSafeElectronVeto.clear();
  if ( doPremixContent_ ) {
    event_.gedReco_TPmetric.clear();
    event_.gedReco_TPid.clear();
    event_.gedReco_TPpt.clear();
    event_.gedReco_TPeta.clear();
    event_.gedReco_TPphi.clear();
    event_.gedReco_TPevent.clear();
  }
  for(size_t iPho=0; iPho<gedPhotonsH->size(); ++iPho) {
    const auto& pho = gedPhotonsH->at(iPho);
    if ( not gedRecoCut_(pho) ) continue;

    event_.gedReco_pt.push_back(pho.pt());
    event_.gedReco_eta.push_back(pho.eta());
    event_.gedReco_phi.push_back(pho.phi());
    for(auto&& b : event_.gedReco_misc) b.push_back(pho);
    for(auto&& b : event_.gedReco_valuemaps) b.push_back(iPho);

    event_.gedReco_conversionSafeElectronVeto.push_back( !ConversionTools::hasMatchedPromptElectron(pho.superCluster(), gedGsfElectronsH, conversionsH, beamspotH->position()) );

    if ( doPremixContent_ ) {
      std::cout << "New photon pt=" << pho.pt() << ", all calohits:" << std::endl;
      for(auto hit_frac : pho.superCluster()->hitsAndFractions()) {
        std::map<int, float> track_e;
        for(auto ch : *caloHitsH) {
          if ( ch.id() == hit_frac.first ) {
            track_e[ch.geantTrackId()] += ch.energy();
          }
        }
        for(auto te : track_e) {
          for (const auto& tp : *trackingParticlesH) {
            if ( tp.g4Tracks().size() > 0 and tp.g4Tracks()[0].trackId() == (unsigned) te.first ) {
              std::cout << "G4 track in TP, sum e =" << te.second << std::endl;
              std::cout << tp <<std::endl;
            }
          }
        }
      }
      std::cout << std::endl;

      const TrackingParticle * matchedTp{nullptr};
      float metric = std::numeric_limits<float>::infinity();
      for (const auto& tp : *trackingParticlesH) {
        // if ( tp.parentVertex()->nSourceTracks() > 0 ) continue;
        // float mtmp = std::hypot(reco::deltaR(tp, pho), std::log10(tp.energy()/pho.energy()));
        float mtmp = std::hypot((tp.vertex()-pho.caloPosition()).R(), 10*std::log10(tp.energy()/pho.energy()));
        if ( mtmp < metric ) {
          matchedTp = &tp;
          metric = mtmp;
        }
      }
      event_.gedReco_TPmetric.push_back(metric);
      event_.gedReco_TPid.push_back(matchedTp->pdgId());
      event_.gedReco_TPpt.push_back(matchedTp->pt());
      event_.gedReco_TPeta.push_back(matchedTp->eta());
      event_.gedReco_TPphi.push_back(matchedTp->phi());
      event_.gedReco_TPevent.push_back(matchedTp->eventId().bunchCrossing()*1000+matchedTp->eventId().event());
    }
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
  event_.gen_conversionZ.clear();
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

    if ( doPremixContent_ ) {
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
      // event_.gen_conversionRho.push_back(gen_conversionRho);
      (void) gen_conversionRho;
      event_.gen_nGeantTracks.push_back(gen_nGeantTracks);
    }

    float conversionRho = -1.;
    float conversionZ = 0.;
    for(auto& pmc : mcPhotons) {
      // auto simTrack = std::find_if(simTracks->begin(), simTracks->end(), [pmc](const SimTrack& t) { return t.trackId() == (size_t) pmc.trackId(); });
      // size_t iGenPart = simTrack->genpartIndex();
      // (but this is genParticles not prunedGenParticles)
      // This is easier than the more correct alternative, first one is always initial G4 track
      if ( reco::deltaR(pmc.fourMomentum(), gp.p4()) < 0.001 ) {
        // PhotonMCTruth::vertex() is the conversion vertex (not mother vertex)
        // So the first one will always have isAConversion() true, but where it
        // interacted tells us if it matters, since all photons will at least interact
        // by the time they hit ECAL
        conversionRho = pmc.vertex().perp();
        conversionZ = pmc.vertex().z();
        break;
      }
    }
    event_.gen_conversionRho.push_back(conversionRho);
    event_.gen_conversionZ.push_back(conversionZ);

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
  
  // Now do reco -> gen match
  // which theoretically is not the same as gen -> reco
  event_.localReco_iGen.clear();
  for (size_t iReco=0; iReco<event_.localReco_pt.size(); ++iReco) {
    int local_iGen = -1;
    float minDrLocalGen = 999.;
    for (size_t iGen=0; iGen<event_.gen_pt.size(); ++iGen) {
      float drTemp = reco::deltaR(event_.gen_eta[iGen], event_.gen_phi[iGen], event_.localReco_eta[iReco], event_.localReco_phi[iReco]);
      if ( drTemp < minDrLocalGen ) {
        local_iGen = iGen;
        minDrLocalGen = drTemp;
      }
    }
    // Define 'no match'
    if ( minDrLocalGen > 0.1 ) local_iGen = -1;
    event_.localReco_iGen.push_back(local_iGen);
  }
  event_.gedReco_iGen.clear();
  for (size_t iReco=0; iReco<event_.gedReco_pt.size(); ++iReco) {
    int ged_iGen = -1;
    float minDrGedGen = 999.;
    for (size_t iGen=0; iGen<event_.gen_pt.size(); ++iGen) {
      float drTemp = reco::deltaR(event_.gen_eta[iGen], event_.gen_phi[iGen], event_.gedReco_eta[iReco], event_.gedReco_phi[iReco]);
      if ( drTemp < minDrGedGen ) {
        ged_iGen = iGen;
        minDrGedGen = drTemp;
      }
    }
    // Define 'no match'
    if ( minDrGedGen > 0.1 ) ged_iGen = -1;
    event_.gedReco_iGen.push_back(ged_iGen);
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
