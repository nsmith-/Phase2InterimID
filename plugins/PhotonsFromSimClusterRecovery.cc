// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

class PhotonsFromSimClusterRecovery : public edm::stream::EDProducer<> {
  public:
    explicit PhotonsFromSimClusterRecovery(const edm::ParameterSet&);
    ~PhotonsFromSimClusterRecovery();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginStream(edm::StreamID) override;
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endStream() override;

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    edm::EDGetTokenT<reco::PhotonCollection> barrelPhotonToken_;
    StringCutObjectSelector<reco::Photon> barrelCut_;
    edm::EDGetTokenT<reco::PFCandidateCollection> pfCandsToken_;
    edm::EDGetTokenT<reco::SuperClusterCollection> hgcalSuperClustersToken_;
    StringCutObjectSelector<reco::SuperCluster> scCut_;
    edm::EDGetTokenT<reco::VertexCollection> primaryVerticesToken_;
};

PhotonsFromSimClusterRecovery::PhotonsFromSimClusterRecovery(const edm::ParameterSet& iConfig):
  barrelPhotonToken_(consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("barrelPhotons"))),
  barrelCut_(iConfig.getParameter<std::string>("barrelCut")),
  pfCandsToken_(consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandidates"))),
  hgcalSuperClustersToken_(consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("pfSuperClustersHGCal"))),
  scCut_(iConfig.getParameter<std::string>("endcapSuperClusterCut")),
  primaryVerticesToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices")))
{
  produces<reco::PhotonCollection>();
  produces<reco::PhotonCoreCollection>("cores");
}


PhotonsFromSimClusterRecovery::~PhotonsFromSimClusterRecovery()
{
}

void
PhotonsFromSimClusterRecovery::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle<reco::PhotonCollection> barrelPhotonHandle;
  iEvent.getByToken(barrelPhotonToken_, barrelPhotonHandle);

  Handle<reco::PFCandidateCollection> pfCandsHandle;
  iEvent.getByToken(pfCandsToken_, pfCandsHandle);

  Handle<reco::SuperClusterCollection> scHandle;
  iEvent.getByToken(hgcalSuperClustersToken_, scHandle);

  Handle<reco::VertexCollection> verticesH;
  iEvent.getByToken(primaryVerticesToken_, verticesH);


  std::unique_ptr<reco::PhotonCollection> photons(new reco::PhotonCollection());
  std::unique_ptr<reco::PhotonCoreCollection> photonCores(new reco::PhotonCoreCollection());
  auto coreRefProd = iEvent.getRefBeforePut<reco::PhotonCoreCollection>("cores");

  for(const auto& pho : *barrelPhotonHandle) {
    if ( barrelCut_(pho) ) {
      photons->emplace_back(pho);
    }
  }

  std::vector<const reco::PFCandidate*> simPhotonsHGCal;
  for(const auto& pf : *pfCandsHandle) {
    // simPFProducer doesn't bother to fill in ecalEnergy
    // so easy way to tell if from HGCal or barrel
    if ( pf.pdgId() == 22 and pf.ecalEnergy() == 0. ) {
      simPhotonsHGCal.push_back(&pf);
    }
  }

  std::vector<Ref<reco::SuperClusterCollection>> matchedSuperClusters(simPhotonsHGCal.size());
  for(size_t iSc=0; iSc<scHandle->size(); ++iSc) {
    const auto& sc = scHandle->at(iSc);
    for(size_t iPf=0; iPf<simPhotonsHGCal.size(); ++iPf) {
      // PF constructed from sim cluster collection, but superClusterRef() is invalid :(
      // In parallel, superclusters built from same collection, hopefully a seed matches
      if ( reco::deltaR(*sc.seed().get(), *simPhotonsHGCal[iPf]) < 0.001 ) {
        matchedSuperClusters[iPf] = Ref<reco::SuperClusterCollection>(scHandle, iSc);
      }
    }
  }

  const math::XYZPoint vertex = (verticesH->size()>0) ? verticesH->at(0).position() : math::XYZPoint(0., 0., 0.);

  for(size_t iPho=0; iPho<simPhotonsHGCal.size(); ++iPho) {
    if ( matchedSuperClusters[iPho].isNull() ) {
      continue;
    }
    if ( not scCut_(*matchedSuperClusters[iPho].get()) ) continue;

    photonCores->emplace_back();
    auto& core = photonCores->back();
    core.setSuperCluster(matchedSuperClusters[iPho]);
    core.setPFlowPhoton(true);
    core.setStandardPhoton(true);

    // Use seed or full supercluster for default kinematics?
    const auto* cluster = matchedSuperClusters[iPho].get();
    // const auto* cluster = matchedSuperClusters[iPho]->seed().get();

    const auto p3 = (cluster->position() - vertex).unit() * cluster->energy();
    reco::Candidate::LorentzVector p4(p3.x(), p3.y(), p3.z(), cluster->energy());
    photons->emplace_back(p4, cluster->position(), reco::PhotonCoreRef(coreRefProd, photonCores->size()-1), vertex);
  }

  iEvent.put(std::move(photons));
  iEvent.put(std::move(photonCores), "cores");
}


void
PhotonsFromSimClusterRecovery::beginStream(edm::StreamID)
{
}


void
PhotonsFromSimClusterRecovery::endStream() {
}


/*
void
PhotonsFromSimClusterRecovery::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 

/*
void
PhotonsFromSimClusterRecovery::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 

/*
void
PhotonsFromSimClusterRecovery::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 

/*
void
PhotonsFromSimClusterRecovery::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 

void
PhotonsFromSimClusterRecovery::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(PhotonsFromSimClusterRecovery);
