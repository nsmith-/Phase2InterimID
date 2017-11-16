// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/PatCandidates/interface/Photon.h"

class Phase2PhotonMerger : public edm::stream::EDProducer<> {
  public:
    explicit Phase2PhotonMerger(const edm::ParameterSet&);
    ~Phase2PhotonMerger();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginStream(edm::StreamID) override;
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endStream() override;

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    edm::EDGetTokenT<pat::PhotonCollection> barrelPhotonToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> barrelIDToken_;
    StringCutObjectSelector<pat::Photon> barrelCut_;
    edm::EDGetTokenT<pat::PhotonCollection> endcapPhotonToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> endcapIDToken_;
    StringCutObjectSelector<pat::Photon> endcapCut_;
};

Phase2PhotonMerger::Phase2PhotonMerger(const edm::ParameterSet& iConfig):
  barrelPhotonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("barrelPhotons"))),
  barrelIDToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("barrelID"))),
  barrelCut_(iConfig.getParameter<std::string>("barrelCut")),
  endcapPhotonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("endcapPhotons"))),
  endcapIDToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("endcapID"))),
  endcapCut_(iConfig.getParameter<std::string>("endcapCut"))
{
  produces<pat::PhotonCollection>();
}


Phase2PhotonMerger::~Phase2PhotonMerger()
{
}

void
Phase2PhotonMerger::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle<pat::PhotonCollection> barrelPhotonHandle;
  iEvent.getByToken(barrelPhotonToken_, barrelPhotonHandle);
  const auto barrelInputId = barrelPhotonHandle.id();

  Handle<ValueMap<float>> barrelIDHandle;
  iEvent.getByToken(barrelIDToken_, barrelIDHandle);

  Handle<pat::PhotonCollection> endcapPhotonHandle;
  iEvent.getByToken(endcapPhotonToken_, endcapPhotonHandle);
  const auto endcapInputId = endcapPhotonHandle.id();

  Handle<ValueMap<float>> endcapIDHandle;
  iEvent.getByToken(endcapIDToken_, endcapIDHandle);

  std::unique_ptr<pat::PhotonCollection> photons(new pat::PhotonCollection());

  for(size_t iPho=0; iPho<barrelPhotonHandle->size(); ++iPho) {
    pat::Photon phoNew(barrelPhotonHandle->at(iPho));
    phoNew.addUserFloat("mvaValue", barrelIDHandle->get(barrelInputId , iPho));
    if ( barrelCut_(phoNew) ) photons->push_back(phoNew);
  }
  for(size_t iPho=0; iPho<endcapPhotonHandle->size(); ++iPho) {
    pat::Photon phoNew(endcapPhotonHandle->at(iPho));
    // Common access method
    phoNew.addUserFloat("mvaValue", endcapIDHandle->get(endcapInputId , iPho));
    if ( endcapCut_(phoNew) ) photons->push_back(phoNew);
  }

  iEvent.put(std::move(photons));
}


void
Phase2PhotonMerger::beginStream(edm::StreamID)
{
}


void
Phase2PhotonMerger::endStream() {
}


/*
void
Phase2PhotonMerger::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 

/*
void
Phase2PhotonMerger::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 

/*
void
Phase2PhotonMerger::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 

/*
void
Phase2PhotonMerger::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 

void
Phase2PhotonMerger::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(Phase2PhotonMerger);
