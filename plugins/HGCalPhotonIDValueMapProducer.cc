// -*- C++ -*-
//
// Package:   RecoEgamma/HGCalPhotonIDValueMapProducer
// Class:    HGCalPhotonIDValueMapProducer
// 
/**\class HGCalPhotonIDValueMapProducer HGCalPhotonIDValueMapProducer.cc RecoEgamma/HGCalPhotonIDValueMapProducer/plugins/HGCalPhotonIDValueMapProducer.cc

 Description: [one line class summary]

 Implementation:
    [Notes on implementation]
*/
//
// Original Author:  Nicholas Charles Smith
//      Created:  Wed, 05 Apr 2017 12:17:43 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"


class HGCalPhotonIDValueMapProducer : public edm::stream::EDProducer<> {
  public:
    explicit HGCalPhotonIDValueMapProducer(const edm::ParameterSet&);
    ~HGCalPhotonIDValueMapProducer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginStream(edm::StreamID) override;
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endStream() override;

    // ----------member data ---------------------------
    edm::EDGetTokenT<reco::PhotonCollection> photonsToken_;
    std::map<const std::string, std::vector<float>> maps_;
};

HGCalPhotonIDValueMapProducer::HGCalPhotonIDValueMapProducer(const edm::ParameterSet& iConfig) :
  photonsToken_(consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons")))
{
  // Define here all the ValueMap names to output
  maps_["dummy"] = {};
  maps_["dummy2"] = {};

  for(auto&& kv : maps_) {
    produces<edm::ValueMap<float>>(kv.first);
  }
}


HGCalPhotonIDValueMapProducer::~HGCalPhotonIDValueMapProducer()
{
}


// ------------ method called to produce the data  ------------
void
HGCalPhotonIDValueMapProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle<reco::PhotonCollection> photonsH;
  iEvent.getByToken(photonsToken_, photonsH);
 
  // Clear previous map
  for(auto&& kv : maps_) kv.second.clear();

  for(size_t iPho=0; iPho<photonsH->size(); ++iPho) {
    const auto& pho = photonsH->at(iPho);

    // Fill here all the ValueMaps from their appropriate functions
    float var = pho.pt();
    maps_["dummy"].push_back(var);

    float var2 = pho.eta();
    maps_["dummy2"].push_back(var2);
  }

  for(auto&& kv : maps_) {
    // Check we didn't forget any values
    if ( kv.second.size() != photonsH->size() ) {
      throw cms::Exception("HGCalPhotonIDValueMapProducer") << "We have a miscoded value map producer, since the variable " << kv.first << " wasn't filled.";
    }
    // Do the filling
    auto out = std::make_unique<edm::ValueMap<float>>();
    edm::ValueMap<float>::Filler filler(*out);
    filler.insert(photonsH, kv.second.begin(), kv.second.end());
    filler.fill();
    // and put it into the event
    iEvent.put(std::move(out), kv.first);
  }
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
HGCalPhotonIDValueMapProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
HGCalPhotonIDValueMapProducer::endStream() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HGCalPhotonIDValueMapProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalPhotonIDValueMapProducer);
