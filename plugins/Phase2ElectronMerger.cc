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

#include "DataFormats/PatCandidates/interface/Electron.h"

class Phase2ElectronMerger : public edm::stream::EDProducer<> {
  public:
    explicit Phase2ElectronMerger(const edm::ParameterSet&);
    ~Phase2ElectronMerger();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginStream(edm::StreamID) override;
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endStream() override;

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    edm::EDGetTokenT<pat::ElectronCollection> barrelElectronToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> barrelIDToken_;
    StringCutObjectSelector<pat::Electron> barrelCut_;
    edm::EDGetTokenT<pat::ElectronCollection> endcapElectronToken_;
    edm::EDGetTokenT<edm::ValueMap<float>> endcapIDToken_;
    StringCutObjectSelector<pat::Electron> endcapCut_;
};

Phase2ElectronMerger::Phase2ElectronMerger(const edm::ParameterSet& iConfig):
  barrelElectronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("barrelElectrons"))),
  barrelIDToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("barrelID"))),
  barrelCut_(iConfig.getParameter<std::string>("barrelCut")),
  endcapElectronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("endcapElectrons"))),
  endcapIDToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("endcapID"))),
  endcapCut_(iConfig.getParameter<std::string>("endcapCut"))
{
  produces<pat::ElectronCollection>();
}


Phase2ElectronMerger::~Phase2ElectronMerger()
{
}

void
Phase2ElectronMerger::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle<pat::ElectronCollection> barrelElectronHandle;
  iEvent.getByToken(barrelElectronToken_, barrelElectronHandle);
  const auto barrelInputId = barrelElectronHandle.id();

  Handle<ValueMap<float>> barrelIDHandle;
  iEvent.getByToken(barrelIDToken_, barrelIDHandle);

  Handle<pat::ElectronCollection> endcapElectronHandle;
  iEvent.getByToken(endcapElectronToken_, endcapElectronHandle);
  const auto endcapInputId = endcapElectronHandle.id();

  Handle<ValueMap<float>> endcapIDHandle;
  iEvent.getByToken(endcapIDToken_, endcapIDHandle);

  std::unique_ptr<pat::ElectronCollection> electrons(new pat::ElectronCollection());

  for(size_t iEle=0; iEle<barrelElectronHandle->size(); ++iEle) {
    pat::Electron eleNew(barrelElectronHandle->at(iEle));
    eleNew.addUserFloat("mvaValue", barrelIDHandle->get(barrelInputId , iEle));
    if ( barrelCut_(eleNew) ) electrons->push_back(eleNew);
  }
  for(size_t iEle=0; iEle<endcapElectronHandle->size(); ++iEle) {
    pat::Electron eleNew(endcapElectronHandle->at(iEle));
    // Common access method
    eleNew.addUserFloat("mvaValue", endcapIDHandle->get(endcapInputId , iEle));
    if ( endcapCut_(eleNew) ) electrons->push_back(eleNew);
  }

  iEvent.put(std::move(electrons));
}


void
Phase2ElectronMerger::beginStream(edm::StreamID)
{
}


void
Phase2ElectronMerger::endStream() {
}


/*
void
Phase2ElectronMerger::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 

/*
void
Phase2ElectronMerger::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 

/*
void
Phase2ElectronMerger::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 

/*
void
Phase2ElectronMerger::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 

void
Phase2ElectronMerger::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(Phase2ElectronMerger);
