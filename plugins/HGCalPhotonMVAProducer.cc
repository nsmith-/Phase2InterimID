// Original Author:  Nicholas Charles Smith

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
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "TMVA/Reader.h"

class HGCalPhotonMVAProducer : public edm::stream::EDProducer<> {
  public:
    explicit HGCalPhotonMVAProducer(const edm::ParameterSet&);
    ~HGCalPhotonMVAProducer() override;

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    void beginStream(edm::StreamID) override;
    void produce(edm::Event&, const edm::EventSetup&) override;
    void endStream() override;

    // ----------member data ---------------------------
    edm::EDGetTokenT<edm::View<reco::Photon>> photonsToken_;
    edm::EDGetTokenT<double> rhoToken_;
    edm::EDGetTokenT<edm::View<reco::GsfElectron>> electronsFromMultiClToken_;
    bool usePat_;
    // We'll need these if not using PAT
    edm::EDGetTokenT<reco::GsfElectronCollection> gedGsfElectronsToken_;
    edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
    edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
    std::map<const std::string, edm::EDGetTokenT<edm::ValueMap<float>>> maps_;

    static constexpr size_t NBarrelVars{17u};
    static constexpr size_t NEndcapVars{18u};
    std::array<float, NBarrelVars> barrelVars_;
    std::array<float, NEndcapVars> endcapVars_;

    std::unique_ptr<TMVA::Reader> barrelReader_;
    std::unique_ptr<TMVA::Reader> endcapReader_;
};

HGCalPhotonMVAProducer::HGCalPhotonMVAProducer(const edm::ParameterSet& iConfig) :
  photonsToken_(consumes<edm::View<reco::Photon>>(iConfig.getParameter<edm::InputTag>("photons"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoSrc"))),
  electronsFromMultiClToken_(consumes<edm::View<reco::GsfElectron>>(iConfig.getParameter<edm::InputTag>("electronsFromMultiCl"))),
  usePat_(iConfig.getParameter<bool>("usePAT")),
  gedGsfElectronsToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("gedGsfElectrons"))),
  conversionsToken_(consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conversions"))),
  beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot")))
{
  if ( not usePat_ ) {
    auto externalVars = iConfig.getParameter<edm::ParameterSet>("externalVariables");
    for (auto name : externalVars.getParameterNames()) {
      if ( externalVars.existsAs<edm::InputTag>(name) ) {
        maps_[name] = consumes<edm::ValueMap<float>>(externalVars.getParameter<edm::InputTag>(name));
      }
    }
  }

  barrelReader_.reset(new TMVA::Reader("!Color:Silent"));
  barrelReader_->AddVariable("full5x5_sigmaIetaIeta", &barrelVars_[0]);
  barrelReader_->AddVariable("full5x5_sigmaIetaIphi", &barrelVars_[1]);
  barrelReader_->AddVariable("full5x5_sigmaIphiIphi", &barrelVars_[2]);
  barrelReader_->AddVariable("etaWidth",              &barrelVars_[3]);
  barrelReader_->AddVariable("phiWidth",              &barrelVars_[4]);
  barrelReader_->AddVariable("full5x5_r9",            &barrelVars_[5]);
  barrelReader_->AddVariable("full5x5_s4",            &barrelVars_[6]);
  barrelReader_->AddVariable("hadronicOverEm",        &barrelVars_[7]);
  barrelReader_->AddVariable("chargedHadronIso",      &barrelVars_[8]);
  barrelReader_->AddVariable("neutralHadronIso",      &barrelVars_[9]);
  barrelReader_->AddVariable("photonIso",             &barrelVars_[10]);
  barrelReader_->AddVariable("hasPixelSeed",          &barrelVars_[11]);
  barrelReader_->AddVariable("scRawEnergy",           &barrelVars_[12]);
  barrelReader_->AddVariable("scEta",                 &barrelVars_[13]);
  barrelReader_->AddVariable("trkSumPt",              &barrelVars_[14]);
  barrelReader_->AddVariable("eVeto",                 &barrelVars_[15]);
  barrelReader_->AddVariable("rho",                   &barrelVars_[16]);
  barrelReader_->BookMVA("onlyone", iConfig.getParameter<edm::FileInPath>("barrelTrainingFile").fullPath());

  endcapReader_.reset(new TMVA::Reader("!Color:Silent"));
  endcapReader_->AddVariable("sigmaUU",              &endcapVars_[0]);
  endcapReader_->AddVariable("sigmaVV",              &endcapVars_[1]);
  endcapReader_->AddVariable("e4oEtot",              &endcapVars_[2]);
  endcapReader_->AddVariable("layerEfrac10",         &endcapVars_[3]);
  endcapReader_->AddVariable("layerEfrac90",         &endcapVars_[4]);
  endcapReader_->AddVariable("FHoverE",              &endcapVars_[5]);
  endcapReader_->AddVariable("measuredDepth",        &endcapVars_[6]);
  endcapReader_->AddVariable("depthCompatibility",   &endcapVars_[7]);
  endcapReader_->AddVariable("isoRing0",             &endcapVars_[8]);
  endcapReader_->AddVariable("isoRing1",             &endcapVars_[9]);
  endcapReader_->AddVariable("isoRing2",             &endcapVars_[10]);
  endcapReader_->AddVariable("isoRing3",             &endcapVars_[11]);
  endcapReader_->AddVariable("isoRing4",             &endcapVars_[12]);
  endcapReader_->AddVariable("scEnergy",             &endcapVars_[13]);
  endcapReader_->AddVariable("matchedTrackChi2",     &endcapVars_[14]);
  endcapReader_->AddVariable("matchedTrackHits",     &endcapVars_[15]);
  endcapReader_->AddVariable("matchedTrackLostHits", &endcapVars_[16]);
  endcapReader_->AddVariable("rho",                  &endcapVars_[17]);
  endcapReader_->BookMVA("onlyone", iConfig.getParameter<edm::FileInPath>("endcapTrainingFile").fullPath());

  produces<edm::ValueMap<float>>();
}


HGCalPhotonMVAProducer::~HGCalPhotonMVAProducer()
{
}


// ------------ method called to produce the data  ------------
void
HGCalPhotonMVAProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle<edm::View<reco::Photon>> photonsH;
  iEvent.getByToken(photonsToken_, photonsH);

  Handle<double> rhoH;
  iEvent.getByToken(rhoToken_, rhoH);

  Handle<edm::View<reco::GsfElectron>> ecalDrivenElectronsH;
  iEvent.getByToken(electronsFromMultiClToken_, ecalDrivenElectronsH);

  Handle<reco::GsfElectronCollection> gedGsfElectronsH;
  Handle<reco::ConversionCollection> conversionsH;
  Handle<reco::BeamSpot> beamspotH;
  if ( not usePat_ ) {
    iEvent.getByToken(gedGsfElectronsToken_, gedGsfElectronsH);
    iEvent.getByToken(conversionsToken_, conversionsH);
    iEvent.getByToken(beamspotToken_, beamspotH);
  }

  std::map<const std::string, Handle<ValueMap<float>>> handles;
  for(auto&& kv : maps_) {
    iEvent.getByToken(kv.second, handles[kv.first]);
  }

  std::vector<float> mvaValues(photonsH->size(), -1.);

  for(size_t iPho=0; iPho<photonsH->size(); ++iPho) {
    const auto& pho = photonsH->at(iPho);
    if ( usePat_ && reinterpret_cast<const pat::Photon*>(&pho) == nullptr ) {
      throw cms::Exception("HGCalPhotonMVAProducer") << "You lied to me, this isn't a pat::Photon";
    }
    auto externalvar = [&](const std::string& key) -> float {
      if ( usePat_ ) {
        return reinterpret_cast<const pat::Photon*>(&pho)->userFloat("hgcPhotonID:"+key);
      }
      if ( handles.find(key) != handles.end() and handles[key].isValid() ) {
        return handles[key]->get(photonsH.id(), iPho);
      }
      throw cms::Exception("HGCalPhotonMVAProducer") << "No such external var exists! " << key;
      return 0.f;
    };

    if(pho.isEB()) {
      barrelVars_[0]  = pho.full5x5_sigmaIetaIeta(); // full5x5_sigmaIetaIeta
      barrelVars_[1]  = pho.full5x5_showerShapeVariables().sigmaIetaIphi; // full5x5_sigmaIetaIphi
      barrelVars_[2]  = pho.full5x5_showerShapeVariables().sigmaIphiIphi; // full5x5_sigmaIphiIphi
      barrelVars_[3]  = pho.superCluster()->etaWidth(); // etaWidth
      barrelVars_[4]  = pho.superCluster()->phiWidth(); // phiWidth
      barrelVars_[5]  = pho.r9(); // full5x5_r9
      barrelVars_[6]  = pho.showerShapeVariables().e2x2/pho.showerShapeVariables().e5x5; // full5x5_s4
      barrelVars_[7]  = pho.hadronicOverEm(); // hadronicOverEm
      barrelVars_[8]  = pho.chargedHadronIso(); // chargedHadronIso
      barrelVars_[9]  = pho.neutralHadronIso(); // neutralHadronIso
      barrelVars_[10] = pho.photonIso(); //photonIso
      barrelVars_[11] = pho.hasPixelSeed(); //hasPixelSeed
      barrelVars_[12] = pho.superCluster()->rawEnergy(); //scRawEnergy
      barrelVars_[13] = pho.superCluster()->eta(); //scEta
      barrelVars_[14] = pho.trkSumPtSolidConeDR04(); //trkSumPtSolidConeDR04
      if ( usePat_ ) {
        barrelVars_[15] = reinterpret_cast<const pat::Photon*>(&pho)->passElectronVeto(); //conversionSafeElectronVeto
      } else {
        barrelVars_[15] = !ConversionTools::hasMatchedPromptElectron(pho.superCluster(), gedGsfElectronsH, conversionsH, beamspotH->position()); //conversionSafeElectronVeto
      }
      barrelVars_[16] = *rhoH; //rho

      mvaValues[iPho] = barrelReader_->EvaluateMVA("onlyone");
    }
    else {
      endcapVars_[0]  = externalvar("sigmaUU");
      endcapVars_[1]  = externalvar("sigmaVV");
      endcapVars_[2]  = externalvar("e4oEtot");
      endcapVars_[3]  = externalvar("layerEfrac10");
      endcapVars_[4]  = externalvar("layerEfrac90");
      endcapVars_[5]  = externalvar("seedEnergyFH") / externalvar("seedEnergyEE");
      endcapVars_[6]  = externalvar("measuredDepth");
      endcapVars_[7]  = externalvar("depthCompatibility");
      endcapVars_[8]  = externalvar("caloIsoRing0");
      endcapVars_[9]  = externalvar("caloIsoRing1");
      endcapVars_[10] = externalvar("caloIsoRing2");
      endcapVars_[11] = externalvar("caloIsoRing3");
      endcapVars_[12] = externalvar("caloIsoRing4");
      endcapVars_[13] = pho.superCluster()->rawEnergy();
      endcapVars_[14] = 1.e5; // matchedGsfChi2
      endcapVars_[15] = -1.; // matchedGsfHits
      endcapVars_[16] = 0.; // matchedGsfLostHits
      for(const auto& el : *ecalDrivenElectronsH) {
        if ( reco::deltaR(*el.superCluster(), *pho.superCluster()) > 0.01 ) continue;
        endcapVars_[14] = el.gsfTrack()->chi2(); // matchedGsfChi2
        endcapVars_[15] = el.gsfTrack()->hitPattern().trackerLayersWithMeasurement(); // matchedGsfHits
        endcapVars_[16] = el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS); // matchedGsfLostHits
        break;
      }
      endcapVars_[17] = *rhoH; // rho
      mvaValues[iPho] = endcapReader_->EvaluateMVA("onlyone");
    }
  }

  // Do the filling
  auto out = std::make_unique<edm::ValueMap<float>>();
  edm::ValueMap<float>::Filler filler(*out);
  filler.insert(photonsH, mvaValues.begin(), mvaValues.end());
  filler.fill();
  // and put it into the event
  iEvent.put(std::move(out), "");
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
HGCalPhotonMVAProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
HGCalPhotonMVAProducer::endStream() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HGCalPhotonMVAProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalPhotonMVAProducer);
