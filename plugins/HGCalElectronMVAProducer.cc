// Authors: F. Beaudette, A. Lobanov, C. Ochando
// Ported by N. Smith
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

#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "TMVA/Reader.h"

class HGCalElectronMVAProducer : public edm::stream::EDProducer<> {
  public:
    explicit HGCalElectronMVAProducer(const edm::ParameterSet&);
    ~HGCalElectronMVAProducer() override;

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    void beginStream(edm::StreamID) override;
    void produce(edm::Event&, const edm::EventSetup&) override;
    void endStream() override;

    void initReader(TMVA::Reader& reader, const std::string& file, bool barrel);

    // ----------member data ---------------------------
    edm::EDGetTokenT<edm::View<reco::GsfElectron>> electronsToken_;
    bool usePat_;
    std::map<const std::string, edm::EDGetTokenT<edm::ValueMap<float>>> maps_;

    std::unique_ptr<TMVA::Reader> barrelLowPtReader_;
    std::unique_ptr<TMVA::Reader> barrelHighPtReader_;
    std::unique_ptr<TMVA::Reader> endcapLowPtReader_;
    std::unique_ptr<TMVA::Reader> endcapHighPtReader_;

    // Variables
    float tmva_ele_kfhits;
    float tmva_ele_gsfhits;
    float tmva_ele_kfchi2;
    float tmva_ele_gsfchi2;
    float tmva_ele_fbrem;

    // E-p matching
    float tmva_ele_eelepout;
    float tmva_ele_deltaetaele;
    float tmva_ele_deltaphiele;

    // HGC-specific
    float tmva_ele_sigmauu;
    float tmva_ele_sigmavv;
    float tmva_ele_hgc_eigenvalues1;
    float tmva_ele_hgc_eigenvalues2;
    float tmva_ele_hgc_eigenvalues3;
    float tmva_ele_hgc_sigmas1;
    float tmva_ele_hgc_sigmas2;
    float tmva_ele_hgc_sigmas3;
    float tmva_ele_hgc_pcaAxisX;
    float tmva_ele_hgc_pcaAxisY;
    float tmva_ele_hgc_pcaAxisZ;
    float tmva_ele_hgc_pcaPosX;
    float tmva_ele_hgc_pcaPosY;
    float tmva_ele_hgc_pcaPosZ;
    float tmva_ele_hgc_FHoverEE;
    float tmva_ele_hgc_nlay;
    float tmva_ele_hgc_firstlay;
    float tmva_ele_hgc_lastlay;
    float tmva_ele_hgc_EE4overEE;
    float tmva_ele_hgc_layEfrac10;
    float tmva_ele_hgc_layEfrac90;
    float tmva_ele_hgc_depthCompat;

    float tmva_ele_oldsigmaietaieta;
    float tmva_ele_oldsigmaiphiiphi;
    float tmva_ele_oldcircularity;
    float tmva_ele_oldr9;
    float tmva_ele_scletawidth;
    float tmva_ele_sclphiwidth;
    float tmva_ele_he;
    float tmva_ele_ep;
    float tmva_ele_eseedpout;
    float tmva_ele_deltaetaseed;
    float tmva_ele_deltaphiseed;

    float tmva_dummy_spectator;

    float pTLimit_;
};

HGCalElectronMVAProducer::HGCalElectronMVAProducer(const edm::ParameterSet& iConfig) :
  electronsToken_(consumes<edm::View<reco::GsfElectron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
  usePat_(iConfig.getParameter<bool>("usePAT"))
{
  if ( not usePat_ ) {
    auto externalVars = iConfig.getParameter<edm::ParameterSet>("externalVariables");
    for (auto name : externalVars.getParameterNames()) {
      if ( externalVars.existsAs<edm::InputTag>(name) ) {
        maps_[name] = consumes<edm::ValueMap<float>>(externalVars.getParameter<edm::InputTag>(name));
      }
    }
  }

  barrelLowPtReader_= std::make_unique<TMVA::Reader>("!Color:Silent");
  barrelHighPtReader_= std::make_unique<TMVA::Reader>("!Color:Silent");
  endcapLowPtReader_= std::make_unique<TMVA::Reader>("!Color:Silent");
  endcapHighPtReader_= std::make_unique<TMVA::Reader>("!Color:Silent");
  std::string barrelLowPtFile(iConfig.getParameter<edm::FileInPath>("barrelLowPt").fullPath());
  std::string barrelHighPtFile(iConfig.getParameter<edm::FileInPath>("barrelHighPt").fullPath());
  std::string endcapLowPtFile(iConfig.getParameter<edm::FileInPath>("endcapLowPt").fullPath());
  std::string endcapHighPtFile(iConfig.getParameter<edm::FileInPath>("endcapHighPt").fullPath());
  initReader(*barrelLowPtReader_, barrelLowPtFile,true);
  initReader(*barrelHighPtReader_,barrelHighPtFile,true);
  initReader(*endcapLowPtReader_, endcapLowPtFile,false);
  initReader(*endcapHighPtReader_,endcapHighPtFile,false);

  pTLimit_ = 20.;

  produces<edm::ValueMap<float>>();
}


HGCalElectronMVAProducer::~HGCalElectronMVAProducer()
{
}


// ------------ method called to produce the data  ------------
void
HGCalElectronMVAProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle<edm::View<reco::GsfElectron>> electronsH;
  iEvent.getByToken(electronsToken_, electronsH);

  std::map<const std::string, Handle<ValueMap<float>>> handles;
  for(auto&& kv : maps_) {
    iEvent.getByToken(kv.second, handles[kv.first]);
  }

  std::vector<float> mvaValues(electronsH->size(), -1.);

  for(size_t iEl=0; iEl<electronsH->size(); ++iEl) {
    const auto& electron = electronsH->at(iEl);
    if ( usePat_ && reinterpret_cast<const pat::Electron*>(&electron) == nullptr ) {
      throw cms::Exception("HGCalElectronMVAProducer") << "You lied to me, this isn't a pat::Electron";
    }
    auto externalvar = [&](const std::string& key) -> float {
      if ( usePat_ ) {
        return reinterpret_cast<const pat::Electron*>(&electron)->userFloat("hgcElectronID:"+key);
      }
      if ( handles.find(key) != handles.end() and handles[key].isValid() ) {
        return handles[key]->get(electronsH.id(), iEl);
      }
      throw cms::Exception("HGCalElectronMVAProducer") << "No such external var exists! " << key;
      return 0.f;
    };

    float bdtValue = -1.;

    // Track-based
    tmva_ele_fbrem = electron.fbrem();
    tmva_ele_gsfhits = electron.gsfTrack()->hitPattern().trackerLayersWithMeasurement();
    tmva_ele_gsfchi2 = electron.gsfTrack()->normalizedChi2();

    reco::TrackRef myTrackRef = electron.closestCtfTrackRef();
    bool validKF = myTrackRef.isAvailable() && myTrackRef.isNonnull();
    tmva_ele_kfhits = (validKF) ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : -1.;
    tmva_ele_kfchi2 = (validKF) ? myTrackRef->normalizedChi2() : -1.;

    if(electron.isEB()) {


      tmva_ele_oldsigmaietaieta = electron.full5x5_sigmaIetaIeta();
      tmva_ele_oldsigmaiphiiphi = electron.full5x5_sigmaIphiIphi();
      tmva_ele_oldcircularity = ( electron.showerShape().e5x5!=0.) ? 1.-electron.showerShape().e1x5/electron.showerShape().e5x5 : -1. ;
      tmva_ele_oldr9 = electron.full5x5_r9();
      tmva_ele_scletawidth = electron.superCluster()->etaWidth();
      tmva_ele_sclphiwidth = electron.superCluster()->phiWidth();
      tmva_ele_he = electron.hcalOverEcal();
      tmva_ele_ep = electron.eSuperClusterOverP();
      tmva_ele_eseedpout = electron.eSeedClusterOverPout();
      tmva_ele_deltaetaseed = electron.deltaEtaSeedClusterTrackAtCalo();
      tmva_ele_deltaphiseed = electron.deltaPhiSeedClusterTrackAtCalo();

      bdtValue = (electron.pt() < pTLimit_ ) ? barrelLowPtReader_->EvaluateMVA("BDTSimpleCat") : barrelHighPtReader_->EvaluateMVA("BDTSimpleCat");
    }
    else if (externalvar("sigmaUU") == 0.) {
      bdtValue = -2.;
    }
    else {
      tmva_ele_hgc_depthCompat = externalvar("depthCompatibility");



      // Track-matching
      tmva_ele_eelepout = electron.eEleClusterOverPout();
      tmva_ele_deltaetaele = electron.deltaEtaEleClusterTrackAtCalo();
      tmva_ele_deltaphiele =  electron.deltaPhiEleClusterTrackAtCalo();

      // Cluster shapes
      // PCA related
      bool goodEle = externalvar("sigmaVV") != 0.;
      tmva_ele_hgc_eigenvalues1 = (goodEle) ? externalvar("pcaEig1") : -1.;
      tmva_ele_hgc_eigenvalues2 = (goodEle) ? externalvar("pcaEig2") : -1.;
      tmva_ele_hgc_eigenvalues3 = (goodEle) ? externalvar("pcaEig3") : -1.;

      tmva_ele_hgc_sigmas1 = (goodEle) ?  externalvar("pcaSig1") : -1. ;
      tmva_ele_hgc_sigmas2 = (goodEle) ?  externalvar("pcaSig2") : -1. ;
      tmva_ele_hgc_sigmas3 = (goodEle) ?  externalvar("pcaSig3") : -1. ;
      tmva_ele_hgc_pcaAxisX = (goodEle) ?  externalvar("pcaAxisX") : -1.;
      tmva_ele_hgc_pcaAxisY = (goodEle) ?  externalvar("pcaAxisY") : -1.;
      tmva_ele_hgc_pcaAxisZ = (goodEle) ?  externalvar("pcaAxisZ") : -1.;
      tmva_ele_hgc_pcaPosX = (goodEle) ?  externalvar("pcaPositionX") : -1.;
      tmva_ele_hgc_pcaPosY = (goodEle) ?  externalvar("pcaPositionY") : -1.;
      tmva_ele_hgc_pcaPosZ = (goodEle) ?  externalvar("pcaPositionZ") : -1.;

      // transverse shapes
      tmva_ele_sigmauu = externalvar("sigmaUU");
      tmva_ele_sigmavv = externalvar("sigmaVV");

      // long profile
      tmva_ele_hgc_nlay = externalvar("nLayers");
      tmva_ele_hgc_firstlay = externalvar("firstLayer");
      tmva_ele_hgc_lastlay = externalvar("lastLayer");
      tmva_ele_hgc_EE4overEE = externalvar("e4oEtot");
      tmva_ele_hgc_layEfrac10 = externalvar("layerEfrac10");
      tmva_ele_hgc_layEfrac90 = externalvar("layerEfrac90");
      tmva_ele_hgc_FHoverEE = ( externalvar("ecEnergyEE") !=0. ) ? externalvar("ecEnergyFH") / externalvar("ecEnergyEE") : -1.;

      bdtValue = (electron.pt() < pTLimit_ ) ? endcapLowPtReader_->EvaluateMVA("BDTSimpleCat") : endcapHighPtReader_->EvaluateMVA("BDTSimpleCat");
    }

    mvaValues[iEl] = bdtValue;
  }

  // Do the filling
  auto out = std::make_unique<edm::ValueMap<float>>();
  edm::ValueMap<float>::Filler filler(*out);
  filler.insert(electronsH, mvaValues.begin(), mvaValues.end());
  filler.fill();
  // and put it into the event
  iEvent.put(std::move(out), "");
}

void HGCalElectronMVAProducer::initReader(TMVA::Reader& reader, const std::string& file, bool barrel) {
  reader.AddVariable("ele_kfhits",&tmva_ele_kfhits);
  reader.AddVariable("ele_gsfhits",&tmva_ele_gsfhits);
  reader.AddVariable("ele_kfchi2",&tmva_ele_kfchi2);
  reader.AddVariable("ele_gsfchi2",&tmva_ele_gsfchi2);
  reader.AddVariable("ele_fbrem",&tmva_ele_fbrem);

  if (!barrel) {
    // E-p matching
    reader.AddVariable("ele_eelepout",&tmva_ele_eelepout);
    reader.AddVariable("ele_deltaetaele",&tmva_ele_deltaetaele);
    reader.AddVariable("ele_deltaphiele",&tmva_ele_deltaetaele);
    // Shower Shapes
    reader.AddVariable("ele_sigmauu",&tmva_ele_sigmauu);
    reader.AddVariable("ele_sigmavv",&tmva_ele_sigmavv);
    reader.AddVariable("ele_hgc_eigenvalues1",&tmva_ele_hgc_eigenvalues1);
    reader.AddVariable("ele_hgc_eigenvalues2",&tmva_ele_hgc_eigenvalues2);
    reader.AddVariable("ele_hgc_eigenvalues3",&tmva_ele_hgc_eigenvalues3);
    reader.AddVariable("ele_hgc_sigmas1",&tmva_ele_hgc_sigmas1);
    reader.AddVariable("ele_hgc_sigmas2",&tmva_ele_hgc_sigmas2);
    reader.AddVariable("ele_hgc_sigmas3",&tmva_ele_hgc_sigmas3);
    reader.AddVariable("ele_hgc_pcaAxisX",&tmva_ele_hgc_pcaAxisX);
    reader.AddVariable("ele_hgc_pcaAxisY",&tmva_ele_hgc_pcaAxisY);
    reader.AddVariable("ele_hgc_pcaAxisZ",&tmva_ele_hgc_pcaAxisZ);
    reader.AddVariable("ele_hgc_pcaPosX",&tmva_ele_hgc_pcaPosX);
    reader.AddVariable("ele_hgc_pcaPosY",&tmva_ele_hgc_pcaPosY);
    reader.AddVariable("ele_hgc_pcaPosZ",&tmva_ele_hgc_pcaPosZ);
    reader.AddVariable("ele_hgc_FHoverEE",&tmva_ele_hgc_FHoverEE);
    reader.AddVariable("ele_hgc_nlay",&tmva_ele_hgc_nlay);
    reader.AddVariable("ele_hgc_firstlay",&tmva_ele_hgc_firstlay);
    reader.AddVariable("ele_hgc_lastlay",&tmva_ele_hgc_lastlay);
    reader.AddVariable("ele_hgc_EE4overEE",&tmva_ele_hgc_EE4overEE);
    reader.AddVariable("ele_hgc_layEfrac10",&tmva_ele_hgc_layEfrac10);
    reader.AddVariable("ele_hgc_layEfrac90",&tmva_ele_hgc_layEfrac90);
    reader.AddVariable("ele_hgc_depthCompat",&tmva_ele_hgc_depthCompat);

    reader.AddSpectator("ele_pT",&tmva_dummy_spectator);
    reader.AddSpectator("ele_eta",&tmva_dummy_spectator);
    reader.AddSpectator("Nvtx",&tmva_dummy_spectator);
    reader.AddSpectator("ele_ET",&tmva_dummy_spectator);
    reader.AddSpectator("NPU",&tmva_dummy_spectator);
    reader.AddSpectator("PU_density",&tmva_dummy_spectator);
  } else {
    reader.AddVariable("ele_oldsigmaietaieta",&tmva_ele_oldsigmaietaieta);
    reader.AddVariable("ele_oldsigmaiphiiphi",&tmva_ele_oldsigmaiphiiphi);
    reader.AddVariable("ele_oldcircularity",&tmva_ele_oldcircularity);
    reader.AddVariable("ele_oldr9",&tmva_ele_oldr9);
    reader.AddVariable("ele_scletawidth",&tmva_ele_scletawidth);
    reader.AddVariable("ele_sclphiwidth",&tmva_ele_sclphiwidth);
    reader.AddVariable("ele_he",&tmva_ele_he);
    reader.AddVariable("ele_ep",&tmva_ele_ep);
    reader.AddVariable("ele_eseedpout",&tmva_ele_eseedpout);
    reader.AddVariable("ele_deltaetaseed",&tmva_ele_deltaetaseed);
    reader.AddVariable("ele_deltaphiseed",&tmva_ele_deltaphiseed);

    reader.AddSpectator("ele_pT",&tmva_dummy_spectator);
    reader.AddSpectator("ele_eta",&tmva_dummy_spectator);
    reader.AddSpectator("scl_eta",&tmva_dummy_spectator);
    reader.AddSpectator("Nvtx",&tmva_dummy_spectator);
    reader.AddSpectator("PU_density",&tmva_dummy_spectator);
    reader.AddSpectator("NPU",&tmva_dummy_spectator);
  }
  reader.BookMVA("BDTSimpleCat", file);
}
// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
HGCalElectronMVAProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
HGCalElectronMVAProducer::endStream() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HGCalElectronMVAProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalElectronMVAProducer);
