#include "RecoEgamma/Phase2InterimID/interface/HGCalIDTool.h"

#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"

#include "TPrincipal.h"
#include "FWCore/Utilities/interface/isFinite.h"

HGCalIDTool::HGCalIDTool(const edm::ParameterSet& conf, edm::ConsumesCollector& cc):
  eetok( cc.consumes<HGCRecHitCollection>(conf.getParameter<edm::InputTag>("HGCEEInput")) ),
  fhtok( cc.consumes<HGCRecHitCollection>(conf.getParameter<edm::InputTag>("HGCFHInput")) ),
  bhtok( cc.consumes<HGCRecHitCollection>(conf.getParameter<edm::InputTag>("HGCBHInput")) ),
  hgcPftok( cc.consumes<reco::PFRecHitCollection>(conf.getParameter<edm::InputTag>("HGCPFRecHits")) )
{
  withPileup_ = conf.getParameter<bool>("withPileup");

  // parameters
  mip_ = 0.0000551;
  minenergy_ = 4.;
  rmax_ = 100.; // no transverse limitation for no PU case
  if (withPileup_) rmax_ = 1.5*2.27;
  hovereConesize_ = 0.05;
    
  // HGCAL average medium
  criticalEnergy_ = 0.00536; // in GeV
  radiationLength_ = 0.968; // in cm
    
  // longitudinal parameters
  // mean values
  // shower max <T> = t0 + t1*lny
  // <alpha> = alpha0 + alpha1*lny
  // shower average = alpha/beta  
  meant0_ = -1.396;
  meant1_ = 1.007;
  meanalpha0_ = -0.0433;
  meanalpha1_ = 0.540;
  // sigmas
  // sigma(lnT) = 1 /sigmalnt0 + sigmalnt1*lny; 
  // sigma(lnalpha) = 1 /sigmalnt0 + sigmalnt1*lny; 
  sigmalnt0_ = -2.506;
  sigmalnt1_ = 1.245;
  sigmalnalpha0_ = -0.08442;
  sigmalnalpha1_ = 0.7904;
  // corr(lnalpha,lnt) = corrlnalpha0_+corrlnalphalnt1_*y
  corrlnalphalnt0_ = 0.7858;
  corrlnalphalnt1_ = -0.0232;
    
  // cut values, to be moved as configurable parameters
  cutStartPosition_ = 325.5;
  cutSigmaetaeta_ = 0.0055;
  if (withPileup_) cutSigmaetaeta_ = 0.00480;
  cutHoverem_ = 0.003;
  if (withPileup_) cutHoverem_ = 0.065;
  cutLengthCompatibility_ = 4.0;
}

void HGCalIDTool::getEvent(const edm::Event& ev) {
  rhtools_.getEvent(ev);
  edm::Handle<HGCRecHitCollection> temp;
  ev.getByToken(eetok, temp);
  eerh_ = temp.product();
  ev.getByToken(fhtok, temp);
  fhrh_ = temp.product();
  ev.getByToken(bhtok, temp);
  bhrh_ = temp.product();
  edm::Handle<reco::PFRecHitCollection> temp2;
  ev.getByToken(hgcPftok, temp2);
  pfRecHits_ = temp2.product();
}

void HGCalIDTool::getEventSetup(const edm::EventSetup& es) {
  rhtools_.getEventSetup(es);
}

bool HGCalIDTool::setClusterPtr(const reco::CaloCluster * cluster)
{
  // Check that actually HGCal
  int det_group = cluster->hitsAndFractions()[0].first.det();
  int detector = cluster->hitsAndFractions()[0].first.subdetId();
  if ( (detector==HcalEndcap || det_group == DetId::Forward) ) {
    cluster_ = cluster;
    showerPos_ = math::XYZPoint(0.,0.,0.);
    showerDir_ = math::XYZVector(0.,0.,0.);
    rebuildRecHitFractions();
    calculateShowerPositionAndAxis();
    return true;
  }

  cluster_ = nullptr;
  return false;
}

double HGCalIDTool::getHadronFraction()
{
  if ( cluster_ == nullptr ) {
    throw cms::Exception("HGCalIDTool") << "Please call setClusterPtr() first! (and make sure it returns true, otherwise cluster is not from HGCAL)" << std::endl;
  }

  float energy=0.f, energyHad=0.f;
  const auto& hits = cluster_->hitsAndFractions();
  for( const auto& hit : hits ) {
    const auto& id = hit.first;
    const float fraction = hit.second;
    if( id.det() == DetId::Forward ) {
      switch( id.subdetId() ) {
      case HGCEE:
        energy += eerh_->find(id)->energy()*fraction;
        break;
      case HGCHEF:
        {
          const float temp = fhrh_->find(id)->energy();
          energy += temp*fraction;
          energyHad += temp*fraction;
        }
        break;
      default:
        throw cms::Exception("HGCalHGCalIDTool")
          << " Cluster contains hits that are not from HGCal! " << std::endl;
      }
    } else if ( id.det() == DetId::Hcal && id.subdetId() == HcalEndcap ) {
      const float temp = bhrh_->find(id)->energy();
      energy += temp*fraction;
      energyHad += temp*fraction;
    } else {
      throw cms::Exception("HGCalHGCalIDTool")
        << " Cluster contains hits that are not from HGCal! " << std::endl;
    }    
  }
  float fraction = -1.f;
  if( energy > 0.f ) {
    fraction = energyHad/energy;
  }
  return fraction;
}

math::XYZPoint HGCalIDTool::getStartPosition()
{
  if ( cluster_ == nullptr ) {
    throw cms::Exception("HGCalIDTool") << "Please call setClusterPtr() first! (and make sure it returns true, otherwise cluster is not from HGCAL)" << std::endl;
  }

  math::XYZPoint firstPos;
  double zmin = 10000.0;  
  for (unsigned int ih=0;ih<recHitFractions_.size();++ih) {
    const auto& refhit = recHitFractions_[ih].recHitRef();
    const auto& pos = refhit->position();
    const DetId & id_(refhit->detId()) ;
    if (id_.det()==DetId::Forward) {      
      if (std::abs(pos.z())<zmin) {
	firstPos = pos;
	zmin = std::abs(pos.z());
      }
    }
  }

  // refine the first position estimation, taking the max energy in the first layer 
  double maxfirstenergy=0.; 
  for (unsigned int ih=0;ih<recHitFractions_.size();++ih) {
    const auto& refhit = recHitFractions_[ih].recHitRef();
    const auto& pos = refhit->position();
    const DetId & id_(refhit->detId());
    if (id_.det()==DetId::Forward) {
      if (std::abs(pos.z()) != zmin) continue;
      if (refhit->energy() > maxfirstenergy) {
	firstPos = pos;
	maxfirstenergy = refhit->energy();
      }
    }
  }
    
  // finally refine firstPos x and y using the meaured direction 
  // 
  double lambda = (firstPos-showerPos_).z()/showerDir_.z();
  math::XYZPoint extraPos = showerPos_ + lambda*showerDir_;	 
  firstPos = extraPos;

  return firstPos;

}

double HGCalIDTool::getSigmaEtaEta()
{
  if ( cluster_ == nullptr ) {
    throw cms::Exception("HGCalIDTool") << "Please call setClusterPtr() first! (and make sure it returns true, otherwise cluster is not from HGCAL)" << std::endl;
  }

  double sigmaetaeta=0., sumnrj=0.;
  math::XYZPoint firstPos = getStartPosition();
    
  for (unsigned int ih=0;ih<recHitFractions_.size();++ih) {
    const auto& refhit = recHitFractions_[ih].recHitRef();
    const DetId & id_ = refhit->detId() ;       
    if (id_.det()==DetId::Forward) {
      math::XYZPoint cellPos(refhit->position());
      math::XYZVector radius, longitudinal, transverse;
      radius = cellPos - firstPos;
      // distances in local coordinates
      longitudinal =  (radius.Dot(showerDir_))*showerDir_.unit()/showerDir_.R();
      transverse = radius - longitudinal;
      // apply energy cut cut
      if (!withPileup_ || refhit->energy()>minenergy_*mip_) {
	// simple transversal cut, later can refine as function of depth
	if (!withPileup_ || transverse.R() < rmax_) {
	  const double deta = (cellPos.eta()-showerPos_.eta());
	  sigmaetaeta += deta*deta*refhit->energy();
	  sumnrj += refhit->energy();
	}
      }
    }
  }

  sigmaetaeta /= sumnrj;
  sigmaetaeta = sqrt(sigmaetaeta);

  // now correct the eta dependency
  double feta;
  constexpr double feta_0 = 0.00964148 - 0.01078431*1.5 + 0.00495703*1.5*1.5;
  const double clu_eta = std::abs(cluster_->eta()); 
  feta = 0.00964148 - clu_eta*(0.0107843 - 0.00495703*clu_eta);
  sigmaetaeta *= feta_0 / feta ;

  return sigmaetaeta;

}

double HGCalIDTool::getLengthCompatibility()
{
  if ( cluster_ == nullptr ) {
    throw cms::Exception("HGCalIDTool") << "Please call setClusterPtr() first! (and make sure it returns true, otherwise cluster is not from HGCAL)" << std::endl;
  }

  double lengthCompatibility=0., predictedLength=0., predictedSigma=0.;
	
  // shower length	 
  const double length =  (showerPos_ - getStartPosition()).R();
  const double cluster_emEnergy = cluster_->energy() * (1.- getHadronFraction());
  const double lny = cluster_emEnergy/criticalEnergy_>1. ? std::log(cluster_emEnergy/criticalEnergy_) : 0.;

  // inject here parametrization results
  const double meantmax = meant0_ + meant1_*lny;
  const double meanalpha = meanalpha0_ + meanalpha1_*lny;
  const double sigmalntmax = 1.0 / (sigmalnt0_+sigmalnt1_*lny);
  const double sigmalnalpha = 1.0 / (sigmalnalpha0_+sigmalnalpha1_*lny);
  const double corrlnalphalntmax = corrlnalphalnt0_+corrlnalphalnt1_*lny;
  
  const double invbeta = meantmax/(meanalpha-1.);
  predictedLength = meanalpha*invbeta;
  predictedLength *= radiationLength_;
  
  double sigmaalpha = meanalpha*sigmalnalpha;
  if (sigmaalpha<0.) sigmaalpha = 1.;
  double sigmatmax = meantmax*sigmalntmax;
  if (sigmatmax<0.) sigmatmax = 1.;
  
  predictedSigma = sigmalnalpha*sigmalnalpha/((meanalpha-1.)*(meanalpha-1.));
  predictedSigma += sigmalntmax*sigmalntmax;
  predictedSigma -= 2*sigmalnalpha*sigmalntmax*corrlnalphalntmax/(meanalpha-1.);
  predictedSigma = predictedLength*sqrt(predictedSigma);
  
  lengthCompatibility = (predictedLength-length)/predictedSigma;
  
  return lengthCompatibility;
  
}

void HGCalIDTool::rebuildRecHitFractions()
{
  recHitFractions_.clear();

  const auto& hits = cluster_->hitsAndFractions();
  for( const auto& hit : hits ) {
    const auto& detId = hit.first;
    const float fraction = hit.second;
    size_t pos{0};
    for ( const auto& recHit : *pfRecHits_ ) {
      if ( recHit.detId() == detId ) break;
      pos++;
    }
    if ( pos < pfRecHits_->size() ) {
      recHitFractions_.emplace_back(reco::PFRecHitRef(pfRecHits_, pos), fraction);
    }
    else {
      throw cms::Exception("MissingRecHit") << "Cluster has missing PFRecHit. cluster eta " << cluster_->eta() << " rechits size " << pfRecHits_->size();
    }
  }
}

void HGCalIDTool::calculateShowerPositionAndAxis()
{ 
  if( !cluster_->seed() ) {
    throw cms::Exception("ClusterWithNoSeed") << "Found a cluster with no seed: " << cluster_;
  }
  auto pca_ = std::make_unique<TPrincipal>(3, "D");

  double cl_energy = 0;  
  double max_e = 0.0;
  double avg_time = 0.0;
  double time_norm = 0.0;
  PFLayer::Layer max_e_layer = PFLayer::NONE;
  reco::PFRecHitRef refseed;  
  double pcavars[3];  

  for( const reco::PFRecHitFraction& rhf : recHitFractions_ ) {
    const reco::PFRecHitRef& refhit = rhf.recHitRef();
    double rh_energy = refhit->energy();
    double rh_time = refhit->time();
    cl_energy += rh_energy * rhf.fraction();
    if( rh_time > 0.0 ) { // time == -1 means no measurement
      // all times are offset by one nanosecond in digitizer
      // remove that here so all times of flight
      // are with respect to (0,0,0)
      avg_time += (rh_time - 1.0); 
      time_norm += 1.0;
    }
    if( rh_energy > max_e ) {
      max_e = rh_energy;
      max_e_layer = rhf.recHitRef()->layer();
    }  
    if( refhit->detId() == cluster_->seed() ) refseed = refhit;
    const double rh_fraction = rhf.fraction();
    rh_energy = refhit->energy()*rh_fraction;
    if( edm::isNotFinite(rh_energy) ) {
      throw cms::Exception("PFClusterAlgo")
	<<"rechit " << refhit->detId() << " has a NaN energy... " 
	<< "The input of the particle flow clustering seems to be corrupted.";
    }    
    pcavars[0] = refhit->position().x();
    pcavars[1] = refhit->position().y();
    pcavars[2] = refhit->position().z();     
    int nhit = int( rh_energy*100 ); // put rec_hit energy in units of 10 MeV

    for( int i = 0; i < nhit; ++i ) {
      pca_->AddRow(pcavars);
    }
      
  }

  // calculate the position
  pca_->MakePrincipals();
  const TVectorD& means = *(pca_->GetMeanValues());
  const TMatrixD& eigens = *(pca_->GetEigenVectors());
  
  math::XYZPoint  barycenter(means[0],means[1],means[2]);
  math::XYZVector axis(eigens(0,0),eigens(1,0),eigens(2,0));

  if( time_norm > 0.0 ) {
    avg_time = avg_time/time_norm;
  } else {
    avg_time = std::numeric_limits<double>::min();
  }

  if( axis.z()*barycenter.z() < 0.0 ) {
    axis = math::XYZVector(-eigens(0,0),-eigens(1,0),-eigens(2,0));
  }
  
  (void) max_e_layer;
  showerPos_ = barycenter;
  showerDir_ = axis;
}
