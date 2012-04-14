
#ifndef ElectronEnergyCalibrator_H
#define ElectronEnergyCalibrator_H

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "FWCore/Framework/interface/ESHandle.h"

class ElectronEnergyCalibrator
{
 public:

  ElectronEnergyCalibrator(std::string dataset, bool isAOD, bool isMC, bool updateEnergyError, bool debug) : dataset_(dataset),
   isAOD_(isAOD), isMC_(isMC), updateEnergyError_(updateEnergyError), debug_(debug) {}

  void correct(reco::GsfElectron &, const edm::Event&, const edm::EventSetup&);

 private:

  void computeNewEnergy( const reco::GsfElectron &, float r9, int run) ;
  void computeEpCombination( reco::GsfElectron & electron ) ;

  float newEnergy_ ;
  float newEnergyError_ ;
  
  math::XYZTLorentzVector newMomentum_ ;
  float errorTrackMomentum_ ;
  float finalMomentumError_ ;

  unsigned long long cacheIDTopo ;
  edm::ESHandle<CaloTopology> caloTopo ;
  
  std::string dataset_;
  bool isAOD_;
  bool isMC_;
  bool updateEnergyError_;
  bool debug_;
   
};

#endif




