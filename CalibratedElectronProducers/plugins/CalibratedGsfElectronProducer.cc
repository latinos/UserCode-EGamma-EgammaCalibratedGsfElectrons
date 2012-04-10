// -*- C++ -*-
//
// Package:    EgammaElectronProducers
// Class:      CalibratedGsfElectronProducer
//
/**\class CalibratedGsfElectronProducer RecoEgamma/ElectronProducers/src/CalibratedGsfElectronProducer.cc

 Description: EDProducer of GsfElectron objects

 Implementation:
     <Notes on implementation>
*/

#include "EgammaCalibratedGsfElectrons/CalibratedElectronProducers/plugins/CalibratedGsfElectronProducer.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "EgammaCalibratedGsfElectrons/CalibratedElectronAlgos/interface/ElectronEnergyCalibrator.h"

#include <iostream>

using namespace edm ;
using namespace std ;
using namespace reco ;

CalibratedGsfElectronProducer::CalibratedGsfElectronProducer( const edm::ParameterSet & cfg )
// : GsfElectronBaseProducer(cfg)
 {

  produces<GsfElectronCollection>();

  inputGsfElectrons = cfg.getParameter<edm::InputTag>("inputGsfElectronsTag");
  dataset = cfg.getParameter<std::string>("inputDataset");
  isMC = cfg.getParameter<bool>("isMC");
  
  //basic checks
  if (isMC&&dataset!="Summer11")
   { throw cms::Exception("CalibratedgsfElectronProducer|ConfigError")<<"Unknown MC dataset" ; }
  if (!isMC&&(dataset!="Prompt"&&dataset!="ReReco"))
   { throw cms::Exception("CalibratedgsfElectronProducer|ConfigError")<<"Unknown Data dataset" ; }
 }
 
CalibratedGsfElectronProducer::~CalibratedGsfElectronProducer()
 {}

void CalibratedGsfElectronProducer::produce( edm::Event & event, const edm::EventSetup & setup )
 {

  edm::Handle<reco::GsfElectronCollection> gsfElectrons ;
  event.getByLabel(inputGsfElectrons,gsfElectrons) ;

  std::auto_ptr<GsfElectronCollection> electrons( new GsfElectronCollection ) ;
  const GsfElectronCollection * oldElectrons = gsfElectrons.product() ;

  GsfElectronCollection::const_iterator electron ;
  GsfElectronCollection::iterator ele ;
  // first clone the initial collection
  for
   ( electron = oldElectrons->begin() ;
     electron != oldElectrons->end() ;
     ++electron )
   {     
      electrons->push_back(*electron->clone()) ;
   }
  
  ElectronEnergyCalibrator theEnCorrector(dataset, isMC);

  for
   ( ele = electrons->begin() ;
     ele != electrons->end() ;
     ++ele )
   {     
      // energy calibration for ecalDriven electrons
      if (ele->core()->ecalDrivenSeed()) {        
        theEnCorrector.correct(*ele, event, setup);
      }
      else {
        //std::cout << "[CalibratedGsfElectronProducer] is tracker driven only!!" << std::endl;
      }
   }
   event.put(electrons) ;

 }


#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
DEFINE_FWK_MODULE(CalibratedGsfElectronProducer);
