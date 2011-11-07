
#ifndef CalibratedGsfElectronProducer_h
#define CalibratedGsfElectronProducer_h

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

class CalibratedGsfElectronProducer: public edm::EDProducer 
 {
  public:

    //static void fillDescriptions( edm::ConfigurationDescriptions & ) ;

    explicit CalibratedGsfElectronProducer( const edm::ParameterSet & ) ;
    virtual ~CalibratedGsfElectronProducer();
    virtual void produce( edm::Event &, const edm::EventSetup & ) ;

  private:

    edm::InputTag inputGsfElectrons ;
    std::string dataset ;
    bool isMC ;
    
 } ;

#endif
