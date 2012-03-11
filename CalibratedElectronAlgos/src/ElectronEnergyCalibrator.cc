

#include "EgammaCalibratedGsfElectrons/CalibratedElectronAlgos/interface/ElectronEnergyCalibrator.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include <CLHEP/Random/RandGaussQ.h>
#include <CLHEP/Random/Random.h>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Utilities/interface/Exception.h"

/****************************************************************************
 *
 * Propagate SC calibration from Zee fit to the electrons
 *
 ****************************************************************************/

using namespace edm;

void ElectronEnergyCalibrator::correct
 ( reco::GsfElectron & electron, const edm::Event& event, const edm::EventSetup& eventSetup)

 {

  // compute r9
  bool validEcalRecHits=true;
  Handle<EcalRecHitCollection> barrelHitHandle;
  EcalRecHitCollection barrelRecHits;
  InputTag ebhits("ecalRecHit","EcalRecHitsEB");
  event.getByLabel(ebhits, barrelHitHandle);
  if (!barrelHitHandle.isValid()) {
//    edm::LogError("PhotonProducer") << "Error! Can't get the product "<<barrelEcalHits_.label();
    validEcalRecHits=false; 
  }
  if (  validEcalRecHits)  barrelRecHits = *(barrelHitHandle.product());
 
  Handle<EcalRecHitCollection> endcapHitHandle;
  InputTag eehits("ecalRecHit","EcalRecHitsEE");
  event.getByLabel(eehits, endcapHitHandle);
  EcalRecHitCollection endcapRecHits;
  if (!endcapHitHandle.isValid()) {
//    edm::LogError("PhotonProducer") << "Error! Can't get the product "<<endcapEcalHits_.label();
    validEcalRecHits=false; 
  }
  if( validEcalRecHits) endcapRecHits = *(endcapHitHandle.product());

  if (cacheIDTopo!=eventSetup.get<CaloTopologyRecord>().cacheIdentifier()){
    cacheIDTopo=eventSetup.get<CaloTopologyRecord>().cacheIdentifier();
    eventSetup.get<CaloTopologyRecord>().get(caloTopo);
  }
  const CaloTopology * topology = caloTopo.product() ;
  const reco::CaloCluster & seedCluster = *(electron.superCluster()->seed()) ;
  // temporary, till CaloCluster->seed() is made available
  DetId seedXtalId = seedCluster.hitsAndFractions()[0].first ;
  int detector = seedXtalId.subdetId() ;
  float e3x3=-999.;
  if (detector==EcalBarrel) {
    e3x3    =   EcalClusterTools::e3x3(*(electron.superCluster()->seed()), &(barrelRecHits), &(*topology)); 
   } else {
    e3x3    =   EcalClusterTools::e3x3(*(electron.superCluster()->seed()), &(endcapRecHits), &(*topology)); 
   }
  float r9 = e3x3/electron.superCluster()->rawEnergy();

  // apply ECAL calibration scale and smearing factors depending on period and categories
  computeNewEnergy(electron, r9, event.run()) ;
  electron.correctEcalEnergy(newEnergy_,newEnergyError_) ;
  
  // apply E-p combination
  //std::cout << "[ElectronEnergCorrector] old comb momentum " << electron.p4(reco::GsfElectron::P4_COMBINATION).t() << std::endl;
  //std::cout << "[ElectronEnergCorrector] old comb momentum error " << electron.p4Error(reco::GsfElectron::P4_COMBINATION) << std::endl;
  computeEpCombination(electron) ;
  //electron.correctMomentum(newMomentum_,errorTrackMomentum_,finalMomentumError_);
  //std::cout << "[ElectronEnergCorrector] old comb momentum " << electron.p4().t() << std::endl;
  //electron.setP4(newMomentum_) ;
  //electron.setP4Error(finalMomentumError_); // this method does not exist
  electron.correctMomentum(newMomentum_,errorTrackMomentum_,finalMomentumError_);
  //std::cout << "[ElectronEnergCorrector] new comb momentum " << electron.p4(reco::GsfElectron::P4_COMBINATION).t() << std::endl;
  //std::cout << "[ElectronEnergCorrector] new comb momentum error " << electron.p4Error(reco::GsfElectron::P4_COMBINATION) << std::endl;
 }

void ElectronEnergyCalibrator::computeNewEnergy
 ( const reco::GsfElectron & electron, float r9, int run)
 {
  double scEnergy = electron.superCluster()->energy() ;
  float corr=0.;
  float dsigMC=0., corrMC=0.;
  newEnergyError_ = electron.ecalEnergyError() ;

  // Compute correction depending on run, categories and dataset
  // Corrections for the PromptReco from R. Paramattti et al.
  //   https://indico.cern.ch/getFile.py/access?contribId=7&sessionId=1&resId=0&materialId=slides&confId=155805 (Oct03, PromptV6, 05Aug, 05Jul)
  //   https://indico.cern.ch/getFile.py/access?contribId=2&resId=0&materialId=slides&confId=149567 (PromptV5)
  //   https://indico.cern.ch/getFile.py/access?contribId=2&resId=0&materialId=slides&confId=149567 (05Jul)
  //   https://hypernews.cern.ch/HyperNews/CMS/get/AUX/2011/07/06/16:50:04-57776-ScaleAndResolution_20110706.pdf (May10+PromptV4)
  // Correction for the ReReco from R. paramatti et al. (private communication, AN in preparation)
  // Corrections for PromptReco are run and R9 dependant, corrections for the ReReco are categories or EB+/EB-/EE+/EE- dependant
  // Correction for MC is a gaussian smearing for the resolution, averaged from the results over the three periods
   edm::Service<edm::RandomNumberGenerator> rng;
   if ( ! rng.isAvailable()) {
     throw cms::Exception("Configuration")
       << "XXXXXXX requires the RandomNumberGeneratorService\n"
          "which is not present in the configuration file.  You must add the service\n"
          "in the configuration file or remove the modules that require it.";
   }
  
  // data corrections 
  if (!isMC_) {
    // corrections for prompt
    if (dataset_=="Prompt") {
      if (run>=160431 && run<=167784) {
	if (electron.isEB()) {
	  if (run>=160431 && run<=163869) {
            if (r9>=0.94) corr = +0.0047;
            if (r9<0.94) corr = -0.0025;
	  } else if (run>=165071 && run<=165970) {
            if (r9>=0.94) corr = +0.0007;
            if (r9<0.94) corr = -0.0049;
	  } else if (run>=165971 && run<=166502) {
            if (r9>=0.94) corr = -0.0003;
            if (r9<0.94) corr = -0.0067;
	  } else if (run>=166503 && run<=166861) {
            if (r9>=0.94) corr = -0.0011;
            if (r9<0.94) corr = -0.0063;
	  } else if (run>=166862 && run<=167784) {
            if (r9>=0.94) corr = -0.0014;
            if (r9<0.94) corr = -0.0074;
	  } 
	} else if (electron.isEE()) {
	  if (run>=160431 && run<=163869) {
            if (r9>=0.94) corr = -0.0058;
            if (r9<0.94) corr = +0.0010;
	  } else if (run>=165071 && run<=165970) {
            if (r9>=0.94) corr = -0.0249;
            if (r9<0.94) corr = -0.0062;
	  } else if (run>=165971 && run<=166502) {
            if (r9>=0.94) corr = -0.0376;
            if (r9<0.94) corr = -0.0133;
	  } else if (run>=166503 && run<=166861) {
            if (r9>=0.94) corr = -0.0450;
            if (r9<0.94) corr = -0.0178;
	  } else if (run>=166862 && run<=167784) {
            if (r9>=0.94) corr = -0.0561;
            if (r9<0.94) corr = -0.0273;
	  } 
	}    
      } else if (run>=1700053 && run <=172619) {
	if (electron.isEB()) {
	  if (r9>=0.94) corr = -0.0011;
	  if (r9<0.94) corr = -0.0067;
	} else if (electron.isEE()) {
	  if (r9>=0.94) corr = +0.0009;
	  if (r9<0.94) corr = -0.0046;
	}  
      } else if (run>=172620 && run <=175770) {
	if (electron.isEB()) {
	  if (r9>=0.94) corr = -0.0046;
	  if (r9<0.94) corr = -0.0104;
	} else if (electron.isEE()) {
	  if (r9>=0.94) corr = +0.0337;
	  if (r9<0.94) corr = +0.0250;
        }  
      } else if (run>=175860 && run<=177139) {                      // prompt-v1 corrections for 2011B [ 175860 - 177139 ]
        if (electron.isEB() && fabs(electron.superCluster()->eta())>=1 and r9<0.94) corr = -0.0228;
        if (electron.isEB() && fabs(electron.superCluster()->eta())>=1 and r9>=0.94) corr = -0.0118;
        if (electron.isEB() && fabs(electron.superCluster()->eta())<1 and r9<0.94) corr = -0.0075;
        if (electron.isEB() && fabs(electron.superCluster()->eta())<1 and r9>=0.94) corr = -0.0034;
        if (electron.isEE() && fabs(electron.superCluster()->eta())>=2 and r9<0.94) corr = -0.0041;
        if (electron.isEE() && fabs(electron.superCluster()->eta())>=2 and r9>=0.94) corr = +0.0019;
        if (electron.isEE() && fabs(electron.superCluster()->eta())<2 and r9<0.94) corr = +0.0147;
        if (electron.isEE() && fabs(electron.superCluster()->eta())<2 and r9>=0.94) corr = +0.0168;
      } else if (run>=177140 && run<=178421) {                      // prompt-v1 corrections for 2011B [ 177140 - 178421 ]
        if (electron.isEB() && fabs(electron.superCluster()->eta())>=1 and r9<0.94) corr = -0.0239;
        if (electron.isEB() && fabs(electron.superCluster()->eta())>=1 and r9>=0.94) corr = -0.0129;
        if (electron.isEB() && fabs(electron.superCluster()->eta())<1 and r9<0.94) corr = -0.0079;
        if (electron.isEB() && fabs(electron.superCluster()->eta())<1 and r9>=0.94) corr = -0.0038;
        if (electron.isEE() && fabs(electron.superCluster()->eta())>=2 and r9<0.94) corr = -0.0011;
        if (electron.isEE() && fabs(electron.superCluster()->eta())>=2 and r9>=0.94) corr = +0.0049;
        if (electron.isEE() && fabs(electron.superCluster()->eta())<2 and r9<0.94) corr = +0.0236;
        if (electron.isEE() && fabs(electron.superCluster()->eta())<2 and r9>=0.94) corr = +0.0257;
      } else if (run>=178424 && run<=180252) {                      // prompt-v1 corrections for 2011B [ 178424 - 180252 ]
        if (electron.isEB() && fabs(electron.superCluster()->eta())>=1 and r9<0.94) corr = -0.0260;
        if (electron.isEB() && fabs(electron.superCluster()->eta())>=1 and r9>=0.94) corr = -0.0150;
        if (electron.isEB() && fabs(electron.superCluster()->eta())<1 and r9<0.94) corr = -0.0094;
        if (electron.isEB() && fabs(electron.superCluster()->eta())<1 and r9>=0.94) corr = -0.0052;
        if (electron.isEE() && fabs(electron.superCluster()->eta())>=2 and r9<0.94) corr = -0.0050;
        if (electron.isEE() && fabs(electron.superCluster()->eta())>=2 and r9>=0.94) corr = +0.0009;
        if (electron.isEE() && fabs(electron.superCluster()->eta())<2 and r9<0.94) corr = +0.0331;
        if (electron.isEE() && fabs(electron.superCluster()->eta())<2 and r9>=0.94) corr = +0.0353;
      } 
    // corrections for rereco  
    } else if (dataset_=="ReReco") {                     // corrections for ReReco
      // values from https://indico.cern.ch/conferenceDisplay.py?confId=146386
      if (run>=160329 && run <=168437) {                 // Jul05 period 160329-168437
        if (electron.isEB() && fabs(electron.superCluster()->eta())>=1 and r9<0.94) corr = -0.0150;
        if (electron.isEB() && fabs(electron.superCluster()->eta())>=1 and r9>=0.94) corr = -0.0039;
        if (electron.isEB() && fabs(electron.superCluster()->eta())<1 and r9<0.94) corr = -0.0014;
        if (electron.isEB() && fabs(electron.superCluster()->eta())<1 and r9>=0.94) corr = +0.0028;
        if (electron.isEE() && fabs(electron.superCluster()->eta())>=2 and r9<0.94) corr = -0.0050;
        if (electron.isEE() && fabs(electron.superCluster()->eta())>=2 and r9>=0.94) corr = +0.0010;
        if (electron.isEE() && fabs(electron.superCluster()->eta())<2 and r9<0.94) corr = -0.0025;
        if (electron.isEE() && fabs(electron.superCluster()->eta())<2 and r9>=0.94) corr = -0.0005;
      } else if (run>=170053 && run <=172619) {          // Aug05 period 170053-172619
        if (electron.isEB() && fabs(electron.superCluster()->eta())>=1 and r9<0.94) corr = -0.0191;
        if (electron.isEB() && fabs(electron.superCluster()->eta())>=1 and r9>=0.94) corr = -0.0081;
        if (electron.isEB() && fabs(electron.superCluster()->eta())<1 and r9<0.94) corr = -0.0030;
        if (electron.isEB() && fabs(electron.superCluster()->eta())<1 and r9>=0.94) corr = +0.0012;
        if (electron.isEE() && fabs(electron.superCluster()->eta())>=2 and r9<0.94) corr = +0.0052;
        if (electron.isEE() && fabs(electron.superCluster()->eta())>=2 and r9>=0.94) corr = +0.0113;
        if (electron.isEE() && fabs(electron.superCluster()->eta())<2 and r9<0.94) corr = +0.0041;
        if (electron.isEE() && fabs(electron.superCluster()->eta())<2 and r9>=0.94) corr = +0.0062;
      } else if (run>=172620 && run <=175770) {          // Oct03 period
        if (electron.isEB() && fabs(electron.superCluster()->eta())>=1 and r9<0.94) corr = -0.0150;
        if (electron.isEB() && fabs(electron.superCluster()->eta())>=1 and r9>=0.94) corr = -0.0039;
        if (electron.isEB() && fabs(electron.superCluster()->eta())<1 and r9<0.94) corr = +0.0001;
        if (electron.isEB() && fabs(electron.superCluster()->eta())<1 and r9>=0.94) corr = +0.0043;
        if (electron.isEE() && fabs(electron.superCluster()->eta())>=2 and r9<0.94) corr = +0.0001;
        if (electron.isEE() && fabs(electron.superCluster()->eta())>=2 and r9>=0.94) corr = +0.0062;
        if (electron.isEE() && fabs(electron.superCluster()->eta())<2 and r9<0.94) corr = +0.0026;
        if (electron.isEE() && fabs(electron.superCluster()->eta())<2 and r9>=0.94) corr = +0.0048;
      } else if (run>=175860 && run<=177139) {                      // prompt-v1 corrections for 2011B [ 175860 - 177139 ]
        if (electron.isEB() && fabs(electron.superCluster()->eta())>=1 and r9<0.94) corr = -0.0228;
        if (electron.isEB() && fabs(electron.superCluster()->eta())>=1 and r9>=0.94) corr = -0.0118;
        if (electron.isEB() && fabs(electron.superCluster()->eta())<1 and r9<0.94) corr = -0.0075;
        if (electron.isEB() && fabs(electron.superCluster()->eta())<1 and r9>=0.94) corr = -0.0034;
        if (electron.isEE() && fabs(electron.superCluster()->eta())>=2 and r9<0.94) corr = -0.0041;
        if (electron.isEE() && fabs(electron.superCluster()->eta())>=2 and r9>=0.94) corr = +0.0019;
        if (electron.isEE() && fabs(electron.superCluster()->eta())<2 and r9<0.94) corr = +0.0147;
        if (electron.isEE() && fabs(electron.superCluster()->eta())<2 and r9>=0.94) corr = +0.0168;
      } else if (run>=177140 && run<=178421) {                      // prompt-v1 corrections for 2011B [ 177140 - 178421 ]
        if (electron.isEB() && fabs(electron.superCluster()->eta())>=1 and r9<0.94) corr = -0.0239;
        if (electron.isEB() && fabs(electron.superCluster()->eta())>=1 and r9>=0.94) corr = -0.0129;
        if (electron.isEB() && fabs(electron.superCluster()->eta())<1 and r9<0.94) corr = -0.0079;
        if (electron.isEB() && fabs(electron.superCluster()->eta())<1 and r9>=0.94) corr = -0.0038;
        if (electron.isEE() && fabs(electron.superCluster()->eta())>=2 and r9<0.94) corr = -0.0011;
        if (electron.isEE() && fabs(electron.superCluster()->eta())>=2 and r9>=0.94) corr = +0.0049;
        if (electron.isEE() && fabs(electron.superCluster()->eta())<2 and r9<0.94) corr = +0.0236;
        if (electron.isEE() && fabs(electron.superCluster()->eta())<2 and r9>=0.94) corr = +0.0257;
      } else if (run>=178424 && run<=180252) {                      // prompt-v1 corrections for 2011B [ 178424 - 180252 ]
        if (electron.isEB() && fabs(electron.superCluster()->eta())>=1 and r9<0.94) corr = -0.0260;
        if (electron.isEB() && fabs(electron.superCluster()->eta())>=1 and r9>=0.94) corr = -0.0150;
        if (electron.isEB() && fabs(electron.superCluster()->eta())<1 and r9<0.94) corr = -0.0094;
        if (electron.isEB() && fabs(electron.superCluster()->eta())<1 and r9>=0.94) corr = -0.0052;
        if (electron.isEE() && fabs(electron.superCluster()->eta())>=2 and r9<0.94) corr = -0.0050;
        if (electron.isEE() && fabs(electron.superCluster()->eta())>=2 and r9>=0.94) corr = +0.0009;
        if (electron.isEE() && fabs(electron.superCluster()->eta())<2 and r9<0.94) corr = +0.0331;
        if (electron.isEE() && fabs(electron.superCluster()->eta())<2 and r9>=0.94) corr = +0.0353;
      } 
    // corrections for januray 16 rereco  
    } else if (dataset_=="Jan16ReReco") {                     // corrections for january 16 ReReco
      // values from http://indico.cern.ch/getFile.py/access?contribId=2&resId=0&materialId=slides&confId=176520
      if (electron.isEB() && fabs(electron.superCluster()->eta())<1 and r9<0.94) {
        if (run>=160431 && run<=167913) corr = -0.0014; 
	if (run>=170000 && run<=172619) corr = -0.0016;   
	if (run>=172620 && run<=173692) corr = -0.0017;  
	if (run>=175830 && run<=177139) corr = -0.0021;
	if (run>=177140 && run<=178421) corr = -0.0025;
	if (run>=178424 && run<=180252) corr = -0.0024;
      } else if (electron.isEB() && fabs(electron.superCluster()->eta())<1 and r9>=0.94) {      
        if (run>=160431 && run<=167913) corr = 0.0059; 
	if (run>=170000 && run<=172619) corr = 0.0046;   
	if (run>=172620 && run<=173692) corr = 0.0045;  
	if (run>=175830 && run<=177139) corr = 0.0042;
	if (run>=177140 && run<=178421) corr = 0.0038;
	if (run>=178424 && run<=180252) corr = 0.0039;
      } else if (electron.isEB() && fabs(electron.superCluster()->eta())>=1 and r9<0.94) {      
        if (run>=160431 && run<=167913) corr = -0.0045; 
	if (run>=170000 && run<=172619) corr = -0.0066;   
	if (run>=172620 && run<=173692) corr = -0.0058;  
	if (run>=175830 && run<=177139) corr = -0.0073;
	if (run>=177140 && run<=178421) corr = -0.0075;
	if (run>=178424 && run<=180252) corr = -0.0071;
      } else if (electron.isEB() && fabs(electron.superCluster()->eta())>=1 and r9>=0.94) {      
        if (run>=160431 && run<=167913) corr = 0.0084; 
	if (run>=170000 && run<=172619) corr = 0.0063;   
	if (run>=172620 && run<=173692) corr = 0.0071;  
	if (run>=175830 && run<=177139) corr = 0.0056;
	if (run>=177140 && run<=178421) corr = 0.0054;
	if (run>=178424 && run<=180252) corr = 0.0058;
      } else if (electron.isEE() && fabs(electron.superCluster()->eta())<2 and r9<0.94) {
        if (run>=160431 && run<=167913) corr = -0.0082; 
	if (run>=170000 && run<=172619) corr = -0.0025;   
	if (run>=172620 && run<=173692) corr = -0.0035;  
	if (run>=175830 && run<=177139) corr = -0.0017;
	if (run>=177140 && run<=178421) corr = -0.0010;
	if (run>=178424 && run<=180252) corr = 0.0030;
      } else if (electron.isEE() && fabs(electron.superCluster()->eta())<2 and r9>=0.94) {      
        if (run>=160431 && run<=167913) corr = -0.0033; 
	if (run>=170000 && run<=172619) corr = 0.0024;   
	if (run>=172620 && run<=173692) corr = 0.0014;  
	if (run>=175830 && run<=177139) corr = 0.0032;
	if (run>=177140 && run<=178421) corr = 0.0040;
	if (run>=178424 && run<=180252) corr = 0.0079;
      } else if (electron.isEE() && fabs(electron.superCluster()->eta())>=2 and r9<0.94) {      
        if (run>=160431 && run<=167913) corr = -0.0064; 
	if (run>=170000 && run<=172619) corr = -0.0046;   
	if (run>=172620 && run<=173692) corr = -0.0029;  
	if (run>=175830 && run<=177139) corr = -0.0040;
	if (run>=177140 && run<=178421) corr = -0.0050;
	if (run>=178424 && run<=180252) corr = -0.0059;
      } else if (electron.isEE() && fabs(electron.superCluster()->eta())>=2 and r9>=0.94) {      
        if (run>=160431 && run<=167913) corr = 0.0042; 
	if (run>=170000 && run<=172619) corr = 0.0060;   
	if (run>=172620 && run<=173692) corr = 0.0077;  
	if (run>=175830 && run<=177139) corr = 0.0067;
	if (run>=177140 && run<=178421) corr = 0.0056;
	if (run>=178424 && run<=180252) corr = 0.0047;
      } 
    }
  } 
    // MC smearing dsig is needed also for data for theenergy error, take it from the last MC values consistant
    // with the data choice
//  else  { // MC corrections
    if (dataset_=="Summer11"||dataset_=="ReReco") { // values from https://indico.cern.ch/conferenceDisplay.py?confId=146386
      if (electron.isEB() && fabs(electron.superCluster()->eta())<1 && r9<0.94) dsigMC = 0.01;
      if (electron.isEB() && fabs(electron.superCluster()->eta())<1 && r9>=0.94) dsigMC = 0.0099;
      if (electron.isEB() && fabs(electron.superCluster()->eta())>=1 && r9<0.94) dsigMC = 0.0217;
      if (electron.isEB() && fabs(electron.superCluster()->eta())>=1 && r9>=0.94) dsigMC = 0.0157;
      if (electron.isEE() && fabs(electron.superCluster()->eta())<2 && r9<0.94) dsigMC = 0.0326;
      if (electron.isEE() && fabs(electron.superCluster()->eta())<2 && r9>=0.94) dsigMC = 0.0330;
      if (electron.isEE() && fabs(electron.superCluster()->eta())>=2 && r9<0.94) dsigMC = 0.0331;
      if (electron.isEE() && fabs(electron.superCluster()->eta())>=2 && r9>=0.94) dsigMC = 0.0378;
    } else if (dataset_=="Fall11"||dataset_=="Jan16ReReco") { // values from https://hypernews.cern.ch/HyperNews/CMS/get/higgs2g/634.html, consistant with Jan16ReReco corrections
      if (electron.isEB() && fabs(electron.superCluster()->eta())<1 && r9<0.94) dsigMC = 0.0096;
      if (electron.isEB() && fabs(electron.superCluster()->eta())<1 && r9>=0.94) dsigMC = 0.0074;
      if (electron.isEB() && fabs(electron.superCluster()->eta())>=1 && r9<0.94) dsigMC = 0.0196;
      if (electron.isEB() && fabs(electron.superCluster()->eta())>=1 && r9>=0.94) dsigMC = 0.0141;
      if (electron.isEE() && fabs(electron.superCluster()->eta())<2 && r9<0.94) dsigMC = 0.0279;
      if (electron.isEE() && fabs(electron.superCluster()->eta())<2 && r9>=0.94) dsigMC = 0.0268;
      if (electron.isEE() && fabs(electron.superCluster()->eta())>=2 && r9<0.94) dsigMC = 0.0301;
      if (electron.isEE() && fabs(electron.superCluster()->eta())>=2 && r9>=0.94) dsigMC = 0.0293;   
    }
//  }
  
  // now correct the energy
  // correction for data
  if (!isMC_) newEnergy_ = scEnergy/(1+corr);
  // smearing for MC
  if (isMC_) {
    CLHEP::RandGaussQ gaussDistribution(rng->getEngine(), 1.,dsigMC);
    corrMC = gaussDistribution.fire();
    newEnergy_ = scEnergy*corrMC;  
  }  
  // correct energy error for MC and for data as error is obtained from (ideal) MC parametrisation
  if (updateEnergyError_)
   newEnergyError_ = sqrt(newEnergyError_*newEnergyError_ + dsigMC*dsigMC*newEnergy_*newEnergy_) ;
  std::cout << "[ElectronEnergyCalibrator] SC corrected energy " << electron.superCluster()->energy() << " new corrected energy " << newEnergy_ << std::endl;
  std::cout << "[ElectronEnergyCalibrator] SC corrected energy error " << electron.ecalEnergyError() << " new corrected energy error " << newEnergyError_ << std::endl;

 }


void ElectronEnergyCalibrator::computeEpCombination
 ( reco::GsfElectron & electron )
 {

  float scEnergy = electron.ecalEnergy() ;
  int elClass = electron.classification() ;

  float trackMomentum  = electron.trackMomentumAtVtx().R() ;
  errorTrackMomentum_ = 999. ;
  
  // retreive momentum error 
  //MultiGaussianState1D qpState(MultiGaussianStateTransform::multiState1D(vtxTsos,0));
  //GaussianSumUtilities1D qpUtils(qpState);
  errorTrackMomentum_ = electron.trackMomentumError();

  float finalMomentum = electron.p4().t(); // initial
  float finalMomentumError = 999.;
  
  // first check for large errors
 
  if (errorTrackMomentum_/trackMomentum > 0.5 && electron.ecalEnergyError()/scEnergy <= 0.5) {
    finalMomentum = scEnergy;    finalMomentumError = electron.ecalEnergyError();
   }
  else if (errorTrackMomentum_/trackMomentum <= 0.5 && electron.ecalEnergyError()/scEnergy > 0.5){
    finalMomentum = trackMomentum;  finalMomentumError = errorTrackMomentum_;
   }
  else if (errorTrackMomentum_/trackMomentum > 0.5 && electron.ecalEnergyError()/scEnergy > 0.5){
    if (errorTrackMomentum_/trackMomentum < electron.ecalEnergyError()/scEnergy) {
      finalMomentum = trackMomentum; finalMomentumError = errorTrackMomentum_;
     }
    else{
      finalMomentum = scEnergy; finalMomentumError = electron.ecalEnergyError();
     }
  }
  
  // then apply the combination algorithm
  else {

     // calculate E/p and corresponding error
    float eOverP = scEnergy / trackMomentum;
    float errorEOverP = sqrt(
			     (electron.ecalEnergyError()/trackMomentum)*(electron.ecalEnergyError()/trackMomentum) +
			     (scEnergy*errorTrackMomentum_/trackMomentum/trackMomentum)*
			     (scEnergy*errorTrackMomentum_/trackMomentum/trackMomentum));
    
    if ( eOverP  > 1 + 2.5*errorEOverP )
      {
	finalMomentum = scEnergy; finalMomentumError = electron.ecalEnergyError();
	if ((elClass==reco::GsfElectron::GOLDEN) && electron.isEB() && (eOverP<1.15))
	  {
	    if (scEnergy<15) {finalMomentum = trackMomentum ; finalMomentumError = errorTrackMomentum_;}
	  }
      }
    else if ( eOverP < 1 - 2.5*errorEOverP )
      {
	finalMomentum = scEnergy; finalMomentumError = electron.ecalEnergyError();
	if (elClass==reco::GsfElectron::SHOWERING)
	  {
	    if (electron.isEB())
	      {
		if(scEnergy<18) {finalMomentum = trackMomentum; finalMomentumError = errorTrackMomentum_;}
	      }
	    else if (electron.isEE())
	      {
		if(scEnergy<13) {finalMomentum = trackMomentum; finalMomentumError = errorTrackMomentum_;}
	      }
	    else
	      { edm::LogWarning("ElectronMomentumCorrector::correct")<<"nor barrel neither endcap electron ?!" ; }
	  }
	else if (electron.isGap())
	  {
	    if(scEnergy<60) {finalMomentum = trackMomentum; finalMomentumError = errorTrackMomentum_;}
	  }
      }
    else 
      {
	// combination
	finalMomentum = (scEnergy/electron.ecalEnergyError()/electron.ecalEnergyError() + trackMomentum/errorTrackMomentum_/errorTrackMomentum_) /
	  (1/electron.ecalEnergyError()/electron.ecalEnergyError() + 1/errorTrackMomentum_/errorTrackMomentum_);
	float finalMomentumVariance = 1 / (1/electron.ecalEnergyError()/electron.ecalEnergyError() + 1/errorTrackMomentum_/errorTrackMomentum_);
	finalMomentumError = sqrt(finalMomentumVariance);
      }
    
  }
  
  math::XYZTLorentzVector oldMomentum = electron.p4() ;
  newMomentum_ = math::XYZTLorentzVector
   ( oldMomentum.x()*finalMomentum/oldMomentum.t(),
     oldMomentum.y()*finalMomentum/oldMomentum.t(),
     oldMomentum.z()*finalMomentum/oldMomentum.t(),
     finalMomentum ) ;
  finalMomentumError_ =  finalMomentumError;  
  //std::cout << "[ElectronEnergCorrector] old comb momentum " << oldMomentum.t() << " new comb momentum " << newMomentum_.t() << std::endl;

 }
