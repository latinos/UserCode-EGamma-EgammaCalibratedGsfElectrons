
import FWCore.ParameterSet.Config as cms


#==============================================================================
# corrected gsf electrons
#==============================================================================

calibratedGsfElectrons = cms.EDProducer("CalibratedGsfElectronProducer",

    # input collections
    inputGsfElectronsTag = cms.InputTag("gsfElectrons"),

    # data or MC corrections
    # if isMC is false, data corrections are applied
    isMC = cms.bool(False),

    # set to True to read AOD format
    isAOD = cms.bool(True),
    
    # set to True to get debugging printout   
    debug = cms.bool(False),
    
    # set to True to propagate MC extra smearing to the electron momentum error  
    updateEnergyError = cms.bool(True),
    
)


