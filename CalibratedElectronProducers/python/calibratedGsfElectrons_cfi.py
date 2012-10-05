
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
    isAOD = cms.bool(False),
    
    # set to True to get debugging printout   
    debug = cms.bool(False),

    #set to True to apply corrections
    applyCorrections = cms.bool(True),
    
    # input datasets
    # Prompt means May10+Promptv4+Aug05+Promptv6 for 2011
    # ReReco means Jul05+Aug05+Oct03 for 2011
    # Jan16ReReco means Jan16 for 2011
    # Summer11 means summer11 MC..
    # Summer12 means summer12 MC pre-ICHEP.
    # Summer12_DR53X_HCP2012 means summer12 MC for HCP 2012.
    # 2012Jul13ReReco

    #inputDataset = cms.string("ReReco"),
    inputDataset = cms.string("Prompt"),
    
)


