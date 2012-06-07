import FWCore.ParameterSet.Config as cms
import os

from TrackingTools.Configuration.TrackingTools_cff import *

process = cms.Process("NewEScale")

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.EventContent.EventContent_cff")

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedGsfElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(1),
        engineName = cms.untracked.string('TRandom3')
    ),
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring
    (
#    'file:/data_CMS/cms/charlot/Run2011/AllCandidatesEPS11/HZZCandidates.root'
#    '/store/data/Summer11/DYToEE_M-800_TuneZ2_7TeV-pythia6-tauola/AODSIM/PU_S3_START42_V11-v2/0000/0ABF7CD0-8888-E011-8561-1CC1DE051038.root'    
    '/store/data/Run2012A/DoubleElectron/AOD/PromptReco-v1/000/190/645/4812E65D-8D82-E111-B659-003048F024DC.root'
    ),
#    eventsToProcess = cms.untracked.VEventRange('173243:16706390')   
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
#    fileName = cms.untracked.string('CandidateZ_newEscale.root')
#    fileName = cms.untracked.string('testMC.root')
)

process.load("EgammaCalibratedGsfElectrons.CalibratedElectronProducers.calibratedGsfElectrons_cfi")

# dataset to correct
#process.calibratedGsfElectrons.inputDataset = cms.string("Jan16ReReco")
process.calibratedGsfElectrons.inputDataset = cms.string("Prompt2012")
#process.calibratedGsfElectrons.inputDataset = cms.string("Summer12")
#process.calibratedGsfElectrons.inputDataset = cms.string("Summer11")
#process.calibratedGsfElectrons.inputDataset = cms.string("Fall11")
#process.calibratedGsfElectrons.inputDataset = cms.string("Fall11")
process.calibratedGsfElectrons.isMC = cms.bool(False)
process.calibratedGsfElectrons.isAOD = cms.bool(True)
process.calibratedGsfElectrons.updateEnergyError = cms.bool(True)
#process.calibratedGsfElectrons.debug = cms.bool(True)


process.p = cms.Path(process.calibratedGsfElectrons)

#process.outpath = cms.EndPath(process.out)
process.GlobalTag.globaltag = 'GR_P_V32::All'
#process.GlobalTag.globaltag = 'START52_V9::All'




