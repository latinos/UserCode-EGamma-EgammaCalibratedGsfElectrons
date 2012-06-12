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
#    '/store/data/Run2011A/DoubleElectron/AOD/03Oct2011-v1/0000/8696915B-1AEF-E011-91ED-003048678B7C.root'    
#    '/store/data/Run2011A/DoubleElectron/AOD/PromptReco-v6/000/173/659/4645CC17-B8CD-E011-8BA0-001D09F2906A.root'
#     '/store/data/Run2011A/DoubleElectron/AOD/16Jan2012-v1/0000/DC832C9D-AA43-E111-AB14-00261894382A.root' 
     '/store/data/Run2011B/DoubleElectron/AOD/16Jan2012-v1/0000/D2B7044D-8444-E111-9172-00A0D1EE9644.root'    ),
    eventsToProcess = cms.untracked.VEventRange('180250:591651181')   
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
#    fileName = cms.untracked.string('CandidateZ_newEscale.root')
    fileName = cms.untracked.string('testMC.root')
)

process.load("EgammaCalibratedGsfElectrons.CalibratedElectronProducers.calibratedGsfElectrons_cfi")

process.calibratedGsfElectrons.inputDataset = cms.string("ReReco2012")
#process.calibratedGsfElectrons.inputDataset = cms.string("Jan16ReReco")
#process.calibratedGsfElectrons.inputDataset = cms.string("Prompt")
#process.calibratedGsfElectrons.inputDataset = cms.string("ReReco")
#process.calibratedGsfElectrons.inputDataset = cms.string("Summer11")
#process.calibratedGsfElectrons.inputDataset = cms.string("Fall11")
process.calibratedGsfElectrons.isMC = cms.bool(False)
process.calibratedGsfElectrons.isAOD = cms.bool(True)
#process.calibratedGsfElectrons.isAOD = cms.bool(False)
process.calibratedGsfElectrons.updateEnergyError = cms.bool(True)
process.calibratedGsfElectrons.debug = cms.bool(True)


process.p = cms.Path(process.calibratedGsfElectrons)

process.outpath = cms.EndPath(process.out)
process.GlobalTag.globaltag = 'GR_R_52_V9::All'
#process.GlobalTag.globaltag = 'GR_R_42_V18::All'
#process.GlobalTag.globaltag = 'START42_V11::All'




