import FWCore.ParameterSet.Config as cms
import os

from TrackingTools.Configuration.TrackingTools_cff import *

process = cms.Process("NewScaleData")

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.EventContent.EventContent_cff")

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",

    # Include a PSet for each module label that needs a
    # random engine.  The name is the module label.
    # You must supply a seed or seeds.
    # Optionally an engine type can be specified
    #t1 = cms.PSet(
    #    initialSeed = cms.untracked.uint32(81)
    #),
    calibratedGsfElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(1),
        engineName = cms.untracked.string('TRandom3')
    ),
    #t3 = cms.PSet(
    #    engineName = cms.untracked.string('HepJamesRandom'),
    #    initialSeed = cms.untracked.uint32(84)
    #),
    #t4 = cms.PSet(
    #    engineName = cms.untracked.string('RanecuEngine'),
    #    initialSeedSet = cms.untracked.vuint32(1, 2)
    #),
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring
    (
    'file:/data_CMS/cms/charlot/Run2011/AllCandidatesEPS11/HZZCandidates.root'
#    '/store/data/Summer11/DYToEE_M-800_TuneZ2_7TeV-pythia6-tauola/AODSIM/PU_S3_START42_V11-v2/0000/0ABF7CD0-8888-E011-8561-1CC1DE051038.root'    
    )
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    fileName = cms.untracked.string('HZZCandidates_newEscale.root')
#    fileName = cms.untracked.string('testMC.root')
)

process.load("EgammaCalibratedGsfElectrons.CalibratedElectronProducers.calibratedGsfElectrons_cfi")

# dataset to correct
process.calibratedGsfElectrons.inputDataset = cms.string("ReReco")
#process.calibratedGsfElectrons.inputDataset = cms.string("Summer11")
process.calibratedGsfElectrons.isMC = cms.bool(False)

process.p = cms.Path(process.calibratedGsfElectrons)

process.outpath = cms.EndPath(process.out)
process.GlobalTag.globaltag = 'GR_R_42_V18::All'




