import FWCore.ParameterSet.Config as cms
import os

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
#     '/store/data/Run2011A/DoubleElectron/AOD/16Jan2012-v1/0000/DC832C9D-AA43-E111-AB14-00261894382A.root' 
#     '/store/data/Run2011B/DoubleElectron/AOD/16Jan2012-v1/0000/D2B7044D-8444-E111-9172-00A0D1EE9644.root'   
#     '/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V5-v2/0000/FE123555-F27A-E111-8E40-003048D46046.root'
#     '/store/data/Run2012A/DoubleElectron/RECO/PromptReco-v1/000/193/575/F2385979-F399-E111-BE08-003048F1110E.root'     
#      '/store/data/Run2012B/DoubleElectron/AOD/PromptReco-v1/000/196/350/F6DE6C77-4BB8-E111-8B09-0030486780AC.root'
#      'file:/data_CMS/cms/charlot/CMSSW/CMSSW_5_2_3_Escale2012/src/EgammaCalibratedGsfElectrons/CalibratedElectronProducers/test/hzzEvents_Run2012_NEW_AOD.root'    
      'file:/data_CMS/cms/charlot/CMSSW/CMSSW_5_2_3_Escale2012/src/EgammaCalibratedGsfElectrons/CalibratedElectronProducers/test/hzzEvents_Run2012_PromptReco_Jun08JSON_AOD.fromUCSD.root'    
    ) 
    ,eventsToProcess = cms.untracked.VEventRange('194051:6362525')
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    fileName = cms.untracked.string('testData.root')
#    fileName = cms.untracked.string('testMC.root')
)

process.load("EgammaCalibratedGsfElectrons.CalibratedElectronProducers.calibratedGsfElectrons_cfi")

#process.calibratedGsfElectrons.inputDataset = cms.string("Jun08ReReco")
#process.calibratedGsfElectrons.inputDataset = cms.string("May23ReReco")
process.calibratedGsfElectrons.inputDataset = cms.string("Prompt2012")
#process.calibratedGsfElectrons.inputDataset = cms.string("Summer12")
#process.calibratedGsfElectrons.inputDataset = cms.string("Jan16ReReco")
#process.calibratedGsfElectrons.inputDataset = cms.string("Prompt")
#process.calibratedGsfElectrons.inputDataset = cms.string("ReReco")
#process.calibratedGsfElectrons.inputDataset = cms.string("Summer11")
#process.calibratedGsfElectrons.inputDataset = cms.string("Fall11")
process.calibratedGsfElectrons.isMC = cms.bool(False)
process.calibratedGsfElectrons.isAOD = cms.bool(True)
process.calibratedGsfElectrons.updateEnergyError = cms.bool(True)
#process.calibratedGsfElectrons.debug = cms.bool(True)


process.p = cms.Path(process.calibratedGsfElectrons)

#process.outpath = cms.EndPath(process.out)
#process.GlobalTag.globaltag = 'GR_R_52_V8::All'
process.GlobalTag.globaltag = 'GR_P_V37::All'
#process.GlobalTag.globaltag = 'START52_V11::All'




