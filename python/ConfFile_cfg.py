import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.connect = cms.string('frontier://FrontierProd/CMS_COND_31X_GLOBALTAG')
process.GlobalTag.globaltag = cms.string('POSTLS162_V2::All')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    outputCommands = cms.untracked.vstring(
                                           'keep *_*_*'),
    fileNames = cms.untracked.vstring(
        "file:/hdfs/store/user/alevine/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6_Sept12AtUW_Fall13dr-tsg_PU40bx25_POSTLS162_V2-v1/05a81b8d696d27a5c3c2ca036967addd/skim_100_1_35f.root"
    )
)

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff") #loading the configuration
from PhysicsTools.PatAlgos.tools.tauTools import * # to get the most up-to-date discriminators when using patTaus
switchToPFTauHPS(process) # this line and the one above can be ignored when running on reco::PFTauCollection directly

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('tauTree.root')
)

process.tauAnalyzer = cms.EDAnalyzer('tauAnalyzer',
                      recoTau              = cms.InputTag("hpsPFTauProducer"),
                      recoTauDiscriminator = cms.InputTag("hpsPFTauDiscriminationByLooseIsolation")

)

process.p = cms.Path(process.PFTau*process.tauAnalyzer)
