import FWCore.ParameterSet.Config as cms
import os


from FWCore.ParameterSet.VarParsing import VarParsing

#input cmsRun options
options = VarParsing ('analysis')
ptions.inputFiles = '/store/mc/RunIISpring15DR74/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/Asympt50nsRaw_MCRUN2_74_V9A-v3/70000/44B5E6C2-3B08-E511-A88B-0025907253D2.root'
options.outputFile = "MiniAOD_FR_QCD.root"
options.parseArguments()

#name the process
process = cms.Process("TreeProducerFromMiniAOD")

#Make the framework shutup
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100;
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')

#50 ns global tag for MC replace with 'GR_P_V56' for prompt reco. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Prompt_reconstruction_Global_Tag 
process.GlobalTag.globaltag = 'MCRUN2_74_V9A'

#how many events to run over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
#input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
)

#whether running MC or data
options.register(
    'isMC',
    1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    'Set to 1 for simulated samples - updates GT, emulates HCAL TPGs.'
)

##################################################
# Main
process.byLooseCombinedIsolationDeltaBetaCorr3Hits = cms.EDAnalyzer("MiniAODfakeRate_alt",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("byLooseCombinedIsolationDeltaBetaCorr3Hits") 
)
process.byMediumCombinedIsolationDeltaBetaCorr3Hits = cms.EDAnalyzer("MiniAODfakeRate_alt",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("byMediumCombinedIsolationDeltaBetaCorr3Hits")
)
process.byTightCombinedIsolationDeltaBetaCorr3Hits = cms.EDAnalyzer("MiniAODfakeRate_alt",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("byTightCombinedIsolationDeltaBetaCorr3Hits")
)
#It tells me that byCombinedIsolationDeltaBetaCorrRaw3Hits is not in the miniAOD
#process.byCombinedIsolationDeltaBetaCorr3Hits = cms.EDAnalyzer("MiniAODeffi",
 #   vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
  #  taus = cms.InputTag("slimmedTaus"),
   # jets = cms.InputTag("slimmedJets"),
    #tauID = cms.string("byCombinedIsolationDeltaBetaCorrRaw3Hits")
#)
#process.ChargedIsoPtSum = cms.EDAnalyzer("MiniAODeffi",
 #   vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
  #  taus = cms.InputTag("slimmedTaus"),
   # jets = cms.InputTag("slimmedJets"),
    #tauID = cms.string("chargedIsoPtSum")
#)
process.neutralIsoPtSum= cms.EDAnalyzer("MiniAODfakeRate_alt",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("neutralIsoPtSum")
)
process.puCorrPtSum= cms.EDAnalyzer("MiniAODfakeRate_alt",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("puCorrPtSum")
)
process.againstMuonLoose3 = cms.EDAnalyzer("MiniAODfakeRate_alt",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("againstMuonLoose3")
)
process.againstMuonTight3 = cms.EDAnalyzer("MiniAODfakeRate_alt",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("againstMuonTight3")
)
process.againstElectronVLooseMVA5 = cms.EDAnalyzer("MiniAODfakeRate_alt",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("againstElectronVLooseMVA5")
)
process.againstElectronLooseMVA5 = cms.EDAnalyzer("MiniAODfakeRate_alt",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("againstElectronLooseMVA5")
)
process.againstElectronMediumMVA5 = cms.EDAnalyzer("MiniAODfakeRate_alt",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    tauID = cms.string("againstElectronMediumMVA5")
)

###################################################
#Global sequence

process.p = cms.Path(
                     process.byLooseCombinedIsolationDeltaBetaCorr3Hits*
		     process.byMediumCombinedIsolationDeltaBetaCorr3Hits*
		     process.byTightCombinedIsolationDeltaBetaCorr3Hits*
 	             #process.byCombinedIsolationDeltaBetaCorrRaw3Hits*
		     #process.chargedIsoPtSum*
		     process.neutralIsoPtSum*
	 	     process.puCorrPtSum*
		     process.againstMuonLoose3*
	 	     process.againstMuonTight3*
		     process.againstElectronVLooseMVA5*
		     process.againstElectronLooseMVA5*
		     process.againstElectronMediumMVA5
                     )

#output file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
)

#print out all processes used when running- useful check to see if module ran
#UNCOMMENT BELOW
#dump_file = open('dump.py','w')
#dump_file.write(process.dumpPython())
