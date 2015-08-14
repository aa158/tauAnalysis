import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

#input cmsRun options
options = VarParsing ('analysis')
options.inputFiles = '/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/244/00000/68275270-7C27-E511-B1F0-02163E011A46.root'

#'/store/mc/RunIISpring15DR74/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/00000/2C0279E7-9B16-E511-B5FF-90B11C2ACF20.root'

#'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/162/00000/160C08A3-4227-E511-B829-02163E01259F.root' 

#'/store/data/Run2015B/SingleMuon/MINIAOD/PromptReco-v1/000/251/162/00000/160C08A3-4227-E511-B829-02163E01259F.root'

#'/store/mc/RunIISpring15DR74/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/00000/2C0279E7-9B16-E511-B5FF-90B11C2ACF20.root'

#'/store/mc/RunIISpring15DR74/ZprimeToTauTau_M_500_TuneCUETP8M1_tauola_13TeV_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/60000/CAE0DF38-7D19-E511-9DA4-0025904A8EC8.root'

#options.inputFiles = '/store/mc/RunIISpring15DR74/SUSYGluGluToHToTauTau_M-250_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/80000/6EBE4F55-4E03-E511-9FE6-0CC47A13D416.root'
options.outputFile = "dataMCcomp.root"

#whether running MC or data
options.register(
    'isMC',
    0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    'Set to 1 for simulated samples - updates GT, emulates HCAL TPGs.'
)

options.register(
    'crossSection',
    -1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Process cross section"
)

options.parseArguments()

#name the process
process = cms.Process("TreeProducerFromMiniAOD")
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100;
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#50 ns global tag for MC replace with 'GR_P_V56' for prompt reco. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Prompt_reconstruction_Global_Tag 

#how many events to run over
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
)

if options.isMC:
    isMC = cms.bool(True)
    process.GlobalTag.globaltag = 'MCRUN2_74_V9A'
else:
    isMC = cms.bool(False)
    process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v0'
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt').getVLuminosityBlockRange()

##################################################
# Main
process.byLooseCombinedIsolationDeltaBetaCorr3Hits = cms.EDAnalyzer("dataMCcomp",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    muons = cms.InputTag('slimmedMuons'),
    tauID = cms.string("byLooseCombinedIsolationDeltaBetaCorr3Hits"),
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    isMC = isMC
)
process.byMediumCombinedIsolationDeltaBetaCorr3Hits = cms.EDAnalyzer("dataMCcomp",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    muons = cms.InputTag('slimmedMuons'),
    tauID = cms.string("byMediumCombinedIsolationDeltaBetaCorr3Hits"),
    isMC = isMC
)
process.byTightCombinedIsolationDeltaBetaCorr3Hits = cms.EDAnalyzer("dataMCcomp",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    muons = cms.InputTag('slimmedMuons'),
    tauID = cms.string("byTightCombinedIsolationDeltaBetaCorr3Hits"),
    isMC = isMC
)
#It tells me that byCombinedIsolationDeltaBetaCorrRaw3Hits is not in the miniAOD
#process.byCombinedIsolationDeltaBetaCorr3Hits = cms.EDAnalyzer("dataMCcomp",
 #   vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
  #  taus = cms.InputTag("slimmedTaus"),
   # jets = cms.InputTag("slimmedJets"),
    #tauID = cms.string("byCombinedIsolationDeltaBetaCorrRaw3Hits")
#)
#process.ChargedIsoPtSum = cms.EDAnalyzer("dataMCcomp",
 #   vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
  #  taus = cms.InputTag("slimmedTaus"),
   # jets = cms.InputTag("slimmedJets"),
    #tauID = cms.string("chargedIsoPtSum")
#)
process.neutralIsoPtSum= cms.EDAnalyzer("dataMCcomp",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    muons = cms.InputTag('slimmedMuons'),
    tauID = cms.string("neutralIsoPtSum"),
    xSec = cms.untracked.double(options.crossSection),
    isMC = isMC
)
process.puCorrPtSum= cms.EDAnalyzer("dataMCcomp",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    muons = cms.InputTag('slimmedMuons'),
    tauID = cms.string("puCorrPtSum"),
    xSec = cms.untracked.double(options.crossSection),
    isMC = isMC
)
process.againstMuonLoose3 = cms.EDAnalyzer("dataMCcomp",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    muons = cms.InputTag('slimmedMuons'),
    tauID = cms.string("againstMuonLoose3"),
    xSec = cms.untracked.double(options.crossSection),
    isMC = isMC
)
process.againstMuonTight3 = cms.EDAnalyzer("dataMCcomp",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    muons = cms.InputTag('slimmedMuons'),
    tauID = cms.string("againstMuonTight3"),
    xSec = cms.untracked.double(options.crossSection),
    isMC = isMC
)
process.againstElectronVLooseMVA5 = cms.EDAnalyzer("dataMCcomp",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    muons = cms.InputTag('slimmedMuons'),
    tauID = cms.string("againstElectronVLooseMVA5"),
    xSec = cms.untracked.double(options.crossSection),
    isMC = isMC
)
process.againstElectronLooseMVA5 = cms.EDAnalyzer("dataMCcomp",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    muons = cms.InputTag('slimmedMuons'),
    tauID = cms.string("againstElectronLooseMVA5"),
    xSec = cms.untracked.double(options.crossSection),
    isMC = isMC
)
process.againstElectronMediumMVA5 = cms.EDAnalyzer("dataMCcomp",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    muons = cms.InputTag('slimmedMuons'),
    tauID = cms.string("againstElectronMediumMVA5"),
    xSec = cms.untracked.double(options.crossSection),
    isMC = isMC
)

process.decayModeFinding = cms.EDAnalyzer("dataMCcomp",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    muons = cms.InputTag('slimmedMuons'),
    tauID = cms.string("decayModeFinding"),
    xSec = cms.untracked.double(options.crossSection),
    isMC = isMC
)

process.decayModeFindingNewDMs = cms.EDAnalyzer("dataMCcomp",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    muons = cms.InputTag('slimmedMuons'),
    tauID = cms.string("decayModeFindingNewDMs"),
    xSec = cms.untracked.double(options.crossSection),
    isMC = isMC
)

process.kOneProng0PiZero = cms.EDAnalyzer("dataMCcomp",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    muons = cms.InputTag('slimmedMuons'),
    tauID = cms.string("kOneProng0PiZero"),
    xSec = cms.untracked.double(options.crossSection),
    isMC = isMC
)

process.kOneProng1PiZero = cms.EDAnalyzer("dataMCcomp",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    muons = cms.InputTag('slimmedMuons'),
    tauID = cms.string("kOneProng1PiZero"),
    xSec = cms.untracked.double(options.crossSection),
    isMC = isMC
)

process.kOneProng2PiZero = cms.EDAnalyzer("dataMCcomp",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    muons = cms.InputTag('slimmedMuons'),
    tauID = cms.string("kOneProng2PiZero"),
    xSec = cms.untracked.double(options.crossSection),
    isMC = isMC
)
###################################################
#Global sequence

process.p = cms.Path(
                     process.byLooseCombinedIsolationDeltaBetaCorr3Hits*
		     process.byMediumCombinedIsolationDeltaBetaCorr3Hits*
		     process.byTightCombinedIsolationDeltaBetaCorr3Hits*
 		     #process.byCombinedIsolationDeltaBetaCorrRaw3Hits*
		     #process.chargedIsoPtSum*
		     #process.neutralIsoPtSum*
	 	     #process.puCorrPtSum*
		     #process.againstMuonLoose3*
	 	     #process.againstMuonTight3*
		     #process.againstElectronVLooseMVA5*
		     #process.againstElectronLooseMVA5*
		     #process.againstElectronMediumMVA5*
		     process.decayModeFinding*
		     process.decayModeFindingNewDMs
		     #process.kOneProng0PiZero*
		     #process.kOneProng1PiZero*
		     #process.kOneProng2PiZero
                     )

process.TFileService = cms.Service("TFileService",
        fileName = cms.string(options.outputFile)
)

#print out all processes used when running- useful check to see if module ran
#UNCOMMENT BELOW
#dump_file = open('dump.py','w')
#dump_file.write(process.dumpPython())
