import FWCore.ParameterSet.Config as cms

OUTFILE='out_L1TEmulatorTEST.root'

process = cms.Process('L1TEMULATION')

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

# Select the Message Logger output you would like to see:

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.load('L1Trigger/L1TYellow/l1t_debug_messages_cfi')
#process.load('L1Trigger/L1TYellow/l1t_info_messages_cfi')

process.load('L1Trigger/L1TCalorimeter/l1tStage1CaloParams_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-10)
    )

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    #fileNames = cms.untracked.vstring("/store/relval/CMSSW_7_0_0_pre8/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/START70_V1-v1/00000/262AA156-744A-E311-9829-002618943945.root")
    # fileNames = cms.untracked.vstring("/store/RelVal/CMSSW_7_0_0_pre4/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_ST62_V8-v1/00000/22610530-FC24-E311-AF35-003048FFD7C2.root")
    #fileNames = cms.untracked.vstring("/store/relval/CMSSW_7_0_0_pre4/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_ST62_V8-v1/00000/22610530-FC24-E311-AF35-003048FFD7C2.root")
    fileNames = cms.untracked.vstring("file:/uscmst1b_scratch/lpc1/lpctrig/apana/L1Upgrade/262AA156-744A-E311-9829-002618943945.root")
    #fileNames = cms.untracked.vstring("root://eoscms.cern.ch//eos/cms/store/group/comm_trigger/L1Trigger/apana/262AA156-744A-E311-9829-002618943945.root")   
    # fileNames = cms.untracked.vstring("root://eoscms.cern.ch//eos/cms/store/group/comm_trigger/L1Trigger/apana/DYJetsToLL_M-50_13TeV-pythia6_Fall13dr-tsg_PU40bx25__skim_150_1_6UN.root")
    #fileNames = cms.untracked.vstring("file:test.root")
    )


process.output = cms.OutputModule(
    "PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    # outputCommands = cms.untracked.vstring('keep *'),
    outputCommands = cms.untracked.vstring('drop *',
                                           'keep *_*_*_L1TEMULATION'),
    fileName = cms.untracked.string(OUTFILE),
    dataset = cms.untracked.PSet(
    filterName = cms.untracked.string(''),
    dataTier = cms.untracked.string('')
    )
                                           )
process.options = cms.untracked.PSet()

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
## process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS1', '')
process.GlobalTag.connect = cms.string('frontier://FrontierProd/CMS_COND_31X_GLOBALTAG')
process.GlobalTag.globaltag = cms.string('POSTLS162_V2::All')

#currently runs on gctDigis present in file already
#see L1TriggerEmulatorHI_newLayer2.py for example of 
#how to run on simRctDigis
process.rctLayer2Format = cms.EDProducer(
    "l1t::L1TCaloRCTToUpgradeConverter",
    regionTag = cms.InputTag("gctDigis"),
    emTag = cms.InputTag("gctDigis"))

process.Layer2HW = cms.EDProducer(
    "l1t::Stage1Layer2Producer",
    CaloRegions = cms.InputTag("rctLayer2Format"),
    CaloEmCands = cms.InputTag("rctLayer2Format"),
    FirmwareVersion = cms.uint32(2)  ## 1=HI algo, 2= pp algo
    )

process.Layer2Phys = cms.EDProducer("l1t::PhysicalEtAdder",
                                    InputCollection = cms.InputTag("Layer2HW")
)


process.Layer2gctFormat = cms.EDProducer("l1t::L1TCaloUpgradeToGCTConverter",
                                         InputCollection = cms.InputTag("Layer2Phys")
)

process.L1Packer = cms.EDProducer("l1t::L1TDigiToRaw",
        InputLabel = cms.InputTag("caloStage1"),
        FedId = cms.int32(100),
        FWId = cms.uint32(1),
        CaloTowers = cms.InputTag(""),
        EGammas = cms.InputTag("Layer2HW"),
        EtSums = cms.InputTag("Layer2HW"),
        Jets = cms.InputTag("Layer2HW"),
        Taus = cms.InputTag("Layer2HW"))

process.L1Unpacker = cms.EDProducer("l1t::L1TRawToDigi",
        InputLabel = cms.InputTag("L1Packer"),
        FedId = cms.int32(100))

process.gctDigis = cms.EDProducer("GctRawToDigi",
    unpackSharedRegions = cms.bool(False),
    numberOfGctSamplesToUnpack = cms.uint32(1),
    verbose = cms.untracked.bool(False),
    numberOfRctSamplesToUnpack = cms.uint32(1),
    inputLabel = cms.InputTag("rawDataCollector"),
    unpackerVersion = cms.uint32(0),
    gctFedId = cms.untracked.int32(745),
    hltMode = cms.bool(False)
)

process.Layer2 = cms.Sequence(
        process.rctLayer2Format
        *process.Layer2HW
        *process.Layer2Phys
        *process.Layer2gctFormat
        ## *process.L1Packer
        ## *process.L1Unpacker
        )


process.p1 = cms.Path(
    process.gctDigis*
    process.Layer2
    )

process.output_step = cms.EndPath(process.output)

process.schedule = cms.Schedule(
    process.p1, process.output_step
    )

# Spit out filter efficiency at the end.
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
