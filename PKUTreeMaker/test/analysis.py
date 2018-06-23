import FWCore.ParameterSet.Config as cms

process = cms.Process( "TEST" )
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True),
				     SkipEvent = cms.untracked.vstring('ProductNotFound'))
corrJetsOnTheFly = True
runOnMC = True
chsorpuppi = True  # AK4Chs or AK4Puppi
#****************************************************************************************************#
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
if runOnMC:
   process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'
elif not(runOnMC):
   process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7'

##########			                                                             
hltFiltersProcessName = 'RECO'
if runOnMC:
   hltFiltersProcessName = 'PAT' #'RECO'
reducedConversionsName = 'RECO'
if runOnMC:
   reducedConversionsName= 'PAT' #'RECO'

process.load("VAJets.PKUCommon.goodMuons_cff")
process.load("VAJets.PKUCommon.goodElectrons_cff")
process.load("VAJets.PKUCommon.goodPhotons_cff")
process.load("VAJets.PKUCommon.leptonicW_cff")

# If Update
process.goodMuons.src = "slimmedMuons"
process.goodElectrons.src = "slimmedElectrons"
process.goodPhotons.src = "slimmedPhotons"
process.Wtoenu.MET  = "slimmedMETs"
process.Wtomunu.MET = "slimmedMETs"

process.load("VAJets.PKUCommon.goodJets_cff") 
if chsorpuppi:
      process.goodAK4Jets.src = "slimmedJets"
else:
      process.goodAK4Jets.src = "slimmedJetsPuppi"
 

#process.goodOfflinePrimaryVertex = cms.EDFilter("VertexSelector",
#                                       src = cms.InputTag("offlineSlimmedPrimaryVertices"),
#                                       cut = cms.string("chi2!=0 && ndof >= 4.0 && abs(z) <= 24.0 && abs(position.Rho) <= 2.0"),
#                                       filter = cms.bool(False)
#                                       )

WBOSONCUT = "pt > 0.0"

process.leptonicVSelector = cms.EDFilter("CandViewSelector",
                                       src = cms.InputTag("leptonicV"),
                                       cut = cms.string( WBOSONCUT ), 
                                       filter = cms.bool(False)
                                       )

process.leptonicVFilter = cms.EDFilter("CandViewCountFilter",
                                       src = cms.InputTag("leptonicV"),
                                       minNumber = cms.uint32(0),
                                       filter = cms.bool(False)
                                       )


process.leptonSequence = cms.Sequence(process.muSequence +
                                      process.eleSequence +
                                      process.leptonicVSequence +
                                      process.leptonicVSelector +
                                      process.leptonicVFilter )

process.jetSequence = cms.Sequence(process.NJetsSequence)


process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.load("RecoMET.METFilters.BadChargedCandidateFilter_cfi")
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.metfilterSequence = cms.Sequence(process.BadPFMuonFilter+process.BadChargedCandidateFilter)
#begin------------JEC on the fly--------
if runOnMC:
   jecLevelsAK4chs = [
          'Summer16_23Sep2016V4_MC_L1FastJet_AK4PFchs.txt',
          'Summer16_23Sep2016V4_MC_L2Relative_AK4PFchs.txt',
          'Summer16_23Sep2016V4_MC_L3Absolute_AK4PFchs.txt'
    ]
   jecLevelsAK4puppi = [
          'Summer16_23Sep2016V4_MC_L1FastJet_AK4PFPuppi.txt',
          'Summer16_23Sep2016V4_MC_L2Relative_AK4PFPuppi.txt',
          'Summer16_23Sep2016V4_MC_L3Absolute_AK4PFPuppi.txt'
    ]
else:
   jecLevelsAK4chs = [
          'Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFchs.txt',
          'Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFchs.txt',
          'Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFchs.txt',
	  'Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFchs.txt'
    ]
   jecLevelsAK4puppi = [
          'Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFPuppi.txt',
          'Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFPuppi.txt',
          'Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFPuppi.txt',
	  'Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFPuppi.txt'
    ]
#end------------JEC on the fly--------

if chsorpuppi:
      ak4jecsrc = jecLevelsAK4chs
else:
      ak4jecsrc = jecLevelsAK4puppi
 
process.load("RecoEgamma/PhotonIdentification/PhotonIDValueMapProducer_cfi")
   
process.treeDumper = cms.EDAnalyzer("PKUTreeMaker",
                                    originalNEvents = cms.int32(1),
                                    crossSectionPb = cms.double(1),
                                    targetLumiInvPb = cms.double(1.0),
                                    PKUChannel = cms.string("VW_CHANNEL"),
                                    isGen = cms.bool(False),
				    RunOnMC = cms.bool(runOnMC), 
                                    generator =  cms.InputTag("generator"),
#                                    lhe =  cms.InputTag("externalLHEProducer"),
                                    pileup  =   cms.InputTag("slimmedAddPileupInfo"),  
                                    leptonicVSrc = cms.InputTag("leptonicV"),
                                    rho = cms.InputTag("fixedGridRhoFastjetAll"),   
                                    ak4jetsSrc = cms.InputTag("cleanAK4Jets"),      
#                                    photonSrc = cms.InputTag("goodPhotons"),
                                    photonSrc = cms.InputTag("slimmedPhotons"),
                                    genSrc =  cms.InputTag("prunedGenParticles"),  
                                    jecAK4chsPayloadNames = cms.vstring( jecLevelsAK4chs ),
                                    jecAK4PayloadNames = cms.vstring( ak4jecsrc ),
                                    metSrc = cms.InputTag("slimmedMETs"),
                                    vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),  
                                    t1jetSrc = cms.InputTag("slimmedJets"),      
                                    t1muSrc = cms.InputTag("slimmedMuons"),       
                                    looseelectronSrc = cms.InputTag("vetoElectrons"),
                                    electrons = cms.InputTag("slimmedElectrons"),
                                    conversions = cms.InputTag("reducedEgamma","reducedConversions",reducedConversionsName),
                                    beamSpot = cms.InputTag("offlineBeamSpot","","RECO"),
                                    loosemuonSrc = cms.InputTag("looseMuons"),
                                    hltToken    = cms.InputTag("TriggerResults","","HLT"),
                                    elPaths1     = cms.vstring("HLT_Ele23_WPTight_Gsf_v*"),
                                    elPaths2     = cms.vstring("HLT_Ele27_WPTight_Gsf_v*"),
                                    muPaths1     = cms.vstring("HLT_IsoMu20_v*","HLT_IsoTkMu20_v*"),
#                                    muPaths2     = cms.vstring("HLT_IsoMu22_v*","HLT_IsoTkMu22_v*"),
                                    muPaths2     = cms.vstring("HLT_IsoMu24_v*","HLT_IsoTkMu24_v*"),
                                    muPaths3     = cms.vstring("HLT_IsoMu27_v*","HLT_IsoTkMu27_v*"),
				    noiseFilter = cms.InputTag('TriggerResults','', hltFiltersProcessName),
				    noiseFilterSelection_HBHENoiseFilter = cms.string('Flag_HBHENoiseFilter'),
                                    noiseFilterSelection_HBHENoiseIsoFilter = cms.string("Flag_HBHENoiseIsoFilter"),
				    noiseFilterSelection_globalTightHaloFilter = cms.string('Flag_globalTightHalo2016Filter'),
                                    noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter = cms.string('Flag_EcalDeadCellTriggerPrimitiveFilter'),
				    noiseFilterSelection_goodVertices = cms.string('Flag_goodVertices'),
				    noiseFilterSelection_eeBadScFilter = cms.string('Flag_eeBadScFilter'),
                                    noiseFilterSelection_badMuon = cms.InputTag('BadPFMuonFilter'),
                                    noiseFilterSelection_badChargedHadron = cms.InputTag('BadChargedCandidateFilter'),
                                    full5x5SigmaIEtaIEtaMap   = cms.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta"),
                                    phoChargedIsolation = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
                                    phoNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
                                    phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
                                    effAreaChHadFile = cms.FileInPath("RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfChargedHadrons_90percentBased.txt"),
                                    effAreaNeuHadFile= cms.FileInPath("RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased.txt"),
                                    effAreaPhoFile   = cms.FileInPath("RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfPhotons_90percentBased.txt")
#                                    effAreaChHadFile = cms.FileInPath("RecoEgamma/PhotonIdentification/data/Spring15/effAreaPhotons_cone03_pfChargedHadrons_25ns_NULLcorrection.txt"),
#                                    effAreaNeuHadFile= cms.FileInPath("RecoEgamma/PhotonIdentification/data/Spring15/effAreaPhotons_cone03_pfNeutralHadrons_25ns_90percentBased.txt"),
#                                    effAreaPhoFile   = cms.FileInPath("RecoEgamma/PhotonIdentification/data/Spring15/effAreaPhotons_cone03_pfPhotons_25ns_90percentBased.txt")
                                    )


process.analysis = cms.Path(
#                            process.goodOfflinePrimaryVertex +
                            process.leptonSequence +
                            process.jetSequence +
                            process.metfilterSequence +
#                           process.photonSequence +
                            process.photonIDValueMapProducer*process.treeDumper)

### Source
process.load("VAJets.PKUCommon.data.RSGravitonToWW_kMpl01_M_1000_Tune4C_13TeV_pythia8")
process.source.fileNames = [
#"/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/04E807A8-9EBE-E611-ABE7-FA163EC9F0BF.root"
"/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/165F54A0-A3BE-E611-B3F7-0025905A606A.root"
]
                       
process.maxEvents.input =10000

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.MessageLogger.cerr.FwkReport.limit = 99999999

process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string("treePKU.root")
                                   )