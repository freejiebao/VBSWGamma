from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName   = 'SMu16D-v1'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.inputFiles =['Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFchs.txt','Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFchs.txt','Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFchs.txt','Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFchs.txt','Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFPuppi.txt','Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFPuppi.txt','Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFPuppi.txt','Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFPuppi.txt']
# Name of the CMSSW configuration file
config.JobType.psetName    = 'analysis_data.py'
config.JobType.allowUndistributedCMSSW = True

config.section_("Data")
config.Data.inputDataset = '/SingleMuon/Run2016D-03Feb2017-v1/MINIAOD'
config.Data.inputDBS = 'global'

config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'Automatic'     #added by myself

config.Data.unitsPerJob = 40
config.Data.lumiMask = 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
#config.Data.runRange = '246908-258750'
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.outLFNDirBase = '/store/group/phys_jetmet/qihuang/'
#config.Data.outLFNDirBase = '/store/user/qihuang/'
config.Data.publication = False
config.Data.outputDatasetTag = 'SMu16D-v1'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'



