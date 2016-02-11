from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'DoubleMuon_Run2015B'
#config.General.workArea = ''

config.JobType.pluginName = 'Analysis'

config.JobType.psetName = 'RootTreeDA.py'

config.Data.inputDataset = '/DoubleMuon/Run2015B-PromptReco-v1/MINIAOD'

config.Data.lumiMask = ''

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/ekennedy/smh2mu/data/DoubleMuon_Run2015B-PromptReco-v1'
config.Data.publication = False

config.Site.storageSite = 'T2_CH_CERN'
