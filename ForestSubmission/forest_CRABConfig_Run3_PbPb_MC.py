from CRABClient.UserUtilities import config
from CRABClient.UserUtilities import getUsername
username = getUsername()

###############################################################################
# INPUT/OUTPUT SETTINGS

jobTag = 'Run3_PbPb_2025MC_ZToEE'
# JPsiToEE
#input = '/JpsiDielectron_pTHatMin4_HydjetEmbedded_Pythia8_TuneCP5_1510pre6/fdamas-PAT_151X_mcRun3_2025_realistic_HI_v1-5249b5d2d214ceff3b59bef72e572410/USER'
# ZToEE
input = '/JpsiDielectron_pTHatMin4_HydjetEmbedded_Pythia8_TuneCP5_1510pre6/fdamas-ZtoEE_PAT_151X_mcRun3_2025_realistic_HI_v1-81feeb1db1a253aa05da8c225ee26acc/USER'
# QCD Photon
#input = '/PhotonQCD_pTMin20_HydjetEmbedded_Pythia8_TuneCP5_1510pre6/fdamas-PAT_151X_mcRun3_2025_realistic_HI_v1-5249b5d2d214ceff3b59bef72e572410/USER'
inputDatabase = 'phys03'
output = '/store/group/phys_heavyions/' + username + '/Run3_PbPb_2025MC/'
outputServer = 'T2_CH_CERN'

###############################################################################

config = config()

config.General.requestName = jobTag
config.General.workArea = 'CrabWorkArea'
config.General.transferOutputs = True

config.JobType.psetName = 'forest_CMSSWConfig_Run3_PbPb_MC_miniAOD.py'
config.JobType.pluginName = 'Analysis'
config.JobType.maxMemoryMB = 3000 # CRAB gets crabby if this is >3000
config.JobType.pyCfgParams = ['noprint']
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = input
config.Data.inputDBS = inputDatabase
config.Data.outLFNDirBase = output
config.Data.splitting = 'FileBased' # 'EventAwareLumiBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = -1
config.Data.publication = False
config.Data.allowNonValidInputDataset = True

config.Site.storageSite = outputServer
