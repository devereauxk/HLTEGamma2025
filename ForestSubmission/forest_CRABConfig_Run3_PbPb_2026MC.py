from CRABClient.UserUtilities import config
from CRABClient.UserUtilities import getUsername
username = getUsername()

###############################################################################
# INPUT/OUTPUT SETTINGS

# 2026 PbPb MC samples

# QCD Photon
#jobTag = 'Run3_PbPb_2026MC_ZToEE'
#input = '/PhotonQCD_pTHatMin20_HydjetEmbedded_1610pre3/fdamas-PATwith161pre3_151X_mcRun3_2025_realistic_HI_v5-5c048b2d868b98d3d5d21d398292d46c/USER'
#inputDatabase = 'phys03'
#output = '/store/group/phys_heavyions/' + username + '/Run3_PbPb_2026MC/'
#outputServer = 'T2_CH_CERN'

# J/psi -> EE
#jobTag = 'Run3_PbPb_2026MC_JpsiToEE'
#input = '/JpsiDielectron_pTHatMin4_HydjetEmbedded_Pythia8_TuneCP5_1610pre3/fdamas-PATwith161pre3_151X_mcRun3_2025_realistic_HI_v5-5c048b2d868b98d3d5d21d398292d46c/USER'
#inputDatabase = 'phys03'
#output = '/store/group/phys_heavyions/' + username + '/Run3_PbPb_2026MC/'
#outputServer = 'T2_CH_CERN'

# DY
jobTag = 'Run3_PbPb_2026MC_DY'
input = '/DrellYan_HighMass_MadGraph_HydjetEmbedded_1610pre3/fdamas-PATwith161pre4_151X_mcRun3_2025_realistic_HI_v5-eaa0399b9218a690ee453ab5f1aeb831/USER'
inputDatabase = 'phys03'
output = '/store/group/phys_heavyions/' + username + '/Run3_PbPb_2026MC/'
outputServer = 'T2_CH_CERN'

###############################################################################

config = config()

config.General.requestName = jobTag
config.General.workArea = 'CrabWorkArea'
config.General.transferOutputs = True

config.JobType.psetName = 'forest_CMSSWConfig_Run3_PbPb_2026MC_miniAOD.py'
config.JobType.pluginName = 'Analysis'
config.JobType.maxMemoryMB = 3000 # CRAB gets crabby if this is >3000
config.JobType.pyCfgParams = ['noprint']
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = input
config.Data.inputDBS = inputDatabase
config.Data.outLFNDirBase = output
config.Data.splitting = 'Automatic' # 'EventAwareLumiBased' 
config.Data.unitsPerJob = 180
config.Data.totalUnits = -1
config.Data.publication = False
config.Data.allowNonValidInputDataset = True

config.Site.storageSite = outputServer
