#!/bin/sh
g++ HLTperformance_HI2023.cpp -o Execute_HLTperformance `root-config --cflags --libs` -lASImage

#./Execute_HLTperformance --trigType 0 > logs/PbPb0.log 2> logs/PbPb0.err < /dev/null &
#./Execute_HLTperformance --trigType 1 > logs/PbPb1.log 2> logs/PbPb1.err < /dev/null &
#./Execute_HLTperformance --trigType 2 > logs/PbPb2.log 2> logs/PbPb2.err < /dev/null &
#./Execute_HLTperformance --trigType 3 > logs/PbPb3.log 2> logs/PbPb3.err < /dev/null &
#./Execute_HLTperformance --trigType 4 > logs/PbPb4.log 2> logs/PbPb4.err < /dev/null &
#./Execute_HLTperformance --trigType 5 > logs/PbPb5.log 2> logs/PbPb5.err < /dev/null &
#./Execute_HLTperformance --trigType 6 > logs/PbPb6.log 2> logs/PbPb6.err < /dev/null &
#./Execute_HLTperformance --trigType 7 > logs/PbPb7.log 2> logs/PbPb7.err < /dev/null &
#./Execute_HLTperformance --trigType 8 > logs/PbPb8.log 2> logs/PbPb8.err < /dev/null &

#./Execute_HLTperformance --trigType 0 --isPbPb false > logs/pp0.log 2> logs/pp0.err < /dev/null &
#./Execute_HLTperformance --trigType 1 --isPbPb false > logs/pp1.log 2> logs/pp1.err < /dev/null &
#./Execute_HLTperformance --trigType 2 --isPbPb false > logs/pp2.log 2> logs/pp2.err < /dev/null &
#./Execute_HLTperformance --trigType 3 --isPbPb false > logs/pp3.log 2> logs/pp3.err < /dev/null &
#./Execute_HLTperformance --trigType 4 --isPbPb false > logs/pp4.log 2> logs/pp4.err < /dev/null &
#./Execute_HLTperformance --trigType 5 --isPbPb false > logs/pp5.log 2> logs/pp5.err < /dev/null &
#./Execute_HLTperformance --trigType 6 --isPbPb false > logs/pp6.log 2> logs/pp6.err < /dev/null &
#./Execute_HLTperformance --trigType 7 --isPbPb false > logs/pp7.log 2> logs/pp7.err < /dev/null &

#./Execute_HLTperformance --trigType 0 --noL1 true > logs/PbPb0_noL1.log 2> logs/PbPb0_noL1.err < /dev/null &
#./Execute_HLTperformance --trigType 1 --noL1 true > logs/PbPb1_noL1.log 2> logs/PbPb1_noL1.err < /dev/null &
#./Execute_HLTperformance --trigType 2 --noL1 true > logs/PbPb2_noL1.log 2> logs/PbPb2_noL1.err < /dev/null &
#./Execute_HLTperformance --trigType 3 --noL1 true > logs/PbPb3_noL1.log 2> logs/PbPb3_noL1.err < /dev/null &
#./Execute_HLTperformance --trigType 8 --noL1 true > logs/PbPb8_noL1.log 2> logs/PbPb8_noL1.err < /dev/null &
#./Execute_HLTperformance --trigType 0 --isPbPb false --noL1 true > logs/pp0_noL1.log 2> logs/pp0_noL1.err < /dev/null &
#./Execute_HLTperformance --trigType 1 --isPbPb false --noL1 true > logs/pp1_noL1.log 2> logs/pp1_noL1.err < /dev/null &
#./Execute_HLTperformance --trigType 2 --isPbPb false --noL1 true > logs/pp2_noL1.log 2> logs/pp2_noL1.err < /dev/null &
#./Execute_HLTperformance --trigType 3 --isPbPb false --noL1 true > logs/pp3_noL1.log 2> logs/pp3_noL1.err < /dev/null &

#dataDir=/eos/cms/store/group/phys_heavyions/jviinika/run3RapidValidation/PbPb2023_run374289_HIPhysicsRawPrime0_withDFinder_2023-09-26/0000/
#dataDir=/eos/cms/store/group/phys_heavyions/jmijusko/run3RapidValidation/PbPb2023_run374322_PhysicsHIPhysicsRawPrime0_withDFinder_2023-09-28/0000/
#dataDir=/eos/cms/store/group/phys_heavyions/jviinika/run3RapidValidation/PbPb2023_run374345_HIPhysicsRawPrime0_triggerObjects_2023-09-29/CRAB_UserFiles/crab_PbPb2023_run374345_HIPhysicsRawPrime0_triggerObjects_2023-09-29/230930_011704/0000/
#dataDir="/eos/cms/store/group/phys_heavyions/wangj/Forest2023/HIPhysicsRawPrime[3,4,5,6,8,9]/crab_HiForestMINIAOD_230930_HIPhysicsRawPrime[3,4,5,6,8,9]_374354_Dpt2trk1/*/0000/"
#dataDir="/eos/cms/store/group/phys_heavyions/jviinika/run3RapidValidation/HIPhysicsRawPrime0_HIRun2023A-PromptReco-v1_run374322_2023-09-30/0000/"

# KD: commented this out
#dataDir="/eos/cms/store/group/phys_heavyions/wangj/Forest2023/HIPhysicsRawPrime9/crab_HiForestMINIAOD_230930_HIPhysicsRawPrime9_374354_Dpt2trk1/231001_011507/0000/*.root",
#dataDir+="/eos/cms/store/group/phys_heavyions/wangj/Forest2023/HIPhysicsRawPrime8/crab_HiForestMINIAOD_230930_HIPhysicsRawPrime8_374354_Dpt2trk1/231001_011439/0000/*.root",
#dataDir+="/eos/cms/store/group/phys_heavyions/wangj/Forest2023/HIPhysicsRawPrime6/crab_HiForestMINIAOD_230930_HIPhysicsRawPrime6_374354_Dpt2trk1/231001_011111/0000/*.root",
#dataDir+="/eos/cms/store/group/phys_heavyions/wangj/Forest2023/HIPhysicsRawPrime5/crab_HiForestMINIAOD_230930_HIPhysicsRawPrime5_374354_Dpt2trk1/231001_010540/0000/*.root",
#dataDir+="/eos/cms/store/group/phys_heavyions/wangj/Forest2023/HIPhysicsRawPrime4/crab_HiForestMINIAOD_230930_HIPhysicsRawPrime4_374354_Dpt2trk1/231001_010449/0000/*.root",
#dataDir+="/eos/cms/store/group/phys_heavyions/wangj/Forest2023/HIPhysicsRawPrime3/crab_HiForestMINIAOD_230930_HIPhysicsRawPrime3_374354_Dpt2trk1/231001_010436/0000/*.root"

#runtext="Run 374322,374345,374354 (5.36 TeV)"
runtext=""
folder="figs"
logpath="logs/log"
suffixText=""
suffixText1=""
drMax=0.5

######################### 2024 MC Zee sample #########################
dataDir="/eos/cms/store/group/phys_heavyions/prdas/EGamma/Run3_PbPb_2024_MC/Ze10e10/Embedded/Pythia8_Embedded_Ze10e10_TuneCP5_2024/crab_20241026_145435/241026_125438/0000/*.root"
tag="2024PbPb_Zee"
plotLabel="Pinchun's code, 2024 PbPb MC ZtoEE"

# MC single ele low pT
#./Execute_HLTperformance --trigType 4 --isL1denom true --nocut false --isHLTObj false --isMC true --isPbPb true --trigsuf "13" --PSvec 1,1,1     --LPSvec 1,1,1 --L1ID 2,2,3 --dataDir "$dataDir" --plotLabel "${plotLabel}" --suffix "${tag}_SingleE_lowPt" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_${tag}_isL1denom_0.log 2>${logpath}_${tag}_isL1denom_0.err < /dev/null &

# MC single ele high pT
./Execute_HLTperformance --trigType 5 --isL1denom true --nocut false --isHLTObj false --isMC true --isPbPb true --trigsuf "13" --PSvec 1,1,1     --LPSvec 1,1,1 --L1ID 3,4,4 --dataDir "$dataDir" --plotLabel "${plotLabel}" --suffix "${tag}_SingleE_highPt" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_${tag}_isL1denom_0.log 2>${logpath}_${tag}_isL1denom_0.err < /dev/null &

######################### 2024 MC Zee sample, no embdedding ####################
dataDir="/eos/cms/store/group/phys_heavyions/prdas/EGamma/Run3_PbPb_2024_MC/Ze10e10/NoEmbedding/Pythia8_Embedded_Ze10e10_TuneCP5_2024/crab_20241026_145054/241026_125057/0000/*.root"
tag="2024PbPb_Zee_noembedding"
plotLabel="Pinchun's code, 2024 PbPb MC ZtoEE, no embedding"

# MC single ele low pT
#./Execute_HLTperformance --trigType 4 --isL1denom true --nocut false --isHLTObj false --isMC true --isPbPb true --trigsuf "13" --PSvec 1,1,1     --LPSvec 1,1,1 --L1ID 2,2,3 --dataDir "$dataDir" --plotLabel "${plotLabel}" --suffix "${tag}_SingleE_lowPt" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_${tag}_isL1denom_0.log 2>${logpath}_${tag}_isL1denom_0.err < /dev/null &

# MC single ele high pT
#./Execute_HLTperformance --trigType 5 --isL1denom true --nocut false --isHLTObj false --isMC true --isPbPb true --trigsuf "13" --PSvec 1,1,1     --LPSvec 1,1,1 --L1ID 3,4,4 --dataDir "$dataDir" --plotLabel "${plotLabel}" --suffix "${tag}_SingleE_highPt" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_${tag}_isL1denom_0.log 2>${logpath}_${tag}_isL1denom_0.err < /dev/null &



######################### 2025 MC Zee sample #########################
dataDir="/eos/cms/store/group/phys_heavyions/kdeverea/Run3_PbPb_2025MC/JpsiDielectron_pTHatMin4_HydjetEmbedded_Pythia8_TuneCP5_1510pre6/crab_Run3_PbPb_2025MC_ZToEE/251115_102359/0000/merged.root"
tag="2025PbPb_Zee"
plotLabel="Pinchun's code, 2025 PbPb MC ZtoEE"

# MC single ele low pT
#./Execute_HLTperformance --trigType 4 --isL1denom true --nocut false --isHLTObj false --isMC true --isPbPb true --trigsuf "16" --PSvec 1,1,1     --LPSvec 1,1,1 --L1ID 2,2,3 --dataDir "$dataDir" --plotLabel "${plotLabel}" --suffix "${tag}_SingleE_lowPt" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_${tag}_isL1denom_0.log 2>${logpath}_${tag}_isL1denom_0.err < /dev/null &

# MC single ele high pT
#./Execute_HLTperformance --trigType 5 --isL1denom true --nocut false --isHLTObj false --isMC true --isPbPb true --trigsuf "16" --PSvec 1,1,1     --LPSvec 1,1,1 --L1ID 3,4,4 --dataDir "$dataDir" --plotLabel "${plotLabel}" --suffix "${tag}_SingleE_highPt" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_${tag}_isL1denom_0.log 2>${logpath}_${tag}_isL1denom_0.err < /dev/null &


###################### 2025 Quick reco PbPb sample ######################
fileList="../run399499_forests.txt"
tag="2025PbPb_quickreco399499"
plotLabel="Pinchun's code, 2025 PbPb run 399499"

# Data single ele low pT
#./Execute_HLTperformance --trigType 4 --isL1denom true --nocut false --isHLTObj false --isMC false --isPbPb true --trigsuf "17" --PSvec 1,1,1     --LPSvec 1,1,1 --L1ID 2,2,3 --fileList "$fileList" --dataDir "$dataDir" --plotLabel "${plotLabel}" --suffix "${tag}_SingleE_lowPt" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_${tag}_isL1denom_0.log 2>${logpath}_${tag}_isL1denom_0.err < /dev/null &

# MC single ele high pT
#./Execute_HLTperformance --trigType 5 --isL1denom true --nocut false --isHLTObj false --isMC false --isPbPb true --trigsuf "16" --PSvec 1,1,1     --LPSvec 1,1,1 --L1ID 3,4,4 --fileList "$fileList" --dataDir "$dataDir" --plotLabel "${plotLabel}" --suffix "${tag}_SingleE_highPt" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_${tag}_isL1denom_0.log 2>${logpath}_${tag}_isL1denom_0.err < /dev/null &



#./Execute_HLTperformance --trigType 0 --isL1denom true --nocut false --isHLTObj true --isMC false --trigsuf "9" --PSvec 1,1,1     --LPSvec 1,1,1 --L1ID 4,4,5 --dataDir "$dataDir" --suffix "$suffixText" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_isL1denom_0.log 2>${logpath}_isL1denom_0.err < /dev/null &
#./Execute_HLTperformance --trigType 0 --isL1denom false --nocut false --isHLTObj true --isMC false --trigsuf "9" --PSvec 1,1,1     --LPSvec 1,1,1 --L1ID 4,4,5 --dataDir "$dataDir" --suffix "$suffixText" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_0.log 2>${logpath}_0.err < /dev/null &
#./Execute_HLTperformance --trigType 1 --isL1denom true --nocut false --isHLTObj true --isMC false --trigsuf "8" --PSvec 1700,10,1 --LPSvec 3,3,3 --L1ID 2,2,2 --dataDir "$dataDir" --suffix "$suffixText" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_1.log 2>${logpath}_1.err < /dev/null &
#./Execute_HLTperformance --trigType 2 --isL1denom true --nocut false --isHLTObj true --isMC false --trigsuf "9" --PSvec 1,1,1     --LPSvec 1,1,1 --L1ID 4,4,5 --dataDir "$dataDir" --suffix "$suffixText" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_isL1denom_2.log 2>${logpath}_isL1denom_2.err < /dev/null &
#./Execute_HLTperformance --trigType 2 --isL1denom false --nocut false --isHLTObj true --isMC false --trigsuf "9" --PSvec 1,1,1     --LPSvec 1,1,1 --L1ID 4,4,5 --dataDir "$dataDir" --suffix "$suffixText" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_2.log 2>${logpath}_2.err < /dev/null &
#./Execute_HLTperformance --trigType 3 --isL1denom true --nocut false --isHLTObj true --isMC false --trigsuf "8" --PSvec 650,5,1   --LPSvec 3,3,3 --L1ID 2,2,2 --dataDir "$dataDir" --suffix "$suffixText" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_3.log 2>${logpath}_3.err < /dev/null &
#./Execute_HLTperformance --trigType 4 --isL1denom true --nocut false --isHLTObj true --isMC false --trigsuf "8" --PSvec 20,7,1    --LPSvec 3,3,1 --L1ID 2,2,3 --dataDir "$dataDir" --suffix "$suffixText" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_4.log 2>${logpath}_4.err < /dev/null &
#./Execute_HLTperformance --trigType 5 --isL1denom true --nocut false --isHLTObj true --isMC false --trigsuf "8" --PSvec 1,1,1     --LPSvec 1,1,1 --L1ID 3,4,4 --dataDir "$dataDir" --suffix "$suffixText" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_5.log 2>${logpath}_5.err < /dev/null &
#./Execute_HLTperformance --trigType 6 --isL1denom true --nocut false --isHLTObj true --isMC false --trigsuf "8" --PSvec 1,1,1     --LPSvec 3,3,3 --L1ID 6,6,6 --dataDir "$dataDir" --suffix "$suffixText" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_6.log 2>${logpath}_6.err < /dev/null &
#./Execute_HLTperformance --trigType 7 --isL1denom true --nocut false --isHLTObj true --isMC false --trigsuf "8" --PSvec 1,1,1     --LPSvec 3,3,3 --L1ID 6,6,6 --dataDir "$dataDir" --suffix "$suffixText" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_7.log 2>${logpath}_7.err < /dev/null &
#./Execute_HLTperformance --trigType 8 --isL1denom true --nocut false --isHLTObj true --isMC false --trigsuf "1" --PSvec 1,1,1     --LPSvec 3,3,3 --L1ID 6,6,6 --dataDir "$dataDir" --suffix "$suffixText" --runtext "$runtext" --folder $folder --dRmax $drMax  > ${logpath}_8.log 2>${logpath}_8.err < /dev/null &