# HLT EGamma studies for 2025 PbPb run

Kyle Devereaux, MIT

## Foresting
Forest production happens in `./Forest Submission`. For the 2025 PbPb run we use `CMSSW_15_1_0`, for your current run use which ever CMSSW version is appropriate. The setup looks like this
```
cmsrel CMSSW_15_1_0
cd CMSSW_15_1_0/src
cmsenv
git cms-merge-topic CmsHI:forest_CMSSW_15_1_X
scram build -j4
voms-proxy-init --voms cms --valid 168:00
```
CRAB file to submit foresting jobs is at `forest_CRABConfig_Run3_PbPb_MC.py`. Change the input minAOD path and output parameters. Also change the Global Tag (GT) and Configuration in `forest_CMSSWConfig_Run3_PbPb_MC_miniAOD.py` to the most up-to-date choices. For HLT studies ensure the forests to have the following information:
```
hltanalysis : info on L1 and HLT triggers that fire
hiEvtAnalyzer : event-level information
ggHiNtuplizer : electron and photon track information
hltobject : info on objects that fire the HLT trigger
```

2025 PbPb MC samples: https://docs.google.com/spreadsheets/d/1hnxPcfO3gXvGtVTpd9hHhg-AhlyECpGRMtRAf6NRXco/edit?gid=1414989721#gid=1414989721

2025 OO MC samples: https://docs.google.com/spreadsheets/d/1hnxPcfO3gXvGtVTpd9hHhg-AhlyECpGRMtRAf6NRXco/edit?gid=222281008#gid=222281008

For foresting guide see:
https://github.com/jdlang/MITHIGAnalysis2024/tree/mainOO/ForestSubmission/20250624_ForestSubmission_OO2025_FastPrivateReco

For current foresting scripts see: (click through to your current CMSSW release)
https://github.com/CmsHI/cmssw/tree/forest_CMSSW_15_0_X/HeavyIonsAnalysis/Configuration/test

Recipe to go from Gen-Sim to miniAOD are here: https://codimd.web.cern.ch/URi6opRRTX2z2LmK7XTpNg#


## HLT emulation
YOU DO NOT NEED TO EMULATE IF THE FOREST IS UPTO DATE WITH THE HLT MENU OR YOU ARE STUDYING DATA.

Emulation is only required if you are studying MC, and the forests you produced do not include HLT paths you want to study. Such is the case when you are testing out a new menu with new paths. Emulation is performed on RECO-level files and reproduces information on L1/HLT trigger descions, as well as (pT, eta, phi) information on the track which fired the trigger. Forests are produced downstream from the same RECO files. The HLT emulation events are matched to forest events  with a special strategy (hash table, etc.) which match on their (run, lumi section, event) information.

Not implemented here.

For implementation follow:
https://github.com/claytoniousfunk/HLT_emulation/tree/main

General emulation info: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGlobalHLT


## Trigger efficiency studies
Trigger analysis scripts are all named `triggerAnalysis_[ele/pho]_[mc/prompt/PbPb].C`.
Run with `root -l -q -b scriptname.C`.
The scripts have inputs
```
inputForest : path to forests you made in the Foresting step
inputHLT : path to HLT emulation files for matching to forests [not implemented right now]
output_base : name tag to distinguish outputs
nfiles : total number of forest files to process
minHiBin : minimum centrality cut, multiply percentage by 2, for example 30% -> minHiBin=60
maxHiBin : maximum centrality cut  
```
The scripts should produce turn-on curve plots for electron/photon HLT performance and L1 performance in the given centrality bin, storing plots to `plots/` and histograms to `output/`.

The HLT paths studied, L1 paths,  quality cuts, Gen-matching criteria, HLT emulation matching criteria, etc. can be modified within the scripts.

Scripts borrow from: https://github.com/claytoniousfunk/triggerAnalysis

And from: https://github.com/pinchunchou/HLT2024


## Making/editing a HLT menus

Not handled in this repo, but just for general know-how ...

### ConfDB
Where you create a new menu or make edits to an existing one. If you just need to view a menu without editing, use https://cmshltcfg.app.cern.ch/

### Setup
https://indico.cern.ch/event/1131290/contributions/4747094/attachments/2394836/4094536/HI_HLT_Tutorial_Feb-21.pdf

### Running
(on local machine)

```
ssh -XY <username>@lxplus929.cern.ch
cd hlt-confdb/
ant gui
vncserver -geometry 1880x970 :2
```

Generally can be ran on any lxplusXXX terminal, but I request 929 so it is consistent across login attempts.

(on local machine, new terminal)

```
ssh -fN -Y -L 5902:localhost:5902 <username>@lxplus929.cern.ch
```

(in VNC Connect application on local)

Open VNC Connect, connect to `localhost:2` . Login using your VNC password when you set during setup. It should open a remote desktop on lxplus.

(in remote desktop)

Open Applications/System Tools/Terminal.

```
cd hlt-confdb
./start
```

(in ConfDBGUI in remote desktop )

Establish database connection via

```
Setup: Run 3 Offline
Password: --------
```