# MyRootMaker
MyRootMaker ntuplizer 

https://twiki.cern.ch/twiki/bin/view/Sandbox/MyRootMakerFrom70XTo72X
https://twiki.cern.ch/twiki/bin/view/Sandbox/MyRootMakerFrom72XTo74X

Recipe:

    cmsrel CMSSW_7_4_6_patch2
    cd CMSSW_7_4_6_patch2/src/
    cmsenv
    git cms-init

    git clone https://github.com/ekenn003/PFIsolation.git
    git clone https://github.com/ekenn003/MyRootMaker.git

    scram b -j 8
    
    cmsRun MyRootMaker/MyRootMaker/RootTreeMC_mini.py    
