# MyRootMaker
MyRootMaker ntuplizer 

Note: electron ID doesn't work if used in a release less than 747.

https://twiki.cern.ch/twiki/bin/view/Sandbox/MyRootMakerFrom72XTo74X

Recipe:

    cmsrel CMSSW_7_4_7
    cd CMSSW_7_4_7/src/
    cmsenv
    git cms-init

    git clone https://github.com/ekenn003/PFIsolation.git
    git clone https://github.com/ekenn003/MyRootMaker.git

    scram b -j 8
    
    cmsRun MyRootMaker/MyRootMaker/RootTreeMC_mini.py (monte carlo)  
    cmsRun MyRootMaker/MyRootMaker/RootTreeDA_mini.py (data)
