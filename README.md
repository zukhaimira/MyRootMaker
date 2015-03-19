# MyRootMaker
MyRootMaker ntuplizer 

https://twiki.cern.ch/twiki/bin/view/Sandbox/MyRootMakerFrom70XTo72X

Recipe:

    cmsrel CMSSW_7_2_4
    cd CMSSW_7_2_4/src/
    cmsenv
    git cms-init

    git cms-merge-topic HuguesBrun:trigElecIdInCommonIsoSelection720
    cp /afs/cern.ch/user/i/ikrav/public/EGMCode/GsfEleFull5x5SigmaIEtaIEtaCut72X.cc RecoEgamma/ElectronIdentification/plugins/cuts/
    cp /afs/cern.ch/user/i/ikrav/public/EGMCode/cutBasedElectronID_PHYS14_PU20bx25_V0_miniAOD_cff.py RecoEgamma/ElectronIdentification/python/Identification/

    git clone https://github.com:ekenn003/PFIsolation.git
    git clone https://github.com:ekenn003/MyRootMaker.git

    scram b -j 8
    
    cmsRun MyRootMaker/MyRootMaker/RootTreeMC_mini.py    
