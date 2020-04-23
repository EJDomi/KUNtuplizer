# KUNtuplizer
To check out code:
```
setenv SCRAM_ARCH slc7_amd64_gcc700 ; ###chs/ tcsh 
export SCRAM_ARCH=slc7_amd64_gcc700 ; ### bash
cmsrel CMSSW_10_2_0
cd CMSSW_10_2_0/src
cmsenv
git cms-init
git clone -b performance_tuples git@github.com:EJDomi/KUNtuplizer.git Framework/KUNtuplizer
scram b -j4
```
To produce ntuple:
```
cd Framework/KUNtuplizer/test
cmsRun NtuplerMaker_cfg.py
```
