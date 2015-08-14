This package can be installed in a recent CMSSW area e.g. CMSSW_7_4_5. inside RecoTagTools 

Example installation:
```
cmsrel CMSSW_7_4_5
cd CMSSW_7_4_5/src
cmsenv
git cms-init
git cms-addpkg RecoTauTag/RecoTau
git clone https://github.com/lmdodd/tauAnalysis.git RecoTauTag/tauAnalysis
scram b -j8
```


```
cd test 
cmsRun runMINIAOD_FR.py
cmsRun runMINIAOD_FR_QCD.py
cmsRun runAODSIM_effi.py
python plot_Mini_effi.py [effi].root label[optional]
python plot_Mini_FR_Jets.py [fakeRate].root label[optional]
python plot_Mini_FR_QCD.py [fakeRate_QCD].root label[optional]
```

note: in plot_def.py change the save directory to your favorite place. e.g. not outputs 

At Wisconsin for farmout- the files submitted are in the various .txt files in /test. The files will appear in /hdfs/store/user/[user]

```
cd test
source farmout.sh [JOB_NAME] #for standard isolation cones
python plot.py [effi].root label[optional]
python plotFR.py [fakeRate].root label[optional]
```

Please note, the dataMCcomp plugin is running but the corresponding plotting script is under construction.
