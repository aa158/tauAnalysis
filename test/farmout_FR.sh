#!/bin/bash

EXPECTED_ARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: $0 JOB_NAME"
fi

farmoutAnalysisJobs $1-WJets_HTrange \
  --input-files-per-job=1 \
  --job-generates-output-name \
  --infer-cmssw-path \
  --input-file-list=WJets_HTrange.txt \
  --input-dir=root://cmsxrootd.fnal.gov/ \
  --assume-input-files-exist \
  ./runMINIAOD_FR.py  \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName'

farmoutAnalysisJobs $1-zp_M500 \
  --input-files-per-job=1 \
  --job-generates-output-name \
  --infer-cmssw-path \
  --input-file-list=Ztautau_M500.txt \
  --input-dir=root://cmsxrootd.fnal.gov/ \
  --assume-input-files-exist \
  ./runMINIAOD_FR.py  \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName'

farmoutAnalysisJobs $1-DYJets \
  --input-files-per-job=1 \
  --job-generates-output-name \
  --infer-cmssw-path \
  --input-file-list=DYJets_M50.txt \
  --input-dir=root://cmsxrootd.fnal.gov/ \
  --assume-input-files-exist \
  ./runMINIAOD_FR.py  \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName'
