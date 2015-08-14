#!/bin/bash

EXPECTED_ARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: $0 JOB_NAME"
fi

farmoutAnalysisJobs $1-Ztautau_M500 \
  --input-files-per-job=1 \
  --job-generates-output-name \
  --infer-cmssw-path \
  --input-file-list=Ztautau_M500.txt \
  --input-dir=root://cmsxrootd.fnal.gov/ \
  --assume-input-files-exist \
  ./runMINIAOD_effi.py  \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName'

farmoutAnalysisJobs $1-Ztautau_M1500 \
  --input-files-per-job=1 \
  --job-generates-output-name \
  --infer-cmssw-path \
  --input-file-list=Ztautau_M1500.txt \
  --input-dir=root://cmsxrootd.fnal.gov/ \
  --assume-input-files-exist \
  ./runMINIAOD_effi.py  \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName'
