#!/bin/bash

EXPECTED_ARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: $0 JOB_NAME"
fi

farmoutAnalysisJobs $1-DYJets \
  --input-files-per-job=1 \
  --job-generates-output-name \
  --infer-cmssw-path \
  --input-file-list=DYJets_M50.txt \
  --input-dir=root://cmsxrootd.fnal.gov/ \
  --assume-input-files-exist \
  ./runMINIAOD_effi.py  \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName'
