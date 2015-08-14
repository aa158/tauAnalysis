#!/bin/bash

EXPECTED_ARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: $0 JOB_NAME"
fi

farmoutAnalysisJobs $1-WJets_orig \
  --input-files-per-job=5 \
  --job-generates-output-name \
  --infer-cmssw-path \
  --input-file-list=WJets.txt \
  --input-dir=root://cmsxrootd.fnal.gov/ \
  --assume-input-files-exist \
  ./runMINIAOD_FR.py  \
  'inputFiles=$inputFileNames' 'outputFile=$outputFileName'
