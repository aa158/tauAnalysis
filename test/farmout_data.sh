#!/bin/bash

EXPECTED_ARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: $0 JOB_NAME"
fi

farmoutAnalysisJobs $1-WJets_Tune \
  --input-files-per-job=1 \
  --job-generates-output-name \
  --infer-cmssw-path \
  --input-file-list=WJets_Tune.txt \
  --input-dir=root://cmsxrootd.fnal.gov/ \
  --assume-input-files-exist \
  ./runDataMCcomp_singleMu.py \
  'inputFiles=$inputFileNames' \
  'outputFile=$outputFileName' \
  'isMC=1' \
  'crossSection=61526.7' \

farmoutAnalysisJobs $1-DYJets_M50 \
  --input-files-per-job=1 \
  --job-generates-output-name \
  --infer-cmssw-path \
  --input-file-list=DYJets_M50.txt \
  --input-dir=root://cmsxrootd.fnal.gov/ \
  --assume-input-files-exist \
  ./runDataMCcomp_singleMu.py \
  'inputFiles=$inputFileNames' \
  'outputFile=$outputFileName' \
  'isMC=1' \
  'crossSection=6025' \

farmoutAnalysisJobs $1-singleMuon \
  --input-files-per-job=1 \
  --job-generates-output-name \
  --infer-cmssw-path \
  --input-file-list=singleMuon.txt \
  --input-dir=root://cmsxrootd.fnal.gov/ \
  --assume-input-files-exist \
  ./runDataMCcomp_singleMu.py \
  'inputFiles=$inputFileNames' \
  'outputFile=$outputFileName' \
  'isMC=0' \
  'crossSection=-1' \
