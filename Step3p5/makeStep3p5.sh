#!/bin/bash

macroDir=${2}
inputDir=${3}
inputDirNom=${4}
filename=${1}

scratch=${PWD}
scramv1 project CMSSW CMSSW_5_3_6
cd CMSSW_5_3_6/src
eval `scramv1 runtime -sh`

cp -p $macroDir/*.cc .
cp -p $macroDir/*.so .
cp -p $macroDir/*.d .
export PATH=$PATH:$macroDir

XRDpath=root://cmsxrootd-site.fnal.gov//eos/uscms/store/${inputDir##/eos/uscms/store/}
XRDpathNom=root://cmsxrootd-site.fnal.gov//eos/uscms/store/${inputDirNom##/eos/uscms/store/}
root -l -b -q $macroDir/makeStep3p5.C\(\"$macroDir\",\"$XRDpath/$filename\",\"$XRDpathNom/$filename\",\"$scratch/$filename\"\)
