from copy import deepcopy
from glob import glob
from ROOT import *
from utils import *

class Sample:

    def __init__(self,name,sampleType=None,inputDir=None,fileIdentifier='',altName=None,channel=None,systematic=''):
        self.name=name;  self.inputDir=inputDir; self.fileIdentifier=fileIdentifier; self.altName=altName; self.channel=channel

        if name==sampleType: print "WARNING: sample name equals sample type:", sample.name

        self.setType(sampleType)

        self.systematic=systematic
        
        #what's wrong
        
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def setType(self,sampleType):
        if type(sampleType)!=type(''): return

        self.type=sampleType

        self.isData=False; self.isMC=False; self.isSignal=False; self.isBackground=False;
        if isEqual(self.type,'data'):
            self.isData=True
        elif isEqual(self.type,'signal'):
            self.isMC=True
            self.isSignal=True
        else:
            self.isMC=True
            self.isBackground=True

    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def setInputList(self,path):
        if not self.inputDir:
            print "No inputDir specified for sample:",self.name
            return

        path+='/'+self.inputDir+'/*'+self.fileIdentifier+'*.root'
        self.inputList=glob(path)

        if len(self.inputList)==0:
            print "No input files found for sample:",self.name,"matching",path
            return

    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def makeTChain(self,treeName):
        if not self.inputList:
            print "No inputList specified for sample:",self.name

        self.chain=TChain(treeName)
        for file in self.inputList:
            self.chain.Add(file)

    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    def clone(self,name):
        theClone=deepcopy(self)
        theClone.name=name
        return theClone

################################################################################################################################################

#name,sampleType,inputDir,fileIdentifier,altName,channel):

Wh_125p6_0P=Sample('Wh_125p6_0P','signal','Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v19/nominal','WHiggs0P_M-125p6','Wh (CP = 0^{+})')
Wh_125p6_0PH=Sample('Wh_125p6_0PH','signal','Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v19/nominal','WHiggs0PH_M-125p6')
Wh_125p6_0M=Sample('Wh_125p6_0M','signal','Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v19/nominal','WHiggs0M_M-125p6','Wh (CP = 0^{-})')
Wh_125p6_0Mf05ph0=Sample('Wh_125p6_0Mf05ph0','signal','Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v19/nominal','WHiggs0Mf05ph0_M-125p6','Wh (CP mixed 50/50)')

signals=[Wh_125p6_0P,Wh_125p6_0PH,Wh_125p6_0M,Wh_125p6_0Mf05ph0]

WZ=Sample('WZ','VZ','Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v19/nominal/WZ_lumiWeighted')
ZZ=Sample('ZZ','VZ','Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v19/nominal/ZZ_lumiWeighted')
WW=Sample('WW','VV','Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v19/nominal/WW_lumiWeighted')

diboson=[WZ,ZZ,WW]

WJets=Sample('WJetsInclusive','WJets','Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v19/nominal/WJets_boostedMG_lumiWeighted')

ZJets=Sample('ZJetsInclusive','ZJets','Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v19/nominal/ZJets_boostedMG_lumiWeighted')

TTbar_FullLept=Sample('TTbar_FullLept','ttbar','Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v19/nominal/TTbar_FullLept_lumiWeighted')
TTbar_SemiLept=Sample('TTbar_SemiLept','ttbar','Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v19/nominal/TTbar_SemiLept_lumiWeighted')
TTbar_Hadronic=Sample('TTbar_Hadronic','ttbar','Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v19/nominal/TTbar_Hadronic_lumiWeighted')

ttbar=[TTbar_FullLept,TTbar_SemiLept,TTbar_Hadronic]

singleTop_sChannel=Sample('singleTop_sChannel','singleTop','Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v19/nominal/SingleTop_sChannel_lumiWeighted')
singleTop_tChannel=Sample('singleTop_tChannel','singleTop','Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v19/nominal/SingleTop_tChannel_lumiWeighted')
singleTop_tW=Sample('singleTop_tW','singleTop','Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v19/nominal/SingleTop_tW_lumiWeighted')

singleTop=[singleTop_sChannel,singleTop_tChannel,singleTop_tW]

QCD_BCtoE=Sample('QCD_BCtoE','QCD','Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v19/nominal/QCD_BCtoE_80toInf_lumiWeighted',channel='el')
QCD_Inclusive=Sample('QCD_Inclusive','QCD','Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v19/nominal/QCD_Inclusive_50toInf_lumiWeighted',channel='el')   #FIXME - Is this correct?
QCD_MuEnriched=Sample('QCD_MuEnriched','QCD','Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v19/nominal/QCD_MuEnrichedPt5_80toInf_lumiWeighted',channel='mu')

QCD=[QCD_BCtoE,QCD_MuEnriched]

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

dataEl=Sample('DataEl','Data','Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_DATA_split_varsBDTsAdded_v19','SingleElectron',channel='el')
dataMu=Sample('DataMu','Data','Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_DATA_split_varsBDTsAdded_v19','SingleMu',channel='mu')

data=[dataEl,dataMu]

#-----------------------------------------------------------------------------------------------------------------------------------------------

samplesForPlotting=[Wh_125p6_0P,Wh_125p6_0M,Wh_125p6_0Mf05ph0]+diboson+[WJets,ZJets]+ttbar+singleTop+QCD+data
#samplesForPlotting=[Wh_125p6_0P,Wh_125p6_0M,Wh_125p6_0Mf05ph0]+diboson+[ZJets]+ttbar+singleTop+QCD+data   #no W+Jets

allSamples=samplesForPlotting+[Wh_125p6_0PH]

#-----------------------------------------------------------------------------------------------------------------------------------------------

WJetsHW=Sample('WJets_shapeSys','WJets','Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v19/nominal/WJets_boostedHW_lumiWeighted',systematic='WJetsShapeUp')
ttbarMCatNLO=Sample('ttbar_shapeSys','ttbar','Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v19/nominal/TTbar_mcAtNLO_lumiWeighted',systematic='ttbarShapeUp')
systematicSamples=[WJetsHW,ttbarMCatNLO]

systematics=['JECDown','JECUp','JERDown','JERUp','btagDown','btagUp','mistagDown','mistagUp']
for sample in allSamples:
    if sample.isMC and not sample.systematic:
        for systematic in systematics:
            s=sample.clone(sample.name+'_'+systematic)
            s.systematic=systematic
            s.inputDir=s.inputDir.replace('nominal',systematic.replace('Up','_up').replace('Down','_down'))
            if s.isSignal: s.inputDir+='/'+s.fileIdentifier+'_lumiWeighted'
            systematicSamples.append(s)

allSamples+=systematicSamples

