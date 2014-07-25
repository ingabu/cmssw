#!/usr/bin/python

if __name__=='__main__':
    from Sample import *
    from Plot import Plot
    from cuts import *
    from DataCard import *

import pickle

from ROOT import *
from tdrStyle import *
from sys import argv
from array import array
from numpy import linspace
import os

setTDRStyle()
gStyle.SetOptStat(False)
gROOT.SetBatch(1)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if len(argv)>1:
    outputDir=argv[1]
else:
    outputDir='plots'

if len(argv)>2:
    inputDir=argv[2]
else:
    inputDir='/eos/uscms/store/user/lpcmbja/noreplica/ssagir/step4/2014_7_3'

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#SETTINGS

DEBUG=False

fillEmptyBins=False
blind=True
applyNormSFs=True
unroll2D=True

doBDT=True

doCutTable=False
doTheta=False
makeDataCard=True

doAllSys=True
doJECSys=False
doJERSys=False
doBTagSys=False
doMisTagSys=False
doWJetsShapeSys=False
doTTbarShapeSys=False
doStatSys=False

doCuts=[
    'bdt',
    #'mjj',
    #'WLF',
    #'WHF',
    #'ttbar'
    ]

doVtypes=[
    2,
    3
    ]

doBoosts=[
    #'low',
    'med',
    'high'
    ]

offset=1; cutTTbar=-0.35; cutWJet=0.0; cutVV=0.2
BDTMin=-1+offset
BDTMax=-1+offset+8
BDTStitching="""(({3}+BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4)*(BDT_8TeV_H125Sig_NewTTbarBkg_newCuts4<{0}))+\
((2+{3}+BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4)*((BDT_8TeV_H125Sig_NewTTbarBkg_newCuts4>{0})*(BDT_8TeV_H125Sig_0b1b2bWjetsBkg_newCuts4<{1})))+\
((4+{3}+BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4)*((BDT_8TeV_H125Sig_NewTTbarBkg_newCuts4>{0})*(BDT_8TeV_H125Sig_0b1b2bWjetsBkg_newCuts4>{1})*(BDT_8TeV_H125Sig_VVBkg_newCuts4<{2})))+\
((6+{3}+BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4)*((BDT_8TeV_H125Sig_NewTTbarBkg_newCuts4>{0})*(BDT_8TeV_H125Sig_0b1b2bWjetsBkg_newCuts4>{1})*(BDT_8TeV_H125Sig_VVBkg_newCuts4>{2})))\
""".format(cutTTbar,cutWJet,cutVV,offset)

medBoostBDTBins =linspace(-1+offset, -1+offset+0.7, 15).tolist() + [-1+offset+2] + linspace(-1+offset+2+0.45,-1+offset+2+0.9, 10).tolist() +[-1+offset+4] + linspace(-1+offset+4+0.55,-1+offset+4+1, 10).tolist() + [-1+offset+6] + linspace(-1+offset+6+0.6, -1+offset+6+1.1, 11).tolist() + [BDTMax]
highBoostBDTBins=linspace(-1+offset, -1+offset+0.7, 15).tolist() + [-1+offset+2] + linspace(-1+offset+2+0.45,-1+offset+2+1, 12).tolist() +[-1+offset+4] + linspace(-1+offset+4+0.55,-1+offset+4+1.05, 11).tolist() + [-1+offset+6] + linspace(-1+offset+6+0.6, -1+offset+6+1.2, 13).tolist() + [BDTMax]

kPS=13499645.7033 #After averaging over e/mu and 0P/0M
kHO=3883419.80677 #After averaging over e/mu and 0P/0M
MELA_SMvPS='MELA_PS/(({0}*MELA_SM) + MELA_PS)'.format(kPS)  #This arbitrary factor needs to be optimized
MELA_SMvHO='MELA_HO/(({0}*MELA_SM) + MELA_HO)'.format(kHO)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

showOverflow=True

doShapeComparison=False   #FIXME - not updated

elLumi=19040
muLumi=19040

sigmaFracUnc={}
sigmaFracUnc['VZ']=0.2
sigmaFracUnc['VV']=0.2
sigmaFracUnc['WJets']=0.2
sigmaFracUnc['ZJets']=0.2
sigmaFracUnc['ttbar']=0.15
sigmaFracUnc['singleTop']=0.15
sigmaFracUnc['QCD']=0.25
lumiFracUnc=.022   #FIXME?

signalMagFrac=20

plotBackgrounds=['QCD','ZJets','WJets','singleTop','ttbar','VV','VZ']
backgroundFillColors={'QCD':ROOT.kMagenta,'ZJets':ROOT.kYellow-7,'WJets':ROOT.kGreen-3,'singleTop':ROOT.kCyan-7,'ttbar':ROOT.kBlue-7,'VV':ROOT.kGray+2,'VZ':ROOT.kRed-7}
backgroundLineColors={'QCD':ROOT.kMagenta+1,'ZJets':ROOT.kYellow-4,'WJets':ROOT.kGreen-2,'singleTop':ROOT.kCyan-3,'ttbar':ROOT.kBlue-3,'VV':ROOT.kGray+3,'VZ':ROOT.kRed-4}

treeName='tree'
                
#---------------------------------------------------------------------------------------------------------------------------------------------
"""
if doTheta or makeDataCard:
    doVtypes=[2,3]
    doBoosts=['med','high']
    doCuts=['bdt']
    doCutTable=False
    
if doCutTable:
    doTheta=False
    makeDataCard=False
"""
##################################################################################################################################################################

if __name__=='__main__':
    
    print 'Welcome to doPostProc!'
    logFile = open(outputDir + '/log.txt','w')
    logFile.close()

    if makeDataCard:
        dataCard=DataCard()
        cardFile=outputDir+'/dataCard.txt'

    yields={}
    plots=[]
    for cuts in doCuts:
        yields[cuts]={}
        for Vtype in doVtypes:
            yields[cuts][Vtype]={}
            for boost in doBoosts:
                yields[cuts][Vtype][boost]={}

    if makeDataCard:
        for Vtype in doVtypes:
            plots+=[
                #Plot(name='h_mass',distribution='H.mass',binsX=[0]+linspace(40,240,21).tolist()+[500],xTitle='m(h) [GeV]',yLog=True,cuts=cuts,Vtype=Vtype,boost='med'),  
                #Plot(name='h_mass',distribution='H.mass',binsX=[0]+linspace(60,240,19).tolist()+[500],xTitle='m(h) [GeV]',yLog=True,cuts=cuts,Vtype=Vtype,boost='high'),
                
                #Plot(name='BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg',distribution='BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4',binsX=linspace(-1,-0.76,4).tolist()+linspace(-0.72,0.04,20).tolist()+[1],xTitle='BDT',yLog=True,cuts=cuts,Vtype=Vtype,boost='med'),
                #Plot(name='BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg',distribution='BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4',binsX=linspace(-1,-0.76,4).tolist()+linspace(-0.72,0.2,24).tolist()+[1],xTitle='BDT',yLog=True,cuts=cuts,Vtype=Vtype,boost='high'),

                #Plot(name='allBDTs',distribution=BDTStitching,binsX=medBoostBDTBins,xTitle='BDT',yLog=True,cuts=cuts,Vtype=Vtype,boost='med'), 
                #Plot(name='allBDTs',distribution=BDTStitching,binsX=highBoostBDTBins,xTitle='BDT',yLog=True,cuts=cuts,Vtype=Vtype,boost='high'),
                

                Plot(name='mainBDT_v_SMvPS', distribution='BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4:('+MELA_SMvPS+')',binsX=[0]+linspace(0.1,0.8,11).tolist()+[0.85,1],xTitle='L(0^{-})/(L(0^{+})+L(0^{-}))',binsY=[-1]+linspace(-0.75, -0.15, 7).tolist()+[1],yTitle='BDT',yLog=False,cuts=cuts,Vtype=Vtype,boost='med'),
                Plot(name='mainBDT_v_SMvPS', distribution='BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4:('+MELA_SMvPS+')',binsX=[0]+linspace(0.1,0.8,11).tolist()+[0.85,1],xTitle='L(0^{-})/(L(0^{+})+L(0^{-}))',binsY=[-1]+linspace(-0.75, -0.15, 7).tolist()+[1],yTitle='BDT',yLog=False,cuts=cuts,Vtype=Vtype,boost='high'),

                #nominal prime
                #Plot(name='mainBDT_v_VstarMass', distribution='BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4:x_mVH',binsX=[0]+linspace(325,450,6).tolist()+[550,1200],xTitle='m(Vh) [GeV]',binsY=[-1]+linspace(-0.75, -0.15, 7).tolist()+[1],yTitle='BDT',yLog=False,cuts=cuts,Vtype=Vtype,boost='med'),
                #Plot(name='mainBDT_v_VstarMass',distribution='BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4:x_mVH',binsX=[0]+linspace(400,900,11).tolist()+[1050,1200],xTitle='m(Vh) [GeV]',binsY=[-1]+linspace(-0.75, 0.05, 9).tolist()+[1],yTitle='BDT',yLog=False,cuts=cuts,Vtype=Vtype,boost='high'),

                #nominal
                #Plot(name='mainBDT_v_VstarMass', distribution='BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4:x_mVH',binsX=[0]+linspace(266+2./3,733+1./3,8).tolist()+[1200],xTitle='m(Vh) [GeV]',binsY=[-1]+linspace(-0.76, -0.2, 8).tolist()+[1],yTitle='BDT',yLog=False,cuts=cuts,Vtype=Vtype,boost='med'),
                #Plot(name='mainBDT_v_VstarMass',distribution='BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4:x_mVH',binsX=[0]+linspace(400,933+1./3,9).tolist()+[1066+2./3,1200],xTitle='m(Vh) [GeV]',binsY=[-1]+linspace(-0.76, 0.04, 11).tolist()+[1],yTitle='BDT',yLog=False,cuts=cuts,Vtype=Vtype,boost='high'),   

                #larger floor
                #Plot(name='mainBDT_v_VstarMass', distribution='BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4:x_mVH',binsX=[0]+linspace(333+1./3,666+2./3,6).tolist()+[1200],xTitle='m(Vh) [GeV]',binsY=[-1]+linspace(-0.68, -0.28, 5).tolist()+[1],yTitle='BDT',yLog=False,cuts=cuts,Vtype=Vtype,boost='med'),
                #Plot(name='mainBDT_v_VstarMass',distribution='BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4:x_mVH',binsX=[0]+linspace(466+2./3,1000,9).tolist()+[1200],xTitle='m(Vh) [GeV]',binsY=[-1]+linspace(-0.68, -0.04, 9).tolist()+[1],yTitle='BDT',yLog=False,cuts=cuts,Vtype=Vtype,boost='high'),   

                #course
                #Plot(name='mainBDT_v_VstarMass', distribution='BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4:x_mVH',binsX=[0]+linspace(266+2./3,733+1./3,4).tolist()+[1200],xTitle='m(Vh) [GeV]',binsY=[-1]+linspace(-0.76, -0.2, 4).tolist()+[1],yTitle='BDT',yLog=False,cuts=cuts,Vtype=Vtype,boost='med'),
                #Plot(name='mainBDT_v_VstarMass',distribution='BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4:x_mVH',binsX=[0]+linspace(400,933+1./3,5).tolist()+[1066+2./3,1200],xTitle='m(Vh) [GeV]',binsY=[-1]+linspace(-0.76, 0.04, 6).tolist()+[1],yTitle='BDT',yLog=False,cuts=cuts,Vtype=Vtype,boost='high'), 
                
                ]
    else:
        for Vtype in doVtypes:
            for cuts in doCuts:
                yields[cuts][Vtype]={}
                for boost in doBoosts:
                    yields[cuts][Vtype][boost]={}

                    if doCutTable:
                        plots+=[Plot(name='dummy',distribution='H.pt',nBinsX=1,xMin=0,xMax=500000,cuts=cuts,Vtype=Vtype,boost=boost)]
                    else:
                        plots+=[
                            #Plot(name='h_mass',distribution='H.mass',nBinsX=25,xMin=0,xMax=500,xTitle='m(h) [GeV]',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='x_mVH',distribution='x_mVH',nBinsX=40,xMin=0,xMax=1200,xTitle='m(VH) [GeV]',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),

                            #Plot(name='hMass_v_VstarMass',distribution='H.mass:x_mVH',binsX=[0]+range(300,700,25)+range(700,1001,100)+[1200],xTitle='m(Vh) [GeV]',binsY=[0,50]+range(75,175,10)+range(175,251,25),yTitle='m(h) [GeV]',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='hMass_v_VstarMass',distribution='H.mass:x_mVH',nBinsX=20,xMin=0,xMax=1000,xTitle='m(Vh) [GeV]',nBinsY=20,yMin=0,yMax=500,yTitle='m(h) [GeV]',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost)                    

                            #Plot(name='BDT_8TeV_H125Sig_NewTTbarBkg',distribution='BDT_8TeV_H125Sig_NewTTbarBkg_newCuts4',nBinsX=25,xMin=-1,xMax=1,xTitle='BDT',yLog=True,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='BDT_8TeV_H125Sig_0b1b2bWjetsBkg',distribution='BDT_8TeV_H125Sig_0b1b2bWjetsBkg_newCuts4',nBinsX=25,xMin=-1,xMax=1,xTitle='BDT',yLog=True,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='BDT_8TeV_H125Sig_VVBkg',distribution='BDT_8TeV_H125Sig_VVBkg_newCuts4',nBinsX=25,xMin=-1,xMax=1,xTitle='BDT',yLog=True,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg',distribution='BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4',nBinsX=25,xMin=-1,xMax=1,xTitle='BDT',yLog=True,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='allBDTs',distribution=BDTStitching,nBinsX=60,xMin=BDTMin,xMax=BDTMax,xTitle='BDT',yLog=True,cuts=cuts,Vtype=Vtype,boost=boost),

                            Plot(name='mainBDT_v_VstarMass',distribution='BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4:x_mVH',binsX=[200,400,450,500,550,600,700,850,1200],xTitle='m(Vh) [GeV]',binsY=[-1,-0.4,-0.2,0,0.1,0.2,1],yTitle='BDT',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='mainBDT_v_VstarMass',distribution='BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4:x_mVH',nBinsX=15,xMin=200,xMax=1200,xTitle='m(Vh) [GeV]',nBinsY=25,yMin=-1,yMax=1,yTitle='BDT',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='allBDTs_v_VstarMass',distribution=BDTStitching+':x_mVH',                                   nBinsX=15,xMin=200,xMax=1200,xTitle='m(Vh) [GeV]',nBinsY=40,yMin=0,yMax=BDTMax,yTitle='BDT',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='BDT_8TeV_H125Sig_NewTTbarBkg',distribution='BDT_8TeV_H125Sig_NewTTbarBkg_newCuts4',nBinsX=25,xMin=-1,xMax=1,xTitle='BDT',yLog=True,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='BDT_8TeV_H125Sig_0b1b2bWjetsBkg',distribution='BDT_8TeV_H125Sig_0b1b2bWjetsBkg_newCuts4',nBinsX=25,xMin=-1,xMax=1,xTitle='BDT',yLog=True,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='BDT_8TeV_H125Sig_VVBkg',distribution='BDT_8TeV_H125Sig_VVBkg_newCuts4',nBinsX=25,xMin=-1,xMax=1,xTitle='BDT',yLog=True,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg',distribution='BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4',nBinsX=25,xMin=-1,xMax=1,xTitle='BDT',yLog=True,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='allBDTs',distribution=BDTStitching,nBinsX=60,xMin=BDTMin,xMax=BDTMax,xTitle='BDT',yLog=True,cuts=cuts,Vtype=Vtype,boost=boost),

                            #Plot(name='MELA_SM',          distribution='MELA_SM',nBinsX=25,xMin=0,xMax=.015, xTitle='L(0^{+})',     yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='MELA_SM_smallXMax',distribution='MELA_SM',nBinsX=25,xMin=0,xMax=.003, xTitle='L(0^{+})',     yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='MELA_PS',          distribution='MELA_PS',nBinsX=25,xMin=0,xMax=50000,xTitle='L(0^{-})',     yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='MELA_HO_smallXMax',distribution='MELA_HO',nBinsX=25,xMin=0,xMax=10000,xTitle='L(0^{+}_{HO})',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='MELA_HO',          distribution='MELA_HO',nBinsX=25,xMin=0,xMax=20000,xTitle='L(0^{+}_{HO})',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),

                            #Plot(name='MELA_SMvPS',distribution=MELA_SMvPS,nBinsX=25,xMin=0,xMax=1, xTitle='L(0^{-})/(L(0^{+})+L(0^{-}))',          yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='MELA_SMvHO',distribution=MELA_SMvHO,nBinsX=25,xMin=0,xMax=1, xTitle='L(0^{+}_{HO})/(L(0^{+})+L(0^{+}_{HO}))',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),

                            #Plot(name='hMass_v_VstarMass',distribution='H.mass:x_mVH',nBinsX=25,xMin=0,xMax=1000,xTitle='m(Vh) [GeV]',nBinsY=25,yMin=0,yMax=500,yTitle='m(h) [GeV]',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='hMass_v_VstarMass',distribution='H.mass:x_mVH',binsX=[0]+range(300,700,25)+range(700,1001,100),xTitle='m(Vh) [GeV]',binsY=[0,50]+range(75,175,10)+range(175,251,25),yTitle='m(h) [GeV]',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),

                            Plot(name='hMass_v_VstarMass',distribution='H.mass:x_mVH',nBinsX=20,xMin=0,xMax=1000,xTitle='m(Vh) [GeV]',nBinsY=20,yMin=0,yMax=500,yTitle='m(h) [GeV]',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost)
                            #Plot(name='mainBDT_v_VstarMass',distribution='BDT_8TeV_H125Sig_LFHFWjetsNewTTbarVVBkg_newCuts4:x_mVH',nBinsX=15,xMin=200,xMax=1200,xTitle='m(Vh) [GeV]',nBinsY=25,yMin=-1,yMax=1,yTitle='BDT',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='allBDTs_v_VstarMass',distribution=BDTStitching+':x_mVH',                                   nBinsX=15,xMin=200,xMax=1200,xTitle='m(Vh) [GeV]',nBinsY=40,yMin=0,yMax=BDTMax,yTitle='BDT',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='allBDTs_v_MELA_SMvPS',distribution=BDTStitching+':'+MELA_SMvPS,nBinsX=15,xMin=0,xMax=1, xTitle='L(0^{-})/(L(0^{+})+L(0^{-}))',nBinsY=40,yMin=0,yMax=BDTMax,yTitle='BDT',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),

                            #Plot(name='nPVs',distribution='nPVs',nBinsX=60,xMin=0,xMax=60,xTitle='nPVs',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='h_pT',distribution='H.pt',nBinsX=25,xMin=0,xMax=500,xTitle='p_{T}(h) [GeV]',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='h_mass',distribution='H.mass',nBinsX=25,xMin=0,xMax=500,xTitle='m(h) [GeV]',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='h_eta',distribution='H.eta',nBinsX=20,xMin=-4,xMax=4,xTitle='#eta(h)',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='h_dRjj',distribution='H.dR',nBinsX=20,xMin=0,xMax=10,xTitle='#deltaR(j_{1},j_{2})',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='hJet1_ptCorr',distribution='hJet_ptCorr[0]',nBinsX=25,xMin=0,xMax=500,xTitle='p_{T}(j_{1}) [GeV]',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='hJet2_ptCorr',distribution='hJet_ptCorr[1]',nBinsX=25,xMin=0,xMax=500,xTitle='p_{T}(j_{2}) [GeV]',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='hJet1_csv',distribution='hJet_csv[0]',nBinsX=25,xMin=0,xMax=1,xTitle='csv(j_{1}) [GeV]',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='hJet2_csv',distribution='hJet_csv[1]',nBinsX=25,xMin=0,xMax=1,xTitle='csv(j_{2}) [GeV]',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='V_pT',distribution='V.pt',nBinsX=25,xMin=0,xMax=500,xTitle='p_{T}(V) [GeV]',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='METtype1corr_et',distribution='METtype1corr.et',nBinsX=25,xMin=0,xMax=500,xTitle='E_{T}^{miss} [GeV]',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='METtype1corr_sig',distribution='METtype1corr.sig',nBinsX=25,xMin=0,xMax=10,xTitle='E_{T}^{miss} significance',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='naJets',distribution='naJets',nBinsX=20,xMin=0,xMax=30,xTitle='N_{aj}',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='nalep',distribution='nalep',nBinsX=20,xMin=0,xMax=20,xTitle='N_{al}',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='hVdPhi',distribution='HVdPhi',nBinsX=20,xMin=0,xMax=3.3,xTitle='#Delta#phi(V,h)',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='lMETdPhi',distribution='lMETdPhi',nBinsX=20,xMin=-3.3,xMax=3.3,xTitle='#Delta#phi(E_{T}^{miss},l)',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='x_costheta1',distribution='x_costheta1',nBinsX=20,xMin=-1,xMax=1,xTitle='Cos(#theta_{1})',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='x_costheta2',distribution='x_costheta2',nBinsX=20,xMin=-1,xMax=1,xTitle='Cos(#theta_{2})',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='x_phi',distribution='x_phi',nBinsX=20,xMin=-3.3,xMax=3.3,xTitle='#phi',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='x_costhetastar',distribution='x_costhetastar',nBinsX=20,xMin=0,xMax=1,xTitle='Cos(#theta*)',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='x_phi1',distribution='x_phi1',nBinsX=20,xMin=-3.3,xMax=3.3,xTitle='#phi_{1}',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='x_mVh',distribution='x_mVH',nBinsX=25,xMin=0,xMax=1200,xTitle='m(Vh) [GeV]',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            #Plot(name='x_rapidityVh',distribution='x_rapidityVH',nBinsX=20,xMin=-3,xMax=3,xTitle='y(Vh)',yLog=False,cuts=cuts,Vtype=Vtype,boost=boost),
                            ]

    for plot in plots:
        y=plot.Prepare()
        try: yields[plot.cuts][plot.Vtype][plot.boost]=y
        except: pass
        plot.Draw()
        #plot.Write()
        
        if makeDataCard:
            dataCard.appendChannel('Vtype'+str(plot.Vtype)+'_'+plot.boost+'Boost',y)
    pickle.dump(yields,open(outputDir+'/yields.p','wb'))

    stdout_old = sys.stdout
    logFile = open(outputDir + '/log.txt','a')
    sys.stdout = logFile

    cWidth=15; nameWidth=30
    for Vtype in doVtypes:
        for cuts in doCuts:
            print ("VTYPE:"+str(Vtype)).ljust(nameWidth+(2*cWidth+1))
            print ("CUTS:"+cuts).ljust(nameWidth+(2*cWidth+1))
            print "".ljust(nameWidth),
            for boost in doBoosts: print boost.ljust(cWidth),
            print
            for histName in ['W_light','W_b','W_bb','ZJets','ttbar','singleTop','QCD','VZ','VV']:
                print histName.ljust(nameWidth),
                try:
                    for boost in doBoosts:
                        print str(round(yields[cuts][Vtype][boost][histName],2)).ljust(cWidth),
                except: pass
                print
            for sample in allSamples:
                #if sample.isData: continue
                print sample.name.ljust(nameWidth),
                try:
                    for boost in doBoosts:
                        print str(round(yields[cuts][Vtype][boost][sample.name],2)).ljust(cWidth),
                except: pass
                print
            for sample in ['Total Background','Data']:
                print sample.replace(' ','_').ljust(nameWidth),
                try:
                    for boost in doBoosts:
                        print str(round(yields[cuts][Vtype][boost][sample],2)).ljust(cWidth),
                    print
                except: pass
            print 'Data/Background'.ljust(nameWidth),
            for boost in doBoosts:
                try: print str(round(yields[cuts][Vtype][boost]['Data']/yields[cuts][Vtype][boost]['Total Background'],3)).ljust(cWidth),
                except: pass
            print 3*'\n'

    sys.stdout = stdout_old
    logFile.close()

    if makeDataCard:
        #sys.stdout=open(cardFile,'w')
        #dataCard.construct(sys.stdout)
        dataCard.distribution='_'.join(plots[0].name.split('_')[:-2])
        dataCard.construct()
        dataCard.toTxt(cardFile)
             



    
