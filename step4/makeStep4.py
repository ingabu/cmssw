import os,shutil,datetime
import getpass
from ROOT import *

#IO directories must be full paths

#inputDir='/eos/uscms/store/user/sethzenz/fromdcache/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v19/nominal'   #nominal MC
inputDir='/eos/uscms/store/user/sethzenz/fromdcache/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v19'   #all MC
#inputDir='/eos/uscms/store/user/sethzenz/fromdcache/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_DATA_split_varsBDTsAdded_v19'   #data

#outputDir='/eos/uscms/store/user/jstupak/Vh/step4'
outputDir='/eos/uscms/store/user/lpcmbja/noreplica/jstupak/step4'

#########################################################################################################################
#Helper function for shutil.copytree

def files(dir, files):
    return [f for f in files if os.path.isfile(os.path.join(dir, f))]

#########################################################################################################################

runDir=os.getcwd()

if not os.path.isdir('condor'): os.mkdir('condor')

gROOT.ProcessLine('.x compileStep4.C')

cTime=datetime.datetime.now()
date='%i_%i_%i'%(cTime.year,cTime.month,cTime.day)

condorDir='/uscmst1b_scratch/lpc1/3DayLifetime/'+getpass.getuser()+'/condorLogs/%s/%s'%(date,inputDir.split('/')[-1])
outputDir+='/%s/%s'%(date,inputDir.split('/')[-1])

shutil.copytree(inputDir,outputDir,ignore=files)
shutil.copytree(inputDir,condorDir,ignore=files) 

os.system('voms-proxy-init -valid 168:00')
proxyPath=os.popen('voms-proxy-info -path')
proxyPath=proxyPath.readline().strip()

for directory, subDirectories, files in os.walk(inputDir):
    if files:
        relPath=directory.replace(inputDir,'')
        
        dict={'RUNDIR':runDir, 'RELPATH':relPath, 'CONDORDIR':condorDir, 'INPUTDIR':inputDir, 'PROXY':proxyPath}
        jdfName='%(CONDORDIR)s/%(RELPATH)s/job'%dict

        jdf=open(jdfName,'w')
        jdf.write(
"""x509userproxy = %(PROXY)s
universe = vanilla
Executable = %(RUNDIR)s/makeStep4.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Output = %(CONDORDIR)s/%(RELPATH)s/out
Error = %(CONDORDIR)s/%(RELPATH)s/err
Log = %(CONDORDIR)s/%(RELPATH)s/log
Notification = Never
Arguments = %(RUNDIR)s %(INPUTDIR)s/%(RELPATH)s

Queue 1"""%dict)
        jdf.close()
        
        os.chdir('%s/%s'%(outputDir,relPath))
        os.system('condor_submit '+jdfName)
                 
