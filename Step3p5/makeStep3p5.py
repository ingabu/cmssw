import os,shutil,datetime
import getpass
from ROOT import *

#IO directories must be full paths

inputDir='/eos/uscms/store/user/sethzenz/fromdcache/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V6_MC_varsAddedSummed_v19'   #All MC
outputDir='/eos/uscms/store/user/lpcmbja/noreplica/ssagir/step3p5'

#########################################################################################################################
#Helper function for shutil.copytree

def files(dir, files):
    return [f for f in files if os.path.isfile(os.path.join(dir, f))]

#########################################################################################################################

runDir=os.getcwd()

if not os.path.isdir('condor'): os.mkdir('condor')

gROOT.ProcessLine('.x compileStep3p5.C')

cTime=datetime.datetime.now()
date='%i_%i_%i'%(cTime.year,cTime.month,cTime.day)

condorDir='/uscmst1b_scratch/lpc1/3DayLifetime/'+getpass.getuser()+'/condorLogs/step3p5/%s/%s'%(date,inputDir.split('/')[-1])#'%s/condor/%s/%s'%(runDir,date,inputDir.split('/')[-1])
outputDir+='/%s/%s'%(date,inputDir.split('/')[-1])

shutil.copytree(inputDir,outputDir,ignore=files)
shutil.copytree(inputDir,condorDir,ignore=files)

os.system('voms-proxy-init -valid 168:00')
proxyPath=os.popen('voms-proxy-info -path')
proxyPath=proxyPath.readline().strip() 

for directory, subDirectories, files in os.walk(inputDir):
	relPath=directory.replace(inputDir,'')
	if files and relPath.startswith('/JER_'):
		if relPath.split('/')[2].startswith('WHiggs0'):relPathNom='/nominal'+relPath[len(relPath.split('/')[1])+1:-len(relPath.split('/')[-1])-1]
		else: relPathNom='/nominal'+relPath[len(relPath.split('/')[1])+1:]
		for file in files:
			if file.endswith(".root"):
				dict={'RUNDIR':runDir, 'RELPATH':relPath, 'RELPATHNOM':relPathNom, 'CONDORDIR':condorDir, 'INPUTDIR':inputDir, 'FILENAME':file, 'PROXY':proxyPath}
				jdfName='%(CONDORDIR)s/%(RELPATH)s/%(FILENAME)s.job'%dict

				jdf=open(jdfName,'w')
				jdf.write(
"""x509userproxy = %(PROXY)s
universe = vanilla
Executable = %(RUNDIR)s/makeStep3p5.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Output = %(CONDORDIR)s/%(RELPATH)s/%(FILENAME)s.out
Error = %(CONDORDIR)s/%(RELPATH)s/%(FILENAME)s.err
Log = %(CONDORDIR)s/%(RELPATH)s/%(FILENAME)s.log
Notification = Never
Arguments = %(FILENAME)s %(RUNDIR)s %(INPUTDIR)s/%(RELPATH)s %(INPUTDIR)s/%(RELPATHNOM)s

Queue 1"""%dict)
				jdf.close()
		
				os.chdir('%s/%s'%(outputDir,relPath))
				os.system('condor_submit '+jdfName)
                 
