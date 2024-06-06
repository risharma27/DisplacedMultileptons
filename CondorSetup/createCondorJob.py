###############################################################################################
# This is a python script that creates a condor job for only one sample (eg. DYJetsToLL_M50)
# It uses runana.C, which is kept in the same directory as the condor scripts.
# It has to be run in a For loop to create condor jobs for all the samples
# This is done using 'submitCondorJobsInBulk.py'
###############################################################################################

# Note : full paths are denoted by all caps.

import sys
import os,argparse
from datetime import datetime

timestamp = datetime.now().strftime("%Y-%m-%d")


#input arguments:
jobname    = sys.argv[1] 
samplepath = sys.argv[2] # input path
dumpdir    = sys.argv[3]
samplename = sys.argv[4] # input name (eg: DYJetsToLL_M50)
data       = sys.argv[5] # "0" = mC, "1" = data
year       = sys.argv[6]
lep        = sys.argv[7] #"el" = electron dataset, "mu" = muon dataset
flag       = sys.argv[8] #default = "flag", if EGamma or SingleElectron "electron_dataset"
scriptname = sys.argv[9]
skim_str   = sys.argv[10]
debug_str  = sys.argv[11]
#For the old setup:
#era = 'Z'
#mc = 'wz'

#_____________________________________________________________
#
#                     DO NOT TOUCH BELOW
#_____________________________________________________________

if debug_str == "True" :   debug = True
elif debug_str == "False": debug = False
else : print("Error: Give correct argument to the debug string!")

if skim_str == "True" :   skim = True
elif skim_str == "False": skim = False
else : print("Error: Give correct argument to the debug string!")

if debug == True : print('going in the createCondorJob.py script ...')
print('debug mode = '+str(debug))
print('skim mode = '+str(skim))

outputtag = 'hst'
if skim : outputtag = 'skim'

#Setting output directories:
INDIR = samplepath
CMSDIR="/home/work/risharma/work/CMSSW_10_3_1/src"
DUMP = dumpdir
OUTDIR=DUMP+"/output/"+jobname+"/"+samplename+"_"+timestamp
SUMDIR=DUMP+"/sumoutput/"+jobname+"/"+samplename+"_"+timestamp
LOGDIR=DUMP+"/log/"+jobname+"/"+samplename+"_"+timestamp
if debug == False:
    os.system("mkdir -p "+OUTDIR)
    os.system("mkdir -p "+LOGDIR)
    os.system("mkdir -p "+SUMDIR)
    os.system("chmod -R 777 "+OUTDIR)
    os.system("chmod -R 777 "+LOGDIR)
    os.system("chmod -R 777 "+SUMDIR)


########################################################################################
# Preparing the shell script containing the root command and the arguments
# This shell script is called for each file in the input directory via the .condor file.
########################################################################################

SCRIPTDIR = os.getcwd()            #Where the condor scripts are located
os.system('cd ..')                 #go back to the main directory
MAINDIR = os.getcwd()              #path of the main directory
USER = os.environ.get('USER')

runScript=SCRIPTDIR+"/runJobsCondor_simulation.sh"
#os.system(f"[ -f +"runScript+" ] && rm -rf "+runScript) #if the file already exists, remove it
os.system("touch "+runScript)
os.system("chmod a+x "+runScript)

#arguments = r'\(\"$1\",\"$2\",\"$3\",\"$4\",\"$5\"\)' #raw string
arguments = r'\(\"$1\",\"$2\",\"$3\",\"$4\",\"$5\",\"$6\"\,\"$7\"\)'
#arguments = r'\(\"$1\",\"$2\",\"$3\",\"$4\",\"$5\",\"$6\",\"$7\"\)' #for the old setup
processline = 'root -q -b -l '+scriptname+arguments

list_of_instructions=[
    "#!/bin/bash",
    "export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch",
    "export COIN_FULL_INDIRECT_RENDERING=1",
    "echo $VO_CMS_SW_DIR",
    "source $VO_CMS_SW_DIR/cmsset_default.sh",
    "cd "+CMSDIR,
    "eval `scramv1 runtime -sh`",
    processline
]
with open(runScript, 'w') as file:
    for string in list_of_instructions:
        file.write(string + '\n')


        
##################################################################
#Preparing the condor file that goes as an argument to 'condor_q'
##################################################################

condorfile=SCRIPTDIR+"/submitJobsCondor_"+USER+"_"+samplename+"_"+timestamp+".condor"
#os.system(f"[ -f +"condorfile+" ] && rm -rf "+condorfile) #if the file already exists, remove it
os.system("touch "+condorfile)
os.system("chmod a+x "+condorfile)

#The following information goes in the beginning of the condor file.
condorfile_info=[
    "universe = vanilla",
    'Requirements = (Arch == "X86_64")',
    "Executable = "+runScript,
    "should_transfer_files = YES",
    "Notification=never",
    "",
    ""
]

#We need to add the runana.C script and the  input parameters for each files in this condor file.
#For that, we need to loop over all files in the input directory.
input_files = os.listdir(INDIR)
filecount=0
for filename in input_files:
    if filename.endswith(".root"):
        filecount+=1
        
        if INDIR.split("/")[-3].startswith('HT'):
            ofile = outputtag+"_"+INDIR.split("/")[-3]+"_"+str(filecount)+".root"
            sfile = outputtag+"_"+INDIR.split("/")[-3]+"_"+str(filecount)+".txt"

        else:
            ofile = outputtag+"_"+INDIR.split("/")[-1]+"_"+str(filecount)+".root"
            sfile = outputtag+"_"+INDIR.split("/")[-1]+"_"+str(filecount)+".txt"
        
        if debug == True :
            print('file no : '+str(filecount))
            print('Input = '+INDIR+'/'+filename)
            print('Output = '+OUTDIR+"//$(Cluster)_"+ofile+"_"+str(filecount)+".root\n")
            print('SumOutput = '+SUMDIR+"//$(Cluster)_"+sfile+"_"+str(filecount)+".txt\n")

        #The following list (for each file) needs to be added to the condor file. 
        info_for_each_file = [
            "output = "+LOGDIR+"//$(Cluster)_data_"+filename+".out",
            "error = "+LOGDIR+"//$(Cluster)_data_"+filename+".err",
            "log = "+LOGDIR+"//$(Cluster)_data_"+filename+".log",
            #"arguments = "+INDIR+"/"+filename+" "+OUTDIR+"//$(Cluster)_"+ofile+"_data.root "+str(data)+" "+str(year)+" "+str(lep),
            "arguments = "+INDIR+"/"+filename+" "+OUTDIR+"//$(Cluster)_"+ofile+" "+SUMDIR+"//$(Cluster)_"+sfile+" "+str(data)+" "+str(year)+" "+str(lep)+" "+str(flag),
            #"arguments = "+INDIR+"/"+filename+" "+OUTDIR+"//$(Cluster)_"+ofile+"_"+str(filecount)+".root "+str(data)+" "+str(year)+" "+str(lep)+" "+str(era)+" "+str(mc),#for the old setup
            "queue",
            ""
        ]
        condorfile_info += info_for_each_file
        if debug == True : break

#Let's write this list to the condor file:
with open(condorfile, 'w') as file:
    for string in condorfile_info:
        file.write(string + '\n')


        
#####################
# Submitting the job
#####################

#Once the condor file is written properly, the following line will submit the job.
condorfile = "submitJobsCondor_"+USER+"_"+samplename+"_"+timestamp+".condor"
main_processline = "condor_submit "+condorfile

print("\n______"+samplename+"______")
print("Input directory : "+INDIR)
print("Output directory: "+OUTDIR)
print("Summary directory: "+SUMDIR)
print("Script :"+scriptname)

if debug == True : print('main processline = '+main_processline)
elif debug == False : os.system("eval "+main_processline)
else : print("Error: Give correct argument to the debug string!")
    
######################
#cleaning up the mess
######################

#All the condor files are dumped in the following folder after they are used.
dumpjobs = "previous_jobs/"+timestamp+"/"+samplename
print('condor dump : '+dumpjobs+"/"+condorfile)
os.system("mkdir -p "+dumpjobs)
os.system("mv "+condorfile+" "+dumpjobs+"/.")
if debug == True : print('... done with the condor script.\n')
