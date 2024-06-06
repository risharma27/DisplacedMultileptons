###################################################################################################
# This takes all the Physics processes mentioned in the list and submits condor job for them
# using createCondorJob.py. The way I am reading the input files from the specific directories
# depends on how they are kept. So, adjust your json file accordingly.
# Make sure to make chnages to runana.C also, if you want to use codes from a different directory.
###################################################################################################
#Make sure to use python3 for this.

import sys
import os,argparse
import json
import time #for calculating how much time taken
start_time = time.time()

#Picking the arguments and flags:
parser=argparse.ArgumentParser()
parser.add_argument('--jobname',type=str, required=True, help='AnalysisName: Such as VLL2018_Mar1_v0')
parser.add_argument('--dryrun' ,type=bool,required=False,help='Check If everything is correct before submitting')
parser.add_argument('--test'   ,type=bool,required=False,help='Check one job')
parser.add_argument('--debug'  ,type=bool,required=False,help='For printing information', default=False)
args=parser.parse_args()

dryrun = args.dryrun #Don't submit any job
jobname= args.jobname
test   = args.test  #submit from only one folder.
debug  = args.debug #for debugging the condor-script

#Global parameters:
year = 2016
process_signal = True
skim = False
#condorsamples = ["DYJetsToLL","HTbinnedWJets","QCD_MuEnriched","QCD_EMEnriched","SingleTop","TTBar","WW","WZ","ZZ","WGamma","ZGamma","SingleElectron","SingleMuon"]
#condorsamples = ["SingleMuon", "SingleElectron", "DYJetsToLL"]
#condorsamples = ["SingleMuon","SingleElectron","DYJetsToLL","HTbinnedWJets","TTBar","QCD_MuEnriched","WGamma","ZGamma"]
#condorsamples  = ["SingleMuon", "QCD_MuEnriched", "DYJetsToLL"]
condorsamples  = ["ZZ"]
dumpdir = "/home/work/risharma/work/git/temp/test/histFiles"
script = "/home/work/risharma/work/git/DisplacedMultileptons/CondorSetup/runana.C"

#_____________________________________________________________
#
#                     DO NOT TOUCH BELOW
#_____________________________________________________________

if dryrun == False: os.system('mkdir -p '+dumpdir)

#initializing a list for summary:
list_processed=[]

if process_signal == True : jsonfile = '../InputJsons/skim_sample2016.json'
with open(jsonfile,'r') as infile:
    samplelist = json.load(infile)

list_processed=[]
print(f'\n\033[93mSubmitting condor jobs ...\033[0m')

#Creating condor jobs for each process mentioned in the list:
for item in condorsamples:
    for samplegroup, groupvalue in samplelist.items():
        if samplegroup == item: #The entry from the list must be identical to the key in the json file.

            list_processed.append(samplegroup)
            print('\n----------------------------------------------')
            print(f'\033[91m Submitting jobs for {samplegroup} \033[0m')
            print('----------------------------------------------')

            #The following loop runs for individual samples.
            #For example, samplegroup = 'DYJetsToLL', sample = 'M50' 
            #It runs 'createCondorJobs.py' for each sub-sample
            for sample, val in groupvalue.items():
                #if samplegroup.startswith('VLL'): sample = samplegroup+'_'+sample
                indir = val['samplepath']
                data = str(val['data'])
                lep = 'mu'
                flag = 'flag'
                if samplegroup == 'EGamma' or samplegroup == 'SingleElectron':
                    lep  = 'el'
                    flag = 'electron_dataset'

                arguments = f'{jobname} {indir} {dumpdir} {sample} {data} {year} {lep} {flag} {script} {skim} {debug}'
                processline = 'python3 createCondorJob.py '+arguments

                if dryrun == True : print(processline)
                else: os.system(processline)

                if test == True: break #only for one subsample

        if test==True : break #only for one sample


end_time = time.time()
time_taken = end_time-start_time

print('\n\033[93mDone!\033[0m\n')
print('------------------------------------------------------------------------------')
print('\033[93mSummary:\033[0m')
print(f'\033[93mSamples : {list_processed}\033[0m')
print(f'\033[93mjobname : {jobname}\033[0m')
print(f'\033[93mOutput directory (hst files): {dumpdir}\033[0m')
print(f'\033[93mtime taken : {time_taken:.2f} seconds\033[0m')
print('-------------------------------------------------------------------------------\n')
