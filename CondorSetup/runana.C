//##########################################################################################
//This picks the appropriate libraries and run the anaCond.C script using gROOT.
//The inputs to runana can be given from the condor submission scripts.
//While anaCond.C is used for only one root file (locally), runana.C is used to run on loops.
//While anaCond.C sits at the directory where the library is, runana.C can sit anywhere.
//##########################################################################################

//test command:
//.x runana.C("/home/work/alaha1/public/RunII_ULSamples/2018/DYJetsToLL/M50/VLL_DYJetsToLL_M50_98.root", "/home/work/phazarik1/work/General-Analysis-template/AnaCodes/BasicEvtSelection/test_outputs/hst_dy_runana.root", "0", "2018", "1")

// I will add more parameters, if needed later.

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

void runana(TString ifname="inputfile" , TString ofname="outputfile" , TString sfname="sumfile" , TString data="data", TString year="year", TString lep="lep", TString flag="flag"){

  TString path_to_anaCond = "/home/work/risharma/work/git/DisplacedMultileptons";
  //TString path_to_anaCond = "/home/work/risharma/work/git/DisplacedMultileptons/Skimmer";
  TString path_to_library = path_to_anaCond;
  
  //TString anastring =".x " + path_to_anaCond+"/anaCond.C(\""+ifname+"\",\""+ofname+"\",\""+data+"\",\""+year+"\",\""+lep+"\")"; //This is the line that is being executed.
  TString anastring = ".x " + path_to_anaCond + "/anaCond.C(\"" + ifname + "\",\"" + ofname + "\",\"" + sfname + "\",\"" + data + "\",\"" + year + "\",\"" + lep + "\",\"" + flag + "\")";

  //cout<<"executing the following command:"<<endl;
  //cout<<anastring<<endl;
  
  gSystem->Load(path_to_library+"/disp_ml_C.so");
  //gSystem->Load(path_to_library+"/SingleTopAna_C.so");
  gROOT->ProcessLine(anastring);
}
