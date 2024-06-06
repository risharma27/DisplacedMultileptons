
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
/*
This is a driver script.
It decides which code to run over which sample, the names
of output files and so on.
*/


void ana(int sample=0){
  const char *hstfilename, *sumfilename;
  //Declare a chain for input files.
  TChain *chain = new TChain("Events"); //"Events"
  //Declare an instance of our code class
  
  disp_ml m_selec;
  //SingleTopAna m_selec;
  
  if(sample==0){
    //Add one file to chain. This is the input file.
    chain->Add("inputs/DYJetsToLL_M-50.root");
    //Set Names of outputfiles
    hstfilename = "hst_output/hst_DYToLL.root";
    sumfilename = "sum_output/sum_DYToLL.txt";
    //Set some options
    m_selec.SetData(0); //MC=0, data=1
    m_selec.SetYear(2016);
  }

  if(sample==1){
    //Add one file to chain. This is the input file.
    chain->Add("inputs/SingleElectron_2016postVFP_F/*");
    //Set Names of outputfiles
    hstfilename = "hst_output/SE_2016postVFP_F.root";
    sumfilename = "sum_output/SE_2016postVFP_F.txt";
    //Set some options
    m_selec.SetData(1); //MC=0, data=1
    m_selec.SetYear(2016);
    m_selec.SetFlag("electron_dataset");
  }

   if(sample==2){
    //Add one file to chain. This is the input file.
    chain->Add("inputs/SingleMuon_2016postVFP_F/*");
    //Set Names of outputfiles
    hstfilename = "hst_output/SM_2016postVFP_F.root";
    sumfilename = "sum_output/SM_2016postVFP_F.txt";
    //Set some options
    m_selec.SetData(1); //MC=0, data=1
    m_selec.SetYear(2016);
  }

   if(sample==3){
     //Add one file to chain. This is the input file.
     chain->Add("../Skimmer/skim_output/skim_SE.root");
     //Set Names of outputfiles
     hstfilename = "hst_output/hst_SEskim.root";
     sumfilename = "sum_output/sum_SEskim.txt";
     //Set some options
     m_selec.SetData(1); //MC=0, data=1
     m_selec.SetYear(2016);
   }

    if(sample==4){
     //Add one file to chain. This is the input file.
     chain->Add("/home/work/alaha1/public/RunII_ULSamples/2016/DYJetsToLL/postVFP/M50/*");
     //Set Names of outputfiles
     hstfilename = "hst_output/DY_M50_2016.root";
     sumfilename = "sum_output/DY_M50_2016.txt";
     //Set some options
     m_selec.SetData(1); //MC=0, data=1
     m_selec.SetYear(2016);
   }
   
  std::cout<<"Output files are "<<hstfilename<<" and "<<sumfilename<<std::endl;
  // Set some more options.. set the output file names.
  m_selec.SetHstFileName(hstfilename);
  m_selec.SetSumFileName(sumfilename);
  m_selec.SetVerbose(10);//set verbosity level for output.
  // Call the process function which runs the code.
  chain->Process(&m_selec);

}
