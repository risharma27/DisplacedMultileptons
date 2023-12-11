
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
/*
  This is a driver script.
  It decides which code to run over which sample, the names
  of output files and so on.
*/


void testana(int sample=0){
  const char *hstfilename, *sumfilename, *skimfilename;
  //Declare a chain for input files.
  TChain *chain = new TChain("Events"); //"Events"
  //Declare an instance of our code class
  
  disp_ml m_selec;
  
  if(sample==0){
    //Add one file to chain. This is the input file.
    chain->Add("../inputs/DYJetsToLL_M-50.root");
    //Set Names of outputfiles
    hstfilename = "../hst_output/hst_DYToLL.root";
    sumfilename = "../sum_output/sum_DYToLL.txt";
    skimfilename = "../skim_output/skim_DYToLL.root";
    //Set some options
    m_selec.SetData(0); //MC=0, data=1
    m_selec.SetYear(2016);
  }

  if(sample==1){
    //Add one file to chain. This is the input file.
    chain->Add("../inputs/SingleElectron_2016postVFP_F/*");
    //Set Names of outputfiles
    hstfilename = "../hst_output/hst_SE.root";
    sumfilename = "../sum_output/sum_SE.txt";
    skimfilename = "../skim_output/skim_SE.root";
    //Set some options
    m_selec.SetData(1); //MC=0, data=1
    m_selec.SetYear(2016);
  }

  if(sample==2){
    //Add one file to chain. This is the input file.
    chain->Add("../inputs/SingleMuon_2016postVFP_F/*");
    //Set Names of outputfiles
    hstfilename = "../hst_output/hst_SM.root";
    sumfilename = "../sum_output/sum_SM.txt";
    skimfilename = "../skim_output/skim_SM.root";
    //Set some options
    m_selec.SetData(1); //MC=0, data=1
    m_selec.SetYear(2016);
  }

  std::cout<<"Output file :"<<skimfilename<<std::endl;
  // Set some more options.. set the output file names.
  m_selec.SetHstFileName(hstfilename);
  m_selec.SetSumFileName(sumfilename);
  m_selec.SetSkimFileName(skimfilename);
  m_selec.SetVerbose(10);//set verbosity level for output.
  // Call the process function which runs the code.
  chain->Process(&m_selec);

}
