//#########################################################################################################
//Instead of the regular ana.C, this anaCond.C is sent to the condor job.
//It requires the arguments like input and output filenames and some other parameters.
//This driver script is run using runana.C, which uses gROOT to run this from outside the root terminal.
//#########################################################################################################

//Execution line : .x anaCond.C("/home/work/alaha1/public/RunII_ULSamples/2018/DYJetsToLL/M50/VLL_DYJetsToLL_M50_98.root", "test_outputs/test.root", "0", "2018", "mu")

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
//#include <boost/lexical_cast.hpp>// for lexical_cast()

void anaCond( TString ifname , TString ofname , TString sfname, TString data, TString year, TString lep, TString flag)
{
  //Description of the parameters:
  //1. ifname : Input file name with full path.
  //2. ofname : Output file name with full path.
  //3. data   : 0, or 1 depedning on whether 'ifname' is MC or data.
  //4. year   : Which year of data taking. (affects trigger paths and efficiencies)
  //5. lep    : If muon dataset, lep=1, if electron dataset, lep=0

  gROOT->Time();
  const char *hstfilename;
  const char *sumfilename;
  TChain *chain = new TChain("Events");
  disp_ml m_selec;
  //SingleTopAna m_selec;

  //Setting parameters from the input arguments:
  TString input = ifname;
  bool manual = input.EndsWith(".root");
  bool single=false;
  if(!manual) input += "/*.root"; //This makes sure that the input filenames always end with .root
  chain->Add(input);
  
  hstfilename = ofname;
  sumfilename = sfname;

  //SetHstFileName:
  m_selec.SetHstFileName(hstfilename);
  m_selec.SetSumFileName(sumfilename);

  //SetVerbose:
  m_selec.SetVerbose(1);

  //SetData:
  if(data=="0") m_selec.SetData(0);
  if(data=="1") m_selec.SetData(1);
 
  //SetYear:
  if(year=="2016") m_selec.SetYear(2016);
  if(year=="2017") m_selec.SetYear(2017);
  if(year=="2018") m_selec.SetYear(2018);

  //SetEra
  //m_selec.SetEra(1);

  //SetLep:
  //This parameter decides whether we are working with muon or electron dataset:
  m_selec.SetLep(1);
  if(lep=="el") m_selec.SetLep(0);
  if(lep=="mu") m_selec.SetLep(1);

  m_selec.SetFlag(flag);
  
  chain->Process(&m_selec);
  gROOT->Time();
}
