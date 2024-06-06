#define disp_ml_cxx

#include "disp_ml.h"
#include <TH2.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

using namespace std;

//Including the header files
#include "Setup/EventSelection.h"
#include "Setup/dispml_evsel_plots.h"
#include "Setup/other_evsel_plots.h"
#include "Setup/CustomFunctions.h"
#include "Setup/ProduceGenCollection.h"
#include "Setup/ProduceRecoCollection.h"
#include "Setup/BookHistograms.h"

//Corrections
#include "Setup/Corrections/ApplyCorrections.h"
#include "Setup/Corrections/TriggerEfficiency.h"
#include "Setup/Corrections/ScaleFactors/ScaleFactors_2016UL_preVFP.h"
#include "Setup/Corrections/ScaleFactors/ScaleFactors_2016UL_postVFP.h"
#include "Setup/Corrections/ScaleFactors/ScaleFactors_2017UL.h"
#include "Setup/Corrections/ScaleFactors/ScaleFactors_2018UL.h"

void disp_ml::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
}

void disp_ml::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
  
  //Initialization of the counters:
 
  nEvtTotal    = 0;
  nEvtGood     = 0;
  nEvtTrigger  = 0;
  nEvtPass     = 0;
  
  //Other custom counters can be initialized here.
  
  _HstFile = new TFile(_HstFileName,"recreate");
  BookHistograms();
}

void disp_ml::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  
  _HstFile->Write();
  _HstFile->Close();

  float goodevtfrac = ((float)nEvtGood)/((float)nEvtTotal);
  float trigevtfrac = ((float)nEvtTrigger)/((float)nEvtTotal);
  float passevtfrac = ((float)nEvtPass)/((float)nEvtTotal);

  //The following lines are displayed on the root prompt.
  
  cout<<"---------------------------------------------"<<endl;
  cout<<"Summary:"<<endl;
  cout<<"nEvtTotal = "<<nEvtTotal<<endl;
  cout<<"nEvtGood = "<<nEvtGood<<" ("<<goodevtfrac*100<<" %)"<<endl;
  cout<<"nEvtTrigger = "<<nEvtTrigger<<" ("<<trigevtfrac*100<<" %)"<<endl;
  cout<<"nEvtPass = "<<nEvtPass<<" ("<<passevtfrac*100<<" %)"<<endl;
  cout<<"---------------------------------------------"<<endl;

  //The following lines are written on the sum_<process name>.txt file
  ofstream fout(_SumFileName);
  
  fout<<"Total events ran = "<<nEvtTotal<<endl;
  fout<<"Total good events  = "<<nEvtGood<<endl;
  fout<<" "<<endl;

  /*
  fout<<"2l1d event list"<<endl;
  for(int i=0; i<(int)evt_2l1d.size(); i++){
    fout<<evt_2l1d.at(i)<<" ";
  }
  fout<<"\n \n";
  fout<<"1l2d event list"<<endl;
  for(int i=0; i<(int)evt_1l2d.size(); i++){
    fout<<evt_1l2d.at(i)<<" ";
  }
  fout<<"\n \n";
  fout<<"3d event list"<<endl;
  for(int i=0; i<(int)evt_3d.size(); i++){
    fout<<evt_3d.at(i)<<" ";
  }
  fout<<"\n \n";
    
  //fout<<"2LonZ_triggeff"<<endl;
  for(int i=0; i<(int)singlemuon_triggeff.size(); i++){
    fout<<singlemuon_mu0pt.at(i)<<" "<<singlemuon_triggeff.at(i)<<endl;    
  }
  */
  
}

void disp_ml::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file. 
}

Bool_t disp_ml::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  fReader.SetLocalEntry(entry);
  if(_data == 0)
    fReader_MC.SetLocalEntry(entry);
  if(_data == 1)
    fReader_Data.SetLocalEntry(entry);
  
  //Verbosity determines the number of processed events after which the root prompt is supposed to display a status update.
  if(_verbosity==0 && nEvtTotal%1000000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;      
  else if(_verbosity>0 && nEvtTotal%1000000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;

  nEvtTotal++;         //Total number of events containing everything (including the trash events).
 
  h.nevt->Fill(0);
  
  //The following flags throws away some events based on unwanted properties (such as detector problems)
  GoodEvt2018 = (_year==2018 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  GoodEvt2017 = (_year==2017 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  GoodEvt2016 = (_year==2016 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  
  GoodEvt = GoodEvt2018 && GoodEvt2017 && GoodEvt2016;                          
  
  if(GoodEvt){
    nEvtGood++;                         //Total number of events containing goodEvents. The analysis is done for these good events.

    h.nevt->Fill(1);

    triggerRes=true; //Always true for MC


    if(_data==1){
      triggerRes = false;
      bool muon_trigger = false;
      bool electron_trigger = false;
      if     (_year==2016) {muon_trigger = (*HLT_IsoMu24==1); electron_trigger = (*HLT_Ele27_WPTight_Gsf==1);}
      else if(_year==2017) {muon_trigger = (*HLT_IsoMu27==1); electron_trigger = (*HLT_Ele32_WPTight_Gsf==1);}
      else if(_year==2018) {muon_trigger = (*HLT_IsoMu24==1); electron_trigger = (*HLT_Ele27_WPTight_Gsf==1);}

      //Muons are preferrred over electrons.
      //For the electron dataset, pick up only those events which do not fire a Muon trigger.
      //Otherwise there will be overcounting.
      
      triggerRes = muon_trigger || (!muon_trigger && electron_trigger);
      //triggerRes = electron_trigger || (!electron_trigger && muon_trigger);
      if(_flag == "electron_dataset" && muon_trigger) triggerRes = false;
      
    }


    if(triggerRes){
      nEvtTrigger++; //Total number of events that pass the trigger
      h.nevt->Fill(2);

      //#################################################################//
      //                         RecoParticle Block                      //
      //#################################################################//

      recoMuon.clear();
      recoElectron.clear();
      recoLepton.clear();
      recoJet.clear();
      bJet.clear();
    
      RecoLeptonArray();
      RecoJetArray();

      Sortpt(recoMuon);
      Sortpt(recoElectron);
      Sortpt(recoLepton);
      Sortpt(recoJet);
      Sortpt(bJet);

      
      //#################################################################//
      //                         GenParticle Block                      //
      //#################################################################//

   
      Muon.clear();
      Electron.clear();
      lightLep.clear();
      promptMuon.clear();
      displacedMuon.clear();
      promptElectron.clear();
      displacedElectron.clear();
      promptLepton.clear();
      displacedLepton.clear();
      

      if(_data==0){
	genMuon.clear();
	genElectron.clear();
	genLightlep.clear();

	GenLeptonArray();

	Sortpt(genMuon);
	Sortpt(genElectron);
	Sortpt(genLightlep);
      
	/*

	//####################### MC genmatching #############################//

	std::pair<vector<int>, vector<float>> mu_result = dR_matching(recoMuon, genMuon, 0.05);
	vector<int> mu_matchto_genmu = mu_result.first;
	vector<float> mu_delRmin_genmu = mu_result.second;
   
	for(int i=0; i<(int)mu_matchto_genmu.size(); i++){
	int mumatch=mu_matchto_genmu.at(i);
	float mumatchdR=mu_delRmin_genmu.at(i);
	if(mumatch>-1 && mumatchdR<0.05){
	//int mumomid = genMuon.at(mumatch).momid;
	//h.motherID[1]->Fill(mumomid);
	//h.mu_dr[1]->Fill(mumatchdR);
	Muon.push_back(recoMuon.at(i));
	lightLep.push_back(recoMuon.at(i));
	bool promptmuon = fabs(recoMuon.at(i).dxy)<0.05 && fabs(recoMuon.at(i).dz)<0.1;
	bool displacedmuon = fabs(recoMuon.at(i).dxy)>0.05 && fabs(recoMuon.at(i).dz)<10; //the dz<10 cut for displaced muons is to reduce cosmic muons backgrounds
	  
	if(promptmuon){
	promptMuon.push_back(recoMuon.at(i));
	promptLepton.push_back(recoMuon.at(i));
	}
	else if(displacedmuon){
	displacedMuon.push_back(recoMuon.at(i));
	displacedLepton.push_back(recoMuon.at(i));
	}
	}
	}

	

	std::pair<vector<int>, vector<float>> el_result = dR_matching(recoElectron, genElectron, 0.05);
	vector<int> el_matchto_genel = el_result.first;
	vector<float> el_delRmin_genel = el_result.second;

	for(int i=0; i<(int)el_matchto_genel.size(); i++){
	int elmatch=el_matchto_genel.at(i);
	float elmatchdR=el_delRmin_genel.at(i);
	if(elmatch>-1 && elmatchdR<0.05){
	//int elmomid = genElon.at(elmatch).momid;
	//h.motherID[1]->Fill(elmomid);
	//h.el_dr[1]->Fill(elmatchdR);
	Electron.push_back(recoElectron.at(i));
	lightLep.push_back(recoElectron.at(i));
	bool promptelectron = fabs(recoElectron.at(i).dxy)<0.05 && fabs(recoElectron.at(i).dz)<0.1;
	bool displacedelectron = fabs(recoElectron.at(i).dxy)>0.05;         //no cut on dz for displaced electrons
	    
	//cout<<"el_dz"<<recoElectron.at(i).dz<<endl;
	    
	if(promptelectron){
	promptElectron.push_back(recoElectron.at(i));
	promptLepton.push_back(recoElectron.at(i));
	}
	else if(displacedelectron){
	displacedElectron.push_back(recoElectron.at(i));
	displacedLepton.push_back(recoElectron.at(i));
	}
	}
	}

	*/

      } //if(_data==0)


	// if(_data==1){
	
      for(int i=0; i<(int)recoMuon.size(); i++){
	Muon.push_back(recoMuon.at(i));
	lightLep.push_back(recoMuon.at(i));
	bool promptmuon = fabs(recoMuon.at(i).dxy)<0.05 && fabs(recoMuon.at(i).dz)<0.1;
	bool displacedmuon = fabs(recoMuon.at(i).dxy)>0.05 && fabs(recoMuon.at(i).dz)<10; 
	if(promptmuon){
	  promptMuon.push_back(recoMuon.at(i));
	  promptLepton.push_back(recoMuon.at(i));
	}
	else if(displacedmuon){
	  displacedMuon.push_back(recoMuon.at(i));
	  displacedLepton.push_back(recoMuon.at(i));
	}

      }
      
      for(int j=0; j<(int)recoElectron.size(); j++){
	Electron.push_back(recoElectron.at(j));
	lightLep.push_back(recoElectron.at(j));
	bool promptelectron = fabs(recoElectron.at(j).dxy)<0.05 && fabs(recoElectron.at(j).dz)<0.1;
	bool displacedelectron = fabs(recoElectron.at(j).dxy)>0.05;
	if(promptelectron){
	  promptElectron.push_back(recoElectron.at(j));
	  promptLepton.push_back(recoElectron.at(j));
	}
	else if(displacedelectron){
	  displacedElectron.push_back(recoElectron.at(j));
	  displacedLepton.push_back(recoElectron.at(j));
	}
      }
      
      //  } //if(_data==1)

    
      Sortpt(Muon);
      Sortpt(Electron);
      Sortpt(lightLep);
      Sortpt(promptMuon);
      Sortpt(displacedMuon);
      Sortpt(promptElectron);
      Sortpt(displacedElectron);
      Sortpt(promptLepton);
      Sortpt(displacedLepton);
       	  
      metpt  = *MET_pt;
      metphi = *MET_phi;
      
      scalefactor = 1.0;
      triggeff = 1.0;
      evtwt = 1.0; //default value
      ttcr_sf = 1.0; //default value
      
      //cout<<scalefactor<<" "<<triggeff<<" "<<evtwt<<endl;
           
      //Applying trigger to MC
      
      bool triggered_events = false;
      if(abs(recoLepton.at(0).id)==11 && recoLepton.at(0).v.Pt()>27) triggered_events = true;
      if(abs(recoLepton.at(0).id)==13 && recoLepton.at(0).v.Pt()>24) triggered_events = true;	
    
      if(triggered_events){

	//----------------------------------------------------------------
	//Event-selection is done right after creating the object arrays.
	//evt_wt is also calculated alongwith.
	//This is done before any plotting.
	
	EventSelection();

	h.nevsel->Fill(0);

	for(int i=0; i<(int)Muon.size(); i++){
	  h.dxy[0]->Fill(Muon.at(i).dxy);
	  h.dxy[1]->Fill(fabs(Muon.at(i).dxy));
	  h.dz[0]->Fill(Muon.at(i).dz);
	  h.dz[1]->Fill(fabs(Muon.at(i).dz));
	  h.ip3d[0]->Fill(Muon.at(i).ip3d);
	  h.sip3d[0]->Fill(Muon.at(i).sip3d);
	}

	for(int i=0; i<(int)Electron.size(); i++){
	  h.dxy[2]->Fill(Electron.at(i).dxy);
	  h.dxy[3]->Fill(fabs(Electron.at(i).dxy));	
	  h.dz[2]->Fill(Electron.at(i).dz);
	  h.dz[3]->Fill(fabs(Electron.at(i).dz));	
	  h.ip3d[1]->Fill(Electron.at(i).ip3d);
	  h.sip3d[1]->Fill(Electron.at(i).sip3d);
	}

	for(int i=0; i<(int)lightLep.size(); i++){
	  h.dxy[4]->Fill(lightLep.at(i).dxy);
	  h.dxy[5]->Fill(fabs(lightLep.at(i).dxy));	
	  h.dz[4]->Fill(lightLep.at(i).dz);
	  h.dz[5]->Fill(fabs(lightLep.at(i).dz));
	  h.ip3d[2]->Fill(lightLep.at(i).ip3d);
	  h.sip3d[2]->Fill(lightLep.at(i).sip3d);
	}
    
	//##################### ANALYSIS BLOCK  ####################//
	  
	dispml_evsel_plots();
	other_evsel_plots();
	
      }//triggered_events

    }//triggerRes
  
  }//GoodEvt

  return kTRUE;
}

 
