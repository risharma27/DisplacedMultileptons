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
#include "Setup/GetEvtWeight.h"
#include "Setup/CustomFunctions.h"
#include "Setup/ProduceGenCollection.h"
#include "Setup/ProduceRecoCollection.h"

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

	if(_data==0){
    
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

      //----------------------------------------------------------------
      //Event-selection is done right after creating the object arrays.
      //evt_wt is also calculated alongwith.
      //This is done before any plotting.
      EventSelection();
      evtwt = GetEvtWeight();

      //Applying trigger to MC  
      bool single_muon = false;
      bool single_electron = false;
 
      if((int)Muon.size()>0 && Muon.at(0).v.Pt()>24)            single_muon = true;
      if((int)Electron.size()>0 && Electron.at(0).v.Pt()>27)    single_electron = true;
       
      bool triggered_events = false;
      //If the event has a single muon passing the trigger, keep it.
      if(single_muon) triggered_events=true;
      //If the event does not pass the single muon trigger then check for the single electron trigger, if it does then keep the event.
      else if(!single_muon && single_electron) triggered_events=true;    

      /*
	if(single_electron) triggered_events=true;
	else if(!single_electron && single_muon) triggered_events=true;
      */

    
      if(triggered_events){

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

	for(int i=0; i<(int)Electron.size(); i++){
	  h.mediumlep_iso[0]->Fill(Electron.at(i).reliso03);
	}
	for(int i=0; i<(int)Muon.size(); i++){
	  h.mediumlep_iso[1]->Fill(Muon.at(i).reliso03);
	}

      
    
	//##################### ANALYSIS BLOCK  ####################//

      
	if(evt_dispml){
		
	  dispml_evsel_plots(evtwt);
	
	}//evt_dispml
    
      }//triggered_events

    }//triggerRes
  
  }//GoodEvt

  return kTRUE;
}

//######################################
//        USER DEFINED FUNCTIONS
//######################################




void disp_ml::BookHistograms()
{

  h.nevt = new TH1F("nEvents", "0-nEvtTotal, 1-nEvtGood, 2-nEvtTrigger, 3-nEvtPass",5,0,5);

  h.dxy[0]   = new TH1F("mu_dxy", "mu_dxy", 1000, -50, 50);
  h.dz[0]    = new TH1F("mu_dz", "mu_dz", 1000, -50, 50);
  h.dxy[1]   = new TH1F("mu_|dxy|", "mu_|dxy|", 5000, 0, 50);
  h.dz[1]    = new TH1F("mu_|dz|", "mu_|dz|", 1000, 0, 100);
  h.ip3d[0]  = new TH1F("mu_ip3d", "mu_ip3d", 1000, 0, 100);
  h.sip3d[0] = new TH1F("mu_sip3d", "mu_sip3d", 1000, 0, 100);

  h.dxy[2]   = new TH1F("el_dxy", "el_dxy", 1000, -50, 50);
  h.dz[2]    = new TH1F("el_dz", "el_dz", 1000, -50, 50);
  h.dxy[3]   = new TH1F("el_|dxy|", "el_|dxy|", 5000, 0, 50);
  h.dz[3]    = new TH1F("el_|dz|", "el_|dz|", 1000, 0, 100);
  h.ip3d[1]  = new TH1F("el_ip3d", "el_ip3d", 1000, 0, 100);
  h.sip3d[1] = new TH1F("el_sip3d", "el_sip3d", 1000, 0, 100);

  h.dxy[4]   = new TH1F("lep_dxy", "lep_dxy", 1000, -50, 50);
  h.dz[4]    = new TH1F("lep_dz", "lep_dz", 1000, -50, 50);
  h.dxy[5]   = new TH1F("lep_|dxy|", "lep_|dxy|", 5000, 0, 50);
  h.dz[5]    = new TH1F("lep_|dz|", "lep_|dz|", 1000, 0, 100);
  h.ip3d[2]  = new TH1F("lep_ip3d", "lep_ip3d", 1000, 0, 100);
  h.sip3d[2] = new TH1F("lep_sip3d", "lep_sip3d", 1000, 0, 100);
  
  h.mediumlep_iso[0] = new TH1F("iso_mvamedium_el", "", 150, 0, 15);
  h.mediumlep_iso[1] = new TH1F("iso_medium_mu", "", 150, 0, 15);
  //h.elBitmap = new TH1F("el_bitmap", "", 10, 0, 10);

  /*
    h.evtweight[0][0] = new TH1F("2l1d_sf", "2l1d_sf", 50, 0, 5);
    h.evtweight[0][1] = new TH1F("2l1d_trigeff", "2l1d_trigeff", 10, 0, 1);
    h.evtweight[0][2] = new TH1F("2l1d_evtwt", "2l1d_evtwt", 50, 0, 5);
    h.evtweight[1][0] = new TH1F("1l2d_sf", "1l2d_sf", 100, 0, 100);
    h.evtweight[1][1] = new TH1F("1l2d_trigeff", "1l2d_trigeff", 10, 0, 1);
    h.evtweight[1][2] = new TH1F("1l2d_evtwt", "1l2d_evtwt", 100, 0, 100);
    h.evtweight[2][0] = new TH1F("3d_sf", "3d_sf", 100, 0, 100);
    h.evtweight[2][1] = new TH1F("3d_trigeff", "3d_trigeff", 10, 0, 1);
    h.evtweight[2][2] = new TH1F("3d_evtwt", "3d_evtwt", 100, 0, 100);
  */
 
  h._2LonZ[0] = new TH1F("zcr_invmass", "zcr_invmass", 200, 0, 200);
  h._2LonZ[1] = new TH1F("zcr_met", "zcr_met", 200, 0, 200);

  h._3L[0]  = new TH1F("3L_invmass_3l", "3L_invmass_3l", 500, 0, 500);
  h._3L[1]  = new TH1F("3L_invmass_l0l1", "3L_invmass_l0l1", 200, 0, 200);
  h._3L[2]  = new TH1F("3L_invmass_l1l2", "3L_invmass_l1l2", 200, 0, 200);
  h._3L[3]  = new TH1F("3L_invmass_l2l0", "3L_invmass_l2l0", 200, 0, 200);
  h._3L[4]  = new TH1F("3L_met", "3L_met", 200, 0, 200);
  h._3L[5]  = new TH1F("3L_pt0", "3L_pt0", 200, 0, 200);
  h._3L[6]  = new TH1F("3L_pt1", "3L_pt1", 200, 0, 200);
  h._3L[7]  = new TH1F("3L_pt2", "3L_pt2", 200, 0, 200);
  h._3L[8]  = new TH1F("3L_lt", "3L_lt", 500, 0, 500);
  h._3L[9]  = new TH1F("3L_njet", "3L_njet", 10, 0, 10);
  h._3L[10] = new TH1F("3L_ht", "3L_ht", 500, 0, 500);
  h._3L[11] = new TH1F("3L_st", "3L_st", 500, 0, 500);

  h.mumud[0] = new TH1F("mumud_invmass_ll", "mumud_invmass_ll", 200, 0, 200);
  h.mumud[1] = new TH1F("mumud_invmass_3l", "mumud_invmass_3l", 500, 0, 500);
  h.mumud[2] = new TH1F("mumud_met", "mumud_met", 200, 0, 200);

  h.eed[0] = new TH1F("eed_invmass_ll", "eed_invmass_ll", 200, 0, 200);
  h.eed[1] = new TH1F("eed_invmass_3l", "eed_invmass_3l", 500, 0, 500);
  h.eed[2] = new TH1F("eed_met", "eed_met", 200, 0, 200);
 
  h.nevsel = new TH1F("nEvSel", "1: 2l1d, 2: 1l2d, 3: 3d", 5,0,5);
  TString evsel_name[3] = {"2l1d_", "1l2d_", "3d_"};
  TString plotname[45] = {"met","pt_3l","imass_3l","pt0","pt1","pt2","pt_l0l1","delR_l0l1","delPhi_l0l1","delPhi_l0met","imass_l0l1","mt0","pt_l1l2","delR_l1l2","delPhi_l1l2","delPhi_l1met","imass_l1l2","mt1","pt_l2l0","delR_l2l0","delPhi_l2l0","delPhi_l2met","imass_l2l0","mt2","HT","njet","dRmin_l0j","dRmin_l1j","dRmin_l2j","l0_dxy","l0_dz","l0_ip3d","l0_sip3d","l0_reliso03","l1_dxy","l1_dz","l1_ip3d","l1_sip3d","l1_reliso03","l2_dxy","l2_dz","l2_ip3d","l2_sip3d","l2_reliso03","bjets"};
  int nbins[45] = {200,500,500,200,200,200,500,100,32,32,500,200,500,100,32,32,500,200,500,100,32,32,500,200,200,10,100,100,100,2000,2000,200,500,1500,2000,2000,200,500,1500,2000,2000,200,1000,1500,20};
  float blo[45] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-10,-10,0,0,0,-10,-10,0,0,0,-10,-10,0,0,0,0};
  float bhi[45] = {200,500,500,200,200,200,500,10,3.2,3.2,500,200,500,10,3.2,3.2,500,200,500,10,3.2,3.2,500,200,200,10,100,100,100,10,10,10,50,15.0,10,10,10,50,15.0,10,10,10,100,15.0,20};
  for(int ievsel=0; ievsel<3; ievsel++){
    TString name1 = evsel_name[ievsel] + "flavor";
    h.flavor[ievsel] = new TH1F(name1,"0:#mu#mu#mu, 1:#mu#mue, 2:#mue#mu, 3:#muee, 4:eee, 5:e#mue, 6:ee#mu, 7:e#mu#mu",10,0.0,10.0);
    for(int iplot=0; iplot<45; iplot++){      
      TString name2 = evsel_name[ievsel] + plotname[iplot];
      //cout << "Creating histogram " << name2 << " with nbins = " << nbins[iplot] << ", blo = " << blo[iplot] << ", bhi = " << bhi[iplot] << endl;
      h.dispml_h[ievsel][iplot] = new TH1F(name2,name2,nbins[iplot],blo[iplot],bhi[iplot]);
     
      //h.bb_h[icr][iplot]->Sumw2();
    
    }
  }

  int n_bins[15] = {500,200,200,100,200,200,64,20,500,64,20,500,64,20,500};
  float b_lo[15] = {0,0,-10.0,0.0,0.0,0.0,-3.2,0,0,-3.2,0,0,-3.2,0,0};
  float b_hi[15] = {500,200,10.0,10.0,200.0,200,3.2,10,500,3.2,10,500,3.2,10,500};
  TString flav_type[2] = {"e", "mu"};
  TString plotnames[15] = {"M_3l", "met", "l2_dxy", "l2_ip3d", "l2_sip3d", "mt2", "dphi_l0l1", "dR_l0l1", "M_l0l1", "dphi_l1l2", "dR_l1l2", "M_l1l2", "dphi_l2l0", "dR_l2l0", "M_l2l0"};
  int p=0;
  for(int flav=0; flav<2; flav++){
    for(int plot=0; plot<15; plot++){
      TString name = "2l1d_" + flav_type[flav] + "_" + plotnames[plot];
      h._2l1d[plot+p] = new TH1F(name,name,n_bins[plot],b_lo[plot],b_hi[plot]);
    }
    p=15;
  }

 
  //############################################################################################################################
  
}//BookHistograms()
