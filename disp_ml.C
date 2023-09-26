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
#include "CustomFunctions.h"
#include "ProduceGenCollection.h"
#include "ProduceRecoCollection.h"

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
    promptLepton.clear();
    displacedLepton.clear();

    if(_data==0){
      genMuon.clear();
      genElectron.clear();

      GenLeptonArray();

      Sortpt(genMuon);
      Sortpt(genElectron);

      //####################### MC genmatching #############################//
    
      std::pair<vector<int>, vector<float>> mu_result = dR_matching(recoMuon, genMuon);
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
	  if(promptmuon) promptLepton.push_back(recoMuon.at(i));
	  else if(displacedmuon) displacedLepton.push_back(recoMuon.at(i));
	}
      }

      std::pair<vector<int>, vector<float>> el_result = dR_matching(recoElectron, genElectron);
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
	  if(promptelectron) promptLepton.push_back(recoElectron.at(i));
	  else if(displacedelectron) displacedLepton.push_back(recoElectron.at(i));
	}
      }

    } //if(_data==0)


    if(_data==1){
      for(int i=0; i<(int)recoMuon.size(); i++){
	Muon.push_back(recoMuon.at(i));
	lightLep.push_back(recoMuon.at(i));
	bool promptmuon = fabs(recoMuon.at(i).dxy)<0.05 && fabs(recoMuon.at(i).dz)<0.1;
	bool displacedmuon = fabs(recoMuon.at(i).dxy)>0.05 && fabs(recoMuon.at(i).dz)<10; 
	if(promptmuon) promptLepton.push_back(recoMuon.at(i));
	else if(displacedmuon) displacedLepton.push_back(recoMuon.at(i));
      }
      
      for(int j=0; j<(int)recoElectron.size(); j++){
	Electron.push_back(recoElectron.at(j));
	lightLep.push_back(recoElectron.at(j));
	bool promptelectron = fabs(recoElectron.at(j).dxy)<0.05 && fabs(recoElectron.at(j).dz)<0.1;
	bool displacedelectron = fabs(recoElectron.at(j).dxy)>0.05;
	if(promptelectron) promptLepton.push_back(recoElectron.at(j));
	else if(displacedelectron) displacedLepton.push_back(recoElectron.at(j)); 
      }
      
    } //if(_data==1)

    
    Sortpt(Muon);
    Sortpt(Electron);
    Sortpt(lightLep);
    Sortpt(promptLepton);
    Sortpt(displacedLepton);
    
    
    //####################### ANALYSIS STARTS HERE ######################//
    
    float metpt = *MET_pt;
    float metphi = *MET_phi;

    /*
    float pv_x = *PV_x;
    float pv_y = *PV_y;
    float pv_2D = sqrt(pow(pv_x,2)+pow(pv_y,2));

    SV2D.clear();
    Delta2D.clear();
    
    for(unsigned int i=0; i<(*nSV); i++){
      float sv_dxy = SV_dxy[i];
      float sv_2D = fabs(sv_dxy);
      float delta_2D = sv_2D - pv_2D;  //this variable is the distance between the primary vertex and the zeroth secondary vertex.
      SV2D.push_back(sv_2D);
      Delta2D.push_back(delta_2D);
    }
    */

    
    //##################### EVENT SELECTION ####################//


    //Applying trigger
      
    bool single_muon = false;
    bool single_electron = false;

    
    if((int)Muon.size()>0){
      if(Muon.at(0).v.Pt()>24) single_muon = true;
    }
      
    else if((int)Electron.size()>0){
      if(Electron.at(0).v.Pt()>27) single_electron = true;
    }
    
    bool triggered_events = false;
    //If the event has a single muon passing the trigger, keep it.
    if(single_muon) triggered_events=true;
    //If the event does not pass the single muon trigger then check for the single electron trigger, if it does then keep the event.
    else if(!single_muon && single_electron) triggered_events=true;    

    myLep.clear();
    
    if(triggered_events){
      bool _2l1d = false;
      bool _1l2d = false;
      bool _3d = false;

      if((int)promptLepton.size()==2 && (int)displacedLepton.size()>0) _2l1d = true;
      else if((int)promptLepton.size()==1 && (int)displacedLepton.size()>1) _1l2d = true;
      else if((int)promptLepton.size()==0 && (int)displacedLepton.size()>2) _3d = true;
          
    int evsel = -1;
    if(_2l1d){
      evsel=0;
      myLep.push_back(promptLepton.at(0));
      myLep.push_back(promptLepton.at(1));
      myLep.push_back(displacedLepton.at(0));
    }
    
    else if(_1l2d){
      evsel=1;
      myLep.push_back(promptLepton.at(0));
      myLep.push_back(displacedLepton.at(0));
      myLep.push_back(displacedLepton.at(1));
    }
    
    else if(_3d){
      evsel=2;
      myLep.push_back(displacedLepton.at(0));
      myLep.push_back(displacedLepton.at(1));
      myLep.push_back(displacedLepton.at(2));
    }
    
      
    //if(evsel==-1) return 0;
    if(evsel!=-1){
      nEvtPass++;
      h.nevt->Fill(3);
      h.dispml_h[evsel][0]->Fill(metpt);
      float sum_pt = 0.0;
      for(int i=0; i<(int)myLep.size(); i++){
	sum_pt = sum_pt + myLep.at(i).v.Pt();
      }
      h.dispml_h[evsel][1]->Fill(sum_pt);
      float imass = ((myLep.at(0).v + myLep.at(1).v) + myLep.at(2).v).M();
      h.dispml_h[evsel][2]->Fill(imass);
      for(int j=3; j<6; j++){
	h.dispml_h[evsel][j]->Fill(myLep.at(j-3).v.Pt());
      }
      float pt_ll[3], delR_ll[3], delPhi_ll[3], M_ll[3];
      for(int i=0; i<3; i++){
	for(int j=i+1; j<3; j++){
	  int index = -1;
	  if(i==0 && j==1) index=0;
	  else if(i==1 && j==2) index=1;
	  else if(i == 0 && j == 2) index=2;
	  pt_ll[index]=myLep.at(i).v.Pt()+myLep.at(j).v.Pt();
	  delR_ll[index]=myLep.at(i).v.DeltaR(myLep.at(j).v);
	  delPhi_ll[index]=delta_phi(myLep.at(i).v.Phi(), myLep.at(j).v.Phi());
	  //delPhi_ll[index]=myLep.at(i).v.DeltaPhi(myLep.at(j).v);
	  M_ll[index]=(myLep.at(i).v+myLep.at(j).v).M();
	}
      }

      float delphi_lmet[3];
      for(int k=0; k<3; k++){
	delphi_lmet[k] = delta_phi(metphi, myLep.at(k).v.Phi());
      }
    
      for(int index=0; index<3; index++){
	int p=6;
	h.dispml_h[evsel][index+p]->Fill(pt_ll[index]);
	h.dispml_h[evsel][index+p+1]->Fill(delR_ll[index]);
	h.dispml_h[evsel][index+p+2]->Fill(delPhi_ll[index]);
	h.dispml_h[evsel][index+p+3]->Fill(delphi_lmet[index]);
	h.dispml_h[evsel][index+p+4]->Fill(M_ll[index]);
	p=p+4;
      }
    }//evsel events

    
    }//triggered_events

   


    //########### ANALYSIS ENDS HERE ##############

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

  TString evsel_name[3] = {"2l1d", "1l2d", "3d"};
  TString plotname[25] = {"met","pt_3l","imass_3l","pt0","pt1","pt2","pt_l0l1","delR_l0l1","delPhi_l0l1","delPhi_l0met","imass_l0l1","pt_l1l2","delR_l1l2","delPhi_l1l2","delPhi_l1met","imass_l1l2","pt_l2l0","delR_l2l0","delPhi_l2l0","delPhi_l2met","imass_l2l0"};
  int nbins[25] = {200,500,500,200,200,200,500,20,64,64,500,500,20,64,64,500,500,20,64,64,500};
  float blo[25] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  float bhi[25] = {200,500,500,200,200,200,500,10,3.2,3.2,500,500,10,3.2,3.2,500,500,10,3.2,3.2,500};
  for(int ievsel=0; ievsel<3; ievsel++){
    for(int iplot=0; iplot<6; iplot++){
      TString name = evsel_name[ievsel] + plotname[iplot];
      h.dispml_h[ievsel][iplot] = new TH1F(name,name,nbins[iplot],blo[iplot],bhi[iplot]);
     
      //h.bb_h[icr][iplot]->Sumw2();
    }
    for(int iplot=6; iplot<21; iplot++){
      TString name = evsel_name[ievsel] + plotname[iplot];
      h.dispml_h[ievsel][iplot] = new TH1F(name,name,nbins[iplot],blo[iplot],bhi[iplot]);
    }
  }
 

 
  //############################################################################################################################
  
  
}

