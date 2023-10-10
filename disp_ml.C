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
#include "Setup/CustomFunctions.h"
#include "Setup/ProduceGenCollection.h"
#include "Setup/ProduceRecoCollection.h"

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

    } //if(_data==0)


    if(_data==1){
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
      
    } //if(_data==1)

    
    Sortpt(Muon);
    Sortpt(Electron);
    Sortpt(lightLep);
    Sortpt(promptMuon);
    Sortpt(displacedMuon);
    Sortpt(promptElectron);
    Sortpt(displacedElectron);
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

      
    for(int i=0; i<3; i++){
      myLep[i].clear();         //clearing myLep[evsel] for each evsel.
    }
    
    if(triggered_events){

      h.nevsel->Fill(0);

      //h.n_dispL->Fill(displacedlepton.size());
	
      bool _2l1d = false;
      bool _1l2d = false;
      bool _3d = false;

      //an event can pass all three selections, and all the three analysis will be orthogonal

      if((int)promptLepton.size()>1 && (int)displacedLepton.size()>0) _2l1d = true;
      if((int)promptLepton.size()>0 && (int)displacedLepton.size()>1) _1l2d = true;
      if((int)promptLepton.size()>=0 && (int)displacedLepton.size()>2) _3d = true;
          
      int evsel = -1;
      if(_2l1d){
	evsel=0;
	h.nevsel->Fill(1);
	myLep[evsel].push_back(promptLepton.at(0));
	myLep[evsel].push_back(promptLepton.at(1));
	myLep[evsel].push_back(displacedLepton.at(0));
      }
    
      if(_1l2d){
	evsel=1;
	h.nevsel->Fill(2);
	myLep[evsel].push_back(promptLepton.at(0));
	myLep[evsel].push_back(displacedLepton.at(0));
	myLep[evsel].push_back(displacedLepton.at(1));
      }
    
      if(_3d){
	evsel=2;
	h.nevsel->Fill(3);
	myLep[evsel].push_back(displacedLepton.at(0));
	myLep[evsel].push_back(displacedLepton.at(1));
	myLep[evsel].push_back(displacedLepton.at(2));
      }

    
      
      //if(evsel==-1) return 0;
      if(evsel!=-1){
	nEvtPass++;
	h.nevt->Fill(3);

	//***************************************************** Flavor Classification *****************************************************************//

	if(abs(myLep[evsel].at(0).id)==13){
	  if(abs(myLep[evsel].at(1).id)==13 && abs(myLep[evsel].at(2).id)==13)                  h.flavor[evsel]->Fill(0);       //mumumu
	  else if(abs(myLep[evsel].at(1).id)==13 && abs(myLep[evsel].at(2).id)==11)             h.flavor[evsel]->Fill(1);       //mumue
	  else if(abs(myLep[evsel].at(1).id)==11 && abs(myLep[evsel].at(2).id)==13)             h.flavor[evsel]->Fill(2);       //muemu
	  else if(abs(myLep[evsel].at(1).id)==11 && abs(myLep[evsel].at(2).id)==11)             h.flavor[evsel]->Fill(3);       //muee
	}	

	else if(abs(myLep[evsel].at(0).id)==11){
	  if(abs(myLep[evsel].at(1).id)==11 && abs(myLep[evsel].at(2).id)==11)                  h.flavor[evsel]->Fill(4);       //eee
	  else if(abs(myLep[evsel].at(1).id)==13 && abs(myLep[evsel].at(2).id)==11)             h.flavor[evsel]->Fill(5);       //emue
	  else if(abs(myLep[evsel].at(1).id)==11 && abs(myLep[evsel].at(2).id)==13)             h.flavor[evsel]->Fill(6);       //eemu
	  else if(abs(myLep[evsel].at(1).id)==13 && abs(myLep[evsel].at(2).id)==13)             h.flavor[evsel]->Fill(7);       //emumu	
	}

	//********************************************************************************************************************************************//

	h.dispml_h[evsel][0]->Fill(metpt);
	float sum_pt = 0.0;
	for(int i=0; i<(int)myLep[evsel].size(); i++){
	  sum_pt = sum_pt + myLep[evsel].at(i).v.Pt();
	}
	h.dispml_h[evsel][1]->Fill(sum_pt);
	float imass = ((myLep[evsel].at(0).v + myLep[evsel].at(1).v) + myLep[evsel].at(2).v).M();
	h.dispml_h[evsel][2]->Fill(imass);
	for(int j=3; j<6; j++){
	  h.dispml_h[evsel][j]->Fill(myLep[evsel].at(j-3).v.Pt());
	}
	float pt_ll[3], delR_ll[3], delPhi_ll[3], M_ll[3];
	for(int i=0; i<3; i++){
	  for(int j=i+1; j<3; j++){
	    int index = -1;
	    if(i==0 && j==1) index=0;
	    else if(i==1 && j==2) index=1;
	    else if(i == 0 && j == 2) index=2;
	  
	    pt_ll[index]=myLep[evsel].at(i).v.Pt()+myLep[evsel].at(j).v.Pt();
	    delR_ll[index]=myLep[evsel].at(i).v.DeltaR(myLep[evsel].at(j).v);
	    delPhi_ll[index]=delta_phi(myLep[evsel].at(i).v.Phi(), myLep[evsel].at(j).v.Phi());
	    //delPhi_ll[index]=myLep[evsel].at(i).v.DeltaPhi(myLep[evsel].at(j).v);
	    M_ll[index]=(myLep[evsel].at(i).v+myLep[evsel].at(j).v).M();
	  }
	}

	float delphi_lmet[3],transvmass[3];
	for(int k=0; k<3; k++){
	  delphi_lmet[k] = delta_phi(metphi, myLep[evsel].at(k).v.Phi());
	  transvmass[k] = transv_mass(myLep[evsel].at(k).v.Pt(), metpt, delphi_lmet[k]);
	}

	int p=6;
	for(int index=0; index<3; index++){
	  h.dispml_h[evsel][index+p]->Fill(pt_ll[index]);
	  h.dispml_h[evsel][index+p+1]->Fill(delR_ll[index]);
	  h.dispml_h[evsel][index+p+2]->Fill(delPhi_ll[index]);
	  h.dispml_h[evsel][index+p+3]->Fill(delphi_lmet[index]);
	  h.dispml_h[evsel][index+p+4]->Fill(M_ll[index]);
	  h.dispml_h[evsel][index+p+5]->Fill(transvmass[index]);
	  p=p+5;
	}

	
	float jet_pt = 0.0;
	for(int i=0; i<(int)recoJet.size(); i++){
	  jet_pt = jet_pt + recoJet.at(i).v.Pt();	
	}

	h.dispml_h[evsel][24]->Fill(jet_pt);

      

	std::pair<vector<int>, vector<float>> result = dR_matching(myLep[evsel], recoJet);
	vector<int> myLep_matchto_recoJet = result.first;
	vector<float> myLep_delRmin_recoJet = result.second;

	for(int i=0; i<(int)myLep_matchto_recoJet.size(); i++){
	  int matchind=myLep_matchto_recoJet.at(i);
	  float matchdR=myLep_delRmin_recoJet.at(i);
	  if(matchind>-1){
	    h.dispml_h[evsel][i+25]->Fill(matchdR);
	  }
	  else{
	    h.dispml_h[evsel][i+25]->Fill(100);
	  }
	}

	
	int q=28;
	for(int i=0; i<(int)myLep[evsel].size(); i++){
	  h.dispml_h[evsel][i+q]->Fill(myLep[evsel].at(i).dxy);
	  h.dispml_h[evsel][i+q+1]->Fill(myLep[evsel].at(i).dz);
	  h.dispml_h[evsel][i+q+2]->Fill(myLep[evsel].at(i).ip3d);
	  h.dispml_h[evsel][i+q+3]->Fill(myLep[evsel].at(i).sip3d);
	  h.dispml_h[evsel][i+q+4]->Fill(myLep[evsel].at(i).reliso03);
	  q=q+4;
	}
	
      
      }//evsel events

    
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
  h.nevsel = new TH1F("nEvSel", "1: 2l1d, 2: 1l2d, 3: 3d", 5,0,5);
  
  TString evsel_name[3] = {"2l1d_", "1l2d_", "3d_"};
  TString plotname[43] = {"met","pt_3l","imass_3l","pt0","pt1","pt2","pt_l0l1","delR_l0l1","delPhi_l0l1","delPhi_l0met","imass_l0l1","mt0","pt_l1l2","delR_l1l2","delPhi_l1l2","delPhi_l1met","imass_l1l2","mt1","pt_l2l0","delR_l2l0","delPhi_l2l0","delPhi_l2met","imass_l2l0","mt2","HT","dRmin_l0j","dRmin_l1j","dRmin_l2j","l0_dxy","l0_dz","l0_ip3d","l0_sip3d","l0_reliso03","l1_dxy","l1_dz","l1_ip3d","l1_sip3d","l1_reliso03","l2_dxy","l2_dz","l2_ip3d","l2_sip3d","l2_reliso03"};
  int nbins[43] = {200,500,500,200,200,200,500,100,32,32,500,200,500,100,32,32,500,200,500,100,32,32,500,200,200,100,100,100,2000,2000,200,200,15,2000,2000,200,200,15,2000,2000,200,200,15};
  float blo[43] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-10,-10,0,0,0,-10,-10,0,0,0,-10,-10,0,0,0};
  float bhi[43] = {200,500,500,200,200,200,500,10,3.2,3.2,500,200,500,10,3.2,3.2,500,200,500,10,3.2,3.2,500,200,200,100,100,100,10,10,10,10,0.15,10,10,10,10,0.15,10,10,10,10,0.15};
  for(int ievsel=0; ievsel<3; ievsel++){
    TString name1 = evsel_name[ievsel] + "flavor";
    h.flavor[ievsel] = new TH1F(name1,"0:#mu#mu#mu, 1:#mu#mue, 2:#mue#mu, 3:#muee, 4:eee, 5:e#mue, 6:ee#mu, 7:e#mu#mu",10,0,10);
    for(int iplot=0; iplot<43; iplot++){      
      TString name2 = evsel_name[ievsel] + plotname[iplot];
      //cout << "Creating histogram " << name2 << " with nbins = " << nbins[iplot] << ", blo = " << blo[iplot] << ", bhi = " << bhi[iplot] << endl;
      h.dispml_h[ievsel][iplot] = new TH1F(name2,name2,nbins[iplot],blo[iplot],bhi[iplot]);
     
      //h.bb_h[icr][iplot]->Sumw2();
    
    }
  }
  
 

 
  //############################################################################################################################
  
  
}
