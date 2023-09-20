#define disp_ml_cxx

#include "disp_ml.h"
#include <TH2.h>
#include <TStyle.h>

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
  nEvtRan        = 0;
  nEvtTotal      = 0;
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

  //The following lines are displayed on the root prompt.
  cout<<"Total events ran = "<<nEvtRan<<endl;
  cout<<"Total good events = "<<nEvtTotal<<endl;

  //The following lines are written on the sum_<process name>.txt file
  ofstream fout(_SumFileName);
  fout<<"Total events ran = "<<nEvtRan<<endl;
  fout<<"Total good events  = "<<nEvtTotal<<endl;
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
  
  //The following flags throws away some events based on unwanted properties (such as detector problems)
  GoodEvt2018 = (_year==2018 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  GoodEvt2017 = (_year==2017 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  GoodEvt2016 = (_year==2016 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  
  GoodEvt = GoodEvt2018 && GoodEvt2017 && GoodEvt2016;
  
  nEvtRan++;                             //Total number of events containing everything (including the trash events).
  
  if(GoodEvt){
    nEvtTotal++;                         //Total number of events containing goodEvents
                                         //The analysis is done for these good events.



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

    /*
    if(_data==1){
      //trigger2018 = (_year==2018 ? (_lep==1 ? *HLT_IsoMu24==1 : _lep==0 && *HLT_Ele32_WPTight_Gsf) : 1);
      //trigger2017 = (_year==2017 ? (_lep==1 ? *HLT_IsoMu27==1 : _lep==0 && (*HLT_Ele32_WPTight_Gsf)) : 1);
      trigger2016 = (_year==2016 ? (_lep==1 ? (*HLT_IsoMu24==1) : _lep==0 && *HLT_Ele27_WPTight_Gsf) : 1);
     
      //triggerRes = trigger2018 && trigger2017 && trigger2016;
      triggerRes = trigger2016;
    }
    */
    

    //Construction of the arrays:
      
    //Muons
    
    int nmu = 0;                       
    recoMuon.clear();
    recoLepton.clear();
    promptLepton.clear();
    displacedLepton.clear();
    
    for(unsigned int i=0; i<(*nMuon); i++){
      Lepton temp;                      
      temp.v.SetPtEtaPhiM(Muon_pt[i],Muon_eta[i],Muon_phi[i],0.105); //the muon mass in GeV is 0.105
      temp.id = -13*Muon_charge[i];    //pdgID for mu- = 13, pdgID for mu+ = -13  
      temp.ind = i;
      temp.tightID = Muon_tightId[i];
      temp.mediumID = Muon_mediumId[i];
      temp.looseID = Muon_looseId[i];
      temp.reliso03 = Muon_pfRelIso03_all[i];
      temp.dxy = Muon_dxy[i];
      temp.dz = Muon_dz[i];
    
      bool passCuts = temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4;
      passCuts = passCuts && Muon_pfRelIso03_all[i]<0.15 && Muon_mediumId[i];
      
      if(passCuts){
	recoMuon.push_back(temp);
	recoLepton.push_back(temp);
      }

      if(passCuts && temp.v.Pt()>25 && fabs(Muon_dxy[i])<0.05 && fabs(Muon_dz[i])<0.1){
	promptLepton.push_back(temp);
      }

      if(passCuts && fabs(Muon_dxy[i])>0.01){
	displacedLepton.push_back(temp);
      }
    }
    
    Sortpt(recoMuon);

    genmatchedMu.clear();
    genmatchedLep.clear();

    if(_data==0){
      genMuon.clear();
      for(unsigned int i=0; i<(*nGenPart); i++){
	Lepton temp;
	temp.status = GenPart_status[i];
	if(temp.status==1){
	  temp.v.SetPtEtaPhiM(GenPart_pt[i],GenPart_eta[i],GenPart_phi[i],GenPart_mass[i]);
	  temp.ind=i; temp.pdgid=GenPart_pdgId[i]; temp.momid=MotherID(i,GenPart_genPartIdxMother[i]);
	  bool passCut = abs(temp.pdgid)==13 && temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4;
	  	
	  if(passCut){
	    genMuon.push_back(temp);
	  }
	}
      }
      
      Sortpt(genMuon);
     

      std::pair<vector<int>, vector<float>> mu_result = dR_matching(recoMuon, genMuon);
      vector<int> mu_matchto_genmu = mu_result.first;
      vector<float> mu_delRmin_genmu = mu_result.second;

      for(int i=0; i<(int)mu_matchto_genmu.size(); i++){
	int mumatch=mu_matchto_genmu.at(i);
	float mumatchdR=mu_delRmin_genmu.at(i);
	//h.mu_dr[0]->Fill(mumatchdR);
	if(mumatch>-1 && mumatchdR<0.05){
	  int mumomid = genMuon.at(mumatch).momid;
	  //h.motherID[1]->Fill(mumomid);
	  //h.mu_dr[1]->Fill(mumatchdR);
	  genmatchedMu.push_back(recoMuon.at(i));
	  genmatchedLep.push_back(recoMuon.at(i));
	}
      }
    } //if(_data==0)

    if(_data==1){
      for(int i=0; i<(int)recoMuon.size(); i++){
	genmatchedMu.push_back(recoMuon.at(i));
	genmatchedLep.push_back(recoMuon.at(i));
      }
    }

    Sortpt(genmatchedMu);

    

    //##########################################################################################################################################################//


    //Low Pt Electrons

    lowptElectron.clear();
    for(unsigned int i=0; i<(*nLowPtElectron); i++){

      Lepton temp;
      temp.v.SetPtEtaPhiM(LowPtElectron_pt[i],LowPtElectron_eta[i],LowPtElectron_phi[i],LowPtElectron_mass[i]);
      temp.id = -11*LowPtElectron_charge[i];
      temp.ind = i;
      temp.reliso = LowPtElectron_miniPFRelIso_all[i];

      bool passCuts = fabs(temp.v.Eta())<2.4 && temp.v.Pt()>10;
     
      if(passCuts){
	lowptElectron.push_back(temp);
      }
    }


    
    //Electrons
    
    recoElectron.clear();
    
    for(unsigned int i=0; i<(*nElectron); i++){
      
      Lepton temp;
      temp.v.SetPtEtaPhiM(Electron_pt[i],Electron_eta[i],Electron_phi[i],Electron_mass[i]);
      temp.id = -11*Electron_charge[i];
      temp.ind = i; 
      temp.reliso03 = Electron_pfRelIso03_all[i]; //h.el_reliso03->Fill(temp.reliso03);
      temp.cutBased = Electron_cutBased[i];
      temp.dxy = Electron_dxy[i];
      temp.dz = Electron_dz[i];
   
      bool passCuts = fabs(temp.v.Eta())<2.4 && temp.v.Pt()>10;
      passCuts = passCuts && Electron_pfRelIso03_all[i]<0.15 && Electron_cutBased[i]>2;
      					
      if(passCuts){
	recoElectron.push_back(temp);
	recoLepton.push_back(temp);
      }

      if(passCuts && temp.v.Pt()>25 && fabs(Electron_dxy[i])<0.05 && fabs(Electron_dz[i])<0.1){
	promptLepton.push_back(temp);
      }

      if(passCuts && fabs(Electron_dxy[i])>0.01){
	displacedLepton.push_back(temp);
      }
    }
    
    Sortpt(recoElectron);
    Sortpt(recoLepton);
    Sortpt(promptLepton);
    Sortpt(displacedLepton);
    
    genmatchedEle.clear();
    
    if(_data==0){
      genElectron.clear();  //truth electrons
      for(unsigned int i=0; i<(*nGenPart); i++){
	Lepton temp;

	/*
	  if(nEvtTotal==12){
	  cout<<i<<"  "<<GenPart_pdgId[i]<<"  "<<GenPart_genPartIdxMother[i]<<"  "<<GenPart_pdgId[GenPart_genPartIdxMother[i]]<<endl;
	  }
	*/
	
	temp.status = GenPart_status[i];
	if(temp.status==1){
	  temp.v.SetPtEtaPhiM(GenPart_pt[i],GenPart_eta[i],GenPart_phi[i],GenPart_mass[i]);
	  temp.ind=i; temp.pdgid=GenPart_pdgId[i]; temp.momid=MotherID(i,GenPart_genPartIdxMother[i]);
	  bool passCut = fabs(temp.v.Eta())<2.4 && temp.v.Pt()>10;
	  passCut = passCut && abs(temp.pdgid)==11;

	  
	  if(passCut){
	    genElectron.push_back(temp);
	  }
	}
      }
      Sortpt(genElectron);

   
      std::pair<vector<int>, vector<float>> el_result = dR_matching(recoElectron, genElectron);
      vector<int> el_matchto_genel = el_result.first;
      vector<float> el_delRmin_genel = el_result.second;

      // cout<<(int)recoElectron.size()<<" "<<(int)el_matchto_genel.size()<<" "<<el_delRmin_genel.size()<<endl;
    
      for(int i=0; i<(int)el_matchto_genel.size(); i++){
	int elmatch=el_matchto_genel.at(i);
	float elmatchdR=el_delRmin_genel.at(i);
	//h.el_dr[0]->Fill(elmatchdR);
	if(elmatch>-1 && elmatchdR<0.05){
	  int elmomid = genElectron.at(elmatch).momid;
	  //h.motherID[0]->Fill(elmomid);
	  //h.el_dr[1]->Fill(elmatchdR);
	  genmatchedEle.push_back(recoElectron.at(i));
	  genmatchedLep.push_back(recoElectron.at(i));
	}
      }
    } //if(_data==0)

    
    if(_data==1){
      for(int i=0; i<(int)recoElectron.size(); i++){
	genmatchedEle.push_back(recoElectron.at(i));
	genmatchedLep.push_back(recoElectron.at(i));
      }
    }

    Sortpt(genmatchedEle);
    Sortpt(genmatchedLep);
    

    genmatched_lowptEle.clear();
    if(_data==0){
      std::pair<vector<int>, vector<float>> lowptel_result = dR_matching(lowptElectron, genElectron);
      vector<int> lowptel_matchto_genel = lowptel_result.first;
      vector<float> lowptel_delRmin_genel = lowptel_result.second;

      for(int i=0; i<(int)lowptel_matchto_genel.size(); i++){
	int lowptelmatch=lowptel_matchto_genel.at(i);
	float lowptelmatchdR=lowptel_delRmin_genel.at(i);
     
	if(lowptelmatch>-1 && lowptelmatchdR<0.05){
	  genmatched_lowptEle.push_back(lowptElectron.at(i));
	}
      } 
    } //if(_data==0)


    if(_data==1){
      for(int i=0; i<(int)lowptElectron.size(); i++){
	genmatched_lowptEle.push_back(lowptElectron.at(i));
      }
    }

    Sortpt(genmatched_lowptEle);

    //##############################################################################################################################################//
    
    float metpt = *MET_pt;
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

    
    //#################################################### EVENT SELECTION #########################################################################//


    bool _2l1d = (int)promptLepton.size()==2 && (int)displacedLepton.size()>0;
    bool _1l2d = (int)promptLepton.size()==1 && (int)displacedLepton.size()>1;
    bool _3d = (int)promptLepton.size()==0 && (int)displacedLepton.size()>2;

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

    if(evsel==-1) return 0;
    else{
      h.dispml_h[evsel][0]->Fill(metpt);
      float sum_pt = 0.0;
      for(int i=0; i<(int)myLep.size(); i++){
	sum_pt = sum_pt + myLep.at(i).v.Pt();
      }
      h.dispml_h[evsel][1]->Fill(sum_pt);
      float imass = (myLep.at(0).v + myLep.at(1).v);
      imass = (imass + myLep.at(2).v).M();
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
	  delPhi_ll[index]=myLep.at(i).v.DeltaPhi(myLep.at(j).v);
	  M_ll[index]=(myLep.at(i).v+myLep.at(j).v).M();
	}
      }
    
      for(int index=0; index<3; index++){
	int p=6;
	h.dispml_h[evsel][index+p]->Fill(pt_ll[index]);
	h.dispml_h[evsel][index+p+1]->Fill(delR_ll[index]);
	h.dispml_h[evsel][index+p+2]->Fill(delPhi_ll[index]);
	h.dispml_h[evsel][index+p+3]->Fill(M_ll[index]);
	p=p+3;
      }
    }
    
    
  



    /*

    
    //1 prompt, 2 displaced leptons
    
    if((int)promptLepton.size()>0 && (int)displacedLepton.size()>1){
     
      h.met[0]->Fill(metpt);

      h.prompt_pt[0]->Fill(promptLepton.at(0).v.Pt());
      h.disp_pt[0]->Fill(displacedLepton.at(0).v.Pt());
      h.disp_pt[1]->Fill(displacedLepton.at(1).v.Pt());
      h.disp_pt[2]->Fill(displacedLepton.at(0).v.Pt()+displacedLepton.at(1).v.Pt());
      
      float imass_d0d1 = (displacedLepton.at(0).v+displacedLepton.at(1).v).M();
      h.dispLep_invmass[0]->Fill(imass_d0d1);

      float delphi_l0d0 = delta_phi(promptLepton.at(0).v.Phi(), displacedLepton.at(0).v.Phi());
      float delphi_l0d1 = delta_phi(promptLepton.at(0).v.Phi(), displacedLepton.at(1).v.Phi());
      float delphi_d0d1 = delta_phi(displacedLepton.at(0).v.Phi(), displacedLepton.at(1).v.Phi());

      h.delphi_ll[0]->Fill(delphi_l0d0);
      h.delphi_ll[1]->Fill(delphi_l0d1);
      h.delphi_ll[2]->Fill(delphi_d0d1);

      float delR_l0d0 = (promptLepton.at(0).v).DeltaR(displacedLepton.at(0).v);
      float delR_l0d1 = (promptLepton.at(0).v).DeltaR(displacedLepton.at(1).v);
      float delR_d0d1 = (displacedLepton.at(0).v).DeltaR(displacedLepton.at(1).v);
      h.delR_ll[0]->Fill(delR_l0d0);
      h.delR_ll[1]->Fill(delR_l0d1);
      h.delR_ll[2]->Fill(delR_d0d1);
      
      h.PV_2D[0]->Fill(pv_2D);
      if(SV2D.size()>0){
	h.SV_2D[0]->Fill(SV2D.at(0));
	h.delta2D[0]->Fill(Delta2D.at(0));
      }
    }

    

    //3 displaced leptons

    if((int)displacedLepton.size()>2){

      h.met[1]->Fill(metpt);

      h.disp_pt[3]->Fill(displacedLepton.at(0).v.Pt());
      h.disp_pt[4]->Fill(displacedLepton.at(1).v.Pt());
      h.disp_pt[5]->Fill(displacedLepton.at(2).v.Pt());
      
      float imass_d0d1d2 = ((displacedLepton.at(0).v+displacedLepton.at(1).v)+displacedLepton.at(2).v).M();
      h.dispLep_invmass[1]->Fill(imass_d0d1d2);

      float delphi_d0d1 = delta_phi(displacedLepton.at(0).v.Phi(), displacedLepton.at(1).v.Phi());
      float delphi_d0d2 = delta_phi(displacedLepton.at(0).v.Phi(), displacedLepton.at(2).v.Phi());
      float delphi_d1d2 = delta_phi(displacedLepton.at(1).v.Phi(), displacedLepton.at(2).v.Phi());

      h.delphi_ll[3]->Fill(delphi_d0d1);
      h.delphi_ll[4]->Fill(delphi_d0d2);
      h.delphi_ll[5]->Fill(delphi_d1d2);

      float delR_d0d1 = (displacedLepton.at(0).v).DeltaR(displacedLepton.at(1).v);
      float delR_d0d2 = (displacedLepton.at(0).v).DeltaR(displacedLepton.at(2).v);
      float delR_d1d2 = (displacedLepton.at(1).v).DeltaR(displacedLepton.at(2).v);

      h.delR_ll[3]->Fill(delR_d0d1);
      h.delR_ll[4]->Fill(delR_d0d2);
      h.delR_ll[5]->Fill(delR_d1d2);

      h.PV_2D[1]->Fill(pv_2D);
      if(SV2D.size()>0){
	h.SV_2D[1]->Fill(SV2D.at(0));
	h.delta2D[1]->Fill(Delta2D.at(0));
      }
    }
    
    */


    //########### ANALYSIS ENDS HERE ##############
  
  }//GoodEvt

  return kTRUE;
}

//######################################
//        USER DEFINED FUNCTIONS
//######################################


void disp_ml::Sortpt(vector<Lepton> vec)
{
  
  for(int i=0; i<(int)vec.size()-1; i++){
    for(int j=i+1; j<(int)vec.size(); j++){
      if( vec[i].v.Pt() < vec[j].v.Pt() ) swap(vec.at(i), vec.at(j));
    }
  }
}

int disp_ml::MotherID(int partindex, int momindex)
{
  int parid = GenPart_pdgId[partindex];
  int momid = GenPart_pdgId[momindex];
  while(parid==momid){
    partindex=momindex;
    momindex=GenPart_genPartIdxMother[momindex];
    parid =GenPart_pdgId[partindex];
    momid = GenPart_pdgId[momindex];
  }
  return momid;
}


std::pair<vector<int>, vector<float>> disp_ml::dR_matching(vector<Lepton> vec1, vector<Lepton> vec2)
{
  float delR_min = 999; int match = -1;
  vector<int> foundMatch;
  vector<float> delRmin;
  for(int i=0; i<(int)vec1.size(); i++){
    for(int j=0; j<(int)vec2.size(); j++){
      float delR = (vec1.at(i).v).DeltaR(vec2.at(j).v);
      if(delR<delR_min){
	delR_min=delR; match=j;
      }
    }
    foundMatch.push_back(match);
    delRmin.push_back(delR_min);
  }
  
  return std::make_pair(foundMatch, delRmin);
}


float disp_ml::delta_phi(float phi1, float phi2)
{
  //The correct deltaPhi falls in the interval [0 , pi]
  phi1 = TVector2::Phi_0_2pi(phi1);
  phi2 = TVector2::Phi_0_2pi(phi2);
  float dphi = fabs(phi1 - phi2);
  if(dphi>TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
  return dphi;
}

float disp_ml::transv_mass(float E_lep, float MET, float dphi)
{
  //The inputs are the Energy of the lepton, MET and dPhi between the lepton and MET
  float mT = sqrt(2* E_lep * MET *(1-cos(dphi)));
  return mT;
}


vector<int> disp_ml::pt_binning_count(vector<Lepton> vec)
{
  float j=10.0; int count=0;
  vector<int> pt_binned_count; 
  for(int k=0; k<9; k++){
    for(int i=0; i<(int)vec.size(); i++){
      if(vec.at(i).v.Pt()>j && vec.at(i).v.Pt()<=(j+10.0)){
	count=count+1;
      }
    }
    pt_binned_count.push_back(count);
    j=j+10.0;
    count=0;
  }
  
  return pt_binned_count;
}



void disp_ml::BookHistograms()
{
  
  TString evsel_name[3] = {"2l1d", "1l2d", "3d"};
  TString plotname[7] = {"met","pt_3l","imass_3l","pt0","pt1","pt2","delR_ll","delPhi_ll","imass_ll"};
  int nbins[5] = {200,200,20,64,200};
  float blo[5] = {0,0,0,-3.2,0};
  float bhi[5] = {200,200,10,3.2,200};
  for(int ievsel=0; ievsel<3; ievsel++){
    for(int iplot=0; iplot<2; iplot++){
      TString name = evsel_name[ievsel] + plotname[iplot];
      h.dispml_h[ievsel][iplot] = new TH1F(name,name,nbins[iplot],blo[iplot],bhi[iplot]);
      h.dispml_h[ievsel][iplot] = new TH1F(name,name,nbins[iplot],blo[iplot],bhi[iplot]);
      

      //h.bb_h[icr][iplot]->Sumw2();
    }
    for(int iplot=1; iplot<7; iplot++){
      for(int k=0; k<3; k++){
	TString name = evsel_name[evsel] + plotname[iplot] + TString(Form("%d", ilep));
	h.dispml_h[ievsel][iplot] = new TH1F(name,name,nbins[iplot],blo[iplot],bhi[iplot]);
      }
    }
  }

  

  /*

  //1l2d channel
  h.met[0] = new TH1F("1l2d_metpt", "", 200, 0, 200);
  h.prompt_pt[0] = new TH1F("1l2d_pt_l0", "", 200, 0, 200);
  h.disp_pt[0] = new TH1F("1l2d_pt_d0", "", 200, 0, 200);
  h.disp_pt[1] = new TH1F("1l2d_pt_d1", "", 200, 0, 200);
  h.disp_pt[2] = new TH1F("1l2d_pt_d0d1", "", 200, 0, 200);
  h.dispLep_invmass[0] = new TH1F("1l2d_imass_d0d1", "", 200, 0, 200);
  h.delphi_ll[0] = new TH1F("1l2d_delphi_l0d0", "", 64, 0, 3.2);
  h.delphi_ll[1] = new TH1F("1l2d_delphi_l0d1", "", 64, 0, 3.2);
  h.delphi_ll[2] = new TH1F("1l2d_delphi_d0d1", "", 64, 0, 3.2);
  h.delR_ll[0] = new TH1F("1l2d_delR_l0d0", "", 20, 0, 10);
  h.delR_ll[1] = new TH1F("1l2d_delR_l0d1", "", 20, 0, 10);
  h.delR_ll[2] = new TH1F("1l2d_delR_d0d1", "", 20, 0, 10);
  h.PV_2D[0] = new TH1F("1l2d_pv2D", "", 20, 0, 10);
  h.SV_2D[0] = new TH1F("1l2d_sv2D", "", 40, 0, 20);
  h.delta2D[0] = new TH1F("1l2d_delta2D", "", 40, 0, 20);

  //3d channel
  h.met[1] = new TH1F("3d_metpt", "", 200, 0, 200);
  h.disp_pt[3] = new TH1F("3d_pt_d0", "", 200, 0, 200);
  h.disp_pt[4] = new TH1F("3d_pt_d1", "", 200, 0, 200);
  h.disp_pt[5] = new TH1F("3d_pt_d2", "", 200, 0, 200);
  h.dispLep_invmass[1] = new TH1F("3d_imass_d0d1d2", "", 200, 0, 200);
  h.delphi_ll[3] = new TH1F("3d_delphi_d0d1", "", 64, 0, 3.2);
  h.delphi_ll[4] = new TH1F("3d_delphi_d0d2", "", 64, 0, 3.2);
  h.delphi_ll[5] = new TH1F("3d_delphi_d1d2", "", 64, 0, 3.2);
  h.delR_ll[3] = new TH1F("3d_delR_d0d1", "", 20, 0, 10);
  h.delR_ll[4] = new TH1F("3d_delR_d0d2", "", 20, 0, 10);
  h.delR_ll[5] = new TH1F("3d_delR_d1d2", "", 20, 0, 10);
  h.PV_2D[1] = new TH1F("3d_pv2D", "", 20, 0, 10);
  h.SV_2D[1] = new TH1F("3d_sv2D", "", 40, 0, 20);
  h.delta2D[1] = new TH1F("3d_delta2D", "", 40, 0, 20);

  */
  //############################################################################################################################
  
  
}

