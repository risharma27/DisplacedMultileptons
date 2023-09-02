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
      trigger2018 = (_year==2018 ? (_lep==1 ? *HLT_IsoMu24==1 : _lep==0 && *HLT_Ele32_WPTight_Gsf) : 1);
      trigger2017 = (_year==2017 ? (_lep==1 ? *HLT_IsoMu27==1 : _lep==0 && (*HLT_Ele32_WPTight_Gsf)) : 1);
      trigger2016 = (_year==2016 ? (_lep==1 ? (*HLT_IsoMu24==1) : _lep==0 && *HLT_Ele27_WPTight_Gsf) : 1);
      
      triggerRes = trigger2018 && trigger2017 && trigger2016;
    }

    

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

      if(passCuts && fabs(Muon_dxy[i])<0.05 && fabs(Muon_dz[i])<0.1){
	promptLepton.push_back(temp);
      }

      if(passCuts && fabs(Muon_dxy[i]<0.01)){
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

      if(passCuts && fabs(Electron_dxy[i])<0.05 && fabs(Electron_dz[i])<0.1){
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





    
    //#################################################### EVENT SELECTION #########################################################################//


    //1 prompt, 2 displaced leptons
    
    if(promptLepton.size()==1 && displacedLepton.size()==2){
      float dispLep_imass = (displacedLepton.at(0).v+displacedLepton.at(1).v).M();
      h.dispLep_invmass->Fill(dispLep_imass);
    }
    

   
  
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

  h.dispLep_invmass = new TH1F("imass_D0D1", "", 200, 0, 200);
  
  //############################################################################################################################
  
  
}

