#define nano9Ana_cxx

#include "nano9Ana.h"
#include <TH2.h>
#include <TStyle.h>

void nano9Ana::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  TString option = GetOption();
}

void nano9Ana::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  
  TString option = GetOption();
  
  //Initialization of the counters:
  nEvtRan        = 0;
  nEvtTotal      = 0;

  nEvtSel = 0;
  //Other custom counters can be initialized here.
  
  _HstFile = new TFile(_HstFileName,"recreate");
  BookHistograms();
}

void nano9Ana::SlaveTerminate()
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
  fout << " " << endl;
  fout<<"2l1d event list"<<endl;
  for(unsigned int i=0; i<arr_2l1d.size(); i++){
    fout<<arr_2l1d[i]<<" ";
  } //2l1d array 
  fout <<"\n \n";
  fout<<"1l2d event list"<<endl;  
  for(unsigned int i=0; i<arr_1l2d.size(); i++){
    fout <<arr_1l2d[i]<<" ";
  } //1l2d array  
  fout <<"\n \n";
  fout<<"3d event list"<<endl;	 
  for(unsigned int i=0; i<arr_3d.size(); i++){
     fout<<arr_3d[i]<<" ";
    }//3d array
 }
  
  
void nano9Ana::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file. 
}

Bool_t nano9Ana::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
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
  if(_verbosity==0 && nEvtTotal%100000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;      
  else if(_verbosity>0 && nEvtTotal%100000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;
  
  //The following flags throws away some events based on unwanted properties (such as detector problems)
  GoodEvt2018 = (_year==2018 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  GoodEvt2017 = (_year==2017 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  GoodEvt2016 = (_year==2016 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  
   GoodEvt = GoodEvt2018 && GoodEvt2017 && GoodEvt2016;
  
  nEvtRan++;                             //Total number of events containing everything (including the trash events).
  
  if(GoodEvt){
    nEvtTotal++;     //Total number of events containing goodEvents
                                         //The analysis is done for these good events.

    triggerRes=true; //Always true for MC


    if(_data==1){
      triggerRes = false;
      bool muon_trigger = false;
      bool electron_trigger = false;
      if     (_year==2016) {muon_trigger = (*HLT_IsoMu24==1); electron_trigger = (*HLT_Ele27_WPTight_Gsf==1);}
      //else if(_year==2017) {muon_trigger = (HLT_IsoMu27==1); electron_trigger = (HLT_Ele32_WPTight_Gsf==1);}
      //else if(_year==2018) {muon_trigger = (HLT_IsoMu24==1); electron_trigger = (HLT_Ele27_WPTight_Gsf==1);}

      //Muons are preferrred over electrons.
      //For the electron dataset, pick up only those events which do not fire a Muon trigger.
      //Otherwise there will be overcounting.
      triggerRes = muon_trigger || (!muon_trigger && electron_trigger);      
    }

    if(triggerRes){
    //Construction of the arrays:
      
    //goodMu array :
    int nmu = 0;                         // This counts the number of muons in each event.
    goodMu.clear();                      // Make sure to empty the array from previous event.
    for(unsigned int i=0; i<(*nMuon); i++){
                                         // This loop runs over all the muon candidates. Some of them will pass our selection criteria.
                                         // These will be stored in the goodMu array.
      Lepton temp;                       // 'temp' is the i-th candidate.
      temp.v.SetPtEtaPhiM(Muon_pt[i],Muon_eta[i],Muon_phi[i],0.105); //the muon mass in GeV is 0.105
      temp.id = -13*Muon_charge[i];    //pdgID for mu- = 13, pdgID for mu+ = -13  
      temp.ind = i; 

      //These are the flags the 'temp' object i.e. the muon candidate has to pass.
      bool passCuts = temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4 && Muon_mediumId[i];
      //passCuts = passCuts && Muon_pfRelIso04_all[i]<0.15;
      //passCuts = passCuts && fabs(Muon_dxy[i])<0.05 && fabs(Muon_dz[i])<0.1;
      
      if(passCuts){
	goodMu.push_back(temp); // If 'temp' satisfies all the conditions, it is pushed back into goodMu
      }
    }                                  // This 'for' loop has created a goodMu array.



    //goodEle array:
    int nele = 0;
    goodEle.clear();
    for(unsigned int i=0; i<(*nElectron); i++){
      Lepton temp;
      temp.v.SetPtEtaPhiM(Electron_pt[i], Electron_eta[i], Electron_phi[i], 0.000511);
      temp.id = -11*Electron_charge[i];
      temp.ind = i;

      bool passCuts = temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4 && Electron_cutBased[i]>2;

      if(passCuts){
	goodEle.push_back(temp);
      }
    }


    
    //proLep
    int pLep=0;
    ProLep.clear();
    for(unsigned int i=0; i<(*nMuon); i++){
      Lepton temp;
      temp.v.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
      temp.id = -13*Muon_charge[i];
      temp.ind = i;

      bool PassCuts = temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4 && Muon_mediumId[i] && /*Muon_pfRelIso03_all[i]<0.15 &&*/ fabs(Muon_dxy[i])<0.05 && fabs(Muon_dz[i])<0.1;

      if(PassCuts){
	ProLep.push_back(temp);
      }
    }
    for(unsigned int i=0; i<(*nElectron); i++){
      Lepton temp;
      temp.v.SetPtEtaPhiM(Electron_pt[i], Electron_eta[i], Electron_phi[i], 0.000511);
      temp.id = -11*Electron_charge[i];
      temp.ind = i;
      
      bool PassCuts = temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4 && /*Electron_pfRelIso03_all[i]<0.15 &&*/ fabs(Electron_dxy[i])<0.05 && fabs(Electron_dz[i])<0.1 && Electron_cutBased[i]>1;
      
      if(PassCuts){
	ProLep.push_back(temp);
      }
    }


    
    //DisLep array:
    int dLep=0;
    DisLep.clear();
    for(unsigned int i=0; i<(*nMuon); i++){
      Lepton temp;
      temp.v.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
      temp.id = -13*Muon_charge[i];
      temp.ind = i;

      bool PassCuts = temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4 && Muon_mediumId[i] && /*Muon_pfRelIso03_all[i]<0.15 &&*/ fabs(Muon_dxy[i])>0.05 && fabs(Muon_dz[i])<10;

      if(PassCuts){
	DisLep.push_back(temp);
      }
    }
    for(unsigned int i=0; i<(*nElectron); i++){
      Lepton temp;
      temp.v.SetPtEtaPhiM(Electron_pt[i], Electron_eta[i], Electron_phi[i], 0.000511);
      temp.id = -11*Electron_charge[i];
      temp.ind = i;
      
      bool PassCuts = temp.v.Pt()>10 && fabs(temp.v.Eta())<2.5 && /*Electron_pfRelIso03_all[i]<0.15 &&*/ fabs(Electron_dxy[i])>0.05 && Electron_cutBased[i]>1; //&& fabs(Electron_dz[i])<0.1;
      
      if(PassCuts){
	DisLep.push_back(temp);
      }
    }

    
    //goodJet array
    int njet=0;
    goodJet.clear();
    for(unsigned int i=0; i<(*nJet); i++){
      Lepton temp;
      temp.v.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);
      temp.id = 0;
      temp.ind = i;

      bool cleanjet = true;
      for(unsigned int j=0; j<ProLep.size(); j++){
	if(temp.v.DeltaR(ProLep.at(j).v)<0.4){
	  cleanjet = false;
	  break;
	}
      }

      bool cleanjet1 = true;
      for(unsigned int j=0; j<DisLep.size(); j++){
	if(temp.v.DeltaR(DisLep.at(j).v)<0.4){
	  cleanjet1 = false;
	  break;
	}
      }

      bool PassCuts = temp.v.Pt()>30 && fabs(temp.v.Eta())<2.4 && Jet_jetId[i] >=1;

      if(PassCuts && cleanjet && cleanjet1){
	goodJet.push_back(temp);
      }
    }
    

    if(_data==0){
      //genMu coming from Z
      int gmu=0;
      genMu.clear();
      for(unsigned int i=0; i<(*nGenPart); i++){
	Lepton temp;
	temp.v.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
	temp.id = GenPart_pdgId[i];
	temp.ind = i;
	
	bool PassCuts = temp.v.Pt()>5 && fabs(temp.v.Eta())<5 && GenPart_status[i] == 1 && abs(GenPart_pdgId[i]) == 13;
	
	int momind = GenPart_genPartIdxMother[temp.ind];
	int momId = GenPart_pdgId[momind];
	while(momind == temp.id){
	  momind = GenPart_genPartIdxMother[momind];
	  momId = GenPart_pdgId[momind];
	}
	
	if(PassCuts && momId == 23){
	  genMu.push_back(temp);
	}
      } //genMu for ends

      int matIdx=-1;
      if((int)genMu.size()>0 && (int)goodMu.size()>0){
	for(unsigned int i=0; i<genMu.size(); i++){
	  for(unsigned int j=0; j<goodMu.size(); j++){
	    float dR = genMu.at(i).v.DeltaR(goodMu.at(j).v);
	    if(dR<0.05){
	      matIdx = j;
	    }
	  }
	  if(matIdx > -1){ //matIdx>-1 to consider for the case where it doesnt find any match, i.e dR!<0.05
	    h.Mu_dxy->Fill(Muon_dxy[goodMu.at(matIdx).ind]);
	    h.Mu_dz->Fill(Muon_dz[goodMu.at(matIdx).ind]);
	  }
	}
      }
      
      //genEle array
      int gEle=0;
      genEle.clear();
      for(unsigned int i=0; i<(*nGenPart); i++){
	Lepton temp;
	temp.v.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
	temp.id=GenPart_pdgId[i];
	temp.ind=i;
	
	bool PassCuts = temp.v.Pt()>5 && fabs(temp.v.Eta())<5 && GenPart_status[i]==1 && abs(GenPart_pdgId[i])==11;
	
	if(PassCuts){
	  genEle.push_back(temp);
	}
      } //genEle for loop ends

      int MatIdx=-1;
      if((int)genEle.size()>0 && (int)goodEle.size()>0){
	for(unsigned int i=0; i<genEle.size(); i++){
	  for(unsigned int j=0; j<goodEle.size(); j++){
	    float dR = genEle.at(i).v.DeltaR(goodEle.at(j).v);
	    if(dR<0.05){
	      MatIdx = j;
	    }
	  }
	  if(MatIdx > -1){ //matIdx>-1 to consider for the case where it doesnt find any match, i.e dR!<0.05
	    h.Ele_dxy->Fill(Electron_dxy[goodEle.at(MatIdx).ind]);
	    h.Ele_dz->Fill(Electron_dz[goodEle.at(MatIdx).ind]);
	  }
	}
      }
    } // if(_data==0) ends

 
    //Now we sort the goodMu in decreasing pT order
    Sortpt(goodMu);
    Sortpt(goodEle);
    Sortpt(goodJet);
    Sortpt(ProLep);
    Sortpt(DisLep);
    
    if(_data==0){
      Sortpt(genMu);
      Sortpt(genEle);
    }
    
    
    //##############
    // Analysis
    //##############

   
    //For muons
    
    h.nmuons->Fill((int)goodMu.size());   //Fill the size of goodMu array

    //Plot pT of all muons 
    for(int i=0; i<(int)goodMu.size(); i++){
      h.mupt->Fill(goodMu.at(i).v.Pt());
    }
    
    //Plotting the leading muon pT in each event.
    if((int)goodMu.size()>0){             
      h.muprop[0]->Fill(goodMu.at(0).v.Pt());
      h.muprop[1]->Fill(goodMu.at(0).v.Phi());
    }

    if((int)goodMu.size()>0){
      h.Mu0_dxy->Fill(Muon_dxy[goodMu.at(0).ind]);
      h.Mu0_dz->Fill(Muon_dz[goodMu.at(0).ind]);
      
      float a = Muon_dxy[goodMu.at(0).ind];
      float b = Muon_dz[goodMu.at(0).ind];
      float c = sqrt((a*a) + (b*b));
      h.Mu0_c->Fill(c);
    }
    
    if((int)goodMu.size()>1){
      h.Mu1_dxy->Fill(Muon_dxy[goodMu.at(1).ind]);
      h.Mu1_dz->Fill(Muon_dz[goodMu.at(1).ind]);
      
      float a = Muon_dxy[goodMu.at(1).ind];
      float b = Muon_dz[goodMu.at(1).ind];
      float c = sqrt((a*a) + (b*b));
      h.Mu1_c->Fill(c);
    }
    
    if((int)goodMu.size()>2){
      h.Mu2_dxy->Fill(Muon_dxy[goodMu.at(2).ind]);
      h.Mu2_dz->Fill(Muon_dz[goodMu.at(2).ind]);
      
      float a = Muon_dxy[goodMu.at(2).ind];
      float b = Muon_dz[goodMu.at(2).ind];
      float c = sqrt((a*a) + (b*b));
      h.Mu2_c->Fill(c);
    }
    
    if((int)goodMu.size()>3){
      h.Mu3_dxy->Fill(Muon_dxy[goodMu.at(3).ind]);
      h.Mu3_dz->Fill(Muon_dz[goodMu.at(3).ind]);
      
      float a = Muon_dxy[goodMu.at(3).ind];
      float b = Muon_dz[goodMu.at(3).ind];
      float c = sqrt((a*a) + (b*b));
      h.Mu3_c->Fill(c);
    }

    if((int)goodEle.size()>0){
      h.Ele0_dxy->Fill(Electron_dxy[goodEle.at(0).ind]);
      h.Ele0_dz->Fill(Electron_dz[goodEle.at(0).ind]);
      
      float a = Electron_dxy[goodEle.at(0).ind];
      float b = Electron_dz[goodEle.at(0).ind];
      float c = sqrt((a*a) + (b*b));
      h.Ele0_c->Fill(c);
    }
    
    if((int)goodEle.size()>1){
      h.Ele1_dxy->Fill(Electron_dxy[goodEle.at(1).ind]);
      h.Ele1_dz->Fill(Electron_dz[goodEle.at(1).ind]);
      
      float a = Electron_dxy[goodEle.at(1).ind];
      float b = Electron_dz[goodEle.at(1).ind];
      float c = sqrt((a*a) + (b*b));
      h.Ele1_c->Fill(c);
    }
    
    if((int)goodEle.size()>2){
      h.Ele2_dxy->Fill(Electron_dxy[goodEle.at(2).ind]);
      h.Ele2_dz->Fill(Electron_dz[goodEle.at(2).ind]);
      
      float a = Electron_dxy[goodEle.at(2).ind];
      float b = Electron_dz[goodEle.at(2).ind];
      float c = sqrt((a*a) + (b*b));
      h.Ele2_c->Fill(c);
    }
    
    if((int)goodEle.size()>3){
      h.Ele3_dxy->Fill(Electron_dxy[goodEle.at(3).ind]);
      h.Ele3_dz->Fill(Electron_dz[goodEle.at(3).ind]);
      
      float a = Electron_dxy[goodEle.at(3).ind];
      float b = Electron_dz[goodEle.at(3).ind];
      float c = sqrt((a*a) + (b*b));
      h.Ele3_c->Fill(c);
    }


      
    //Applying trigger
    
    bool single_muon = false;
    bool single_electron = false;
    
    if((int)goodMu.size()>0 && goodMu.at(0).v.Pt()>24)            single_muon = true;
    if((int)goodEle.size()>0 && goodEle.at(0).v.Pt()>27)    single_electron = true;
       
    bool triggered_events = false;
    //If the event has a single muon passing the trigger, keep it.
    if(single_muon) triggered_events=true;
    //If the event does not pass the single muon trigger then check for the single electron trigger, if it does then keep the event.
    else if(!single_muon && single_electron) triggered_events=true;    

    if(triggered_events){
      

      if((int)ProLep.size()>0 && (int)DisLep.size()>1){
	arr_1l2d.push_back((int)nEvtTotal);
      } //1l2d array

      if((int)ProLep.size()>1 && (int)DisLep.size()>0){
	arr_2l1d.push_back((int)nEvtTotal);
      } //2l1d array

      if((int)ProLep.size()>=0 && (int)DisLep.size()>2){
	arr_3d.push_back((int)nEvtTotal);
      } //3d array
      
      
	  
      


      //############################################################
      // ####################### 1l2d channel ############################

      

      //1l2d channel
      if((int)ProLep.size()>0 && (int)DisLep.size()>1){ //atleast 1 prompt and atleast 2 displaced

	//mass plots
	h.m0_l0d0->Fill((ProLep.at(0).v + DisLep.at(0).v).M());
	h.m0_d0d1->Fill((DisLep.at(0).v + DisLep.at(1).v).M());
	h.m0_l0d1->Fill((ProLep.at(0).v + DisLep.at(1).v).M());
	h.m0_l0d0d1->Fill((ProLep.at(0).v + DisLep.at(0).v+DisLep.at(1).v).M());
	
	//dPhi plots
	h.dPhi0_l0d0->Fill(ProLep.at(0).v.DeltaPhi(DisLep.at(0).v));
	h.dPhi0_d0d1->Fill(DisLep.at(0).v.DeltaPhi(DisLep.at(1).v));
	h.dPhi0_l0d1->Fill(ProLep.at(0).v.DeltaPhi(DisLep.at(1).v));

	//dPhi_MET_lepton
        float dPhi0_l0MET = delta_phi(ProLep.at(0).v.Pt(),*MET_pt);
	h.dPhi0_l0MET->Fill(dPhi0_l0MET);
        float dPhi0_d0MET = delta_phi(DisLep.at(0).v.Pt(),*MET_pt);
	h.dPhi0_d0MET->Fill(dPhi0_d0MET);
        float dPhi0_d1MET = delta_phi(DisLep.at(1).v.Pt(),*MET_pt);
	h.dPhi0_d1MET->Fill(dPhi0_d1MET);

	//dR plots
	h.dR0_l0d0->Fill(ProLep.at(0).v.DeltaR(DisLep.at(0).v));
	h.dR0_d0d1->Fill(DisLep.at(0).v.DeltaR(DisLep.at(1).v));
	h.dR0_l0d1->Fill(ProLep.at(0).v.DeltaR(DisLep.at(1).v));

	//pT plots
	h.pT0_l0->Fill(ProLep.at(0).v.Pt());
	h.pT0_d0->Fill(DisLep.at(0).v.Pt());
	h.pT0_d1->Fill(DisLep.at(1).v.Pt());
	h.pT0_l0d0->Fill((ProLep.at(0).v + DisLep.at(0).v).Pt());
	h.pT0_d0d1->Fill((DisLep.at(0).v + DisLep.at(1).v).Pt());
	h.pT0_l0d1->Fill((ProLep.at(0).v + DisLep.at(1).v).Pt());
	h.pT0_l0d0d1->Fill((ProLep.at(0).v + DisLep.at(0).v + DisLep.at(1).v).Pt());
       
      	//Isolation, dxy, dz, sqrt(dxy^2 + dz^2)
	if(abs(ProLep.at(0).id)==11){ //for ProLep.at(0)  //ProLep.at(0) is electron
	  h.Iso0_l0->Fill(Electron_pfRelIso03_all[ProLep.at(0).ind]);
	  h.dxy0_l0->Fill(Electron_dxy[ProLep.at(0).ind]);
	  h.dz0_l0->Fill(Electron_dz[ProLep.at(0).ind]);
	  float a = Electron_dxy[ProLep.at(0).ind];
	  float b = Electron_dz[ProLep.at(0).ind];
	  float c = sqrt((a*a)+(b*b));
	  h.c0_l0->Fill(c);
	  h.ip3d0_l0->Fill(Electron_ip3d[ProLep.at(0).ind]);
	  h.sip3d0_l0->Fill(Electron_sip3d[ProLep.at(0).ind]);
	}
	else if(abs(ProLep.at(0).id)==13){  //ProLep.at(0) is muon
	  h.Iso0_l0->Fill(Muon_pfRelIso03_all[ProLep.at(0).ind]);
	  h.dxy0_l0->Fill(Muon_dxy[ProLep.at(0).ind]);
	  h.dz0_l0->Fill(Muon_dz[ProLep.at(0).ind]);
	  float a = Muon_dxy[ProLep.at(0).ind];
	  float b = Muon_dz[ProLep.at(0).ind];
	  float c = sqrt((a*a)+(b*b));
	  h.c0_l0->Fill(c);
	  h.ip3d0_l0->Fill(Muon_ip3d[ProLep.at(0).ind]);
	  h.sip3d0_l0->Fill(Muon_sip3d[ProLep.at(0).ind]);
	}
	if(abs(DisLep.at(0).id)==11){ //for DisLep.at(0)  //DisLep.at(0) is electron
	  h.Iso0_d0->Fill(Electron_pfRelIso03_all[DisLep.at(0).ind]);
	  h.dxy0_d0->Fill(Electron_dxy[DisLep.at(0).ind]);
	  h.dz0_d0->Fill(Electron_dz[DisLep.at(0).ind]);
	  float a = Electron_dxy[DisLep.at(0).ind];
	  float b = Electron_dz[DisLep.at(0).ind];
	  float c = sqrt((a*a)+(b*b));
	  h.c0_d0->Fill(c);
	  h.ip3d0_d0->Fill(Electron_ip3d[DisLep.at(0).ind]);
	  h.sip3d0_d0->Fill(Electron_sip3d[DisLep.at(0).ind]);
	}
	else if(abs(DisLep.at(0).id)==13){  //DisLep.at(0) is muon
	  h.Iso0_d0->Fill(Muon_pfRelIso03_all[DisLep.at(0).ind]);
	  h.dxy0_d0->Fill(Muon_dxy[DisLep.at(0).ind]);
	  h.dz0_d0->Fill(Muon_dz[DisLep.at(0).ind]);
	  float a = Muon_dxy[DisLep.at(0).ind];
	  float b = Muon_dz[DisLep.at(0).ind];
	  float c = sqrt((a*a)+(b*b));
	  h.c0_d0->Fill(c);
	  h.ip3d0_d0->Fill(Muon_ip3d[DisLep.at(0).ind]);
	  h.sip3d0_d0->Fill(Muon_sip3d[DisLep.at(0).ind]);
	}
	if(abs(DisLep.at(1).id)==11){ //for DisLep.at(1)  //DisLep.at(1) is electron
	  h.Iso0_d1->Fill(Electron_pfRelIso03_all[DisLep.at(1).ind]);
	  h.dxy0_d1->Fill(Electron_dxy[DisLep.at(1).ind]);
	  h.dz0_d1->Fill(Electron_dz[DisLep.at(1).ind]);
	  float a = Electron_dxy[DisLep.at(1).ind];
	  float b = Electron_dz[DisLep.at(1).ind];
	  float c = sqrt((a*a)+(b*b));
	  h.c0_d1->Fill(c);
	  h.ip3d0_d1->Fill(Electron_ip3d[DisLep.at(1).ind]);
	  h.sip3d0_d1->Fill(Electron_sip3d[DisLep.at(1).ind]);
	}
	else if(abs(DisLep.at(1).id)==13){  //DisLep.at(1) is muon
	  h.Iso0_d1->Fill(Muon_pfRelIso03_all[DisLep.at(1).ind]);
	  h.dxy0_d1->Fill(Muon_dxy[DisLep.at(1).ind]);
	  h.dz0_d1->Fill(Muon_dz[DisLep.at(1).ind]);
	  float a = Muon_dxy[DisLep.at(1).ind];
	  float b = Muon_dz[DisLep.at(1).ind];
	  float c = sqrt((a*a)+(b*b));
	  h.c0_d1->Fill(c);
	  h.ip3d0_d1->Fill(Muon_ip3d[DisLep.at(1).ind]);
	  h.sip3d0_d1->Fill(Muon_sip3d[DisLep.at(1).ind]);
	}
	
	//MET
	h.MET0->Fill(*MET_pt);
	
	//Transverse mass-l0
	float p = ProLep.at(0).v.Pt();
	float q = *MET_pt;
	float dphi = delta_phi(p,q);
	float mT = transv_mass(p,q,dphi);
	h.mT0_l0->Fill(mT);
	//Transverse mass-d0
	float a = DisLep.at(0).v.Pt();
	float b = *MET_pt;
	float Dphi = delta_phi(a,b);
	float mt = transv_mass(a,b,Dphi);
	h.mT0_d0->Fill(mt);
	//Transverse mass-d1
	float c = DisLep.at(1).v.Pt();
	float d = *MET_pt;
	float dPhi = delta_phi(c,d);
	float Mt = transv_mass(c,d,dPhi);
	h.mT0_d1->Fill(Mt);
	
	//dRmin between jet and l0,d0,d1
	float dRmin_l0j = 10000;
	float dRmin_d0j = 10000;
	float dRmin_d1j = 10000;
	if((int)goodJet.size()>0){
	  for(unsigned int i=0; i<goodJet.size(); i++){
	    float dR_l0j = goodJet.at(i).v.DeltaR(ProLep.at(0).v);
	    float dR_d0j = goodJet.at(i).v.DeltaR(DisLep.at(0).v);
	    float dR_d1j = goodJet.at(i).v.DeltaR(DisLep.at(1).v);
	    if(dR_l0j < dRmin_l0j){
	      dRmin_l0j = dR_l0j;
	    }
	    if(dR_d0j < dRmin_d0j){
	      dRmin_d0j = dR_d0j;
	    }
	    if(dR_d1j < dRmin_d1j){
	      dRmin_d1j = dR_d1j;
	    }
	  } 
	  h.dRmin0_l0j->Fill(dRmin_l0j);
	  h.dRmin0_d0j->Fill(dRmin_d0j);
	  h.dRmin0_d1j->Fill(dRmin_d1j);

	  //HT
	  float HT=0;
	  for(unsigned int i=0; i<goodJet.size(); i++){
	    HT = HT + goodJet.at(i).v.Pt();
	  }
	  h.HT0->Fill(HT);
	} //if goodJet>0 ends

	if((int)goodJet.size()>=0){
	  h.njets0->Fill(goodJet.size());
	}
	
	//flavor classification
	int flav =10;
	if(abs(ProLep.at(0).id)==11 && abs(DisLep.at(0).id)==11 && abs(DisLep.at(1).id)==11){ 
	  flav = 0; //e e e
	}
	if(abs(ProLep.at(0).id)==11 && abs(DisLep.at(0).id)==11 && abs(DisLep.at(1).id)==13){
	  flav = 1; //e e mu 
	}
	if(abs(ProLep.at(0).id)==11 && abs(DisLep.at(0).id)==13 && abs(DisLep.at(1).id)==13){
	  flav = 2; //e mu mu
	}
	if(abs(ProLep.at(0).id)==13 && abs(DisLep.at(0).id)==13 && abs(DisLep.at(1).id)==13){
	  flav = 3; //mu mu mu
	}
	if(abs(ProLep.at(0).id)==13 && abs(DisLep.at(0).id)==13 && abs(DisLep.at(1).id)==11){
	  flav = 4; //mu mu e
	}
	if(abs(ProLep.at(0).id)==13 && abs(DisLep.at(0).id)==11 && abs(DisLep.at(1).id)==11){
	  flav = 5; //mu e e
	}
	if(abs(ProLep.at(0).id)==11 && abs(DisLep.at(0).id)==13 && abs(DisLep.at(1).id)==11){
	  flav = 6; //e mu e
	}
	if(abs(ProLep.at(0).id)==13 && abs(DisLep.at(0).id)==11 && abs(DisLep.at(1).id)==13){
	  flav = 7; //mu e mu
	}
	h.flav0->Fill(flav);
      } //if ProLep>0 && DisLep>1 ends
      
      
      
      //############################################################
      
      
      
      //#################### 1l2d channel - J/Psi veto #########################
      
      
      
      //1l2d channel
      if((int)ProLep.size()>0 && (int)DisLep.size()>1){  //atleast 1 prompt and atleast 2 displaced
	if(((ProLep.at(0).v + DisLep.at(0).v).M() < 2.6 or (ProLep.at(0).v + DisLep.at(0).v).M() > 4.2) &&  //m_l0d0 < 2.6 and m_l0d0 > 4.2 
	   ((DisLep.at(0).v + DisLep.at(1).v).M() < 2.6 or (DisLep.at(0).v + DisLep.at(1).v).M() > 4.2) &&  //m_d0d1 < 2.6 and m_d0d1 > 4.2
	   ((DisLep.at(0).v.DeltaR(DisLep.at(1).v)) > 0.3)){  //dR_d0d1 > 0.3
	

	  //mass plots
	  h.m0_Jveto_l0d0->Fill((ProLep.at(0).v + DisLep.at(0).v).M());
	  h.m0_Jveto_d0d1->Fill((DisLep.at(0).v + DisLep.at(1).v).M());
	  h.m0_Jveto_l0d1->Fill((ProLep.at(0).v + DisLep.at(1).v).M());
	  h.m0_Jveto_l0d0d1->Fill((ProLep.at(0).v + DisLep.at(0).v+DisLep.at(1).v).M());
	  
	  //dPhi plots
	  h.dPhi0_Jveto_l0d0->Fill(ProLep.at(0).v.DeltaPhi(DisLep.at(0).v));
	  h.dPhi0_Jveto_d0d1->Fill(DisLep.at(0).v.DeltaPhi(DisLep.at(1).v));
	  h.dPhi0_Jveto_l0d1->Fill(ProLep.at(0).v.DeltaPhi(DisLep.at(1).v));
	  
	  //dPhi_MET_lepton
	  float dPhi0_Jveto_l0MET = delta_phi(ProLep.at(0).v.Pt(),*MET_pt);
	  h.dPhi0_Jveto_l0MET->Fill(dPhi0_Jveto_l0MET);
	  float dPhi0_Jveto_d0MET = delta_phi(DisLep.at(0).v.Pt(),*MET_pt);
	  h.dPhi0_Jveto_d0MET->Fill(dPhi0_Jveto_d0MET);
	  float dPhi0_Jveto_d1MET = delta_phi(DisLep.at(1).v.Pt(),*MET_pt);
	  h.dPhi0_Jveto_d1MET->Fill(dPhi0_Jveto_d1MET);
	  
	  //dR plots
	  h.dR0_Jveto_l0d0->Fill(ProLep.at(0).v.DeltaR(DisLep.at(0).v));
	  h.dR0_Jveto_d0d1->Fill(DisLep.at(0).v.DeltaR(DisLep.at(1).v));
	  h.dR0_Jveto_l0d1->Fill(ProLep.at(0).v.DeltaR(DisLep.at(1).v));
	  
	  //pT plots
	  h.pT0_Jveto_l0->Fill(ProLep.at(0).v.Pt());
	  h.pT0_Jveto_d0->Fill(DisLep.at(0).v.Pt());
	  h.pT0_Jveto_d1->Fill(DisLep.at(1).v.Pt());
	  h.pT0_Jveto_l0d0->Fill((ProLep.at(0).v + DisLep.at(0).v).Pt());
	  h.pT0_Jveto_d0d1->Fill((DisLep.at(0).v + DisLep.at(1).v).Pt());
	  h.pT0_Jveto_l0d1->Fill((ProLep.at(0).v + DisLep.at(1).v).Pt());
	  h.pT0_Jveto_l0d0d1->Fill((ProLep.at(0).v + DisLep.at(0).v + DisLep.at(1).v).Pt());
	  
	  //Isolation, dxy, dz, sqrt(dxy^2 + dz^2)
	  if(abs(ProLep.at(0).id)==11){ //for ProLep.at(0)  //ProLep.at(0) is electron
	    h.Iso0_Jveto_l0->Fill(Electron_pfRelIso03_all[ProLep.at(0).ind]);
	    h.dxy0_Jveto_l0->Fill(Electron_dxy[ProLep.at(0).ind]);
	    h.dz0_Jveto_l0->Fill(Electron_dz[ProLep.at(0).ind]);
	    float a = Electron_dxy[ProLep.at(0).ind];
	    float b = Electron_dz[ProLep.at(0).ind];
	    float c = sqrt((a*a)+(b*b));
	    h.c0_Jveto_l0->Fill(c);
	    h.ip3d0_Jveto_l0->Fill(Electron_ip3d[ProLep.at(0).ind]);
	    h.sip3d0_Jveto_l0->Fill(Electron_sip3d[ProLep.at(0).ind]);
	  }
	  else if(abs(ProLep.at(0).id)==13){  //ProLep.at(0) is muon
	    h.Iso0_Jveto_l0->Fill(Muon_pfRelIso03_all[ProLep.at(0).ind]);
	    h.dxy0_Jveto_l0->Fill(Muon_dxy[ProLep.at(0).ind]);
	    h.dz0_Jveto_l0->Fill(Muon_dz[ProLep.at(0).ind]);
	    float a = Muon_dxy[ProLep.at(0).ind];
	    float b = Muon_dz[ProLep.at(0).ind];
	    float c = sqrt((a*a)+(b*b));
	    h.c0_Jveto_l0->Fill(c);
	    h.ip3d0_Jveto_l0->Fill(Muon_ip3d[ProLep.at(0).ind]);
	    h.sip3d0_Jveto_l0->Fill(Muon_sip3d[ProLep.at(0).ind]);
	  }
	  if(abs(DisLep.at(0).id)==11){ //for DisLep.at(0)  //DisLep.at(0) is electron
	    h.Iso0_Jveto_d0->Fill(Electron_pfRelIso03_all[DisLep.at(0).ind]);
	    h.dxy0_Jveto_d0->Fill(Electron_dxy[DisLep.at(0).ind]);
	    h.dz0_Jveto_d0->Fill(Electron_dz[DisLep.at(0).ind]);
	    float a = Electron_dxy[DisLep.at(0).ind];
	    float b = Electron_dz[DisLep.at(0).ind];
	    float c = sqrt((a*a)+(b*b));
	    h.c0_Jveto_d0->Fill(c);
	    h.ip3d0_Jveto_d0->Fill(Electron_ip3d[DisLep.at(0).ind]);
	    h.sip3d0_Jveto_d0->Fill(Electron_sip3d[DisLep.at(0).ind]);
	  }
	  else if(abs(DisLep.at(0).id)==13){  //DisLep.at(0) is muon
	    h.Iso0_Jveto_d0->Fill(Muon_pfRelIso03_all[DisLep.at(0).ind]);
	    h.dxy0_Jveto_d0->Fill(Muon_dxy[DisLep.at(0).ind]);
	    h.dz0_Jveto_d0->Fill(Muon_dz[DisLep.at(0).ind]);
	    float a = Muon_dxy[DisLep.at(0).ind];
	    float b = Muon_dz[DisLep.at(0).ind];
	    float c = sqrt((a*a)+(b*b));
	    h.c0_Jveto_d0->Fill(c);
	    h.ip3d0_Jveto_d0->Fill(Muon_ip3d[DisLep.at(0).ind]);
	    h.sip3d0_Jveto_d0->Fill(Muon_sip3d[DisLep.at(0).ind]);
	  }
	  if(abs(DisLep.at(1).id)==11){ //for DisLep.at(1)  //DisLep.at(1) is electron
	    h.Iso0_Jveto_d1->Fill(Electron_pfRelIso03_all[DisLep.at(1).ind]);
	    h.dxy0_Jveto_d1->Fill(Electron_dxy[DisLep.at(1).ind]);
	    h.dz0_Jveto_d1->Fill(Electron_dz[DisLep.at(1).ind]);
	    float a = Electron_dxy[DisLep.at(1).ind];
	    float b = Electron_dz[DisLep.at(1).ind];
	    float c = sqrt((a*a)+(b*b));
	    h.c0_Jveto_d1->Fill(c);
	    h.ip3d0_Jveto_d1->Fill(Electron_ip3d[DisLep.at(1).ind]);
	    h.sip3d0_Jveto_d1->Fill(Electron_sip3d[DisLep.at(1).ind]);
	  }
	  else if(abs(DisLep.at(1).id)==13){  //DisLep.at(1) is muon
	    h.Iso0_Jveto_d1->Fill(Muon_pfRelIso03_all[DisLep.at(1).ind]);
	    h.dxy0_Jveto_d1->Fill(Muon_dxy[DisLep.at(1).ind]);
	    h.dz0_Jveto_d1->Fill(Muon_dz[DisLep.at(1).ind]);
	    float a = Muon_dxy[DisLep.at(1).ind];
	    float b = Muon_dz[DisLep.at(1).ind];
	    float c = sqrt((a*a)+(b*b));
	    h.c0_Jveto_d1->Fill(c);
	    h.ip3d0_Jveto_d1->Fill(Muon_ip3d[DisLep.at(1).ind]);
	    h.sip3d0_Jveto_d1->Fill(Muon_sip3d[DisLep.at(1).ind]);
	  }
	  
	  //MET
	  h.MET0_Jveto->Fill(*MET_pt);
	  
	  //Transverse mass-l0
	  float p = ProLep.at(0).v.Pt();
	  float q = *MET_pt;
	  float dphi = delta_phi(p,q);
	  float mT = transv_mass(p,q,dphi);
	  h.mT0_Jveto_l0->Fill(mT);
	  //Transverse mass-d0
	  float a = DisLep.at(0).v.Pt();
	  float b = *MET_pt;
	  float Dphi = delta_phi(a,b);
	  float mt = transv_mass(a,b,Dphi);
	  h.mT0_Jveto_d0->Fill(mt);
	  //Transverse mass-d1
	  float c = DisLep.at(1).v.Pt();
	  float d = *MET_pt;
	  float dPhi = delta_phi(c,d);
	  float Mt = transv_mass(c,d,dPhi);
	  h.mT0_Jveto_d1->Fill(Mt);
	  
	  //dRmin between jet and l0,d0,d1
	  float dRmin_l0j = 10000;
	  float dRmin_d0j = 10000;
	  float dRmin_d1j = 10000;
	  if((int)goodJet.size()>0){
	    for(unsigned int i=0; i<goodJet.size(); i++){
	      float dR_l0j = goodJet.at(i).v.DeltaR(ProLep.at(0).v);
	      float dR_d0j = goodJet.at(i).v.DeltaR(DisLep.at(0).v);
	      float dR_d1j = goodJet.at(i).v.DeltaR(DisLep.at(1).v);
	      if(dR_l0j < dRmin_l0j){
		dRmin_l0j = dR_l0j;
	      }
	      if(dR_d0j < dRmin_d0j){
		dRmin_d0j = dR_d0j;
	      }
	      if(dR_d1j < dRmin_d1j){
		dRmin_d1j = dR_d1j;
	      }
	    } 
	    h.dRmin0_Jveto_l0j->Fill(dRmin_l0j);
	    h.dRmin0_Jveto_d0j->Fill(dRmin_d0j);
	    h.dRmin0_Jveto_d1j->Fill(dRmin_d1j);
	    
	    //HT
	    float HT=0;
	    for(unsigned int i=0; i<goodJet.size(); i++){
	      HT = HT + goodJet.at(i).v.Pt();
	    }
	    h.HT0_Jveto->Fill(HT);
	  } //if goodJet>0 ends
	  
	  if((int)goodJet.size()>=0){
	    h.njets0_Jveto->Fill(goodJet.size());
	  }
	  
	  //flavor classification
	  int flav =10;
	  if(abs(ProLep.at(0).id)==11 && abs(DisLep.at(0).id)==11 && abs(DisLep.at(1).id)==11){ 
	    flav = 0; //e e e
	  }
	  if(abs(ProLep.at(0).id)==11 && abs(DisLep.at(0).id)==11 && abs(DisLep.at(1).id)==13){
	    flav = 1; //e e mu 
	  }
	  if(abs(ProLep.at(0).id)==11 && abs(DisLep.at(0).id)==13 && abs(DisLep.at(1).id)==13){
	    flav = 2; //e mu mu
	  }
	  if(abs(ProLep.at(0).id)==13 && abs(DisLep.at(0).id)==13 && abs(DisLep.at(1).id)==13){
	  flav = 3; //mu mu mu
	  }
	  if(abs(ProLep.at(0).id)==13 && abs(DisLep.at(0).id)==13 && abs(DisLep.at(1).id)==11){
	    flav = 4; //mu mu e
	  }
	  if(abs(ProLep.at(0).id)==13 && abs(DisLep.at(0).id)==11 && abs(DisLep.at(1).id)==11){
	    flav = 5; //mu e e
	  }
	  if(abs(ProLep.at(0).id)==11 && abs(DisLep.at(0).id)==13 && abs(DisLep.at(1).id)==11){
	    flav = 6; //e mu e
	  }
	  if(abs(ProLep.at(0).id)==13 && abs(DisLep.at(0).id)==11 && abs(DisLep.at(1).id)==13){
	    flav = 7; //mu e mu
	  }
	  h.flav0_Jveto->Fill(flav);
	}
      } //if ProLep>0 && DisLep>1 ends
           
      
      
      //############################################################
      
      
      
      //#################### 1l2d channel - d0, d1 = e ########################
      
      
      
      //1l2d channel
      if((int)ProLep.size()>0 && (int)DisLep.size()>1){  //atleast 1 prompt and atleast 2 displaced
	if(((ProLep.at(0).v + DisLep.at(0).v).M() < 2.6 or (ProLep.at(0).v + DisLep.at(0).v).M() > 4.2) &&  //m_l0d0 < 2.6 and m_l0d0 > 4.2 
	   ((DisLep.at(0).v + DisLep.at(1).v).M() < 2.6 or (DisLep.at(0).v + DisLep.at(1).v).M() > 4.2) &&  //m_d0d1 < 2.6 and m_d0d1 > 4.2
	   ((DisLep.at(0).v.DeltaR(DisLep.at(1).v)) > 0.3)){  //dR_d0d1 > 0.3
	  if((abs(DisLep.at(0).id) == 11) && (abs(DisLep.at(1).id) == 11)){
	    
	    //mass plots
	    h.m0_ee_l0d0->Fill((ProLep.at(0).v + DisLep.at(0).v).M());
	    h.m0_ee_l0d1->Fill((ProLep.at(0).v + DisLep.at(1).v).M());
	    h.m0_ee_d0d1->Fill((DisLep.at(0).v + DisLep.at(1).v).M());
	    h.m0_ee_l0d0d1->Fill((ProLep.at(0).v + DisLep.at(0).v + DisLep.at(1).v).M());
	    
	    //dPhi plots
	    h.dPhi0_ee_l0d0->Fill(ProLep.at(0).v.DeltaPhi(DisLep.at(0).v));
	    h.dPhi0_ee_l0d1->Fill(ProLep.at(0).v.DeltaPhi(DisLep.at(1).v));
	    h.dPhi0_ee_d0d1->Fill(DisLep.at(0).v.DeltaPhi(DisLep.at(1).v));
	    
	    //dR plots
	    h.dR0_ee_l0d0->Fill(ProLep.at(0).v.DeltaR(DisLep.at(0).v));
	    h.dR0_ee_l0d1->Fill(ProLep.at(0).v.DeltaR(DisLep.at(1).v));
	    h.dR0_ee_d0d1->Fill(DisLep.at(0).v.DeltaR(DisLep.at(1).v));
	    
	    //pT plots
	    h.pT0_ee_l0->Fill(ProLep.at(0).v.Pt());
	    h.pT0_ee_d0->Fill(DisLep.at(0).v.Pt());
	    h.pT0_ee_d1->Fill(DisLep.at(1).v.Pt());
	    
	    //MET
	    h.MET0_ee->Fill(*MET_pt);

	    //Isolation plots
	    if(abs(ProLep.at(0).id)==11){
	      h.Iso0_ee_l0->Fill(Electron_pfRelIso03_all[ProLep.at(0).ind]);
	    }
	    else if(abs(ProLep.at(0).id)==13){
	      h.Iso0_ee_l0->Fill(Muon_pfRelIso03_all[ProLep.at(0).ind]);
	    }
	    h.Iso0_ee_d0->Fill(Electron_pfRelIso03_all[DisLep.at(0).ind]);
	    h.Iso0_ee_d1->Fill(Electron_pfRelIso03_all[DisLep.at(1).ind]);
	    
	    //dxy, dz, ip3d, sip3d
	    h.dxy0_ee_d0->Fill(Electron_dxy[DisLep.at(0).ind]);
	    h.dz0_ee_d0->Fill(Electron_dz[DisLep.at(0).ind]);
	    h.ip3d0_ee_d0->Fill(Electron_ip3d[DisLep.at(0).ind]);
	    h.sip3d0_ee_d0->Fill(Electron_sip3d[DisLep.at(0).ind]);
	    
	    //dxy, dz, ip3d, sip3d
	    h.dxy0_ee_d1->Fill(Electron_dxy[DisLep.at(1).ind]);
	    h.dz0_ee_d1->Fill(Electron_dz[DisLep.at(1).ind]);
	    h.ip3d0_ee_d1->Fill(Electron_ip3d[DisLep.at(1).ind]);
	    h.sip3d0_ee_d1->Fill(Electron_sip3d[DisLep.at(1).ind]);
	    
	    //Transverse mass-l0
	    float p = ProLep.at(0).v.Pt();
	    float q = *MET_pt;
	    float dphi = delta_phi(p,q);
	    float mT = transv_mass(p,q,dphi);
	    h.mT0_ee_l0->Fill(mT);
	    //Transverse mass-d0
	    float a = DisLep.at(0).v.Pt();
	    float b = *MET_pt;
	    float Dphi = delta_phi(a,b);
	    float mt = transv_mass(a,b,Dphi);
	    h.mT0_ee_d0->Fill(mt);
	    //Transverse mass-d1
	    float c = DisLep.at(1).v.Pt();
	    float d = *MET_pt;
	    float dPhi = delta_phi(c,d);
	    float Mt = transv_mass(c,d,dPhi);
	    h.mT0_ee_d1->Fill(Mt);
	    
	    if((int)goodJet.size()>=0){
	      h.njets0_ee->Fill(goodJet.size());
	    }
	    
	    //flavor classification
	    int flav =10;
	    if(abs(ProLep.at(0).id)==11 && abs(DisLep.at(0).id)==11 && abs(DisLep.at(1).id)==11){ 
	      flav = 0; //e e e
	    }
	    if(abs(ProLep.at(0).id)==11 && abs(DisLep.at(0).id)==11 && abs(DisLep.at(1).id)==13){
	      flav = 1; //e e mu 
	    }
	    if(abs(ProLep.at(0).id)==11 && abs(DisLep.at(0).id)==13 && abs(DisLep.at(1).id)==13){
	      flav = 2; //e mu mu
	    }
	    if(abs(ProLep.at(0).id)==13 && abs(DisLep.at(0).id)==13 && abs(DisLep.at(1).id)==13){
	      flav = 3; //mu mu mu
	    }
	    if(abs(ProLep.at(0).id)==13 && abs(DisLep.at(0).id)==13 && abs(DisLep.at(1).id)==11){
	      flav = 4; //mu mu e
	    }
	    if(abs(ProLep.at(0).id)==13 && abs(DisLep.at(0).id)==11 && abs(DisLep.at(1).id)==11){
	      flav = 5; //mu e e
	    }
	    if(abs(ProLep.at(0).id)==11 && abs(DisLep.at(0).id)==13 && abs(DisLep.at(1).id)==11){
	      flav = 6; //e mu e
	    }
	    if(abs(ProLep.at(0).id)==13 && abs(DisLep.at(0).id)==11 && abs(DisLep.at(1).id)==13){
	      flav = 7; //mu e mu
	    }
	    h.flav0_ee->Fill(flav);
	  }  //if abs(DisLep.at(0)) == 11 && abs(DisLep.at(1)) == 11 ends
	}  //if not J/Psi (mass cuts to exclude J/Psi ends here)
      } //if ProLep>1 && DisLep>0 ends

      
      
      //############################################################
      
      
      
      //#################### 1l2d channel - d0-mu, d1-e #######################
      
      
      
      //1l2d channel
      if((int)ProLep.size()>0 && (int)DisLep.size()>1){  //atleast 1 prompt and atleast 2 displaced
	if(((ProLep.at(0).v + DisLep.at(0).v).M() < 2.6 or (ProLep.at(0).v + DisLep.at(0).v).M() > 4.2) &&  //m_l0d0 < 2.6 and m_l0d0 > 4.2 
	   ((DisLep.at(0).v + DisLep.at(1).v).M() < 2.6 or (DisLep.at(0).v + DisLep.at(1).v).M() > 4.2) &&  //m_d0d1 < 2.6 and m_d0d1 > 4.2
	   ((DisLep.at(0).v.DeltaR(DisLep.at(1).v)) > 0.3)){  //dR_d0d1 > 0.3
	  if((abs(DisLep.at(0).id) == 13) && (abs(DisLep.at(1).id) == 11)){
	    
	    //mass plots
	    h.m0_MuE_l0d0->Fill((ProLep.at(0).v + DisLep.at(0).v).M());
	    h.m0_MuE_l0d1->Fill((ProLep.at(0).v + DisLep.at(1).v).M());
	    h.m0_MuE_d0d1->Fill((DisLep.at(0).v + DisLep.at(1).v).M());
	    h.m0_MuE_l0d0d1->Fill((ProLep.at(0).v + DisLep.at(0).v + DisLep.at(1).v).M());
	    
	    //dPhi plots
	    h.dPhi0_MuE_l0d0->Fill(ProLep.at(0).v.DeltaPhi(DisLep.at(0).v));
	    h.dPhi0_MuE_l0d1->Fill(ProLep.at(0).v.DeltaPhi(DisLep.at(1).v));
	    h.dPhi0_MuE_d0d1->Fill(DisLep.at(0).v.DeltaPhi(DisLep.at(1).v));
	    
	    //dR plots
	    h.dR0_MuE_l0d0->Fill(ProLep.at(0).v.DeltaR(DisLep.at(0).v));
	    h.dR0_MuE_l0d1->Fill(ProLep.at(0).v.DeltaR(DisLep.at(1).v));
	    h.dR0_MuE_d0d1->Fill(DisLep.at(0).v.DeltaR(DisLep.at(1).v));
	    
	    //pT plots
	    h.pT0_MuE_l0->Fill(ProLep.at(0).v.Pt());
	    h.pT0_MuE_d0->Fill(DisLep.at(0).v.Pt());
	    h.pT0_MuE_d1->Fill(DisLep.at(1).v.Pt());
	    
	    //MET
	    h.MET0_MuE->Fill(*MET_pt);
	    
	    //Isolation plots
	    if(abs(ProLep.at(0).id)==11){
	      h.Iso0_MuE_l0->Fill(Electron_pfRelIso03_all[ProLep.at(0).ind]);
	    }
	    else if(abs(ProLep.at(0).id)==13){
	      h.Iso0_MuE_l0->Fill(Muon_pfRelIso03_all[ProLep.at(0).ind]);
	    }
	    h.Iso0_MuE_d0->Fill(Muon_pfRelIso03_all[DisLep.at(0).ind]);
	    h.Iso0_MuE_d1->Fill(Electron_pfRelIso03_all[DisLep.at(1).ind]);
	    
	    //dxy, dz, ip3d, sip3d_d0
	    h.dxy0_MuE_d0->Fill(Muon_dxy[DisLep.at(0).ind]);
	    h.dz0_MuE_d0->Fill(Muon_dz[DisLep.at(0).ind]);
	    h.ip3d0_MuE_d0->Fill(Muon_ip3d[DisLep.at(0).ind]);
	    h.sip3d0_MuE_d0->Fill(Muon_sip3d[DisLep.at(0).ind]);
	    
	    //dxy, dz, ip3d, sip3d_d1
	    h.dxy0_MuE_d1->Fill(Electron_dxy[DisLep.at(1).ind]);
	    h.dz0_MuE_d1->Fill(Electron_dz[DisLep.at(1).ind]);
	    h.ip3d0_MuE_d1->Fill(Electron_ip3d[DisLep.at(1).ind]);
	    h.sip3d0_MuE_d1->Fill(Electron_sip3d[DisLep.at(1).ind]);
	    
	    //Transverse mass-l0
	    float p = ProLep.at(0).v.Pt();
	    float q = *MET_pt;
	    float dphi = delta_phi(p,q);
	    float mT = transv_mass(p,q,dphi);
	    h.mT0_MuE_l0->Fill(mT);
	    //Transverse mass-d0
	    float a = DisLep.at(0).v.Pt();
	    float b = *MET_pt;
	    float Dphi = delta_phi(a,b);
	    float mt = transv_mass(a,b,Dphi);
	    h.mT0_MuE_d0->Fill(mt);
	    //Transverse mass-d1
	    float c = DisLep.at(1).v.Pt();
	    float d = *MET_pt;
	    float dPhi = delta_phi(c,d);
	    float Mt = transv_mass(c,d,dPhi);
	    h.mT0_MuE_d1->Fill(Mt);
	    
	    if((int)goodJet.size()>=0){
	      h.njets0_MuE->Fill(goodJet.size());
	    }
	    
	    //flavor classification
	    int flav =10;
	    if(abs(ProLep.at(0).id)==11 && abs(DisLep.at(0).id)==11 && abs(DisLep.at(1).id)==11){ 
	      flav = 0; //e e e
	    }
	    if(abs(ProLep.at(0).id)==11 && abs(DisLep.at(0).id)==11 && abs(DisLep.at(1).id)==13){
	      flav = 1; //e e mu 
	    }
	    if(abs(ProLep.at(0).id)==11 && abs(DisLep.at(0).id)==13 && abs(DisLep.at(1).id)==13){
	      flav = 2; //e mu mu
	    }
	    if(abs(ProLep.at(0).id)==13 && abs(DisLep.at(0).id)==13 && abs(DisLep.at(1).id)==13){
	      flav = 3; //mu mu mu
	    }
	    if(abs(ProLep.at(0).id)==13 && abs(DisLep.at(0).id)==13 && abs(DisLep.at(1).id)==11){
	      flav = 4; //mu mu e
	    }
	    if(abs(ProLep.at(0).id)==13 && abs(DisLep.at(0).id)==11 && abs(DisLep.at(1).id)==11){
	      flav = 5; //mu e e
	    }
	    if(abs(ProLep.at(0).id)==11 && abs(DisLep.at(0).id)==13 && abs(DisLep.at(1).id)==11){
	      flav = 6; //e mu e
	    }
	    if(abs(ProLep.at(0).id)==13 && abs(DisLep.at(0).id)==11 && abs(DisLep.at(1).id)==13){
	      flav = 7; //mu e mu
	    }
	    h.flav0_MuE->Fill(flav);
	  }  //if abs(DisLep.at(0)) == 13 && abs(DisLep.at(1)) == 11 ends
	}  //if not J/Psi (mass cuts to exclude J/Psi ends here)
      } //if ProLep>1 && DisLep>0 ends
                 
      
      
      //############################################################
      
      
      
      //#################### 1l2d channel - d0, d1 = mu ########################
      
      
      
      //1l2d channel
      if((int)ProLep.size()>0 && (int)DisLep.size()>1){  //atleast 1 prompt and atleast 2 displaced
	if(((ProLep.at(0).v + DisLep.at(0).v).M() < 2.6 or (ProLep.at(0).v + DisLep.at(0).v).M() > 4.2) &&  //m_l0d0 < 2.6 and m_l0d0 > 4.2 
	   ((DisLep.at(0).v + DisLep.at(1).v).M() < 2.6 or (DisLep.at(0).v + DisLep.at(1).v).M() > 4.2) &&  //m_d0d1 < 2.6 and m_d0d1 > 4.2
	   ((DisLep.at(0).v.DeltaR(DisLep.at(1).v)) > 0.3)){  //dR_d0d1 > 0.3
	  if((abs(DisLep.at(0).id) == 13) && (abs(DisLep.at(1).id) == 13)){
	    
	    //mass plots
	    h.m0_MuMu_l0d0->Fill((ProLep.at(0).v + DisLep.at(0).v).M());
	    h.m0_MuMu_l0d1->Fill((ProLep.at(0).v + DisLep.at(1).v).M());
	    h.m0_MuMu_d0d1->Fill((DisLep.at(0).v + DisLep.at(1).v).M());
	    h.m0_MuMu_l0d0d1->Fill((ProLep.at(0).v + DisLep.at(0).v + DisLep.at(1).v).M());
	    
	    //dPhi plots
	    h.dPhi0_MuMu_l0d0->Fill(ProLep.at(0).v.DeltaPhi(DisLep.at(0).v));
	    h.dPhi0_MuMu_l0d1->Fill(ProLep.at(0).v.DeltaPhi(DisLep.at(1).v));
	    h.dPhi0_MuMu_d0d1->Fill(DisLep.at(0).v.DeltaPhi(DisLep.at(1).v));
	    
	    //dR plots
	    h.dR0_MuMu_l0d0->Fill(ProLep.at(0).v.DeltaR(DisLep.at(0).v));
	    h.dR0_MuMu_l0d1->Fill(ProLep.at(0).v.DeltaR(DisLep.at(1).v));
	    h.dR0_MuMu_d0d1->Fill(DisLep.at(0).v.DeltaR(DisLep.at(1).v));
	    
	    //pT plots
	    h.pT0_MuMu_l0->Fill(ProLep.at(0).v.Pt());
	    h.pT0_MuMu_d0->Fill(DisLep.at(0).v.Pt());
	    h.pT0_MuMu_d1->Fill(DisLep.at(1).v.Pt());
	    
	    //MET
	    h.MET0_MuMu->Fill(*MET_pt);
	    
	    //Isolation plots
	    if(abs(ProLep.at(0).id)==11){
	      h.Iso0_MuMu_l0->Fill(Electron_pfRelIso03_all[ProLep.at(0).ind]);
	    }
	    else if(abs(ProLep.at(0).id)==13){
	      h.Iso0_MuMu_l0->Fill(Muon_pfRelIso03_all[ProLep.at(0).ind]);
	    }
	    h.Iso0_MuMu_d0->Fill(Muon_pfRelIso03_all[DisLep.at(0).ind]);
	    h.Iso0_MuMu_d1->Fill(Muon_pfRelIso03_all[DisLep.at(1).ind]);
	    
	    //dxy, dz, ip3d, sip3d
	    h.dxy0_MuMu_d0->Fill(Electron_dxy[DisLep.at(0).ind]);
	    h.dz0_MuMu_d0->Fill(Electron_dz[DisLep.at(0).ind]);
	    h.ip3d0_MuMu_d0->Fill(Electron_ip3d[DisLep.at(0).ind]);
	    h.sip3d0_MuMu_d0->Fill(Electron_sip3d[DisLep.at(0).ind]);
	    
	    //dxy, dz, ip3d, sip3d
	    h.dxy0_MuMu_d1->Fill(Electron_dxy[DisLep.at(1).ind]);
	    h.dz0_MuMu_d1->Fill(Electron_dz[DisLep.at(1).ind]);
	    h.ip3d0_MuMu_d1->Fill(Electron_ip3d[DisLep.at(1).ind]);
	    h.sip3d0_MuMu_d1->Fill(Electron_sip3d[DisLep.at(1).ind]);
	    
	    //Transverse mass-l0
	    float p = ProLep.at(0).v.Pt();
	    float q = *MET_pt;
	    float dphi = delta_phi(p,q);
	    float mT = transv_mass(p,q,dphi);
	    h.mT0_MuMu_l0->Fill(mT);
	    //Transverse mass-d0
	    float a = DisLep.at(0).v.Pt();
	    float b = *MET_pt;
	    float Dphi = delta_phi(a,b);
	    float mt = transv_mass(a,b,Dphi);
	    h.mT0_MuMu_d0->Fill(mt);
	    //Transverse mass-d1
	    float c = DisLep.at(1).v.Pt();
	    float d = *MET_pt;
	    float dPhi = delta_phi(c,d);
	    float Mt = transv_mass(c,d,dPhi);
	    h.mT0_MuMu_d1->Fill(Mt);
	    
	    if((int)goodJet.size()>=0){
	      h.njets0_MuMu->Fill(goodJet.size());
	    }
	    
	    //flavor classification
	    int flav =10;
	    if(abs(ProLep.at(0).id)==11 && abs(DisLep.at(0).id)==11 && abs(DisLep.at(1).id)==11){ 
	      flav = 0; //e e e
	    }
	    if(abs(ProLep.at(0).id)==11 && abs(DisLep.at(0).id)==11 && abs(DisLep.at(1).id)==13){
	      flav = 1; //e e mu 
	    }
	    if(abs(ProLep.at(0).id)==11 && abs(DisLep.at(0).id)==13 && abs(DisLep.at(1).id)==13){
	      flav = 2; //e mu mu
	    }
	    if(abs(ProLep.at(0).id)==13 && abs(DisLep.at(0).id)==13 && abs(DisLep.at(1).id)==13){
	      flav = 3; //mu mu mu
	    }
	    if(abs(ProLep.at(0).id)==13 && abs(DisLep.at(0).id)==13 && abs(DisLep.at(1).id)==11){
	      flav = 4; //mu mu e
	    }
	    if(abs(ProLep.at(0).id)==13 && abs(DisLep.at(0).id)==11 && abs(DisLep.at(1).id)==11){
	      flav = 5; //mu e e
	    }
	    if(abs(ProLep.at(0).id)==11 && abs(DisLep.at(0).id)==13 && abs(DisLep.at(1).id)==11){
	      flav = 6; //e mu e
	    }
	    if(abs(ProLep.at(0).id)==13 && abs(DisLep.at(0).id)==11 && abs(DisLep.at(1).id)==13){
	      flav = 7; //mu e mu
	    }
	    h.flav0_MuMu->Fill(flav);
	  }  //if abs(DisLep.at(0)) == 11 && abs(DisLep.at(1)) == 11 ends
	}  //if not J/Psi (mass cuts to exclude J/Psi ends here)
      } //if ProLep>1 && DisLep>0 ends
      
      
      
      //############################################################
      //############################################################
      
      
      
      //####################### 0l3d channel ############################
      
      
      
      //0l3d channel
      if((int)ProLep.size()>=0 && (int)DisLep.size()>2){ //atleast 0 prompt and atleast 3 displaced

	
	//mass plots
	h.m1_d0d1->Fill((DisLep.at(0).v + DisLep.at(1).v).M());
	h.m1_d1d2->Fill((DisLep.at(1).v + DisLep.at(2).v).M());
	h.m1_d0d2->Fill((DisLep.at(0).v + DisLep.at(2).v).M());
	h.m1_d0d1d2->Fill((DisLep.at(0).v + DisLep.at(1).v + DisLep.at(2).v).M());
	
	//dPhi plots
	h.dPhi1_d0d1->Fill(DisLep.at(0).v.DeltaPhi(DisLep.at(1).v));
	h.dPhi1_d1d2->Fill(DisLep.at(1).v.DeltaPhi(DisLep.at(2).v));
	h.dPhi1_d0d2->Fill(DisLep.at(0).v.DeltaPhi(DisLep.at(2).v));
	
	//dPhi_MET_lepton
        float dPhi1_d0MET = delta_phi(DisLep.at(0).v.Pt(),*MET_pt);
	h.dPhi1_d0MET->Fill(dPhi1_d0MET);
        float dPhi1_d1MET = delta_phi(DisLep.at(1).v.Pt(),*MET_pt);
	h.dPhi1_d1MET->Fill(dPhi1_d1MET);
        float dPhi1_d2MET = delta_phi(DisLep.at(2).v.Pt(),*MET_pt);
	h.dPhi1_d2MET->Fill(dPhi1_d2MET);
	
	//dR plots
	h.dR1_d0d1->Fill(DisLep.at(0).v.DeltaR(DisLep.at(1).v));
	h.dR1_d1d2->Fill(DisLep.at(1).v.DeltaR(DisLep.at(2).v));
	h.dR1_d0d2->Fill(DisLep.at(0).v.DeltaR(DisLep.at(2).v));
	
	//pT plots
	h.pT1_d0->Fill(DisLep.at(0).v.Pt());
	h.pT1_d1->Fill(DisLep.at(1).v.Pt());
	h.pT1_d2->Fill(DisLep.at(2).v.Pt());
	h.pT1_d0d1->Fill((DisLep.at(0).v + DisLep.at(1).v).Pt());
	h.pT1_d1d2->Fill((DisLep.at(1).v + DisLep.at(2).v).Pt());
	h.pT1_d0d2->Fill((DisLep.at(0).v + DisLep.at(2).v).Pt());
	h.pT1_d0d1d2->Fill((DisLep.at(0).v + DisLep.at(1).v + DisLep.at(2).v).Pt());
	
	//Isolation, dxy, dz, sqrt(dxy^2 + dz^2)
	if(abs(DisLep.at(0).id)==11){ //for DisLep.at(0)  //DisLep.at(0) is electron
	  h.Iso1_d0->Fill(Electron_pfRelIso03_all[DisLep.at(0).ind]);
	  h.dxy1_d0->Fill(Electron_dxy[DisLep.at(0).ind]);
	  h.dz1_d0->Fill(Electron_dz[DisLep.at(0).ind]);
	  float a = Electron_dxy[DisLep.at(0).ind];
	  float b = Electron_dz[DisLep.at(0).ind];
	  float c = sqrt((a*a)+(b*b));
	  h.c1_d0->Fill(c);
	  h.ip3d1_d0->Fill(Electron_ip3d[DisLep.at(0).ind]);
	  h.sip3d1_d0->Fill(Electron_sip3d[DisLep.at(0).ind]);
	}
	else if(abs(DisLep.at(0).id)==13){  //DisLep.at(0) is muon
	  h.Iso1_d0->Fill(Muon_pfRelIso03_all[DisLep.at(0).ind]);
	  h.dxy1_d0->Fill(Muon_dxy[DisLep.at(0).ind]);
	  h.dz1_d0->Fill(Muon_dz[DisLep.at(0).ind]);
	  float a = Muon_dxy[DisLep.at(0).ind];
	  float b = Muon_dz[DisLep.at(0).ind];
	  float c = sqrt((a*a)+(b*b));
	  h.c1_d0->Fill(c);
	  h.ip3d1_d0->Fill(Muon_ip3d[DisLep.at(0).ind]);
	  h.sip3d1_d0->Fill(Muon_sip3d[DisLep.at(0).ind]);
	}
	if(abs(DisLep.at(1).id)==11){ //for DisLep.at(1)  //DisLep.at(1) is electron
	  h.Iso1_d1->Fill(Electron_pfRelIso03_all[DisLep.at(1).ind]);
	  h.dxy1_d1->Fill(Electron_dxy[DisLep.at(1).ind]);
	  h.dz1_d1->Fill(Electron_dz[DisLep.at(1).ind]);
	  float a = Electron_dxy[DisLep.at(1).ind];
	  float b = Electron_dz[DisLep.at(1).ind];
	  float c = sqrt((a*a)+(b*b));
	  h.c1_d1->Fill(c);
	  h.ip3d1_d1->Fill(Electron_ip3d[DisLep.at(1).ind]);
	  h.sip3d1_d1->Fill(Electron_sip3d[DisLep.at(1).ind]);
	}
	else if(abs(DisLep.at(1).id)==13){  //DisLep.at(1) is muon
	  h.Iso1_d1->Fill(Muon_pfRelIso03_all[DisLep.at(1).ind]);
	  h.dxy1_d1->Fill(Muon_dxy[DisLep.at(1).ind]);
	  h.dz1_d1->Fill(Muon_dz[DisLep.at(1).ind]);
	  float a = Muon_dxy[DisLep.at(1).ind];
	  float b = Muon_dz[DisLep.at(1).ind];
	  float c = sqrt((a*a)+(b*b));
	  h.c1_d1->Fill(c);
	  h.ip3d1_d1->Fill(Muon_ip3d[DisLep.at(1).ind]);
	  h.sip3d1_d1->Fill(Muon_sip3d[DisLep.at(1).ind]);
	}
	if(abs(DisLep.at(2).id)==11){ //for DisLep.at(2)  //DisLep.at(2) is electron
	  h.Iso1_d2->Fill(Electron_pfRelIso03_all[DisLep.at(2).ind]);
	  h.dxy1_d2->Fill(Electron_dxy[DisLep.at(2).ind]);
	  h.dz1_d2->Fill(Electron_dz[DisLep.at(2).ind]);
	  float a = Electron_dxy[DisLep.at(2).ind];
	  float b = Electron_dz[DisLep.at(2).ind];
	  float c = sqrt((a*a)+(b*b));
	  h.c1_d2->Fill(c);
	  h.ip3d1_d2->Fill(Electron_ip3d[DisLep.at(2).ind]);
	  h.sip3d1_d2->Fill(Electron_sip3d[DisLep.at(2).ind]);
	}
	else if(abs(DisLep.at(2).id)==13){  //DisLep.at(2) is muon
	  h.Iso1_d2->Fill(Muon_pfRelIso03_all[DisLep.at(2).ind]);
	  h.dxy1_d2->Fill(Muon_dxy[DisLep.at(2).ind]);
	  h.dz1_d2->Fill(Muon_dz[DisLep.at(2).ind]);
	  float a = Muon_dxy[DisLep.at(2).ind];
	  float b = Muon_dz[DisLep.at(2).ind];
	  float c = sqrt((a*a)+(b*b));
	  h.c1_d2->Fill(c);
	  h.ip3d1_d2->Fill(Muon_ip3d[DisLep.at(2).ind]);
	  h.sip3d1_d2->Fill(Muon_sip3d[DisLep.at(2).ind]);
	}
	
	//MET
	h.MET1->Fill(*MET_pt);
	
	//Transverse mass-d0
	float p = DisLep.at(0).v.Pt();
	float q = *MET_pt;
	float dphi = delta_phi(p,q);
	float mT = transv_mass(p,q,dphi);
	h.mT1_d0->Fill(mT);
	//Transverse mass-d1
	float a = DisLep.at(1).v.Pt();
	float b = *MET_pt;
	float Dphi = delta_phi(a,b);
	float mt = transv_mass(a,b,Dphi);
	h.mT1_d1->Fill(mt);
	//Transverse mass-d2
	float c = DisLep.at(2).v.Pt();
	float d = *MET_pt;
	float dPhi = delta_phi(c,d);
	float Mt = transv_mass(c,d,dPhi);
	h.mT1_d2->Fill(Mt);
	
	//dRmin between jet and d0,d1,d2
	float dRmin_d0j = 10000;
	float dRmin_d1j = 10000;
	float dRmin_d2j = 10000;
	if((int)goodJet.size()>0){
	  for(unsigned int i=0; i<goodJet.size(); i++){
	    float dR_d0j = goodJet.at(i).v.DeltaR(DisLep.at(0).v);
	    float dR_d1j = goodJet.at(i).v.DeltaR(DisLep.at(1).v);
	    float dR_d2j = goodJet.at(i).v.DeltaR(DisLep.at(2).v);
	    if(dR_d0j < dRmin_d0j){
	      dRmin_d0j = dR_d0j;
	    }
	    if(dR_d1j < dRmin_d1j){
	      dRmin_d1j = dR_d1j;
	    }
	    if(dR_d2j < dRmin_d2j){
	      dRmin_d2j = dR_d2j;
	    }
	  } 
	  h.dRmin1_d0j->Fill(dRmin_d0j);
	  h.dRmin1_d1j->Fill(dRmin_d1j);
	  h.dRmin1_d2j->Fill(dRmin_d2j);
	  
	  //HT
	  float HT=0;
	  for(unsigned int i=0; i<goodJet.size(); i++){
	    HT = HT + goodJet.at(i).v.Pt();
	  }
	  h.HT1->Fill(HT);
	} //if goodJet>0 ends

	if((int)goodJet.size()>=0){
	  h.njets1->Fill(goodJet.size());
	}
	
	//flavor classification
	int flav =10;
	if(abs(DisLep.at(0).id)==11 && abs(DisLep.at(1).id)==11 && abs(DisLep.at(2).id)==11){ 
	  flav = 0; //e e e
	}
	if(abs(DisLep.at(0).id)==11 && abs(DisLep.at(1).id)==11 && abs(DisLep.at(2).id)==13){
	  flav = 1; //e e mu 
	}
	if(abs(DisLep.at(0).id)==11 && abs(DisLep.at(1).id)==13 && abs(DisLep.at(2).id)==13){
	  flav = 2; //e mu mu
	}
	if(abs(DisLep.at(0).id)==13 && abs(DisLep.at(1).id)==13 && abs(DisLep.at(2).id)==13){
	  flav = 3; //mu mu mu
	}
	if(abs(DisLep.at(0).id)==13 && abs(DisLep.at(1).id)==13 && abs(DisLep.at(2).id)==11){
	  flav = 4; //mu mu e
	}
	if(abs(DisLep.at(0).id)==13 && abs(DisLep.at(1).id)==11 && abs(DisLep.at(2).id)==11){
	  flav = 5; //mu e e
	}
	if(abs(DisLep.at(0).id)==11 && abs(DisLep.at(1).id)==13 && abs(DisLep.at(2).id)==11){
	  flav = 6; //e mu e
	}
	if(abs(DisLep.at(0).id)==13 && abs(DisLep.at(1).id)==11 && abs(DisLep.at(2).id)==13){
	  flav = 7; //mu e mu
	}
	h.flav1->Fill(flav);
      } //if ProLep>=0 && DisLep>2 ends
      
      
      
      //############################################################
      //############################################################
      
      
      
      //####################### 2l1d channel ############################
      
      
      
      //2l1d channel
      if((int)ProLep.size()>1 && (int)DisLep.size()>0 && ((ProLep.at(0).v + ProLep.at(1).v).M()>12)){ //atleast 2 prompt and atleast 1 displaced
	
	//mass plots
	h.m2_l0l1->Fill((ProLep.at(0).v + ProLep.at(1).v).M());
	h.m2_l0d0->Fill((ProLep.at(0).v + DisLep.at(0).v).M());
	h.m2_l1d0->Fill((ProLep.at(1).v + DisLep.at(0).v).M());
	h.m2_l0l1d0->Fill((ProLep.at(0).v + ProLep.at(1).v + DisLep.at(0).v).M());

	//dPhi plots
	h.dPhi2_l0l1->Fill(ProLep.at(0).v.DeltaPhi(ProLep.at(1).v));
	h.dPhi2_l0d0->Fill(ProLep.at(0).v.DeltaPhi(DisLep.at(0).v));
	h.dPhi2_l1d0->Fill(ProLep.at(1).v.DeltaPhi(DisLep.at(0).v));

	//dPhi_MET_lepton
        float dPhi2_l0MET = delta_phi(ProLep.at(0).v.Pt(),*MET_pt);
	h.dPhi2_l0MET->Fill(dPhi2_l0MET);
        float dPhi2_l1MET = delta_phi(ProLep.at(1).v.Pt(),*MET_pt);
	h.dPhi2_l1MET->Fill(dPhi2_l1MET);
        float dPhi2_d0MET = delta_phi(DisLep.at(0).v.Pt(),*MET_pt);
	h.dPhi2_d0MET->Fill(dPhi2_d0MET);

	//dR plots
	h.dR2_l0l1->Fill(ProLep.at(0).v.DeltaR(ProLep.at(1).v));
	h.dR2_l0d0->Fill(ProLep.at(0).v.DeltaR(DisLep.at(0).v));
	h.dR2_l1d0->Fill(ProLep.at(1).v.DeltaR(DisLep.at(0).v));

	//pT plots
	h.pT2_l0->Fill(ProLep.at(0).v.Pt());
	h.pT2_l1->Fill(ProLep.at(1).v.Pt());
	h.pT2_d0->Fill(DisLep.at(0).v.Pt());
	h.pT2_l0l1->Fill((ProLep.at(0).v + ProLep.at(1).v).Pt());
	h.pT2_l0d0->Fill((ProLep.at(0).v + DisLep.at(0).v).Pt());
	h.pT2_l1d0->Fill((ProLep.at(1).v + DisLep.at(0).v).Pt());
	h.pT2_l0l1d0->Fill((ProLep.at(0).v + ProLep.at(1).v + DisLep.at(0).v).Pt());

	//(Isolation, dxy, dz, c, ip3d, sip3d)
	if(abs(ProLep.at(0).id)==11){ //for ProLep.at(0)  //ProLep.at(0) is electron
	  h.Iso2_l0->Fill(Electron_pfRelIso03_all[ProLep.at(0).ind]);
	  h.dxy2_l0->Fill(Electron_dxy[ProLep.at(0).ind]);
	  h.dz2_l0->Fill(Electron_dz[ProLep.at(0).ind]);
	  float a = Electron_dxy[ProLep.at(0).ind];
	  float b = Electron_dz[ProLep.at(0).ind];
	  float c = ((a*a)+(b*b));
	  h.c2_l0->Fill(c);
	  h.ip3d2_l0->Fill(Electron_ip3d[ProLep.at(0).ind]);
	  h.sip3d2_l0->Fill(Electron_sip3d[ProLep.at(0).ind]);
	}
	else if(abs(ProLep.at(0).id)==13){  //ProLep.at(0) is muon
	  h.Iso2_l0->Fill(Muon_pfRelIso03_all[ProLep.at(0).ind]);
	  h.dxy2_l0->Fill(Muon_dxy[ProLep.at(0).ind]);
	  h.dz2_l0->Fill(Muon_dz[ProLep.at(0).ind]);
	  float a = Muon_dxy[ProLep.at(0).ind];
	  float b = Muon_dz[ProLep.at(0).ind];
	  float c = ((a*a)+(b*b));
	  h.c2_l0->Fill(c);
	  h.ip3d2_l0->Fill(Muon_ip3d[ProLep.at(0).ind]);
	  h.sip3d2_l0->Fill(Muon_sip3d[ProLep.at(0).ind]);
	}
	if(abs(ProLep.at(1).id)==11){ //for ProLep.at(1)  //ProLep.at(1) is electron
	  h.Iso2_l1->Fill(Electron_pfRelIso03_all[ProLep.at(1).ind]);
	  h.dxy2_l1->Fill(Electron_dxy[ProLep.at(1).ind]);
	  h.dz2_l1->Fill(Electron_dz[ProLep.at(1).ind]);
	  float a = Electron_dxy[ProLep.at(1).ind];
	  float b = Electron_dz[ProLep.at(1).ind];
	  float c = ((a*a)+(b*b));
	  h.c2_l1->Fill(c);
	  h.ip3d2_l1->Fill(Electron_ip3d[ProLep.at(1).ind]);
	  h.sip3d2_l1->Fill(Electron_sip3d[ProLep.at(1).ind]);
	}
	else if(abs(ProLep.at(1).id)==13){  //ProLep.at(1) is muon
	  h.Iso2_l1->Fill(Muon_pfRelIso03_all[ProLep.at(1).ind]);
	  h.dxy2_l1->Fill(Muon_dxy[ProLep.at(1).ind]);
	  h.dz2_l1->Fill(Muon_dz[ProLep.at(1).ind]);
	  float a = Muon_dxy[ProLep.at(1).ind];
	  float b = Muon_dz[ProLep.at(1).ind];
	  float c = ((a*a)+(b*b));
	  h.c2_l1->Fill(c);
	  h.ip3d2_l1->Fill(Muon_ip3d[ProLep.at(1).ind]);
	  h.sip3d2_l1->Fill(Muon_sip3d[ProLep.at(1).ind]);
	}
	if(abs(DisLep.at(0).id)==11){ //for DisLep.at(0)  //DisLep.at(0) is electron
	  h.Iso2_d0->Fill(Electron_pfRelIso03_all[DisLep.at(0).ind]);
	  h.dxy2_d0->Fill(Electron_dxy[DisLep.at(0).ind]);
	  h.dz2_d0->Fill(Electron_dz[DisLep.at(0).ind]);
	  float a = Electron_dxy[DisLep.at(0).ind];
	  float b = Electron_dz[DisLep.at(0).ind];
	  float c = ((a*a)+(b*b));
	  h.c2_d0->Fill(c);
	  h.ip3d2_d0->Fill(Electron_ip3d[DisLep.at(0).ind]);
	  h.sip3d2_d0->Fill(Electron_sip3d[DisLep.at(0).ind]);
	}
	else if(abs(DisLep.at(0).id)==13){  //DisLep.at(0) is muon
	  h.Iso2_d0->Fill(Muon_pfRelIso03_all[DisLep.at(0).ind]);
	  h.dxy2_d0->Fill(Muon_dxy[DisLep.at(0).ind]);
	  h.dz2_d0->Fill(Muon_dz[DisLep.at(0).ind]);
	  float a = Muon_dxy[DisLep.at(0).ind];
	  float b = Muon_dz[DisLep.at(0).ind];
	  float c = ((a*a)+(b*b));
	  h.c2_d0->Fill(c);
	  h.ip3d2_d0->Fill(Muon_ip3d[DisLep.at(0).ind]);
	  h.sip3d2_d0->Fill(Muon_sip3d[DisLep.at(0).ind]);
	}

	//MET
	h.MET2->Fill(*MET_pt);

	//Transverse mass-l0
	float p = ProLep.at(0).v.Pt();
	float q = *MET_pt;
	float dphi = delta_phi(p,q);
	float mT = transv_mass(p,q,dphi);
	h.mT2_l0->Fill(mT);
	//Transverse mass-l1
	float a = ProLep.at(1).v.Pt();
	float b = *MET_pt;
	float Dphi = delta_phi(a,b);
	float mt = transv_mass(a,b,Dphi);
	h.mT2_l1->Fill(mt);
	//Transverse mass-d0
	float c = DisLep.at(0).v.Pt();
	float d = *MET_pt;
	float dPhi = delta_phi(c,d);
	float Mt = transv_mass(c,d,dPhi);
	h.mT2_d0->Fill(Mt);

	//dRmin between jet and d0,d1,d2
	float dRmin_l0j = 10000;
	float dRmin_l1j = 10000;
	float dRmin_d0j = 10000;
	if((int)goodJet.size()>0){
	  for(unsigned int i=0; i<goodJet.size(); i++){
	    float dR_l0j = goodJet.at(i).v.DeltaR(ProLep.at(0).v);
	    float dR_l1j = goodJet.at(i).v.DeltaR(ProLep.at(1).v);
	    float dR_d0j = goodJet.at(i).v.DeltaR(DisLep.at(0).v);
	    if(dR_l0j < dRmin_l0j){
	      dRmin_l0j = dR_l0j;
	    }
	    if(dR_l1j < dRmin_l1j){
	      dRmin_l1j = dR_l1j;
	    }
	    if(dR_d0j < dRmin_d0j){
	      dRmin_d0j = dR_d0j;
	    }
	  } 
	  h.dRmin2_l0j->Fill(dRmin_l0j);
	  h.dRmin2_l1j->Fill(dRmin_l1j);
	  h.dRmin2_d0j->Fill(dRmin_d0j);
	  
	  //HT
	  float HT=0;
	  for(unsigned int i=0; i<goodJet.size(); i++){
	    HT = HT + goodJet.at(i).v.Pt();
	  }
	  h.HT2->Fill(HT);
	} //if goodJet>0 ends
	
	if((int)goodJet.size()>=0){
	  h.njets2->Fill(goodJet.size());
	}

	//flavor classification
	int flav =10;
	if(abs(ProLep.at(0).id)==11 && abs(ProLep.at(1).id)==11 && abs(DisLep.at(0).id)==11){ 
	  flav = 0; //e e e
	}
	if(abs(ProLep.at(0).id)==11 && abs(ProLep.at(1).id)==11 && abs(DisLep.at(0).id)==13){
	  flav = 1; //e e mu 
	}
	if(abs(ProLep.at(0).id)==11 && abs(ProLep.at(1).id)==13 && abs(DisLep.at(0).id)==13){
	  flav = 2; //e mu mu
	}
	if(abs(ProLep.at(0).id)==13 && abs(ProLep.at(1).id)==13 && abs(DisLep.at(0).id)==13){
	  flav = 3; //mu mu mu
	}
	if(abs(ProLep.at(0).id)==13 && abs(ProLep.at(1).id)==13 && abs(DisLep.at(0).id)==11){
	  flav = 4; //mu mu e
	}
	if(abs(ProLep.at(0).id)==13 && abs(ProLep.at(1).id)==11 && abs(DisLep.at(0).id)==11){
	  flav = 5; //mu e e
	}
	if(abs(ProLep.at(0).id)==11 && abs(ProLep.at(1).id)==13 && abs(DisLep.at(0).id)==11){
	  flav = 6; //e mu e
	}
	if(abs(ProLep.at(0).id)==13 && abs(ProLep.at(1).id)==11 && abs(DisLep.at(0).id)==13){
	  flav = 7; //mu e mu
	}
	h.flav2->Fill(flav);
      } //if ProLep>1 && DisLep>0 ends
      
      
      
	//###########################################################
      
      
      
	//#####################2l1d - Z vetoed#############################
      
      
      
      //2l1d channel
      if((int)ProLep.size()>1 && (int)DisLep.size()>0){ //atleast 2 prompt and atleast 1 displaced
	if((((ProLep.at(0).v + ProLep.at(1).v).M()) < 76 or ((ProLep.at(0).v + ProLep.at(1).v).M()) > 106)  //mass_l0l1 <76 and mass_l0l1 > 106
	   && (((ProLep.at(0).v + ProLep.at(1).v + DisLep.at(0).v).M()) < 76 or ((ProLep.at(0).v + ProLep.at(1).v + DisLep.at(0).v).M()) > 106) && //mass_l0l1d0 < 76 and mass_l0l1d0 > 106
	   ((ProLep.at(1).v.DeltaR(DisLep.at(0).v))>0.4)){  //dR_l1d0 >0.4
       
	  //mass plots
	  h.m2_Zveto_l0l1->Fill((ProLep.at(0).v + ProLep.at(1).v).M());
	  h.m2_Zveto_l0d0->Fill((ProLep.at(0).v + DisLep.at(0).v).M());
	  h.m2_Zveto_l1d0->Fill((ProLep.at(1).v + DisLep.at(0).v).M());
	  h.m2_Zveto_l0l1d0->Fill((ProLep.at(0).v + ProLep.at(1).v + DisLep.at(0).v).M());
	  
	  //dPhi plots
	  h.dPhi2_Zveto_l0l1->Fill(ProLep.at(0).v.DeltaPhi(ProLep.at(1).v));
	  h.dPhi2_Zveto_l0d0->Fill(ProLep.at(0).v.DeltaPhi(DisLep.at(0).v));
	  h.dPhi2_Zveto_l1d0->Fill(ProLep.at(1).v.DeltaPhi(DisLep.at(0).v));
	  
	  //dPhi_MET_lepton
	  float dPhi2_Zveto_l0MET = delta_phi(ProLep.at(0).v.Pt(),*MET_pt);
	  h.dPhi2_Zveto_l0MET->Fill(dPhi2_Zveto_l0MET);
	  float dPhi2_Zveto_l1MET = delta_phi(ProLep.at(1).v.Pt(),*MET_pt);
	  h.dPhi2_Zveto_l1MET->Fill(dPhi2_Zveto_l1MET);
	  float dPhi2_Zveto_d0MET = delta_phi(DisLep.at(0).v.Pt(),*MET_pt);
	  h.dPhi2_Zveto_d0MET->Fill(dPhi2_Zveto_d0MET);
	  
	  //dR plots
	  h.dR2_Zveto_l0l1->Fill(ProLep.at(0).v.DeltaR(ProLep.at(1).v));
	  h.dR2_Zveto_l0d0->Fill(ProLep.at(0).v.DeltaR(DisLep.at(0).v));
	  h.dR2_Zveto_l1d0->Fill(ProLep.at(1).v.DeltaR(DisLep.at(0).v));

	  //pT plots
	  h.pT2_Zveto_l0->Fill(ProLep.at(0).v.Pt());
	  h.pT2_Zveto_l1->Fill(ProLep.at(1).v.Pt());
	  h.pT2_Zveto_d0->Fill(DisLep.at(0).v.Pt());
	  h.pT2_Zveto_l0l1->Fill((ProLep.at(0).v + ProLep.at(1).v).Pt());
	  h.pT2_Zveto_l0d0->Fill((ProLep.at(0).v + DisLep.at(0).v).Pt());
	  h.pT2_Zveto_l1d0->Fill((ProLep.at(1).v + DisLep.at(0).v).Pt());
	  h.pT2_Zveto_l0l1d0->Fill((ProLep.at(0).v + ProLep.at(1).v + DisLep.at(0).v).Pt());
	  
	  //(Isolation, dxy, dz, c, ip3d, sip3d)
	  if(abs(ProLep.at(0).id)==11){ //for ProLep.at(0)  //ProLep.at(0) is electron
	    h.Iso2_Zveto_l0->Fill(Electron_pfRelIso03_all[ProLep.at(0).ind]);
	    h.dxy2_Zveto_l0->Fill(Electron_dxy[ProLep.at(0).ind]);
	    h.dz2_Zveto_l0->Fill(Electron_dz[ProLep.at(0).ind]);
	    float A = Electron_dxy[ProLep.at(0).ind];
	    float B = Electron_dz[ProLep.at(0).ind];
	    float C = ((A*A)+(B*B));
	    h.c2_Zveto_l0->Fill(C);
	    h.ip3d2_Zveto_l0->Fill(Electron_ip3d[ProLep.at(0).ind]);
	    h.sip3d2_Zveto_l0->Fill(Electron_sip3d[ProLep.at(0).ind]);
	  }
	  else if(abs(ProLep.at(0).id)==13){  //ProLep.at(0) is muon
	    h.Iso2_Zveto_l0->Fill(Muon_pfRelIso03_all[ProLep.at(0).ind]);
	    h.dxy2_Zveto_l0->Fill(Muon_dxy[ProLep.at(0).ind]);
	    h.dz2_Zveto_l0->Fill(Muon_dz[ProLep.at(0).ind]);
	    float A = Muon_dxy[ProLep.at(0).ind];
	    float B = Muon_dz[ProLep.at(0).ind];
	    float C = ((A*A)+(B*B));
	    h.c2_Zveto_l0->Fill(C);
	    h.ip3d2_Zveto_l0->Fill(Muon_ip3d[ProLep.at(0).ind]);
	    h.sip3d2_Zveto_l0->Fill(Muon_sip3d[ProLep.at(0).ind]);
	  }
	  if(abs(ProLep.at(1).id)==11){ //for ProLep.at(1)  //ProLep.at(1) is electron
	    h.Iso2_Zveto_l1->Fill(Electron_pfRelIso03_all[ProLep.at(1).ind]);
	    h.dxy2_Zveto_l1->Fill(Electron_dxy[ProLep.at(1).ind]);
	    h.dz2_Zveto_l1->Fill(Electron_dz[ProLep.at(1).ind]);
	    float A = Electron_dxy[ProLep.at(1).ind];
	    float B = Electron_dz[ProLep.at(1).ind];
	    float C = ((A*A)+(B*B));
	    h.c2_Zveto_l1->Fill(C);
	    h.ip3d2_Zveto_l1->Fill(Electron_ip3d[ProLep.at(1).ind]);
	    h.sip3d2_Zveto_l1->Fill(Electron_sip3d[ProLep.at(1).ind]);
	  }
	  else if(abs(ProLep.at(1).id)==13){  //ProLep.at(1) is muon
	    h.Iso2_Zveto_l1->Fill(Muon_pfRelIso03_all[ProLep.at(1).ind]);
	    h.dxy2_Zveto_l1->Fill(Muon_dxy[ProLep.at(1).ind]);
	    h.dz2_Zveto_l1->Fill(Muon_dz[ProLep.at(1).ind]);
	    float A = Muon_dxy[ProLep.at(1).ind];
	    float B = Muon_dz[ProLep.at(1).ind];
	    float C = ((A*A)+(B*B));
	    h.c2_Zveto_l1->Fill(C);
	    h.ip3d2_Zveto_l1->Fill(Muon_ip3d[ProLep.at(1).ind]);
	    h.sip3d2_Zveto_l1->Fill(Muon_sip3d[ProLep.at(1).ind]);
	  }
	  if(abs(DisLep.at(0).id)==11){ //for DisLep.at(0)  //DisLep.at(0) is electron
	    h.Iso2_Zveto_d0->Fill(Electron_pfRelIso03_all[DisLep.at(0).ind]);
	    h.dxy2_Zveto_d0->Fill(Electron_dxy[DisLep.at(0).ind]);
	    h.dz2_Zveto_d0->Fill(Electron_dz[DisLep.at(0).ind]);
	    float A = Electron_dxy[DisLep.at(0).ind];
	    float B = Electron_dz[DisLep.at(0).ind];
	    float C = ((A*A)+(B*B));
	    h.c2_Zveto_d0->Fill(C);
	    h.ip3d2_Zveto_d0->Fill(Electron_ip3d[DisLep.at(0).ind]);
	    h.sip3d2_Zveto_d0->Fill(Electron_sip3d[DisLep.at(0).ind]);
	  }
	  else if(abs(DisLep.at(0).id)==13){  //DisLep.at(0) is muon
	    h.Iso2_Zveto_d0->Fill(Muon_pfRelIso03_all[DisLep.at(0).ind]);
	    h.dxy2_Zveto_d0->Fill(Muon_dxy[DisLep.at(0).ind]);
	    h.dz2_Zveto_d0->Fill(Muon_dz[DisLep.at(0).ind]);
	    float A = Muon_dxy[DisLep.at(0).ind];
	    float B = Muon_dz[DisLep.at(0).ind];
	    float C = ((A*A)+(B*B));
	    h.c2_Zveto_d0->Fill(C);
	    h.ip3d2_Zveto_d0->Fill(Muon_ip3d[DisLep.at(0).ind]);
	    h.sip3d2_Zveto_d0->Fill(Muon_sip3d[DisLep.at(0).ind]);
	  }
	  
	  //MET
	  h.MET2_Zveto->Fill(*MET_pt);
	  
	  //Transverse mass-l0
	  float P = ProLep.at(0).v.Pt();
	  float Q = *MET_pt;
	  float delPhi = delta_phi(P,Q);
	  float transvM = transv_mass(P,Q,delPhi);
	  h.mT2_Zveto_l0->Fill(transvM);
	  //Transverse mass-l1
	  float R = ProLep.at(1).v.Pt();
	  float S = *MET_pt;
	  float delPHi= delta_phi(R,S);
	  float TransvM = transv_mass(R,S,delPHi);
	  h.mT2_Zveto_l1->Fill(TransvM);
	  //Transverse mass-d0
	  float T = DisLep.at(0).v.Pt();
	  float U = *MET_pt;
	  float DelPhi = delta_phi(T,U);
	  float Transvm = transv_mass(T,U,DelPhi);
	  h.mT2_Zveto_d0->Fill(Transvm);
	  
	  //dRmin between jet and d0,d1,d2
	  float dRmin_Zveto_l0j = 10000;
	  float dRmin_Zveto_l1j = 10000;
	  float dRmin_Zveto_d0j = 10000;
	  if((int)goodJet.size()>0){
	    for(unsigned int i=0; i<goodJet.size(); i++){
	      float dR_Zveto_l0j = goodJet.at(i).v.DeltaR(ProLep.at(0).v);
	      float dR_Zveto_l1j = goodJet.at(i).v.DeltaR(ProLep.at(1).v);
	      float dR_Zveto_d0j = goodJet.at(i).v.DeltaR(DisLep.at(0).v);
	      if(dR_Zveto_l0j < dRmin_Zveto_l0j){
		dRmin_Zveto_l0j = dR_Zveto_l0j;
	      }
	      if(dR_Zveto_l1j < dRmin_Zveto_l1j){
		dRmin_Zveto_l1j = dR_Zveto_l1j;
	      }
	      if(dR_Zveto_d0j < dRmin_Zveto_d0j){
		dRmin_Zveto_d0j = dR_Zveto_d0j;
	      }
	    } 
	    h.dRmin2_Zveto_l0j->Fill(dRmin_Zveto_l0j);
	    h.dRmin2_Zveto_l1j->Fill(dRmin_Zveto_l1j);
	    h.dRmin2_Zveto_d0j->Fill(dRmin_Zveto_d0j);
	    
	    //HT
	    float HT_Zveto=0;
	    for(unsigned int i=0; i<goodJet.size(); i++){
	      HT_Zveto = HT_Zveto + goodJet.at(i).v.Pt();
	    }
	    h.HT2_Zveto->Fill(HT_Zveto);
	  } //if goodJet>0 ends
	  
	if((int)goodJet.size()>=0){
	  h.njets2_Zveto->Fill(goodJet.size());
	}
	  
	  //flavor classification
	  int Flav =10;
	  if(abs(ProLep.at(0).id)==11 && abs(ProLep.at(1).id)==11 && abs(DisLep.at(0).id)==11){ 
	    Flav = 0; //e e e
	  }
	  if(abs(ProLep.at(0).id)==11 && abs(ProLep.at(1).id)==11 && abs(DisLep.at(0).id)==13){
	    Flav = 1; //e e mu 
	  }
	  if(abs(ProLep.at(0).id)==11 && abs(ProLep.at(1).id)==13 && abs(DisLep.at(0).id)==13){
	    Flav = 2; //e mu mu
	  }
	  if(abs(ProLep.at(0).id)==13 && abs(ProLep.at(1).id)==13 && abs(DisLep.at(0).id)==13){
	    Flav = 3; //mu mu mu
	  }
	  if(abs(ProLep.at(0).id)==13 && abs(ProLep.at(1).id)==13 && abs(DisLep.at(0).id)==11){
	    Flav = 4; //mu mu e
	  }
	  if(abs(ProLep.at(0).id)==13 && abs(ProLep.at(1).id)==11 && abs(DisLep.at(0).id)==11){
	    Flav = 5; //mu e e
	  }
	  if(abs(ProLep.at(0).id)==11 && abs(ProLep.at(1).id)==13 && abs(DisLep.at(0).id)==11){
	    Flav = 6; //e mu e
	  }
	  if(abs(ProLep.at(0).id)==13 && abs(ProLep.at(1).id)==11 && abs(DisLep.at(0).id)==13){
	    Flav = 7; //mu e mu
	  }
	  h.flav2_Zveto->Fill(Flav);
	}	
      } //if ProLep>1 && DisLep>0 ends
      
      
      
	//###########################################################
      
      
      
	//#######################2l1d - CR###############################
      
      
      
      //2l1d channel
      if((int)ProLep.size()>1 && (int)DisLep.size()>0){ //atleast 2 prompt and atleast 1 displaced
	if((((ProLep.at(0).v + ProLep.at(1).v).M()) > 76 && ((ProLep.at(0).v + ProLep.at(1).v).M()) < 106) &&  //DiLep mass betn 76 and 106
	   (((ProLep.at(0).v + ProLep.at(1).v + DisLep.at(0).v).M()) > 76 && ((ProLep.at(0).v + ProLep.at(1).v + DisLep.at(0).v).M()) < 106) &&  //TriLep mass betn 76 and 106
	   (ProLep.at(0).id * ProLep.at(1).id)<0 &&  //oppositely charged prompt leptons
	   ((abs(ProLep.at(0).id) == 11 && abs(ProLep.at(1).id)==11) or (abs(ProLep.at(0).id) ==13 && abs(ProLep.at(1).id) ==13)) &&  //same flavour prompt leptons
	   (*MET_pt<40)) {  //MET <40
	  
	  //mass plots
	  h.m2_Z_CR_l0l1->Fill((ProLep.at(0).v + ProLep.at(1).v).M());
	  h.m2_Z_CR_l0d0->Fill((ProLep.at(0).v + DisLep.at(0).v).M());
	  h.m2_Z_CR_l1d0->Fill((ProLep.at(1).v + DisLep.at(0).v).M());
	  h.m2_Z_CR_l0l1d0->Fill((ProLep.at(0).v + ProLep.at(1).v + DisLep.at(0).v).M());
	  
	  //dPhi plots
	  h.dPhi2_Z_CR_l0l1->Fill(ProLep.at(0).v.DeltaPhi(ProLep.at(1).v));
	  h.dPhi2_Z_CR_l0d0->Fill(ProLep.at(0).v.DeltaPhi(DisLep.at(0).v));
	  h.dPhi2_Z_CR_l1d0->Fill(ProLep.at(1).v.DeltaPhi(DisLep.at(0).v));
	  
	  //dPhi_MET_lepton
	  float dPhi2_CR_l0MET = delta_phi(ProLep.at(0).v.Pt(),*MET_pt);
	  h.dPhi2_Z_CR_l0MET->Fill(dPhi2_CR_l0MET);
	  float dPhi2_CR_l1MET = delta_phi(ProLep.at(1).v.Pt(),*MET_pt);
	  h.dPhi2_Z_CR_l1MET->Fill(dPhi2_CR_l1MET);
	  float dPhi2_CR_d0MET = delta_phi(DisLep.at(0).v.Pt(),*MET_pt);
	  h.dPhi2_Z_CR_d0MET->Fill(dPhi2_CR_d0MET);
	  
	  //dR plots
	  h.dR2_Z_CR_l0l1->Fill(ProLep.at(0).v.DeltaR(ProLep.at(1).v));
	  h.dR2_Z_CR_l0d0->Fill(ProLep.at(0).v.DeltaR(DisLep.at(0).v));
	  h.dR2_Z_CR_l1d0->Fill(ProLep.at(1).v.DeltaR(DisLep.at(0).v));

	  //pT plots
	  h.pT2_Z_CR_l0->Fill(ProLep.at(0).v.Pt());
	  h.pT2_Z_CR_l1->Fill(ProLep.at(1).v.Pt());
	  h.pT2_Z_CR_d0->Fill(DisLep.at(0).v.Pt());
	  h.pT2_Z_CR_l0l1->Fill((ProLep.at(0).v + ProLep.at(1).v).Pt());
	  h.pT2_Z_CR_l0d0->Fill((ProLep.at(0).v + DisLep.at(0).v).Pt());
	  h.pT2_Z_CR_l1d0->Fill((ProLep.at(1).v + DisLep.at(0).v).Pt());
	  h.pT2_Z_CR_l0l1d0->Fill((ProLep.at(0).v + ProLep.at(1).v + DisLep.at(0).v).Pt());
	  
	  //(Isolation, dxy, dz, c, ip3d, sip3d)
	  if(abs(ProLep.at(0).id)==11){ //for ProLep.at(0)  //ProLep.at(0) is electron
	    h.Iso2_Z_CR_l0->Fill(Electron_pfRelIso03_all[ProLep.at(0).ind]);
	    h.dxy2_Z_CR_l0->Fill(Electron_dxy[ProLep.at(0).ind]);
	    h.dz2_Z_CR_l0->Fill(Electron_dz[ProLep.at(0).ind]);
	    float A = Electron_dxy[ProLep.at(0).ind];
	    float B = Electron_dz[ProLep.at(0).ind];
	    float C = ((A*A)+(B*B));
	    h.c2_Z_CR_l0->Fill(C);
	    h.ip3d2_Z_CR_l0->Fill(Electron_ip3d[ProLep.at(0).ind]);
	    h.sip3d2_Z_CR_l0->Fill(Electron_sip3d[ProLep.at(0).ind]);
	  }
	  else if(abs(ProLep.at(0).id)==13){  //ProLep.at(0) is muon
	    h.Iso2_Z_CR_l0->Fill(Muon_pfRelIso03_all[ProLep.at(0).ind]);
	    h.dxy2_Z_CR_l0->Fill(Muon_dxy[ProLep.at(0).ind]);
	    h.dz2_Z_CR_l0->Fill(Muon_dz[ProLep.at(0).ind]);
	    float A = Muon_dxy[ProLep.at(0).ind];
	    float B = Muon_dz[ProLep.at(0).ind];
	    float C = ((A*A)+(B*B));
	    h.c2_Z_CR_l0->Fill(C);
	    h.ip3d2_Z_CR_l0->Fill(Muon_ip3d[ProLep.at(0).ind]);
	    h.sip3d2_Z_CR_l0->Fill(Muon_sip3d[ProLep.at(0).ind]);
	  }
	  if(abs(ProLep.at(1).id)==11){ //for ProLep.at(1)  //ProLep.at(1) is electron
	    h.Iso2_Z_CR_l1->Fill(Electron_pfRelIso03_all[ProLep.at(1).ind]);
	    h.dxy2_Z_CR_l1->Fill(Electron_dxy[ProLep.at(1).ind]);
	    h.dz2_Z_CR_l1->Fill(Electron_dz[ProLep.at(1).ind]);
	    float A = Electron_dxy[ProLep.at(1).ind];
	    float B = Electron_dz[ProLep.at(1).ind];
	    float C = ((A*A)+(B*B));
	    h.c2_Z_CR_l1->Fill(C);
	    h.ip3d2_Z_CR_l1->Fill(Electron_ip3d[ProLep.at(1).ind]);
	    h.sip3d2_Z_CR_l1->Fill(Electron_sip3d[ProLep.at(1).ind]);
	  }
	  else if(abs(ProLep.at(1).id)==13){  //ProLep.at(1) is muon
	    h.Iso2_Z_CR_l1->Fill(Muon_pfRelIso03_all[ProLep.at(1).ind]);
	    h.dxy2_Z_CR_l1->Fill(Muon_dxy[ProLep.at(1).ind]);
	    h.dz2_Z_CR_l1->Fill(Muon_dz[ProLep.at(1).ind]);
	    float A = Muon_dxy[ProLep.at(1).ind];
	    float B = Muon_dz[ProLep.at(1).ind];
	    float C = ((A*A)+(B*B));
	    h.c2_Z_CR_l1->Fill(C);
	    h.ip3d2_Z_CR_l1->Fill(Muon_ip3d[ProLep.at(1).ind]);
	    h.sip3d2_Z_CR_l1->Fill(Muon_sip3d[ProLep.at(1).ind]);
	  }
	  if(abs(DisLep.at(0).id)==11){ //for DisLep.at(0)  //DisLep.at(0) is electron
	    h.Iso2_Z_CR_d0->Fill(Electron_pfRelIso03_all[DisLep.at(0).ind]);
	    h.dxy2_Z_CR_d0->Fill(Electron_dxy[DisLep.at(0).ind]);
	    h.dz2_Z_CR_d0->Fill(Electron_dz[DisLep.at(0).ind]);
	    float A = Electron_dxy[DisLep.at(0).ind];
	    float B = Electron_dz[DisLep.at(0).ind];
	    float C = ((A*A)+(B*B));
	    h.c2_Z_CR_d0->Fill(C);
	    h.ip3d2_Z_CR_d0->Fill(Electron_ip3d[DisLep.at(0).ind]);
	    h.sip3d2_Z_CR_d0->Fill(Electron_sip3d[DisLep.at(0).ind]);
	  }
	  else if(abs(DisLep.at(0).id)==13){  //DisLep.at(0) is muon
	    h.Iso2_Z_CR_d0->Fill(Muon_pfRelIso03_all[DisLep.at(0).ind]);
	    h.dxy2_Z_CR_d0->Fill(Muon_dxy[DisLep.at(0).ind]);
	    h.dz2_Z_CR_d0->Fill(Muon_dz[DisLep.at(0).ind]);
	    float A = Muon_dxy[DisLep.at(0).ind];
	    float B = Muon_dz[DisLep.at(0).ind];
	    float C = ((A*A)+(B*B));
	    h.c2_Z_CR_d0->Fill(C);
	    h.ip3d2_Z_CR_d0->Fill(Muon_ip3d[DisLep.at(0).ind]);
	    h.sip3d2_Z_CR_d0->Fill(Muon_sip3d[DisLep.at(0).ind]);
	  }
	  
	  //MET
	  h.MET2_Z_CR->Fill(*MET_pt);
	  
	  //Transverse mass-l0
	  float P = ProLep.at(0).v.Pt();
	  float Q = *MET_pt;
	  float delPhi = delta_phi(P,Q);
	  float transvM = transv_mass(P,Q,delPhi);
	  h.mT2_Z_CR_l0->Fill(transvM);
	  //Transverse mass-l1
	  float R = ProLep.at(1).v.Pt();
	  float S = *MET_pt;
	  float delPHi= delta_phi(R,S);
	  float TransvM = transv_mass(R,S,delPHi);
	  h.mT2_Z_CR_l1->Fill(TransvM);
	  //Transverse mass-d0
	  float T = DisLep.at(0).v.Pt();
	  float U = *MET_pt;
	  float DelPhi = delta_phi(T,U);
	  float Transvm = transv_mass(T,U,DelPhi);
	  h.mT2_Z_CR_d0->Fill(Transvm);
	  
	  //dRmin between jet and d0,d1,d2
	  float dRmin_CR_l0j = 10000;
	  float dRmin_CR_l1j = 10000;
	  float dRmin_CR_d0j = 10000;
	  if((int)goodJet.size()>0){
	    for(unsigned int i=0; i<goodJet.size(); i++){
	      float dR_CR_l0j = goodJet.at(i).v.DeltaR(ProLep.at(0).v);
	      float dR_CR_l1j = goodJet.at(i).v.DeltaR(ProLep.at(1).v);
	      float dR_CR_d0j = goodJet.at(i).v.DeltaR(DisLep.at(0).v);
	      if(dR_CR_l0j < dRmin_CR_l0j){
		dRmin_CR_l0j = dR_CR_l0j;
	      }
	      if(dR_CR_l1j < dRmin_CR_l1j){
		dRmin_CR_l1j = dR_CR_l1j;
	      }
	      if(dR_CR_d0j < dRmin_CR_d0j){
		dRmin_CR_d0j = dR_CR_d0j;
	      }
	    } 
	    h.dRmin2_Z_CR_l0j->Fill(dRmin_CR_l0j);
	    h.dRmin2_Z_CR_l1j->Fill(dRmin_CR_l1j);
	    h.dRmin2_Z_CR_d0j->Fill(dRmin_CR_d0j);
	    
	    //HT
	    float HT_CR=0;
	    for(unsigned int i=0; i<goodJet.size(); i++){
	      HT_CR = HT_CR + goodJet.at(i).v.Pt();
	    }
	    h.HT2_Z_CR->Fill(HT_CR);
	  } //if goodJet>0 ends

	if((int)goodJet.size()>=0){
	  h.njets2_Z_CR->Fill(goodJet.size());
	}
	  
	  //flavor classification
	  int Flav =10;
	  if(abs(ProLep.at(0).id)==11 && abs(ProLep.at(1).id)==11 && abs(DisLep.at(0).id)==11){ 
	    Flav = 0; //e e e
	  }
	  if(abs(ProLep.at(0).id)==11 && abs(ProLep.at(1).id)==11 && abs(DisLep.at(0).id)==13){
	    Flav = 1; //e e mu 
	  }
	  if(abs(ProLep.at(0).id)==11 && abs(ProLep.at(1).id)==13 && abs(DisLep.at(0).id)==13){
	    Flav = 2; //e mu mu
	  }
	  if(abs(ProLep.at(0).id)==13 && abs(ProLep.at(1).id)==13 && abs(DisLep.at(0).id)==13){
	    Flav = 3; //mu mu mu
	  }
	  if(abs(ProLep.at(0).id)==13 && abs(ProLep.at(1).id)==13 && abs(DisLep.at(0).id)==11){
	    Flav = 4; //mu mu e
	  }
	  if(abs(ProLep.at(0).id)==13 && abs(ProLep.at(1).id)==11 && abs(DisLep.at(0).id)==11){
	    Flav = 5; //mu e e
	  }
	  if(abs(ProLep.at(0).id)==11 && abs(ProLep.at(1).id)==13 && abs(DisLep.at(0).id)==11){
	    Flav = 6; //e mu e
	  }
	  if(abs(ProLep.at(0).id)==13 && abs(ProLep.at(1).id)==11 && abs(DisLep.at(0).id)==13){
	    Flav = 7; //mu e mu
	  }
	  h.flav2_Z_CR->Fill(Flav);
	}	
      } //if ProLep>1 && DisLep>0 ends
      
      
      
      //############################################################
      
      
      
      //################### 2l1d channel - m_l0l1 <12 ########################
      
      
      
      //2l1d channel
      if((int)ProLep.size()>1 && (int)DisLep.size()>0 && ((ProLep.at(0).v + ProLep.at(1).v).M()<12)){ //atleast 2 prompt and atleast 1 displaced
	// " _12_ " => m_l0l1<12 cut

	//mass plots
	h.m2_12_l0l1->Fill((ProLep.at(0).v + ProLep.at(1).v).M());
	h.m2_12_l0d0->Fill((ProLep.at(0).v + DisLep.at(0).v).M());
	h.m2_12_l1d0->Fill((ProLep.at(1).v + DisLep.at(0).v).M());
	h.m2_12_l0l1d0->Fill((ProLep.at(0).v + ProLep.at(1).v + DisLep.at(0).v).M());
	
	//dPhi plots
	h.dPhi2_12_l0l1->Fill(ProLep.at(0).v.DeltaPhi(ProLep.at(1).v));
	h.dPhi2_12_l0d0->Fill(ProLep.at(0).v.DeltaPhi(DisLep.at(0).v));
	h.dPhi2_12_l1d0->Fill(ProLep.at(1).v.DeltaPhi(DisLep.at(0).v));
	
	//dR plots
	h.dR2_12_l0l1->Fill(ProLep.at(0).v.DeltaR(ProLep.at(1).v));
	h.dR2_12_l0d0->Fill(ProLep.at(0).v.DeltaR(DisLep.at(0).v));
	h.dR2_12_l1d0->Fill(ProLep.at(1).v.DeltaR(DisLep.at(0).v));
	
	//pT plots
	h.pT2_12_l0->Fill(ProLep.at(0).v.Pt());
	h.pT2_12_l1->Fill(ProLep.at(1).v.Pt());
	h.pT2_12_d0->Fill(DisLep.at(0).v.Pt());
	
	//MET
	h.MET2_12->Fill(*MET_pt);
	
	//(Isolation, dxy, dz, c, ip3d, sip3d)
	if(abs(ProLep.at(0).id)==11){ //for ProLep.at(0)  //ProLep.at(0) is electron
	  h.Iso2_12_l0->Fill(Electron_pfRelIso03_all[ProLep.at(0).ind]);
	  h.dxy2_12_l0->Fill(Electron_dxy[ProLep.at(0).ind]);
	  h.dz2_12_l0->Fill(Electron_dz[ProLep.at(0).ind]);
	  h.ip3d2_12_l0->Fill(Electron_ip3d[ProLep.at(0).ind]);
	  h.sip3d2_12_l0->Fill(Electron_sip3d[ProLep.at(0).ind]);
	}
	else if(abs(ProLep.at(0).id)==13){  //ProLep.at(0) is muon
	  h.Iso2_12_l0->Fill(Muon_pfRelIso03_all[ProLep.at(0).ind]);
	  h.dxy2_12_l0->Fill(Muon_dxy[ProLep.at(0).ind]);
	  h.dz2_12_l0->Fill(Muon_dz[ProLep.at(0).ind]);
	  h.ip3d2_12_l0->Fill(Muon_ip3d[ProLep.at(0).ind]);
	  h.sip3d2_12_l0->Fill(Muon_sip3d[ProLep.at(0).ind]);
	}
	if(abs(ProLep.at(1).id)==11){ //for ProLep.at(1)  //ProLep.at(1) is electron
	  h.Iso2_12_l1->Fill(Electron_pfRelIso03_all[ProLep.at(1).ind]);
	  h.dxy2_12_l1->Fill(Electron_dxy[ProLep.at(1).ind]);
	  h.dz2_12_l1->Fill(Electron_dz[ProLep.at(1).ind]);
	  h.ip3d2_12_l1->Fill(Electron_ip3d[ProLep.at(1).ind]);
	  h.sip3d2_12_l1->Fill(Electron_sip3d[ProLep.at(1).ind]);
	}
	else if(abs(ProLep.at(1).id)==13){  //ProLep.at(1) is muon
	  h.Iso2_12_l1->Fill(Muon_pfRelIso03_all[ProLep.at(1).ind]);
	  h.dxy2_12_l1->Fill(Muon_dxy[ProLep.at(1).ind]);
	  h.dz2_12_l1->Fill(Muon_dz[ProLep.at(1).ind]);
	  h.ip3d2_12_l1->Fill(Muon_ip3d[ProLep.at(1).ind]);
	  h.sip3d2_12_l1->Fill(Muon_sip3d[ProLep.at(1).ind]);
	}
	if(abs(DisLep.at(0).id)==11){ //for DisLep.at(0)  //DisLep.at(0) is electron
	  h.Iso2_12_d0->Fill(Electron_pfRelIso03_all[DisLep.at(0).ind]);
	  h.dxy2_12_d0->Fill(Electron_dxy[DisLep.at(0).ind]);
	  h.dz2_12_d0->Fill(Electron_dz[DisLep.at(0).ind]);
	  h.ip3d2_12_d0->Fill(Electron_ip3d[DisLep.at(0).ind]);
	  h.sip3d2_12_d0->Fill(Electron_sip3d[DisLep.at(0).ind]);
	}
	else if(abs(DisLep.at(0).id)==13){  //DisLep.at(0) is muon
	  h.Iso2_12_d0->Fill(Muon_pfRelIso03_all[DisLep.at(0).ind]);
	  h.dxy2_12_d0->Fill(Muon_dxy[DisLep.at(0).ind]);
	  h.dz2_12_d0->Fill(Muon_dz[DisLep.at(0).ind]);
	  h.ip3d2_12_d0->Fill(Muon_ip3d[DisLep.at(0).ind]);
	  h.sip3d2_12_d0->Fill(Muon_sip3d[DisLep.at(0).ind]);
	}
	
	//Transverse mass-l0
	float p = ProLep.at(0).v.Pt();
	float q = *MET_pt;
	float dphi = delta_phi(p,q);
	float mT = transv_mass(p,q,dphi);
	h.mT2_12_l0->Fill(mT);
	//Transverse mass-l1
	float a = ProLep.at(1).v.Pt();
	float b = *MET_pt;
	float Dphi = delta_phi(a,b);
	float mt = transv_mass(a,b,Dphi);
	h.mT2_12_l1->Fill(mt);
	//Transverse mass-d0
	float c = DisLep.at(0).v.Pt();
	float d = *MET_pt;
	float dPhi = delta_phi(c,d);
	float Mt = transv_mass(c,d,dPhi);
	h.mT2_12_d0->Fill(Mt);
	
	if((int)goodJet.size()>=0){
	  h.njets2_12->Fill(goodJet.size());
	}
	
	//flavor classification
	int flav =10;
	if(abs(ProLep.at(0).id)==11 && abs(ProLep.at(1).id)==11 && abs(DisLep.at(0).id)==11){ 
	  flav = 0; //e e e
	}
	if(abs(ProLep.at(0).id)==11 && abs(ProLep.at(1).id)==11 && abs(DisLep.at(0).id)==13){
	  flav = 1; //e e mu 
	}
	if(abs(ProLep.at(0).id)==11 && abs(ProLep.at(1).id)==13 && abs(DisLep.at(0).id)==13){
	  flav = 2; //e mu mu
	}
	if(abs(ProLep.at(0).id)==13 && abs(ProLep.at(1).id)==13 && abs(DisLep.at(0).id)==13){
	  flav = 3; //mu mu mu
	}
	if(abs(ProLep.at(0).id)==13 && abs(ProLep.at(1).id)==13 && abs(DisLep.at(0).id)==11){
	  flav = 4; //mu mu e
	}
	if(abs(ProLep.at(0).id)==13 && abs(ProLep.at(1).id)==11 && abs(DisLep.at(0).id)==11){
	  flav = 5; //mu e e
	}
	if(abs(ProLep.at(0).id)==11 && abs(ProLep.at(1).id)==13 && abs(DisLep.at(0).id)==11){
	  flav = 6; //e mu e
	}
	if(abs(ProLep.at(0).id)==13 && abs(ProLep.at(1).id)==11 && abs(DisLep.at(0).id)==13){
	  flav = 7; //mu e mu
	}
	h.flav2_12->Fill(flav);
      } //if ProLep>1 && DisLep>0 ends
      
      
      
      //############################################################
      
      
      
      //################## 2l1d channel - l0 l1 e, m_l0l1>12 ######################
      
      
      
      //2l1d channel
      if((int)ProLep.size()>1 && (int)DisLep.size()>0 && ((ProLep.at(0).v + ProLep.at(1).v).M()>12) && (abs(DisLep.at(0).id) == 11)){ //atleast 2 prompt and atleast 1 displaced
	// "_e_" => Displaced Lepton is electron
	
	//mass plots
	h.m2_e_l0l1->Fill((ProLep.at(0).v + ProLep.at(1).v).M());
	h.m2_e_l0d0->Fill((ProLep.at(0).v + DisLep.at(0).v).M());
	h.m2_e_l1d0->Fill((ProLep.at(1).v + DisLep.at(0).v).M());
	h.m2_e_l0l1d0->Fill((ProLep.at(0).v + ProLep.at(1).v + DisLep.at(0).v).M());
	
	//dPhi plots
	h.dPhi2_e_l0l1->Fill(ProLep.at(0).v.DeltaPhi(ProLep.at(1).v));
	h.dPhi2_e_l0d0->Fill(ProLep.at(0).v.DeltaPhi(DisLep.at(0).v));
	h.dPhi2_e_l1d0->Fill(ProLep.at(1).v.DeltaPhi(DisLep.at(0).v));
	
	//dR plots
	h.dR2_e_l0l1->Fill(ProLep.at(0).v.DeltaR(ProLep.at(1).v));
	h.dR2_e_l0d0->Fill(ProLep.at(0).v.DeltaR(DisLep.at(0).v));
	h.dR2_e_l1d0->Fill(ProLep.at(1).v.DeltaR(DisLep.at(0).v));
	
	//pT plots
	h.pT2_e_l0->Fill(ProLep.at(0).v.Pt());
	h.pT2_e_l1->Fill(ProLep.at(1).v.Pt());
	h.pT2_e_d0->Fill(DisLep.at(0).v.Pt());
	
	//MET
	h.MET2_e->Fill(*MET_pt);

	//isolation
	if(abs(ProLep.at(0).id)==11){
	  h.Iso2_e_l0->Fill(Electron_pfRelIso03_all[ProLep.at(0).ind]);
	}
	else if(abs(ProLep.at(0).id)==13){
	  h.Iso2_e_l0->Fill(Muon_pfRelIso03_all[ProLep.at(0).ind]);
	}
	if(abs(ProLep.at(1).id)==11){
	  h.Iso2_e_l1->Fill(Electron_pfRelIso03_all[ProLep.at(1).ind]);
	}
	else if(abs(ProLep.at(1).id)==13){
	  h.Iso2_e_l1->Fill(Muon_pfRelIso03_all[ProLep.at(1).ind]);
	}
	h.Iso2_e_d0->Fill(Electron_pfRelIso03_all[DisLep.at(0).ind]);

	//dxy, dz, ip3d, sip3d
	h.dxy2_e_d0->Fill(Electron_dxy[DisLep.at(0).ind]);
	h.dz2_e_d0->Fill(Electron_dz[DisLep.at(0).ind]);
	h.ip3d2_e_d0->Fill(Electron_ip3d[DisLep.at(0).ind]);
	h.sip3d2_e_d0->Fill(Electron_sip3d[DisLep.at(0).ind]);
	
	//Transverse mass-l0
	float p = ProLep.at(0).v.Pt();
	float q = *MET_pt;
	float dphi = delta_phi(p,q);
	float mT = transv_mass(p,q,dphi);
	h.mT2_e_l0->Fill(mT);
	//Transverse mass-l1
	float a = ProLep.at(1).v.Pt();
	float b = *MET_pt;
	float Dphi = delta_phi(a,b);
	float mt = transv_mass(a,b,Dphi);
	h.mT2_e_l1->Fill(mt);
	//Transverse mass-d0
	float c = DisLep.at(0).v.Pt();
	float d = *MET_pt;
	float dPhi = delta_phi(c,d);
	float Mt = transv_mass(c,d,dPhi);
	h.mT2_e_d0->Fill(Mt);
	
	if((int)goodJet.size()>=0){
	  h.njets2_e->Fill(goodJet.size());
	}
	
	//flavor classification
	int flav =10;
	if(abs(ProLep.at(0).id)==11 && abs(ProLep.at(1).id)==11 && abs(DisLep.at(0).id)==11){ 
	  flav = 0; //e e e
	}
	if(abs(ProLep.at(0).id)==11 && abs(ProLep.at(1).id)==11 && abs(DisLep.at(0).id)==13){
	  flav = 1; //e e mu 
	}
	if(abs(ProLep.at(0).id)==11 && abs(ProLep.at(1).id)==13 && abs(DisLep.at(0).id)==13){
	  flav = 2; //e mu mu
	}
	if(abs(ProLep.at(0).id)==13 && abs(ProLep.at(1).id)==13 && abs(DisLep.at(0).id)==13){
	  flav = 3; //mu mu mu
	}
	if(abs(ProLep.at(0).id)==13 && abs(ProLep.at(1).id)==13 && abs(DisLep.at(0).id)==11){
	  flav = 4; //mu mu e
	}
	if(abs(ProLep.at(0).id)==13 && abs(ProLep.at(1).id)==11 && abs(DisLep.at(0).id)==11){
	  flav = 5; //mu e e
	}
	if(abs(ProLep.at(0).id)==11 && abs(ProLep.at(1).id)==13 && abs(DisLep.at(0).id)==11){
	  flav = 6; //e mu e
	}
	if(abs(ProLep.at(0).id)==13 && abs(ProLep.at(1).id)==11 && abs(DisLep.at(0).id)==13){
	  flav = 7; //mu e mu
	}
	h.flav2_e->Fill(flav);
      } //if ProLep>1 && DisLep>0 ends
      
      
      
      //############################################################
      
      
      
      //################# 2l1d channel - l0 l1 mu, m_l0l1>12 #####################
      
      
      
      //2l1d channel
      if((int)ProLep.size()>1 && (int)DisLep.size()>0 && ((ProLep.at(0).v + ProLep.at(1).v).M()>12) && (abs(DisLep.at(0).id) == 13)){ //atleast 2 prompt and atleast 1 displaced
	// "_mu_" => Displaced Lepton is muon
	
	//mass plots
	h.m2_mu_l0l1->Fill((ProLep.at(0).v + ProLep.at(1).v).M());
	h.m2_mu_l0d0->Fill((ProLep.at(0).v + DisLep.at(0).v).M());
	h.m2_mu_l1d0->Fill((ProLep.at(1).v + DisLep.at(0).v).M());
	h.m2_mu_l0l1d0->Fill((ProLep.at(0).v + ProLep.at(1).v + DisLep.at(0).v).M());
	
	//dPhi plots
	h.dPhi2_mu_l0l1->Fill(ProLep.at(0).v.DeltaPhi(ProLep.at(1).v));
	h.dPhi2_mu_l0d0->Fill(ProLep.at(0).v.DeltaPhi(DisLep.at(0).v));
	h.dPhi2_mu_l1d0->Fill(ProLep.at(1).v.DeltaPhi(DisLep.at(0).v));
	
	//dR plots
	h.dR2_mu_l0l1->Fill(ProLep.at(0).v.DeltaR(ProLep.at(1).v));
	h.dR2_mu_l0d0->Fill(ProLep.at(0).v.DeltaR(DisLep.at(0).v));
	h.dR2_mu_l1d0->Fill(ProLep.at(1).v.DeltaR(DisLep.at(0).v));
	
	//pT plots
	h.pT2_mu_l0->Fill(ProLep.at(0).v.Pt());
	h.pT2_mu_l1->Fill(ProLep.at(1).v.Pt());
	h.pT2_mu_d0->Fill(DisLep.at(0).v.Pt());
	
	//MET
	h.MET2_mu->Fill(*MET_pt);

	//isolation
	if(abs(ProLep.at(0).id)==11){
	  h.Iso2_mu_l0->Fill(Electron_pfRelIso03_all[ProLep.at(0).ind]);
	}
	else if(abs(ProLep.at(0).id)==13){
	  h.Iso2_mu_l0->Fill(Muon_pfRelIso03_all[ProLep.at(0).ind]);
	}
	if(abs(ProLep.at(1).id)==11){
	  h.Iso2_mu_l1->Fill(Electron_pfRelIso03_all[ProLep.at(1).ind]);
	}
	else if(abs(ProLep.at(1).id)==13){
	  h.Iso2_mu_l1->Fill(Muon_pfRelIso03_all[ProLep.at(1).ind]);
	}
	h.Iso2_mu_d0->Fill(Muon_pfRelIso03_all[DisLep.at(0).ind]);

	//dxy, dz, ip3d, sip3d
	h.dxy2_mu_d0->Fill(Muon_dxy[DisLep.at(0).ind]);
	h.dz2_mu_d0->Fill(Muon_dz[DisLep.at(0).ind]);
	h.ip3d2_mu_d0->Fill(Muon_ip3d[DisLep.at(0).ind]);
	h.sip3d2_mu_d0->Fill(Muon_sip3d[DisLep.at(0).ind]);
	
	//Transverse mass-l0
	float p = ProLep.at(0).v.Pt();
	float q = *MET_pt;
	float dphi = delta_phi(p,q);
	float mT = transv_mass(p,q,dphi);
	h.mT2_mu_l0->Fill(mT);
	//Transverse mass-l1
	float a = ProLep.at(1).v.Pt();
	float b = *MET_pt;
	float Dphi = delta_phi(a,b);
	float mt = transv_mass(a,b,Dphi);
	h.mT2_mu_l1->Fill(mt);
	//Transverse mass-d0
	float c = DisLep.at(0).v.Pt();
	float d = *MET_pt;
	float dPhi = delta_phi(c,d);
	float Mt = transv_mass(c,d,dPhi);
	h.mT2_mu_d0->Fill(Mt);
	
	if((int)goodJet.size()>=0){
	  h.njets2_mu->Fill(goodJet.size());
	}
	
	//flavor classification
	int flav =10;
	if(abs(ProLep.at(0).id)==11 && abs(ProLep.at(1).id)==11 && abs(DisLep.at(0).id)==11){ 
	  flav = 0; //e e e
	}
	if(abs(ProLep.at(0).id)==11 && abs(ProLep.at(1).id)==11 && abs(DisLep.at(0).id)==13){
	  flav = 1; //e e mu 
	}
	if(abs(ProLep.at(0).id)==11 && abs(ProLep.at(1).id)==13 && abs(DisLep.at(0).id)==13){
	  flav = 2; //e mu mu
	}
	if(abs(ProLep.at(0).id)==13 && abs(ProLep.at(1).id)==13 && abs(DisLep.at(0).id)==13){
	  flav = 3; //mu mu mu
	}
	if(abs(ProLep.at(0).id)==13 && abs(ProLep.at(1).id)==13 && abs(DisLep.at(0).id)==11){
	  flav = 4; //mu mu e
	}
	if(abs(ProLep.at(0).id)==13 && abs(ProLep.at(1).id)==11 && abs(DisLep.at(0).id)==11){
	  flav = 5; //mu e e
	}
	if(abs(ProLep.at(0).id)==11 && abs(ProLep.at(1).id)==13 && abs(DisLep.at(0).id)==11){
	  flav = 6; //e mu e
	}
	if(abs(ProLep.at(0).id)==13 && abs(ProLep.at(1).id)==11 && abs(DisLep.at(0).id)==13){
	  flav = 7; //mu e mu
	}
	h.flav2_mu->Fill(flav);
      } //if ProLep>1 && DisLep>0 ends
      
      
      
      
    } //triggered_events ends
    
    
    
    
    
    //########### ANALYSIS ENDS HERE ##############
    } //triggerRes ends
  }//GoodEvt
  
  if(nEvtTotal>0){
    h.nEvt->Fill(nEvtTotal);
  }
 
    
  
  return kTRUE;
}


//######################################
//        USER DEFINED FUNCTIONS
//######################################


void nano9Ana::Sortpt(vector<Lepton> vec)
{
  
  for(int i=0; i<(int)vec.size()-1; i++){
    for(int j=i+1; j<(int)vec.size(); j++){
      if( vec[i].v.Pt() < vec[j].v.Pt() ) swap(vec.at(i), vec.at(j));
    }
  }
}

int nano9Ana::get_mother(int i)
{
  int pid = GenPart_pdgId[i];
  int motherid = GenPart_pdgId[GenPart_genPartIdxMother[i]];
  
  return motherid;
}


int nano9Ana::GenMother(int ind, int mom_ind)
{
  int p_id = GenPart_pdgId[ind];
  int m_id = GenPart_pdgId[mom_ind];
  while(p_id==m_id){
    ind = mom_ind;
    mom_ind = GenPart_genPartIdxMother[ind];
    p_id = GenPart_pdgId[ind];
    m_id = GenPart_pdgId[mom_ind];
  }
  return m_id;
}


float nano9Ana::delta_phi(float phi1, float phi2)
{
  //The correct deltaPhi falls in the interval [0 , pi]
  phi1 = TVector2::Phi_0_2pi(phi1);
  phi2 = TVector2::Phi_0_2pi(phi2);
  float dphi = fabs(phi1 - phi2);
  if(dphi>TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
  return dphi;
}

float nano9Ana::transv_mass(float E_lep, float MET, float dphi)
{
  //The inputs are the Energy of the lepton, MET and dPhi between the lepton and MET
  float mT = sqrt(2* E_lep * MET *(1-cos(dphi)));
  return mT;
}


void nano9Ana::BookHistograms()
{
  //The histograms are booked here.
  //Binning etc are done here.
  //These histograms are stored in the hst_<process name>.root file in the same order.

  //Example : new TH1F ("hst_name", "hst title", NBins, startVal, EndVal);

  h.nEvt = new TH1F("nEvt","nEvents",500, 0, 500);
  
  h.nmuons = new TH1F("nmuons", "number of muons", 10, 0, 10);
  h.mupt= new TH1F("mu_pt", "Muon p_{T}", 200, 0, 200);
  h.muprop[0] = new TH1F("mu0_pt","Leading muon p_{T}",200,0,200); h.muprop[0]->Sumw2();
  h.muprop[1] = new TH1F("mu0_phi","Leading muon phi",64,-3.2,3.2);h.muprop[1]->Sumw2();
  h.Mu0_dxy = new TH1F("Mu0_dxy","Mu0_dxy", 200, -1, 1);
  h.Mu0_dz = new TH1F("Mu0_dz","Mu0_dz", 200, -1, 1);
  h.Mu0_c = new TH1F("Mu0_c","Mu0_c", 200, -1, 1);
  h.Mu1_dxy = new TH1F("Mu1_dxy","Mu1_dxy", 200, -1, 1);
  h.Mu1_dz = new TH1F("Mu1_dz","Mu1_dz", 200, -1, 1);
  h.Mu1_c = new TH1F("Mu1_c","Mu1_c", 200, -1, 1);
  h.Mu2_dxy = new TH1F("Mu2_dxy","Mu2_dxy", 200, -1, 1);
  h.Mu2_dz = new TH1F("Mu2_dz","Mu2_dz", 200, -1, 1);
  h.Mu2_c = new TH1F("Mu2_c","Mu2_c", 200, -1, 1);
  h.Mu3_dxy = new TH1F("Mu3_dxy","Mu3_dxy", 200, -1, 1);
  h.Mu3_dz = new TH1F("Mu3_dz","Mu3_dz", 200, -1, 1);
  h.Mu3_c = new TH1F("Mu3_c","Mu3_c", 200, -1, 1);

  h.Ele0_dxy = new TH1F("Ele0_dxy","Ele0_dxy", 200, -1, 1);
  h.Ele0_dz = new TH1F("Ele0_dz","Ele0_dz", 200, -1, 1);
  h.Ele0_c = new TH1F("Ele0_c","Ele0_c", 200, -1, 1);
  h.Ele1_dxy = new TH1F("Ele1_dxy","Ele1_dxy", 200, -1, 1);
  h.Ele1_dz = new TH1F("Ele1_dz","Ele1_dz", 200, -1, 1);
  h.Ele1_c = new TH1F("Ele1_c","Ele1_c", 200, -1, 1);
  h.Ele2_dxy = new TH1F("Ele2_dxy","Ele2_dxy", 200, -1, 1);
  h.Ele2_dz = new TH1F("Ele2_dz","Ele2_dz", 200, -1, 1);
  h.Ele2_c = new TH1F("Ele2_c","Ele2_c", 200, -1, 1);
  h.Ele3_dxy = new TH1F("Ele3_dxy","Ele3_dxy", 200, -1, 1);
  h.Ele3_dz = new TH1F("Ele3_dz","Ele3_dz", 200, -1, 1);
  h.Ele3_c = new TH1F("Ele3_c","Ele3_c", 200, -1, 1);
  
  if(_data==0){
    h.Mu_dxy = new TH1F("Mu_dxy","Mu_dxy",800, -1, 1);
    h.Mu_dz = new TH1F("Mu_dz","Mu_dz",800, -1, 1);
    h.Ele_dxy = new TH1F("Ele_dxy","Ele_dxy", 800, -1, 1);
    h.Ele_dz = new TH1F("Ele_dz","Ele_dz",800, -1, 1);
  }


  
  //#############################################################
  //#############################################################
  
  
  
  //#######################  1l2d channel  #############################
  
  
  
  //mass
  h.m0_l0d0 = new TH1F("1l2d_m_l0d0","1l2d_m_l0d0",200, 0, 200);  
  h.m0_d0d1 = new TH1F("1l2d_m_d0d1","1l2d_m_d0d1",200, 0, 200);
  h.m0_l0d1 = new TH1F("1l2d_m_l0d1","1l2d_m_l0d1",200, 0, 200);
  h.m0_l0d0d1 = new TH1F("1l2d_m_l0d0d1","1l2d_m_l0d0d1",200, 0, 200);

  //dPhi
  h.dPhi0_l0d0 = new TH1F("1l2d_dPhi_l0d0","1l2d_dPhi_l0d0",64, -3.2, 3.2);  
  h.dPhi0_d0d1 = new TH1F("1l2d_dPhi_d0d1","1l2d_dPhi_d0d1",64, -3.2, 3.2);
  h.dPhi0_l0d1 = new TH1F("1l2d_dPhi_l0d1","1l2d_dPhi_l0d1",64, -3.2, 3.2);

  //dPhi_lep_MET
  h.dPhi0_l0MET = new TH1F("1l2d_dPhi_l0MET","1l2d_dPhi_l0_MET",32, 0, 3.2);  
  h.dPhi0_d0MET = new TH1F("1l2d_dPhi_d0MET","1l2d_dPhi_d0_MET",32, 0, 3.2);
  h.dPhi0_d1MET = new TH1F("1l2d_dPhi_d1MET","1l2d_dPhi_d1_MET",32, 0, 3.2);

  //dR
  h.dR0_l0d0 = new TH1F("1l2d_dR_l0d0","1l2d_dR_l0d0",100, 0, 10);  
  h.dR0_d0d1 = new TH1F("1l2d_dR_d0d1","1l2d_dR_d0d1",100, 0, 10);
  h.dR0_l0d1 = new TH1F("1l2d_dR_l0d1","1l2d_dR_l0d1",100, 0, 10);

  //pT
  h.pT0_l0 = new TH1F("1l2d_pT_l0","1l2d_pT_l0",200, 0, 200);  
  h.pT0_d0 = new TH1F("1l2d_pT_d0","1l2d_pT_d0",200, 0, 200);
  h.pT0_d1 = new TH1F("1l2d_pT_d1","1l2d_pT_d1",200, 0, 200);
  h.pT0_l0d0 = new TH1F("1l2d_pT_l0d0","1l2d_pT_l0d0",200, 0, 200);
  h.pT0_d0d1 = new TH1F("1l2d_pT_d0d1","1l2d_pT_d0d1",200, 0, 200);
  h.pT0_l0d1 = new TH1F("1l2d_pT_l0d1","1l2d_pT_l0d1",200, 0, 200);
  h.pT0_l0d0d1 = new TH1F("1l2d_pT_l0d0d1","1l2d_pT_l0d0d1",200, 0, 200);

  //Isolation (Iso03)
  h.Iso0_l0 = new TH1F("1l2d_Iso03_l0","1l2d_Iso03_l0",800, 0, 2);  
  h.Iso0_d0 = new TH1F("1l2d_Iso03_d0","1l2d_Iso03_d0",800, 0, 2);
  h.Iso0_d1 = new TH1F("1l2d_Iso03_d1","1l2d_Iso03_d1",800, 0, 2);

  //(dxy, dz, c, ip3d, sip3d)_l0
  h.dxy0_l0 = new TH1F("1l2d_dxy_l0","1l2d_dxy_l0",200, -1, 1);
  h.dz0_l0 = new TH1F("1l2d_dz_l0","1l2d_dz_l0",200, -1, 1);
  h.c0_l0 = new TH1F("1l2d_c_l0","1l2d_sqrt(dxy^2 + dz^2)_l0",200, -1, 1);
  h.ip3d0_l0 = new TH1F("1l2d_ip3d_l0","1l2d_ip3d_l0",200, 0, 20);
  h.sip3d0_l0 = new TH1F("1l2d_sip3d_l0","1l2d_sip3d_l0",200, 0, 200);

  //(dxy, dz, c, ip3d, sip3d)_d0
  h.dxy0_d0 = new TH1F("1l2d_dxy_d0","1l2d_dxy_d0",200, -1, 1);
  h.dz0_d0 = new TH1F("1l2d_dz_d0","1l2d_dz_d0",200, -1, 1);
  h.c0_d0 = new TH1F("1l2d_c_d0","1l2d_sqrt(dxy^2 + dz^2)_d0",200, -1, 1);
  h.ip3d0_d0 = new TH1F("1l2d_ip3d_d0","1l2d_ip3d_d0",200, 0, 20);
  h.sip3d0_d0 = new TH1F("1l2d_sip3d_d0","1l2d_sip3d_d0",200, 0, 200);

  //(dxy, dz, c, ip3d, sip3d)_d1
  h.dxy0_d1 = new TH1F("1l2d_dxy_d1","1l2d_dxy_d1",200, -1, 1);
  h.dz0_d1 = new TH1F("1l2d_dz_d1","1l2d_dz_d1",200, -1, 1);
  h.c0_d1 = new TH1F("1l2d_c_d1","1l2d_sqrt(dxy^2 + dz^2)_d1",200, -1, 1);
  h.ip3d0_d1 = new TH1F("1l2d_ip3d_d1","1l2d_ip3d_d1",200, 0, 20);
  h.sip3d0_d1 = new TH1F("1l2d_sip3d_d1","1l2d_sip3d_d1",200, 0, 200);

  //MET
  h.MET0 = new TH1F("1l2d_MET","1l2d_MET", 200, 0, 200);  

  //Transverse mass (mT)
  h.mT0_l0 = new TH1F("1l2d_mT_l0","1l2d_mT_l0",200, 0, 200);  
  h.mT0_d0 = new TH1F("1l2d_mT_d0","1l2d_mT_d0",200, 0, 200);
  h.mT0_d1 = new TH1F("1l2d_mT_d1","1l2d_mT_d1",200, 0, 200);

  //dRmin_lep_jet
  h.dRmin0_l0j = new TH1F("1l2d_dRmin_l0j","1l2d_dRmin_l0j",500, 0, 5); 
  h.dRmin0_d0j = new TH1F("1l2d_dRmin_d0j","1l2d_dRmin_d0j",500, 0, 5);
  h.dRmin0_d1j = new TH1F("1l2d_dRmin_d1j","1l2d_dRmin_d1j",500, 0, 5);

  //HT
  h.HT0 = new TH1F("1l2d_HT","1l2d_HT",200,0,200);  

  //njets
  h.njets0 = new TH1F("1l2d_njets","1l2d_njets",10, 0, 10); 

  //flavor classification
  h.flav0 = new TH1F("1l2d_flav","1l2d_flav, 0:eee 1:ee#mu 2:e#mu#mu 3:#mu#mu#mu 4:#mu#mue 5:#muee 6:e#mue 7:#mue#mu",8, 0, 8); 


  
  //#############################################################
  

  
  //####################  1l2d channel - J/Psi veto #########################
  
  
  
  //mass
  h.m0_Jveto_l0d0 = new TH1F("1l2d_Jveto_m_l0d0","1l2d_Jveto_m_l0d0",200, 0, 200);  
  h.m0_Jveto_d0d1 = new TH1F("1l2d_Jveto_m_d0d1","1l2d_Jveto_m_d0d1",200, 0, 200);
  h.m0_Jveto_l0d1 = new TH1F("1l2d_Jveto_m_l0d1","1l2d_Jveto_m_l0d1",200, 0, 200);
  h.m0_Jveto_l0d0d1 = new TH1F("1l2d_Jveto_m_l0d0d1","1l2d_Jveto_m_l0d0d1",200, 0, 200);

  //dPhi
  h.dPhi0_Jveto_l0d0 = new TH1F("1l2d_Jveto_dPhi_l0d0","1l2d_Jveto_dPhi_l0d0",64, -3.2, 3.2);  
  h.dPhi0_Jveto_d0d1 = new TH1F("1l2d_Jveto_dPhi_d0d1","1l2d_Jveto_dPhi_d0d1",64, -3.2, 3.2);
  h.dPhi0_Jveto_l0d1 = new TH1F("1l2d_Jveto_dPhi_l0d1","1l2d_Jveto_dPhi_l0d1",64, -3.2, 3.2);

  //dPhi_lep_MET
  h.dPhi0_Jveto_l0MET = new TH1F("1l2d_Jveto_dPhi_l0MET","1l2d_Jveto_dPhi_l0_MET",32, 0, 3.2);  
  h.dPhi0_Jveto_d0MET = new TH1F("1l2d_Jveto_dPhi_d0MET","1l2d_Jveto_dPhi_d0_MET",32, 0, 3.2);
  h.dPhi0_Jveto_d1MET = new TH1F("1l2d_Jveto_dPhi_d1MET","1l2d_Jveto_dPhi_d1_MET",32, 0, 3.2);

  //dR
  h.dR0_Jveto_l0d0 = new TH1F("1l2d_Jveto_dR_l0d0","1l2d_Jveto_dR_l0d0",100, 0, 10);  
  h.dR0_Jveto_d0d1 = new TH1F("1l2d_Jveto_dR_d0d1","1l2d_Jveto_dR_d0d1",100, 0, 10);
  h.dR0_Jveto_l0d1 = new TH1F("1l2d_Jveto_dR_l0d1","1l2d_Jveto_dR_l0d1",100, 0, 10);

  //pT
  h.pT0_Jveto_l0 = new TH1F("1l2d_Jveto_pT_l0","1l2d_Jveto_pT_l0",200, 0, 200);  
  h.pT0_Jveto_d0 = new TH1F("1l2d_Jveto_pT_d0","1l2d_Jveto_pT_d0",200, 0, 200);
  h.pT0_Jveto_d1 = new TH1F("1l2d_Jveto_pT_d1","1l2d_Jveto_pT_d1",200, 0, 200);
  h.pT0_Jveto_l0d0 = new TH1F("1l2d_Jveto_pT_l0d0","1l2d_Jveto_pT_l0d0",200, 0, 200);
  h.pT0_Jveto_d0d1 = new TH1F("1l2d_Jveto_pT_d0d1","1l2d_Jveto_pT_d0d1",200, 0, 200);
  h.pT0_Jveto_l0d1 = new TH1F("1l2d_Jveto_pT_l0d1","1l2d_Jveto_pT_l0d1",200, 0, 200);
  h.pT0_Jveto_l0d0d1 = new TH1F("1l2d_Jveto_pT_l0d0d1","1l2d_Jveto_pT_l0d0d1",200, 0, 200);

  //Isolation (Iso03)
  h.Iso0_Jveto_l0 = new TH1F("1l2d_Jveto_Iso03_l0","1l2d_Jveto_Iso03_l0",800, 0, 2);  
  h.Iso0_Jveto_d0 = new TH1F("1l2d_Jveto_Iso03_d0","1l2d_Jveto_Iso03_d0",800, 0, 2);
  h.Iso0_Jveto_d1 = new TH1F("1l2d_Jveto_Iso03_d1","1l2d_Jveto_Iso03_d1",800, 0, 2);

  //(dxy, dz, c, ip3d, sip3d)_l0
  h.dxy0_Jveto_l0 = new TH1F("1l2d_Jveto_dxy_l0","1l2d_Jveto_dxy_l0",200, -1, 1);
  h.dz0_Jveto_l0 = new TH1F("1l2d_Jveto_dz_l0","1l2d_Jveto_dz_l0",200, -1, 1);
  h.c0_Jveto_l0 = new TH1F("1l2d_Jveto_c_l0","1l2d_Jveto_sqrt(dxy^2 + dz^2)_l0",200, -1, 1);
  h.ip3d0_Jveto_l0 = new TH1F("1l2d_Jveto_ip3d_l0","1l2d_Jveto_ip3d_l0",200, 0, 20);
  h.sip3d0_Jveto_l0 = new TH1F("1l2d_Jveto_sip3d_l0","1l2d_Jveto_sip3d_l0",200, 0, 200);

  //(dxy, dz, c, ip3d, sip3d)_d0
  h.dxy0_Jveto_d0 = new TH1F("1l2d_Jveto_dxy_d0","1l2d_Jveto_dxy_d0",200, -1, 1);
  h.dz0_Jveto_d0 = new TH1F("1l2d_Jveto_dz_d0","1l2d_Jveto_dz_d0",200, -1, 1);
  h.c0_Jveto_d0 = new TH1F("1l2d_Jveto_c_d0","1l2d_Jveto_sqrt(dxy^2 + dz^2)_d0",200, -1, 1);
  h.ip3d0_Jveto_d0 = new TH1F("1l2d_Jveto_ip3d_d0","1l2d_Jveto_ip3d_d0",200, 0, 20);
  h.sip3d0_Jveto_d0 = new TH1F("1l2d_Jveto_sip3d_d0","1l2d_Jveto_sip3d_d0",200, 0, 200);

  //(dxy, dz, c, ip3d, sip3d)_d1
  h.dxy0_Jveto_d1 = new TH1F("1l2d_Jveto_dxy_d1","1l2d_Jveto_dxy_d1",200, -1, 1);
  h.dz0_Jveto_d1 = new TH1F("1l2d_Jveto_dz_d1","1l2d_Jveto_dz_d1",200, -1, 1);
  h.c0_Jveto_d1 = new TH1F("1l2d_Jveto_c_d1","1l2d_Jveto_sqrt(dxy^2 + dz^2)_d1",200, -1, 1);
  h.ip3d0_Jveto_d1 = new TH1F("1l2d_Jveto_ip3d_d1","1l2d_Jveto_ip3d_d1",200, 0, 20);
  h.sip3d0_Jveto_d1 = new TH1F("1l2d_Jveto_sip3d_d1","1l2d_Jveto_sip3d_d1",200, 0, 200);

  //MET
  h.MET0_Jveto = new TH1F("1l2d_Jveto_MET","1l2d_Jveto_MET", 200, 0, 200);  

  //Transverse mass (mT)
  h.mT0_Jveto_l0 = new TH1F("1l2d_Jveto_mT_l0","1l2d_Jveto_mT_l0",200, 0, 200);  
  h.mT0_Jveto_d0 = new TH1F("1l2d_Jveto_mT_d0","1l2d_Jveto_mT_d0",200, 0, 200);
  h.mT0_Jveto_d1 = new TH1F("1l2d_Jveto_mT_d1","1l2d_Jveto_mT_d1",200, 0, 200);

  //dRmin_lep_jet
  h.dRmin0_Jveto_l0j = new TH1F("1l2d_Jveto_dRmin_l0j","1l2d_Jveto_dRmin_l0j",500, 0, 5); 
  h.dRmin0_Jveto_d0j = new TH1F("1l2d_Jveto_dRmin_d0j","1l2d_Jveto_dRmin_d0j",500, 0, 5);
  h.dRmin0_Jveto_d1j = new TH1F("1l2d_Jveto_dRmin_d1j","1l2d_Jveto_dRmin_d1j",500, 0, 5);

  //HT
  h.HT0_Jveto = new TH1F("1l2d_Jveto_HT","1l2d_Jveto_HT",200,0,200);  

  //njets
  h.njets0_Jveto = new TH1F("1l2d_Jveto_njets","1l2d_Jveto_njets",10, 0, 10); 

  //flavor classification
  h.flav0_Jveto = new TH1F("1l2d_Jveto_flav","1l2d_Jveto_flav, 0:eee 1:ee#mu 2:e#mu#mu 3:#mu#mu#mu 4:#mu#mue 5:#muee 6:e#mue 7:#mue#mu",8, 0, 8);
  
  
  
  //#############################################################
  
  
  
  //####################  1l2d channel - l0 e e  ##########################
  // "_ee_" => the two displaced leptons are electrons
  
  
  
  //mass
  h.m0_ee_l0d0 = new TH1F("1l2d_ee_m_l0d0","1l2d_m_l0d0 (DisLep.at(0) and DisLep.at(1) are electrons)", 500, 0, 500);  
  h.m0_ee_l0d1 = new TH1F("1l2d_ee_m_l0d1","1l2d_m_l0d1 (DisLep.at(0) and DisLep.at(1) are electrons)", 500, 0, 500);
  h.m0_ee_d0d1 = new TH1F("1l2d_ee_m_d0d1","1l2d_m_d0d1 (DisLep.at(0) and DisLep.at(1) are electrons)", 500, 0, 500);
  h.m0_ee_l0d0d1 = new TH1F("1l2d_ee_m_l0d0d1","1l2d_m_l0d0d1 (DisLep.at(0) and DisLep.at(1) are electrons)", 500, 0, 500);

  //dPhi
  h.dPhi0_ee_l0d0 = new TH1F("1l2d_ee_dPhi_l0d0","1l2d_dPhi_l0d0 (DisLep.at(0) and DisLep.at(1) are electrons)", 64, -3.2, 3.2);  
  h.dPhi0_ee_l0d1 = new TH1F("1l2d_ee_dPhi_l0d1","1l2d_dPhi_l0d1 (DisLep.at(0) and DisLep.at(1) are electrons)", 64, -3.2, 3.2);
  h.dPhi0_ee_d0d1 = new TH1F("1l2d_ee_dPhi_d0d1","1l2d_dPhi_d0d1 (DisLep.at(0) and DisLep.at(1) are electrons)", 64, -3.2, 3.2);
 
  //dR
  h.dR0_ee_l0d0 = new TH1F("1l2d_ee_dR_l0d0","1l2d_dR_l0d0 (DisLep.at(0) and DisLep.at(1) are electrons)", 100, 0, 10); 
  h.dR0_ee_l0d1 = new TH1F("1l2d_ee_dR_l0d1","1l2d_dR_l0d1 (DisLep.at(0) and DisLep.at(1) are electrons)", 100, 0, 10);
  h.dR0_ee_d0d1 = new TH1F("1l2d_ee_dR_d0d1","1l2d_dR_d0d1 (DisLep.at(0) and DisLep.at(1) are electrons)", 100, 0, 10);

  //pT
  h.pT0_ee_l0 = new TH1F("1l2d_ee_pT_l0","1l2d_pT_l0 (DisLep.at(0) and DisLep.at(1) are electrons)", 200, 0, 200);  
  h.pT0_ee_d0 = new TH1F("1l2d_ee_pT_d0","1l2d_pT_d0 (DisLep.at(0) and DisLep.at(1) are electrons)", 200, 0, 200);
  h.pT0_ee_d1 = new TH1F("1l2d_ee_pT_d1","1l2d_pT_d1 (DisLep.at(0) and DisLep.at(1) are electrons)", 200, 0, 200);
  
  //MET
  h.MET0_ee = new TH1F("1l2d_ee_MET","1l2d_MET (DisLep.at(0) and DisLep.at(1) are electrons)", 200, 0, 200);
  
  //Isolation (Iso03)
  h.Iso0_ee_l0 = new TH1F("1l2d_ee_Iso03_l0","1l2d_ee_Iso03_l0(d0=e, d1=e)",800, 0, 2);  
  h.Iso0_ee_d0 = new TH1F("1l2d_ee_Iso03_d0","1l2d_ee_Iso03_d0(d0=e, d1=e)",800, 0, 2);
  h.Iso0_ee_d1 = new TH1F("1l2d_ee_Iso03_d1","1l2d_ee_Iso03_d1(d0=e, d1=e)",800, 0, 2);

  //dxy, dz, ip3d, sip3d
  h.dxy0_ee_d0 = new TH1F("1l2d_dxy_ee_d0","1l2d_dxy_d0 (DisLep.at(0) and DisLep.at(1) are electrons)", 200, -1, 1);
  h.dz0_ee_d0 = new TH1F("1l2d_dz_ee_d0","1l2d_dz_d0 (DisLep.at(0) and DisLep.at(1) are electrons)", 200, -1, 1);
  h.ip3d0_ee_d0 = new TH1F("1l2d_ip3d_ee_d0","1l2d_ip3d_d0 (DisLep.at(0) and DisLep.at(1) are electrons)", 200, 0, 20);
  h.sip3d0_ee_d0 = new TH1F("1l2d_sip3d_ee_d0","1l2d_sip3d_d0 (DisLep.at(0) and DisLep.at(1) are electrons)", 200, 0, 200);
  
  //dxy, dz, ip3d, sip3d
  h.dxy0_ee_d1 = new TH1F("1l2d_dxy_ee_d1","1l2d_dxy_d1 (DisLep.at(0) and DisLep.at(1) are electrons)", 200, -1, 1);
  h.dz0_ee_d1 = new TH1F("1l2d_dz_ee_d1","1l2d_dz_d1 (DisLep.at(0) and DisLep.at(1) are electrons)", 200, -1, 1);
  h.ip3d0_ee_d1 = new TH1F("1l2d_ip3d_ee_d1","1l2d_ip3d_d1 (DisLep.at(0) and DisLep.at(1) are electrons)", 200, 0, 20);
  h.sip3d0_ee_d1 = new TH1F("1l2d_sip3d_ee_d1","1l2d_sip3d_d1 (DisLep.at(0) and DisLep.at(1) are electrons)", 200, 0, 200); 
  
  //Transverse mass (mT)
  h.mT0_ee_l0 = new TH1F("1l2d_ee_mT_l0","1l2d_mT_l0 (DisLep.at(0) and DisLep.at(1) are electrons)",200, 0, 200);  
  h.mT0_ee_d0 = new TH1F("1l2d_ee_mT_d0","1l2d_mT_d0 (DisLep.at(0) and DisLep.at(1) are electrons)",200, 0, 200);
  h.mT0_ee_d1 = new TH1F("1l2d_ee_mT_d1","1l2d_mT_d1 (DisLep.at(0) and DisLep.at(1) are electrons)",200, 0, 200);
  
  //njets
  h.njets0_ee = new TH1F("1l2d_ee_njets","1l2d_njets (DisLep.at(0) and DisLep.at(1) are electrons)",10, 0, 10); 

  //flavor classification
  h.flav0_ee = new TH1F("1l2d_ee_flav","1l2d_flav (DisLep.at(0) and DisLep.at(1) are electrons), 0:eee 1:ee#mu 2:e#mu#mu 3:#mu#mu#mu 4:#mu#mue 5:#muee 6:e#mue 7:#mue#mu",8, 0, 8);

  
  
  //#############################################################
  
  
  
  //####################  1l2d channel - l0 mu e  ##########################
  // "_MuE_" => DisLep.at(0) == Mu and DisLep.at(1) == E
  
  
  
  //mass
  h.m0_MuE_l0d0 = new TH1F("1l2d_MuE_m_l0d0","1l2d_m_l0d0 (DisLep.at(0) = mu and DisLep.at(1) = e)", 500, 0, 500);  
  h.m0_MuE_l0d1 = new TH1F("1l2d_MuE_m_l0d1","1l2d_m_l0d1 (DisLep.at(0) = mu and DisLep.at(1) = e )", 500, 0, 500);
  h.m0_MuE_d0d1 = new TH1F("1l2d_MuE_m_d0d1","1l2d_m_d0d1 (DisLep.at(0) = mu and DisLep.at(1)= e)", 500, 0, 500);
  h.m0_MuE_l0d0d1 = new TH1F("1l2d_MuE_m_l0d0d1","1l2d_m_l0d0d1 (DisLep.at(0) = mu and DisLep.at(1) = e)", 500, 0, 500);

  //dPhi
  h.dPhi0_MuE_l0d0 = new TH1F("1l2d_MuE_dPhi_l0d0","1l2d_dPhi_l0d0 (DisLep.at(0) = mu and DisLep.at(1) = e)", 64, -3.2, 3.2);  
  h.dPhi0_MuE_l0d1 = new TH1F("1l2d_MuE_dPhi_l0d1","1l2d_dPhi_l0d1 (DisLep.at(0) = mu and DisLep.at(1) = e)", 64, -3.2, 3.2);
  h.dPhi0_MuE_d0d1 = new TH1F("1l2d_MuE_dPhi_d0d1","1l2d_dPhi_d0d1 (DisLep.at(0) = mu and DisLep.at(1) = e)", 64, -3.2, 3.2);
 
  //dR
  h.dR0_MuE_l0d0 = new TH1F("1l2d_MuE_dR_l0d0","1l2d_dR_l0d0 (DisLep.at(0) = mu and DisLep.at(1) = e)", 100, 0, 10); 
  h.dR0_MuE_l0d1 = new TH1F("1l2d_MuE_dR_l0d1","1l2d_dR_l0d1 (DisLep.at(0) = mu and DisLep.at(1) = e)", 100, 0, 10);
  h.dR0_MuE_d0d1 = new TH1F("1l2d_MuE_dR_d0d1","1l2d_dR_d0d1 (DisLep.at(0) = mu and DisLep.at(1) = e)", 100, 0, 10);

  //pT
  h.pT0_MuE_l0 = new TH1F("1l2d_MuE_pT_l0","1l2d_pT_l0 (DisLep.at(0) = mu and DisLep.at(1) = e)", 200, 0, 200);  
  h.pT0_MuE_d0 = new TH1F("1l2d_MuE_pT_d0","1l2d_pT_d0 (DisLep.at(0) = mu and DisLep.at(1) = e)", 200, 0, 200);
  h.pT0_MuE_d1 = new TH1F("1l2d_MuE_pT_d1","1l2d_pT_d1 (DisLep.at(0) = mu and DisLep.at(1) = e)", 200, 0, 200);
  
  //MET
  h.MET0_MuE = new TH1F("1l2d_MuE_MET","1l2d_MET (DisLep.at(0) = mu and DisLep.at(1) = e)", 200, 0, 200);
  
  //Isolation (Iso03)
  h.Iso0_MuE_l0 = new TH1F("1l2d_MuE_Iso03_l0","1l2d_MuE_Iso03_l0(d0=Mu, d1=E)",800, 0, 2);  
  h.Iso0_MuE_d0 = new TH1F("1l2d_MuE_Iso03_d0","1l2d_MuE_Iso03_d0(d0=Mu, d1=E)",800, 0, 2);
  h.Iso0_MuE_d1 = new TH1F("1l2d_MuE_Iso03_d1","1l2d_MuE_Iso03_d1(d0=Mu, d1=E)",800, 0, 2);

  //dxy, dz, ip3d, sip3d
  h.dxy0_MuE_d0 = new TH1F("1l2d_dxy_MuE_d0","1l2d_dxy_d0 (DisLep.at(0) = mu and DisLep.at(1) = e)", 200, -1, 1);
  h.dz0_MuE_d0 = new TH1F("1l2d_dz_MuE_d0","1l2d_dz_d0 (DisLep.at(0) = mu and DisLep.at(1) = e)", 200, -1, 1);
  h.ip3d0_MuE_d0 = new TH1F("1l2d_ip3d_MuE_d0","1l2d_ip3d_d0 (DisLep.at(0) = mu and DisLep.at(1) = e)", 200, 0, 20);
  h.sip3d0_MuE_d0 = new TH1F("1l2d_sip3d_MuE_d0","1l2d_sip3d_d0 (DisLep.at(0) = mu and DisLep.at(1) = e )", 200, 0, 200);
  
  //dxy, dz, ip3d, sip3d
  h.dxy0_MuE_d1 = new TH1F("1l2d_dxy_MuE_d1","1l2d_dxy_d1 (DisLep.at(0) = mu and DisLep.at(1) = e)", 200, -1, 1);
  h.dz0_MuE_d1 = new TH1F("1l2d_dz_MuE_d1","1l2d_dz_d1 (DisLep.at(0) = mu and DisLep.at(1) = e)", 200, -1, 1);
  h.ip3d0_MuE_d1 = new TH1F("1l2d_ip3d_MuE_d1","1l2d_ip3d_d1 (DisLep.at(0) = mu and DisLep.at(1) = e)", 200, 0, 20);
  h.sip3d0_MuE_d1 = new TH1F("1l2d_sip3d_MuE_d1","1l2d_sip3d_d1 (DisLep.at(0) = muand DisLep.at(1) = e)", 200, 0, 200); 
  
  //Transverse mass (mT)
  h.mT0_MuE_l0 = new TH1F("1l2d_MuE_mT_l0","1l2d_mT_l0 (DisLep.at(0) = mu and DisLep.at(1) = e)",200, 0, 200);  
  h.mT0_MuE_d0 = new TH1F("1l2d_MuE_mT_d0","1l2d_mT_d0 (DisLep.at(0) = mu and DisLep.at(1) = e)",200, 0, 200);
  h.mT0_MuE_d1 = new TH1F("1l2d_MuE_mT_d1","1l2d_mT_d1 (DisLep.at(0) = mu and DisLep.at(1) = e)",200, 0, 200);
  
  //njets
  h.njets0_MuE = new TH1F("1l2d_MuE_njets","1l2d_njets (DisLep.at(0) = mu and DisLep.at(1) = e)",10, 0, 10); 

  //flavor classification
  h.flav0_MuE = new TH1F("1l2d_MuE_flav","1l2d_flav (DisLep.at(0) = mu and DisLep.at(1) = e), 0:eee 1:ee#mu 2:e#mu#mu 3:#mu#mu#mu 4:#mu#mue 5:#muee 6:e#mue 7:#mue#mu",8, 0, 8);
    
  
  
  //#############################################################
  
  
  
  //####################  1l2d channel - l0 mu mu ##########################
  // "_MuMu_" => the two displaced leptons are electrons
  
  
  
  //mass
  h.m0_MuMu_l0d0 = new TH1F("1l2d_MuMu_m_l0d0","1l2d_m_l0d0 (DisLep.at(0) and DisLep.at(1) are muons)", 500, 0, 500);  
  h.m0_MuMu_l0d1 = new TH1F("1l2d_MuMu_m_l0d1","1l2d_m_l0d1 (DisLep.at(0) and DisLep.at(1) are muons)", 500, 0, 500);
  h.m0_MuMu_d0d1 = new TH1F("1l2d_MuMu_m_d0d1","1l2d_m_d0d1 (DisLep.at(0) and DisLep.at(1) are muons)", 500, 0, 500);
  h.m0_MuMu_l0d0d1 = new TH1F("1l2d_MuMu_m_l0d0d1","1l2d_m_l0d0d1 (DisLep.at(0) and DisLep.at(1) are muons)", 500, 0, 500);

  //dPhi
  h.dPhi0_MuMu_l0d0 = new TH1F("1l2d_MuMu_dPhi_l0d0","1l2d_dPhi_l0d0 (DisLep.at(0) and DisLep.at(1) are muons)", 64, -3.2, 3.2);  
  h.dPhi0_MuMu_l0d1 = new TH1F("1l2d_MuMu_dPhi_l0d1","1l2d_dPhi_l0d1 (DisLep.at(0) and DisLep.at(1) are muons)", 64, -3.2, 3.2);
  h.dPhi0_MuMu_d0d1 = new TH1F("1l2d_MuMu_dPhi_d0d1","1l2d_dPhi_d0d1 (DisLep.at(0) and DisLep.at(1) are muons)", 64, -3.2, 3.2);
 
  //dR
  h.dR0_MuMu_l0d0 = new TH1F("1l2d_MuMu_dR_l0d0","1l2d_dR_l0d0 (DisLep.at(0) and DisLep.at(1) are muons)", 100, 0, 10); 
  h.dR0_MuMu_l0d1 = new TH1F("1l2d_MuMu_dR_l0d1","1l2d_dR_l0d1 (DisLep.at(0) and DisLep.at(1) are muons)", 100, 0, 10);
  h.dR0_MuMu_d0d1 = new TH1F("1l2d_MuMu_dR_d0d1","1l2d_dR_d0d1 (DisLep.at(0) and DisLep.at(1) are muons)", 100, 0, 10);

  //pT
  h.pT0_MuMu_l0 = new TH1F("1l2d_MuMu_pT_l0","1l2d_pT_l0 (DisLep.at(0) and DisLep.at(1) are muons)", 200, 0, 200);  
  h.pT0_MuMu_d0 = new TH1F("1l2d_MuMu_pT_d0","1l2d_pT_d0 (DisLep.at(0) and DisLep.at(1) are muons)", 200, 0, 200);
  h.pT0_MuMu_d1 = new TH1F("1l2d_MuMu_pT_d1","1l2d_pT_d1 (DisLep.at(0) and DisLep.at(1) are muons)", 200, 0, 200);
  
  //MET
  h.MET0_MuMu = new TH1F("1l2d_MuMu_MET","1l2d_MET (DisLep.at(0) and DisLep.at(1) are muons)", 200, 0, 200);
  
  //Isolation (Iso03)
  h.Iso0_MuMu_l0 = new TH1F("1l2d_MuMu_Iso03_l0","1l2d_MuMu_Iso03_l0(d0=Mu, d1=Mu)",800, 0, 2);  
  h.Iso0_MuMu_d0 = new TH1F("1l2d_MuMu_Iso03_d0","1l2d_MuMu_Iso03_d0(d0=Mu, d1=Mu)",800, 0, 2);
  h.Iso0_MuMu_d1 = new TH1F("1l2d_MuMu_Iso03_d1","1l2d_MuMu_Iso03_d1(d0=Mu, d1=Mu)",800, 0, 2);

  //dxy, dz, ip3d, sip3d
  h.dxy0_MuMu_d0 = new TH1F("1l2d_dxy_MuMu_d0","1l2d_dxy_d0 (DisLep.at(0) and DisLep.at(1) are muons)", 200, -1, 1);
  h.dz0_MuMu_d0 = new TH1F("1l2d_dz_MuMu_d0","1l2d_dz_d0 (DisLep.at(0) and DisLep.at(1) are muons)", 200, -1, 1);
  h.ip3d0_MuMu_d0 = new TH1F("1l2d_ip3d_MuMu_d0","1l2d_ip3d_d0 (DisLep.at(0) and DisLep.at(1) are muons)", 200, 0, 20);
  h.sip3d0_MuMu_d0 = new TH1F("1l2d_sip3d_MuMu_d0","1l2d_sip3d_d0 (DisLep.at(0) and DisLep.at(1) are muons)", 200, 0, 200);
  
  //dxy, dz, ip3d, sip3d
  h.dxy0_MuMu_d1 = new TH1F("1l2d_dxy_MuMu_d1","1l2d_dxy_d1 (DisLep.at(0) and DisLep.at(1) are muons)", 200, -1, 1);
  h.dz0_MuMu_d1 = new TH1F("1l2d_dz_MuMu_d1","1l2d_dz_d1 (DisLep.at(0) and DisLep.at(1) are muons)", 200, -1, 1);
  h.ip3d0_MuMu_d1 = new TH1F("1l2d_ip3d_MuMu_d1","1l2d_ip3d_d1 (DisLep.at(0) and DisLep.at(1) are muons)", 200, 0, 20);
  h.sip3d0_MuMu_d1 = new TH1F("1l2d_sip3d_MuMu_d1","1l2d_sip3d_d1 (DisLep.at(0) and DisLep.at(1) are muons)", 200, 0, 200); 
  
  //Transverse mass (mT)
  h.mT0_MuMu_l0 = new TH1F("1l2d_MuMu_mT_l0","1l2d_mT_l0 (DisLep.at(0) and DisLep.at(1) are muons)",200, 0, 200);  
  h.mT0_MuMu_d0 = new TH1F("1l2d_MuMu_mT_d0","1l2d_mT_d0 (DisLep.at(0) and DisLep.at(1) are muons)",200, 0, 200);
  h.mT0_MuMu_d1 = new TH1F("1l2d_MuMu_mT_d1","1l2d_mT_d1 (DisLep.at(0) and DisLep.at(1) are muons)",200, 0, 200);
  
  //njets
  h.njets0_MuMu = new TH1F("1l2d_MuMu_njets","1l2d_njets (DisLep.at(0) and DisLep.at(1) are muons)",10, 0, 10); 

  //flavor classification
  h.flav0_MuMu = new TH1F("1l2d_MuMu_flav","1l2d_flav (DisLep.at(0) and DisLep.at(1) are muons), 0:eee 1:ee#mu 2:e#mu#mu 3:#mu#mu#mu 4:#mu#mue 5:#muee 6:e#mue 7:#mue#mu",8, 0, 8);
  
  
  
  //#############################################################
  //#############################################################
  
  
  
  //#######################  0l3d channel  #############################
  
  
  
  //mass
  h.m1_d0d1 = new TH1F("0l3d_m_d0d1","0l3d_m_d0d1", 200, 0, 200);  //mass
  h.m1_d1d2 = new TH1F("0l3d_m_d1d2","0l3d_m_d1d2", 200, 0, 200);
  h.m1_d0d2 = new TH1F("0l3d_m_d0d2","0l3d_m_d0d2", 200, 0, 200);
  h.m1_d0d1d2 = new TH1F("0l3d_m_d0d1d2","0l3d_m_d0d1d2", 400, 0, 400);
  
  //dPhi
  h.dPhi1_d0d1 = new TH1F("0l3d_dPhi_d0d1","0l3d_dPhi_d0d1", 64, -3.2, 3.2);  //dPhi
  h.dPhi1_d1d2 = new TH1F("0l3d_dPhi_d1d2","0l3d_dPhi_d1d2", 64, -3.2, 3.2);
  h.dPhi1_d0d2 = new TH1F("0l3d_dPhi_d0d2","0l3d_dPhi_d0d2", 64, -3.2, 3.2);

  //dPhi_lep_MET
  h.dPhi1_d0MET = new TH1F("0l3d_dPhi_d0MET","0l3d_dPhi_d0_MET", 32, 0, 3.2);  //dPhi_lep_MET
  h.dPhi1_d1MET = new TH1F("0l3d_dPhi_d1MET","0l3d_dPhi_d1_MET", 32, 0, 3.2);
  h.dPhi1_d2MET = new TH1F("0l3d_dPhi_d2MET","0l3d_dPhi_d2_MET", 32, 0, 3.2);

  //dR
  h.dR1_d0d1 = new TH1F("0l3d_dR_d0d1","0l3d_dR_d0d1", 100, 0, 10);  //dR
  h.dR1_d1d2 = new TH1F("0l3d_dR_d1d2","0l3d_dR_d1d2", 100, 0, 10);
  h.dR1_d0d2 = new TH1F("0l3d_dR_d0d2","0l3d_dR_d0d2", 100, 0, 10);

  //pT
  h.pT1_d0 = new TH1F("0l3d_pT_d0","0l3d_pT_d0", 200, 0, 200);  //pT
  h.pT1_d1 = new TH1F("0l3d_pT_d1","0l3d_pT_d1", 200, 0, 200);
  h.pT1_d2 = new TH1F("0l3d_pT_d2","0l3d_pT_d2", 200, 0, 200);
  h.pT1_d0d1 = new TH1F("0l3d_pT_d0d1","0l3d_pT_d0d1", 200, 0, 200);
  h.pT1_d1d2 = new TH1F("0l3d_pT_d1d2","0l3d_pT_d1d2", 200, 0, 200);
  h.pT1_d0d2 = new TH1F("0l3d_pT_d0d2","0l3d_pT_d0d2", 200, 0, 200);
  h.pT1_d0d1d2 = new TH1F("0l3d_pT_d0d1d2","0l3d_pT_d0d1d2", 200, 0, 200);
  
  //Isolation (Iso03)
  h.Iso1_d0 = new TH1F("0l3d_Iso03_d0","0l3d_Iso03_d0",800, 0, 2);  
  h.Iso1_d1 = new TH1F("0l3d_Iso03_d1","0l3d_Iso03_d1",800, 0, 2);
  h.Iso1_d2 = new TH1F("0l3d_Iso03_d2","0l3d_Iso03_d2",800, 0, 2);

  //(dxy, dz, c, ip3d, sip3d)_d0
  h.dxy1_d0 = new TH1F("0l3d_dxy_d0","0l3d_dxy_d0",200, -1, 1);
  h.dz1_d0 = new TH1F("0l3d_dz_d0","0l3d_dz_d0",200, -1, 1);
  h.c1_d0 = new TH1F("0l3d_c_d0","0l3d_sqrt(dxy^2 + dz^2)_d0",200, -1, 1);
  h.ip3d1_d0 = new TH1F("0l3d_ip3d_d0","0l3d_ip3d_d0",200, 0, 20);
  h.sip3d1_d0 = new TH1F("0l3d_sip3d_d0","0l3d_sip3d_d0",200, 0, 200);

  //(dxy, dz, c, ip3d, sip3d)_d1
  h.dxy1_d1 = new TH1F("0l3d_dxy_d1","0l3d_dxy_d1",200, -1, 1);
  h.dz1_d1 = new TH1F("0l3d_dz_d1","0l3d_dz_d1",200, -1, 1);
  h.c1_d1 = new TH1F("0l3d_c_d1","0l3d_sqrt(dxy^2 + dz^2)_d1",200, -1, 1);
  h.ip3d1_d1 = new TH1F("0l3d_ip3d_d1","0l3d_ip3d_d1",200, 0, 20);
  h.sip3d1_d1 = new TH1F("0l3d_sip3d_d1","0l3d_sip3d_d1",200, 0, 200);

  //(dxy, dz, c, ip3d, sip3d)_d2
  h.dxy1_d2 = new TH1F("0l3d_dxy_d2","0l3d_dxy_d2",200, -1, 1);
  h.dz1_d2 = new TH1F("0l3d_dz_d2","0l3d_dz_d2",200, -1, 1);
  h.c1_d2 = new TH1F("0l3d_c_d2","0l3d_sqrt(dxy^2 + dz^2)_d2",200, -1, 1);
  h.ip3d1_d2 = new TH1F("0l3d_ip3d_d2","0l3d_ip3d_d2",200, 0, 20);
  h.sip3d1_d2 = new TH1F("0l3d_sip3d_d2","0l3d_sip3d_d2",200, 0, 200);

  //MET
  h.MET1 = new TH1F("0l3d_MET","0l3d_MET1", 200, 0, 200);  //MET

  //Transverse mass (mT)
  h.mT1_d0 = new TH1F("0l3d_mT_d0","0l3d_mT_d0",200, 0, 200);  //mT
  h.mT1_d1 = new TH1F("0l3d_mT_d1","0l3d_mT_d1",200, 0, 200);
  h.mT1_d2 = new TH1F("0l3d_mT_d2","0l3d_mT_d2",200, 0, 200);

  //dRmin_lep_jet
  h.dRmin1_d0j = new TH1F("0l3d_dRmin_d0j","0l3d_dRmin_d0j",500, 0, 5);  //dRmin_lep_jet
  h.dRmin1_d1j = new TH1F("0l3d_dRmin_d1j","0l3d_dRmin_d1j",500, 0, 5);
  h.dRmin1_d2j = new TH1F("0l3d_dRmin_d2j","0l3d_dRmin_d2j",500, 0, 5);

  //HT
  h.HT1 = new TH1F("0l3d_HT","0l3d_HT",200,0,200);  //HT

  //njets
  h.njets1 = new TH1F("0l3d_njets","0l3d_njets",10, 0, 10);  //njets

  //flavor classification
  h.flav1 = new TH1F("0l3d_flav","0l3d_flav, 0:eee 1:ee#mu 2:e#mu#mu 3:#mu#mu#mu 4:#mu#mue 5:#muee 6:e#mue 7:#mue#mu",8, 0, 8);  //flav



  
  //#############################################################
  //#############################################################

  
  
  //#######################  2l1d channel  #############################


  
  //mass
  h.m2_l0l1 = new TH1F("2l1d_m_l0l1","2l1d_m_l0l1", 500, 0, 500);  
  h.m2_l0d0 = new TH1F("2l1d_m_l0d0","2l1d_m_l0d0", 500, 0, 500);
  h.m2_l1d0 = new TH1F("2l1d_m_l1d0","2l1d_m_l1d0", 500, 0, 500);
  h.m2_l0l1d0 = new TH1F("2l1d_m_l0l1d0","2l1d_m_l0l1d0", 500, 0, 500);

  //dPhi
  h.dPhi2_l0l1 = new TH1F("2l1d_dPhi_l0l1","2l1d_dPhi_l0l1", 64, -3.2, 3.2);  
  h.dPhi2_l0d0 = new TH1F("2l1d_dPhi_l0d0","2l1d_dPhi_l0d0", 64, -3.2, 3.2);
  h.dPhi2_l1d0 = new TH1F("2l1d_dPhi_l1d0","2l1d_dPhi_l1d0", 64, -3.2, 3.2);

  //dPhi_lep_MET
  h.dPhi2_l0MET = new TH1F("2l1d_dPhi_l0MET","2l1d_dPhi_l0_MET", 32, 0, 3.2); 
  h.dPhi2_l1MET = new TH1F("2l1d_dPhi_l1MET","2l1d_dPhi_l1_MET", 32, 0, 3.2);
  h.dPhi2_d0MET = new TH1F("2l1d_dPhi_d0MET","2l1d_dPhi_d0_MET", 32, 0, 3.2);

  //dR
  h.dR2_l0l1 = new TH1F("2l1d_dR_l0l1","2l1d_dR_l0l1", 100, 0, 10); 
  h.dR2_l0d0 = new TH1F("2l1d_dR_l0d0","2l1d_dR_l0d0", 100, 0, 10);
  h.dR2_l1d0 = new TH1F("2l1d_dR_l1d0","2l1d_dR_l1d0", 100, 0, 10);

  //pT
  h.pT2_l0 = new TH1F("2l1d_pT_l0","2l1d_pT_l0", 200, 0, 200);  
  h.pT2_l1 = new TH1F("2l1d_pT_l1","2l1d_pT_l1", 200, 0, 200);
  h.pT2_d0 = new TH1F("2l1d_pT_d0","2l1d_pT_d0", 200, 0, 200);
  h.pT2_l0l1 = new TH1F("2l1d_pT_l0l1","2l1d_pT_l0l1", 200, 0, 200);
  h.pT2_l0d0 = new TH1F("2l1d_pT_l0d0","2l1d_pT_l0d0", 200, 0, 200);
  h.pT2_l1d0 = new TH1F("2l1d_pT_l1d0","2l1d_pT_l1d0", 200, 0, 200);
  h.pT2_l0l1d0 = new TH1F("2l1d_pT_l0l1d0","2l1d_pT_l0l1d0", 200, 0, 200);

  //Isolation (Iso03)
  h.Iso2_l0 = new TH1F("2l1d_Iso03_l0","2l1d_Iso03_l0",800, 0, 2);  
  h.Iso2_l1 = new TH1F("2l1d_Iso03_l1","2l1d_Iso03_l1",800, 0, 2);
  h.Iso2_d0 = new TH1F("2l1d_Iso03_d0","2l1d_Iso03_d0",800, 0, 2);
  
  //(dxy, dz, c, ip3d, sip3d)_l0
  h.dxy2_l0 = new TH1F("2l1d_dxy_l0","2l1d_dxy_l0",200, -1, 1);
  h.dz2_l0 = new TH1F("2l1d_dz_l0","2l1d_dz_l0",200, -1, 1);
  h.c2_l0 = new TH1F("2l1d_c_l0","2l1d_sqrt(dxy^2 + dz^2)_l0",200, -1, 1);
  h.ip3d2_l0 = new TH1F("2l1d_ip3d_l0","2l1d_ip3d_l0",200, 0, 20);
  h.sip3d2_l0 = new TH1F("2l1d_sip3d_l0","2l1d_sip3d_l0",200, 0, 200);
  
  //(dxy, dz, c, ip3d, sip3d)_l1
  h.dxy2_l1 = new TH1F("2l1d_dxy_l1","2l1d_dxy_l1",200, -1, 1);
  h.dz2_l1 = new TH1F("2l1d_dz_l1","2l1d_dz_l1",200, -1, 1);
  h.c2_l1 = new TH1F("2l1d_c_l1","2l1d_sqrt(dxy^2 + dz^2)_l1",200, -1, 1);
  h.ip3d2_l1 = new TH1F("2l1d_ip3d_l1","2l1d_ip3d_l1",200, 0, 20);
  h.sip3d2_l1 = new TH1F("2l1d_sip3d_l1","2l1d_sip3d_l1",200, 0, 200);
  
  //(dxy, dz, c, ip3d, sip3d)_d0
  h.dxy2_d0 = new TH1F("2l1d_dxy_d0","2l1d_dxy_d0",200, -1, 1);
  h.dz2_d0 = new TH1F("2l1d_dz_d0","2l1d_dz_d0",200, -1, 1);
  h.c2_d0 = new TH1F("2l1d_c_d0","2l1d_sqrt(dxy^2 + dz^2)_d0",200, -1, 1);
  h.ip3d2_d0 = new TH1F("2l1d_ip3d_d0","2l1d_ip3d_d0",200, 0, 20);
  h.sip3d2_d0 = new TH1F("2l1d_sip3d_d0","2l1d_sip3d_d0",200, 0, 200);

  //MET
  h.MET2 = new TH1F("2l1d_MET","2l1d_MET2", 200, 0, 200); 

  //Transverse mass (mT)
  h.mT2_l0 = new TH1F("2l1d_mT_l0","2l1d_mT_l0",200, 0, 200);  
  h.mT2_l1 = new TH1F("2l1d_mT_l1","2l1d_mT_l1",200, 0, 200);
  h.mT2_d0 = new TH1F("2l1d_mT_d0","2l1d_mT_d0",200, 0, 200);

  //dRmin_lep_jet
  h.dRmin2_l0j = new TH1F("2l1d_dRmin_l0j","2l1d_dRmin_l0j",500, 0, 5); 
  h.dRmin2_l1j = new TH1F("2l1d_dRmin_l1j","2l1d_dRmin_l1j",500, 0, 5);
  h.dRmin2_d0j = new TH1F("2l1d_dRmin_d0j","2l1d_dRmin_d0j",500, 0, 5);

  //HT
  h.HT2 = new TH1F("2l1d_HT","2l1d_HT",200,0,200); 

  //njets
  h.njets2 = new TH1F("2l1d_njets","2l1d_njets",10, 0, 10); 

  //flavor classification
  h.flav2 = new TH1F("2l1d_flav","2l1d_flav, 0:eee 1:ee#mu 2:e#mu#mu 3:#mu#mu#mu 4:#mu#mue 5:#muee 6:e#mue 7:#mue#mu",8, 0, 8);



  //#############################################################

  
  
  //#######################  2l1d channel - Z vetoed #############################


 
  //mass
  h.m2_Zveto_l0l1 = new TH1F("2l1d_Zveto_m_l0l1","2l1d_Zveto_m_l0l1", 500, 0, 500);  
  h.m2_Zveto_l0d0 = new TH1F("2l1d_Zveto_m_l0d0","2l1d_Zveto_m_l0d0", 500, 0, 500);
  h.m2_Zveto_l1d0 = new TH1F("2l1d_Zveto_m_l1d0","2l1d_Zveto_m_l1d0", 500, 0, 500);
  h.m2_Zveto_l0l1d0 = new TH1F("2l1d_Zveto_m_l0l1d0","2l1d_Zveto_m_l0l1d0", 500, 0, 500);

  //dPhi
  h.dPhi2_Zveto_l0l1 = new TH1F("2l1d_Zveto_dPhi_l0l1","2l1d_Zveto_dPhi_l0l1", 64, -3.2, 3.2);  
  h.dPhi2_Zveto_l0d0 = new TH1F("2l1d_Zveto_dPhi_l0d0","2l1d_Zveto_dPhi_l0d0", 64, -3.2, 3.2);
  h.dPhi2_Zveto_l1d0 = new TH1F("2l1d_Zveto_dPhi_l1d0","2l1d_Zveto_dPhi_l1d0", 64, -3.2, 3.2);

  //dPhi_lep_MET
  h.dPhi2_Zveto_l0MET = new TH1F("2l1d_Zveto_dPhi_l0MET","2l1d_Zveto_dPhi_l0_MET", 32, 0, 3.2); 
  h.dPhi2_Zveto_l1MET = new TH1F("2l1d_Zveto_dPhi_l1MET","2l1d_Zveto_dPhi_l1_MET", 32, 0, 3.2);
  h.dPhi2_Zveto_d0MET = new TH1F("2l1d_Zveto_dPhi_d0MET","2l1d_Zveto_dPhi_d0_MET", 32, 0, 3.2);

  //dR
  h.dR2_Zveto_l0l1 = new TH1F("2l1d_Zveto_dR_l0l1","2l1d_Zveto_dR_l0l1", 100, 0, 10); 
  h.dR2_Zveto_l0d0 = new TH1F("2l1d_Zveto_dR_l0d0","2l1d_Zveto_dR_l0d0", 100, 0, 10);
  h.dR2_Zveto_l1d0 = new TH1F("2l1d_Zveto_dR_l1d0","2l1d_Zveto_dR_l1d0", 100, 0, 10);

  //pT
  h.pT2_Zveto_l0 = new TH1F("2l1d_Zveto_pT_l0","2l1d_Zveto_pT_l0", 200, 0, 200);  
  h.pT2_Zveto_l1 = new TH1F("2l1d_Zveto_pT_l1","2l1d_Zveto_pT_l1", 200, 0, 200);
  h.pT2_Zveto_d0 = new TH1F("2l1d_Zveto_pT_d0","2l1d_Zveto_pT_d0", 200, 0, 200);
  h.pT2_Zveto_l0l1 = new TH1F("2l1d_Zveto_pT_l0l1","2l1d_Zveto_pT_l0l1", 200, 0, 200);
  h.pT2_Zveto_l0d0 = new TH1F("2l1d_Zveto_pT_l0d0","2l1d_Zveto_pT_l0d0", 200, 0, 200);
  h.pT2_Zveto_l1d0 = new TH1F("2l1d_Zveto_pT_l1d0","2l1d_Zveto_pT_l1d0", 200, 0, 200);
  h.pT2_Zveto_l0l1d0 = new TH1F("2l1d_Zveto_pT_l0l1d0","2l1d_Zveto_pT_l0l1d0", 200, 0, 200);

  //Isolation (Iso03)
  h.Iso2_Zveto_l0 = new TH1F("2l1d_Zveto_Iso03_l0","2l1d_Zveto_Iso03_l0",800, 0, 2);  
  h.Iso2_Zveto_l1 = new TH1F("2l1d_Zveto_Iso03_l1","2l1d_Zveto_Iso03_l1",800, 0, 2);
  h.Iso2_Zveto_d0 = new TH1F("2l1d_Zveto_Iso03_d0","2l1d_Zveto_Iso03_d0",800, 0, 2);
  
  //(dxy, dz, c, ip3d, sip3d)_l0
  h.dxy2_Zveto_l0 = new TH1F("2l1d_Zveto_dxy_l0","2l1d_Zveto_dxy_l0",200, -1, 1);
  h.dz2_Zveto_l0 = new TH1F("2l1d_Zveto_dz_l0","2l1d_Zveto_dz_l0",200, -1, 1);
  h.c2_Zveto_l0 = new TH1F("2l1d_Zveto_c_l0","2l1d_Zveto_sqrt(dxy^2 + dz^2)_l0",200, -1, 1);
  h.ip3d2_Zveto_l0 = new TH1F("2l1d_Zveto_ip3d_l0","2l1d_Zveto_ip3d_l0",200, 0, 20);
  h.sip3d2_Zveto_l0 = new TH1F("2l1d_Zveto_sip3d_l0","2l1d_Zveto_sip3d_l0",200, 0, 200);
  
  //(dxy, dz, c, ip3d, sip3d)_l1
  h.dxy2_Zveto_l1 = new TH1F("2l1d_Zveto_dxy_l1","2l1d_Zveto_dxy_l1",200, -1, 1);
  h.dz2_Zveto_l1 = new TH1F("2l1d_Zveto_dz_l1","2l1d_Zveto_dz_l1",200, -1, 1);
  h.c2_Zveto_l1 = new TH1F("2l1d_Zveto_c_l1","2l1d_Zveto_sqrt(dxy^2 + dz^2)_l1",200, -1, 1);
  h.ip3d2_Zveto_l1 = new TH1F("2l1d_Zveto_ip3d_l1","2l1d_Zveto_ip3d_l1",200, 0, 20);
  h.sip3d2_Zveto_l1 = new TH1F("2l1d_Zveto_sip3d_l1","2l1d_Zveto_sip3d_l1",200, 0, 200);
  
  //(dxy, dz, c, ip3d, sip3d)_d0
  h.dxy2_Zveto_d0 = new TH1F("2l1d_Zveto_dxy_d0","2l1d_Zveto_dxy_d0",200, -1, 1);
  h.dz2_Zveto_d0 = new TH1F("2l1d_Zveto_dz_d0","2l1d_Zveto_dz_d0",200, -1, 1);
  h.c2_Zveto_d0 = new TH1F("2l1d_Zveto_c_d0","2l1d_Zveto_sqrt(dxy^2 + dz^2)_d0",200, -1, 1);
  h.ip3d2_Zveto_d0 = new TH1F("2l1d_Zveto_ip3d_d0","2l1d_Zveto_ip3d_d0",200, 0, 20);
  h.sip3d2_Zveto_d0 = new TH1F("2l1d_Zveto_sip3d_d0","2l1d_Zveto_sip3d_d0",200, 0, 200);

  //MET
  h.MET2_Zveto = new TH1F("2l1d_Zveto_MET","2l1d_Zveto_MET2", 200, 0, 200); 

  //Transverse mass (mT)
  h.mT2_Zveto_l0 = new TH1F("2l1d_Zveto_mT_l0","2l1d_Zveto_mT_l0",200, 0, 200);  
  h.mT2_Zveto_l1 = new TH1F("2l1d_Zveto_mT_l1","2l1d_Zveto_mT_l1",200, 0, 200);
  h.mT2_Zveto_d0 = new TH1F("2l1d_Zveto_mT_d0","2l1d_Zveto_mT_d0",200, 0, 200);

  //dRmin_lep_jet
  h.dRmin2_Zveto_l0j = new TH1F("2l1d_Zveto_dRmin_l0j","2l1d_dRmin_l0j",500, 0, 5); 
  h.dRmin2_Zveto_l1j = new TH1F("2l1d_Zveto_dRmin_l1j","2l1d_Zveto_dRmin_l1j",500, 0, 5);
  h.dRmin2_Zveto_d0j = new TH1F("2l1d_Zveto_dRmin_d0j","2l1d_Zveto_dRmin_d0j",500, 0, 5);

  //HT
  h.HT2_Zveto = new TH1F("2l1d_Zveto_HT","2l1d_Zveto_HT",200,0,200); 

  //njets
  h.njets2_Zveto = new TH1F("2l1d_Zveto_njets","2l1d_Zveto_njets",10, 0, 10); 

  //flavor classification
  h.flav2_Zveto = new TH1F("2l1d_Zveto_flav","2l1d_Zveto_flav, 0:eee 1:ee#mu 2:e#mu#mu 3:#mu#mu#mu 4:#mu#mue 5:#muee 6:e#mue 7:#mue#mu",8, 0, 8);



  //#############################################################

  
  
  //#######################  2l1d channel - CR #############################


 
  //mass
  h.m2_Z_CR_l0l1 = new TH1F("2l1d_Z_CR_m_l0l1","2l1d_Z_CR_m_l0l1", 200, 0, 200);  
  h.m2_Z_CR_l0d0 = new TH1F("2l1d_Z_CR_m_l0d0","2l1d_Z_CR_m_l0d0", 200, 0, 200);
  h.m2_Z_CR_l1d0 = new TH1F("2l1d_Z_CR_m_l1d0","2l1d_Z_CR_m_l1d0", 200, 0, 200);
  h.m2_Z_CR_l0l1d0 = new TH1F("2l1d_Z_CR_m_l0l1d0","2l1d_Z_CR_m_l0l1d0", 200, 0, 200);

  //dPhi
  h.dPhi2_Z_CR_l0l1 = new TH1F("2l1d_Z_CR_dPhi_l0l1","2l1d_Z_CR_dPhi_l0l1", 64, -3.2, 3.2);  
  h.dPhi2_Z_CR_l0d0 = new TH1F("2l1d_Z_CR_dPhi_l0d0","2l1d_Z_CR_dPhi_l0d0", 64, -3.2, 3.2);
  h.dPhi2_Z_CR_l1d0 = new TH1F("2l1d_Z_CR_dPhi_l1d0","2l1d_Z_CR_dPhi_l1d0", 64, -3.2, 3.2);

  //dPhi_lep_MET
  h.dPhi2_Z_CR_l0MET = new TH1F("2l1d_Z_CR_dPhi_l0MET","2l1d_Z_CR_dPhi_l0_MET", 32, 0, 3.2); 
  h.dPhi2_Z_CR_l1MET = new TH1F("2l1d_Z_CR_dPhi_l1MET","2l1d_Z_CR_dPhi_l1_MET", 32, 0, 3.2);
  h.dPhi2_Z_CR_d0MET = new TH1F("2l1d_Z_CR_dPhi_d0MET","2l1d_Z_CR_dPhi_d0_MET", 32, 0, 3.2);

  //dR
  h.dR2_Z_CR_l0l1 = new TH1F("2l1d_Z_CR_dR_l0l1","2l1d_Z_CR_dR_l0l1", 100, 0, 10); 
  h.dR2_Z_CR_l0d0 = new TH1F("2l1d_Z_CR_dR_l0d0","2l1d_Z_CR_dR_l0d0", 100, 0, 10);
  h.dR2_Z_CR_l1d0 = new TH1F("2l1d_Z_CR_dR_l1d0","2l1d_Z_CR_dR_l1d0", 100, 0, 10);

  //pT
  h.pT2_Z_CR_l0 = new TH1F("2l1d_Z_CR_pT_l0","2l1d_Z_CR_pT_l0", 200, 0, 200);  
  h.pT2_Z_CR_l1 = new TH1F("2l1d_Z_CR_pT_l1","2l1d_Z_CR_pT_l1", 200, 0, 200);
  h.pT2_Z_CR_d0 = new TH1F("2l1d_Z_CR_pT_d0","2l1d_Z_CR_pT_d0", 200, 0, 200);
  h.pT2_Z_CR_l0l1 = new TH1F("2l1d_Z_CR_pT_l0l1","2l1d_Z_CR_pT_l0l1", 200, 0, 200);
  h.pT2_Z_CR_l0d0 = new TH1F("2l1d_Z_CR_pT_l0d0","2l1d_Z_CR_pT_l0d0", 200, 0, 200);
  h.pT2_Z_CR_l1d0 = new TH1F("2l1d_Z_CR_pT_l1d0","2l1d_Z_CR_pT_l1d0", 200, 0, 200);
  h.pT2_Z_CR_l0l1d0 = new TH1F("2l1d_Z_CR_pT_l0l1d0","2l1d_Z_CR_pT_l0l1d0", 200, 0, 200);

  //Isolation (Iso03)
  h.Iso2_Z_CR_l0 = new TH1F("2l1d_Z_CR_Iso03_l0","2l1d_Z_CR_Iso03_l0",800, 0, 2);  
  h.Iso2_Z_CR_l1 = new TH1F("2l1d_Z_CR_Iso03_l1","2l1d_Z_CR_Iso03_l1",800, 0, 2);
  h.Iso2_Z_CR_d0 = new TH1F("2l1d_Z_CR_Iso03_d0","2l1d_Z_CR_Iso03_d0",800, 0, 2);
  
  //(dxy, dz, c, ip3d, sip3d)_l0
  h.dxy2_Z_CR_l0 = new TH1F("2l1d_Z_CR_dxy_l0","2l1d_Z_CR_dxy_l0",200, -1, 1);
  h.dz2_Z_CR_l0 = new TH1F("2l1d_Z_CR_dz_l0","2l1d_Z_CR_dz_l0",200, -1, 1);
  h.c2_Z_CR_l0 = new TH1F("2l1d_Z_CR_c_l0","2l1d_Z_CR_sqrt(dxy^2 + dz^2)_l0",200, -1, 1);
  h.ip3d2_Z_CR_l0 = new TH1F("2l1d_Z_CR_ip3d_l0","2l1d_Z_CR_ip3d_l0",200, 0, 20);
  h.sip3d2_Z_CR_l0 = new TH1F("2l1d_Z_CR_sip3d_l0","2l1d_Z_CR_sip3d_l0",200, 0, 200);
  
  //(dxy, dz, c, ip3d, sip3d)_l1
  h.dxy2_Z_CR_l1 = new TH1F("2l1d_Z_CR_dxy_l1","2l1d_Z_CR_dxy_l1",200, -1, 1);
  h.dz2_Z_CR_l1 = new TH1F("2l1d_Z_CR_dz_l1","2l1d_Z_CR_dz_l1",200, -1, 1);
  h.c2_Z_CR_l1 = new TH1F("2l1d_Z_CR_c_l1","2l1d_Z_CR_sqrt(dxy^2 + dz^2)_l1",200, -1, 1);
  h.ip3d2_Z_CR_l1 = new TH1F("2l1d_Z_CR_ip3d_l1","2l1d_Z_CR_ip3d_l1",200, 0, 20);
  h.sip3d2_Z_CR_l1 = new TH1F("2l1d_Z_CR_sip3d_l1","2l1d_Z_CR_sip3d_l1",200, 0, 200);
  
  //(dxy, dz, c, ip3d, sip3d)_d0
  h.dxy2_Z_CR_d0 = new TH1F("2l1d_Z_CR_dxy_d0","2l1d_Z_CR_dxy_d0",200, -1, 1);
  h.dz2_Z_CR_d0 = new TH1F("2l1d_Z_CR_dz_d0","2l1d_Z_CR_dz_d0",200, -1, 1);
  h.c2_Z_CR_d0 = new TH1F("2l1d_Z_CR_c_d0","2l1d_Z_CR_sqrt(dxy^2 + dz^2)_d0",200, -1, 1);
  h.ip3d2_Z_CR_d0 = new TH1F("2l1d_Z_CR_ip3d_d0","2l1d_Z_CR_ip3d_d0",200, 0, 20);
  h.sip3d2_Z_CR_d0 = new TH1F("2l1d_Z_CR_sip3d_d0","2l1d_Z_CR_sip3d_d0",200, 0, 200);

  //MET
  h.MET2_Z_CR = new TH1F("2l1d_Z_CR_MET","2l1d_Z_CR_MET2", 200, 0, 200); 

  //Transverse mass (mT)
  h.mT2_Z_CR_l0 = new TH1F("2l1d_Z_CR_mT_l0","2l1d_Z_CR_mT_l0",200, 0, 200);  
  h.mT2_Z_CR_l1 = new TH1F("2l1d_Z_CR_mT_l1","2l1d_Z_CR_mT_l1",200, 0, 200);
  h.mT2_Z_CR_d0 = new TH1F("2l1d_Z_CR_mT_d0","2l1d_Z_CR_mT_d0",200, 0, 200);

  //dRmin_lep_jet
  h.dRmin2_Z_CR_l0j = new TH1F("2l1d_Z_CR_dRmin_l0j","2l1d_Z_CR_dRmin_l0j",500, 0, 5); 
  h.dRmin2_Z_CR_l1j = new TH1F("2l1d_Z_CR_dRmin_l1j","2l1d_Z_CR_dRmin_l1j",500, 0, 5);
  h.dRmin2_Z_CR_d0j = new TH1F("2l1d_Z_CR_dRmin_d0j","2l1d_Z_CR_dRmin_d0j",500, 0, 5);

  //HT
  h.HT2_Z_CR = new TH1F("2l1d_Z_CR_HT","2l1d_Z_CR_HT",200,0,200); 

  //njets
  h.njets2_Z_CR = new TH1F("2l1d_Z_CR_njets","2l1d_Z_CR_njets",10, 0, 10); 

  //flavor classification
  h.flav2_Z_CR = new TH1F("2l1d_Z_CR_flav","2l1d_Z_CR_flav, 0:eee 1:ee#mu 2:e#mu#mu 3:#mu#mu#mu 4:#mu#mue 5:#muee 6:e#mue 7:#mue#mu",8, 0, 8);
  
  
  
  //#############################################################
  
  
  
  //###################  2l1d channel - m_l0l1<12  #########################
  
  
  
  //mass
  h.m2_12_l0l1 = new TH1F("2l1d_12_m_l0l1","2l1d_m_l0l1(m_l0l1<12)", 500, 0, 500);  
  h.m2_12_l0d0 = new TH1F("2l1d_12_m_l0d0","2l1d_m_l0d0(m_l0l1<12)", 500, 0, 500);
  h.m2_12_l1d0 = new TH1F("2l1d_12_m_l1d0","2l1d_m_l1d0(m_l0l1<12)", 500, 0, 500);
  h.m2_12_l0l1d0 = new TH1F("2l1d_12_m_l0l1d0","2l1d_m_l0l1d0(m_l0l1<12)", 500, 0, 500);
  
  //dPhi
  h.dPhi2_12_l0l1 = new TH1F("2l1d_12_dPhi_l0l1","2l1d_dPhi_l0l1(m_l0l1<12)", 64, -3.2, 3.2);  
  h.dPhi2_12_l0d0 = new TH1F("2l1d_12_dPhi_l0d0","2l1d_dPhi_l0d0(m_l0l1<12)", 64, -3.2, 3.2);
  h.dPhi2_12_l1d0 = new TH1F("2l1d_12_dPhi_l1d0","2l1d_dPhi_l1d0(m_l0l1<12)", 64, -3.2, 3.2);

  //dR
  h.dR2_12_l0l1 = new TH1F("2l1d_12_dR_l0l1","2l1d_dR_l0l1(m_l0l1<12)", 100, 0, 10); 
  h.dR2_12_l0d0 = new TH1F("2l1d_12_dR_l0d0","2l1d_dR_l0d0(m_l0l1<12)", 100, 0, 10);
  h.dR2_12_l1d0 = new TH1F("2l1d_12_dR_l1d0","2l1d_dR_l1d0(m_l0l1<12)", 100, 0, 10);

  //pT
  h.pT2_12_l0 = new TH1F("2l1d_12_pT_l0","2l1d_pT_l0(m_l0l1<12)", 200, 0, 200);  
  h.pT2_12_l1 = new TH1F("2l1d_12_pT_l1","2l1d_pT_l1(m_l0l1<12)", 200, 0, 200);
  h.pT2_12_d0 = new TH1F("2l1d_12_pT_d0","2l1d_pT_d0(m_l0l1<12)", 200, 0, 200);
  
  //MET
  h.MET2_12 = new TH1F("2l1d_12_MET","2l1d_MET2(m_l0l1<12)", 200, 0, 200);
  
  //Isolation (Iso03)
  h.Iso2_12_l0 = new TH1F("2l1d_12_Iso03_l0","2l1d_12_Iso03_l0(m_l0l1<12)",800, 0, 2);  
  h.Iso2_12_l1 = new TH1F("2l1d_12_Iso03_l1","2l1d_12_Iso03_l1(m_l0l1<12)",800, 0, 2);
  h.Iso2_12_d0 = new TH1F("2l1d_12_Iso03_d0","2l1d_12_Iso03_d0(m_l0l1<12)",800, 0, 2);

  //dxy, dz, ip3d, sip3d_d0
  h.dxy2_12_l0 = new TH1F("2l1d_12_dxy_l0","2l1d_12_dxy_l0(m_l0l1<12)",200, -1, 1);
  h.dz2_12_l0 = new TH1F("2l1d_12_dz_l0","2l1d_12_dz_l0(m_l0l1<12)",200, -1, 1);
  h.ip3d2_12_l0 = new TH1F("2l1d_12_ip3d_l0","2l1d_12_ip3d_l0(m_l0l1<12)",200, 0, 20);
  h.sip3d2_12_l0 = new TH1F("2l1d_12_sip3d_l0","2l1d_12_sip3d_l0(m_l0l1<12)",200, 0, 200);

  //dxy, dz, ip3d, sip3d_l1
  h.dxy2_12_l1 = new TH1F("2l1d_12_dxy_l1","2l1d_12_dxy_l1(m_l0l1<12)",200, -1, 1);
  h.dz2_12_l1 = new TH1F("2l1d_12_dz_l1","2l1d_12_dz_l1(m_l0l1<12)",200, -1, 1);
  h.ip3d2_12_l1 = new TH1F("2l1d_12_ip3d_l1","2l1d_12_ip3d_l1(m_l0l1<12)",200, 0, 20);
  h.sip3d2_12_l1 = new TH1F("2l1d_12_sip3d_l1","2l1d_12_sip3d_l1(m_l0l1<12)",200, 0, 200);
  
  //dxy, dz, ip3d, sip3d_d0
  h.dxy2_12_d0 = new TH1F("2l1d_12_dxy_d0","2l1d_12_dxy_d0(m_l0l1<12)",200, -1, 1);
  h.dz2_12_d0 = new TH1F("2l1d_12_dz_d0","2l1d_12_dz_d0(m_l0l1<12)",200, -1, 1);
  h.ip3d2_12_d0 = new TH1F("2l1d_12_ip3d_d0","2l1d_12_ip3d_d0(m_l0l1<12)",200, 0, 20);
  h.sip3d2_12_d0 = new TH1F("2l1d_12_sip3d_d0","2l1d_12_sip3d_d0(m_l0l1<12)",200, 0, 200);
  
  //Transverse mass (mT)
  h.mT2_12_l0 = new TH1F("2l1d_12_mT_l0","2l1d_mT_l0(m_l0l1<12)",200, 0, 200);  
  h.mT2_12_l1 = new TH1F("2l1d_12_mT_l1","2l1d_mT_l1(m_l0l1<12)",200, 0, 200);
  h.mT2_12_d0 = new TH1F("2l1d_12_mT_d0","2l1d_mT_d0(m_l0l1<12)",200, 0, 200);

  //njets
  h.njets2_12 = new TH1F("2l1d_12_njets","2l1d_njets(m_l0l1<12)",10, 0, 10); 

  //flavor classification
  h.flav2_12 = new TH1F("2l1d_12_flav","2l1d_flav(m_l0l1<12), 0:eee 1:ee#mu 2:e#mu#mu 3:#mu#mu#mu 4:#mu#mue 5:#muee 6:e#mue 7:#mue#mu",8, 0, 8);
  
  
  
  //#############################################################
  
  
  
  //####################  2l1d channel - l0 l1 e  ###########################
  
  
  
  //mass
  h.m2_e_l0l1 = new TH1F("2l1d_e_m_l0l1","2l1d_m_l0l1 (DisLep.at(0) is electron)", 500, 0, 500);  
  h.m2_e_l0d0 = new TH1F("2l1d_e_m_l0d0","2l1d_m_l0d0 (DisLep.at(0) is electron)", 500, 0, 500);
  h.m2_e_l1d0 = new TH1F("2l1d_e_m_l1d0","2l1d_m_l1d0 (DisLep.at(0) is electron)", 500, 0, 500);
  h.m2_e_l0l1d0 = new TH1F("2l1d_e_m_l0l1d0","2l1d_m_l0l1d0 (DisLep.at(0) is electron)", 500, 0, 500);
  
  //dPhi
  h.dPhi2_e_l0l1 = new TH1F("2l1d_e_dPhi_l0l1","2l1d_dPhi_l0l1 (DisLep.at(0) is electron)", 64, -3.2, 3.2);  
  h.dPhi2_e_l0d0 = new TH1F("2l1d_e_dPhi_l0d0","2l1d_dPhi_l0d0 (DisLep.at(0) is electron)", 64, -3.2, 3.2);
  h.dPhi2_e_l1d0 = new TH1F("2l1d_e_dPhi_l1d0","2l1d_dPhi_l1d0 (DisLep.at(0) is electron)", 64, -3.2, 3.2);

  //dR
  h.dR2_e_l0l1 = new TH1F("2l1d_e_dR_l0l1","2l1d_dR_l0l1 (DisLep.at(0) is electron)", 100, 0, 10); 
  h.dR2_e_l0d0 = new TH1F("2l1d_e_dR_l0d0","2l1d_dR_l0d0 (DisLep.at(0) is electron)", 100, 0, 10);
  h.dR2_e_l1d0 = new TH1F("2l1d_e_dR_l1d0","2l1d_dR_l1d0 (DisLep.at(0) is electron)", 100, 0, 10);

  //pT
  h.pT2_e_l0 = new TH1F("2l1d_e_pT_l0","2l1d_pT_l0 (DisLep.at(0) is electron)", 200, 0, 200);  
  h.pT2_e_l1 = new TH1F("2l1d_e_pT_l1","2l1d_pT_l1 (DisLep.at(0) is electron)", 200, 0, 200);
  h.pT2_e_d0 = new TH1F("2l1d_e_pT_d0","2l1d_pT_d0 (DisLep.at(0) is electron)", 200, 0, 200);
  
  //MET
  h.MET2_e = new TH1F("2l1d_e_MET","2l1d_MET2 (DisLep.at(0) is electron)", 200, 0, 200);
  
  //Isolation (Iso03)
  h.Iso2_e_l0 = new TH1F("2l1d_e_Iso03_l0","2l1d_e_Iso03_l0(d0=e)",800, 0, 2);  
  h.Iso2_e_l1 = new TH1F("2l1d_e_Iso03_l1","2l1d_e_Iso03_l1(d0=e)",800, 0, 2);
  h.Iso2_e_d0 = new TH1F("2l1d_e_Iso03_d0","2l1d_e_Iso03_d0(d0=e)",800, 0, 2);

  //dxy, dz, ip3d, sip3d
  h.dxy2_e_d0 = new TH1F("2l1d_dxy_e","2l1d_dxy_e (DisLep.at(0) is electron)", 200, -1, 1);
  h.dz2_e_d0 = new TH1F("2l1d_dz_e","2l1d_dz_e (DisLep.at(0) is electron)", 200, -1, 1);
  h.ip3d2_e_d0 = new TH1F("2l1d_ip3d_e","2l1d_ip3d_e (DisLep.at(0) is electron)", 200, 0, 20);
  h.sip3d2_e_d0 = new TH1F("2l1d_sip3d_e","2l1d_sip3d_e (DisLep.at(0) is electron)", 200, 0, 200); 

  //Transverse mass (mT)
  h.mT2_e_l0 = new TH1F("2l1d_e_mT_l0","2l1d_mT_l0 (DisLep.at(0) is electron)",200, 0, 200);  
  h.mT2_e_l1 = new TH1F("2l1d_e_mT_l1","2l1d_mT_l1 (DisLep.at(0) is electron)",200, 0, 200);
  h.mT2_e_d0 = new TH1F("2l1d_e_mT_d0","2l1d_mT_d0 (DisLep.at(0) is electron)",200, 0, 200);

  //njets
  h.njets2_e = new TH1F("2l1d_e_njets","2l1d_njets (DisLep.at(0) is electron)",10, 0, 10); 

  //flavor classification
  h.flav2_e = new TH1F("2l1d_e_flav","2l1d_flav (DisLep.at(0) is electron), 0:eee 1:ee#mu 2:e#mu#mu 3:#mu#mu#mu 4:#mu#mue 5:#muee 6:e#mue 7:#mue#mu",8, 0, 8);
  
  
  
  //#############################################################
  
  
  
  //####################  2l1d channel - l0 l1 mu  ##########################
  
  
  
  //mass
  h.m2_mu_l0l1 = new TH1F("2l1d_mu_m_l0l1","2l1d_m_l0l1 (DisLep.at(0) is muon)", 500, 0, 500);  
  h.m2_mu_l0d0 = new TH1F("2l1d_mu_m_l0d0","2l1d_m_l0d0 (DisLep.at(0) is muon)", 500, 0, 500);
  h.m2_mu_l1d0 = new TH1F("2l1d_mu_m_l1d0","2l1d_m_l1d0 (DisLep.at(0) is muon)", 500, 0, 500);
  h.m2_mu_l0l1d0 = new TH1F("2l1d_mu_m_l0l1d0","2l1d_m_l0l1d0 (DisLep.at(0) is muon)", 500, 0, 500);
  
  //dPhi
  h.dPhi2_mu_l0l1 = new TH1F("2l1d_mu_dPhi_l0l1","2l1d_dPhi_l0l1 (DisLep.at(0) is muon)", 64, -3.2, 3.2);  
  h.dPhi2_mu_l0d0 = new TH1F("2l1d_mu_dPhi_l0d0","2l1d_dPhi_l0d0 (DisLep.at(0) is muon)", 64, -3.2, 3.2);
  h.dPhi2_mu_l1d0 = new TH1F("2l1d_mu_dPhi_l1d0","2l1d_dPhi_l1d0 (DisLep.at(0) is muon)", 64, -3.2, 3.2);

  //dR
  h.dR2_mu_l0l1 = new TH1F("2l1d_mu_dR_l0l1","2l1d_dR_l0l1 (DisLep.at(0) is muon)", 100, 0, 10); 
  h.dR2_mu_l0d0 = new TH1F("2l1d_mu_dR_l0d0","2l1d_dR_l0d0 (DisLep.at(0) is muon)", 100, 0, 10);
  h.dR2_mu_l1d0 = new TH1F("2l1d_mu_dR_l1d0","2l1d_dR_l1d0 (DisLep.at(0) is muon)", 100, 0, 10);

  //pT
  h.pT2_mu_l0 = new TH1F("2l1d_mu_pT_l0","2l1d_pT_l0 (DisLep.at(0) is muon)", 200, 0, 200);  
  h.pT2_mu_l1 = new TH1F("2l1d_mu_pT_l1","2l1d_pT_l1 (DisLep.at(0) is muon)", 200, 0, 200);
  h.pT2_mu_d0 = new TH1F("2l1d_mu_pT_d0","2l1d_pT_d0 (DisLep.at(0) is muon)", 200, 0, 200);
  
  //MET
  h.MET2_mu = new TH1F("2l1d_mu_MET","2l1d_MET2 (DisLep.at(0) is muon)", 200, 0, 200);
  
  //Isolation (Iso03)
  h.Iso2_mu_l0 = new TH1F("2l1d_mu_Iso03_l0","2l1d_mu_Iso03_l0(d0=mu)",800, 0, 2);  
  h.Iso2_mu_l1 = new TH1F("2l1d_mu_Iso03_l1","2l1d_mu_Iso03_l1(d0=mu)",800, 0, 2);
  h.Iso2_mu_d0 = new TH1F("2l1d_mu_Iso03_d0","2l1d_mu_Iso03_d0(d0=mu)",800, 0, 2);
  
  //dxy, dz, ip3d, sip3d
  h.dxy2_mu_d0 = new TH1F("2l1d_dxy_mu","2l1d_dxy_e (DisLep.at(0) is muon)", 200, -1, 1);
  h.dz2_mu_d0 = new TH1F("2l1d_dz_mu","2l1d_dz_e (DisLep.at(0) is muon)", 200, -1, 1);
  h.ip3d2_mu_d0 = new TH1F("2l1d_ip3d_mu","2l1d_ip3d_e (DisLep.at(0) is muon)", 200, 0, 20);
  h.sip3d2_mu_d0 = new TH1F("2l1d_sip3d_mu","2l1d_sip3d_e (DisLep.at(0) is muon)", 200, 0, 200); 

  //Transverse mass (mT)
  h.mT2_mu_l0 = new TH1F("2l1d_mu_mT_l0","2l1d_mT_l0 (DisLep.at(0) is muon)",200, 0, 200);  
  h.mT2_mu_l1 = new TH1F("2l1d_mu_mT_l1","2l1d_mT_l1 (DisLep.at(0) is muon)",200, 0, 200);
  h.mT2_mu_d0 = new TH1F("2l1d_mu_mT_d0","2l1d_mT_d0 (DisLep.at(0) is muon)",200, 0, 200);

  //njets
  h.njets2_mu = new TH1F("2l1d_mu_njets","2l1d_njets (DisLep.at(0) is muon)",10, 0, 10); 

  //flavor classification
  h.flav2_mu = new TH1F("2l1d_mu_flav","2l1d_flav (DisLep.at(0) is muon), 0:eee 1:ee#mu 2:e#mu#mu 3:#mu#mu#mu 4:#mu#mue 5:#muee 6:e#mue 7:#mue#mu",8, 0, 8);
  
  
  
  //############################################################################################################################
  
  
}


