void disp_ml::EventSelection(){
 
  //defining the selections:
  _2l1d = false, _1l2d = false, _3d = false;
  evt_dispml = false;
  evt_mumud  = false;
  evt_eed    = false;
  evt_llmu   = false;
  evt_lle    = false;
  evt_emud   = false;
  evt_lee    = false;
  evt_lmumu  = false;
  evt_lemu   = false;
  evt_2LonZ  = false;
  evt_2LSS   = false;
  evt_mumu = false, evt_ee = false, evt_emu = false, evt_mue = false;
  evt_3L     = false;

  //Control Regions
  cr_ttbar_2l1d = false;  vr_ttbar_2l1d = false;
  cr_qcd_2l1d   = false;  cr_qcd_1l2d = false; 
  cr_wjets_1l2d = false;
  
  bool all_lep_isolated = true;
  for(int i=0; i<(int)lightLep.size(); i++){
    if(lightLep.at(i).reliso03>0.15){
      all_lep_isolated = false;
      break;
    }
  }

  //HT (scalar sum of all jet pT) of the event
  float evt_ht = 0.0;
  for(int i=0; i<(int)recoJet.size(); i++){
    evt_ht = evt_ht + recoJet.at(i).v.Pt();
  }
  
 
  //2LonZ (L=mu only)
  if((int)promptLepton.size()==2 && all_lep_isolated && abs(promptLepton.at(0).id)==13 && abs(promptLepton.at(1).id)==13){
    if(promptLepton.at(0).id == -promptLepton.at(1).id){
      float ll_invmass = (promptLepton.at(0).v+promptLepton.at(1).v).M();
      if(76.0<ll_invmass && ll_invmass<106.0 && metpt<40.0) evt_2LonZ = true;
    }
  }

  //2LSS (exclusive)
  if((int)promptLepton.size()==2 && all_lep_isolated){
    if((promptLepton.at(0).charge == promptLepton.at(1).charge) && evt_ht<500.0){
      h._2LSS[6]->Fill(evt_ht);
      evt_2LSS = true;
      if(abs(promptLepton.at(0).id)==13 && abs(promptLepton.at(1).id)==13)       {evt_mumu = true;}
      else if(abs(promptLepton.at(0).id)==11 && abs(promptLepton.at(1).id)==11)  {evt_ee   = true;}
      else if(abs(promptLepton.at(0).id)==11 && abs(promptLepton.at(1).id)==13)  {evt_emu  = true;}
      else if(abs(promptLepton.at(0).id)==13 && abs(promptLepton.at(1).id)==11)  {evt_mue  = true;}
    }
  }

  
  //3L (exclusive)
  if((int)promptLepton.size()==3 && all_lep_isolated){
    if(promptLepton.at(0).v.Pt()>15.0 && promptLepton.at(1).v.Pt()>15.0 && promptLepton.at(2).v.Pt()>15.0) evt_3L = true;
  }
  
  
  //*****************************************************************************
  //all 3 displaced final states
  //*****************************************************************************
  if((int)promptLepton.size()>1 && (int)displacedLepton.size()>0 && all_lep_isolated)                 _2l1d = true;  //2l1d
  if((int)promptLepton.size()>0 && (int)displacedLepton.size()>1 && all_lep_isolated && metpt>20.0)   _1l2d = true;  //1l2d
  if((int)promptLepton.size()>=0 && (int)displacedLepton.size()>2)                                    _3d   = true;  //3d

  vec_evsel.clear();
     
  for(int i=0; i<3; i++){
    myLep[i].clear();  //clearing myLep[evsel] for each evsel.
  }
	
  if(_2l1d){
    evt_2l1d.push_back(nEvtTotal);
    //h.nevsel->Fill(1);
    vec_evsel.push_back(0);
    myLep[0].push_back(promptLepton.at(0));
    myLep[0].push_back(promptLepton.at(1));
    myLep[0].push_back(displacedLepton.at(0));
  }
    
  if(_1l2d){
    evt_1l2d.push_back(nEvtTotal);
    //h.count_1l2d->Fill(1);
    //h.nevsel->Fill(2);
    vec_evsel.push_back(1);
    myLep[1].push_back(promptLepton.at(0));
    myLep[1].push_back(displacedLepton.at(0));
    myLep[1].push_back(displacedLepton.at(1));
  }
    
  if(_3d){
    evt_3d.push_back(nEvtTotal);
    //h.nevsel->Fill(3);
    vec_evsel.push_back(2);
    myLep[2].push_back(displacedLepton.at(0));
    myLep[2].push_back(displacedLepton.at(1));
    myLep[2].push_back(displacedLepton.at(2));
  }
  
  if((int)vec_evsel.size()>0){
    evt_dispml = true;
    nEvtPass++;
    h.nevt->Fill(3);
  }

  //******************************************************************************

  //flavor classified selections in 2l1d
  if(_2l1d){
    if(abs(promptLepton.at(0).id)==13 && abs(promptLepton.at(1).id)==13)      evt_mumud = true; //mumud
    if(abs(promptLepton.at(0).id)==11 && abs(promptLepton.at(1).id)==11)      evt_eed   = true; //eed
    if((abs(promptLepton.at(0).id)==11 && abs(promptLepton.at(1).id)==13) || (abs(promptLepton.at(0).id)==13 && abs(promptLepton.at(1).id)==11)) evt_emud = true; //emud || mued
    if(abs(displacedLepton.at(0).id)==13) evt_llmu = true;
    if(abs(displacedLepton.at(0).id)==11) evt_lle = true;
  }

  //flavor classified selections in 1l2d
  if(_1l2d){
    h.count_1l2d->Fill(0);
    if(abs(displacedLepton.at(0).id)==11 && abs(displacedLepton.at(1).id)==11) evt_lee = true; //lee
    else if(abs(displacedLepton.at(0).id)==13 && abs(displacedLepton.at(1).id)==13) evt_lmumu = true; //lmumu
    else if((abs(displacedLepton.at(0).id)==13 && abs(displacedLepton.at(1).id)==11) || (abs(displacedLepton.at(0).id)==11 && abs(displacedLepton.at(1).id)==13)) evt_lemu = true; //lmue or lemu
  }

  
  if(evt_lemu)  h.count_1l2d->Fill(1);
  if(evt_lee)   h.count_1l2d->Fill(2);
  if(evt_lmumu) h.count_1l2d->Fill(3);
  
  //************* Control Regions *****************************************//

  //TTBar Control Region for 2l1d event selection

  //targeting TTBarToDileptonic final state
  if(_2l1d){
    if((int)recoJet.size()>1 && metpt>30.0){
      
      if(promptLepton.at(0).v.Pt()>25 && promptLepton.at(1).v.Pt()>20 && (promptLepton.at(0).charge == -promptLepton.at(1).charge) && (int)bJet.size()==1){
	float invmass_l0l1 = (promptLepton.at(0).v+promptLepton.at(1).v).M();
	if(invmass_l0l1>20.0){
	  //cout<<promptLepton.at(0).id<<" "<<promptLepton.at(1).id<<endl;
	  if(abs(promptLepton.at(0).id) == abs(promptLepton.at(1).id)){
	    if(invmass_l0l1<76.0 || invmass_l0l1>106.0) cr_ttbar_2l1d = true;  //rejecting Drell-Yan background
	  }
	  else cr_ttbar_2l1d = true;
	}
      }
    }
  }//if(_2l1d)
  

  //TTBar Validation Region
  if(_2l1d){ 
    if((int)recoJet.size()>1 && metpt>30.0){
      if(promptLepton.at(0).v.Pt()>25 && promptLepton.at(1).v.Pt()>20 && (promptLepton.at(0).charge == -promptLepton.at(1).charge) && (int)bJet.size()>1){
	float invmass_l0l1 = (promptLepton.at(0).v+promptLepton.at(1).v).M();
	if(invmass_l0l1>20.0){
	  //cout<<promptLepton.at(0).id<<" "<<promptLepton.at(1).id<<endl;
	  if(abs(promptLepton.at(0).id) == abs(promptLepton.at(1).id)){
	    if(invmass_l0l1<76.0 || invmass_l0l1>106.0) vr_ttbar_2l1d = true;  //rejecting Drell-Yan background
	  }
	  else vr_ttbar_2l1d = true;
	}
      }
    }
  }//if(_2l1d)
  
  
  //QCD Control Region
  
  /*
  if(_1l2d){
    if(promptLepton.at(0).reliso03>1.0 && displacedLepton.at(0).reliso03>1.0 && displacedLepton.at(1).reliso03>1.0){
      if(5.0<displacedLepton.at(0).sip3d && displacedLepton.at(0).sip3d<=25.0)  cr_1l2d_qcd = true;
      else if(displacedLepton.at(0).sip3d>25.0)                                 vr_1l2d_qcd = true;
    }
  } 
  */

  if(_2l1d){
    if(metpt<40.0 && displacedLepton.at(0).reliso03>0.15 && promptLepton.at(0).sip3d>10)      cr_qcd_2l1d   = true;
  }
 
  if(_1l2d){
    if(metpt<40.0 && displacedLepton.at(0).reliso03>0.15)      cr_qcd_1l2d   = true;
    else if(metpt>40.0)                                        cr_wjets_1l2d = true;
  }
 

}//EventSelection()
