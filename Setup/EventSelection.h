void disp_ml::EventSelection(){
 
  //defining the selections:
  _2l1d = false, _1l2d = false, _3d = false;
  evt_dispml = false;
  evt_mumud  = false;
  evt_eed    = false;  
  evt_2LonZ  = false;
  evt_3L     = false;

  bool all_lep_isolated = true;
  for(int i=0; i<(int)lightLep.size(); i++){
    if(lightLep.at(i).reliso03>0.15) all_lep_isolated = false;
  }
  
    
  //2LonZ
  if((int)promptLepton.size()==2 && all_lep_isolated){
    float ll_invmass = (promptLepton.at(0).v+promptLepton.at(1).v).M();
    if(76.0<ll_invmass && ll_invmass<106.0 && metpt<40.0) evt_2LonZ = true;
  }

  
  //3L (exclusive)
  if((int)promptLepton.size()==3 && all_lep_isolated) evt_3L = true;

  //*****************************************************************************
  //all 3 displaced final states
  //*****************************************************************************
  if((int)promptLepton.size()>1 && (int)displacedLepton.size()>0 && all_lep_isolated)   _2l1d = true;
  if((int)promptLepton.size()>0 && (int)displacedLepton.size()>1 && all_lep_isolated)   _1l2d = true;
  if((int)promptLepton.size()>=0 && (int)displacedLepton.size()>2 && all_lep_isolated)  _3d   = true;

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
  
}//EventSelection()
