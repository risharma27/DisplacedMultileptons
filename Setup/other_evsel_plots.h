void disp_ml::other_evsel_plots(){

  //2LonZ
  float invmass_ll = -1.0;
  if(evt_2LonZ){
    invmass_ll = (promptLepton.at(0).v+promptLepton.at(1).v).M();
    
    //********************************* Get Event Weight *********************************************************//
    //corrections are only applied on MC
    if(_data==0){
      float lep0SF = LeptonIDSF(promptLepton.at(0).id, promptLepton.at(0).v.Pt(), promptLepton.at(0).v.Eta());
      float lep1SF = LeptonIDSF(promptLepton.at(1).id, promptLepton.at(1).v.Pt(), promptLepton.at(1).v.Eta());
      scalefactor = lep0SF * lep1SF;
	
      float e1=SingleLepTrigger_eff(promptLepton.at(0).id, promptLepton.at(0).v.Pt(), promptLepton.at(0).v.Eta());
      float e2=SingleLepTrigger_eff(promptLepton.at(1).id, promptLepton.at(1).v.Pt(), promptLepton.at(1).v.Eta());	  
      triggeff=1-((1-e1)*(1-e2));
      
      evtwt = scalefactor * triggeff;
      //**********************************************************************************************************//
    }
    
    h._2LonZ[0]->Fill(invmass_ll, evtwt);
    h._2LonZ[1]->Fill(metpt, evtwt);
    
  }//if(evt_2LonZ)

  
  //2LSS
  if(evt_2LSS){
    h._2LSS[0]->Fill(0);
    invmass_ll = (promptLepton.at(0).v+promptLepton.at(1).v).M();
    if(evt_mumu)      {h._2LSS[0]->Fill(1);}
    else if(evt_ee)   {h._2LSS[0]->Fill(2);}
    else if(evt_emu)  {h._2LSS[0]->Fill(3);}
    else if(evt_mue)  {h._2LSS[0]->Fill(4);}
    
    if(_data==0){
      float lep0SF = LeptonIDSF(promptLepton.at(0).id, promptLepton.at(0).v.Pt(), promptLepton.at(0).v.Eta());
      float lep1SF = LeptonIDSF(promptLepton.at(1).id, promptLepton.at(1).v.Pt(), promptLepton.at(1).v.Eta());
      scalefactor = lep0SF * lep1SF;
	
      float e1=SingleLepTrigger_eff(promptLepton.at(0).id, promptLepton.at(0).v.Pt(), promptLepton.at(0).v.Eta());
      float e2=SingleLepTrigger_eff(promptLepton.at(1).id, promptLepton.at(1).v.Pt(), promptLepton.at(1).v.Eta());	  
      triggeff=1-((1-e1)*(1-e2));
      
      evtwt = scalefactor * triggeff;
    }
    
    h._2LSS[1]->Fill(invmass_ll, evtwt);
    h._2LSS[2]->Fill(metpt, evtwt);

    if(evt_mumu)                {h._2LSS[3]->Fill(invmass_ll, evtwt);}
    else if(evt_ee)             {h._2LSS[4]->Fill(invmass_ll, evtwt);}
    else if(evt_emu || evt_mue) {h._2LSS[5]->Fill(invmass_ll, evtwt);}
    
  }//if(evt_2LSS)
  
  
  //3L
  float _3L_invmass_lll = -1.0, _3L_invmass_l0l1 = -1.0, _3L_invmass_l1l2 = -1.0,  _3L_invmass_l2l0 = -1.0, _3L_ht = 0.0, _3L_lt = 0.0, _3L_st = 0.0;
  if(evt_3L){
    //*************************************************************************************************************
    if(_data==0){
      float lep0SF = LeptonIDSF(promptLepton.at(0).id, promptLepton.at(0).v.Pt(), promptLepton.at(0).v.Eta());
      float lep1SF = LeptonIDSF(promptLepton.at(1).id, promptLepton.at(1).v.Pt(), promptLepton.at(1).v.Eta());
      float lep2SF = LeptonIDSF(promptLepton.at(2).id, promptLepton.at(2).v.Pt(), promptLepton.at(2).v.Eta());
      scalefactor = lep0SF * lep1SF * lep2SF;
	
      float e1=SingleLepTrigger_eff(promptLepton.at(0).id, promptLepton.at(0).v.Pt(), promptLepton.at(0).v.Eta());
      float e2=SingleLepTrigger_eff(promptLepton.at(1).id, promptLepton.at(1).v.Pt(), promptLepton.at(1).v.Eta());
      float e3=SingleLepTrigger_eff(promptLepton.at(2).id, promptLepton.at(2).v.Pt(), promptLepton.at(2).v.Eta());
      triggeff=1-((1-e1)*(1-e2)*(1-e3));
      
      evtwt = scalefactor * triggeff;
    //*************************************************************************************************************//
    }
     
    _3L_invmass_lll = (promptLepton.at(0).v+promptLepton.at(1).v+promptLepton.at(2).v).M();
    _3L_invmass_l0l1 = (promptLepton.at(0).v+promptLepton.at(1).v).M();
    _3L_invmass_l1l2 = (promptLepton.at(1).v+promptLepton.at(2).v).M();
    _3L_invmass_l2l0 = (promptLepton.at(2).v+promptLepton.at(0).v).M();
    h._3L[0]->Fill(_3L_invmass_lll, evtwt);
    h._3L[1]->Fill(_3L_invmass_l0l1, evtwt);
    h._3L[2]->Fill(_3L_invmass_l1l2, evtwt);
    h._3L[3]->Fill(_3L_invmass_l2l0, evtwt);
    h._3L[4]->Fill(metpt, evtwt);
    h._3L[5]->Fill(promptLepton.at(0).v.Pt(), evtwt);
    h._3L[6]->Fill(promptLepton.at(1).v.Pt(), evtwt);
    h._3L[7]->Fill(promptLepton.at(2).v.Pt(), evtwt);
    for(int i=0; i<(int)promptLepton.size(); i++){
      _3L_lt = _3L_lt + promptLepton.at(i).v.Pt();
    }
    h._3L[8]->Fill(_3L_lt, evtwt);
    h._3L[9]->Fill((int)recoJet.size(), evtwt);
    for(int i=0; i<(int)recoJet.size(); i++){
      _3L_ht = _3L_ht + recoJet.at(i).v.Pt();
    }
    h._3L[10]->Fill(_3L_ht, evtwt);
    _3L_st = _3L_ht + _3L_lt + metpt;
    h._3L[11]->Fill(_3L_st, evtwt);
    
  }//if(evt_3L)



  //cr_ttbar
  float cr_ttbar_lt = 0.0, cr_ttbar_ht = 0.0;
  if(cr_ttbar){
    h.cr_ttbar[0]->Fill(metpt);
    h.cr_ttbar[1]->Fill(promptLepton.at(0).v.Pt());
    h.cr_ttbar[2]->Fill(promptLepton.at(1).v.Pt());
    for(int i=0; i<(int)lightLep.size(); i++){
      cr_ttbar_lt = cr_ttbar_lt + lightLep.at(i).v.Pt();
    }
    h.cr_ttbar[3]->Fill(cr_ttbar_lt);
    h.cr_ttbar[4]->Fill((int)recoJet.size());
    h.cr_ttbar[5]->Fill((int)bJet.size());  
    for(int i=0; i<(int)recoJet.size(); i++){
      cr_ttbar_ht = cr_ttbar_ht + recoJet.at(i).v.Pt();
    }
    h.cr_ttbar[6]->Fill(cr_ttbar_ht);
    h.cr_ttbar_ht_rebinned->Fill(cr_ttbar_ht, evtwt);
    h.cr_ttbar[7]->Fill(transv_mass(promptLepton.at(0).v.Pt(), metpt, delta_phi(promptLepton.at(0).v.Phi(), metphi)));
    h.cr_ttbar[8]->Fill(transv_mass(promptLepton.at(1).v.Pt(), metpt, delta_phi(promptLepton.at(1).v.Phi(), metphi)));
    h.cr_ttbar[9]->Fill(promptLepton.at(0).reliso03);
    h.cr_ttbar[10]->Fill(promptLepton.at(1).reliso03);
    h.cr_ttbar[11]->Fill(delta_phi(promptLepton.at(0).v.Phi(), promptLepton.at(1).v.Phi()));
    h.cr_ttbar[12]->Fill((promptLepton.at(0).v+promptLepton.at(1).v).M());
    h.cr_ttbar[13]->Fill((recoJet.at(0).v+recoJet.at(1).v).M());
    h.cr_ttbar[14]->Fill(recoJet.at(0).v.Pt());
    h.cr_ttbar[15]->Fill(recoJet.at(1).v.Pt());
    
  }//if(cr_ttbar)

  //vr_ttbar
  float vr_ttbar_lt = 0.0, vr_ttbar_ht = 0.0;
  if(vr_ttbar){
    h.vr_ttbar[0]->Fill(metpt);
    h.vr_ttbar[1]->Fill(promptLepton.at(0).v.Pt());
    h.vr_ttbar[2]->Fill(promptLepton.at(1).v.Pt());
    for(int i=0; i<(int)lightLep.size(); i++){
      vr_ttbar_lt = vr_ttbar_lt + lightLep.at(i).v.Pt();
    }
    h.vr_ttbar[3]->Fill(vr_ttbar_lt);
    h.vr_ttbar[4]->Fill((int)recoJet.size());
    h.vr_ttbar[5]->Fill((int)bJet.size());  
    for(int i=0; i<(int)recoJet.size(); i++){
      vr_ttbar_ht = vr_ttbar_ht + recoJet.at(i).v.Pt();
    }
    h.vr_ttbar[6]->Fill(vr_ttbar_ht);
    h.vr_ttbar[7]->Fill(transv_mass(promptLepton.at(0).v.Pt(), metpt, delta_phi(promptLepton.at(0).v.Phi(), metphi)));
    h.vr_ttbar[8]->Fill(transv_mass(promptLepton.at(1).v.Pt(), metpt, delta_phi(promptLepton.at(1).v.Phi(), metphi)));
    h.vr_ttbar[9]->Fill(promptLepton.at(0).reliso03);
    h.vr_ttbar[10]->Fill(promptLepton.at(1).reliso03);
    h.vr_ttbar[11]->Fill(delta_phi(promptLepton.at(0).v.Phi(), promptLepton.at(1).v.Phi()));
    h.vr_ttbar[12]->Fill((promptLepton.at(0).v+promptLepton.at(1).v).M());
    h.vr_ttbar[13]->Fill((recoJet.at(0).v+recoJet.at(1).v).M());
    h.vr_ttbar[14]->Fill(recoJet.at(0).v.Pt());
    h.vr_ttbar[15]->Fill(recoJet.at(1).v.Pt());
    
  }//if(vr_ttbar)


  //cr_ttbar_2l1d
  float cr_ttbar_2l1d_lt = 0.0, cr_ttbar_2l1d_ht = 0.0;
  if(cr_ttbar_2l1d){

    //*************************************************************************************************************
    if(_data==0){
      float lep0SF = LeptonIDSF(promptLepton.at(0).id, promptLepton.at(0).v.Pt(), promptLepton.at(0).v.Eta());
      float lep1SF = LeptonIDSF(promptLepton.at(1).id, promptLepton.at(1).v.Pt(), promptLepton.at(1).v.Eta());
      float lep2SF = LeptonIDSF(displacedLepton.at(0).id, displacedLepton.at(0).v.Pt(), displacedLepton.at(0).v.Eta());
      scalefactor = lep0SF * lep1SF * lep2SF;
	
      float e1=SingleLepTrigger_eff(promptLepton.at(0).id, promptLepton.at(0).v.Pt(), promptLepton.at(0).v.Eta());
      float e2=SingleLepTrigger_eff(promptLepton.at(1).id, promptLepton.at(1).v.Pt(), promptLepton.at(1).v.Eta());
      float e3=SingleLepTrigger_eff(displacedLepton.at(0).id, displacedLepton.at(0).v.Pt(), displacedLepton.at(0).v.Eta());
      triggeff=1-((1-e1)*(1-e2)*(1-e3));
      
      evtwt = scalefactor * triggeff;
    }
    //*************************************************************************************************************//
    
    h.cr_ttbar_2l1d[0]->Fill(metpt, evtwt);
    h.cr_ttbar_2l1d[1]->Fill(promptLepton.at(0).v.Pt(), evtwt);
    h.cr_ttbar_2l1d[2]->Fill(promptLepton.at(1).v.Pt(), evtwt);
    for(int i=0; i<(int)lightLep.size(); i++){
      cr_ttbar_2l1d_lt = cr_ttbar_2l1d_lt + lightLep.at(i).v.Pt();
    }
    h.cr_ttbar_2l1d[3]->Fill(cr_ttbar_2l1d_lt, evtwt);
    h.cr_ttbar_2l1d[4]->Fill((int)recoJet.size(), evtwt);
    h.cr_ttbar_2l1d[5]->Fill((int)bJet.size(), evtwt);  
    for(int i=0; i<(int)recoJet.size(); i++){
      cr_ttbar_2l1d_ht = cr_ttbar_2l1d_ht + recoJet.at(i).v.Pt();
    }
    h.cr_ttbar_2l1d[6]->Fill(cr_ttbar_2l1d_ht, evtwt);
    h.cr_ttbar_2l1d_ht_rebinned->Fill(cr_ttbar_2l1d_ht, evtwt);
    h.cr_ttbar_2l1d[7]->Fill(transv_mass(promptLepton.at(0).v.Pt(), metpt, delta_phi(promptLepton.at(0).v.Phi(), metphi)), evtwt);
    h.cr_ttbar_2l1d[8]->Fill(transv_mass(promptLepton.at(1).v.Pt(), metpt, delta_phi(promptLepton.at(1).v.Phi(), metphi)), evtwt);
    h.cr_ttbar_2l1d[9]->Fill(promptLepton.at(0).reliso03, evtwt);
    h.cr_ttbar_2l1d[10]->Fill(promptLepton.at(1).reliso03, evtwt);
    h.cr_ttbar_2l1d[11]->Fill(delta_phi(promptLepton.at(0).v.Phi(), promptLepton.at(1).v.Phi()), evtwt);
    h.cr_ttbar_2l1d[12]->Fill((promptLepton.at(0).v+promptLepton.at(1).v).M(), evtwt);
    h.cr_ttbar_2l1d[13]->Fill((recoJet.at(0).v+recoJet.at(1).v).M(), evtwt);
    h.cr_ttbar_2l1d[14]->Fill(recoJet.at(0).v.Pt(), evtwt);
    h.cr_ttbar_2l1d[15]->Fill(recoJet.at(1).v.Pt(), evtwt);   
			      
  }//if(cr_ttbar_2l1d)

   //vr_ttbar_2l1d
  float vr_ttbar_2l1d_lt = 0.0, vr_ttbar_2l1d_ht = 0.0;
  if(vr_ttbar_2l1d){
    
    //*************************************************************************************************************
    if(_data==0){
      float lep0SF = LeptonIDSF(promptLepton.at(0).id, promptLepton.at(0).v.Pt(), promptLepton.at(0).v.Eta());
      float lep1SF = LeptonIDSF(promptLepton.at(1).id, promptLepton.at(1).v.Pt(), promptLepton.at(1).v.Eta());
      float lep2SF = LeptonIDSF(displacedLepton.at(0).id, displacedLepton.at(0).v.Pt(), displacedLepton.at(0).v.Eta());
      scalefactor = lep0SF * lep1SF * lep2SF;
	
      float e1=SingleLepTrigger_eff(promptLepton.at(0).id, promptLepton.at(0).v.Pt(), promptLepton.at(0).v.Eta());
      float e2=SingleLepTrigger_eff(promptLepton.at(1).id, promptLepton.at(1).v.Pt(), promptLepton.at(1).v.Eta());
      float e3=SingleLepTrigger_eff(displacedLepton.at(0).id, displacedLepton.at(0).v.Pt(), displacedLepton.at(0).v.Eta());
      triggeff=1-((1-e1)*(1-e2)*(1-e3));
      
      evtwt = scalefactor * triggeff;
    }
    //*************************************************************************************************************//
    
    h.vr_ttbar_2l1d[0]->Fill(metpt, evtwt);
    h.vr_ttbar_2l1d[1]->Fill(promptLepton.at(0).v.Pt(), evtwt);
    h.vr_ttbar_2l1d[2]->Fill(promptLepton.at(1).v.Pt(), evtwt);
    for(int i=0; i<(int)lightLep.size(); i++){
      vr_ttbar_2l1d_lt = vr_ttbar_2l1d_lt + lightLep.at(i).v.Pt();
    }
    h.vr_ttbar_2l1d[3]->Fill(vr_ttbar_2l1d_lt, evtwt);
    h.vr_ttbar_2l1d[4]->Fill((int)recoJet.size(), evtwt);
    h.vr_ttbar_2l1d[5]->Fill((int)bJet.size(), evtwt);  
    for(int i=0; i<(int)recoJet.size(); i++){
      vr_ttbar_2l1d_ht = vr_ttbar_2l1d_ht + recoJet.at(i).v.Pt();
    }
    h.vr_ttbar_2l1d[6]->Fill(vr_ttbar_2l1d_ht, evtwt);
    h.vr_ttbar_2l1d[7]->Fill(transv_mass(promptLepton.at(0).v.Pt(), metpt, delta_phi(promptLepton.at(0).v.Phi(), metphi)), evtwt);
    h.vr_ttbar_2l1d[8]->Fill(transv_mass(promptLepton.at(1).v.Pt(), metpt, delta_phi(promptLepton.at(1).v.Phi(), metphi)), evtwt);
    h.vr_ttbar_2l1d[9]->Fill(promptLepton.at(0).reliso03, evtwt);
    h.vr_ttbar_2l1d[10]->Fill(promptLepton.at(1).reliso03, evtwt);
    h.vr_ttbar_2l1d[11]->Fill(delta_phi(promptLepton.at(0).v.Phi(), promptLepton.at(1).v.Phi()), evtwt);
    h.vr_ttbar_2l1d[12]->Fill((promptLepton.at(0).v+promptLepton.at(1).v).M(), evtwt);
    h.vr_ttbar_2l1d[13]->Fill((recoJet.at(0).v+recoJet.at(1).v).M(), evtwt);
    h.vr_ttbar_2l1d[14]->Fill(recoJet.at(0).v.Pt(), evtwt);
    h.vr_ttbar_2l1d[15]->Fill(recoJet.at(1).v.Pt(), evtwt);   
			      
  }//if(vr_ttbar_2l1d)


  //cr_wjets
  float cr_wjets_lt = 0.0, cr_wjets_ht = 0.0;
  if(cr_wjets){
    h.cr_wjets[0]->Fill(metpt);
    h.cr_wjets[1]->Fill(promptLepton.at(0).v.Pt());
    for(int i=0; i<(int)lightLep.size(); i++){
      cr_wjets_lt = cr_wjets_lt + lightLep.at(i).v.Pt();
    }
    h.cr_wjets[2]->Fill(cr_wjets_lt);
    h.cr_wjets[3]->Fill((int)recoJet.size());
    for(int i=0; i<(int)recoJet.size(); i++){
      cr_wjets_ht = cr_wjets_ht + recoJet.at(i).v.Pt();
    }
    h.cr_wjets[4]->Fill(cr_wjets_ht);
    h.cr_wjets[5]->Fill(transv_mass(promptLepton.at(0).v.Pt(), metpt, delta_phi(promptLepton.at(0).v.Phi(), metphi)));
    h.cr_wjets[6]->Fill(promptLepton.at(0).reliso03);
    h.cr_wjets[7]->Fill(delta_phi(promptLepton.at(0).v.Phi(), recoJet.at(0).v.Phi()));
    h.cr_wjets[8]->Fill(recoJet.at(0).v.Pt());
  }//if(cr_wjets)

   //cr_wjets_1l2d
  float cr_wjets_1l2d_lt = 0.0, cr_wjets_1l2d_ht = 0.0;
  if(cr_wjets_1l2d){
    h.cr_wjets_1l2d[0]->Fill(metpt);
    h.cr_wjets_1l2d[1]->Fill(promptLepton.at(0).v.Pt());
    for(int i=0; i<(int)lightLep.size(); i++){
      cr_wjets_1l2d_lt = cr_wjets_1l2d_lt + lightLep.at(i).v.Pt();
    }
    h.cr_wjets_1l2d[2]->Fill(cr_wjets_1l2d_lt);
    h.cr_wjets_1l2d[3]->Fill((int)recoJet.size());
    for(int i=0; i<(int)recoJet.size(); i++){
      cr_wjets_1l2d_ht = cr_wjets_1l2d_ht + recoJet.at(i).v.Pt();
    }
    h.cr_wjets_1l2d[4]->Fill(cr_wjets_1l2d_ht);
    h.cr_wjets_1l2d[5]->Fill(transv_mass(promptLepton.at(0).v.Pt(), metpt, delta_phi(promptLepton.at(0).v.Phi(), metphi)));
    h.cr_wjets_1l2d[6]->Fill(promptLepton.at(0).reliso03);
    h.cr_wjets_1l2d[7]->Fill(delta_phi(promptLepton.at(0).v.Phi(), recoJet.at(0).v.Phi()));
    h.cr_wjets_1l2d[8]->Fill(recoJet.at(0).v.Pt());
  }//if(cr_wjets_1l2d)
  
}//other_evsel_plots


