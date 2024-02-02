void disp_ml::other_evsel_plots(float wt){

  //2LonZ
  float invmass_ll = -1.0; float invmass_lll = -1.0;
  if((int)promptLepton.size()>1 && (promptLepton.at(0).id)*(promptLepton.at(1).id)==-169){  //considering Z->mumu events for now
    invmass_ll = (promptLepton.at(0).v+promptLepton.at(1).v).M();
    if(76.0<invmass_ll && invmass_ll<106.0 && metpt<40.0){
      h._2LonZ[0]->Fill(invmass_ll);
      h._2LonZ[1]->Fill(metpt);
    }
  }


  //mumu+d selection
  if((int)promptLepton.size()>1 && (int)displacedLepton.size()>0){
    if(fabs(promptLepton.at(0).id)==13 && fabs(promptLepton.at(1).id)==13){
      invmass_ll = (promptLepton.at(0).v+promptLepton.at(1).v).M();
      invmass_lll = (promptLepton.at(0).v+promptLepton.at(1).v+displacedLepton.at(0).v).M();
      h.mumud[0]->Fill(invmass_ll);
      h.mumud[1]->Fill(invmass_lll);
      h.mumud[2]->Fill(metpt);
    }
  }

  //ee+d
  if((int)promptLepton.size()>1 && (int)displacedLepton.size()>0){
    if(fabs(promptLepton.at(0).id)==11 && fabs(promptLepton.at(1).id)==11){
      invmass_ll = (promptLepton.at(0).v+promptLepton.at(1).v).M();
      invmass_lll = (promptLepton.at(0).v+promptLepton.at(1).v+displacedLepton.at(0).v).M();
      h.eed[0]->Fill(invmass_ll);
      h.eed[1]->Fill(invmass_lll);
      h.eed[2]->Fill(metpt);
    }
  }

  
  //3L
  float _3L_invmass_lll = -1.0, _3L_invmass_l0l1 = -1.0, _3L_invmass_l1l2 = -1.0,  _3L_invmass_l2l0 = -1.0, _3L_ht = 0.0, _3L_lt = 0.0, _3L_st = 0.0;
  if((int)promptLepton.size()==3){
    _3L_invmass_lll = (promptLepton.at(0).v+promptLepton.at(1).v+promptLepton.at(2).v).M();
    _3L_invmass_l0l1 = (promptLepton.at(0).v+promptLepton.at(1).v).M();
    _3L_invmass_l1l2 = (promptLepton.at(1).v+promptLepton.at(2).v).M();
    _3L_invmass_l2l0 = (promptLepton.at(2).v+promptLepton.at(0).v).M();
    h._3L[0]->Fill(_3L_invmass_lll);
    h._3L[1]->Fill(_3L_invmass_l0l1);
    h._3L[2]->Fill(_3L_invmass_l1l2);
    h._3L[3]->Fill(_3L_invmass_l2l0);
    h._3L[4]->Fill(metpt);
    h._3L[5]->Fill(promptLepton.at(0).v.Pt());
    h._3L[6]->Fill(promptLepton.at(1).v.Pt());
    h._3L[7]->Fill(promptLepton.at(2).v.Pt());
    for(int i=0; i<(int)promptLepton.size(); i++){
      _3L_lt = _3L_lt + promptLepton.at(i).v.Pt();
    }
    h._3L[8]->Fill(_3L_lt);
    h._3L[9]->Fill((int)recoJet.size());
    for(int i=0; i<(int)recoJet.size(); i++){
      _3L_ht = _3L_ht + recoJet.at(i).v.Pt();
    }
    h._3L[10]->Fill(_3L_ht);
    _3L_st = _3L_ht + _3L_lt + metpt;
    h._3L[11]->Fill(_3L_st);
  }
      


}//other_evsel_plots
