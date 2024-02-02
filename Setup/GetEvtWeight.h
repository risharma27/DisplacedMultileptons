float disp_ml::GetEvtWeight(){

  float wt = 1.0; //default value

  if(_data==0){  //these corrections are applied to MC only.
    float scalefactor = 1.0;
    float triggeff = 1.0;

    if(evt_dispml){
      float lep0SF = LeptonIDSF(myLep[evsel].at(0).id, myLep[evsel].at(0).v.Pt(), myLep[evsel].at(0).v.Eta());
      float lep1SF = LeptonIDSF(myLep[evsel].at(1).id, myLep[evsel].at(1).v.Pt(), myLep[evsel].at(1).v.Eta());
      float lep2SF = LeptonIDSF(myLep[evsel].at(2).id, myLep[evsel].at(2).v.Pt(), myLep[evsel].at(2).v.Eta());
      scalefactor = lep0SF * lep1SF * lep2SF;

      float e1=SingleLepTrigger_eff(myLep[evsel].at(0).id, myLep[evsel].at(0).v.Pt(), myLep[evsel].at(0).v.Eta());
      float e2=SingleLepTrigger_eff(myLep[evsel].at(1).id, myLep[evsel].at(1).v.Pt(), myLep[evsel].at(1).v.Eta());
      float e3=SingleLepTrigger_eff(myLep[evsel].at(2).id, myLep[evsel].at(2).v.Pt(), myLep[evsel].at(2).v.Eta());	  
      triggeff=1-((1-e1)*(1-e2)*(1-e3));
    }
  
    wt = scalefactor * triggeff;
  
    //h.evtweight[evsel][0]->Fill(scalefactor);
    //h.evtweight[evsel][1]->Fill(triggeff);
    //h.evtweight[evsel][2]->Fill(evtwt);
 
  }//if(_data=0)

  return wt;
  
}
