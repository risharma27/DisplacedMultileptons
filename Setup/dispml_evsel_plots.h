void disp_ml::dispml_evsel_plots(){

  for(int ev=0; ev<(int)vec_evsel.size(); ev++){
    evsel = vec_evsel.at(ev);

    //********************************* Get Event Weight *********************************************************//
    //corrections are only applied on MC
    if(_data==0){
      float lep0SF = LeptonIDSF(myLep[evsel].at(0).id, myLep[evsel].at(0).v.Pt(), myLep[evsel].at(0).v.Eta());
      float lep1SF = LeptonIDSF(myLep[evsel].at(1).id, myLep[evsel].at(1).v.Pt(), myLep[evsel].at(1).v.Eta());
      float lep2SF = LeptonIDSF(myLep[evsel].at(2).id, myLep[evsel].at(2).v.Pt(), myLep[evsel].at(2).v.Eta());
      scalefactor = lep0SF * lep1SF * lep2SF;
	
      float e1=SingleLepTrigger_eff(myLep[evsel].at(0).id, myLep[evsel].at(0).v.Pt(), myLep[evsel].at(0).v.Eta());
      float e2=SingleLepTrigger_eff(myLep[evsel].at(1).id, myLep[evsel].at(1).v.Pt(), myLep[evsel].at(1).v.Eta());
      float e3=SingleLepTrigger_eff(myLep[evsel].at(2).id, myLep[evsel].at(2).v.Pt(), myLep[evsel].at(2).v.Eta());	  
      triggeff=1-((1-e1)*(1-e2)*(1-e3));
      
      evtwt = scalefactor * triggeff;
     
    }    
    //************************************************************************************************************//
    

    //************** Flavor Classification *************//

    if(abs(myLep[evsel].at(0).id)==13){
      if(abs(myLep[evsel].at(1).id)==13 && abs(myLep[evsel].at(2).id)==13)             h.flavor[evsel]->Fill(0.0, evtwt);       //mumumu
      else if(abs(myLep[evsel].at(1).id)==13 && abs(myLep[evsel].at(2).id)==11)        h.flavor[evsel]->Fill(1.0, evtwt);       //mumue
      else if(abs(myLep[evsel].at(1).id)==11 && abs(myLep[evsel].at(2).id)==13)        h.flavor[evsel]->Fill(2.0, evtwt);       //muemu
      else if(abs(myLep[evsel].at(1).id)==11 && abs(myLep[evsel].at(2).id)==11)        h.flavor[evsel]->Fill(3.0, evtwt);       //muee
    }	

    else if(abs(myLep[evsel].at(0).id)==11){
      if(abs(myLep[evsel].at(1).id)==11 && abs(myLep[evsel].at(2).id)==11)             h.flavor[evsel]->Fill(4.0, evtwt);       //eee
      else if(abs(myLep[evsel].at(1).id)==13 && abs(myLep[evsel].at(2).id)==11)        h.flavor[evsel]->Fill(5.0, evtwt);       //emue
      else if(abs(myLep[evsel].at(1).id)==11 && abs(myLep[evsel].at(2).id)==13)        h.flavor[evsel]->Fill(6.0, evtwt);       //eemu
      else if(abs(myLep[evsel].at(1).id)==13 && abs(myLep[evsel].at(2).id)==13)        h.flavor[evsel]->Fill(7.0, evtwt);       //emumu	
    }

    //***************************************************//

    
    float invmassl0l1 = (myLep[evsel].at(0).v+myLep[evsel].at(1).v).M();
    
    if((evsel==0 && invmassl0l1>12) || evsel==1 || evsel==2){
      h.dispml_h[evsel][0]->Fill(metpt, evtwt);
      float sum_pt = 0.0;
      for(int i=0; i<(int)myLep[evsel].size(); i++){
	sum_pt = sum_pt + myLep[evsel].at(i).v.Pt();
      }
      h.dispml_h[evsel][1]->Fill(sum_pt, evtwt);
      float imass = ((myLep[evsel].at(0).v + myLep[evsel].at(1).v) + myLep[evsel].at(2).v).M();
      h.dispml_h[evsel][2]->Fill(imass, evtwt);
      for(int j=3; j<6; j++){
	h.dispml_h[evsel][j]->Fill(myLep[evsel].at(j-3).v.Pt(), evtwt);
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
	h.dispml_h[evsel][index+p]->Fill(pt_ll[index], evtwt);
	h.dispml_h[evsel][index+p+1]->Fill(delR_ll[index], evtwt);
	h.dispml_h[evsel][index+p+2]->Fill(delPhi_ll[index], evtwt);
	h.dispml_h[evsel][index+p+3]->Fill(delphi_lmet[index], evtwt);
	h.dispml_h[evsel][index+p+4]->Fill(M_ll[index], evtwt);
	h.dispml_h[evsel][index+p+5]->Fill(transvmass[index], evtwt);
	p=p+5;
      }

	
      float jet_pt = 0.0;
      for(int i=0; i<(int)recoJet.size(); i++){
	jet_pt = jet_pt + recoJet.at(i).v.Pt();	
      }

      h.dispml_h[evsel][24]->Fill(jet_pt, evtwt);
      h.dispml_h[evsel][25]->Fill((int)recoJet.size(), evtwt);
      
      std::pair<vector<int>, vector<float>> result = dR_matching(myLep[evsel], recoJet, 0.05);
      vector<int> myLep_matchto_recoJet = result.first;
      vector<float> myLep_delRmin_recoJet = result.second;

      for(int i=0; i<(int)myLep_matchto_recoJet.size(); i++){
	int matchind=myLep_matchto_recoJet.at(i);
	float matchdR=myLep_delRmin_recoJet.at(i);
	if(matchind>-1){
	  h.dispml_h[evsel][i+26]->Fill(matchdR, evtwt);
	}
	else{
	  h.dispml_h[evsel][i+26]->Fill(99, evtwt);
	}
      }
	
      int q=29;
      for(int i=0; i<(int)myLep[evsel].size(); i++){
	h.dispml_h[evsel][i+q]->Fill(myLep[evsel].at(i).dxy, evtwt);
	h.dispml_h[evsel][i+q+1]->Fill(myLep[evsel].at(i).dz, evtwt);
	h.dispml_h[evsel][i+q+2]->Fill(myLep[evsel].at(i).ip3d, evtwt);
	h.dispml_h[evsel][i+q+3]->Fill(myLep[evsel].at(i).sip3d, evtwt);
	h.dispml_h[evsel][i+q+4]->Fill(myLep[evsel].at(i).reliso03, evtwt);
	q=q+4;
      }

      h.dispml_h[evsel][44]->Fill((int)bJet.size(), evtwt);
	  

      //******************************* 2l1d "flavor classified" analysis *******************************//
		
      if(evt_llmu){ //llmu selection
	h.dispml_h[0][45]->Fill((myLep[0].at(0).v+myLep[0].at(1).v+myLep[0].at(2).v).M());
	h.dispml_h[0][46]->Fill(metpt, evtwt);
	h.dispml_h[0][47]->Fill(myLep[0].at(2).dxy, evtwt);
	h.dispml_h[0][48]->Fill(myLep[0].at(2).ip3d, evtwt);
	h.dispml_h[0][49]->Fill(myLep[0].at(2).sip3d, evtwt);
	float delphi_l2met = delta_phi(metphi, myLep[0].at(2).v.Phi());
	h.dispml_h[0][50]->Fill(transv_mass(myLep[0].at(2).v.Pt(), metpt, delphi_l2met), evtwt);
	h.dispml_h[0][51]->Fill(myLep[0].at(0).v.DeltaPhi(myLep[0].at(1).v), evtwt);
	h.dispml_h[0][52]->Fill(myLep[0].at(0).v.DeltaR(myLep[0].at(1).v), evtwt);
	h.dispml_h[0][53]->Fill((myLep[0].at(0).v+myLep[0].at(1).v).M(), evtwt);
	h.dispml_h[0][54]->Fill(myLep[0].at(1).v.DeltaPhi(myLep[0].at(2).v), evtwt);
	h.dispml_h[0][55]->Fill(myLep[0].at(1).v.DeltaR(myLep[0].at(2).v), evtwt);
	h.dispml_h[0][56]->Fill((myLep[0].at(1).v+myLep[0].at(2).v).M(), evtwt);
	h.dispml_h[0][57]->Fill(myLep[0].at(2).v.DeltaPhi(myLep[0].at(0).v), evtwt);
	h.dispml_h[0][58]->Fill(myLep[0].at(2).v.DeltaR(myLep[0].at(0).v), evtwt);
	h.dispml_h[0][59]->Fill((myLep[0].at(2).v+myLep[0].at(0).v).M(), evtwt);	
      }

      else if(evt_lle){ //lle selection
	h.dispml_h[0][60]->Fill((myLep[0].at(0).v+myLep[0].at(1).v+myLep[0].at(2).v).M(), evtwt);
	h.dispml_h[0][61]->Fill(metpt, evtwt);
	h.dispml_h[0][62]->Fill(myLep[0].at(2).dxy, evtwt);
	h.dispml_h[0][63]->Fill(myLep[0].at(2).ip3d, evtwt);
	h.dispml_h[0][64]->Fill(myLep[0].at(2).sip3d, evtwt);
	float delphi_l2met = delta_phi(metphi, myLep[0].at(2).v.Phi());
	h.dispml_h[0][65]->Fill(transv_mass(myLep[0].at(2).v.Pt(), metpt, delphi_l2met), evtwt);
	h.dispml_h[0][66]->Fill(myLep[0].at(0).v.DeltaPhi(myLep[0].at(1).v), evtwt);
	h.dispml_h[0][67]->Fill(myLep[0].at(0).v.DeltaR(myLep[0].at(1).v), evtwt);
	h.dispml_h[0][68]->Fill((myLep[0].at(0).v+myLep[0].at(1).v).M(), evtwt);
	h.dispml_h[0][69]->Fill(myLep[0].at(1).v.DeltaPhi(myLep[0].at(2).v), evtwt);
	h.dispml_h[0][70]->Fill(myLep[0].at(1).v.DeltaR(myLep[0].at(2).v), evtwt);
	h.dispml_h[0][71]->Fill((myLep[0].at(1).v+myLep[0].at(2).v).M(), evtwt);
	h.dispml_h[0][72]->Fill(myLep[0].at(2).v.DeltaPhi(myLep[0].at(0).v), evtwt);
	h.dispml_h[0][73]->Fill(myLep[0].at(2).v.DeltaR(myLep[0].at(0).v), evtwt);
	h.dispml_h[0][74]->Fill((myLep[0].at(2).v+myLep[0].at(0).v).M(), evtwt);
      }

      //if(evt_mumud)  //mumud selection
      //if(evt_eed)    //eed selection


      //***************************************************************************************************//

    }//if((evsel==0 && invmassl0l1>12) || evsel==1 || evsel==2)
    
  }//for ev<vec_evsel.size()
	
}//dispml_evsel_plots
