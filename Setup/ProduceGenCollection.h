// Gen Leptons (Electrons and Muons) Collection //

void disp_ml::GenLeptonArray(){

  for(unsigned int i=0; i<(*nGenPart); i++){
    Lepton temp;
    temp.status = GenPart_status[i];
    if(temp.status==1){
      temp.v.SetPtEtaPhiM(GenPart_pt[i],GenPart_eta[i],GenPart_phi[i],GenPart_mass[i]);
      temp.ind=i;
      temp.pdgid=GenPart_pdgId[i];
      temp.momid=MotherID(i,GenPart_genPartIdxMother[i]);
      
      int lepcharge=0;
      if(GenPart_pdgId[i]>0) lepcharge=-1;
      else lepcharge=1;
      temp.charge = lepcharge;
	
      bool passcutmu = abs(temp.pdgid)==13 && temp.v.Pt()>10 && abs(temp.v.Eta())<2.4;
      bool passcutel = abs(temp.pdgid)==11 && temp.v.Pt()>10 && abs(temp.v.Eta())<2.4;
      if(passcutmu){
	genMuon.push_back(temp);
      }

      if(passcutel){
	genElectron.push_back(temp);
      }
    }
  }

}
