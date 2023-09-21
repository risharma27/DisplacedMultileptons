//Reco Leptons (Electrons and Muons) Collection//

void disp_ml::RecoLeptonArray(){
  
  //Muons
  
  int nmu = 0;                       
  for(unsigned int i=0; i<(*nMuon); i++){
    Lepton temp;                      
    temp.v.SetPtEtaPhiM(Muon_pt[i],Muon_eta[i],Muon_phi[i],0.105); //the muon mass in GeV is 0.105
    temp.id = -13*Muon_charge[i];    //pdgID for mu- = 13, pdgID for mu+ = -13  
    temp.ind = i;
    temp.reliso03 = Muon_pfRelIso03_all[i];
    temp.reliso04 = Muon_pfRelIso04_all[i];
    temp.dxy = Muon_dxy[i];
    temp.dz = Muon_dz[i];
    
    bool ptetacut = temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4;
    bool passcut_mediummu =  ptetacut && Muon_pfRelIso03_all[i]<0.15 && Muon_mediumId[i];
    bool promptmuon = passcut_mediummu && fabs(Muon_dxy[i])<0.05 && fabs(Muon_dz[i])<0.1;
    bool displacedmuon = passcut_mediummu && fabs(Muon_dxy[i])>0.05 && fabs(Muon_dz[i])<10;
      
    if(passcut_mediummu){
      recoMuon.push_back(temp);
      recoLepton.push_back(temp);
    }

    if(promptmuon){
      promptLepton.push_back(temp);
    }

    if(displacedmuon){
      displacedLepton.push_back(temp);
    }
  } //Muons


  //Electrons

  int nel=0;
  for(unsigned int i=0; i<(*nElectron); i++){
    Lepton temp;
    temp.v.SetPtEtaPhiM(Electron_pt[i],Electron_eta[i],Electron_phi[i],Electron_mass[i]);
    temp.id = -11*Electron_charge[i];
    temp.ind = i; 
    temp.reliso03 = Electron_pfRelIso03_all[i]; 
    temp.dxy = Electron_dxy[i];
    temp.dz = Electron_dz[i];
   
    bool ptetacut = temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4;
    passcut_mediumel = ptetacut && Electron_pfRelIso03_all[i]<0.15 && Electron_cutBased[i]>2;
    bool promptelectron = passcut_mediumel && fabs(Electron_dxy[i])<0.05 && fabs(Electron_dz[i])<0.1;
     bool displacedelectron = passcut_mediumel && fabs(Electron_dxy[i]).0.05;
      					
    if(passcut_mediumel){
      recoElectron.push_back(temp);
      recoLepton.push_back(temp);
    }

    if(promptelectron){
      promptLepton.push_back(temp);
    }

    if(displacedelectron){
      displacedLepton.push_back(temp);
    }
  } //Electrons
    
}
