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

