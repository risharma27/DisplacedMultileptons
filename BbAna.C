#define BbAna_cxx
#include "BbAna.h"
#include "SmearingClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>

void BbAna::Begin()
{
  GEV = 1000.;
  MEV2GEV = 0.001;

  _HstFile = new TFile(_HstFileName,"recreate");
  BookHistograms();
  for(int i=0; i<20; i++) npasscuts[i]=0;
  for(int i=0;i<4;i++){
    bb_nev[i] = 0; bb_wgt[i] = 0.;
  }
}
void BbAna::End()
{
  cout << endl << " Total events processed =  " << nEvtTotal << endl;

  
  

  _HstFile->Write();
  _HstFile->Close();
 
  ofstream fout(_SumFileName);
  fout<<"Total_Events                      = "<<npasscuts[0]<<endl;
  fout<<"Total_Event_Weight                = "<<nEvtWeight<<endl;
  fout<<" ============================================  "<<endl;
  fout<<"Event counts and weights in Regions 0,1,2,3=A,B,C,D "<<endl;
  for(int i=0; i<4; i++){
    fout<<"Region "<<i<<": NEvents = "<<bb_nev[i]<<"  EventWeight = "<<bb_wgt[i]<<endl;
    cout<<"Region "<<i<<": NEvents = "<<bb_nev[i]<<"  EventWeight = "<<bb_wgt[i]<<endl;
  }
}

void BbAna::Loop(int nevents)
{
//   In a ROOT session, you can do:
//      Root > .L BbAna.C
//      Root > BbAna t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   if(nevents>0) nentries = nevents;
   nEvtTotal = 0;
   nEvtWeight = 0;
   cout << endl << " Total events    in     = " << nentries  << endl;

   //Initialize the smearing class once.
   SmearingClass mcp_smear;
   mcp_smear.UseScale(1);

   //fChain->SetCacheSize(10000000);
   //fChain->AddBranchToCache("*");
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      nEvtTotal++; npasscuts[0]++;
      nEvtWeight+=Weight;
      if(nEvtTotal%100000==0) 
	cout << nEvtTotal << "  events processed so far " << endl;

      //START MY CODE

      nlepsel = 0; nelsel = 0; nmusel = 0;       NTRK_SEL = 0;
      if(NVertex > 0){
	//h.nvtxtrack->Fill(VertexNtrk->at(0));
	if(VertexNtrk->at(0)<5) continue;
      }
      else if(NVertex <= 0)
	continue;
      npasscuts[1]++;

      //if(NBadJets > 0 && isData) continue;
      int nele_pass = 0;
      for(int i=0; i<NElectrons; i++){
	double elet = 0.;
	if( ElectronSctHits->at(i)+ElectronPixHits->at(i) < 4 ) elet = ElectronClusEt->at(i);
	else elet = ElectronClusE->at(i)/cosh(ElectronEta->at(i));
	bool passele = pass_electron_cut(i) && elet>15*GEV;
	if(passele){
	  elsel[nelsel].mom.SetPtEtaPhiE(ElectronPt->at(i),ElectronEta->at(i),ElectronPhi->at(i),ElectronClusE->at(i));
	  elsel[nelsel].id = 11*ElectronCharge->at(i);
	  elsel[nelsel].ind = i;
	  lepsel[nlepsel] = elsel[nelsel];
	  nelsel++;
	  nlepsel++;
	}

	if(passele){ nele_pass++;
	}
      }
      int nmu_pass = 0;
      for(int i=0; i<NMuons; i++){
	bool passmu = pass_muon_cut(i) && MuonPt->at(i)>15*GEV;
	/* get new charge by accounting for mixing */
	int new_charge = MuonCharge->at(i);
	//cout<<"Mixing On = "<<isMixingOn<<endl;
	h.mumotherid->Fill(abs(MuonMotherID->at(i)));
	if(isMixingOn){
	  //  if parent is B0_s, then mix 0.49927% of events. PDGID=531
	  //  if parent is B0,   then mix 0.1863% of events.  PDGID=511
	  TRandom3 randgen1; 
	  TRandom3 randgen2; 
	  if(abs(MuonMotherID->at(i))==531){
	    randgen1.SetSeed(Event+(i+1)*10); if(randgen1.Rndm() < 0.49927){ new_charge = -1*new_charge; }
	  }
	  else if(abs(MuonMotherID->at(i))==511){
	    randgen1.SetSeed(Event+(i+2)*11); if(randgen1.Rndm() < 0.1863){ new_charge = -1*new_charge; }
	  }
	}
	double myPT = MuonPt->at(i);         double myE = MuonE->at(i);
	if(passmu){
	  // Smear the muon pT if this is MC.
	  if(!isData){
	      double eta1 = MuonEta->at(i);         double ptcb1 = MuonPt->at(i);
	      double ptms1 = MuonPtmsextrap->at(i); double ptid1 = MuonPtid->at(i);
	      mcp_smear.SetSeed(EventNumber);
	      mcp_smear.Event(ptms1,ptid1,ptcb1,eta1);
	      double ptcb1_smear = mcp_smear.pTCB();
	      myPT = ptcb1_smear;
	      myE = sqrt(pow(myPT*cosh(MuonEta->at(i)),2)+pow(105.7,2));
	      // ----------- myPT is smeared
	  }
	  musel[nmusel].mom.SetPtEtaPhiE(myPT,MuonEta->at(i),MuonPhi->at(i),myE);
	  //musel[nmusel].id = 13*MuonCharge->at(i);
	  musel[nmusel].id = 13*new_charge;
	  musel[nmusel].ind = i;
	  lepsel[nlepsel] = musel[nmusel];
	  nmusel++;
	  nlepsel++;
	}
	if(passmu){ nmu_pass++;
	}
      }
      int ntracks_pass = 0;
      for(int i=0; i<NTracks; i++){
	h.trackpt->Fill(TrkPt->at(i)*MEV2GEV);
	bool passtrack = pass_track_cut(i);
	if(passtrack){ ntracks_pass++; NTRK_SEL++; }
      }

      //sort the arrays lepsel, elsel and musel according to pt/et
      Sort();
      //doTagAndProbeIsolation();
      //doEvtCountIsolation();
      //int tmp = doWFakeRate();
      int tmp = dobbPlots();
   }//Loop over events
}
// =======================================================================================================================
void BbAna::Sort()
{
  //sort lepsel;
  Lepton temp;
  for(int i=0; i<nlepsel-1; i++){
    for(int j=i+1; j<nlepsel; j++){
      if(lepsel[i].mom.Pt()<lepsel[j].mom.Pt()){
	temp = lepsel[i];
	lepsel[i] = lepsel[j];
	lepsel[j] = temp;
      }
    }
  }
  for(int i=0; i<nelsel-1; i++){
    for(int j=i+1; j<nelsel; j++){
      if(elsel[i].mom.Pt()<elsel[j].mom.Pt()){
	temp = elsel[i];
	elsel[i] = elsel[j];
	elsel[j] = temp;
      }
    }
  }
  for(int i=0; i<nmusel-1; i++){
    for(int j=i+1; j<nmusel; j++){
      if(musel[i].mom.Pt()<musel[j].mom.Pt()){
	temp = musel[i];
	musel[i] = musel[j];
	musel[j] = temp;
      }
    }
  }

}
bool BbAna::pass_isolation(int ind, int id, double pt)
{
  bool result = false;
  if(id==13){//test muon isolation
    if(Muonptcone20->at(ind)/pt<0.2) result = true;
    //if(Muonptcone20->at(ind)<2*GEV) result = true;
  }
  else if(id==11){//test electron isolation
    double elet = 0.;
    if( ElectronSctHits->at(ind)+ElectronPixHits->at(ind) < 4 ) elet = ElectronClusEt->at(ind);
    else elet = ElectronClusE->at(ind)/cosh(ElectronEta->at(ind));
    if(Electronetcone20->at(ind)/elet<0.2) result = true;
  }
  return result;
}
float BbAna::delta_phi(float phi1, float phi2)
{
  //phi1 = TVector2::Phi_0_2pi(phi1);
  //phi2 = TVector2::Phi_0_2pi(phi2);
    float dphi=fabs(phi1-phi2);
    if ( dphi>M_PI ) dphi = 2.*M_PI - dphi;
    //dphi = dphi*180/M_PI;
    return dphi;
}

float BbAna::delta_eta(float eta1, float eta2)
{
    return fabs(eta1-eta2);
}

float BbAna::delta_r(float eta1, float phi1, float eta2, float phi2)
{
    float dphi=delta_phi(phi1,phi2);
    return sqrt(pow((eta1-eta2),2)+pow(dphi,2));
}
float BbAna::inv_mass(const float p1[4], const float p2[4])
{
    float en = p1[0] + p2[0];
    float px = p1[1] + p2[1];
    float py = p1[2] + p2[2];
    float pz = p1[3] + p2[3];

    float mass2 = pow(en,2)-pow(px,2)-pow(py,2)-pow(pz,2);
    return sqrt(mass2);

}

float BbAna::inv_mass(const vector<TLorentzVector> &v)
{
    TLorentzVector vsum;
    for ( vector<TLorentzVector>::const_iterator i=v.begin(); i!=v.end(); i++ )
        vsum += *i;
    return vsum.M();
}
float BbAna::transverse_mass_pt_phi(float pt1, float phi1, float pt2, float phi2)
{
    return sqrt(2.*pt1*pt2*delta_phi(phi1,phi2));
}
bool BbAna::pass_track_cut(int i)
{
  bool result = false;
  result = TrkPt->at(i)>10*GEV
    && TrkSctHits->at(i)+TrkSctDeadSensors->at(i)>5
    && TrkPixHits->at(i)+TrkPixDeadSensors->at(i)>1
    && TrkPixHoles->at(i)+TrkSctHoles->at(i)<3
    && ( (TrkExpBLayerHit->at(i)==0 || TrkBLayerHits->at(i)>0) );
  result = result && fabs(TrkEta->at(i))<2.4   && fabs(Trkz0->at(i)*sin(Trktheta->at(i)))<1.5;
  int ntrt = TrkTrtHits->at(i) + TrkTrtOutliers->at(i);
  if(fabs(TrkEta->at(i))<1.9)
    result = result && ntrt > 5 && TrkTrtOutliers->at(i)< 0.9*ntrt;
  else {
    if(ntrt>5) result = result && TrkTrtOutliers->at(i) < 0.9*ntrt;
  }

  return result;
}
bool BbAna::pass_muon_cut(int i)
{
  bool result = false;
  result = 
    ( (MuonExpBLayerHit->at(i)==0 || MuonBLayerHits->at(i)>0) )
    && MuonPixHits->at(i)+MuonPixDeadSensors->at(i)>1
    && MuonSctHits->at(i)+MuonSctDeadSensors->at(i)>5 
    && MuonPixHoles->at(i)+MuonSctHoles->at(i)<3
    && fabs(MuonEta->at(i))<2.4 && fabs(Muonz0->at(i)*sin(Muontheta->at(i)))<1.5;
  int ntrt = MuonTrtHits->at(i) + MuonTrtOutliers->at(i);
  if(fabs(MuonEta->at(i))<1.9)
    result = result && ntrt > 5 && MuonTrtOutliers->at(i)< 0.9*ntrt;
  else {
    if(ntrt>5) result = result && MuonTrtOutliers->at(i) < 0.9*ntrt;
  }

//   if(MuonPtms->at(i)<50*GEV){
//     TLorentzVector vms,vid;
//     vms.SetPtEtaPhiE(MuonPtmsextrap->at(i),MuonEta->at(i),MuonPhi->at(i),sqrt(pow(MuonPtmsextrap->at(i)*cosh(MuonEta->at(0)),2)+pow(105.7,2)));
//     vid.SetPtEtaPhiE(MuonPtid->at(i),MuonEta->at(i),MuonPhi->at(i),sqrt(pow(MuonPtid->at(i)*cosh(MuonEta->at(0)),2)+pow(105.7,2)));
//     double p_ms = vms.P();
//     double p_id = vid.P();
//     result = result && (p_ms - p_id)> (-0.4*p_id);
//   }

  return result;  
}
bool BbAna::pass_electron_cut(int i)
{
  bool result = false;

  double elet = 0.;
  if( ElectronSctHits->at(i)+ElectronPixHits->at(i) < 4 ) elet = ElectronClusEt->at(i);
  else elet = ElectronClusE->at(i)/cosh(ElectronEta->at(i));

  result = ElectronisTight->at(i)>9.;
  result = result && ElectronGoodOQ->at(i)>0;
  result = result && elet>20*GEV;
  result = result && ElectronSctHits->at(i)>5 && ElectronPixHits->at(i)>1
    && fabs(Electronz0->at(i)*sin(Electrontheta->at(i)))<2.;


  return result;
}
double BbAna::get_mu_trig_eff(float pt, float eta, float &err)
{
  //return eff for EF_mu20
  //see https://phyweb.lbl.gov:8090/Multiobject/17 for details
  double trigeff = 0.; err = 0.0;

  float p0 = 0.7512;  float ep0 = 0.00662;// These are numbers
  float p1 = 18.84;   float ep1 = 0.8053;// for barrel
  float p2 = 0.5598;  float ep2 = 5.212;//
  if(fabs(eta) > 1.05){
    p0 = 0.8864; ep0 = 0.0049;
    p1 = 19.21;  ep1 = 0.208;
    p2 = 2.08;   ep2 = 0.8225;
  }
  trigeff = p0 / (1 + pow(81,((p1-pt)/p2)) );
  double err0 = ep0 / (1 + pow(81,((p1-pt)/p2)) ); // d(trigeff)/dp0
  double err1 = -1 * p0 * pow((1 + pow(81,((p1-pt)/p2)) ),-2) * pow(81,((p1-pt)/p2)) * log(81) * ep1/p2; // d(trigeff)/dp1
  double err2 = p0 * pow((1 + pow(81,((p1-pt)/p2)) ),-2) * pow(81,((p1-pt)/p2)) * log(81) * (p1-pt) * pow(p2,-2) * ep2;// d(trigeff)/dp2
  err = sqrt(err0*err0 + err1*err1 + err2*err2);


  return trigeff;
}
int BbAna::dobbPlots()
{
  //if(NBadJets>0 && isData) return 0;
  if(nmusel<2) return 0;
  //At least two muons in the event.
  int mu1 = 0;
  int mu2 = 1;

  if(delta_r(MuonEta->at(musel[mu1].ind),MuonPhi->at(musel[mu1].ind),MuonEta->at(musel[mu2].ind),MuonPhi->at(musel[mu2].ind))<0.2) return 0;

  bool pass_TRIG = MuonMatchTrig->at(musel[mu1].ind)==4 ;
  bool hid01 = fabs(Muond0->at(musel[mu1].ind)/Muond0err->at(musel[mu1].ind)) > 3;
  bool hid02 = fabs(Muond0->at(musel[mu2].ind)/Muond0err->at(musel[mu2].ind)) > 3;
  double ptrat = (MuonPtid->at(musel[mu2].ind)-MuonPtms->at(musel[mu2].ind))/MuonPtid->at(musel[mu2].ind);
  bool iso1= pass_isolation(musel[mu1].ind,13,musel[mu1].mom.Pt());//Muonptcone20->at(musel[mu1].ind)/musel[mu1].mom.Pt()<0.2;
  bool iso2= pass_isolation(musel[mu2].ind,13,musel[mu2].mom.Pt());
  bool passPt= musel[mu1].mom.Pt()>25*GEV && musel[mu2].mom.Pt()>15*GEV;

  if(passPt) h.whatsign->Fill(musel[mu1].id*musel[mu2].id/169);
    
  // in the bb MC, plot the pt of muons if both from c, vs either from b.
  if(musel[mu1].id*musel[mu2].id != 169) return 0;

  if(passPt) h.mom12->Fill(abs(MuonMotherID->at(musel[mu1].ind)),abs(MuonMotherID->at(musel[mu2].ind)));

  if(!isData){
    if(fabs(MuonMotherID->at(musel[mu1].ind))>400 && fabs(MuonMotherID->at(musel[mu1].ind))<500
       && fabs(MuonMotherID->at(musel[mu2].ind))>400 && fabs(MuonMotherID->at(musel[mu2].ind))<500 ){//both from c
      h.muonpt_bbcc[0][0]->Fill(musel[mu1].mom.Pt()*MEV2GEV);
      h.muonpt_bbcc[0][1]->Fill(musel[mu2].mom.Pt()*MEV2GEV);
    }
    else if( ( fabs(MuonMotherID->at(musel[mu1].ind))>500 && fabs(MuonMotherID->at(musel[mu1].ind))<600 ) 
	     ||  ( fabs(MuonMotherID->at(musel[mu2].ind))>500 && fabs(MuonMotherID->at(musel[mu2].ind))<600 ) ){
      h.muonpt_bbcc[1][0]->Fill(musel[mu1].mom.Pt()*MEV2GEV);
      h.muonpt_bbcc[1][1]->Fill(musel[mu2].mom.Pt()*MEV2GEV);
    }
  }



  double thisEventsWeight = 1.0;
  if(!isData){
      //weight for the trigger efficiency and ID scalefactors.
      // ID SF
     double mu1id = StacoCBSCF.scaleFactor(musel[mu1].mom);
     double mu2id = StacoCBSCF.scaleFactor(musel[mu2].mom);
     // Trigger Eff
     float err[2] = {0.};
     double mu1trig = get_mu_trig_eff(musel[mu1].mom.Pt()*MEV2GEV,MuonEta->at(musel[mu1].ind),err[0]);
     double mu2trig = get_mu_trig_eff(musel[mu2].mom.Pt()*MEV2GEV,MuonEta->at(musel[mu2].ind),err[1]);
     double mu_trig_weight = 1 - (1-mu1trig)*(1-mu2trig);
     thisEventsWeight = mu_trig_weight*mu1id*mu2id;
  }

  if(isData){
    if(musel[mu2].mom.Pt()<20*GEV) {
      passPt = passPt && pass_TRIG;
    }
  }
 

  bool is_HF = (1 - MuonPtms->at(musel[mu1].ind)/MuonPtid->at(musel[mu1].ind))<0.2 
    && (1 - MuonPtms->at(musel[mu2].ind)/MuonPtid->at(musel[mu2].ind))<0.2;
  bool notSignal= NTRK_SEL < 10;
  //bool notSignal= true;
  float dphi= delta_phi(MuonPhi->at(musel[mu1].ind),MuonPhi->at(musel[mu2].ind));
  //if(passPt && notSignal ) h.bb_ptrat->Fill(ptrat,thisEventsWeight);
  if(passPt){
    h.d0isol->Fill(Muonptcone20->at(musel[mu1].ind)/musel[mu1].mom.Pt(),fabs(Muond0->at(musel[mu1].ind)/Muond0err->at(musel[mu1].ind)));
    int icr = -1;
    if(!hid01 && iso1)     icr=0;
    else if(hid01 && iso1) icr=1;
    else if(hid01 && !iso1 && is_HF)  icr=2; //region C with extra HF enhancing selection is_HF.
    else if(!hid01 && !iso1) icr=3;
//     if(iso1) icr=0;
//     else icr=1;
    if(icr==-1) return 0;
    if(icr==2 && notSignal) h.bb_ptrat->Fill(ptrat,thisEventsWeight);
    if(notSignal){
      h.bb_h[icr][0]->Fill(Muond0->at(musel[mu1].ind)/Muond0err->at(musel[mu1].ind),thisEventsWeight);
      h.bb_h[icr][1]->Fill(Muond0->at(musel[mu2].ind)/Muond0err->at(musel[mu2].ind),thisEventsWeight);
      h.bb_h[icr][2]->Fill(Muonptcone20->at(musel[mu1].ind)/musel[mu1].mom.Pt(),thisEventsWeight);
      h.bb_h[icr][3]->Fill(Muonptcone20->at(musel[mu2].ind)/musel[mu2].mom.Pt(),thisEventsWeight);
      h.bb_h[icr][4]->Fill(musel[mu1].mom.Pt()*MEV2GEV,thisEventsWeight);
      h.bb_h[icr][5]->Fill(musel[mu2].mom.Pt()*MEV2GEV,thisEventsWeight);
      h.bb_h[icr][6]->Fill(MissingET*MEV2GEV,thisEventsWeight);
      h.bb_h[icr][7]->Fill(dphi,thisEventsWeight);
      bb_nev[icr]++; bb_wgt[icr] += thisEventsWeight;
    }
    //h.bb_h[icr][8]->Fill(NTRK_SEL,thisEventsWeight);
    h.bb_h[icr][8]->Fill(NTRK_SEL);

  }
  
  return 1;
}

void BbAna::BookHistograms()
{
  h.whatsign = new TH1F("whatsign","Mu1Q*Mu2Q",5,-2.5,2.5);
  h.mumotherid = new TH1F("mumotherid","Muon Mother PdgId",6000,0,6000);
  h.mom12 = new TH2F("mom12","Mom1 vs Mom2",600,0,6000,600,0,6000);
  h.trackpt = new TH1F("trackpt","trackpt",200,0,200);
  //bb_h[4][9] {d0sig1,d0sig2,isol1,isol2,pt1,pt2,met,dphi,ntrk_sel}
  TString crname[4] = {"rA","rB","rC","rD"};
  TString plotname[9] = {"d0sig1","d0sig2","isol1","isol2","pt1","pt2","met","dphi","NTRK"};
  int nbins[9] = {200,200,200,200,200,200,200,32,50};
  float blo[9] = {-100,-100,0,0,0,0,0,0,0};
  float bhi[9] = {100,100,10,10,200,200,200,3.2,50};
  for(int icr=0; icr<4; icr++){
    for(int iplot=0; iplot<9; iplot++){
      TString name = crname[icr] + plotname[iplot];
      h.bb_h[icr][iplot] = new TH1F(name,name,nbins[iplot],blo[iplot],bhi[iplot]); h.bb_h[icr][iplot]->Sumw2();
    }
  }
  h.bb_ptrat    = new TH1F("bb_ptrat","mu2 ptrat",100,-5,5);
  h.d0isol      = new TH2F("bb_d0isol","d0sig vs Isol (#mu_{1})",200,0,10,40,0,20);

  h.muonpt_bbcc[0][0] = new TH1F("mupt_cc_0","mu1 pT, both c",100,0,100);
  h.muonpt_bbcc[0][1] = new TH1F("mupt_cc_1","mu2 pT, both c",100,0,100);
  h.muonpt_bbcc[1][0] = new TH1F("mupt_bc_0","mu1 pT, one b",100,0,100);
  h.muonpt_bbcc[1][1] = new TH1F("mupt_bc_1","mu2 pT, one b",100,0,100);
  
  //histograms for bbar measurement.
  /*
  h.bb_d0sig[0] = new TH1F("bb_d0sig1","mu1 d0sig",100,-50,50);
  h.bb_d0sig[1] = new TH1F("bb_d0sig2","mu2 d0sig",100,-50,50);
  h.bb_isol[0][0] = new TH1F("bb1_ptcone20","ptcone20 leading muon,high d0 #mu_{1}",100,0,100);
  h.bb_isol[0][1] = new TH1F("bb1_ptcone20pt","ptcone20/p_{T} leading muon,high d0 #mu_{1}",200,0,10);
  h.bb_isol[0][2] = new TH1F("bb1_ptcone20_2","ptcone20 second muon,high d0 #mu_{1}, noisol",100,0,100);
  h.bb_isol[0][3] = new TH1F("bb1_ptcone20pt_2","ptcone20/p_{T} second muon,high d0 #mu_{1}, noisol",200,0,10);
  h.bb_isol[1][0] = new TH1F("bb2_ptcone20","ptcone20 leading muon,high d0 #mu_{1,2}",100,0,100);
  h.bb_isol[1][1] = new TH1F("bb2_ptcone20pt","ptcone20/p_{T} leading muon,high d0 #mu_{1,2}",200,0,10);
  h.bb_isol[1][2] = new TH1F("bb2_ptcone20_2","ptcone20 second muon,high d0 #mu_{1,2}",100,0,100);
  h.bb_isol[1][3] = new TH1F("bb2_ptcone20pt_2","ptcone20/p_{T} second muon,high d0 #mu_{1,2}",200,0,10);
  h.bb_ntrk[0]    = new TH1F("bb_ntrk_niso","NTracks, non-isolated leading muon",50,0,50);
  h.bb_ntrk[1]    = new TH1F("bb_ntrk_nisod0","NTracks, non-isolated leading muon, high d0 #mu_{1}",50,0,50);
  h.bb_ntrk[2]    = new TH1F("bb_ntrk","NTracks, isolated leading muon",50,0,50);
  */

}







// ======================================================================================================================================
int EtaPhiBinning::getSector(double phi) const {
  int sector = -9;
  if (phi>2.905 || phi<=-2.905) sector = 9;
  else if (phi>2.59 && phi<=2.905) sector = 8;
  else if (phi>2.12 && phi<=2.59) sector = 7;
  else if (phi>1.805 && phi<=2.12) sector = 6;
  else if (phi>1.335 && phi<=1.805) sector = 5;
  else if (phi>1.02 && phi<=1.335) sector = 4;
  else if (phi>0.55 && phi<=1.02) sector = 3;
  else if (phi>0.235 && phi<=0.55) sector = 2;
  else if (phi>-0.235 && phi<=0.235) sector = 1;
  else if (phi>-0.55 && phi<=-0.235) sector = 16;
  else if (phi>-1.02 && phi<=-0.55) sector = 15;
  else if (phi>-1.335 && phi<=-1.02) sector = 14;
  else if (phi>-1.805 && phi<=-1.335) sector = 13;
  else if (phi>-2.12 && phi<=-1.805) sector = 12;
  else if (phi>-2.59 && phi<=-2.12) sector = 11;
  else if (phi>-2.905 && phi<=-2.59) sector = 10;
  return sector;
}
int EtaPhiBinning::getECSector(double phi) const {
  int sector = -9;
  if (phi>3.011 && phi<=-3.011) sector = 9;
  else if (phi>2.487 && phi<=3.011) sector = 8;
  else if (phi>2.225 && phi<=2.487) sector = 7;
  else if (phi>1.702 && phi<=2.225) sector = 6;
  else if (phi>1.440 && phi<=1.702) sector = 5;
  else if (phi>0.916 && phi<=1.440) sector = 4;
  else if (phi>0.655 && phi<=0.916) sector = 3;
  else if (phi>0.131 && phi<=0.655) sector = 2;
  else if (phi>-0.131 && phi<=0.131) sector = 1;
  else if (phi>-0.655 && phi<=-0.131) sector = 16;
  else if (phi>-0.916 && phi<=-0.655) sector = 15;
  else if (phi>-1.440 && phi<=-0.916) sector = 14;
  else if (phi>-1.702 && phi<=-1.440) sector = 13;
  else if (phi>-2.225 && phi<=-1.702) sector = 12;
  else if (phi>-2.487 && phi<=-2.225) sector = 11;
  else if (phi>-3.011 && phi<=-2.487) sector = 10;
  return sector;
}
int EtaPhiBinning::getCoarseNSector(double phi) const {
    if ( (fabs(phi)>=0.18 && fabs(phi)<0.285)
            || (fabs(phi)>=0.5 && fabs(phi)<0.605)
            || (fabs(phi)>=0.965 && fabs(phi)<1.07)
            || (fabs(phi)>=1.285 && fabs(phi)<1.39)
            || (fabs(phi)>=1.75 && fabs(phi)<1.855)
            || (fabs(phi)>=2.07 && fabs(phi)<2.175)
            || (fabs(phi)>=2.535 && fabs(phi)<2.64)
            || (fabs(phi)>=2.855 && fabs(phi)<2.96) )
        return 2;
    else return 1;
}
int EtaPhiBinning::symmetricBin(const TLorentzVector* mst) const
{
  //Region bin based on eta and phi
  double mu_phi = mst->Phi();
  double mu_eta = mst->Eta();
  int mu_sector = getSector(mu_phi);
  bool isSmSect = mu_sector%2==0 ? true : false;
  int mu_nsect = getCoarseNSector(mu_phi);
  int mu_ec_sector = getECSector(mu_phi);
  bool isSmECSect = mu_ec_sector%2==0 ? true : false;
  
  //Feet
  if ( fabs(mu_eta)<0.97
       && ((mu_phi<-1.01 && mu_phi>-1.36)
	   || (mu_phi<-1.79 && mu_phi>-2.14)) )
    return binFEET;
  if ( ( fabs(mu_eta)>0.51 && fabs(mu_eta)<0.81)
       && (mu_phi<-1.36 && mu_phi>-1.79) )
    return binFEET;
  //To put tiny regions of overlap into feet 
  //if ( mu_sector==13 && mst->barrelSectors()>1) return binFEET;

  //Non-Feet Barrel
  if ( fabs(mu_eta)<0.97 ) {
    if (mu_nsect==2) return bin2BARREL;
    else if (mu_nsect==1) {
      if (isSmSect) return bin1BARRELSM;
      else return bin1BARRELLG;
    }
  }
  
  //Transition and BIS78
  if ( fabs(mu_eta)>=0.97 && fabs(mu_eta)<1.11 )
    return binTRANSITION;
  
  else if ( fabs(mu_eta)>=1.11 && fabs(mu_eta)<1.19 ) {
    if ( isSmSect || mu_nsect==2 ) return binTRANSITION;
    else {
      if (isSmECSect) return binENDCAPSM;
      else return binENDCAPLG;
    }
    
  } else if ( fabs(mu_eta)>=1.19 && fabs(mu_eta)<1.25 ) {
    if ( isSmSect && mu_nsect==1 ) return binTRANSITION;
    else if ( isSmECSect ) return binENDCAPSM;
    else return binENDCAPLG;
    
  }
  
  //BEE
  if ( fabs(mu_eta)>=1.42 && fabs(mu_eta)<1.72
       && mu_nsect==1 && isSmSect )
    return binBEE;
  
  //Endcap
  if ( fabs(mu_eta)>=1.25 && fabs(mu_eta)<1.97 ) {
    if (isSmECSect) return binENDCAPSM;
    else return binENDCAPLG;
  }
  
  //Forward
  if ( fabs(mu_eta)>=1.97 && fabs(mu_eta)<2.01 ) {
    if ( mu_nsect==1 && !isSmSect ) {
      if (isSmECSect) return binENDCAPSM;
      else return binENDCAPLG;
    } else return binFORWARDSM;
    }
  
    if ( fabs(mu_eta)>=2.01 ) {
    if ( mu_nsect>1 || isSmSect ) return binFORWARDSM;
    else {
    if (isSmECSect) return binFORWARDSM;
    else return binFORWARDLG;
    }
    }
    
    //control flow should now be here!!
    //std::cout<<"Control flow should not be here: MuSelectionToolsLine="<<__LINE__<<std::endl;
    return binUNKNOWN;
    
    }

// 
// MuonEfficiencyCorrections-00-03-08  
// Dated:June 23, 2011
// 
StacoCBScaleEffFactors::StacoCBScaleEffFactors(void) {

    m_last_run_periodB = 178109;

////////////////////////////////////////////////////////////////
// FILL THE SCALE FACTOR AND SCALE FACTOR UNCERTAINTY VECTORS //
////////////////////////////////////////////////////////////////

    m_scale_factor_A = std::vector<double>(10);
    m_scale_factor_uncertainty_A = std::vector<double>(10);
    m_scale_factor_C_high_pt = std::vector<double>(10);
    m_scale_factor_uncertainty_C_high_pt = std::vector<double>(10);
    m_scale_factor_C_low_pt = std::vector<double>(10);
    m_scale_factor_uncertainty_C_low_pt = std::vector<double>(10);

    m_scale_factor_A[0] = 0.99679;
    m_scale_factor_A[1] = 1.00910;
    m_scale_factor_A[2] = 0.99299;
    m_scale_factor_A[3] = 0.96569;
    m_scale_factor_A[4] = 0.91706;
    m_scale_factor_A[5] = 0.99053;
    m_scale_factor_A[6] = 0.99441;
    m_scale_factor_A[7] = 0.98640;
    m_scale_factor_A[8] = 1.00846;
    m_scale_factor_A[9] = 0.99240;

    m_scale_factor_uncertainty_A[0] = 0.00325;
    m_scale_factor_uncertainty_A[1] = 0.00302;
    m_scale_factor_uncertainty_A[2] = 0.00379;
    m_scale_factor_uncertainty_A[3] = 0.00652;
    m_scale_factor_uncertainty_A[4] = 0.00646;
    m_scale_factor_uncertainty_A[5] = 0.00317;
    m_scale_factor_uncertainty_A[6] = 0.00306;
    m_scale_factor_uncertainty_A[7] = 0.00506;
    m_scale_factor_uncertainty_A[8] = 0.00490;
    m_scale_factor_uncertainty_A[9] = 0.00340;

    m_scale_factor_C_high_pt[0] = 0.99385;
    m_scale_factor_C_high_pt[1] = 1.00289;
    m_scale_factor_C_high_pt[2] = 0.99826;
    m_scale_factor_C_high_pt[3] = 0.96668;
    m_scale_factor_C_high_pt[4] = 0.90399;
    m_scale_factor_C_high_pt[5] = 0.98067;
    m_scale_factor_C_high_pt[6] = 0.99212;
    m_scale_factor_C_high_pt[7] = 0.95546;
    m_scale_factor_C_high_pt[8] = 0.99287;
    m_scale_factor_C_high_pt[9] = 0.98921;

    m_scale_factor_uncertainty_C_high_pt[0] = 0.00457;
    m_scale_factor_uncertainty_C_high_pt[1] = 0.00414;
    m_scale_factor_uncertainty_C_high_pt[2] = 0.00512;
    m_scale_factor_uncertainty_C_high_pt[3] = 0.00957;
    m_scale_factor_uncertainty_C_high_pt[4] = 0.00880;
    m_scale_factor_uncertainty_C_high_pt[5] = 0.00420;
    m_scale_factor_uncertainty_C_high_pt[6] = 0.00397;
    m_scale_factor_uncertainty_C_high_pt[7] = 0.00858;
    m_scale_factor_uncertainty_C_high_pt[8] = 0.00681;
    m_scale_factor_uncertainty_C_high_pt[9] = 0.00473;

    m_scale_factor_C_low_pt[0] = 0.98571;
    m_scale_factor_C_low_pt[1] = 1.00240;
    m_scale_factor_C_low_pt[2] = 1.00070;
    m_scale_factor_C_low_pt[3] = 0.98692;
    m_scale_factor_C_low_pt[4] = 0.89969;
    m_scale_factor_C_low_pt[5] = 0.99031;
    m_scale_factor_C_low_pt[6] = 0.99262;
    m_scale_factor_C_low_pt[7] = 0.96997;
    m_scale_factor_C_low_pt[8] = 0.99599;
    m_scale_factor_C_low_pt[9] = 0.99305;

    m_scale_factor_uncertainty_C_low_pt[0] = 0.00501;
    m_scale_factor_uncertainty_C_low_pt[1] = 0.00443;
    m_scale_factor_uncertainty_C_low_pt[2] = 0.00527;
    m_scale_factor_uncertainty_C_low_pt[3] = 0.00943;
    m_scale_factor_uncertainty_C_low_pt[4] = 0.00990;
    m_scale_factor_uncertainty_C_low_pt[5] = 0.00465;
    m_scale_factor_uncertainty_C_low_pt[6] = 0.00453;
    m_scale_factor_uncertainty_C_low_pt[7] = 0.00910;
    m_scale_factor_uncertainty_C_low_pt[8] = 0.00667;
    m_scale_factor_uncertainty_C_low_pt[9] = 0.00448;

// separate values for period B //
    m_scale_factor_A_periodB = std::vector<double>(10);
    m_scale_factor_uncertainty_A_periodB = std::vector<double>(10);
    m_scale_factor_C_periodB = std::vector<double>(10);
    m_scale_factor_uncertainty_C_periodB = std::vector<double>(10);

    m_scale_factor_A_periodB[0] = 0.9837;
    m_scale_factor_A_periodB[1] = 0.9879;
    m_scale_factor_A_periodB[2] = 0.9652;
    m_scale_factor_A_periodB[3] = 0.9179;
    m_scale_factor_A_periodB[4] = 0.9343;
    m_scale_factor_A_periodB[5] = 0.9286;
    m_scale_factor_A_periodB[6] = 0.9824;
    m_scale_factor_A_periodB[7] = 0.9720;
    m_scale_factor_A_periodB[8] = 0.9809;
    m_scale_factor_A_periodB[9] = 1.0080;

    m_scale_factor_uncertainty_A_periodB[0] = 0.0189;
    m_scale_factor_uncertainty_A_periodB[1] = 0.0202;
    m_scale_factor_uncertainty_A_periodB[2] = 0.0221;
    m_scale_factor_uncertainty_A_periodB[3] = 0.0406;
    m_scale_factor_uncertainty_A_periodB[4] = 0.0364;
    m_scale_factor_uncertainty_A_periodB[5] = 0.0249;
    m_scale_factor_uncertainty_A_periodB[6] = 0.0176;
    m_scale_factor_uncertainty_A_periodB[7] = 0.0367;
    m_scale_factor_uncertainty_A_periodB[8] = 0.0351;
    m_scale_factor_uncertainty_A_periodB[9] = 0.0184;

    m_scale_factor_C_periodB[0] = 0.9768;
    m_scale_factor_C_periodB[1] = 0.9583;
    m_scale_factor_C_periodB[2] = 0.9725;
    m_scale_factor_C_periodB[3] = 0.9861;
    m_scale_factor_C_periodB[4] = 0.9085;
    m_scale_factor_C_periodB[5] = 0.9906;
    m_scale_factor_C_periodB[6] = 1.0093;
    m_scale_factor_C_periodB[7] = 0.9339;
    m_scale_factor_C_periodB[8] = 0.9905;
    m_scale_factor_C_periodB[9] = 0.9767;

    m_scale_factor_uncertainty_C_periodB[0] = 0.0214;
    m_scale_factor_uncertainty_C_periodB[1] = 0.0210;
    m_scale_factor_uncertainty_C_periodB[2] = 0.0233;
    m_scale_factor_uncertainty_C_periodB[3] = 0.0388;
    m_scale_factor_uncertainty_C_periodB[4] = 0.0367;
    m_scale_factor_uncertainty_C_periodB[5] = 0.0166;
    m_scale_factor_uncertainty_C_periodB[6] = 0.0137;
    m_scale_factor_uncertainty_C_periodB[7] = 0.0550;
    m_scale_factor_uncertainty_C_periodB[8] = 0.0298;
    m_scale_factor_uncertainty_C_periodB[9] = 0.0201;

}

//*****************************************************************************

////////////////////////
// METHOD scaleFactor //
////////////////////////

double StacoCBScaleEffFactors::scaleFactor(
                                            const TLorentzVector & tlv) const {

    return scaleFactor(tlv, m_last_run_periodB+1);

}

//*****************************************************************************

///////////////////////////////////
// METHOD scaleFactorUncertainty //
///////////////////////////////////

double StacoCBScaleEffFactors::scaleFactorUncertainty(
                                            const TLorentzVector & tlv) const {

    return scaleFactorUncertainty(tlv, m_last_run_periodB+1);

}

//*****************************************************************************

//////////////////////////////
// METHOD scaleFactor(., .) //
//////////////////////////////

double StacoCBScaleEffFactors::scaleFactor(
                                        const TLorentzVector & tlv,
                                        const unsigned int & run_nb) const {

    int bin(m_eta_phi_binning.symmetricBin(&tlv)-1);
    if(bin<0) return 1.; //UNKNOWN BIN
    if (tlv.Eta()>0) {
        if (run_nb>m_last_run_periodB) {
            return m_scale_factor_A[bin];
        } else {
            return m_scale_factor_A_periodB[bin];
        }
    }
    if (run_nb>m_last_run_periodB) {
        if (tlv.Pt()>40000) {
            return m_scale_factor_C_high_pt[bin];
        } else {
            return m_scale_factor_C_low_pt[bin];   
        }
    } else {
        return m_scale_factor_C_periodB[bin];
    }

}

//*****************************************************************************

/////////////////////////////////////////
// METHOD scaleFactorUncertainty(., .) //
/////////////////////////////////////////

double StacoCBScaleEffFactors::scaleFactorUncertainty(
                                            const TLorentzVector & tlv,
                                            const unsigned int & run_nb) const {

    int bin(m_eta_phi_binning.symmetricBin(&tlv)-1);
    if(bin<0) return 1.; //UNKNOWN BIN
    if (tlv.Eta()>0) {
        if (run_nb>m_last_run_periodB) {
            return m_scale_factor_uncertainty_A[bin];
        } else {
            return m_scale_factor_uncertainty_A_periodB[bin];
        }
    }
    if (run_nb>m_last_run_periodB) {
        if (tlv.Pt()>40000) {
            return m_scale_factor_uncertainty_C_high_pt[bin];
        } else {
            return m_scale_factor_uncertainty_C_low_pt[bin];
        }
    } else {
        return m_scale_factor_uncertainty_C_periodB[bin];
    }

}

//*****************************************************************************

/////////////////////////////////////////////
// METHOD scaleFactorSystematicUncertainty //
/////////////////////////////////////////////

double StacoCBScaleEffFactors::scaleFactorSystematicUncertainty(
                                            const TLorentzVector & tlv) const {

    return 0.002+4.2e-6*tlv.E()/1000.0;

}


