void disp_ml::BookHistograms(){

 h.nevt = new TH1F("nEvents", "0-nEvtTotal, 1-nEvtGood, 2-nEvtTrigger, 3-nEvtPass",5,0,5);

  h.charge[0] = new TH1F("muon_charge", "Muon charge", 10, -5, 5);
  h.charge[1] = new TH1F("electron_charge", "Electron charge", 10, -5, 5);

  //double xbins1[7] = {-50, -10, -5, 0, 5, 10, 50};
  //h.hnew1 = (TH1F*)h.dxy[0]->Rebin(4,"hnew1",xbins1);
  // h.dxy[0]   = new TH1F("mu_dxy", "mu_dxy", 6, xbins1);

  h.dxy[0]   = new TH1F("mu_dxy", "mu_dxy", 1000, -50, 50);
  h.dz[0]    = new TH1F("mu_dz", "mu_dz", 1000, -50, 50);
  h.dxy[1]   = new TH1F("mu_|dxy|", "mu_|dxy|", 5000, 0, 50);
  h.dz[1]    = new TH1F("mu_|dz|", "mu_|dz|", 1000, 0, 100);
  h.ip3d[0]  = new TH1F("mu_ip3d", "mu_ip3d", 1000, 0, 100);
  h.sip3d[0] = new TH1F("mu_sip3d", "mu_sip3d", 1000, 0, 100);

  h.dxy[2]   = new TH1F("el_dxy", "el_dxy", 1000, -50, 50);
  h.dz[2]    = new TH1F("el_dz", "el_dz", 1000, -50, 50);
  h.dxy[3]   = new TH1F("el_|dxy|", "el_|dxy|", 5000, 0, 50);
  h.dz[3]    = new TH1F("el_|dz|", "el_|dz|", 1000, 0, 100);
  h.ip3d[1]  = new TH1F("el_ip3d", "el_ip3d", 1000, 0, 100);
  h.sip3d[1] = new TH1F("el_sip3d", "el_sip3d", 1000, 0, 100);

  h.dxy[4]   = new TH1F("lep_dxy", "lep_dxy", 1000, -50, 50);
  h.dz[4]    = new TH1F("lep_dz", "lep_dz", 1000, -50, 50);
  h.dxy[5]   = new TH1F("lep_|dxy|", "lep_|dxy|", 5000, 0, 50);
  h.dz[5]    = new TH1F("lep_|dz|", "lep_|dz|", 1000, 0, 100);
  h.ip3d[2]  = new TH1F("lep_ip3d", "lep_ip3d", 1000, 0, 100);
  h.sip3d[2] = new TH1F("lep_sip3d", "lep_sip3d", 1000, 0, 100);

  /*
    h.evtweight[0][0] = new TH1F("2l1d_sf", "2l1d_sf", 50, 0, 5);
    h.evtweight[0][1] = new TH1F("2l1d_trigeff", "2l1d_trigeff", 10, 0, 1);
    h.evtweight[0][2] = new TH1F("2l1d_evtwt", "2l1d_evtwt", 50, 0, 5);
    h.evtweight[1][0] = new TH1F("1l2d_sf", "1l2d_sf", 100, 0, 100);
    h.evtweight[1][1] = new TH1F("1l2d_trigeff", "1l2d_trigeff", 10, 0, 1);
    h.evtweight[1][2] = new TH1F("1l2d_evtwt", "1l2d_evtwt", 100, 0, 100);
    h.evtweight[2][0] = new TH1F("3d_sf", "3d_sf", 100, 0, 100);
    h.evtweight[2][1] = new TH1F("3d_trigeff", "3d_trigeff", 10, 0, 1);
    h.evtweight[2][2] = new TH1F("3d_evtwt", "3d_evtwt", 100, 0, 100);
  */

  
  h._2LonZ[0] = new TH1F("2LonZ_invmass", "2LonZ_invmass", 200, 0, 200);
  h._2LonZ[1] = new TH1F("2LonZ_met", "2LonZ_met", 200, 0, 200);
    
  h._2LSS[0] = new TH1F("2LSS_flavor", "0: all 2LSS, 1: mumu, 2: ee, 3: emu, 4: mue", 5, 0, 5);
  h._2LSS[1] = new TH1F("2LSS_invmass_ll", "2LSS M_{ll}", 200, 0, 200);
  h._2LSS[2] = new TH1F("2LSS_met", "2LSS_met", 200, 0, 200);
  h._2LSS[3] = new TH1F("mumu_invmass_ll", "mumu M_{#mu#mu}", 200, 0, 200);
  h._2LSS[4] = new TH1F("ee_invmass_ll", "ee M_{ee}", 200, 0, 200);
  h._2LSS[5] = new TH1F("emu_or_mue_invmass_ll", "M_{e#mu or #mue}", 200, 0, 200);
  h._2LSS[6] = new TH1F("2LSS_ht", "2LSS_ht", 500, 0, 500);
 
  h._3L[0]  = new TH1F("3L_invmass_3l", "3L_invmass_3l", 500, 0, 500);
  h._3L[1]  = new TH1F("3L_invmass_l0l1", "3L_invmass_l0l1", 200, 0, 200);
  h._3L[2]  = new TH1F("3L_invmass_l1l2", "3L_invmass_l1l2", 200, 0, 200);
  h._3L[3]  = new TH1F("3L_invmass_l2l0", "3L_invmass_l2l0", 200, 0, 200);
  h._3L[4]  = new TH1F("3L_met", "3L_met", 200, 0, 200);
  h._3L[5]  = new TH1F("3L_pt0", "3L_pt0", 200, 0, 200);
  h._3L[6]  = new TH1F("3L_pt1", "3L_pt1", 200, 0, 200);
  h._3L[7]  = new TH1F("3L_pt2", "3L_pt2", 200, 0, 200);
  h._3L[8]  = new TH1F("3L_lt", "3L_lt", 500, 0, 500);
  h._3L[9]  = new TH1F("3L_njet", "3L_njet", 10, 0, 10);
  h._3L[10] = new TH1F("3L_ht", "3L_ht", 500, 0, 500);
  h._3L[11] = new TH1F("3L_st", "3L_st", 500, 0, 500);

  h.cr_ttbar[0]  = new TH1F("cr_ttbar_met", "met", 200, 0, 200);
  h.cr_ttbar[1]  = new TH1F("cr_ttbar_pt0", "pt0", 200, 0, 200);
  h.cr_ttbar[2]  = new TH1F("cr_ttbar_pt1", "pt1", 200, 0, 200);
  h.cr_ttbar[3]  = new TH1F("cr_ttbar_lt", "lt", 500, 0, 500);
  h.cr_ttbar[4]  = new TH1F("cr_ttbar_njet", "njet", 10, 0, 10);
  h.cr_ttbar[5]  = new TH1F("cr_ttbar_bjet", "bjet", 10, 0, 10);
  h.cr_ttbar[6]  = new TH1F("cr_ttbar_ht", "ht", 500, 0, 500);

  double xbins1[6] = {0, 50, 150, 250, 350, 500};
  h.cr_ttbar_ht_rebinned = new TH1F("cr_ttbar_ht_rebinned", "cr_ttbar_ht_rebinned", 5, xbins1);
  
  h.cr_ttbar[7]  = new TH1F("cr_ttbar_mt0", "mt0", 200, 0, 200);
  h.cr_ttbar[8]  = new TH1F("cr_ttbar_mt1", "mt1", 200, 0, 200);
  h.cr_ttbar[9]  = new TH1F("cr_ttbar_l0iso", "l0iso", 150, 0, 0.15);
  h.cr_ttbar[10] = new TH1F("cr_ttbar_l1iso", "l1iso", 150, 0, 0.15);
  h.cr_ttbar[11] = new TH1F("cr_ttbar_dphil0l1", "#Delta#phi_{ll}", 32, 0, 3.2);
  h.cr_ttbar[12] = new TH1F("cr_ttbar_invmassl0l1", "M_{ll}", 200, 0, 200);
  h.cr_ttbar[13] = new TH1F("cr_ttbar_invmassj0j1", "M_{jj}", 200, 0, 200);
  h.cr_ttbar[14] = new TH1F("cr_ttbar_j0pt", "j0 pt", 200, 0, 200);
  h.cr_ttbar[15] = new TH1F("cr_ttbar_j1pt", "j1 pt", 200, 0, 200);

  h.vr_ttbar[0]  = new TH1F("vr_ttbar_met", "met", 200, 0, 200);
  h.vr_ttbar[1]  = new TH1F("vr_ttbar_pt0", "pt0", 200, 0, 200);
  h.vr_ttbar[2]  = new TH1F("vr_ttbar_pt1", "pt1", 200, 0, 200);
  h.vr_ttbar[3]  = new TH1F("vr_ttbar_lt", "lt", 500, 0, 500);
  h.vr_ttbar[4]  = new TH1F("vr_ttbar_njet", "njet", 10, 0, 10);
  h.vr_ttbar[5]  = new TH1F("vr_ttbar_bjet", "bjet", 10, 0, 10);
  h.vr_ttbar[6]  = new TH1F("vr_ttbar_ht", "ht", 500, 0, 500);
  h.vr_ttbar[7]  = new TH1F("vr_ttbar_mt0", "mt0", 200, 0, 200);
  h.vr_ttbar[8]  = new TH1F("vr_ttbar_mt1", "mt1", 200, 0, 200);
  h.vr_ttbar[9]  = new TH1F("vr_ttbar_l0iso", "l0iso", 150, 0, 0.15);
  h.vr_ttbar[10] = new TH1F("vr_ttbar_l1iso", "l1iso", 150, 0, 0.15);
  h.vr_ttbar[11] = new TH1F("vr_ttbar_dphil0l1", "#Delta#phi_{ll}", 32, 0, 3.2);
  h.vr_ttbar[12] = new TH1F("vr_ttbar_invmassl0l1", "M_{ll}", 200, 0, 200);
  h.vr_ttbar[13] = new TH1F("vr_ttbar_invmassj0j1", "M_{jj}", 200, 0, 200);
  h.vr_ttbar[14] = new TH1F("vr_ttbar_j0pt", "j0 pt", 200, 0, 200);
  h.vr_ttbar[15] = new TH1F("vr_ttbar_j1pt", "j1 pt", 200, 0, 200);

  h.cr_ttbar_2l1d[0]  = new TH1F("cr_ttbar_2l1d_met", "met", 200, 0, 200);
  h.cr_ttbar_2l1d[1]  = new TH1F("cr_ttbar_2l1d_pt0", "pt0", 200, 0, 200);
  h.cr_ttbar_2l1d[2]  = new TH1F("cr_ttbar_2l1d_pt1", "pt1", 200, 0, 200);
  h.cr_ttbar_2l1d[3]  = new TH1F("cr_ttbar_2l1d_lt", "lt", 500, 0, 500);
  h.cr_ttbar_2l1d[4]  = new TH1F("cr_ttbar_2l1d_njet", "njet", 10, 0, 10);
  h.cr_ttbar_2l1d[5]  = new TH1F("cr_ttbar_2l1d_bjet", "bjet", 10, 0, 10);
  h.cr_ttbar_2l1d[6]  = new TH1F("cr_ttbar_2l1d_ht", "ht", 500, 0, 500);

  double xbins2[6] = {0, 50, 150, 250, 350, 500};
  h.cr_ttbar_2l1d_ht_rebinned = new TH1F("cr_ttbar_2l1d_ht_rebinned", "cr_ttbar_2l1d_ht_rebinned", 5, xbins2);

  h.cr_ttbar_2l1d[7]  = new TH1F("cr_ttbar_2l1d_mt0", "mt0", 200, 0, 200);
  h.cr_ttbar_2l1d[8]  = new TH1F("cr_ttbar_2l1d_mt1", "mt1", 200, 0, 200);
  h.cr_ttbar_2l1d[9]  = new TH1F("cr_ttbar_2l1d_l0iso", "l0iso", 150, 0, 0.15);
  h.cr_ttbar_2l1d[10] = new TH1F("cr_ttbar_2l1d_l1iso", "l1iso", 150, 0, 0.15);
  h.cr_ttbar_2l1d[11] = new TH1F("cr_ttbar_2l1d_dphil0l1", "#Delta#phi_{ll}", 32, 0, 3.2);
  h.cr_ttbar_2l1d[12] = new TH1F("cr_ttbar_2l1d_invmassl0l1", "M_{ll}", 200, 0, 200);
  h.cr_ttbar_2l1d[13] = new TH1F("cr_ttbar_2l1d_invmassj0j1", "M_{jj}", 200, 0, 200);
  h.cr_ttbar_2l1d[14] = new TH1F("cr_ttbar_2l1d_j0pt", "j0 pt", 200, 0, 200);
  h.cr_ttbar_2l1d[15] = new TH1F("cr_ttbar_2l1d_j1pt", "j1 pt", 200, 0, 200);

  h.vr_ttbar_2l1d[0]  = new TH1F("vr_ttbar_2l1d_met", "met", 200, 0, 200);
  h.vr_ttbar_2l1d[1]  = new TH1F("vr_ttbar_2l1d_pt0", "pt0", 200, 0, 200);
  h.vr_ttbar_2l1d[2]  = new TH1F("vr_ttbar_2l1d_pt1", "pt1", 200, 0, 200);
  h.vr_ttbar_2l1d[3]  = new TH1F("vr_ttbar_2l1d_lt", "lt", 500, 0, 500);
  h.vr_ttbar_2l1d[4]  = new TH1F("vr_ttbar_2l1d_njet", "njet", 10, 0, 10);
  h.vr_ttbar_2l1d[5]  = new TH1F("vr_ttbar_2l1d_bjet", "bjet", 10, 0, 10);
  h.vr_ttbar_2l1d[6]  = new TH1F("vr_ttbar_2l1d_ht", "ht", 500, 0, 500);
  h.vr_ttbar_2l1d[7]  = new TH1F("vr_ttbar_2l1d_mt0", "mt0", 200, 0, 200);
  h.vr_ttbar_2l1d[8]  = new TH1F("vr_ttbar_2l1d_mt1", "mt1", 200, 0, 200);
  h.vr_ttbar_2l1d[9]  = new TH1F("vr_ttbar_2l1d_l0iso", "l0iso", 150, 0, 0.15);
  h.vr_ttbar_2l1d[10] = new TH1F("vr_ttbar_2l1d_l1iso", "l1iso", 150, 0, 0.15);
  h.vr_ttbar_2l1d[11] = new TH1F("vr_ttbar_2l1d_dphil0l1", "#Delta#phi_{ll}", 32, 0, 3.2);
  h.vr_ttbar_2l1d[12] = new TH1F("vr_ttbar_2l1d_invmassl0l1", "M_{ll}", 200, 0, 200);
  h.vr_ttbar_2l1d[13] = new TH1F("vr_ttbar_2l1d_invmassj0j1", "M_{jj}", 200, 0, 200);
  h.vr_ttbar_2l1d[14] = new TH1F("vr_ttbar_2l1d_j0pt", "j0 pt", 200, 0, 200);
  h.vr_ttbar_2l1d[15] = new TH1F("vr_ttbar_2l1d _j1pt", "j1 pt", 200, 0, 200);

  h.cr_wjets[0]  = new TH1F("cr_wjets_met", "met", 200, 0, 200);
  h.cr_wjets[1]  = new TH1F("cr_wjets_pt0", "pt0", 200, 0, 200);
  h.cr_wjets[2]  = new TH1F("cr_wjets_lt", "lt", 500 , 0, 500);
  h.cr_wjets[3]  = new TH1F("cr_wjets_njets", "njets", 10 , 0, 10);
  h.cr_wjets[4]  = new TH1F("cr_wjets_ht", "ht", 500 , 0, 500);
  h.cr_wjets[5]  = new TH1F("cr_wjets_mt0", "mt0", 200 , 0, 200);
  h.cr_wjets[6]  = new TH1F("cr_wjets_l0iso", "l0iso", 200 , 0, 200);
  h.cr_wjets[7]  = new TH1F("cr_wjets_dphil0j0", "#Delta#Phi_{l0j0}", 200 , 0, 200);
  h.cr_wjets[8]  = new TH1F("cr_wjets_j0pt", "j0 pt", 200 , 0, 200);

  h.cr_wjets_1l2d[0]  = new TH1F("cr_wjets_1l2d_met", "met", 200, 0, 200);
  h.cr_wjets_1l2d[1]  = new TH1F("cr_wjets_1l2d_pt0", "pt0", 200, 0, 200);
  h.cr_wjets_1l2d[2]  = new TH1F("cr_wjets_1l2d_lt", "lt", 500 , 0, 500);
  h.cr_wjets_1l2d[3]  = new TH1F("cr_wjets_1l2d_njets", "njets", 10 , 0, 10);
  h.cr_wjets_1l2d[4]  = new TH1F("cr_wjets_1l2d_ht", "ht", 500 , 0, 500);
  h.cr_wjets_1l2d[5]  = new TH1F("cr_wjets_1l2d_mt0", "mt0", 200 , 0, 200);
  h.cr_wjets_1l2d[6]  = new TH1F("cr_wjets_1l2d_l0iso", "l0iso", 200 , 0, 200);
  h.cr_wjets_1l2d[7]  = new TH1F("cr_wjets_1l2d_dphil0j0", "#Delta#Phi_{l0j0}", 200 , 0, 200);
  h.cr_wjets_1l2d[8]  = new TH1F("cr_wjets_1l2d_j0pt", "j0 pt", 200 , 0, 200);
 

  h.nevsel = new TH1F("nEvSel", "1: 2l1d, 2: 1l2d, 3: 3d", 5,0,5);
  
  TString evsel_name[3] = {"2l1d_", "1l2d_", "3d_"};
  TString plotname[45] = {"met","pt_3l","imass_3l","pt0","pt1","pt2","pt_l0l1","delR_l0l1","delPhi_l0l1","delPhi_l0met","imass_l0l1","mt0","pt_l1l2","delR_l1l2","delPhi_l1l2","delPhi_l1met","imass_l1l2","mt1","pt_l2l0","delR_l2l0","delPhi_l2l0","delPhi_l2met","imass_l2l0","mt2","HT","njet","dRmin_l0j","dRmin_l1j","dRmin_l2j","l0_dxy","l0_dz","l0_ip3d","l0_sip3d","l0_reliso03","l1_dxy","l1_dz","l1_ip3d","l1_sip3d","l1_reliso03","l2_dxy","l2_dz","l2_ip3d","l2_sip3d","l2_reliso03","bjets"};
  int nbins[45] = {200,500,500,200,200,200,500,100,32,32,500,200,500,100,32,32,500,200,500,100,32,32,500,200,200,10,100,100,100,2000,2000,200,500,1500,2000,2000,200,500,1500,2000,2000,200,1000,1500,20};
  float blo[45] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-10,-10,0,0,0,-10,-10,0,0,0,-10,-10,0,0,0,0};
  float bhi[45] = {200,500,500,200,200,200,500,10,3.2,3.2,500,200,500,10,3.2,3.2,500,200,500,10,3.2,3.2,500,200,200,10,100,100,100,10,10,10,50,15.0,10,10,10,50,15.0,10,10,10,100,15.0,20};
  for(int ievsel=0; ievsel<3; ievsel++){
    TString name1 = evsel_name[ievsel] + "flavor";
    h.flavor[ievsel] = new TH1F(name1,"0:#mu#mu#mu, 1:#mu#mue, 2:#mue#mu, 3:#muee, 4:eee, 5:e#mue, 6:ee#mu, 7:e#mu#mu",10,0.0,10.0);
    for(int iplot=0; iplot<45; iplot++){      
      TString name2 = evsel_name[ievsel] + plotname[iplot];
      //cout << "Creating histogram " << name2 << " with nbins = " << nbins[iplot] << ", blo = " << blo[iplot] << ", bhi = " << bhi[iplot] << endl;
      h.dispml_h[ievsel][iplot] = new TH1F(name2,name2,nbins[iplot],blo[iplot],bhi[iplot]);
     
      //h.bb_h[icr][iplot]->Sumw2();
    
    }
  }

  //2l1d flavor classified plots
  int n_bins[15] = {500,200,200,100,200,200,64,20,500,64,20,500,64,20,500};
  float b_lo[15] = {0,0,-10.0,0.0,0.0,0.0,-3.2,0,0,-3.2,0,0,-3.2,0,0};
  float b_hi[15] = {500,200,10.0,10.0,200.0,200,3.2,10,500,3.2,10,500,3.2,10,500};
  TString flav_type[2] = {"e", "mu"};
  TString plotnames[15] = {"M_3l", "met", "l2_dxy", "l2_ip3d", "l2_sip3d", "mt2", "dphi_l0l1", "dR_l0l1", "M_l0l1", "dphi_l1l2", "dR_l1l2", "M_l1l2", "dphi_l2l0", "dR_l2l0", "M_l2l0"};
  int p=0; int q=0;
  for(int flav=0; flav<2; flav++){
    for(int plot=0; plot<15; plot++){
      p=45;
      TString name = "2l1d_" + flav_type[flav] + "_" + plotnames[plot];
      h.dispml_h[0][plot+p+q] = new TH1F(name,name,n_bins[plot],b_lo[plot],b_hi[plot]);
    }
    q=15;
  }

  TString evsel_type[3] = {"2l1d_", "1l2d_", "3d_"};
  TString prefix[3] = {"l0_", "l1_", "l2_"};
  for(int ievsel=0; ievsel<3; ievsel++){
    for(int iplot = 0; iplot<3; iplot++){
      TString name = evsel_type[ievsel] + prefix[iplot] + "|dxy|";
      h.dispml_h[ievsel][75+iplot] = new TH1F(name,name,100,0,10);
    }
  }
  
}//BookHistograms()

  

