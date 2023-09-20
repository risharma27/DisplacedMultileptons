#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include "TLorentzVector.h"
#include "TVector.h"
#include "TVector2.h"
#include "TVector3.h"
#include <vector>
#include <fstream>
#include <iostream>
#include "MuonEfficiencyScaleFactor.h"



#ifndef ETAPHI_BINNING_CLASS_H
#define ETAPHI_BINNING_CLASS_H
class EtaPhiBinning {

    public:
        EtaPhiBinning(){ };
        virtual ~EtaPhiBinning() {;};

        virtual int bin(const TLorentzVector *m) const { 
            if(m->Eta() > 0)
                return this->symmetricBin(m) + 11;

            return 11 - this->symmetricBin(m);
        }
        virtual int symmetricBin(const TLorentzVector *m) const;

        enum binregion{ binUNKNOWN=0, bin1BARRELLG=1, bin1BARRELSM=2, bin2BARREL=3, binFEET=4,
			binTRANSITION=5, binENDCAPLG=6, binENDCAPSM=7, binBEE=8, binFORWARDLG=9, 
			binFORWARDSM=10};

    private:
        int getCoarseNSector(double phi) const;
        int getSector(double phi) const;
        int getECSector(double phi) const;

};

#endif


#ifndef StacoCBScaleEffFactorH
#define StacoCBScaleEffFactorH

class StacoCBScaleEffFactors : public MuonEfficiencyScaleFactor {
  public:
    //! Constructor
    StacoCBScaleEffFactors(void);
    //! Default constructor
    virtual ~StacoCBScaleEffFactors() {}

    // Methods //
    double scaleFactor(const TLorentzVector & tlv) const;
    ///< Get the efficiency scale factor for the given
    ///< fourmomentum. Scale factors for period after B are provided.
    double scaleFactorUncertainty(const TLorentzVector & tlv) const;
    ///< Get the uncertainty of the efficiency scale
    ///< factor for the given fourmomentum. Scale factors errors 
    ///< for period after B are provided.
    double scaleFactor(const TLorentzVector & tlv,
                                    const unsigned int & run_nb) const;
    ///< Get the efficiency scale factor for the given
    ///< fourmomentum and the given run number.
    double scaleFactorUncertainty(const TLorentzVector & tlv,
                                    const unsigned int & run_nb) const;
    ///< Get the uncertainty of the efficiency scale
    ///< factor for the given fourmomentum and the given run number..

    double scaleFactorSystematicUncertainty(const TLorentzVector & tlv) const;
    ///< Get the systematic uncertainty of the scale factor. The momentum
    ///< is assumed to be given in MeV.

  private:
    std::vector<double> m_scale_factor_A;
    std::vector<double> m_scale_factor_uncertainty_A;
    std::vector<double> m_scale_factor_C_high_pt;
    std::vector<double> m_scale_factor_uncertainty_C_high_pt;
    std::vector<double> m_scale_factor_C_low_pt;
    std::vector<double> m_scale_factor_uncertainty_C_low_pt;
    std::vector<double> m_scale_factor_A_periodB;
    std::vector<double> m_scale_factor_uncertainty_A_periodB;
    std::vector<double> m_scale_factor_C_periodB;
    std::vector<double> m_scale_factor_uncertainty_C_periodB;

    unsigned m_last_run_periodB;

    EtaPhiBinning m_eta_phi_binning; // auxiliary binning class
};
//}

#endif



//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue May 17 08:02:14 2011 by ROOT version 5.28/00b
// from TTree CollectionTree/CollectionTree
// found on file: /eliza18/atlas/sdube/2011_multiObj/user.sdube.rnd001_mc10_106051_pythiazmm.110516120452/user.sdube.000359.AANT._00003.root
//////////////////////////////////////////////////////////

#ifndef BbAna_h
#define BbAna_h



class BbAna {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           RunNumber;
   Int_t           EventNumber;
   Char_t          StreamESD_ref[153];
   Char_t          Stream1_ref;
   Char_t          Token[153];
   Int_t           Run;
   Int_t           Event;
   Int_t           Time;
   Int_t           LumiBlock;
   Int_t           BCID;
   Int_t           LVL1ID;
   Double_t        Weight;
   Int_t           IEvent;
   Int_t           StatusElement;
   Int_t           LVL1TriggerType;
   vector<unsigned int> *LVL1TriggerInfo;
   vector<unsigned int> *LVL2TriggerInfo;
   vector<unsigned int> *EventFilterInfo;
   vector<string>  *StreamTagName;
   vector<string>  *StreamTagType;
   Int_t           NVertex;
   vector<double>  *VertexVx;
   vector<double>  *VertexVy;
   vector<double>  *VertexVz;
   vector<int>     *VertexType;
   vector<int>     *VertexNtrk;
   Double_t        BeamPosX;
   Double_t        BeamPosY;
   Double_t        BeamPosZ;
   Int_t           NTrupart;
   vector<double>  *TrupartEta;
   vector<double>  *TrupartPhi;
   vector<double>  *TrupartPx;
   vector<double>  *TrupartPy;
   vector<double>  *TrupartPz;
   vector<double>  *TrupartPt;
   vector<double>  *TrupartE;
   vector<double>  *TrupartCharge;
   vector<double>  *Trupartptcone20;
   vector<int>     *TrupartPdgId;
   vector<int>     *TrupartMotherID;
   Int_t           NElectrons;
   vector<double>  *ElectronEta;
   vector<double>  *ElectronPhi;
   vector<double>  *ElectronPx;
   vector<double>  *ElectronPy;
   vector<double>  *ElectronPz;
   vector<double>  *ElectronPt;
   vector<double>  *ElectronE;
   vector<double>  *ElectronEt;
   vector<double>  *ElectronClusE;
   vector<double>  *ElectronClusEt;
   vector<double>  *ElectronClusEta;
   vector<double>  *ElectronClusPhi;
   vector<double>  *ElectronCharge;
   vector<double>  *Electrond0;
   vector<double>  *Electrond0err;
   vector<double>  *Electronz0;
   vector<double>  *Electronz0err;
   vector<double>  *Electrontheta;
   vector<double>  *Electronetcone30;
   vector<double>  *Electronptcone30;
   vector<double>  *Electronetcone20;
   vector<double>  *Electronptcone20;
   vector<int>     *ElectronExpBLayerHit;
   vector<int>     *ElectronBLayerHits;
   vector<int>     *ElectronPixHits;
   vector<int>     *ElectronSctHits;
   vector<int>     *ElectronTrtHits;
   vector<int>     *ElectronPixHoles;
   vector<int>     *ElectronSctHoles;
   vector<int>     *ElectronTrtOutliers;
   vector<int>     *ElectronPixDeadSensors;
   vector<int>     *ElectronSctDeadSensors;
   vector<double>  *ElectronSumpt03;
   vector<int>     *ElectronisTight;
   vector<int>     *ElectronPdgId;
   vector<int>     *ElectronMotherID;
   vector<int>     *ElectronMatchTrig;
   vector<int>     *ElectronNconv;
   vector<double>  *ElectronEoverP;
   vector<int>     *ElectronGoodOQ;
   Int_t           NMuons;
   vector<double>  *MuonEta;
   vector<double>  *MuonPhi;
   vector<double>  *MuonPx;
   vector<double>  *MuonPy;
   vector<double>  *MuonPz;
   vector<double>  *MuonPt;
   vector<double>  *MuonE;
   vector<double>  *Muond0;
   vector<double>  *Muond0err;
   vector<double>  *Muonz0;
   vector<double>  *Muonz0err;
   vector<double>  *Muontheta;
   vector<int>     *MuonCharge;
   vector<double>  *Muonetcone30;
   vector<double>  *Muonptcone30;
   vector<double>  *Muonetcone20;
   vector<double>  *Muonptcone20;
   vector<int>     *MuonExpBLayerHit;
   vector<int>     *MuonBLayerHits;
   vector<int>     *MuonSctHits;
   vector<int>     *MuonPixHits;
   vector<int>     *MuonTrtHits;
   vector<int>     *MuonSctHoles;
   vector<int>     *MuonPixHoles;
   vector<int>     *MuonTrtOutliers;
   vector<int>     *MuonPixDeadSensors;
   vector<int>     *MuonSctDeadSensors;
   vector<double>  *MuonSumpt03;
   vector<double>  *MuonPtms;
   vector<double>  *MuonPms;
   vector<double>  *MuonEtams;
   vector<double>  *MuonPhims;
   vector<double>  *MuonQms;
   vector<double>  *MuonPtmsextrap;
   vector<double>  *MuonPmsextrap;
   vector<double>  *MuonEtamsextrap;
   vector<double>  *MuonPhimsextrap;
   vector<double>  *MuonQmsextrap;
   vector<double>  *MuonPtid;
   vector<double>  *MuonPid;
   vector<double>  *MuonEtaid;
   vector<double>  *MuonPhiid;
   vector<double>  *MuonQid;
   vector<int>     *MuonPdgId;
   vector<int>     *MuonMotherID;
   vector<double>  *MuonMatchChi2;
   vector<int>     *MuonMatchTrig;
   Int_t           NBadJets;
   Int_t           NUglyJets;
   Int_t           NJets;
   vector<double>  *JetEta;
   vector<double>  *JetPhi;
   vector<double>  *JetPx;
   vector<double>  *JetPy;
   vector<double>  *JetPz;
   vector<double>  *JetPt;
   vector<double>  *JetE;
   vector<double>  *JetEt;
   vector<double>  *JetFlavorWeight;
   vector<double>  *JetFlavorWeightIp3d;
   vector<double>  *JetFlavorWeightSV1;
   vector<double>  *JetEmf;
   vector<double>  *JetQuality;
   vector<double>  *JetHecf;
   vector<double>  *Jetn90;
   vector<double>  *JetTime;
   vector<double>  *JetfracSamplingMax;
   Int_t           NTracks;
   vector<double>  *TrkEta;
   vector<double>  *TrkPhi;
   vector<double>  *TrkPx;
   vector<double>  *TrkPy;
   vector<double>  *TrkPz;
   vector<double>  *TrkPt;
   vector<double>  *TrkE;
   vector<double>  *Trkd0;
   vector<double>  *Trkd0err;
   vector<double>  *Trkz0;
   vector<double>  *Trkz0err;
   vector<double>  *Trktheta;
   vector<int>     *TrkCharge;
   vector<int>     *TrkExpBLayerHit;
   vector<int>     *TrkBLayerHits;
   vector<int>     *TrkSctHits;
   vector<int>     *TrkPixHits;
   vector<int>     *TrkTrtHits;
   vector<int>     *TrkSctHoles;
   vector<int>     *TrkPixHoles;
   vector<int>     *TrkTrtOutliers;
   vector<int>     *TrkPixDeadSensors;
   vector<int>     *TrkSctDeadSensors;
   vector<double>  *TrkSumpt03;
   Double_t        MissingET;
   Double_t        MissingETphi;
   Double_t        Met_RefFinal;
   Double_t        Met_RefFinalphi;
   Double_t        Met_LocHadTopo;
   Double_t        Met_LocHadTopophi;
   Double_t        Met_MuonBoy;
   Double_t        Met_MuonBoyphi;
   Double_t        Met_RefMuonTrack;
   Double_t        Met_RefMuonTrackphi;
   Int_t           NBJets;
   Double_t        PdfReweight;

   // List of branches
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_StreamESD_ref;   //!
   TBranch        *b_Stream1_ref;   //!
   TBranch        *b_Token;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_Time;   //!
   TBranch        *b_LumiBlock;   //!
   TBranch        *b_BCID;   //!
   TBranch        *b_LVL1ID;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_IEvent;   //!
   TBranch        *b_StatusElement;   //!
   TBranch        *b_LVL1TriggerType;   //!
   TBranch        *b_LVL1TriggerInfo;   //!
   TBranch        *b_LVL2TriggerInfo;   //!
   TBranch        *b_EventFilterInfo;   //!
   TBranch        *b_StreamTagName;   //!
   TBranch        *b_StreamTagType;   //!
   TBranch        *b_NVertex;   //!
   TBranch        *b_VertexVx;   //!
   TBranch        *b_VertexVy;   //!
   TBranch        *b_VertexVz;   //!
   TBranch        *b_VertexType;   //!
   TBranch        *b_VertexNtrk;   //!
   TBranch        *b_BeamPosX;   //!
   TBranch        *b_BeamPosY;   //!
   TBranch        *b_BeamPosZ;   //!
   TBranch        *b_NTrupart;   //!
   TBranch        *b_TrupartEta;   //!
   TBranch        *b_TrupartPhi;   //!
   TBranch        *b_TrupartPx;   //!
   TBranch        *b_TrupartPy;   //!
   TBranch        *b_TrupartPz;   //!
   TBranch        *b_TrupartPt;   //!
   TBranch        *b_TrupartE;   //!
   TBranch        *b_TrupartCharge;   //!
   TBranch        *b_Trupartptcone20;   //!
   TBranch        *b_TrupartPdgId;   //!
   TBranch        *b_TrupartMotherID;   //!
   TBranch        *b_NElectrons;   //!
   TBranch        *b_ElectronEta;   //!
   TBranch        *b_ElectronPhi;   //!
   TBranch        *b_ElectronPx;   //!
   TBranch        *b_ElectronPy;   //!
   TBranch        *b_ElectronPz;   //!
   TBranch        *b_ElectronPt;   //!
   TBranch        *b_ElectronE;   //!
   TBranch        *b_ElectronEt;   //!
   TBranch        *b_ElectronClusE;   //!
   TBranch        *b_ElectronClusEt;   //!
   TBranch        *b_ElectronClusEta;   //!
   TBranch        *b_ElectronClusPhi;   //!
   TBranch        *b_ElectronCharge;   //!
   TBranch        *b_Electrond0;   //!
   TBranch        *b_Electrond0err;   //!
   TBranch        *b_Electronz0;   //!
   TBranch        *b_Electronz0err;   //!
   TBranch        *b_Electrontheta;   //!
   TBranch        *b_Electronetcone30;   //!
   TBranch        *b_Electronptcone30;   //!
   TBranch        *b_Electronetcone20;   //!
   TBranch        *b_Electronptcone20;   //!
   TBranch        *b_ElectronExpBLayerHit;   //!
   TBranch        *b_ElectronBLayerHits;   //!
   TBranch        *b_ElectronPixHits;   //!
   TBranch        *b_ElectronSctHits;   //!
   TBranch        *b_ElectronTrtHits;   //!
   TBranch        *b_ElectronPixHoles;   //!
   TBranch        *b_ElectronSctHoles;   //!
   TBranch        *b_ElectronTrtOutliers;   //!
   TBranch        *b_ElectronPixDeadSensors;   //!
   TBranch        *b_ElectronSctDeadSensors;   //!
   TBranch        *b_ElectronSumpt03;   //!
   TBranch        *b_ElectronisTight;   //!
   TBranch        *b_ElectronPdgId;   //!
   TBranch        *b_ElectronMotherID;   //!
   TBranch        *b_ElectronMatchTrig;   //!
   TBranch        *b_ElectronNconv;   //!
   TBranch        *b_ElectronEoverP;   //!
   TBranch        *b_ElectronGoodOQ;   //!
   TBranch        *b_NMuons;   //!
   TBranch        *b_MuonEta;   //!
   TBranch        *b_MuonPhi;   //!
   TBranch        *b_MuonPx;   //!
   TBranch        *b_MuonPy;   //!
   TBranch        *b_MuonPz;   //!
   TBranch        *b_MuonPt;   //!
   TBranch        *b_MuonE;   //!
   TBranch        *b_Muond0;   //!
   TBranch        *b_Muond0err;   //!
   TBranch        *b_Muonz0;   //!
   TBranch        *b_Muonz0err;   //!
   TBranch        *b_Muontheta;   //!
   TBranch        *b_MuonCharge;   //!
   TBranch        *b_Muonetcone30;   //!
   TBranch        *b_Muonptcone30;   //!
   TBranch        *b_Muonetcone20;   //!
   TBranch        *b_Muonptcone20;   //!
   TBranch        *b_MuonExpBLayerHit;   //!
   TBranch        *b_MuonBLayerHits;   //!
   TBranch        *b_MuonSctHits;   //!
   TBranch        *b_MuonPixHits;   //!
   TBranch        *b_MuonTrtHits;   //!
   TBranch        *b_MuonSctHoles;   //!
   TBranch        *b_MuonPixHoles;   //!
   TBranch        *b_MuonTrtOutliers;   //!
   TBranch        *b_MuonPixDeadSensors;   //!
   TBranch        *b_MuonSctDeadSensors;   //!
   TBranch        *b_MuonSumpt03;   //!
   TBranch        *b_MuonPtms;   //!
   TBranch        *b_MuonPms;   //!
   TBranch        *b_MuonEtams;   //!
   TBranch        *b_MuonPhims;   //!
   TBranch        *b_MuonQms;   //!
   TBranch        *b_MuonPtmsextrap;   //!
   TBranch        *b_MuonPmsextrap;   //!
   TBranch        *b_MuonEtamsextrap;   //!
   TBranch        *b_MuonPhimsextrap;   //!
   TBranch        *b_MuonQmsextrap;   //!
   TBranch        *b_MuonPtid;   //!
   TBranch        *b_MuonPid;   //!
   TBranch        *b_MuonEtaid;   //!
   TBranch        *b_MuonPhiid;   //!
   TBranch        *b_MuonQid;   //!
   TBranch        *b_MuonPdgId;   //!
   TBranch        *b_MuonMotherID;   //!
   TBranch        *b_MuonMatchChi2;   //!
   TBranch        *b_MuonMatchTrig;   //!
   TBranch        *b_NBadJets;   //!
   TBranch        *b_NUglyJets;   //!
   TBranch        *b_NJets;   //!
   TBranch        *b_JetEta;   //!
   TBranch        *b_JetPhi;   //!
   TBranch        *b_JetPx;   //!
   TBranch        *b_JetPy;   //!
   TBranch        *b_JetPz;   //!
   TBranch        *b_JetPt;   //!
   TBranch        *b_JetE;   //!
   TBranch        *b_JetEt;   //!
   TBranch        *b_JetFlavorWeight;   //!
   TBranch        *b_JetFlavorWeightIp3d;   //!
   TBranch        *b_JetFlavorWeightSV1;   //!
   TBranch        *b_JetEmf;   //!
   TBranch        *b_JetQuality;   //!
   TBranch        *b_JetHecf;   //!
   TBranch        *b_Jetn90;   //!
   TBranch        *b_JetTime;   //!
   TBranch        *b_JetfracSamplingMax;   //!
   TBranch        *b_NTracks;   //!
   TBranch        *b_TrkEta;   //!
   TBranch        *b_TrkPhi;   //!
   TBranch        *b_TrkPx;   //!
   TBranch        *b_TrkPy;   //!
   TBranch        *b_TrkPz;   //!
   TBranch        *b_TrkPt;   //!
   TBranch        *b_TrkE;   //!
   TBranch        *b_Trkd0;   //!
   TBranch        *b_Trkd0err;   //!
   TBranch        *b_Trkz0;   //!
   TBranch        *b_Trkz0err;   //!
   TBranch        *b_Trktheta;   //!
   TBranch        *b_TrkCharge;   //!
   TBranch        *b_TrkExpBLayerHit;   //!
   TBranch        *b_TrkBLayerHits;   //!
   TBranch        *b_TrkSctHits;   //!
   TBranch        *b_TrkPixHits;   //!
   TBranch        *b_TrkTrtHits;   //!
   TBranch        *b_TrkSctHoles;   //!
   TBranch        *b_TrkPixHoles;   //!
   TBranch        *b_TrkTrtOutliers;   //!
   TBranch        *b_TrkPixDeadSensors;   //!
   TBranch        *b_TrkSctDeadSensors;   //!
   TBranch        *b_TrkSumpt03;   //!
   TBranch        *b_MissingET;   //!
   TBranch        *b_MissingETphi;   //!
   TBranch        *b_Met_RefFinal;   //!
   TBranch        *b_Met_RefFinalphi;   //!
   TBranch        *b_Met_LocHadTopo;   //!
   TBranch        *b_Met_LocHadTopophi;   //!
   TBranch        *b_Met_MuonBoy;   //!
   TBranch        *b_Met_MuonBoyphi;   //!
   TBranch        *b_Met_RefMuonTrack;   //!
   TBranch        *b_Met_RefMuonTrackphi;   //!
   TBranch        *b_NBJets;   //!
   TBranch        *b_PdfReweight;   //!

   BbAna(TTree *tree=0);
   virtual ~BbAna();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
  //My Methods
   void Begin();
   void End();
   void BookHistograms();
   void SetHstFileName(char *HstFileName) { _HstFileName = HstFileName; }
   void SetSumFileName(char *SumFileName) { _SumFileName = SumFileName; }
   void SetData(bool isthisdata) { isData = isthisdata; }
   void SetMixing(bool mixingon) { isMixingOn = mixingon; }

   float delta_phi(float, float);
   float delta_eta(float, float);
   float delta_r(float, float, float, float);
   float inv_mass(const float[4], const float[4]);
   float inv_mass(const vector<TLorentzVector> &);
   float transverse_mass_pt_phi(float, float, float, float);
   bool pass_track_cut(int);
   bool pass_muon_cut(int);
   bool pass_electron_cut(int);
   bool pass_isolation(int,int, double);
   double get_mu_trig_eff(float pt,float eta, float &);
   void Sort();
   int dobbPlots();

 public:
   struct Hists {

     //hists for bbar comparison
     /*

d0sig  |     |
       |     |
       |  B  |    C
       |_____|____________
       |  A  |    D
       |_____|____________
                      Isol

       Signal region is A+B with high NTracks
       Show that isolation is well modeled by comparing shape in B+C with low NTracks
       Show that NTracks   is well modeled by comparing shape in C+D, C
       f = (A+B)/(C+D) in MC.
       Then in data BB_1 = f * (C+D) with high NTracks.
       or BB_2 = (A+B) in low NTracks * fraction of high NTracks from MC.
     */
     //Region A,B,C,D : Plots d0sig1,d0sig2,isol1,isol2, pt1, pt2, met, dphi
     TH1F *whatsign,*mumotherid;
     TH2F *mom12;
     TH1F *trackpt;
     TH1F *bb_h[4][9];
     TH1F *bb_ptrat;
     TH2F *d0isol;
     TH1F *muonpt_bbcc[2][2];
   };
   struct Lepton {
     TLorentzVector mom;
     int id,ind;
   };
 protected:
   Hists h;
 
private:
   TFile *_HstFile;
   char *_HstFileName;
   char *_SumFileName;
   //char *_SumFileName;
   int nleptons, nelectrons, nmuons, nphotons, njets;
   int nEvtTotal;
   int nEvtWeight;
   int npasscuts[20];
   int bb_nev[4];
   float bb_wgt[4];
   Lepton lepsel[10], musel[10], elsel[10];
   int nlepsel, nmusel, nelsel;
   bool isData,isMixingOn;
   int NTRK_SEL,NBJET_SEL;
   float GEV, MEV2GEV;
   StacoCBScaleEffFactors StacoCBSCF; 
};

#endif

#ifdef BbAna_cxx
BbAna::BbAna(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eliza18/atlas/sdube/2011_multiObj/user.sdube.rnd001_mc10_106051_pythiazmm.110516120452/user.sdube.000359.AANT._00003.root");
      if (!f) {
         f = new TFile("/eliza18/atlas/sdube/2011_multiObj/user.sdube.rnd001_mc10_106051_pythiazmm.110516120452/user.sdube.000359.AANT._00003.root");
      }
      tree = (TTree*)gDirectory->Get("CollectionTree");

   }
   Init(tree);
}

BbAna::~BbAna()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t BbAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t BbAna::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void BbAna::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   LVL1TriggerInfo = 0;
   LVL2TriggerInfo = 0;
   EventFilterInfo = 0;
   StreamTagName = 0;
   StreamTagType = 0;
   VertexVx = 0;
   VertexVy = 0;
   VertexVz = 0;
   VertexType = 0;
   VertexNtrk = 0;
   TrupartEta = 0;
   TrupartPhi = 0;
   TrupartPx = 0;
   TrupartPy = 0;
   TrupartPz = 0;
   TrupartPt = 0;
   TrupartE = 0;
   TrupartCharge = 0;
   Trupartptcone20 = 0;
   TrupartPdgId = 0;
   TrupartMotherID = 0;
   ElectronEta = 0;
   ElectronPhi = 0;
   ElectronPx = 0;
   ElectronPy = 0;
   ElectronPz = 0;
   ElectronPt = 0;
   ElectronE = 0;
   ElectronEt = 0;
   ElectronClusE = 0;
   ElectronClusEt = 0;
   ElectronClusEta = 0;
   ElectronClusPhi = 0;
   ElectronCharge = 0;
   Electrond0 = 0;
   Electrond0err = 0;
   Electronz0 = 0;
   Electronz0err = 0;
   Electrontheta = 0;
   Electronetcone30 = 0;
   Electronptcone30 = 0;
   Electronetcone20 = 0;
   Electronptcone20 = 0;
   ElectronExpBLayerHit = 0;
   ElectronBLayerHits = 0;
   ElectronPixHits = 0;
   ElectronSctHits = 0;
   ElectronTrtHits = 0;
   ElectronPixHoles = 0;
   ElectronSctHoles = 0;
   ElectronTrtOutliers = 0;
   ElectronPixDeadSensors = 0;
   ElectronSctDeadSensors = 0;
   ElectronSumpt03 = 0;
   ElectronisTight = 0;
   ElectronPdgId = 0;
   ElectronMotherID = 0;
   ElectronMatchTrig = 0;
   ElectronNconv = 0;
   ElectronEoverP = 0;
   ElectronGoodOQ = 0;
   MuonEta = 0;
   MuonPhi = 0;
   MuonPx = 0;
   MuonPy = 0;
   MuonPz = 0;
   MuonPt = 0;
   MuonE = 0;
   Muond0 = 0;
   Muond0err = 0;
   Muonz0 = 0;
   Muonz0err = 0;
   Muontheta = 0;
   MuonCharge = 0;
   Muonetcone30 = 0;
   Muonptcone30 = 0;
   Muonetcone20 = 0;
   Muonptcone20 = 0;
   MuonExpBLayerHit = 0;
   MuonBLayerHits = 0;
   MuonSctHits = 0;
   MuonPixHits = 0;
   MuonTrtHits = 0;
   MuonSctHoles = 0;
   MuonPixHoles = 0;
   MuonTrtOutliers = 0;
   MuonPixDeadSensors = 0;
   MuonSctDeadSensors = 0;
   MuonSumpt03 = 0;
   MuonPtms = 0;
   MuonPms = 0;
   MuonEtams = 0;
   MuonPhims = 0;
   MuonQms = 0;
   MuonPtmsextrap = 0;
   MuonPmsextrap = 0;
   MuonEtamsextrap = 0;
   MuonPhimsextrap = 0;
   MuonQmsextrap = 0;
   MuonPtid = 0;
   MuonPid = 0;
   MuonEtaid = 0;
   MuonPhiid = 0;
   MuonQid = 0;
   MuonPdgId = 0;
   MuonMotherID = 0;
   MuonMatchChi2 = 0;
   MuonMatchTrig = 0;
   JetEta = 0;
   JetPhi = 0;
   JetPx = 0;
   JetPy = 0;
   JetPz = 0;
   JetPt = 0;
   JetE = 0;
   JetEt = 0;
   JetFlavorWeight = 0;
   JetFlavorWeightIp3d = 0;
   JetFlavorWeightSV1 = 0;
   JetEmf = 0;
   JetQuality = 0;
   JetHecf = 0;
   Jetn90 = 0;
   JetTime = 0;
   JetfracSamplingMax = 0;
   TrkEta = 0;
   TrkPhi = 0;
   TrkPx = 0;
   TrkPy = 0;
   TrkPz = 0;
   TrkPt = 0;
   TrkE = 0;
   Trkd0 = 0;
   Trkd0err = 0;
   Trkz0 = 0;
   Trkz0err = 0;
   Trktheta = 0;
   TrkCharge = 0;
   TrkExpBLayerHit = 0;
   TrkBLayerHits = 0;
   TrkSctHits = 0;
   TrkPixHits = 0;
   TrkTrtHits = 0;
   TrkSctHoles = 0;
   TrkPixHoles = 0;
   TrkTrtOutliers = 0;
   TrkPixDeadSensors = 0;
   TrkSctDeadSensors = 0;
   TrkSumpt03 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("StreamESD_ref", StreamESD_ref, &b_StreamESD_ref);
   fChain->SetBranchAddress("Stream1_ref", &Stream1_ref, &b_Stream1_ref);
   fChain->SetBranchAddress("Token", Token, &b_Token);
   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("Time", &Time, &b_Time);
   fChain->SetBranchAddress("LumiBlock", &LumiBlock, &b_LumiBlock);
   fChain->SetBranchAddress("BCID", &BCID, &b_BCID);
   fChain->SetBranchAddress("LVL1ID", &LVL1ID, &b_LVL1ID);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("IEvent", &IEvent, &b_IEvent);
   fChain->SetBranchAddress("StatusElement", &StatusElement, &b_StatusElement);
   fChain->SetBranchAddress("LVL1TriggerType", &LVL1TriggerType, &b_LVL1TriggerType);
   fChain->SetBranchAddress("LVL1TriggerInfo", &LVL1TriggerInfo, &b_LVL1TriggerInfo);
   fChain->SetBranchAddress("LVL2TriggerInfo", &LVL2TriggerInfo, &b_LVL2TriggerInfo);
   fChain->SetBranchAddress("EventFilterInfo", &EventFilterInfo, &b_EventFilterInfo);
   fChain->SetBranchAddress("StreamTagName", &StreamTagName, &b_StreamTagName);
   fChain->SetBranchAddress("StreamTagType", &StreamTagType, &b_StreamTagType);
   fChain->SetBranchAddress("NVertex", &NVertex, &b_NVertex);
   fChain->SetBranchAddress("VertexVx", &VertexVx, &b_VertexVx);
   fChain->SetBranchAddress("VertexVy", &VertexVy, &b_VertexVy);
   fChain->SetBranchAddress("VertexVz", &VertexVz, &b_VertexVz);
   fChain->SetBranchAddress("VertexType", &VertexType, &b_VertexType);
   fChain->SetBranchAddress("VertexNtrk", &VertexNtrk, &b_VertexNtrk);
   fChain->SetBranchAddress("BeamPosX", &BeamPosX, &b_BeamPosX);
   fChain->SetBranchAddress("BeamPosY", &BeamPosY, &b_BeamPosY);
   fChain->SetBranchAddress("BeamPosZ", &BeamPosZ, &b_BeamPosZ);
   fChain->SetBranchAddress("NTrupart", &NTrupart, &b_NTrupart);
   fChain->SetBranchAddress("TrupartEta", &TrupartEta, &b_TrupartEta);
   fChain->SetBranchAddress("TrupartPhi", &TrupartPhi, &b_TrupartPhi);
   fChain->SetBranchAddress("TrupartPx", &TrupartPx, &b_TrupartPx);
   fChain->SetBranchAddress("TrupartPy", &TrupartPy, &b_TrupartPy);
   fChain->SetBranchAddress("TrupartPz", &TrupartPz, &b_TrupartPz);
   fChain->SetBranchAddress("TrupartPt", &TrupartPt, &b_TrupartPt);
   fChain->SetBranchAddress("TrupartE", &TrupartE, &b_TrupartE);
   fChain->SetBranchAddress("TrupartCharge", &TrupartCharge, &b_TrupartCharge);
   fChain->SetBranchAddress("Trupartptcone20", &Trupartptcone20, &b_Trupartptcone20);
   fChain->SetBranchAddress("TrupartPdgId", &TrupartPdgId, &b_TrupartPdgId);
   fChain->SetBranchAddress("TrupartMotherID", &TrupartMotherID, &b_TrupartMotherID);
   fChain->SetBranchAddress("NElectrons", &NElectrons, &b_NElectrons);
   fChain->SetBranchAddress("ElectronEta", &ElectronEta, &b_ElectronEta);
   fChain->SetBranchAddress("ElectronPhi", &ElectronPhi, &b_ElectronPhi);
   fChain->SetBranchAddress("ElectronPx", &ElectronPx, &b_ElectronPx);
   fChain->SetBranchAddress("ElectronPy", &ElectronPy, &b_ElectronPy);
   fChain->SetBranchAddress("ElectronPz", &ElectronPz, &b_ElectronPz);
   fChain->SetBranchAddress("ElectronPt", &ElectronPt, &b_ElectronPt);
   fChain->SetBranchAddress("ElectronE", &ElectronE, &b_ElectronE);
   fChain->SetBranchAddress("ElectronEt", &ElectronEt, &b_ElectronEt);
   fChain->SetBranchAddress("ElectronClusE", &ElectronClusE, &b_ElectronClusE);
   fChain->SetBranchAddress("ElectronClusEt", &ElectronClusEt, &b_ElectronClusEt);
   fChain->SetBranchAddress("ElectronClusEta", &ElectronClusEta, &b_ElectronClusEta);
   fChain->SetBranchAddress("ElectronClusPhi", &ElectronClusPhi, &b_ElectronClusPhi);
   fChain->SetBranchAddress("ElectronCharge", &ElectronCharge, &b_ElectronCharge);
   fChain->SetBranchAddress("Electrond0", &Electrond0, &b_Electrond0);
   fChain->SetBranchAddress("Electrond0err", &Electrond0err, &b_Electrond0err);
   fChain->SetBranchAddress("Electronz0", &Electronz0, &b_Electronz0);
   fChain->SetBranchAddress("Electronz0err", &Electronz0err, &b_Electronz0err);
   fChain->SetBranchAddress("Electrontheta", &Electrontheta, &b_Electrontheta);
   fChain->SetBranchAddress("Electronetcone30", &Electronetcone30, &b_Electronetcone30);
   fChain->SetBranchAddress("Electronptcone30", &Electronptcone30, &b_Electronptcone30);
   fChain->SetBranchAddress("Electronetcone20", &Electronetcone20, &b_Electronetcone20);
   fChain->SetBranchAddress("Electronptcone20", &Electronptcone20, &b_Electronptcone20);
   fChain->SetBranchAddress("ElectronExpBLayerHit", &ElectronExpBLayerHit, &b_ElectronExpBLayerHit);
   fChain->SetBranchAddress("ElectronBLayerHits", &ElectronBLayerHits, &b_ElectronBLayerHits);
   fChain->SetBranchAddress("ElectronPixHits", &ElectronPixHits, &b_ElectronPixHits);
   fChain->SetBranchAddress("ElectronSctHits", &ElectronSctHits, &b_ElectronSctHits);
   fChain->SetBranchAddress("ElectronTrtHits", &ElectronTrtHits, &b_ElectronTrtHits);
   fChain->SetBranchAddress("ElectronPixHoles", &ElectronPixHoles, &b_ElectronPixHoles);
   fChain->SetBranchAddress("ElectronSctHoles", &ElectronSctHoles, &b_ElectronSctHoles);
   fChain->SetBranchAddress("ElectronTrtOutliers", &ElectronTrtOutliers, &b_ElectronTrtOutliers);
   fChain->SetBranchAddress("ElectronPixDeadSensors", &ElectronPixDeadSensors, &b_ElectronPixDeadSensors);
   fChain->SetBranchAddress("ElectronSctDeadSensors", &ElectronSctDeadSensors, &b_ElectronSctDeadSensors);
   fChain->SetBranchAddress("ElectronSumpt03", &ElectronSumpt03, &b_ElectronSumpt03);
   fChain->SetBranchAddress("ElectronisTight", &ElectronisTight, &b_ElectronisTight);
   fChain->SetBranchAddress("ElectronPdgId", &ElectronPdgId, &b_ElectronPdgId);
   fChain->SetBranchAddress("ElectronMotherID", &ElectronMotherID, &b_ElectronMotherID);
   fChain->SetBranchAddress("ElectronMatchTrig", &ElectronMatchTrig, &b_ElectronMatchTrig);
   fChain->SetBranchAddress("ElectronNconv", &ElectronNconv, &b_ElectronNconv);
   fChain->SetBranchAddress("ElectronEoverP", &ElectronEoverP, &b_ElectronEoverP);
   fChain->SetBranchAddress("ElectronGoodOQ", &ElectronGoodOQ, &b_ElectronGoodOQ);
   fChain->SetBranchAddress("NMuons", &NMuons, &b_NMuons);
   fChain->SetBranchAddress("MuonEta", &MuonEta, &b_MuonEta);
   fChain->SetBranchAddress("MuonPhi", &MuonPhi, &b_MuonPhi);
   fChain->SetBranchAddress("MuonPx", &MuonPx, &b_MuonPx);
   fChain->SetBranchAddress("MuonPy", &MuonPy, &b_MuonPy);
   fChain->SetBranchAddress("MuonPz", &MuonPz, &b_MuonPz);
   fChain->SetBranchAddress("MuonPt", &MuonPt, &b_MuonPt);
   fChain->SetBranchAddress("MuonE", &MuonE, &b_MuonE);
   fChain->SetBranchAddress("Muond0", &Muond0, &b_Muond0);
   fChain->SetBranchAddress("Muond0err", &Muond0err, &b_Muond0err);
   fChain->SetBranchAddress("Muonz0", &Muonz0, &b_Muonz0);
   fChain->SetBranchAddress("Muonz0err", &Muonz0err, &b_Muonz0err);
   fChain->SetBranchAddress("Muontheta", &Muontheta, &b_Muontheta);
   fChain->SetBranchAddress("MuonCharge", &MuonCharge, &b_MuonCharge);
   fChain->SetBranchAddress("Muonetcone30", &Muonetcone30, &b_Muonetcone30);
   fChain->SetBranchAddress("Muonptcone30", &Muonptcone30, &b_Muonptcone30);
   fChain->SetBranchAddress("Muonetcone20", &Muonetcone20, &b_Muonetcone20);
   fChain->SetBranchAddress("Muonptcone20", &Muonptcone20, &b_Muonptcone20);
   fChain->SetBranchAddress("MuonExpBLayerHit", &MuonExpBLayerHit, &b_MuonExpBLayerHit);
   fChain->SetBranchAddress("MuonBLayerHits", &MuonBLayerHits, &b_MuonBLayerHits);
   fChain->SetBranchAddress("MuonSctHits", &MuonSctHits, &b_MuonSctHits);
   fChain->SetBranchAddress("MuonPixHits", &MuonPixHits, &b_MuonPixHits);
   fChain->SetBranchAddress("MuonTrtHits", &MuonTrtHits, &b_MuonTrtHits);
   fChain->SetBranchAddress("MuonSctHoles", &MuonSctHoles, &b_MuonSctHoles);
   fChain->SetBranchAddress("MuonPixHoles", &MuonPixHoles, &b_MuonPixHoles);
   fChain->SetBranchAddress("MuonTrtOutliers", &MuonTrtOutliers, &b_MuonTrtOutliers);
   fChain->SetBranchAddress("MuonPixDeadSensors", &MuonPixDeadSensors, &b_MuonPixDeadSensors);
   fChain->SetBranchAddress("MuonSctDeadSensors", &MuonSctDeadSensors, &b_MuonSctDeadSensors);
   fChain->SetBranchAddress("MuonSumpt03", &MuonSumpt03, &b_MuonSumpt03);
   fChain->SetBranchAddress("MuonPtms", &MuonPtms, &b_MuonPtms);
   fChain->SetBranchAddress("MuonPms", &MuonPms, &b_MuonPms);
   fChain->SetBranchAddress("MuonEtams", &MuonEtams, &b_MuonEtams);
   fChain->SetBranchAddress("MuonPhims", &MuonPhims, &b_MuonPhims);
   fChain->SetBranchAddress("MuonQms", &MuonQms, &b_MuonQms);
   fChain->SetBranchAddress("MuonPtmsextrap", &MuonPtmsextrap, &b_MuonPtmsextrap);
   fChain->SetBranchAddress("MuonPmsextrap", &MuonPmsextrap, &b_MuonPmsextrap);
   fChain->SetBranchAddress("MuonEtamsextrap", &MuonEtamsextrap, &b_MuonEtamsextrap);
   fChain->SetBranchAddress("MuonPhimsextrap", &MuonPhimsextrap, &b_MuonPhimsextrap);
   fChain->SetBranchAddress("MuonQmsextrap", &MuonQmsextrap, &b_MuonQmsextrap);
   fChain->SetBranchAddress("MuonPtid", &MuonPtid, &b_MuonPtid);
   fChain->SetBranchAddress("MuonPid", &MuonPid, &b_MuonPid);
   fChain->SetBranchAddress("MuonEtaid", &MuonEtaid, &b_MuonEtaid);
   fChain->SetBranchAddress("MuonPhiid", &MuonPhiid, &b_MuonPhiid);
   fChain->SetBranchAddress("MuonQid", &MuonQid, &b_MuonQid);
   fChain->SetBranchAddress("MuonPdgId", &MuonPdgId, &b_MuonPdgId);
   fChain->SetBranchAddress("MuonMotherID", &MuonMotherID, &b_MuonMotherID);
   fChain->SetBranchAddress("MuonMatchChi2", &MuonMatchChi2, &b_MuonMatchChi2);
   fChain->SetBranchAddress("MuonMatchTrig", &MuonMatchTrig, &b_MuonMatchTrig);
   fChain->SetBranchAddress("NBadJets", &NBadJets, &b_NBadJets);
   fChain->SetBranchAddress("NUglyJets", &NUglyJets, &b_NUglyJets);
   fChain->SetBranchAddress("NJets", &NJets, &b_NJets);
   fChain->SetBranchAddress("JetEta", &JetEta, &b_JetEta);
   fChain->SetBranchAddress("JetPhi", &JetPhi, &b_JetPhi);
   fChain->SetBranchAddress("JetPx", &JetPx, &b_JetPx);
   fChain->SetBranchAddress("JetPy", &JetPy, &b_JetPy);
   fChain->SetBranchAddress("JetPz", &JetPz, &b_JetPz);
   fChain->SetBranchAddress("JetPt", &JetPt, &b_JetPt);
   fChain->SetBranchAddress("JetE", &JetE, &b_JetE);
   fChain->SetBranchAddress("JetEt", &JetEt, &b_JetEt);
   fChain->SetBranchAddress("JetFlavorWeight", &JetFlavorWeight, &b_JetFlavorWeight);
   fChain->SetBranchAddress("JetFlavorWeightIp3d", &JetFlavorWeightIp3d, &b_JetFlavorWeightIp3d);
   fChain->SetBranchAddress("JetFlavorWeightSV1", &JetFlavorWeightSV1, &b_JetFlavorWeightSV1);
   fChain->SetBranchAddress("JetEmf", &JetEmf, &b_JetEmf);
   fChain->SetBranchAddress("JetQuality", &JetQuality, &b_JetQuality);
   fChain->SetBranchAddress("JetHecf", &JetHecf, &b_JetHecf);
   fChain->SetBranchAddress("Jetn90", &Jetn90, &b_Jetn90);
   fChain->SetBranchAddress("JetTime", &JetTime, &b_JetTime);
   fChain->SetBranchAddress("JetfracSamplingMax", &JetfracSamplingMax, &b_JetfracSamplingMax);
   fChain->SetBranchAddress("NTracks", &NTracks, &b_NTracks);
   fChain->SetBranchAddress("TrkEta", &TrkEta, &b_TrkEta);
   fChain->SetBranchAddress("TrkPhi", &TrkPhi, &b_TrkPhi);
   fChain->SetBranchAddress("TrkPx", &TrkPx, &b_TrkPx);
   fChain->SetBranchAddress("TrkPy", &TrkPy, &b_TrkPy);
   fChain->SetBranchAddress("TrkPz", &TrkPz, &b_TrkPz);
   fChain->SetBranchAddress("TrkPt", &TrkPt, &b_TrkPt);
   fChain->SetBranchAddress("TrkE", &TrkE, &b_TrkE);
   fChain->SetBranchAddress("Trkd0", &Trkd0, &b_Trkd0);
   fChain->SetBranchAddress("Trkd0err", &Trkd0err, &b_Trkd0err);
   fChain->SetBranchAddress("Trkz0", &Trkz0, &b_Trkz0);
   fChain->SetBranchAddress("Trkz0err", &Trkz0err, &b_Trkz0err);
   fChain->SetBranchAddress("Trktheta", &Trktheta, &b_Trktheta);
   fChain->SetBranchAddress("TrkCharge", &TrkCharge, &b_TrkCharge);
   fChain->SetBranchAddress("TrkExpBLayerHit", &TrkExpBLayerHit, &b_TrkExpBLayerHit);
   fChain->SetBranchAddress("TrkBLayerHits", &TrkBLayerHits, &b_TrkBLayerHits);
   fChain->SetBranchAddress("TrkSctHits", &TrkSctHits, &b_TrkSctHits);
   fChain->SetBranchAddress("TrkPixHits", &TrkPixHits, &b_TrkPixHits);
   fChain->SetBranchAddress("TrkTrtHits", &TrkTrtHits, &b_TrkTrtHits);
   fChain->SetBranchAddress("TrkSctHoles", &TrkSctHoles, &b_TrkSctHoles);
   fChain->SetBranchAddress("TrkPixHoles", &TrkPixHoles, &b_TrkPixHoles);
   fChain->SetBranchAddress("TrkTrtOutliers", &TrkTrtOutliers, &b_TrkTrtOutliers);
   fChain->SetBranchAddress("TrkPixDeadSensors", &TrkPixDeadSensors, &b_TrkPixDeadSensors);
   fChain->SetBranchAddress("TrkSctDeadSensors", &TrkSctDeadSensors, &b_TrkSctDeadSensors);
   fChain->SetBranchAddress("TrkSumpt03", &TrkSumpt03, &b_TrkSumpt03);
   fChain->SetBranchAddress("MissingET", &MissingET, &b_MissingET);
   fChain->SetBranchAddress("MissingETphi", &MissingETphi, &b_MissingETphi);
   fChain->SetBranchAddress("Met_RefFinal", &Met_RefFinal, &b_Met_RefFinal);
   fChain->SetBranchAddress("Met_RefFinalphi", &Met_RefFinalphi, &b_Met_RefFinalphi);
   fChain->SetBranchAddress("Met_LocHadTopo", &Met_LocHadTopo, &b_Met_LocHadTopo);
   fChain->SetBranchAddress("Met_LocHadTopophi", &Met_LocHadTopophi, &b_Met_LocHadTopophi);
   fChain->SetBranchAddress("Met_MuonBoy", &Met_MuonBoy, &b_Met_MuonBoy);
   fChain->SetBranchAddress("Met_MuonBoyphi", &Met_MuonBoyphi, &b_Met_MuonBoyphi);
   fChain->SetBranchAddress("Met_RefMuonTrack", &Met_RefMuonTrack, &b_Met_RefMuonTrack);
   fChain->SetBranchAddress("Met_RefMuonTrackphi", &Met_RefMuonTrackphi, &b_Met_RefMuonTrackphi);
   fChain->SetBranchAddress("NBJets", &NBJets, &b_NBJets);
   fChain->SetBranchAddress("PdfReweight", &PdfReweight, &b_PdfReweight);
   Notify();
}

Bool_t BbAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void BbAna::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t BbAna::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef BbAna_cxx



