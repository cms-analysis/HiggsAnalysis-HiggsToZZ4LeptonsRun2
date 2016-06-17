//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May 18 13:01:53 2016 by ROOT version 6.02/13
// from TTree HZZ4LeptonsAnalysis/HZZ4Leptons Analysis Tree
// found on file: /localdata/Syncr13TeV/roottree_leptons_sync_Fall15_HiggsToZZ_76x_vector.root
//////////////////////////////////////////////////////////

#ifndef MonoHiggsAnalysis4mu_h
#define MonoHiggsAnalysis4mu_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <TLorentzVector.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"

using namespace std;

class MonoHiggsAnalysis4mu {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          Run;
   UInt_t          Event;
   UInt_t          LumiSection;
   Float_t         Avginstlumi;
   Int_t           RECO_NVTX;
   vector<float>   *RECO_VERTEX_x;
   vector<float>   *RECO_VERTEX_y;
   vector<float>   *RECO_VERTEX_z;
   vector<float>   *RECO_VERTEX_ndof;
   vector<float>   *RECO_VERTEX_chi2;
   vector<float>   *RECO_VERTEX_ntracks;
   vector<float>   *RECO_VERTEXPROB;
   vector<bool>    *RECO_VERTEX_isValid;
   vector<float>   *RECO_VERTEX_TRACK_PT;
   Int_t           RECO_NTRACK;
   vector<float>   *RECO_TRACK_PT;
   vector<float>   *RECO_TRACK_ETA;
   vector<float>   *RECO_TRACK_PHI;
   vector<float>   *RECO_TRACK_CHI2;
   vector<float>   *RECO_TRACK_CHI2RED;
   vector<float>   *RECO_TRACK_CHI2PROB;
   vector<int>     *RECO_TRACK_NHITS;
   vector<float>   *RECO_TRACK_DXY;
   vector<float>   *RECO_TRACK_DXYERR;
   vector<float>   *RECO_TRACK_DZ;
   vector<float>   *RECO_TRACK_DZERR;
   Float_t         MC_GENMET;
   Char_t          HLTPathsFired;
   Int_t           RECO_nMuHLTMatch;
   Int_t           RECO_nEleHLTMatch;
   vector<float>   *RECOMU_PT_MuHLTMatch;
   vector<float>   *RECOMU_ETA_MuHLTMatch;
   vector<float>   *RECOMU_PHI_MuHLTMatch;
   vector<float>   *RECOELE_PT_EleHLTMatch;
   vector<float>   *RECOELE_ETA_EleHLTMatch;
   vector<float>   *RECOELE_PHI_EleHLTMatch;
   vector<float>   *MC_LEPT_PT;
   vector<float>   *MC_LEPT_ETA;
   vector<float>   *MC_LEPT_PHI;
   vector<float>   *MC_LEPT_THETA;
   vector<int>     *MC_LEPT_PDGID;
   vector<float>   *MC_Z_PT;
   vector<float>   *MC_Z_ETA;
   vector<float>   *MC_Z_PHI;
   vector<float>   *MC_Z_THETA;
   vector<float>   *MC_Z_MASS;
   vector<int>     *MC_Z_PDGID;
   vector<float>   *MC_fourl_MASS;
   vector<float>   *MC_fourl_PT;
   vector<int>     *MC_fourl_PDGID;
   vector<float>   *MC_ZZ_MASS;
   vector<float>   *MC_ZZ_PT;
   vector<float>   *MC_ZZ_ETA;
   vector<float>   *MC_ZZ_PHI;
   vector<float>   *MC_ZZ_THETA;
   vector<int>     *MC_ZZ_PDGID;
   vector<float>   *MC_E;
   vector<float>   *MC_PT;
   vector<float>   *MC_ETA;
   vector<float>   *MC_THETA;
   vector<float>   *MC_PHI;
   vector<float>   *MC_MASS;
   vector<int>     *MC_PDGID;
   Int_t           RECO_NMU;
   vector<float>   *RECOMU_isPFMu;
   vector<float>   *RECOMU_isGlobalMu;
   vector<float>   *RECOMU_isStandAloneMu;
   vector<float>   *RECOMU_isTrackerMu;
   vector<float>   *RECOMU_isCaloMu;
   vector<float>   *RECOMU_E;
   vector<float>   *RECOMU_PT;
   vector<float>   *RECOMU_P;
   vector<float>   *RECOMU_ETA;
   vector<float>   *RECOMU_THETA;
   vector<float>   *RECOMU_PHI;
   vector<float>   *RECOMU_MASS;
   vector<int>     *RECOMU_CHARGE;
   vector<float>   *RECOMU_COV;
   vector<float>   *RECOMU_TRACKISO;
   vector<float>   *RECOMU_TRACKISO_SUMPT;
   vector<float>   *RECOMU_ECALISO;
   vector<float>   *RECOMU_HCALISO;
   vector<float>   *RECOMU_X;
   vector<float>   *RECOMU_PFchHad;
   vector<float>   *RECOMU_PFneuHad;
   vector<float>   *RECOMU_PFphoton;
   vector<float>   *RECOMU_PFPUchAllPart;
   vector<float>   *RECOMU_PFX_dB;
   vector<float>    RECOMU_PFX_dB_new;
   vector<float>   *RECOMU_PFX_rho;
   vector<float>   *RECOMU_SIP;
   vector<float>   *RECOMU_IP;
   vector<float>   *RECOMU_IPERROR;
   vector<float>   *RECOMU_SIP_KF;
   vector<float>   *RECOMU_IP_KF;
   vector<float>   *RECOMU_IPERROR_KF;
   vector<float>   *RECOMU_STIP;
   vector<float>   *RECOMU_SLIP;
   vector<float>   *RECOMU_TIP;
   vector<float>   *RECOMU_LIP;
   vector<float>   *RECOMU_TIPERROR;
   vector<float>   *RECOMU_LIPERROR;
   vector<float>   *RECOMU_numberOfMatches;
   vector<float>   *RECOMU_numberOfMatchedStations;
   vector<float>   *RECOMU_caloCompatibility;
   vector<float>   *RECOMU_segmentCompatibility;
   vector<float>   *RECOMU_glbmuPromptTight;
   vector<float>   *RECOMU_mubesttrkType;
   vector<float>   *RECOMU_mubesttrkDxy;
   vector<float>   *RECOMU_mubesttrkDxyB;
   vector<float>   *RECOMU_mubesttrkDxyError;
   vector<float>   *RECOMU_mubesttrkDz;
   vector<float>   *RECOMU_mubesttrkDzB;
   vector<float>   *RECOMU_mubesttrkDzError;
   vector<float>   *RECOMU_mubesttrkPTError;
   vector<float>   *RECOMU_mutrkPT;
   vector<float>   *RECOMU_mutrkPTError;
   vector<float>   *RECOMU_mutrkDxy;
   vector<float>   *RECOMU_mutrkDxyError;
   vector<float>   *RECOMU_mutrkDxyB;
   vector<float>   *RECOMU_mutrkDz;
   vector<float>   *RECOMU_mutrkDzError;
   vector<float>   *RECOMU_mutrkDzB;
   vector<float>   *RECOMU_mutrkChi2PerNdof;
   vector<int>     *RECOMU_mutrkCharge;
   vector<int>     *RECOMU_mutrkNHits;
   vector<int>     *RECOMU_mutrkNPixHits;
   vector<int>     *RECOMU_mutrkNStripHits;
   vector<int>     *RECOMU_mutrkNMuonHits;
   vector<int>     *RECOMU_mutrktrackerLayersWithMeasurement;
   vector<float>   *RECOMU_muInnertrkDxy;
   vector<float>   *RECOMU_muInnertrkDxyError;
   vector<float>   *RECOMU_muInnertrkDxyB;
   vector<float>   *RECOMU_muInnertrkDz;
   vector<float>   *RECOMU_muInnertrkDzError;
   vector<float>   *RECOMU_muInnertrkDzB;
   vector<float>   *RECOMU_muInnertrkChi2PerNdof;
   vector<float>   *RECOMU_muInnertrktrackerLayersWithMeasurement;
   vector<float>   *RECOMU_muInnertrkPT;
   vector<float>   *RECOMU_muInnertrkPTError;
   vector<int>     *RECOMU_muInnertrkCharge;
   vector<int>     *RECOMU_muInnertrkNHits;
   vector<int>     *RECOMU_muInnertrkNPixHits;
   vector<int>     *RECOMU_muInnertrkNStripHits;
   vector<float>   *RECOMU_trkmuArbitration;
   vector<float>   *RECOMU_trkmu2DCompatibilityLoose;
   vector<float>   *RECOMU_trkmu2DCompatibilityTight;
   vector<float>   *RECOMU_trkmuOneStationLoose;
   vector<float>   *RECOMU_trkmuOneStationTight;
   vector<float>   *RECOMU_trkmuLastStationLoose;
   vector<float>   *RECOMU_trkmuLastStationTight;
   vector<float>   *RECOMU_trkmuOneStationAngLoose;
   vector<float>   *RECOMU_trkmuOneStationAngTight;
   vector<float>   *RECOMU_trkmuLastStationAngLoose;
   vector<float>   *RECOMU_trkmuLastStationAngTight;
   vector<float>   *RECOMU_trkmuLastStationOptimizedLowPtLoose;
   vector<float>   *RECOMU_trkmuLastStationOptimizedLowPtTight;
   vector<int>     *RECOMU_MatchingMCTruth;
   vector<float>   *RECOMU_MatchingMCpT;
   vector<float>   *RECOMU_MatchingMCEta;
   vector<float>   *RECOMU_MatchingMCPhi;
   Int_t           RECO_NELE;
   vector<int>     *RECOELE_isEcalDriven;
   vector<int>     *RECOELE_isTrackerDriven;
   vector<int>     *RECOELE_CHARGE;
   vector<float>   *RECOELE_E;
   vector<float>   *RECOELE_PT;
   vector<float>   *RECOELE_P;
   vector<float>   *RECOELE_ETA;
   vector<float>   *RECOELE_THETA;
   vector<float>   *RECOELE_PHI;
   vector<float>   *RECOELE_MASS;
   vector<float>   *ele_sclRawE;
   vector<float>   *RECOELE_scl_E;
   vector<float>   *RECOELE_scl_Et;
   vector<float>   *RECOELE_scl_Eta;
   vector<float>   *RECOELE_scl_Phi;
   vector<float>   *ele_sclX;
   vector<float>   *ele_sclY;
   vector<float>   *ele_sclZ;
   vector<float>   *RECOELE_PTError;
   vector<float>   *RECOELE_COV;
   vector<float>   *RECOELE_EGMECALISO;
   vector<float>   *RECOELE_EGMHCALISO;
   vector<float>   *RECOELE_EGMX;
   vector<float>   *RECOELE_PFchAllPart;
   vector<float>   *RECOELE_PFchHad;
   vector<float>   *RECOELE_PFneuHad;
   vector<float>   *RECOELE_PFphoton;
   vector<float>   *RECOELE_PFPUchAllPart;
   vector<float>   *RECOELE_PFX_dB;
   vector<float>   *RECOELE_PFX_rho;
   vector<float>   RECOELE_PFX_rho_new;
   vector<float>   *RECOELE_SIP;
   vector<float>   *RECOELE_IP;
   vector<float>   *RECOELE_IPERROR;
   vector<float>   *RECOELE_SIP_KF;
   vector<float>   *RECOELE_IP_KF;
   vector<float>   *RECOELE_IPERROR_KF;
   vector<float>   *RECOELE_STIP;
   vector<float>   *RECOELE_SLIP;
   vector<float>   *RECOELE_TIP;
   vector<float>   *RECOELE_LIP;
   vector<float>   *RECOELE_TIPERROR;
   vector<float>   *RECOELE_LIPERROR;
   vector<float>   *RECOELE_gsftrack_NPixHits;
   vector<float>   *RECOELE_gsftrack_NStripHits;
   vector<float>   *RECOELE_gsftrack_chi2;
   vector<float>   *RECOELE_gsftrack_dxyB;
   vector<float>   *RECOELE_gsftrack_dxy;
   vector<float>   *RECOELE_gsftrack_dxyError;
   vector<float>   *RECOELE_gsftrack_dzB;
   vector<float>   *RECOELE_gsftrack_dz;
   vector<float>   *RECOELE_gsftrack_dzError;
   vector<int>     *RECOELE_gsftrack_losthits;
   vector<int>     *RECOELE_gsftrack_validhits;
   vector<int>     *RECOELE_gsftrack_expected_inner_hits;
   vector<float>   *RECOELE_ep;
   vector<float>   *RECOELE_eSeedp;
   vector<float>   *RECOELE_eSeedpout;
   vector<float>   *RECOELE_eElepout;
   vector<float>   *RECOELE_deltaEtaIn;
   vector<float>   *RECOELE_deltaEtaSeed;
   vector<float>   *RECOELE_deltaEtaEle;
   vector<float>   *RECOELE_deltaPhiIn;
   vector<float>   *RECOELE_deltaPhiSeed;
   vector<float>   *RECOELE_deltaPhiEle;
   vector<int>     *RECOELE_isbarrel;
   vector<int>     *RECOELE_isendcap;
   vector<int>     *RECOELE_isGap;
   vector<int>     *RECOELE_isEBetaGap;
   vector<int>     *RECOELE_isEBphiGap;
   vector<int>     *RECOELE_isEEdeeGap;
   vector<int>     *RECOELE_isEEringGap;
   vector<float>   *RECOELE_sigmaIetaIeta;
   vector<float>   *RECOELE_sigmaEtaEta;
   vector<float>   *RECOELE_e15;
   vector<float>   *RECOELE_e25max;
   vector<float>   *RECOELE_e55;
   vector<float>   *RECOELE_he;
   vector<float>   *RECOELE_fbrem;
   vector<float>   *RECOELE_nbrems;
   vector<float>   *RECOELE_Class;
   vector<float>   *RECOELE_fbrem_mean;
   vector<float>   *RECOELE_fbrem_mode;
   vector<float>   *RECOELE_ecalEnergy;
   vector<int>     *ele_seedSubdet2;
   vector<double>  *ele_seedDphi2;
   vector<double>  *ele_seedDrz2;
   vector<int>     *ele_seedSubdet1;
   vector<double>  *ele_seedDphi1;
   vector<double>  *ele_seedDrz1;
   vector<float>   *RECOELE_mvaTrigV0;
   vector<float>   *RECOELE_mvaNonTrigV0;
   vector<float>   *ConvMapDist;
   vector<float>   *ConvMapDcot;
   vector<int>     *RECOELE_MatchingMCTruth;
   vector<float>   *RECOELE_MatchingMCpT;
   vector<float>   *RECOELE_MatchingMCEta;
   vector<float>   *RECOELE_MatchingMCPhi;
   Int_t           RECO_NPHOT;
   vector<float>   *RECOPHOT_PT;
   vector<float>   *RECOPHOT_ETA;
   vector<float>   *RECOPHOT_PHI;
   vector<float>   *RECOPHOT_THETA;
   Int_t           RECO_NPFPHOT;
   vector<float>   *RECOPFPHOT_PT;
   vector<float>   *RECOPFPHOT_PTError;
   vector<float>   *RECOPFPHOT_ETA;
   vector<float>   *RECOPFPHOT_PHI;
   vector<float>   *RECOPFPHOT_THETA;
   vector<float>   *RECOPFPHOT_PFchAllPart;
   vector<float>   *RECOPFPHOT_PFchHad;
   vector<float>   *RECOPFPHOT_PFneuHad;
   vector<float>   *RECOPFPHOT_PFphoton;
   vector<float>   *RECOPFPHOT_PFPUchAllPart;
   vector<float>   *RECOPFPHOT_PFX_rho;
   vector<int>     *RECOPHOT_MatchingMCTruth;
   vector<float>   *RECOPHOT_MatchingMCpT;
   vector<float>   *RECOPHOT_MatchingMCEta;
   vector<float>   *RECOPHOT_MatchingMCPhi;
   Float_t         PFMet_pt;
   Float_t         PFMet_eta;
   Float_t         PFMet_phi;
   Float_t         PFMet_en;
   Float_t         PFMet_px;
   Float_t         PFMet_py;
   Float_t         PFMet_pz;
   Float_t         PFMet_sumEt;
   Double_t        METSign;
   Int_t           tCHighEff_nb;
   vector<float>   *tCHighEff_BTagJet_PT;
   vector<float>   *tCHighEff_BTagJet_ETA;
   vector<float>   *tCHighEff_BTagJet_PHI;
   vector<float>   *tCHighEff_BTagJet_DISCR;
   Int_t           tCHighPur_nb;
   vector<float>   *tCHighPur_BTagJet_PT;
   vector<float>   *tCHighPur_BTagJet_ETA;
   vector<float>   *tCHighPur_BTagJet_PHI;
   vector<float>   *tCHighPur_BTagJet_DISCR;
   Int_t           cSV_nb;
   vector<float>   *cSV_BTagJet_PT;
   vector<float>   *cSV_BTagJet_ETA;
   vector<float>   *cSV_BTagJet_PHI;
   vector<float>   *cSV_BTagJet_DISCR;
   vector<float>   *cSV_BTagJet_ET;
   vector<float>   *RECO_ZMM_MASS;
   vector<float>   *RECO_ZEE_MASS;
   vector<float>   *RECO_DiLep_MASS;
   vector<float>   *RECO_ZMM_PT;
   vector<float>   *RECO_ZEE_PT;
   vector<float>   *RECO_DiLep_PT;
   vector<float>   *RECO_ZMM_ETA;
   vector<float>   *RECO_ZEE_ETA;
   vector<float>   *RECO_DiLep_ETA;
   vector<float>   *RECO_ZMM_PHI;
   vector<float>   *RECO_ZEE_PHI;
   vector<float>   *RECO_DiLep_PHI;
   vector<float>   *RECO_ZMMss_MASS;
   vector<float>   *RECO_ZEEss_MASS;
   vector<float>   *RECO_ZEM_MASS;
   vector<float>   *RECO_ZMMss_PT;
   vector<float>   *RECO_ZEEss_PT;
   vector<float>   *RECO_ZEM_PT;
   vector<float>   *RECO_ZMMss_ETA;
   vector<float>   *RECO_ZEEss_ETA;
   vector<float>   *RECO_ZEM_ETA;
   vector<float>   *RECO_ZMMss_PHI;
   vector<float>   *RECO_ZEEss_PHI;
   vector<float>   *RECO_ZEM_PHI;
   vector<float>   *RECO_MMMM_MASS;
   vector<float>   *RECO_MMMM_PT;
   vector<float>   *RECO_MMMM_ETA;
   vector<float>   *RECO_MMMM_PHI;
   vector<float>   *RECO_MMMM_MASS_REFIT;
   vector<float>   *RECO_EEEE_MASS;
   vector<float>   *RECO_EEEE_PT;
   vector<float>   *RECO_EEEE_ETA;
   vector<float>   *RECO_EEEE_PHI;
   vector<float>   *RECO_EEEE_MASS_REFIT;
   vector<float>   *RECO_EEMM_MASS;
   vector<float>   *RECO_EEMM_PT;
   vector<float>   *RECO_EEMM_ETA;
   vector<float>   *RECO_EEMM_PHI;
   vector<float>   *RECO_EEMM_MASS_REFIT;
   vector<float>   *RECO_LLL0_MASS;
   vector<float>   *RECO_LLL1_MASS;
   vector<float>   *RECO_LLL2_MASS;
   vector<float>   *RECO_LLL3_MASS;
   vector<float>   *RECO_LLL0_PT;
   vector<float>   *RECO_LLL1_PT;
   vector<float>   *RECO_LLL2_PT;
   vector<float>   *RECO_LLL3_PT;
   vector<float>   *RECO_LLLl0_MASS;
   vector<float>   *RECO_LLLl1_MASS;
   vector<float>   *RECO_LLLl0_PT;
   vector<float>   *RECO_LLLl1_PT;
   vector<float>   *RECO_LLLL0ss_MASS;
   vector<float>   *RECO_LLLL0ss_PT;
   vector<float>   *RECO_LLLL1ss_MASS;
   vector<float>   *RECO_LLLL1ss_PT;
   vector<float>   *RECO_LLLL2ss_MASS;
   vector<float>   *RECO_LLLL2ss_PT;
   vector<float>   *RECO_LLLL_MASS;
   vector<float>   *RECO_LLLL_PT;
   vector<float>   *RECO_LLLL_ETA;
   vector<float>   *RECO_LLLL_PHI;
   vector<bool>    *RECOzMuMu_MatchingMCTruth;
   vector<float>   *RECOzMuMu_MatchingMCpT;
   vector<float>   *RECOzMuMu_MatchingMCmass;
   vector<float>   *RECOzMuMu_MatchingMCEta;
   vector<float>   *RECOzMuMu_MatchingMCPhi;
   vector<bool>    *RECOzEE_MatchingMCTruth;
   vector<float>   *RECOzEE_MatchingMCpT;
   vector<float>   *RECOzEE_MatchingMCmass;
   vector<float>   *RECOzEE_MatchingMCEta;
   vector<float>   *RECOzEE_MatchingMCPhi;
   vector<bool>    *RECOHzzMMMM_MatchingMCTruth;
   vector<float>   *RECOHzzMMMM_MatchingMCpT;
   vector<float>   *RECOHzzMMMM_MatchingMCmass;
   vector<float>   *RECOHzzMMMM_MatchingMCEta;
   vector<float>   *RECOHzzMMMM_MatchingMCPhi;
   vector<bool>    *RECOHzzEEEE_MatchingMCTruth;
   vector<float>   *RECOHzzEEEE_MatchingMCpT;
   vector<float>   *RECOHzzEEEE_MatchingMCmass;
   vector<float>   *RECOHzzEEEE_MatchingMCEta;
   vector<float>   *RECOHzzEEEE_MatchingMCPhi;
   vector<bool>    *RECOHzzEEMM_MatchingMCTruth;
   vector<float>   *RECOHzzEEMM_MatchingMCpT;
   vector<float>   *RECOHzzEEMM_MatchingMCmass;
   vector<float>   *RECOHzzEEMM_MatchingMCEta;
   vector<float>   *RECOHzzEEMM_MatchingMCPhi;
   Int_t           num_PU_vertices;
   Int_t           PU_BunchCrossing;
   vector<float>   *MC_weighting;
   vector<double>  *StdFitVertexX;
   vector<double>  *StdFitVertexY;
   vector<double>  *StdFitVertexZ;
   vector<double>  *StdFitVertexChi2r;
   vector<double>  *StdFitVertexProb;
   vector<double>  *StdFitVertexTrack_PT;
   vector<double>  *StdFitVertexTrack_ETA;
   vector<double>  *StdFitVertexTrack_PHI;
   vector<double>  *KinFitVertexX;
   vector<double>  *KinFitVertexY;
   vector<double>  *KinFitVertexZ;
   vector<double>  *KinFitVertexChi2r;
   vector<double>  *KinFitVertexProb;
   vector<double>  *StdFitVertexXMMMM;
   vector<double>  *StdFitVertexYMMMM;
   vector<double>  *StdFitVertexZMMMM;
   vector<double>  *StdFitVertexChi2rMMMM;
   vector<double>  *StdFitVertexProbMMMM;
   vector<double>  *StdFitVertexTrackMMMM_PT;
   vector<double>  *StdFitVertexTrackMMMM_ETA;
   vector<double>  *StdFitVertexTrackMMMM_PHI;
   vector<double>  *KinFitVertexXMMMM;
   vector<double>  *KinFitVertexYMMMM;
   vector<double>  *KinFitVertexZMMMM;
   vector<double>  *KinFitVertexChi2rMMMM;
   vector<double>  *KinFitVertexProbMMMM;
   vector<double>  *StdFitVertexXEEEE;
   vector<double>  *StdFitVertexYEEEE;
   vector<double>  *StdFitVertexZEEEE;
   vector<double>  *StdFitVertexChi2rEEEE;
   vector<double>  *StdFitVertexProbEEEE;
   vector<double>  *StdFitVertexTrackEEEE_PT;
   vector<double>  *StdFitVertexTrackEEEE_ETA;
   vector<double>  *StdFitVertexTrackEEEE_PHI;
   vector<double>  *KinFitVertexXEEEE;
   vector<double>  *KinFitVertexYEEEE;
   vector<double>  *KinFitVertexZEEEE;
   vector<double>  *KinFitVertexChi2rEEEE;
   vector<double>  *KinFitVertexProbEEEE;
   vector<double>  *StdFitVertexChi2rDiLep;
   vector<double>  *StdFitVertexProbDiLep;
   vector<double>  *StdFitVertexChi2rMMM;
   vector<double>  *StdFitVertexProbMMM;
   vector<double>  *StdFitVertexChi2rMME;
   vector<double>  *StdFitVertexProbMME;
   vector<double>  *StdFitVertexChi2rEEE;
   vector<double>  *StdFitVertexProbEEE;
   vector<double>  *StdFitVertexChi2rMEE;
   vector<double>  *StdFitVertexProbMEE;
   vector<double>  *ftsigma;
   vector<double>  *gdX;
   vector<double>  *gdY;
   vector<double>  *gdZ;
   vector<double>  *ftsigmalag;
   vector<double>  *gdlagX;
   vector<double>  *gdlagY;
   vector<double>  *gdlagZ;
   vector<double>  *gdlagProb;
   vector<double>  *gdlagNdof;
   vector<double>  *ftsigmaMMMM;
   vector<double>  *gdXMMMM;
   vector<double>  *gdYMMMM;
   vector<double>  *gdZMMMM;
   vector<double>  *ftsigmalagMMMM;
   vector<double>  *gdlagXMMMM;
   vector<double>  *gdlagYMMMM;
   vector<double>  *gdlagZMMMM;
   vector<double>  *gdlagProbMMMM;
   vector<double>  *gdlagNdofMMMM;
   vector<double>  *ftsigmaEEEE;
   vector<double>  *gdXEEEE;
   vector<double>  *gdYEEEE;
   vector<double>  *gdZEEEE;
   vector<double>  *ftsigmalagEEEE;
   vector<double>  *gdlagXEEEE;
   vector<double>  *gdlagYEEEE;
   vector<double>  *gdlagZEEEE;
   vector<double>  *gdlagProbEEEE;
   vector<double>  *gdlagNdofEEEE;
   vector<double>  *RECORF_2e2mu_cosTheta1_spin;
   vector<double>  *RECORF_2e2mu_cosTheta2_spin;
   vector<double>  *RECORF_2e2mu_cosThetaStar_spin;
   vector<double>  *RECORF_2e2mu_Phi_spin;
   vector<double>  *RECORF_2e2mu_Phi1_spin;
   vector<double>  *RECORF_2e2mu_Phi2_spin;
   vector<double>  *RECORF_2e2mu_phi1RF_spin;
   vector<double>  *RECORF_2e2mu_phi2RF_spin;
   vector<double>  *RECORF_2e2mu_MELA;
   vector<double>  *RECORF_4e_cosTheta1_spin;
   vector<double>  *RECORF_4e_cosTheta2_spin;
   vector<double>  *RECORF_4e_cosThetaStar_spin;
   vector<double>  *RECORF_4e_Phi_spin;
   vector<double>  *RECORF_4e_Phi1_spin;
   vector<double>  *RECORF_4e_Phi2_spin;
   vector<double>  *RECORF_4e_phi1RF_spin;
   vector<double>  *RECORF_4e_phi2RF_spin;
   vector<double>  *RECORF_4e_MELA;
   vector<double>  *RECORF_4mu_cosTheta1_spin;
   vector<double>  *RECORF_4mu_cosTheta2_spin;
   vector<double>  *RECORF_4mu_cosThetaStar_spin;
   vector<double>  *RECORF_4mu_Phi_spin;
   vector<double>  *RECORF_4mu_Phi1_spin;
   vector<double>  *RECORF_4mu_Phi2_spin;
   vector<double>  *RECORF_4mu_phi1RF_spin;
   vector<double>  *RECORF_4mu_phi2RF_spin;
   vector<double>  *RECORF_4mu_MELA;
   Double_t        BeamSpot_X;
   Double_t        BeamSpot_Y;
   Double_t        BeamSpot_Z;
   vector<int>     *nbPv;
   vector<int>     *Nbdof;
   vector<float>   *PositionRho;
   vector<float>   *PositionX;
   vector<float>   *PositionY;
   vector<float>   *PositionZ;
   Int_t           RECO_PFJET_N;
   vector<int>     *RECO_PFJET_CHARGE;
   vector<float>   *RECO_PFJET_ET;
   vector<float>   *RECO_PFJET_PT;
   vector<float>   *RECO_PFJET_ETA;
   vector<float>   *RECO_PFJET_PHI;
   vector<int>     *RECO_PFJET_PUID;
   vector<float>   *RECO_PFJET_PUID_MVA;
   Double_t        RHO_ele;
   Double_t        RHO_mu;

   // List of branches
   TBranch        *b_irun;   //!
   TBranch        *b_ievt;   //!
   TBranch        *b_ils;   //!
   TBranch        *b_Avginstlumi;   //!
   TBranch        *b_RECO_NVTX;   //!
   TBranch        *b_RECO_VERTEX_x;   //!
   TBranch        *b_RECO_VERTEX_y;   //!
   TBranch        *b_RECO_VERTEX_z;   //!
   TBranch        *b_RECO_VERTEX_ndof;   //!
   TBranch        *b_RECO_VERTEX_chi2;   //!
   TBranch        *b_RECO_VERTEX_ntracks;   //!
   TBranch        *b_RECO_VERTEXPROB;   //!
   TBranch        *b_RECO_VERTEX_isValid;   //!
   TBranch        *b_RECO_VERTEX_TRACK_PT;   //!
   TBranch        *b_RECO_NTRACK;   //!
   TBranch        *b_RECO_TRACK_PT;   //!
   TBranch        *b_RECO_TRACK_ETA;   //!
   TBranch        *b_RECO_TRACK_PHI;   //!
   TBranch        *b_RECO_TRACK_CHI2;   //!
   TBranch        *b_RECO_TRACK_CHI2RED;   //!
   TBranch        *b_RECO_TRACK_CHI2PROB;   //!
   TBranch        *b_RECO_TRACK_NHITS;   //!
   TBranch        *b_RECO_TRACK_DXY;   //!
   TBranch        *b_RECO_TRACK_DXYERR;   //!
   TBranch        *b_RECO_TRACK_DZ;   //!
   TBranch        *b_RECO_TRACK_DZERR;   //!
   TBranch        *b_MC_GENMET;   //!
   TBranch        *b_HLTPathsFired;   //!
   TBranch        *b_RECO_nMuHLTMatch;   //!
   TBranch        *b_RECO_nEleHLTMatch;   //!
   TBranch        *b_RECOMU_PT_MuHLTMatch;   //!
   TBranch        *b_RECOMU_ETA_MuHLTMatch;   //!
   TBranch        *b_RECOMU_PHI_MuHLTMatch;   //!
   TBranch        *b_RECOELE_PT_EleHLTMatch;   //!
   TBranch        *b_RECOELE_ETA_EleHLTMatch;   //!
   TBranch        *b_RECOELE_PHI_EleHLTMatch;   //!
   TBranch        *b_MC_LEPT_PT;   //!
   TBranch        *b_MC_LEPT_ETA;   //!
   TBranch        *b_MC_LEPT_PHI;   //!
   TBranch        *b_MC_LEPT_THETA;   //!
   TBranch        *b_MC_LEPT_PDGID;   //!
   TBranch        *b_MC_Z_PT;   //!
   TBranch        *b_MC_Z_ETA;   //!
   TBranch        *b_MC_Z_PHI;   //!
   TBranch        *b_MC_Z_THETA;   //!
   TBranch        *b_MC_Z_MASS;   //!
   TBranch        *b_MC_Z_PDGID;   //!
   TBranch        *b_MC_fourl_MASS;   //!
   TBranch        *b_MC_fourl_PT;   //!
   TBranch        *b_MC_fourl_PDGID;   //!
   TBranch        *b_MC_ZZ_MASS;   //!
   TBranch        *b_MC_ZZ_PT;   //!
   TBranch        *b_MC_ZZ_ETA;   //!
   TBranch        *b_MC_ZZ_PHI;   //!
   TBranch        *b_MC_ZZ_THETA;   //!
   TBranch        *b_MC_ZZ_PDGID;   //!
   TBranch        *b_MC_E;   //!
   TBranch        *b_MC_PT;   //!
   TBranch        *b_MC_ETA;   //!
   TBranch        *b_MC_THETA;   //!
   TBranch        *b_MC_PHI;   //!
   TBranch        *b_MC_MASS;   //!
   TBranch        *b_MC_PDGID;   //!
   TBranch        *b_RECO_NMU;   //!
   TBranch        *b_RECOMU_isPFMu;   //!
   TBranch        *b_RECOMU_isGlobalMu;   //!
   TBranch        *b_RECOMU_isStandAloneMu;   //!
   TBranch        *b_RECOMU_isTrackerMu;   //!
   TBranch        *b_RECOMU_isCaloMu;   //!
   TBranch        *b_RECOMU_E;   //!
   TBranch        *b_RECOMU_PT;   //!
   TBranch        *b_RECOMU_P;   //!
   TBranch        *b_RECOMU_ETA;   //!
   TBranch        *b_RECOMU_THETA;   //!
   TBranch        *b_RECOMU_PHI;   //!
   TBranch        *b_RECOMU_MASS;   //!
   TBranch        *b_RECOMU_CHARGE;   //!
   TBranch        *b_RECOMU_COV;   //!
   TBranch        *b_RECOMU_TRACKISO;   //!
   TBranch        *b_RECOMU_TRACKISO_SUMPT;   //!
   TBranch        *b_RECOMU_ECALISO;   //!
   TBranch        *b_RECOMU_HCALISO;   //!
   TBranch        *b_RECOMU_X;   //!
   TBranch        *b_RECOMU_PFchHad;   //!
   TBranch        *b_RECOMU_PFneuHad;   //!
   TBranch        *b_RECOMU_PFphoton;   //!
   TBranch        *b_RECOMU_PFPUchAllPart;   //!
   TBranch        *b_RECOMU_PFX_dB;   //!
   TBranch        *b_RECOMU_PFX_rho;   //!
   TBranch        *b_RECOMU_SIP;   //!
   TBranch        *b_RECOMU_IP;   //!
   TBranch        *b_RECOMU_IPERROR;   //!
   TBranch        *b_RECOMU_SIP_KF;   //!
   TBranch        *b_RECOMU_IP_KF;   //!
   TBranch        *b_RECOMU_IPERROR_KF;   //!
   TBranch        *b_RECOMU_STIP;   //!
   TBranch        *b_RECOMU_SLIP;   //!
   TBranch        *b_RECOMU_TIP;   //!
   TBranch        *b_RECOMU_LIP;   //!
   TBranch        *b_RECOMU_TIPERROR;   //!
   TBranch        *b_RECOMU_LIPERROR;   //!
   TBranch        *b_RECOMU_numberOfMatches;   //!
   TBranch        *b_RECOMU_numberOfMatchedStations;   //!
   TBranch        *b_RECOMU_caloCompatibility;   //!
   TBranch        *b_RECOMU_segmentCompatibility;   //!
   TBranch        *b_RECOMU_glbmuPromptTight;   //!
   TBranch        *b_RECOMU_mubesttrkType;   //!
   TBranch        *b_RECOMU_mubesttrkDxy;   //!
   TBranch        *b_RECOMU_mubesttrkDxyB;   //!
   TBranch        *b_RECOMU_mubesttrkDxyError;   //!
   TBranch        *b_RECOMU_mubesttrkDz;   //!
   TBranch        *b_RECOMU_mubesttrkDzB;   //!
   TBranch        *b_RECOMU_mubesttrkDzError;   //!
   TBranch        *b_RECOMU_mubesttrkPTError;   //!
   TBranch        *b_RECOMU_mutrkPT;   //!
   TBranch        *b_RECOMU_mutrkPTError;   //!
   TBranch        *b_RECOMU_mutrkDxy;   //!
   TBranch        *b_RECOMU_mutrkDxyError;   //!
   TBranch        *b_RECOMU_mutrkDxyB;   //!
   TBranch        *b_RECOMU_mutrkDz;   //!
   TBranch        *b_RECOMU_mutrkDzError;   //!
   TBranch        *b_RECOMU_mutrkDzB;   //!
   TBranch        *b_RECOMU_mutrkChi2PerNdof;   //!
   TBranch        *b_RECOMU_mutrkCharge;   //!
   TBranch        *b_RECOMU_mutrkNHits;   //!
   TBranch        *b_RECOMU_mutrkNPixHits;   //!
   TBranch        *b_RECOMU_mutrkNStripHits;   //!
   TBranch        *b_RECOMU_mutrkNMuonHits;   //!
   TBranch        *b_RECOMU_mutrktrackerLayersWithMeasurement;   //!
   TBranch        *b_RECOMU_muInnertrkDxy;   //!
   TBranch        *b_RECOMU_muInnertrkDxyError;   //!
   TBranch        *b_RECOMU_muInnertrkDxyB;   //!
   TBranch        *b_RECOMU_muInnertrkDz;   //!
   TBranch        *b_RECOMU_muInnertrkDzError;   //!
   TBranch        *b_RECOMU_muInnertrkDzB;   //!
   TBranch        *b_RECOMU_muInnertrkChi2PerNdof;   //!
   TBranch        *b_RECOMU_muInnertrktrackerLayersWithMeasurement;   //!
   TBranch        *b_RECOMU_muInnertrkPT;   //!
   TBranch        *b_RECOMU_muInnertrkPTError;   //!
   TBranch        *b_RECOMU_muInnertrkCharge;   //!
   TBranch        *b_RECOMU_muInnertrkNHits;   //!
   TBranch        *b_RECOMU_muInnertrkNPixHits;   //!
   TBranch        *b_RECOMU_muInnertrkNStripHits;   //!
   TBranch        *b_RECOMU_trkmuArbitration;   //!
   TBranch        *b_RECOMU_trkmu2DCompatibilityLoose;   //!
   TBranch        *b_RECOMU_trkmu2DCompatibilityTight;   //!
   TBranch        *b_RECOMU_trkmuOneStationLoose;   //!
   TBranch        *b_RECOMU_trkmuOneStationTight;   //!
   TBranch        *b_RECOMU_trkmuLastStationLoose;   //!
   TBranch        *b_RECOMU_trkmuLastStationTight;   //!
   TBranch        *b_RECOMU_trkmuOneStationAngLoose;   //!
   TBranch        *b_RECOMU_trkmuOneStationAngTight;   //!
   TBranch        *b_RECOMU_trkmuLastStationAngLoose;   //!
   TBranch        *b_RECOMU_trkmuLastStationAngTight;   //!
   TBranch        *b_RECOMU_trkmuLastStationOptimizedLowPtLoose;   //!
   TBranch        *b_RECOMU_trkmuLastStationOptimizedLowPtTight;   //!
   TBranch        *b_RECOMU_MatchingMCTruth;   //!
   TBranch        *b_RECOMU_MatchingMCpT;   //!
   TBranch        *b_RECOMU_MatchingMCEta;   //!
   TBranch        *b_RECOMU_MatchingMCPhi;   //!
   TBranch        *b_RECO_NELE;   //!
   TBranch        *b_RECOELE_isEcalDriven;   //!
   TBranch        *b_RECOELE_isTrackerDriven;   //!
   TBranch        *b_RECOELE_CHARGE;   //!
   TBranch        *b_RECOELE_E;   //!
   TBranch        *b_RECOELE_PT;   //!
   TBranch        *b_RECOELE_P;   //!
   TBranch        *b_RECOELE_ETA;   //!
   TBranch        *b_RECOELE_THETA;   //!
   TBranch        *b_RECOELE_PHI;   //!
   TBranch        *b_RECOELE_MASS;   //!
   TBranch        *b_ele_sclRawE;   //!
   TBranch        *b_RECOELE_scl_E;   //!
   TBranch        *b_RECOELE_scl_Et;   //!
   TBranch        *b_RECOELE_scl_Eta;   //!
   TBranch        *b_RECOELE_scl_Phi;   //!
   TBranch        *b_ele_sclX;   //!
   TBranch        *b_ele_sclY;   //!
   TBranch        *b_ele_sclZ;   //!
   TBranch        *b_RECOELE_PTError;   //!
   TBranch        *b_RECOELE_COV;   //!
   TBranch        *b_RECOELE_EGMECALISO;   //!
   TBranch        *b_RECOELE_EGMHCALISO;   //!
   TBranch        *b_RECOELE_EGMX;   //!
   TBranch        *b_RECOELE_PFchAllPart;   //!
   TBranch        *b_RECOELE_PFchHad;   //!
   TBranch        *b_RECOELE_PFneuHad;   //!
   TBranch        *b_RECOELE_PFphoton;   //!
   TBranch        *b_RECOELE_PFPUchAllPart;   //!
   TBranch        *b_RECOELE_PFX_dB;   //!
   TBranch        *b_RECOELE_PFX_rho;   //!
   TBranch        *b_RECOELE_SIP;   //!
   TBranch        *b_RECOELE_IP;   //!
   TBranch        *b_RECOELE_IPERROR;   //!
   TBranch        *b_RECOELE_SIP_KF;   //!
   TBranch        *b_RECOELE_IP_KF;   //!
   TBranch        *b_RECOELE_IPERROR_KF;   //!
   TBranch        *b_RECOELE_STIP;   //!
   TBranch        *b_RECOELE_SLIP;   //!
   TBranch        *b_RECOELE_TIP;   //!
   TBranch        *b_RECOELE_LIP;   //!
   TBranch        *b_RECOELE_TIPERROR;   //!
   TBranch        *b_RECOELE_LIPERROR;   //!
   TBranch        *b_RECOELE_gsftrack_NPixHits;   //!
   TBranch        *b_RECOELE_gsftrack_NStripHits;   //!
   TBranch        *b_RECOELE_gsftrack_chi2;   //!
   TBranch        *b_RECOELE_gsftrack_dxyB;   //!
   TBranch        *b_RECOELE_gsftrack_dxy;   //!
   TBranch        *b_RECOELE_gsftrack_dxyError;   //!
   TBranch        *b_RECOELE_gsftrack_dzB;   //!
   TBranch        *b_RECOELE_gsftrack_dz;   //!
   TBranch        *b_RECOELE_gsftrack_dzError;   //!
   TBranch        *b_RECOELE_gsftrack_losthits;   //!
   TBranch        *b_RECOELE_gsftrack_validhits;   //!
   TBranch        *b_RECOELE_gsftrack_expected_inner_hits;   //!
   TBranch        *b_RECOELE_ep;   //!
   TBranch        *b_RECOELE_eSeedp;   //!
   TBranch        *b_RECOELE_eSeedpout;   //!
   TBranch        *b_RECOELE_eElepout;   //!
   TBranch        *b_RECOELE_deltaEtaIn;   //!
   TBranch        *b_RECOELE_deltaEtaSeed;   //!
   TBranch        *b_RECOELE_deltaEtaEle;   //!
   TBranch        *b_RECOELE_deltaPhiIn;   //!
   TBranch        *b_RECOELE_deltaPhiSeed;   //!
   TBranch        *b_RECOELE_deltaPhiEle;   //!
   TBranch        *b_RECOELE_isbarrel;   //!
   TBranch        *b_RECOELE_isendcap;   //!
   TBranch        *b_RECOELE_isGap;   //!
   TBranch        *b_RECOELE_isEBetaGap;   //!
   TBranch        *b_RECOELE_isEBphiGap;   //!
   TBranch        *b_RECOELE_isEEdeeGap;   //!
   TBranch        *b_RECOELE_isEEringGap;   //!
   TBranch        *b_RECOELE_sigmaIetaIeta;   //!
   TBranch        *b_RECOELE_sigmaEtaEta;   //!
   TBranch        *b_RECOELE_e15;   //!
   TBranch        *b_RECOELE_e25max;   //!
   TBranch        *b_RECOELE_e55;   //!
   TBranch        *b_RECOELE_he;   //!
   TBranch        *b_RECOELE_fbrem;   //!
   TBranch        *b_RECOELE_nbrems;   //!
   TBranch        *b_RECOELE_Class;   //!
   TBranch        *b_RECOELE_fbrem_mean;   //!
   TBranch        *b_RECOELE_fbrem_mode;   //!
   TBranch        *b_RECOELE_ecalEnergy;   //!
   TBranch        *b_ele_seedSubdet2;   //!
   TBranch        *b_ele_seedDphi2;   //!
   TBranch        *b_ele_seedDrz2;   //!
   TBranch        *b_ele_seedSubdet1;   //!
   TBranch        *b_ele_seedDphi1;   //!
   TBranch        *b_ele_seedDrz1;   //!
   TBranch        *b_RECOELE_mvaTrigV0;   //!
   TBranch        *b_RECOELE_mvaNonTrigV0;   //!
   TBranch        *b_ConvMapDist;   //!
   TBranch        *b_ConvMapDcot;   //!
   TBranch        *b_RECOELE_MatchingMCTruth;   //!
   TBranch        *b_RECOELE_MatchingMCpT;   //!
   TBranch        *b_RECOELE_MatchingMCEta;   //!
   TBranch        *b_RECOELE_MatchingMCPhi;   //!
   TBranch        *b_RECO_NPHOT;   //!
   TBranch        *b_RECOPHOT_PT;   //!
   TBranch        *b_RECOPHOT_ETA;   //!
   TBranch        *b_RECOPHOT_PHI;   //!
   TBranch        *b_RECOPHOT_THETA;   //!
   TBranch        *b_RECO_NPFPHOT;   //!
   TBranch        *b_RECOPFPHOT_PT;   //!
   TBranch        *b_RECOPFPHOT_PTError;   //!
   TBranch        *b_RECOPFPHOT_ETA;   //!
   TBranch        *b_RECOPFPHOT_PHI;   //!
   TBranch        *b_RECOPFPHOT_THETA;   //!
   TBranch        *b_RECOPFPHOT_PFchAllPart;   //!
   TBranch        *b_RECOPFPHOT_PFchHad;   //!
   TBranch        *b_RECOPFPHOT_PFneuHad;   //!
   TBranch        *b_RECOPFPHOT_PFphoton;   //!
   TBranch        *b_RECOPFPHOT_PFPUchAllPart;   //!
   TBranch        *b_RECOPFPHOT_PFX_rho;   //!
   TBranch        *b_RECOPHOT_MatchingMCTruth;   //!
   TBranch        *b_RECOPHOT_MatchingMCpT;   //!
   TBranch        *b_RECOPHOT_MatchingMCEta;   //!
   TBranch        *b_RECOPHOT_MatchingMCPhi;   //!
   TBranch        *b_PFMet_pt;   //!
   TBranch        *b_PFMet_eta;   //!
   TBranch        *b_PFMet_phi;   //!
   TBranch        *b_PFMet_en;   //!
   TBranch        *b_PFMet_px;   //!
   TBranch        *b_PFMet_py;   //!
   TBranch        *b_PFMet_pz;   //!
   TBranch        *b_PFMet_sumEt;   //!
   TBranch        *b_METSign;   //!
   TBranch        *b_tCHighEff_nb;   //!
   TBranch        *b_tCHighEff_BTagJet_PT;   //!
   TBranch        *b_tCHighEff_BTagJet_ETA;   //!
   TBranch        *b_tCHighEff_BTagJet_PHI;   //!
   TBranch        *b_tCHighEff_BTagJet_DISCR;   //!
   TBranch        *b_tCHighPur_nb;   //!
   TBranch        *b_tCHighPur_BTagJet_PT;   //!
   TBranch        *b_tCHighPur_BTagJet_ETA;   //!
   TBranch        *b_tCHighPur_BTagJet_PHI;   //!
   TBranch        *b_tCHighPur_BTagJet_DISCR;   //!
   TBranch        *b_cSV_nb;   //!
   TBranch        *b_cSV_BTagJet_PT;   //!
   TBranch        *b_cSV_BTagJet_ETA;   //!
   TBranch        *b_cSV_BTagJet_PHI;   //!
   TBranch        *b_cSV_BTagJet_DISCR;   //!
   TBranch        *b_cSV_BTagJet_ET;   //!
   TBranch        *b_RECO_ZMM_MASS;   //!
   TBranch        *b_RECO_ZEE_MASS;   //!
   TBranch        *b_RECO_DiLep_MASS;   //!
   TBranch        *b_RECO_ZMM_PT;   //!
   TBranch        *b_RECO_ZEE_PT;   //!
   TBranch        *b_RECO_DiLep_PT;   //!
   TBranch        *b_RECO_ZMM_ETA;   //!
   TBranch        *b_RECO_ZEE_ETA;   //!
   TBranch        *b_RECO_DiLep_ETA;   //!
   TBranch        *b_RECO_ZMM_PHI;   //!
   TBranch        *b_RECO_ZEE_PHI;   //!
   TBranch        *b_RECO_DiLep_PHI;   //!
   TBranch        *b_RECO_ZMMss_MASS;   //!
   TBranch        *b_RECO_ZEEss_MASS;   //!
   TBranch        *b_RECO_ZEM_MASS;   //!
   TBranch        *b_RECO_ZMMss_PT;   //!
   TBranch        *b_RECO_ZEEss_PT;   //!
   TBranch        *b_RECO_ZEM_PT;   //!
   TBranch        *b_RECO_ZMMss_ETA;   //!
   TBranch        *b_RECO_ZEEss_ETA;   //!
   TBranch        *b_RECO_ZEM_ETA;   //!
   TBranch        *b_RECO_ZMMss_PHI;   //!
   TBranch        *b_RECO_ZEEss_PHI;   //!
   TBranch        *b_RECO_ZEM_PHI;   //!
   TBranch        *b_RECO_MMMM_MASS;   //!
   TBranch        *b_RECO_MMMM_PT;   //!
   TBranch        *b_RECO_MMMM_ETA;   //!
   TBranch        *b_RECO_MMMM_PHI;   //!
   TBranch        *b_RECO_MMMM_MASS_REFIT;   //!
   TBranch        *b_RECO_EEEE_MASS;   //!
   TBranch        *b_RECO_EEEE_PT;   //!
   TBranch        *b_RECO_EEEE_ETA;   //!
   TBranch        *b_RECO_EEEE_PHI;   //!
   TBranch        *b_RECO_EEEE_MASS_REFIT;   //!
   TBranch        *b_RECO_EEMM_MASS;   //!
   TBranch        *b_RECO_EEMM_PT;   //!
   TBranch        *b_RECO_EEMM_ETA;   //!
   TBranch        *b_RECO_EEMM_PHI;   //!
   TBranch        *b_RECO_EEMM_MASS_REFIT;   //!
   TBranch        *b_RECO_LLL0_MASS;   //!
   TBranch        *b_RECO_LLL1_MASS;   //!
   TBranch        *b_RECO_LLL2_MASS;   //!
   TBranch        *b_RECO_LLL3_MASS;   //!
   TBranch        *b_RECO_LLL0_PT;   //!
   TBranch        *b_RECO_LLL1_PT;   //!
   TBranch        *b_RECO_LLL2_PT;   //!
   TBranch        *b_RECO_LLL3_PT;   //!
   TBranch        *b_RECO_LLLl0_MASS;   //!
   TBranch        *b_RECO_LLLl1_MASS;   //!
   TBranch        *b_RECO_LLLl0_PT;   //!
   TBranch        *b_RECO_LLLl1_PT;   //!
   TBranch        *b_RECO_LLLL0ss_MASS;   //!
   TBranch        *b_RECO_LLLL0ss_PT;   //!
   TBranch        *b_RECO_LLLL1ss_MASS;   //!
   TBranch        *b_RECO_LLLL1ss_PT;   //!
   TBranch        *b_RECO_LLLL2ss_MASS;   //!
   TBranch        *b_RECO_LLLL2ss_PT;   //!
   TBranch        *b_RECO_LLLL_MASS;   //!
   TBranch        *b_RECO_LLLL_PT;   //!
   TBranch        *b_RECO_LLLL_ETA;   //!
   TBranch        *b_RECO_LLLL_PHI;   //!
   TBranch        *b_RECOzMuMu_MatchingMCTruth;   //!
   TBranch        *b_RECOzMuMu_MatchingMCpT;   //!
   TBranch        *b_RECOzMuMu_MatchingMCmass;   //!
   TBranch        *b_RECOzMuMu_MatchingMCEta;   //!
   TBranch        *b_RECOzMuMu_MatchingMCPhi;   //!
   TBranch        *b_RECOzEE_MatchingMCTruth;   //!
   TBranch        *b_RECOzEE_MatchingMCpT;   //!
   TBranch        *b_RECOzEE_MatchingMCmass;   //!
   TBranch        *b_RECOzEE_MatchingMCEta;   //!
   TBranch        *b_RECOzEE_MatchingMCPhi;   //!
   TBranch        *b_RECOHzzMMMM_MatchingMCTruth;   //!
   TBranch        *b_RECOHzzMMMM_MatchingMCpT;   //!
   TBranch        *b_RECOHzzMMMM_MatchingMCmass;   //!
   TBranch        *b_RECOHzzMMMM_MatchingMCEta;   //!
   TBranch        *b_RECOHzzMMMM_MatchingMCPhi;   //!
   TBranch        *b_RECOHzzEEEE_MatchingMCTruth;   //!
   TBranch        *b_RECOHzzEEEE_MatchingMCpT;   //!
   TBranch        *b_RECOHzzEEEE_MatchingMCmass;   //!
   TBranch        *b_RECOHzzEEEE_MatchingMCEta;   //!
   TBranch        *b_RECOHzzEEEE_MatchingMCPhi;   //!
   TBranch        *b_RECOHzzEEMM_MatchingMCTruth;   //!
   TBranch        *b_RECOHzzEEMM_MatchingMCpT;   //!
   TBranch        *b_RECOHzzEEMM_MatchingMCmass;   //!
   TBranch        *b_RECOHzzEEMM_MatchingMCEta;   //!
   TBranch        *b_RECOHzzEEMM_MatchingMCPhi;   //!
   TBranch        *b_num_PU_vertices;   //!
   TBranch        *b_PU_BunchCrossing;   //!
   TBranch        *b_MC_weighting;   //!
   TBranch        *b_StdFitVertexX;   //!
   TBranch        *b_StdFitVertexY;   //!
   TBranch        *b_StdFitVertexZ;   //!
   TBranch        *b_StdFitVertexChi2r;   //!
   TBranch        *b_StdFitVertexProb;   //!
   TBranch        *b_StdFitVertexTrack_PT;   //!
   TBranch        *b_StdFitVertexTrack_ETA;   //!
   TBranch        *b_StdFitVertexTrack_PHI;   //!
   TBranch        *b_KinFitVertexX;   //!
   TBranch        *b_KinFitVertexY;   //!
   TBranch        *b_KinFitVertexZ;   //!
   TBranch        *b_KinFitVertexChi2r;   //!
   TBranch        *b_KinFitVertexProb;   //!
   TBranch        *b_StdFitVertexXMMMM;   //!
   TBranch        *b_StdFitVertexYMMMM;   //!
   TBranch        *b_StdFitVertexZMMMM;   //!
   TBranch        *b_StdFitVertexChi2rMMMM;   //!
   TBranch        *b_StdFitVertexProbMMMM;   //!
   TBranch        *b_StdFitVertexTrackMMMM_PT;   //!
   TBranch        *b_StdFitVertexTrackMMMM_ETA;   //!
   TBranch        *b_StdFitVertexTrackMMMM_PHI;   //!
   TBranch        *b_KinFitVertexXMMMM;   //!
   TBranch        *b_KinFitVertexYMMMM;   //!
   TBranch        *b_KinFitVertexZMMMM;   //!
   TBranch        *b_KinFitVertexChi2rMMMM;   //!
   TBranch        *b_KinFitVertexProbMMMM;   //!
   TBranch        *b_StdFitVertexXEEEE;   //!
   TBranch        *b_StdFitVertexYEEEE;   //!
   TBranch        *b_StdFitVertexZEEEE;   //!
   TBranch        *b_StdFitVertexChi2rEEEE;   //!
   TBranch        *b_StdFitVertexProbEEEE;   //!
   TBranch        *b_StdFitVertexTrackEEEE_PT;   //!
   TBranch        *b_StdFitVertexTrackEEEE_ETA;   //!
   TBranch        *b_StdFitVertexTrackEEEE_PHI;   //!
   TBranch        *b_KinFitVertexXEEEE;   //!
   TBranch        *b_KinFitVertexYEEEE;   //!
   TBranch        *b_KinFitVertexZEEEE;   //!
   TBranch        *b_KinFitVertexChi2rEEEE;   //!
   TBranch        *b_KinFitVertexProbEEEE;   //!
   TBranch        *b_StdFitVertexChi2rDiLep;   //!
   TBranch        *b_StdFitVertexProbDiLep;   //!
   TBranch        *b_StdFitVertexChi2rMMM;   //!
   TBranch        *b_StdFitVertexProbMMM;   //!
   TBranch        *b_StdFitVertexChi2rMME;   //!
   TBranch        *b_StdFitVertexProbMME;   //!
   TBranch        *b_StdFitVertexChi2rEEE;   //!
   TBranch        *b_StdFitVertexProbEEE;   //!
   TBranch        *b_StdFitVertexChi2rMEE;   //!
   TBranch        *b_StdFitVertexProbMEE;   //!
   TBranch        *b_ftsigma;   //!
   TBranch        *b_gdX;   //!
   TBranch        *b_gdY;   //!
   TBranch        *b_gdZ;   //!
   TBranch        *b_ftsigmalag;   //!
   TBranch        *b_gdlagX;   //!
   TBranch        *b_gdlagY;   //!
   TBranch        *b_gdlagZ;   //!
   TBranch        *b_gdlagProb;   //!
   TBranch        *b_gdlagNdof;   //!
   TBranch        *b_ftsigmaMMMM;   //!
   TBranch        *b_gdXMMMM;   //!
   TBranch        *b_gdYMMMM;   //!
   TBranch        *b_gdZMMMM;   //!
   TBranch        *b_ftsigmalagMMMM;   //!
   TBranch        *b_gdlagXMMMM;   //!
   TBranch        *b_gdlagYMMMM;   //!
   TBranch        *b_gdlagZMMMM;   //!
   TBranch        *b_gdlagProbMMMM;   //!
   TBranch        *b_gdlagNdofMMMM;   //!
   TBranch        *b_ftsigmaEEEE;   //!
   TBranch        *b_gdXEEEE;   //!
   TBranch        *b_gdYEEEE;   //!
   TBranch        *b_gdZEEEE;   //!
   TBranch        *b_ftsigmalagEEEE;   //!
   TBranch        *b_gdlagXEEEE;   //!
   TBranch        *b_gdlagYEEEE;   //!
   TBranch        *b_gdlagZEEEE;   //!
   TBranch        *b_gdlagProbEEEE;   //!
   TBranch        *b_gdlagNdofEEEE;   //!
   TBranch        *b_RECORF_2e2mu_cosTheta1_spin;   //!
   TBranch        *b_RECORF_2e2mu_cosTheta2_spin;   //!
   TBranch        *b_RECORF_2e2mu_cosThetaStar_spin;   //!
   TBranch        *b_RECORF_2e2mu_Phi_spin;   //!
   TBranch        *b_RECORF_2e2mu_Phi1_spin;   //!
   TBranch        *b_RECORF_2e2mu_Phi2_spin;   //!
   TBranch        *b_RECORF_2e2mu_phi1RF_spin;   //!
   TBranch        *b_RECORF_2e2mu_phi2RF_spin;   //!
   TBranch        *b_RECORF_2e2mu_MELA;   //!
   TBranch        *b_RECORF_4e_cosTheta1_spin;   //!
   TBranch        *b_RECORF_4e_cosTheta2_spin;   //!
   TBranch        *b_RECORF_4e_cosThetaStar_spin;   //!
   TBranch        *b_RECORF_4e_Phi_spin;   //!
   TBranch        *b_RECORF_4e_Phi1_spin;   //!
   TBranch        *b_RECORF_4e_Phi2_spin;   //!
   TBranch        *b_RECORF_4e_phi1RF_spin;   //!
   TBranch        *b_RECORF_4e_phi2RF_spin;   //!
   TBranch        *b_RECORF_4e_MELA;   //!
   TBranch        *b_RECORF_4mu_cosTheta1_spin;   //!
   TBranch        *b_RECORF_4mu_cosTheta2_spin;   //!
   TBranch        *b_RECORF_4mu_cosThetaStar_spin;   //!
   TBranch        *b_RECORF_4mu_Phi_spin;   //!
   TBranch        *b_RECORF_4mu_Phi1_spin;   //!
   TBranch        *b_RECORF_4mu_Phi2_spin;   //!
   TBranch        *b_RECORF_4mu_phi1RF_spin;   //!
   TBranch        *b_RECORF_4mu_phi2RF_spin;   //!
   TBranch        *b_RECORF_4mu_MELA;   //!
   TBranch        *b_BeamSpot_X;   //!
   TBranch        *b_BeamSpot_Y;   //!
   TBranch        *b_BeamSpot_Z;   //!
   TBranch        *b_nbPv;   //!
   TBranch        *b_Nbdof;   //!
   TBranch        *b_PositionRho;   //!
   TBranch        *b_PositionX;   //!
   TBranch        *b_PositionY;   //!
   TBranch        *b_PositionZ;   //!
   TBranch        *b_RECO_PFJET_N;   //!
   TBranch        *b_RECO_PFJET_CHARGE;   //!
   TBranch        *b_RECO_PFJET_ET;   //!
   TBranch        *b_RECO_PFJET_PT;   //!
   TBranch        *b_RECO_PFJET_ETA;   //!
   TBranch        *b_RECO_PFJET_PHI;   //!
   TBranch        *b_RECO_PFJET_PUID;   //!
   TBranch        *b_RECO_PFJET_PUID_MVA;   //!
   TBranch        *b_RHO_ele;   //!
   TBranch        *b_RHO_mu;   //!

   MonoHiggsAnalysis4mu(TTree *tree=0,Double_t weight_=1.,string DATA_type_="DATA",string MC_type_="MC");
   virtual ~MonoHiggsAnalysis4mu();
   Double_t weight;
   string DATA_type,MC_type;
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Char_t *name);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   ofstream output_txt; 
   ofstream bnn_file;
   double DELTAPHI( double phi1, double phi2 );
   double EAele(int ,bool );
   double masserror( std::vector<TLorentzVector> Lep, std::vector<double> pterr );
   void initpuweights2011();
   
   double invmass (float M1, float PT1, float ETA1, float PHI1, float M2, float PT2, float ETA2, float PHI2 );
   void printmubnn(int i);
};

#endif

#ifdef MonoHiggsAnalysis4mu_cxx
MonoHiggsAnalysis4mu::MonoHiggsAnalysis4mu(TTree *tree,Double_t weight_, string DATA_type_, string MC_type_) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  weight = weight_;
  DATA_type = DATA_type_;
  MC_type = MC_type_;
  
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("roottree_leptons_sync_Fall15_HiggsToZZ_76x_vector.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("roottree_leptons_sync_Fall15_HiggsToZZ_76x_vector.root");
      }
      f->GetObject("HZZ4LeptonsAnalysis",tree);

   }
   Init(tree);
}

MonoHiggsAnalysis4mu::~MonoHiggsAnalysis4mu()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MonoHiggsAnalysis4mu::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MonoHiggsAnalysis4mu::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MonoHiggsAnalysis4mu::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   RECO_VERTEX_x = 0;
   RECO_VERTEX_y = 0;
   RECO_VERTEX_z = 0;
   RECO_VERTEX_ndof = 0;
   RECO_VERTEX_chi2 = 0;
   RECO_VERTEX_ntracks = 0;
   RECO_VERTEXPROB = 0;
   RECO_VERTEX_isValid = 0;
   RECO_VERTEX_TRACK_PT = 0;
   RECO_TRACK_PT = 0;
   RECO_TRACK_ETA = 0;
   RECO_TRACK_PHI = 0;
   RECO_TRACK_CHI2 = 0;
   RECO_TRACK_CHI2RED = 0;
   RECO_TRACK_CHI2PROB = 0;
   RECO_TRACK_NHITS = 0;
   RECO_TRACK_DXY = 0;
   RECO_TRACK_DXYERR = 0;
   RECO_TRACK_DZ = 0;
   RECO_TRACK_DZERR = 0;
   RECOMU_PT_MuHLTMatch = 0;
   RECOMU_ETA_MuHLTMatch = 0;
   RECOMU_PHI_MuHLTMatch = 0;
   RECOELE_PT_EleHLTMatch = 0;
   RECOELE_ETA_EleHLTMatch = 0;
   RECOELE_PHI_EleHLTMatch = 0;
   MC_LEPT_PT = 0;
   MC_LEPT_ETA = 0;
   MC_LEPT_PHI = 0;
   MC_LEPT_THETA = 0;
   MC_LEPT_PDGID = 0;
   MC_Z_PT = 0;
   MC_Z_ETA = 0;
   MC_Z_PHI = 0;
   MC_Z_THETA = 0;
   MC_Z_MASS = 0;
   MC_Z_PDGID = 0;
   MC_fourl_MASS = 0;
   MC_fourl_PT = 0;
   MC_fourl_PDGID = 0;
   MC_ZZ_MASS = 0;
   MC_ZZ_PT = 0;
   MC_ZZ_ETA = 0;
   MC_ZZ_PHI = 0;
   MC_ZZ_THETA = 0;
   MC_ZZ_PDGID = 0;
   MC_E = 0;
   MC_PT = 0;
   MC_ETA = 0;
   MC_THETA = 0;
   MC_PHI = 0;
   MC_MASS = 0;
   MC_PDGID = 0;
   RECOMU_isPFMu = 0;
   RECOMU_isGlobalMu = 0;
   RECOMU_isStandAloneMu = 0;
   RECOMU_isTrackerMu = 0;
   RECOMU_isCaloMu = 0;
   RECOMU_E = 0;
   RECOMU_PT = 0;
   RECOMU_P = 0;
   RECOMU_ETA = 0;
   RECOMU_THETA = 0;
   RECOMU_PHI = 0;
   RECOMU_MASS = 0;
   RECOMU_CHARGE = 0;
   RECOMU_COV = 0;
   RECOMU_TRACKISO = 0;
   RECOMU_TRACKISO_SUMPT = 0;
   RECOMU_ECALISO = 0;
   RECOMU_HCALISO = 0;
   RECOMU_X = 0;
   RECOMU_PFchHad = 0;
   RECOMU_PFneuHad = 0;
   RECOMU_PFphoton = 0;
   RECOMU_PFPUchAllPart = 0;
   RECOMU_PFX_dB = 0;
   RECOMU_PFX_rho = 0;
   RECOMU_SIP = 0;
   RECOMU_IP = 0;
   RECOMU_IPERROR = 0;
   RECOMU_SIP_KF = 0;
   RECOMU_IP_KF = 0;
   RECOMU_IPERROR_KF = 0;
   RECOMU_STIP = 0;
   RECOMU_SLIP = 0;
   RECOMU_TIP = 0;
   RECOMU_LIP = 0;
   RECOMU_TIPERROR = 0;
   RECOMU_LIPERROR = 0;
   RECOMU_numberOfMatches = 0;
   RECOMU_numberOfMatchedStations = 0;
   RECOMU_caloCompatibility = 0;
   RECOMU_segmentCompatibility = 0;
   RECOMU_glbmuPromptTight = 0;
   RECOMU_mubesttrkType = 0;
   RECOMU_mubesttrkDxy = 0;
   RECOMU_mubesttrkDxyB = 0;
   RECOMU_mubesttrkDxyError = 0;
   RECOMU_mubesttrkDz = 0;
   RECOMU_mubesttrkDzB = 0;
   RECOMU_mubesttrkDzError = 0;
   RECOMU_mubesttrkPTError = 0;
   RECOMU_mutrkPT = 0;
   RECOMU_mutrkPTError = 0;
   RECOMU_mutrkDxy = 0;
   RECOMU_mutrkDxyError = 0;
   RECOMU_mutrkDxyB = 0;
   RECOMU_mutrkDz = 0;
   RECOMU_mutrkDzError = 0;
   RECOMU_mutrkDzB = 0;
   RECOMU_mutrkChi2PerNdof = 0;
   RECOMU_mutrkCharge = 0;
   RECOMU_mutrkNHits = 0;
   RECOMU_mutrkNPixHits = 0;
   RECOMU_mutrkNStripHits = 0;
   RECOMU_mutrkNMuonHits = 0;
   RECOMU_mutrktrackerLayersWithMeasurement = 0;
   RECOMU_muInnertrkDxy = 0;
   RECOMU_muInnertrkDxyError = 0;
   RECOMU_muInnertrkDxyB = 0;
   RECOMU_muInnertrkDz = 0;
   RECOMU_muInnertrkDzError = 0;
   RECOMU_muInnertrkDzB = 0;
   RECOMU_muInnertrkChi2PerNdof = 0;
   RECOMU_muInnertrktrackerLayersWithMeasurement = 0;
   RECOMU_muInnertrkPT = 0;
   RECOMU_muInnertrkPTError = 0;
   RECOMU_muInnertrkCharge = 0;
   RECOMU_muInnertrkNHits = 0;
   RECOMU_muInnertrkNPixHits = 0;
   RECOMU_muInnertrkNStripHits = 0;
   RECOMU_trkmuArbitration = 0;
   RECOMU_trkmu2DCompatibilityLoose = 0;
   RECOMU_trkmu2DCompatibilityTight = 0;
   RECOMU_trkmuOneStationLoose = 0;
   RECOMU_trkmuOneStationTight = 0;
   RECOMU_trkmuLastStationLoose = 0;
   RECOMU_trkmuLastStationTight = 0;
   RECOMU_trkmuOneStationAngLoose = 0;
   RECOMU_trkmuOneStationAngTight = 0;
   RECOMU_trkmuLastStationAngLoose = 0;
   RECOMU_trkmuLastStationAngTight = 0;
   RECOMU_trkmuLastStationOptimizedLowPtLoose = 0;
   RECOMU_trkmuLastStationOptimizedLowPtTight = 0;
   RECOMU_MatchingMCTruth = 0;
   RECOMU_MatchingMCpT = 0;
   RECOMU_MatchingMCEta = 0;
   RECOMU_MatchingMCPhi = 0;
   RECOELE_isEcalDriven = 0;
   RECOELE_isTrackerDriven = 0;
   RECOELE_CHARGE = 0;
   RECOELE_E = 0;
   RECOELE_PT = 0;
   RECOELE_P = 0;
   RECOELE_ETA = 0;
   RECOELE_THETA = 0;
   RECOELE_PHI = 0;
   RECOELE_MASS = 0;
   ele_sclRawE = 0;
   RECOELE_scl_E = 0;
   RECOELE_scl_Et = 0;
   RECOELE_scl_Eta = 0;
   RECOELE_scl_Phi = 0;
   ele_sclX = 0;
   ele_sclY = 0;
   ele_sclZ = 0;
   RECOELE_PTError = 0;
   RECOELE_COV = 0;
   RECOELE_EGMECALISO = 0;
   RECOELE_EGMHCALISO = 0;
   RECOELE_EGMX = 0;
   RECOELE_PFchAllPart = 0;
   RECOELE_PFchHad = 0;
   RECOELE_PFneuHad = 0;
   RECOELE_PFphoton = 0;
   RECOELE_PFPUchAllPart = 0;
   RECOELE_PFX_dB = 0;
   RECOELE_PFX_rho = 0;
   RECOELE_SIP = 0;
   RECOELE_IP = 0;
   RECOELE_IPERROR = 0;
   RECOELE_SIP_KF = 0;
   RECOELE_IP_KF = 0;
   RECOELE_IPERROR_KF = 0;
   RECOELE_STIP = 0;
   RECOELE_SLIP = 0;
   RECOELE_TIP = 0;
   RECOELE_LIP = 0;
   RECOELE_TIPERROR = 0;
   RECOELE_LIPERROR = 0;
   RECOELE_gsftrack_NPixHits = 0;
   RECOELE_gsftrack_NStripHits = 0;
   RECOELE_gsftrack_chi2 = 0;
   RECOELE_gsftrack_dxyB = 0;
   RECOELE_gsftrack_dxy = 0;
   RECOELE_gsftrack_dxyError = 0;
   RECOELE_gsftrack_dzB = 0;
   RECOELE_gsftrack_dz = 0;
   RECOELE_gsftrack_dzError = 0;
   RECOELE_gsftrack_losthits = 0;
   RECOELE_gsftrack_validhits = 0;
   RECOELE_gsftrack_expected_inner_hits = 0;
   RECOELE_ep = 0;
   RECOELE_eSeedp = 0;
   RECOELE_eSeedpout = 0;
   RECOELE_eElepout = 0;
   RECOELE_deltaEtaIn = 0;
   RECOELE_deltaEtaSeed = 0;
   RECOELE_deltaEtaEle = 0;
   RECOELE_deltaPhiIn = 0;
   RECOELE_deltaPhiSeed = 0;
   RECOELE_deltaPhiEle = 0;
   RECOELE_isbarrel = 0;
   RECOELE_isendcap = 0;
   RECOELE_isGap = 0;
   RECOELE_isEBetaGap = 0;
   RECOELE_isEBphiGap = 0;
   RECOELE_isEEdeeGap = 0;
   RECOELE_isEEringGap = 0;
   RECOELE_sigmaIetaIeta = 0;
   RECOELE_sigmaEtaEta = 0;
   RECOELE_e15 = 0;
   RECOELE_e25max = 0;
   RECOELE_e55 = 0;
   RECOELE_he = 0;
   RECOELE_fbrem = 0;
   RECOELE_nbrems = 0;
   RECOELE_Class = 0;
   RECOELE_fbrem_mean = 0;
   RECOELE_fbrem_mode = 0;
   RECOELE_ecalEnergy = 0;
   ele_seedSubdet2 = 0;
   ele_seedDphi2 = 0;
   ele_seedDrz2 = 0;
   ele_seedSubdet1 = 0;
   ele_seedDphi1 = 0;
   ele_seedDrz1 = 0;
   RECOELE_mvaTrigV0 = 0;
   RECOELE_mvaNonTrigV0 = 0;
   ConvMapDist = 0;
   ConvMapDcot = 0;
   RECOELE_MatchingMCTruth = 0;
   RECOELE_MatchingMCpT = 0;
   RECOELE_MatchingMCEta = 0;
   RECOELE_MatchingMCPhi = 0;
   RECOPHOT_PT = 0;
   RECOPHOT_ETA = 0;
   RECOPHOT_PHI = 0;
   RECOPHOT_THETA = 0;
   RECOPFPHOT_PT = 0;
   RECOPFPHOT_PTError = 0;
   RECOPFPHOT_ETA = 0;
   RECOPFPHOT_PHI = 0;
   RECOPFPHOT_THETA = 0;
   RECOPFPHOT_PFchAllPart = 0;
   RECOPFPHOT_PFchHad = 0;
   RECOPFPHOT_PFneuHad = 0;
   RECOPFPHOT_PFphoton = 0;
   RECOPFPHOT_PFPUchAllPart = 0;
   RECOPFPHOT_PFX_rho = 0;
   RECOPHOT_MatchingMCTruth = 0;
   RECOPHOT_MatchingMCpT = 0;
   RECOPHOT_MatchingMCEta = 0;
   RECOPHOT_MatchingMCPhi = 0;
   tCHighEff_BTagJet_PT = 0;
   tCHighEff_BTagJet_ETA = 0;
   tCHighEff_BTagJet_PHI = 0;
   tCHighEff_BTagJet_DISCR = 0;
   tCHighPur_BTagJet_PT = 0;
   tCHighPur_BTagJet_ETA = 0;
   tCHighPur_BTagJet_PHI = 0;
   tCHighPur_BTagJet_DISCR = 0;
   cSV_BTagJet_PT = 0;
   cSV_BTagJet_ETA = 0;
   cSV_BTagJet_PHI = 0;
   cSV_BTagJet_DISCR = 0;
   cSV_BTagJet_ET = 0;
   RECO_ZMM_MASS = 0;
   RECO_ZEE_MASS = 0;
   RECO_DiLep_MASS = 0;
   RECO_ZMM_PT = 0;
   RECO_ZEE_PT = 0;
   RECO_DiLep_PT = 0;
   RECO_ZMM_ETA = 0;
   RECO_ZEE_ETA = 0;
   RECO_DiLep_ETA = 0;
   RECO_ZMM_PHI = 0;
   RECO_ZEE_PHI = 0;
   RECO_DiLep_PHI = 0;
   RECO_ZMMss_MASS = 0;
   RECO_ZEEss_MASS = 0;
   RECO_ZEM_MASS = 0;
   RECO_ZMMss_PT = 0;
   RECO_ZEEss_PT = 0;
   RECO_ZEM_PT = 0;
   RECO_ZMMss_ETA = 0;
   RECO_ZEEss_ETA = 0;
   RECO_ZEM_ETA = 0;
   RECO_ZMMss_PHI = 0;
   RECO_ZEEss_PHI = 0;
   RECO_ZEM_PHI = 0;
   RECO_MMMM_MASS = 0;
   RECO_MMMM_PT = 0;
   RECO_MMMM_ETA = 0;
   RECO_MMMM_PHI = 0;
   RECO_MMMM_MASS_REFIT = 0;
   RECO_EEEE_MASS = 0;
   RECO_EEEE_PT = 0;
   RECO_EEEE_ETA = 0;
   RECO_EEEE_PHI = 0;
   RECO_EEEE_MASS_REFIT = 0;
   RECO_EEMM_MASS = 0;
   RECO_EEMM_PT = 0;
   RECO_EEMM_ETA = 0;
   RECO_EEMM_PHI = 0;
   RECO_EEMM_MASS_REFIT = 0;
   RECO_LLL0_MASS = 0;
   RECO_LLL1_MASS = 0;
   RECO_LLL2_MASS = 0;
   RECO_LLL3_MASS = 0;
   RECO_LLL0_PT = 0;
   RECO_LLL1_PT = 0;
   RECO_LLL2_PT = 0;
   RECO_LLL3_PT = 0;
   RECO_LLLl0_MASS = 0;
   RECO_LLLl1_MASS = 0;
   RECO_LLLl0_PT = 0;
   RECO_LLLl1_PT = 0;
   RECO_LLLL0ss_MASS = 0;
   RECO_LLLL0ss_PT = 0;
   RECO_LLLL1ss_MASS = 0;
   RECO_LLLL1ss_PT = 0;
   RECO_LLLL2ss_MASS = 0;
   RECO_LLLL2ss_PT = 0;
   RECO_LLLL_MASS = 0;
   RECO_LLLL_PT = 0;
   RECO_LLLL_ETA = 0;
   RECO_LLLL_PHI = 0;
   RECOzMuMu_MatchingMCTruth = 0;
   RECOzMuMu_MatchingMCpT = 0;
   RECOzMuMu_MatchingMCmass = 0;
   RECOzMuMu_MatchingMCEta = 0;
   RECOzMuMu_MatchingMCPhi = 0;
   RECOzEE_MatchingMCTruth = 0;
   RECOzEE_MatchingMCpT = 0;
   RECOzEE_MatchingMCmass = 0;
   RECOzEE_MatchingMCEta = 0;
   RECOzEE_MatchingMCPhi = 0;
   RECOHzzMMMM_MatchingMCTruth = 0;
   RECOHzzMMMM_MatchingMCpT = 0;
   RECOHzzMMMM_MatchingMCmass = 0;
   RECOHzzMMMM_MatchingMCEta = 0;
   RECOHzzMMMM_MatchingMCPhi = 0;
   RECOHzzEEEE_MatchingMCTruth = 0;
   RECOHzzEEEE_MatchingMCpT = 0;
   RECOHzzEEEE_MatchingMCmass = 0;
   RECOHzzEEEE_MatchingMCEta = 0;
   RECOHzzEEEE_MatchingMCPhi = 0;
   RECOHzzEEMM_MatchingMCTruth = 0;
   RECOHzzEEMM_MatchingMCpT = 0;
   RECOHzzEEMM_MatchingMCmass = 0;
   RECOHzzEEMM_MatchingMCEta = 0;
   RECOHzzEEMM_MatchingMCPhi = 0;
   MC_weighting = 0;
   StdFitVertexX = 0;
   StdFitVertexY = 0;
   StdFitVertexZ = 0;
   StdFitVertexChi2r = 0;
   StdFitVertexProb = 0;
   StdFitVertexTrack_PT = 0;
   StdFitVertexTrack_ETA = 0;
   StdFitVertexTrack_PHI = 0;
   KinFitVertexX = 0;
   KinFitVertexY = 0;
   KinFitVertexZ = 0;
   KinFitVertexChi2r = 0;
   KinFitVertexProb = 0;
   StdFitVertexXMMMM = 0;
   StdFitVertexYMMMM = 0;
   StdFitVertexZMMMM = 0;
   StdFitVertexChi2rMMMM = 0;
   StdFitVertexProbMMMM = 0;
   StdFitVertexTrackMMMM_PT = 0;
   StdFitVertexTrackMMMM_ETA = 0;
   StdFitVertexTrackMMMM_PHI = 0;
   KinFitVertexXMMMM = 0;
   KinFitVertexYMMMM = 0;
   KinFitVertexZMMMM = 0;
   KinFitVertexChi2rMMMM = 0;
   KinFitVertexProbMMMM = 0;
   StdFitVertexXEEEE = 0;
   StdFitVertexYEEEE = 0;
   StdFitVertexZEEEE = 0;
   StdFitVertexChi2rEEEE = 0;
   StdFitVertexProbEEEE = 0;
   StdFitVertexTrackEEEE_PT = 0;
   StdFitVertexTrackEEEE_ETA = 0;
   StdFitVertexTrackEEEE_PHI = 0;
   KinFitVertexXEEEE = 0;
   KinFitVertexYEEEE = 0;
   KinFitVertexZEEEE = 0;
   KinFitVertexChi2rEEEE = 0;
   KinFitVertexProbEEEE = 0;
   StdFitVertexChi2rDiLep = 0;
   StdFitVertexProbDiLep = 0;
   StdFitVertexChi2rMMM = 0;
   StdFitVertexProbMMM = 0;
   StdFitVertexChi2rMME = 0;
   StdFitVertexProbMME = 0;
   StdFitVertexChi2rEEE = 0;
   StdFitVertexProbEEE = 0;
   StdFitVertexChi2rMEE = 0;
   StdFitVertexProbMEE = 0;
   ftsigma = 0;
   gdX = 0;
   gdY = 0;
   gdZ = 0;
   ftsigmalag = 0;
   gdlagX = 0;
   gdlagY = 0;
   gdlagZ = 0;
   gdlagProb = 0;
   gdlagNdof = 0;
   ftsigmaMMMM = 0;
   gdXMMMM = 0;
   gdYMMMM = 0;
   gdZMMMM = 0;
   ftsigmalagMMMM = 0;
   gdlagXMMMM = 0;
   gdlagYMMMM = 0;
   gdlagZMMMM = 0;
   gdlagProbMMMM = 0;
   gdlagNdofMMMM = 0;
   ftsigmaEEEE = 0;
   gdXEEEE = 0;
   gdYEEEE = 0;
   gdZEEEE = 0;
   ftsigmalagEEEE = 0;
   gdlagXEEEE = 0;
   gdlagYEEEE = 0;
   gdlagZEEEE = 0;
   gdlagProbEEEE = 0;
   gdlagNdofEEEE = 0;
   RECORF_2e2mu_cosTheta1_spin = 0;
   RECORF_2e2mu_cosTheta2_spin = 0;
   RECORF_2e2mu_cosThetaStar_spin = 0;
   RECORF_2e2mu_Phi_spin = 0;
   RECORF_2e2mu_Phi1_spin = 0;
   RECORF_2e2mu_Phi2_spin = 0;
   RECORF_2e2mu_phi1RF_spin = 0;
   RECORF_2e2mu_phi2RF_spin = 0;
   RECORF_2e2mu_MELA = 0;
   RECORF_4e_cosTheta1_spin = 0;
   RECORF_4e_cosTheta2_spin = 0;
   RECORF_4e_cosThetaStar_spin = 0;
   RECORF_4e_Phi_spin = 0;
   RECORF_4e_Phi1_spin = 0;
   RECORF_4e_Phi2_spin = 0;
   RECORF_4e_phi1RF_spin = 0;
   RECORF_4e_phi2RF_spin = 0;
   RECORF_4e_MELA = 0;
   RECORF_4mu_cosTheta1_spin = 0;
   RECORF_4mu_cosTheta2_spin = 0;
   RECORF_4mu_cosThetaStar_spin = 0;
   RECORF_4mu_Phi_spin = 0;
   RECORF_4mu_Phi1_spin = 0;
   RECORF_4mu_Phi2_spin = 0;
   RECORF_4mu_phi1RF_spin = 0;
   RECORF_4mu_phi2RF_spin = 0;
   RECORF_4mu_MELA = 0;
   nbPv = 0;
   Nbdof = 0;
   PositionRho = 0;
   PositionX = 0;
   PositionY = 0;
   PositionZ = 0;
   RECO_PFJET_CHARGE = 0;
   RECO_PFJET_ET = 0;
   RECO_PFJET_PT = 0;
   RECO_PFJET_ETA = 0;
   RECO_PFJET_PHI = 0;
   RECO_PFJET_PUID = 0;
   RECO_PFJET_PUID_MVA = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run", &Run, &b_irun);
   fChain->SetBranchAddress("Event", &Event, &b_ievt);
   fChain->SetBranchAddress("LumiSection", &LumiSection, &b_ils);
   fChain->SetBranchAddress("Avginstlumi", &Avginstlumi, &b_Avginstlumi);
   fChain->SetBranchAddress("RECO_NVTX", &RECO_NVTX, &b_RECO_NVTX);
   fChain->SetBranchAddress("RECO_VERTEX_x", &RECO_VERTEX_x, &b_RECO_VERTEX_x);
   fChain->SetBranchAddress("RECO_VERTEX_y", &RECO_VERTEX_y, &b_RECO_VERTEX_y);
   fChain->SetBranchAddress("RECO_VERTEX_z", &RECO_VERTEX_z, &b_RECO_VERTEX_z);
   fChain->SetBranchAddress("RECO_VERTEX_ndof", &RECO_VERTEX_ndof, &b_RECO_VERTEX_ndof);
   fChain->SetBranchAddress("RECO_VERTEX_chi2", &RECO_VERTEX_chi2, &b_RECO_VERTEX_chi2);
   fChain->SetBranchAddress("RECO_VERTEX_ntracks", &RECO_VERTEX_ntracks, &b_RECO_VERTEX_ntracks);
   fChain->SetBranchAddress("RECO_VERTEXPROB", &RECO_VERTEXPROB, &b_RECO_VERTEXPROB);
   fChain->SetBranchAddress("RECO_VERTEX_isValid", &RECO_VERTEX_isValid, &b_RECO_VERTEX_isValid);
   fChain->SetBranchAddress("RECO_VERTEX_TRACK_PT", &RECO_VERTEX_TRACK_PT, &b_RECO_VERTEX_TRACK_PT);
   fChain->SetBranchAddress("RECO_NTRACK", &RECO_NTRACK, &b_RECO_NTRACK);
   fChain->SetBranchAddress("RECO_TRACK_PT", &RECO_TRACK_PT, &b_RECO_TRACK_PT);
   fChain->SetBranchAddress("RECO_TRACK_ETA", &RECO_TRACK_ETA, &b_RECO_TRACK_ETA);
   fChain->SetBranchAddress("RECO_TRACK_PHI", &RECO_TRACK_PHI, &b_RECO_TRACK_PHI);
   fChain->SetBranchAddress("RECO_TRACK_CHI2", &RECO_TRACK_CHI2, &b_RECO_TRACK_CHI2);
   fChain->SetBranchAddress("RECO_TRACK_CHI2RED", &RECO_TRACK_CHI2RED, &b_RECO_TRACK_CHI2RED);
   fChain->SetBranchAddress("RECO_TRACK_CHI2PROB", &RECO_TRACK_CHI2PROB, &b_RECO_TRACK_CHI2PROB);
   fChain->SetBranchAddress("RECO_TRACK_NHITS", &RECO_TRACK_NHITS, &b_RECO_TRACK_NHITS);
   fChain->SetBranchAddress("RECO_TRACK_DXY", &RECO_TRACK_DXY, &b_RECO_TRACK_DXY);
   fChain->SetBranchAddress("RECO_TRACK_DXYERR", &RECO_TRACK_DXYERR, &b_RECO_TRACK_DXYERR);
   fChain->SetBranchAddress("RECO_TRACK_DZ", &RECO_TRACK_DZ, &b_RECO_TRACK_DZ);
   fChain->SetBranchAddress("RECO_TRACK_DZERR", &RECO_TRACK_DZERR, &b_RECO_TRACK_DZERR);
   fChain->SetBranchAddress("MC_GENMET", &MC_GENMET, &b_MC_GENMET);
   fChain->SetBranchAddress("HLTPathsFired", &HLTPathsFired, &b_HLTPathsFired);
   fChain->SetBranchAddress("RECO_nMuHLTMatch", &RECO_nMuHLTMatch, &b_RECO_nMuHLTMatch);
   fChain->SetBranchAddress("RECO_nEleHLTMatch", &RECO_nEleHLTMatch, &b_RECO_nEleHLTMatch);
   fChain->SetBranchAddress("RECOMU_PT_MuHLTMatch", &RECOMU_PT_MuHLTMatch, &b_RECOMU_PT_MuHLTMatch);
   fChain->SetBranchAddress("RECOMU_ETA_MuHLTMatch", &RECOMU_ETA_MuHLTMatch, &b_RECOMU_ETA_MuHLTMatch);
   fChain->SetBranchAddress("RECOMU_PHI_MuHLTMatch", &RECOMU_PHI_MuHLTMatch, &b_RECOMU_PHI_MuHLTMatch);
   fChain->SetBranchAddress("RECOELE_PT_EleHLTMatch", &RECOELE_PT_EleHLTMatch, &b_RECOELE_PT_EleHLTMatch);
   fChain->SetBranchAddress("RECOELE_ETA_EleHLTMatch", &RECOELE_ETA_EleHLTMatch, &b_RECOELE_ETA_EleHLTMatch);
   fChain->SetBranchAddress("RECOELE_PHI_EleHLTMatch", &RECOELE_PHI_EleHLTMatch, &b_RECOELE_PHI_EleHLTMatch);
   fChain->SetBranchAddress("MC_LEPT_PT", &MC_LEPT_PT, &b_MC_LEPT_PT);
   fChain->SetBranchAddress("MC_LEPT_ETA", &MC_LEPT_ETA, &b_MC_LEPT_ETA);
   fChain->SetBranchAddress("MC_LEPT_PHI", &MC_LEPT_PHI, &b_MC_LEPT_PHI);
   fChain->SetBranchAddress("MC_LEPT_THETA", &MC_LEPT_THETA, &b_MC_LEPT_THETA);
   fChain->SetBranchAddress("MC_LEPT_PDGID", &MC_LEPT_PDGID, &b_MC_LEPT_PDGID);
   fChain->SetBranchAddress("MC_Z_PT", &MC_Z_PT, &b_MC_Z_PT);
   fChain->SetBranchAddress("MC_Z_ETA", &MC_Z_ETA, &b_MC_Z_ETA);
   fChain->SetBranchAddress("MC_Z_PHI", &MC_Z_PHI, &b_MC_Z_PHI);
   fChain->SetBranchAddress("MC_Z_THETA", &MC_Z_THETA, &b_MC_Z_THETA);
   fChain->SetBranchAddress("MC_Z_MASS", &MC_Z_MASS, &b_MC_Z_MASS);
   fChain->SetBranchAddress("MC_Z_PDGID", &MC_Z_PDGID, &b_MC_Z_PDGID);
   fChain->SetBranchAddress("MC_fourl_MASS", &MC_fourl_MASS, &b_MC_fourl_MASS);
   fChain->SetBranchAddress("MC_fourl_PT", &MC_fourl_PT, &b_MC_fourl_PT);
   fChain->SetBranchAddress("MC_fourl_PDGID", &MC_fourl_PDGID, &b_MC_fourl_PDGID);
   fChain->SetBranchAddress("MC_ZZ_MASS", &MC_ZZ_MASS, &b_MC_ZZ_MASS);
   fChain->SetBranchAddress("MC_ZZ_PT", &MC_ZZ_PT, &b_MC_ZZ_PT);
   fChain->SetBranchAddress("MC_ZZ_ETA", &MC_ZZ_ETA, &b_MC_ZZ_ETA);
   fChain->SetBranchAddress("MC_ZZ_PHI", &MC_ZZ_PHI, &b_MC_ZZ_PHI);
   fChain->SetBranchAddress("MC_ZZ_THETA", &MC_ZZ_THETA, &b_MC_ZZ_THETA);
   fChain->SetBranchAddress("MC_ZZ_PDGID", &MC_ZZ_PDGID, &b_MC_ZZ_PDGID);
   fChain->SetBranchAddress("MC_E", &MC_E, &b_MC_E);
   fChain->SetBranchAddress("MC_PT", &MC_PT, &b_MC_PT);
   fChain->SetBranchAddress("MC_ETA", &MC_ETA, &b_MC_ETA);
   fChain->SetBranchAddress("MC_THETA", &MC_THETA, &b_MC_THETA);
   fChain->SetBranchAddress("MC_PHI", &MC_PHI, &b_MC_PHI);
   fChain->SetBranchAddress("MC_MASS", &MC_MASS, &b_MC_MASS);
   fChain->SetBranchAddress("MC_PDGID", &MC_PDGID, &b_MC_PDGID);
   fChain->SetBranchAddress("RECO_NMU", &RECO_NMU, &b_RECO_NMU);
   fChain->SetBranchAddress("RECOMU_isPFMu", &RECOMU_isPFMu, &b_RECOMU_isPFMu);
   fChain->SetBranchAddress("RECOMU_isGlobalMu", &RECOMU_isGlobalMu, &b_RECOMU_isGlobalMu);
   fChain->SetBranchAddress("RECOMU_isStandAloneMu", &RECOMU_isStandAloneMu, &b_RECOMU_isStandAloneMu);
   fChain->SetBranchAddress("RECOMU_isTrackerMu", &RECOMU_isTrackerMu, &b_RECOMU_isTrackerMu);
   fChain->SetBranchAddress("RECOMU_isCaloMu", &RECOMU_isCaloMu, &b_RECOMU_isCaloMu);
   fChain->SetBranchAddress("RECOMU_E", &RECOMU_E, &b_RECOMU_E);
   fChain->SetBranchAddress("RECOMU_PT", &RECOMU_PT, &b_RECOMU_PT);
   fChain->SetBranchAddress("RECOMU_P", &RECOMU_P, &b_RECOMU_P);
   fChain->SetBranchAddress("RECOMU_ETA", &RECOMU_ETA, &b_RECOMU_ETA);
   fChain->SetBranchAddress("RECOMU_THETA", &RECOMU_THETA, &b_RECOMU_THETA);
   fChain->SetBranchAddress("RECOMU_PHI", &RECOMU_PHI, &b_RECOMU_PHI);
   fChain->SetBranchAddress("RECOMU_MASS", &RECOMU_MASS, &b_RECOMU_MASS);
   fChain->SetBranchAddress("RECOMU_CHARGE", &RECOMU_CHARGE, &b_RECOMU_CHARGE);
   fChain->SetBranchAddress("RECOMU_COV", &RECOMU_COV, &b_RECOMU_COV);
   fChain->SetBranchAddress("RECOMU_TRACKISO", &RECOMU_TRACKISO, &b_RECOMU_TRACKISO);
   fChain->SetBranchAddress("RECOMU_TRACKISO_SUMPT", &RECOMU_TRACKISO_SUMPT, &b_RECOMU_TRACKISO_SUMPT);
   fChain->SetBranchAddress("RECOMU_ECALISO", &RECOMU_ECALISO, &b_RECOMU_ECALISO);
   fChain->SetBranchAddress("RECOMU_HCALISO", &RECOMU_HCALISO, &b_RECOMU_HCALISO);
   fChain->SetBranchAddress("RECOMU_X", &RECOMU_X, &b_RECOMU_X);
   fChain->SetBranchAddress("RECOMU_PFchHad", &RECOMU_PFchHad, &b_RECOMU_PFchHad);
   fChain->SetBranchAddress("RECOMU_PFneuHad", &RECOMU_PFneuHad, &b_RECOMU_PFneuHad);
   fChain->SetBranchAddress("RECOMU_PFphoton", &RECOMU_PFphoton, &b_RECOMU_PFphoton);
   fChain->SetBranchAddress("RECOMU_PFPUchAllPart", &RECOMU_PFPUchAllPart, &b_RECOMU_PFPUchAllPart);
   fChain->SetBranchAddress("RECOMU_PFX_dB", &RECOMU_PFX_dB, &b_RECOMU_PFX_dB);
   fChain->SetBranchAddress("RECOMU_PFX_rho", &RECOMU_PFX_rho, &b_RECOMU_PFX_rho);
   fChain->SetBranchAddress("RECOMU_SIP", &RECOMU_SIP, &b_RECOMU_SIP);
   fChain->SetBranchAddress("RECOMU_IP", &RECOMU_IP, &b_RECOMU_IP);
   fChain->SetBranchAddress("RECOMU_IPERROR", &RECOMU_IPERROR, &b_RECOMU_IPERROR);
   fChain->SetBranchAddress("RECOMU_SIP_KF", &RECOMU_SIP_KF, &b_RECOMU_SIP_KF);
   fChain->SetBranchAddress("RECOMU_IP_KF", &RECOMU_IP_KF, &b_RECOMU_IP_KF);
   fChain->SetBranchAddress("RECOMU_IPERROR_KF", &RECOMU_IPERROR_KF, &b_RECOMU_IPERROR_KF);
   fChain->SetBranchAddress("RECOMU_STIP", &RECOMU_STIP, &b_RECOMU_STIP);
   fChain->SetBranchAddress("RECOMU_SLIP", &RECOMU_SLIP, &b_RECOMU_SLIP);
   fChain->SetBranchAddress("RECOMU_TIP", &RECOMU_TIP, &b_RECOMU_TIP);
   fChain->SetBranchAddress("RECOMU_LIP", &RECOMU_LIP, &b_RECOMU_LIP);
   fChain->SetBranchAddress("RECOMU_TIPERROR", &RECOMU_TIPERROR, &b_RECOMU_TIPERROR);
   fChain->SetBranchAddress("RECOMU_LIPERROR", &RECOMU_LIPERROR, &b_RECOMU_LIPERROR);
   fChain->SetBranchAddress("RECOMU_numberOfMatches", &RECOMU_numberOfMatches, &b_RECOMU_numberOfMatches);
   fChain->SetBranchAddress("RECOMU_numberOfMatchedStations", &RECOMU_numberOfMatchedStations, &b_RECOMU_numberOfMatchedStations);
   fChain->SetBranchAddress("RECOMU_caloCompatibility", &RECOMU_caloCompatibility, &b_RECOMU_caloCompatibility);
   fChain->SetBranchAddress("RECOMU_segmentCompatibility", &RECOMU_segmentCompatibility, &b_RECOMU_segmentCompatibility);
   fChain->SetBranchAddress("RECOMU_glbmuPromptTight", &RECOMU_glbmuPromptTight, &b_RECOMU_glbmuPromptTight);
   fChain->SetBranchAddress("RECOMU_mubesttrkType", &RECOMU_mubesttrkType, &b_RECOMU_mubesttrkType);
   fChain->SetBranchAddress("RECOMU_mubesttrkDxy", &RECOMU_mubesttrkDxy, &b_RECOMU_mubesttrkDxy);
   fChain->SetBranchAddress("RECOMU_mubesttrkDxyB", &RECOMU_mubesttrkDxyB, &b_RECOMU_mubesttrkDxyB);
   fChain->SetBranchAddress("RECOMU_mubesttrkDxyError", &RECOMU_mubesttrkDxyError, &b_RECOMU_mubesttrkDxyError);
   fChain->SetBranchAddress("RECOMU_mubesttrkDz", &RECOMU_mubesttrkDz, &b_RECOMU_mubesttrkDz);
   fChain->SetBranchAddress("RECOMU_mubesttrkDzB", &RECOMU_mubesttrkDzB, &b_RECOMU_mubesttrkDzB);
   fChain->SetBranchAddress("RECOMU_mubesttrkDzError", &RECOMU_mubesttrkDzError, &b_RECOMU_mubesttrkDzError);
   fChain->SetBranchAddress("RECOMU_mubesttrkPTError", &RECOMU_mubesttrkPTError, &b_RECOMU_mubesttrkPTError);
   fChain->SetBranchAddress("RECOMU_mutrkPT", &RECOMU_mutrkPT, &b_RECOMU_mutrkPT);
   fChain->SetBranchAddress("RECOMU_mutrkPTError", &RECOMU_mutrkPTError, &b_RECOMU_mutrkPTError);
   fChain->SetBranchAddress("RECOMU_mutrkDxy", &RECOMU_mutrkDxy, &b_RECOMU_mutrkDxy);
   fChain->SetBranchAddress("RECOMU_mutrkDxyError", &RECOMU_mutrkDxyError, &b_RECOMU_mutrkDxyError);
   fChain->SetBranchAddress("RECOMU_mutrkDxyB", &RECOMU_mutrkDxyB, &b_RECOMU_mutrkDxyB);
   fChain->SetBranchAddress("RECOMU_mutrkDz", &RECOMU_mutrkDz, &b_RECOMU_mutrkDz);
   fChain->SetBranchAddress("RECOMU_mutrkDzError", &RECOMU_mutrkDzError, &b_RECOMU_mutrkDzError);
   fChain->SetBranchAddress("RECOMU_mutrkDzB", &RECOMU_mutrkDzB, &b_RECOMU_mutrkDzB);
   fChain->SetBranchAddress("RECOMU_mutrkChi2PerNdof", &RECOMU_mutrkChi2PerNdof, &b_RECOMU_mutrkChi2PerNdof);
   fChain->SetBranchAddress("RECOMU_mutrkCharge", &RECOMU_mutrkCharge, &b_RECOMU_mutrkCharge);
   fChain->SetBranchAddress("RECOMU_mutrkNHits", &RECOMU_mutrkNHits, &b_RECOMU_mutrkNHits);
   fChain->SetBranchAddress("RECOMU_mutrkNPixHits", &RECOMU_mutrkNPixHits, &b_RECOMU_mutrkNPixHits);
   fChain->SetBranchAddress("RECOMU_mutrkNStripHits", &RECOMU_mutrkNStripHits, &b_RECOMU_mutrkNStripHits);
   fChain->SetBranchAddress("RECOMU_mutrkNMuonHits", &RECOMU_mutrkNMuonHits, &b_RECOMU_mutrkNMuonHits);
   fChain->SetBranchAddress("RECOMU_mutrktrackerLayersWithMeasurement", &RECOMU_mutrktrackerLayersWithMeasurement, &b_RECOMU_mutrktrackerLayersWithMeasurement);
   fChain->SetBranchAddress("RECOMU_muInnertrkDxy", &RECOMU_muInnertrkDxy, &b_RECOMU_muInnertrkDxy);
   fChain->SetBranchAddress("RECOMU_muInnertrkDxyError", &RECOMU_muInnertrkDxyError, &b_RECOMU_muInnertrkDxyError);
   fChain->SetBranchAddress("RECOMU_muInnertrkDxyB", &RECOMU_muInnertrkDxyB, &b_RECOMU_muInnertrkDxyB);
   fChain->SetBranchAddress("RECOMU_muInnertrkDz", &RECOMU_muInnertrkDz, &b_RECOMU_muInnertrkDz);
   fChain->SetBranchAddress("RECOMU_muInnertrkDzError", &RECOMU_muInnertrkDzError, &b_RECOMU_muInnertrkDzError);
   fChain->SetBranchAddress("RECOMU_muInnertrkDzB", &RECOMU_muInnertrkDzB, &b_RECOMU_muInnertrkDzB);
   fChain->SetBranchAddress("RECOMU_muInnertrkChi2PerNdof", &RECOMU_muInnertrkChi2PerNdof, &b_RECOMU_muInnertrkChi2PerNdof);
   fChain->SetBranchAddress("RECOMU_muInnertrktrackerLayersWithMeasurement", &RECOMU_muInnertrktrackerLayersWithMeasurement, &b_RECOMU_muInnertrktrackerLayersWithMeasurement);
   fChain->SetBranchAddress("RECOMU_muInnertrkPT", &RECOMU_muInnertrkPT, &b_RECOMU_muInnertrkPT);
   fChain->SetBranchAddress("RECOMU_muInnertrkPTError", &RECOMU_muInnertrkPTError, &b_RECOMU_muInnertrkPTError);
   fChain->SetBranchAddress("RECOMU_muInnertrkCharge", &RECOMU_muInnertrkCharge, &b_RECOMU_muInnertrkCharge);
   fChain->SetBranchAddress("RECOMU_muInnertrkNHits", &RECOMU_muInnertrkNHits, &b_RECOMU_muInnertrkNHits);
   fChain->SetBranchAddress("RECOMU_muInnertrkNPixHits", &RECOMU_muInnertrkNPixHits, &b_RECOMU_muInnertrkNPixHits);
   fChain->SetBranchAddress("RECOMU_muInnertrkNStripHits", &RECOMU_muInnertrkNStripHits, &b_RECOMU_muInnertrkNStripHits);
   fChain->SetBranchAddress("RECOMU_trkmuArbitration", &RECOMU_trkmuArbitration, &b_RECOMU_trkmuArbitration);
   fChain->SetBranchAddress("RECOMU_trkmu2DCompatibilityLoose", &RECOMU_trkmu2DCompatibilityLoose, &b_RECOMU_trkmu2DCompatibilityLoose);
   fChain->SetBranchAddress("RECOMU_trkmu2DCompatibilityTight", &RECOMU_trkmu2DCompatibilityTight, &b_RECOMU_trkmu2DCompatibilityTight);
   fChain->SetBranchAddress("RECOMU_trkmuOneStationLoose", &RECOMU_trkmuOneStationLoose, &b_RECOMU_trkmuOneStationLoose);
   fChain->SetBranchAddress("RECOMU_trkmuOneStationTight", &RECOMU_trkmuOneStationTight, &b_RECOMU_trkmuOneStationTight);
   fChain->SetBranchAddress("RECOMU_trkmuLastStationLoose", &RECOMU_trkmuLastStationLoose, &b_RECOMU_trkmuLastStationLoose);
   fChain->SetBranchAddress("RECOMU_trkmuLastStationTight", &RECOMU_trkmuLastStationTight, &b_RECOMU_trkmuLastStationTight);
   fChain->SetBranchAddress("RECOMU_trkmuOneStationAngLoose", &RECOMU_trkmuOneStationAngLoose, &b_RECOMU_trkmuOneStationAngLoose);
   fChain->SetBranchAddress("RECOMU_trkmuOneStationAngTight", &RECOMU_trkmuOneStationAngTight, &b_RECOMU_trkmuOneStationAngTight);
   fChain->SetBranchAddress("RECOMU_trkmuLastStationAngLoose", &RECOMU_trkmuLastStationAngLoose, &b_RECOMU_trkmuLastStationAngLoose);
   fChain->SetBranchAddress("RECOMU_trkmuLastStationAngTight", &RECOMU_trkmuLastStationAngTight, &b_RECOMU_trkmuLastStationAngTight);
   fChain->SetBranchAddress("RECOMU_trkmuLastStationOptimizedLowPtLoose", &RECOMU_trkmuLastStationOptimizedLowPtLoose, &b_RECOMU_trkmuLastStationOptimizedLowPtLoose);
   fChain->SetBranchAddress("RECOMU_trkmuLastStationOptimizedLowPtTight", &RECOMU_trkmuLastStationOptimizedLowPtTight, &b_RECOMU_trkmuLastStationOptimizedLowPtTight);
   fChain->SetBranchAddress("RECOMU_MatchingMCTruth", &RECOMU_MatchingMCTruth, &b_RECOMU_MatchingMCTruth);
   fChain->SetBranchAddress("RECOMU_MatchingMCpT", &RECOMU_MatchingMCpT, &b_RECOMU_MatchingMCpT);
   fChain->SetBranchAddress("RECOMU_MatchingMCEta", &RECOMU_MatchingMCEta, &b_RECOMU_MatchingMCEta);
   fChain->SetBranchAddress("RECOMU_MatchingMCPhi", &RECOMU_MatchingMCPhi, &b_RECOMU_MatchingMCPhi);
   fChain->SetBranchAddress("RECO_NELE", &RECO_NELE, &b_RECO_NELE);
   fChain->SetBranchAddress("RECOELE_isEcalDriven", &RECOELE_isEcalDriven, &b_RECOELE_isEcalDriven);
   fChain->SetBranchAddress("RECOELE_isTrackerDriven", &RECOELE_isTrackerDriven, &b_RECOELE_isTrackerDriven);
   fChain->SetBranchAddress("RECOELE_CHARGE", &RECOELE_CHARGE, &b_RECOELE_CHARGE);
   fChain->SetBranchAddress("RECOELE_E", &RECOELE_E, &b_RECOELE_E);
   fChain->SetBranchAddress("RECOELE_PT", &RECOELE_PT, &b_RECOELE_PT);
   fChain->SetBranchAddress("RECOELE_P", &RECOELE_P, &b_RECOELE_P);
   fChain->SetBranchAddress("RECOELE_ETA", &RECOELE_ETA, &b_RECOELE_ETA);
   fChain->SetBranchAddress("RECOELE_THETA", &RECOELE_THETA, &b_RECOELE_THETA);
   fChain->SetBranchAddress("RECOELE_PHI", &RECOELE_PHI, &b_RECOELE_PHI);
   fChain->SetBranchAddress("RECOELE_MASS", &RECOELE_MASS, &b_RECOELE_MASS);
   fChain->SetBranchAddress("ele_sclRawE", &ele_sclRawE, &b_ele_sclRawE);
   fChain->SetBranchAddress("RECOELE_scl_E", &RECOELE_scl_E, &b_RECOELE_scl_E);
   fChain->SetBranchAddress("RECOELE_scl_Et", &RECOELE_scl_Et, &b_RECOELE_scl_Et);
   fChain->SetBranchAddress("RECOELE_scl_Eta", &RECOELE_scl_Eta, &b_RECOELE_scl_Eta);
   fChain->SetBranchAddress("RECOELE_scl_Phi", &RECOELE_scl_Phi, &b_RECOELE_scl_Phi);
   fChain->SetBranchAddress("ele_sclX", &ele_sclX, &b_ele_sclX);
   fChain->SetBranchAddress("ele_sclY", &ele_sclY, &b_ele_sclY);
   fChain->SetBranchAddress("ele_sclZ", &ele_sclZ, &b_ele_sclZ);
   fChain->SetBranchAddress("RECOELE_PTError", &RECOELE_PTError, &b_RECOELE_PTError);
   fChain->SetBranchAddress("RECOELE_COV", &RECOELE_COV, &b_RECOELE_COV);
   fChain->SetBranchAddress("RECOELE_EGMECALISO", &RECOELE_EGMECALISO, &b_RECOELE_EGMECALISO);
   fChain->SetBranchAddress("RECOELE_EGMHCALISO", &RECOELE_EGMHCALISO, &b_RECOELE_EGMHCALISO);
   fChain->SetBranchAddress("RECOELE_EGMX", &RECOELE_EGMX, &b_RECOELE_EGMX);
   fChain->SetBranchAddress("RECOELE_PFchAllPart", &RECOELE_PFchAllPart, &b_RECOELE_PFchAllPart);
   fChain->SetBranchAddress("RECOELE_PFchHad", &RECOELE_PFchHad, &b_RECOELE_PFchHad);
   fChain->SetBranchAddress("RECOELE_PFneuHad", &RECOELE_PFneuHad, &b_RECOELE_PFneuHad);
   fChain->SetBranchAddress("RECOELE_PFphoton", &RECOELE_PFphoton, &b_RECOELE_PFphoton);
   fChain->SetBranchAddress("RECOELE_PFPUchAllPart", &RECOELE_PFPUchAllPart, &b_RECOELE_PFPUchAllPart);
   fChain->SetBranchAddress("RECOELE_PFX_dB", &RECOELE_PFX_dB, &b_RECOELE_PFX_dB);
   fChain->SetBranchAddress("RECOELE_PFX_rho", &RECOELE_PFX_rho, &b_RECOELE_PFX_rho);
   fChain->SetBranchAddress("RECOELE_SIP", &RECOELE_SIP, &b_RECOELE_SIP);
   fChain->SetBranchAddress("RECOELE_IP", &RECOELE_IP, &b_RECOELE_IP);
   fChain->SetBranchAddress("RECOELE_IPERROR", &RECOELE_IPERROR, &b_RECOELE_IPERROR);
   fChain->SetBranchAddress("RECOELE_SIP_KF", &RECOELE_SIP_KF, &b_RECOELE_SIP_KF);
   fChain->SetBranchAddress("RECOELE_IP_KF", &RECOELE_IP_KF, &b_RECOELE_IP_KF);
   fChain->SetBranchAddress("RECOELE_IPERROR_KF", &RECOELE_IPERROR_KF, &b_RECOELE_IPERROR_KF);
   fChain->SetBranchAddress("RECOELE_STIP", &RECOELE_STIP, &b_RECOELE_STIP);
   fChain->SetBranchAddress("RECOELE_SLIP", &RECOELE_SLIP, &b_RECOELE_SLIP);
   fChain->SetBranchAddress("RECOELE_TIP", &RECOELE_TIP, &b_RECOELE_TIP);
   fChain->SetBranchAddress("RECOELE_LIP", &RECOELE_LIP, &b_RECOELE_LIP);
   fChain->SetBranchAddress("RECOELE_TIPERROR", &RECOELE_TIPERROR, &b_RECOELE_TIPERROR);
   fChain->SetBranchAddress("RECOELE_LIPERROR", &RECOELE_LIPERROR, &b_RECOELE_LIPERROR);
   fChain->SetBranchAddress("RECOELE_gsftrack_NPixHits", &RECOELE_gsftrack_NPixHits, &b_RECOELE_gsftrack_NPixHits);
   fChain->SetBranchAddress("RECOELE_gsftrack_NStripHits", &RECOELE_gsftrack_NStripHits, &b_RECOELE_gsftrack_NStripHits);
   fChain->SetBranchAddress("RECOELE_gsftrack_chi2", &RECOELE_gsftrack_chi2, &b_RECOELE_gsftrack_chi2);
   fChain->SetBranchAddress("RECOELE_gsftrack_dxyB", &RECOELE_gsftrack_dxyB, &b_RECOELE_gsftrack_dxyB);
   fChain->SetBranchAddress("RECOELE_gsftrack_dxy", &RECOELE_gsftrack_dxy, &b_RECOELE_gsftrack_dxy);
   fChain->SetBranchAddress("RECOELE_gsftrack_dxyError", &RECOELE_gsftrack_dxyError, &b_RECOELE_gsftrack_dxyError);
   fChain->SetBranchAddress("RECOELE_gsftrack_dzB", &RECOELE_gsftrack_dzB, &b_RECOELE_gsftrack_dzB);
   fChain->SetBranchAddress("RECOELE_gsftrack_dz", &RECOELE_gsftrack_dz, &b_RECOELE_gsftrack_dz);
   fChain->SetBranchAddress("RECOELE_gsftrack_dzError", &RECOELE_gsftrack_dzError, &b_RECOELE_gsftrack_dzError);
   fChain->SetBranchAddress("RECOELE_gsftrack_losthits", &RECOELE_gsftrack_losthits, &b_RECOELE_gsftrack_losthits);
   fChain->SetBranchAddress("RECOELE_gsftrack_validhits", &RECOELE_gsftrack_validhits, &b_RECOELE_gsftrack_validhits);
   fChain->SetBranchAddress("RECOELE_gsftrack_expected_inner_hits", &RECOELE_gsftrack_expected_inner_hits, &b_RECOELE_gsftrack_expected_inner_hits);
   fChain->SetBranchAddress("RECOELE_ep", &RECOELE_ep, &b_RECOELE_ep);
   fChain->SetBranchAddress("RECOELE_eSeedp", &RECOELE_eSeedp, &b_RECOELE_eSeedp);
   fChain->SetBranchAddress("RECOELE_eSeedpout", &RECOELE_eSeedpout, &b_RECOELE_eSeedpout);
   fChain->SetBranchAddress("RECOELE_eElepout", &RECOELE_eElepout, &b_RECOELE_eElepout);
   fChain->SetBranchAddress("RECOELE_deltaEtaIn", &RECOELE_deltaEtaIn, &b_RECOELE_deltaEtaIn);
   fChain->SetBranchAddress("RECOELE_deltaEtaSeed", &RECOELE_deltaEtaSeed, &b_RECOELE_deltaEtaSeed);
   fChain->SetBranchAddress("RECOELE_deltaEtaEle", &RECOELE_deltaEtaEle, &b_RECOELE_deltaEtaEle);
   fChain->SetBranchAddress("RECOELE_deltaPhiIn", &RECOELE_deltaPhiIn, &b_RECOELE_deltaPhiIn);
   fChain->SetBranchAddress("RECOELE_deltaPhiSeed", &RECOELE_deltaPhiSeed, &b_RECOELE_deltaPhiSeed);
   fChain->SetBranchAddress("RECOELE_deltaPhiEle", &RECOELE_deltaPhiEle, &b_RECOELE_deltaPhiEle);
   fChain->SetBranchAddress("RECOELE_isbarrel", &RECOELE_isbarrel, &b_RECOELE_isbarrel);
   fChain->SetBranchAddress("RECOELE_isendcap", &RECOELE_isendcap, &b_RECOELE_isendcap);
   fChain->SetBranchAddress("RECOELE_isGap", &RECOELE_isGap, &b_RECOELE_isGap);
   fChain->SetBranchAddress("RECOELE_isEBetaGap", &RECOELE_isEBetaGap, &b_RECOELE_isEBetaGap);
   fChain->SetBranchAddress("RECOELE_isEBphiGap", &RECOELE_isEBphiGap, &b_RECOELE_isEBphiGap);
   fChain->SetBranchAddress("RECOELE_isEEdeeGap", &RECOELE_isEEdeeGap, &b_RECOELE_isEEdeeGap);
   fChain->SetBranchAddress("RECOELE_isEEringGap", &RECOELE_isEEringGap, &b_RECOELE_isEEringGap);
   fChain->SetBranchAddress("RECOELE_sigmaIetaIeta", &RECOELE_sigmaIetaIeta, &b_RECOELE_sigmaIetaIeta);
   fChain->SetBranchAddress("RECOELE_sigmaEtaEta", &RECOELE_sigmaEtaEta, &b_RECOELE_sigmaEtaEta);
   fChain->SetBranchAddress("RECOELE_e15", &RECOELE_e15, &b_RECOELE_e15);
   fChain->SetBranchAddress("RECOELE_e25max", &RECOELE_e25max, &b_RECOELE_e25max);
   fChain->SetBranchAddress("RECOELE_e55", &RECOELE_e55, &b_RECOELE_e55);
   fChain->SetBranchAddress("RECOELE_he", &RECOELE_he, &b_RECOELE_he);
   fChain->SetBranchAddress("RECOELE_fbrem", &RECOELE_fbrem, &b_RECOELE_fbrem);
   fChain->SetBranchAddress("RECOELE_nbrems", &RECOELE_nbrems, &b_RECOELE_nbrems);
   fChain->SetBranchAddress("RECOELE_Class", &RECOELE_Class, &b_RECOELE_Class);
   fChain->SetBranchAddress("RECOELE_fbrem_mean", &RECOELE_fbrem_mean, &b_RECOELE_fbrem_mean);
   fChain->SetBranchAddress("RECOELE_fbrem_mode", &RECOELE_fbrem_mode, &b_RECOELE_fbrem_mode);
   fChain->SetBranchAddress("RECOELE_ecalEnergy", &RECOELE_ecalEnergy, &b_RECOELE_ecalEnergy);
   fChain->SetBranchAddress("ele_seedSubdet2", &ele_seedSubdet2, &b_ele_seedSubdet2);
   fChain->SetBranchAddress("ele_seedDphi2", &ele_seedDphi2, &b_ele_seedDphi2);
   fChain->SetBranchAddress("ele_seedDrz2", &ele_seedDrz2, &b_ele_seedDrz2);
   fChain->SetBranchAddress("ele_seedSubdet1", &ele_seedSubdet1, &b_ele_seedSubdet1);
   fChain->SetBranchAddress("ele_seedDphi1", &ele_seedDphi1, &b_ele_seedDphi1);
   fChain->SetBranchAddress("ele_seedDrz1", &ele_seedDrz1, &b_ele_seedDrz1);
   fChain->SetBranchAddress("RECOELE_mvaTrigV0", &RECOELE_mvaTrigV0, &b_RECOELE_mvaTrigV0);
   fChain->SetBranchAddress("RECOELE_mvaNonTrigV0", &RECOELE_mvaNonTrigV0, &b_RECOELE_mvaNonTrigV0);
   fChain->SetBranchAddress("ConvMapDist", &ConvMapDist, &b_ConvMapDist);
   fChain->SetBranchAddress("ConvMapDcot", &ConvMapDcot, &b_ConvMapDcot);
   fChain->SetBranchAddress("RECOELE_MatchingMCTruth", &RECOELE_MatchingMCTruth, &b_RECOELE_MatchingMCTruth);
   fChain->SetBranchAddress("RECOELE_MatchingMCpT", &RECOELE_MatchingMCpT, &b_RECOELE_MatchingMCpT);
   fChain->SetBranchAddress("RECOELE_MatchingMCEta", &RECOELE_MatchingMCEta, &b_RECOELE_MatchingMCEta);
   fChain->SetBranchAddress("RECOELE_MatchingMCPhi", &RECOELE_MatchingMCPhi, &b_RECOELE_MatchingMCPhi);
   fChain->SetBranchAddress("RECO_NPHOT", &RECO_NPHOT, &b_RECO_NPHOT);
   fChain->SetBranchAddress("RECOPHOT_PT", &RECOPHOT_PT, &b_RECOPHOT_PT);
   fChain->SetBranchAddress("RECOPHOT_ETA", &RECOPHOT_ETA, &b_RECOPHOT_ETA);
   fChain->SetBranchAddress("RECOPHOT_PHI", &RECOPHOT_PHI, &b_RECOPHOT_PHI);
   fChain->SetBranchAddress("RECOPHOT_THETA", &RECOPHOT_THETA, &b_RECOPHOT_THETA);
   fChain->SetBranchAddress("RECO_NPFPHOT", &RECO_NPFPHOT, &b_RECO_NPFPHOT);
   fChain->SetBranchAddress("RECOPFPHOT_PT", &RECOPFPHOT_PT, &b_RECOPFPHOT_PT);
   fChain->SetBranchAddress("RECOPFPHOT_PTError", &RECOPFPHOT_PTError, &b_RECOPFPHOT_PTError);
   fChain->SetBranchAddress("RECOPFPHOT_ETA", &RECOPFPHOT_ETA, &b_RECOPFPHOT_ETA);
   fChain->SetBranchAddress("RECOPFPHOT_PHI", &RECOPFPHOT_PHI, &b_RECOPFPHOT_PHI);
   fChain->SetBranchAddress("RECOPFPHOT_THETA", &RECOPFPHOT_THETA, &b_RECOPFPHOT_THETA);
   fChain->SetBranchAddress("RECOPFPHOT_PFchAllPart", &RECOPFPHOT_PFchAllPart, &b_RECOPFPHOT_PFchAllPart);
   fChain->SetBranchAddress("RECOPFPHOT_PFchHad", &RECOPFPHOT_PFchHad, &b_RECOPFPHOT_PFchHad);
   fChain->SetBranchAddress("RECOPFPHOT_PFneuHad", &RECOPFPHOT_PFneuHad, &b_RECOPFPHOT_PFneuHad);
   fChain->SetBranchAddress("RECOPFPHOT_PFphoton", &RECOPFPHOT_PFphoton, &b_RECOPFPHOT_PFphoton);
   fChain->SetBranchAddress("RECOPFPHOT_PFPUchAllPart", &RECOPFPHOT_PFPUchAllPart, &b_RECOPFPHOT_PFPUchAllPart);
   fChain->SetBranchAddress("RECOPFPHOT_PFX_rho", &RECOPFPHOT_PFX_rho, &b_RECOPFPHOT_PFX_rho);
   fChain->SetBranchAddress("RECOPHOT_MatchingMCTruth", &RECOPHOT_MatchingMCTruth, &b_RECOPHOT_MatchingMCTruth);
   fChain->SetBranchAddress("RECOPHOT_MatchingMCpT", &RECOPHOT_MatchingMCpT, &b_RECOPHOT_MatchingMCpT);
   fChain->SetBranchAddress("RECOPHOT_MatchingMCEta", &RECOPHOT_MatchingMCEta, &b_RECOPHOT_MatchingMCEta);
   fChain->SetBranchAddress("RECOPHOT_MatchingMCPhi", &RECOPHOT_MatchingMCPhi, &b_RECOPHOT_MatchingMCPhi);
   fChain->SetBranchAddress("PFMet_pt", &PFMet_pt, &b_PFMet_pt);
   fChain->SetBranchAddress("PFMet_eta", &PFMet_eta, &b_PFMet_eta);
   fChain->SetBranchAddress("PFMet_phi", &PFMet_phi, &b_PFMet_phi);
   fChain->SetBranchAddress("PFMet_en", &PFMet_en, &b_PFMet_en);
   fChain->SetBranchAddress("PFMet_px", &PFMet_px, &b_PFMet_px);
   fChain->SetBranchAddress("PFMet_py", &PFMet_py, &b_PFMet_py);
   fChain->SetBranchAddress("PFMet_pz", &PFMet_pz, &b_PFMet_pz);
   fChain->SetBranchAddress("PFMet_sumEt", &PFMet_sumEt, &b_PFMet_sumEt);
   fChain->SetBranchAddress("METSign", &METSign, &b_METSign);
   fChain->SetBranchAddress("tCHighEff_nb", &tCHighEff_nb, &b_tCHighEff_nb);
   fChain->SetBranchAddress("tCHighEff_BTagJet_PT", &tCHighEff_BTagJet_PT, &b_tCHighEff_BTagJet_PT);
   fChain->SetBranchAddress("tCHighEff_BTagJet_ETA", &tCHighEff_BTagJet_ETA, &b_tCHighEff_BTagJet_ETA);
   fChain->SetBranchAddress("tCHighEff_BTagJet_PHI", &tCHighEff_BTagJet_PHI, &b_tCHighEff_BTagJet_PHI);
   fChain->SetBranchAddress("tCHighEff_BTagJet_DISCR", &tCHighEff_BTagJet_DISCR, &b_tCHighEff_BTagJet_DISCR);
   fChain->SetBranchAddress("tCHighPur_nb", &tCHighPur_nb, &b_tCHighPur_nb);
   fChain->SetBranchAddress("tCHighPur_BTagJet_PT", &tCHighPur_BTagJet_PT, &b_tCHighPur_BTagJet_PT);
   fChain->SetBranchAddress("tCHighPur_BTagJet_ETA", &tCHighPur_BTagJet_ETA, &b_tCHighPur_BTagJet_ETA);
   fChain->SetBranchAddress("tCHighPur_BTagJet_PHI", &tCHighPur_BTagJet_PHI, &b_tCHighPur_BTagJet_PHI);
   fChain->SetBranchAddress("tCHighPur_BTagJet_DISCR", &tCHighPur_BTagJet_DISCR, &b_tCHighPur_BTagJet_DISCR);
   fChain->SetBranchAddress("cSV_nb", &cSV_nb, &b_cSV_nb);
   fChain->SetBranchAddress("cSV_BTagJet_PT", &cSV_BTagJet_PT, &b_cSV_BTagJet_PT);
   fChain->SetBranchAddress("cSV_BTagJet_ETA", &cSV_BTagJet_ETA, &b_cSV_BTagJet_ETA);
   fChain->SetBranchAddress("cSV_BTagJet_PHI", &cSV_BTagJet_PHI, &b_cSV_BTagJet_PHI);
   fChain->SetBranchAddress("cSV_BTagJet_DISCR", &cSV_BTagJet_DISCR, &b_cSV_BTagJet_DISCR);
   fChain->SetBranchAddress("cSV_BTagJet_ET", &cSV_BTagJet_ET, &b_cSV_BTagJet_ET);
   fChain->SetBranchAddress("RECO_ZMM_MASS", &RECO_ZMM_MASS, &b_RECO_ZMM_MASS);
   fChain->SetBranchAddress("RECO_ZEE_MASS", &RECO_ZEE_MASS, &b_RECO_ZEE_MASS);
   fChain->SetBranchAddress("RECO_DiLep_MASS", &RECO_DiLep_MASS, &b_RECO_DiLep_MASS);
   fChain->SetBranchAddress("RECO_ZMM_PT", &RECO_ZMM_PT, &b_RECO_ZMM_PT);
   fChain->SetBranchAddress("RECO_ZEE_PT", &RECO_ZEE_PT, &b_RECO_ZEE_PT);
   fChain->SetBranchAddress("RECO_DiLep_PT", &RECO_DiLep_PT, &b_RECO_DiLep_PT);
   fChain->SetBranchAddress("RECO_ZMM_ETA", &RECO_ZMM_ETA, &b_RECO_ZMM_ETA);
   fChain->SetBranchAddress("RECO_ZEE_ETA", &RECO_ZEE_ETA, &b_RECO_ZEE_ETA);
   fChain->SetBranchAddress("RECO_DiLep_ETA", &RECO_DiLep_ETA, &b_RECO_DiLep_ETA);
   fChain->SetBranchAddress("RECO_ZMM_PHI", &RECO_ZMM_PHI, &b_RECO_ZMM_PHI);
   fChain->SetBranchAddress("RECO_ZEE_PHI", &RECO_ZEE_PHI, &b_RECO_ZEE_PHI);
   fChain->SetBranchAddress("RECO_DiLep_PHI", &RECO_DiLep_PHI, &b_RECO_DiLep_PHI);
   fChain->SetBranchAddress("RECO_ZMMss_MASS", &RECO_ZMMss_MASS, &b_RECO_ZMMss_MASS);
   fChain->SetBranchAddress("RECO_ZEEss_MASS", &RECO_ZEEss_MASS, &b_RECO_ZEEss_MASS);
   fChain->SetBranchAddress("RECO_ZEM_MASS", &RECO_ZEM_MASS, &b_RECO_ZEM_MASS);
   fChain->SetBranchAddress("RECO_ZMMss_PT", &RECO_ZMMss_PT, &b_RECO_ZMMss_PT);
   fChain->SetBranchAddress("RECO_ZEEss_PT", &RECO_ZEEss_PT, &b_RECO_ZEEss_PT);
   fChain->SetBranchAddress("RECO_ZEM_PT", &RECO_ZEM_PT, &b_RECO_ZEM_PT);
   fChain->SetBranchAddress("RECO_ZMMss_ETA", &RECO_ZMMss_ETA, &b_RECO_ZMMss_ETA);
   fChain->SetBranchAddress("RECO_ZEEss_ETA", &RECO_ZEEss_ETA, &b_RECO_ZEEss_ETA);
   fChain->SetBranchAddress("RECO_ZEM_ETA", &RECO_ZEM_ETA, &b_RECO_ZEM_ETA);
   fChain->SetBranchAddress("RECO_ZMMss_PHI", &RECO_ZMMss_PHI, &b_RECO_ZMMss_PHI);
   fChain->SetBranchAddress("RECO_ZEEss_PHI", &RECO_ZEEss_PHI, &b_RECO_ZEEss_PHI);
   fChain->SetBranchAddress("RECO_ZEM_PHI", &RECO_ZEM_PHI, &b_RECO_ZEM_PHI);
   fChain->SetBranchAddress("RECO_MMMM_MASS", &RECO_MMMM_MASS, &b_RECO_MMMM_MASS);
   fChain->SetBranchAddress("RECO_MMMM_PT", &RECO_MMMM_PT, &b_RECO_MMMM_PT);
   fChain->SetBranchAddress("RECO_MMMM_ETA", &RECO_MMMM_ETA, &b_RECO_MMMM_ETA);
   fChain->SetBranchAddress("RECO_MMMM_PHI", &RECO_MMMM_PHI, &b_RECO_MMMM_PHI);
   fChain->SetBranchAddress("RECO_MMMM_MASS_REFIT", &RECO_MMMM_MASS_REFIT, &b_RECO_MMMM_MASS_REFIT);
   fChain->SetBranchAddress("RECO_EEEE_MASS", &RECO_EEEE_MASS, &b_RECO_EEEE_MASS);
   fChain->SetBranchAddress("RECO_EEEE_PT", &RECO_EEEE_PT, &b_RECO_EEEE_PT);
   fChain->SetBranchAddress("RECO_EEEE_ETA", &RECO_EEEE_ETA, &b_RECO_EEEE_ETA);
   fChain->SetBranchAddress("RECO_EEEE_PHI", &RECO_EEEE_PHI, &b_RECO_EEEE_PHI);
   fChain->SetBranchAddress("RECO_EEEE_MASS_REFIT", &RECO_EEEE_MASS_REFIT, &b_RECO_EEEE_MASS_REFIT);
   fChain->SetBranchAddress("RECO_EEMM_MASS", &RECO_EEMM_MASS, &b_RECO_EEMM_MASS);
   fChain->SetBranchAddress("RECO_EEMM_PT", &RECO_EEMM_PT, &b_RECO_EEMM_PT);
   fChain->SetBranchAddress("RECO_EEMM_ETA", &RECO_EEMM_ETA, &b_RECO_EEMM_ETA);
   fChain->SetBranchAddress("RECO_EEMM_PHI", &RECO_EEMM_PHI, &b_RECO_EEMM_PHI);
   fChain->SetBranchAddress("RECO_EEMM_MASS_REFIT", &RECO_EEMM_MASS_REFIT, &b_RECO_EEMM_MASS_REFIT);
   fChain->SetBranchAddress("RECO_LLL0_MASS", &RECO_LLL0_MASS, &b_RECO_LLL0_MASS);
   fChain->SetBranchAddress("RECO_LLL1_MASS", &RECO_LLL1_MASS, &b_RECO_LLL1_MASS);
   fChain->SetBranchAddress("RECO_LLL2_MASS", &RECO_LLL2_MASS, &b_RECO_LLL2_MASS);
   fChain->SetBranchAddress("RECO_LLL3_MASS", &RECO_LLL3_MASS, &b_RECO_LLL3_MASS);
   fChain->SetBranchAddress("RECO_LLL0_PT", &RECO_LLL0_PT, &b_RECO_LLL0_PT);
   fChain->SetBranchAddress("RECO_LLL1_PT", &RECO_LLL1_PT, &b_RECO_LLL1_PT);
   fChain->SetBranchAddress("RECO_LLL2_PT", &RECO_LLL2_PT, &b_RECO_LLL2_PT);
   fChain->SetBranchAddress("RECO_LLL3_PT", &RECO_LLL3_PT, &b_RECO_LLL3_PT);
   fChain->SetBranchAddress("RECO_LLLl0_MASS", &RECO_LLLl0_MASS, &b_RECO_LLLl0_MASS);
   fChain->SetBranchAddress("RECO_LLLl1_MASS", &RECO_LLLl1_MASS, &b_RECO_LLLl1_MASS);
   fChain->SetBranchAddress("RECO_LLLl0_PT", &RECO_LLLl0_PT, &b_RECO_LLLl0_PT);
   fChain->SetBranchAddress("RECO_LLLl1_PT", &RECO_LLLl1_PT, &b_RECO_LLLl1_PT);
   fChain->SetBranchAddress("RECO_LLLL0ss_MASS", &RECO_LLLL0ss_MASS, &b_RECO_LLLL0ss_MASS);
   fChain->SetBranchAddress("RECO_LLLL0ss_PT", &RECO_LLLL0ss_PT, &b_RECO_LLLL0ss_PT);
   fChain->SetBranchAddress("RECO_LLLL1ss_MASS", &RECO_LLLL1ss_MASS, &b_RECO_LLLL1ss_MASS);
   fChain->SetBranchAddress("RECO_LLLL1ss_PT", &RECO_LLLL1ss_PT, &b_RECO_LLLL1ss_PT);
   fChain->SetBranchAddress("RECO_LLLL2ss_MASS", &RECO_LLLL2ss_MASS, &b_RECO_LLLL2ss_MASS);
   fChain->SetBranchAddress("RECO_LLLL2ss_PT", &RECO_LLLL2ss_PT, &b_RECO_LLLL2ss_PT);
   fChain->SetBranchAddress("RECO_LLLL_MASS", &RECO_LLLL_MASS, &b_RECO_LLLL_MASS);
   fChain->SetBranchAddress("RECO_LLLL_PT", &RECO_LLLL_PT, &b_RECO_LLLL_PT);
   fChain->SetBranchAddress("RECO_LLLL_ETA", &RECO_LLLL_ETA, &b_RECO_LLLL_ETA);
   fChain->SetBranchAddress("RECO_LLLL_PHI", &RECO_LLLL_PHI, &b_RECO_LLLL_PHI);
   fChain->SetBranchAddress("RECOzMuMu_MatchingMCTruth", &RECOzMuMu_MatchingMCTruth, &b_RECOzMuMu_MatchingMCTruth);
   fChain->SetBranchAddress("RECOzMuMu_MatchingMCpT", &RECOzMuMu_MatchingMCpT, &b_RECOzMuMu_MatchingMCpT);
   fChain->SetBranchAddress("RECOzMuMu_MatchingMCmass", &RECOzMuMu_MatchingMCmass, &b_RECOzMuMu_MatchingMCmass);
   fChain->SetBranchAddress("RECOzMuMu_MatchingMCEta", &RECOzMuMu_MatchingMCEta, &b_RECOzMuMu_MatchingMCEta);
   fChain->SetBranchAddress("RECOzMuMu_MatchingMCPhi", &RECOzMuMu_MatchingMCPhi, &b_RECOzMuMu_MatchingMCPhi);
   fChain->SetBranchAddress("RECOzEE_MatchingMCTruth", &RECOzEE_MatchingMCTruth, &b_RECOzEE_MatchingMCTruth);
   fChain->SetBranchAddress("RECOzEE_MatchingMCpT", &RECOzEE_MatchingMCpT, &b_RECOzEE_MatchingMCpT);
   fChain->SetBranchAddress("RECOzEE_MatchingMCmass", &RECOzEE_MatchingMCmass, &b_RECOzEE_MatchingMCmass);
   fChain->SetBranchAddress("RECOzEE_MatchingMCEta", &RECOzEE_MatchingMCEta, &b_RECOzEE_MatchingMCEta);
   fChain->SetBranchAddress("RECOzEE_MatchingMCPhi", &RECOzEE_MatchingMCPhi, &b_RECOzEE_MatchingMCPhi);
   fChain->SetBranchAddress("RECOHzzMMMM_MatchingMCTruth", &RECOHzzMMMM_MatchingMCTruth, &b_RECOHzzMMMM_MatchingMCTruth);
   fChain->SetBranchAddress("RECOHzzMMMM_MatchingMCpT", &RECOHzzMMMM_MatchingMCpT, &b_RECOHzzMMMM_MatchingMCpT);
   fChain->SetBranchAddress("RECOHzzMMMM_MatchingMCmass", &RECOHzzMMMM_MatchingMCmass, &b_RECOHzzMMMM_MatchingMCmass);
   fChain->SetBranchAddress("RECOHzzMMMM_MatchingMCEta", &RECOHzzMMMM_MatchingMCEta, &b_RECOHzzMMMM_MatchingMCEta);
   fChain->SetBranchAddress("RECOHzzMMMM_MatchingMCPhi", &RECOHzzMMMM_MatchingMCPhi, &b_RECOHzzMMMM_MatchingMCPhi);
   fChain->SetBranchAddress("RECOHzzEEEE_MatchingMCTruth", &RECOHzzEEEE_MatchingMCTruth, &b_RECOHzzEEEE_MatchingMCTruth);
   fChain->SetBranchAddress("RECOHzzEEEE_MatchingMCpT", &RECOHzzEEEE_MatchingMCpT, &b_RECOHzzEEEE_MatchingMCpT);
   fChain->SetBranchAddress("RECOHzzEEEE_MatchingMCmass", &RECOHzzEEEE_MatchingMCmass, &b_RECOHzzEEEE_MatchingMCmass);
   fChain->SetBranchAddress("RECOHzzEEEE_MatchingMCEta", &RECOHzzEEEE_MatchingMCEta, &b_RECOHzzEEEE_MatchingMCEta);
   fChain->SetBranchAddress("RECOHzzEEEE_MatchingMCPhi", &RECOHzzEEEE_MatchingMCPhi, &b_RECOHzzEEEE_MatchingMCPhi);
   fChain->SetBranchAddress("RECOHzzEEMM_MatchingMCTruth", &RECOHzzEEMM_MatchingMCTruth, &b_RECOHzzEEMM_MatchingMCTruth);
   fChain->SetBranchAddress("RECOHzzEEMM_MatchingMCpT", &RECOHzzEEMM_MatchingMCpT, &b_RECOHzzEEMM_MatchingMCpT);
   fChain->SetBranchAddress("RECOHzzEEMM_MatchingMCmass", &RECOHzzEEMM_MatchingMCmass, &b_RECOHzzEEMM_MatchingMCmass);
   fChain->SetBranchAddress("RECOHzzEEMM_MatchingMCEta", &RECOHzzEEMM_MatchingMCEta, &b_RECOHzzEEMM_MatchingMCEta);
   fChain->SetBranchAddress("RECOHzzEEMM_MatchingMCPhi", &RECOHzzEEMM_MatchingMCPhi, &b_RECOHzzEEMM_MatchingMCPhi);
   fChain->SetBranchAddress("num_PU_vertices", &num_PU_vertices, &b_num_PU_vertices);
   fChain->SetBranchAddress("PU_BunchCrossing", &PU_BunchCrossing, &b_PU_BunchCrossing);
   fChain->SetBranchAddress("MC_weighting", &MC_weighting, &b_MC_weighting);
   fChain->SetBranchAddress("StdFitVertexX", &StdFitVertexX, &b_StdFitVertexX);
   fChain->SetBranchAddress("StdFitVertexY", &StdFitVertexY, &b_StdFitVertexY);
   fChain->SetBranchAddress("StdFitVertexZ", &StdFitVertexZ, &b_StdFitVertexZ);
   fChain->SetBranchAddress("StdFitVertexChi2r", &StdFitVertexChi2r, &b_StdFitVertexChi2r);
   fChain->SetBranchAddress("StdFitVertexProb", &StdFitVertexProb, &b_StdFitVertexProb);
   fChain->SetBranchAddress("StdFitVertexTrack_PT", &StdFitVertexTrack_PT, &b_StdFitVertexTrack_PT);
   fChain->SetBranchAddress("StdFitVertexTrack_ETA", &StdFitVertexTrack_ETA, &b_StdFitVertexTrack_ETA);
   fChain->SetBranchAddress("StdFitVertexTrack_PHI", &StdFitVertexTrack_PHI, &b_StdFitVertexTrack_PHI);
   fChain->SetBranchAddress("KinFitVertexX", &KinFitVertexX, &b_KinFitVertexX);
   fChain->SetBranchAddress("KinFitVertexY", &KinFitVertexY, &b_KinFitVertexY);
   fChain->SetBranchAddress("KinFitVertexZ", &KinFitVertexZ, &b_KinFitVertexZ);
   fChain->SetBranchAddress("KinFitVertexChi2r", &KinFitVertexChi2r, &b_KinFitVertexChi2r);
   fChain->SetBranchAddress("KinFitVertexProb", &KinFitVertexProb, &b_KinFitVertexProb);
   fChain->SetBranchAddress("StdFitVertexXMMMM", &StdFitVertexXMMMM, &b_StdFitVertexXMMMM);
   fChain->SetBranchAddress("StdFitVertexYMMMM", &StdFitVertexYMMMM, &b_StdFitVertexYMMMM);
   fChain->SetBranchAddress("StdFitVertexZMMMM", &StdFitVertexZMMMM, &b_StdFitVertexZMMMM);
   fChain->SetBranchAddress("StdFitVertexChi2rMMMM", &StdFitVertexChi2rMMMM, &b_StdFitVertexChi2rMMMM);
   fChain->SetBranchAddress("StdFitVertexProbMMMM", &StdFitVertexProbMMMM, &b_StdFitVertexProbMMMM);
   fChain->SetBranchAddress("StdFitVertexTrackMMMM_PT", &StdFitVertexTrackMMMM_PT, &b_StdFitVertexTrackMMMM_PT);
   fChain->SetBranchAddress("StdFitVertexTrackMMMM_ETA", &StdFitVertexTrackMMMM_ETA, &b_StdFitVertexTrackMMMM_ETA);
   fChain->SetBranchAddress("StdFitVertexTrackMMMM_PHI", &StdFitVertexTrackMMMM_PHI, &b_StdFitVertexTrackMMMM_PHI);
   fChain->SetBranchAddress("KinFitVertexXMMMM", &KinFitVertexXMMMM, &b_KinFitVertexXMMMM);
   fChain->SetBranchAddress("KinFitVertexYMMMM", &KinFitVertexYMMMM, &b_KinFitVertexYMMMM);
   fChain->SetBranchAddress("KinFitVertexZMMMM", &KinFitVertexZMMMM, &b_KinFitVertexZMMMM);
   fChain->SetBranchAddress("KinFitVertexChi2rMMMM", &KinFitVertexChi2rMMMM, &b_KinFitVertexChi2rMMMM);
   fChain->SetBranchAddress("KinFitVertexProbMMMM", &KinFitVertexProbMMMM, &b_KinFitVertexProbMMMM);
   fChain->SetBranchAddress("StdFitVertexXEEEE", &StdFitVertexXEEEE, &b_StdFitVertexXEEEE);
   fChain->SetBranchAddress("StdFitVertexYEEEE", &StdFitVertexYEEEE, &b_StdFitVertexYEEEE);
   fChain->SetBranchAddress("StdFitVertexZEEEE", &StdFitVertexZEEEE, &b_StdFitVertexZEEEE);
   fChain->SetBranchAddress("StdFitVertexChi2rEEEE", &StdFitVertexChi2rEEEE, &b_StdFitVertexChi2rEEEE);
   fChain->SetBranchAddress("StdFitVertexProbEEEE", &StdFitVertexProbEEEE, &b_StdFitVertexProbEEEE);
   fChain->SetBranchAddress("StdFitVertexTrackEEEE_PT", &StdFitVertexTrackEEEE_PT, &b_StdFitVertexTrackEEEE_PT);
   fChain->SetBranchAddress("StdFitVertexTrackEEEE_ETA", &StdFitVertexTrackEEEE_ETA, &b_StdFitVertexTrackEEEE_ETA);
   fChain->SetBranchAddress("StdFitVertexTrackEEEE_PHI", &StdFitVertexTrackEEEE_PHI, &b_StdFitVertexTrackEEEE_PHI);
   fChain->SetBranchAddress("KinFitVertexXEEEE", &KinFitVertexXEEEE, &b_KinFitVertexXEEEE);
   fChain->SetBranchAddress("KinFitVertexYEEEE", &KinFitVertexYEEEE, &b_KinFitVertexYEEEE);
   fChain->SetBranchAddress("KinFitVertexZEEEE", &KinFitVertexZEEEE, &b_KinFitVertexZEEEE);
   fChain->SetBranchAddress("KinFitVertexChi2rEEEE", &KinFitVertexChi2rEEEE, &b_KinFitVertexChi2rEEEE);
   fChain->SetBranchAddress("KinFitVertexProbEEEE", &KinFitVertexProbEEEE, &b_KinFitVertexProbEEEE);
   fChain->SetBranchAddress("StdFitVertexChi2rDiLep", &StdFitVertexChi2rDiLep, &b_StdFitVertexChi2rDiLep);
   fChain->SetBranchAddress("StdFitVertexProbDiLep", &StdFitVertexProbDiLep, &b_StdFitVertexProbDiLep);
   fChain->SetBranchAddress("StdFitVertexChi2rMMM", &StdFitVertexChi2rMMM, &b_StdFitVertexChi2rMMM);
   fChain->SetBranchAddress("StdFitVertexProbMMM", &StdFitVertexProbMMM, &b_StdFitVertexProbMMM);
   fChain->SetBranchAddress("StdFitVertexChi2rMME", &StdFitVertexChi2rMME, &b_StdFitVertexChi2rMME);
   fChain->SetBranchAddress("StdFitVertexProbMME", &StdFitVertexProbMME, &b_StdFitVertexProbMME);
   fChain->SetBranchAddress("StdFitVertexChi2rEEE", &StdFitVertexChi2rEEE, &b_StdFitVertexChi2rEEE);
   fChain->SetBranchAddress("StdFitVertexProbEEE", &StdFitVertexProbEEE, &b_StdFitVertexProbEEE);
   fChain->SetBranchAddress("StdFitVertexChi2rMEE", &StdFitVertexChi2rMEE, &b_StdFitVertexChi2rMEE);
   fChain->SetBranchAddress("StdFitVertexProbMEE", &StdFitVertexProbMEE, &b_StdFitVertexProbMEE);
   fChain->SetBranchAddress("ftsigma", &ftsigma, &b_ftsigma);
   fChain->SetBranchAddress("gdX", &gdX, &b_gdX);
   fChain->SetBranchAddress("gdY", &gdY, &b_gdY);
   fChain->SetBranchAddress("gdZ", &gdZ, &b_gdZ);
   fChain->SetBranchAddress("ftsigmalag", &ftsigmalag, &b_ftsigmalag);
   fChain->SetBranchAddress("gdlagX", &gdlagX, &b_gdlagX);
   fChain->SetBranchAddress("gdlagY", &gdlagY, &b_gdlagY);
   fChain->SetBranchAddress("gdlagZ", &gdlagZ, &b_gdlagZ);
   fChain->SetBranchAddress("gdlagProb", &gdlagProb, &b_gdlagProb);
   fChain->SetBranchAddress("gdlagNdof", &gdlagNdof, &b_gdlagNdof);
   fChain->SetBranchAddress("ftsigmaMMMM", &ftsigmaMMMM, &b_ftsigmaMMMM);
   fChain->SetBranchAddress("gdXMMMM", &gdXMMMM, &b_gdXMMMM);
   fChain->SetBranchAddress("gdYMMMM", &gdYMMMM, &b_gdYMMMM);
   fChain->SetBranchAddress("gdZMMMM", &gdZMMMM, &b_gdZMMMM);
   fChain->SetBranchAddress("ftsigmalagMMMM", &ftsigmalagMMMM, &b_ftsigmalagMMMM);
   fChain->SetBranchAddress("gdlagXMMMM", &gdlagXMMMM, &b_gdlagXMMMM);
   fChain->SetBranchAddress("gdlagYMMMM", &gdlagYMMMM, &b_gdlagYMMMM);
   fChain->SetBranchAddress("gdlagZMMMM", &gdlagZMMMM, &b_gdlagZMMMM);
   fChain->SetBranchAddress("gdlagProbMMMM", &gdlagProbMMMM, &b_gdlagProbMMMM);
   fChain->SetBranchAddress("gdlagNdofMMMM", &gdlagNdofMMMM, &b_gdlagNdofMMMM);
   fChain->SetBranchAddress("ftsigmaEEEE", &ftsigmaEEEE, &b_ftsigmaEEEE);
   fChain->SetBranchAddress("gdXEEEE", &gdXEEEE, &b_gdXEEEE);
   fChain->SetBranchAddress("gdYEEEE", &gdYEEEE, &b_gdYEEEE);
   fChain->SetBranchAddress("gdZEEEE", &gdZEEEE, &b_gdZEEEE);
   fChain->SetBranchAddress("ftsigmalagEEEE", &ftsigmalagEEEE, &b_ftsigmalagEEEE);
   fChain->SetBranchAddress("gdlagXEEEE", &gdlagXEEEE, &b_gdlagXEEEE);
   fChain->SetBranchAddress("gdlagYEEEE", &gdlagYEEEE, &b_gdlagYEEEE);
   fChain->SetBranchAddress("gdlagZEEEE", &gdlagZEEEE, &b_gdlagZEEEE);
   fChain->SetBranchAddress("gdlagProbEEEE", &gdlagProbEEEE, &b_gdlagProbEEEE);
   fChain->SetBranchAddress("gdlagNdofEEEE", &gdlagNdofEEEE, &b_gdlagNdofEEEE);
   fChain->SetBranchAddress("RECORF_2e2mu_cosTheta1_spin", &RECORF_2e2mu_cosTheta1_spin, &b_RECORF_2e2mu_cosTheta1_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_cosTheta2_spin", &RECORF_2e2mu_cosTheta2_spin, &b_RECORF_2e2mu_cosTheta2_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_cosThetaStar_spin", &RECORF_2e2mu_cosThetaStar_spin, &b_RECORF_2e2mu_cosThetaStar_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_Phi_spin", &RECORF_2e2mu_Phi_spin, &b_RECORF_2e2mu_Phi_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_Phi1_spin", &RECORF_2e2mu_Phi1_spin, &b_RECORF_2e2mu_Phi1_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_Phi2_spin", &RECORF_2e2mu_Phi2_spin, &b_RECORF_2e2mu_Phi2_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_phi1RF_spin", &RECORF_2e2mu_phi1RF_spin, &b_RECORF_2e2mu_phi1RF_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_phi2RF_spin", &RECORF_2e2mu_phi2RF_spin, &b_RECORF_2e2mu_phi2RF_spin);
   fChain->SetBranchAddress("RECORF_2e2mu_MELA", &RECORF_2e2mu_MELA, &b_RECORF_2e2mu_MELA);
   fChain->SetBranchAddress("RECORF_4e_cosTheta1_spin", &RECORF_4e_cosTheta1_spin, &b_RECORF_4e_cosTheta1_spin);
   fChain->SetBranchAddress("RECORF_4e_cosTheta2_spin", &RECORF_4e_cosTheta2_spin, &b_RECORF_4e_cosTheta2_spin);
   fChain->SetBranchAddress("RECORF_4e_cosThetaStar_spin", &RECORF_4e_cosThetaStar_spin, &b_RECORF_4e_cosThetaStar_spin);
   fChain->SetBranchAddress("RECORF_4e_Phi_spin", &RECORF_4e_Phi_spin, &b_RECORF_4e_Phi_spin);
   fChain->SetBranchAddress("RECORF_4e_Phi1_spin", &RECORF_4e_Phi1_spin, &b_RECORF_4e_Phi1_spin);
   fChain->SetBranchAddress("RECORF_4e_Phi2_spin", &RECORF_4e_Phi2_spin, &b_RECORF_4e_Phi2_spin);
   fChain->SetBranchAddress("RECORF_4e_phi1RF_spin", &RECORF_4e_phi1RF_spin, &b_RECORF_4e_phi1RF_spin);
   fChain->SetBranchAddress("RECORF_4e_phi2RF_spin", &RECORF_4e_phi2RF_spin, &b_RECORF_4e_phi2RF_spin);
   fChain->SetBranchAddress("RECORF_4e_MELA", &RECORF_4e_MELA, &b_RECORF_4e_MELA);
   fChain->SetBranchAddress("RECORF_4mu_cosTheta1_spin", &RECORF_4mu_cosTheta1_spin, &b_RECORF_4mu_cosTheta1_spin);
   fChain->SetBranchAddress("RECORF_4mu_cosTheta2_spin", &RECORF_4mu_cosTheta2_spin, &b_RECORF_4mu_cosTheta2_spin);
   fChain->SetBranchAddress("RECORF_4mu_cosThetaStar_spin", &RECORF_4mu_cosThetaStar_spin, &b_RECORF_4mu_cosThetaStar_spin);
   fChain->SetBranchAddress("RECORF_4mu_Phi_spin", &RECORF_4mu_Phi_spin, &b_RECORF_4mu_Phi_spin);
   fChain->SetBranchAddress("RECORF_4mu_Phi1_spin", &RECORF_4mu_Phi1_spin, &b_RECORF_4mu_Phi1_spin);
   fChain->SetBranchAddress("RECORF_4mu_Phi2_spin", &RECORF_4mu_Phi2_spin, &b_RECORF_4mu_Phi2_spin);
   fChain->SetBranchAddress("RECORF_4mu_phi1RF_spin", &RECORF_4mu_phi1RF_spin, &b_RECORF_4mu_phi1RF_spin);
   fChain->SetBranchAddress("RECORF_4mu_phi2RF_spin", &RECORF_4mu_phi2RF_spin, &b_RECORF_4mu_phi2RF_spin);
   fChain->SetBranchAddress("RECORF_4mu_MELA", &RECORF_4mu_MELA, &b_RECORF_4mu_MELA);
   fChain->SetBranchAddress("BeamSpot_X", &BeamSpot_X, &b_BeamSpot_X);
   fChain->SetBranchAddress("BeamSpot_Y", &BeamSpot_Y, &b_BeamSpot_Y);
   fChain->SetBranchAddress("BeamSpot_Z", &BeamSpot_Z, &b_BeamSpot_Z);
   fChain->SetBranchAddress("nbPv", &nbPv, &b_nbPv);
   fChain->SetBranchAddress("Nbdof", &Nbdof, &b_Nbdof);
   fChain->SetBranchAddress("PositionRho", &PositionRho, &b_PositionRho);
   fChain->SetBranchAddress("PositionX", &PositionX, &b_PositionX);
   fChain->SetBranchAddress("PositionY", &PositionY, &b_PositionY);
   fChain->SetBranchAddress("PositionZ", &PositionZ, &b_PositionZ);
   fChain->SetBranchAddress("RECO_PFJET_N", &RECO_PFJET_N, &b_RECO_PFJET_N);
   fChain->SetBranchAddress("RECO_PFJET_CHARGE", &RECO_PFJET_CHARGE, &b_RECO_PFJET_CHARGE);
   fChain->SetBranchAddress("RECO_PFJET_ET", &RECO_PFJET_ET, &b_RECO_PFJET_ET);
   fChain->SetBranchAddress("RECO_PFJET_PT", &RECO_PFJET_PT, &b_RECO_PFJET_PT);
   fChain->SetBranchAddress("RECO_PFJET_ETA", &RECO_PFJET_ETA, &b_RECO_PFJET_ETA);
   fChain->SetBranchAddress("RECO_PFJET_PHI", &RECO_PFJET_PHI, &b_RECO_PFJET_PHI);
   fChain->SetBranchAddress("RECO_PFJET_PUID", &RECO_PFJET_PUID, &b_RECO_PFJET_PUID);
   fChain->SetBranchAddress("RECO_PFJET_PUID_MVA", &RECO_PFJET_PUID_MVA, &b_RECO_PFJET_PUID_MVA);
   fChain->SetBranchAddress("RHO_ele", &RHO_ele, &b_RHO_ele);
   fChain->SetBranchAddress("RHO_mu", &RHO_mu, &b_RHO_mu);
   Notify();
}

Bool_t MonoHiggsAnalysis4mu::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MonoHiggsAnalysis4mu::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MonoHiggsAnalysis4mu::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

/* void  MonoHiggsAnalysis4mu::printmubnn(int i){     */
/*   bnn_file  */
/*                 << RECOMU_PT->at(i)  << " "   */
/*                 << RECOMU_ETA->at(i) << " "   */
/*                 << RECOMU_PHI->at(i) << " "   */
/*                 << RECOMU_CHARGE->at(i) << " " */
/*                 << RECOMU_PFX_dB->at(i) << " " */
/*                 << RECOMU_SIP->at(i) << " "; */
/* } */
#endif // #ifdef MonoHiggsAnalysis4mu_cxx
