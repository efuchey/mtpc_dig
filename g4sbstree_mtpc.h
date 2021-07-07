//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 15 12:12:06 2021 by ROOT version 6.14/04
// from TTree T/Geant4 SBS Simulation
// found on file: g4sbs_out_protons_60-400MeV_30-70degrees_1kevents.root
//////////////////////////////////////////////////////////

#ifndef g4sbstree_mtpc_h
#define g4sbstree_mtpc_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class g4sbstree_mtpc {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        ev_count;
   Double_t        ev_rate;
   Double_t        ev_solang;
   Double_t        ev_sigma;
   Double_t        ev_W2;
   Double_t        ev_xbj;
   Double_t        ev_Q2;
   Double_t        ev_th;
   Double_t        ev_ph;
   Double_t        ev_sigmaDIS;
   Double_t        ev_sigmaTDIS;
   Double_t        ev_Aperp;
   Double_t        ev_Apar;
   Double_t        ev_Pt;
   Double_t        ev_Pl;
   Double_t        ev_vx;
   Double_t        ev_vy;
   Double_t        ev_vz;
   Double_t        ev_ep;
   Double_t        ev_np;
   Double_t        ev_p1p;
   Double_t        ev_p2p;
   Double_t        ev_pip;
   Double_t        ev_epx;
   Double_t        ev_epy;
   Double_t        ev_epz;
   Double_t        ev_npx;
   Double_t        ev_npy;
   Double_t        ev_npz;
   Double_t        ev_p1px;
   Double_t        ev_p1py;
   Double_t        ev_p1pz;
   Double_t        ev_p2px;
   Double_t        ev_p2py;
   Double_t        ev_p2pz;
   Double_t        ev_pipx;
   Double_t        ev_pipy;
   Double_t        ev_pipz;
   Double_t        ev_nth;
   Double_t        ev_nph;
   Double_t        ev_p1th;
   Double_t        ev_p1ph;
   Double_t        ev_p2th;
   Double_t        ev_p2ph;
   Double_t        ev_pith;
   Double_t        ev_piph;
   Double_t        ev_pmperp;
   Double_t        ev_pmpar;
   Double_t        ev_pmparsm;
   Double_t        ev_z;
   Double_t        ev_phperp;
   Double_t        ev_phih;
   Double_t        ev_MX2;
   Double_t        ev_xpi;
   Double_t        ev_tpi;
   Double_t        ev_xa;
   Double_t        ev_pt;
   Double_t        ev_nu;
   Double_t        ev_ya;
   Double_t        ev_y;
   Double_t        ev_f2p;
   Double_t        ev_f2pi;
   Double_t        ev_ypi;
   Double_t        ev_Sx;
   Double_t        ev_Sy;
   Double_t        ev_Sz;
   Int_t           ev_nucl;
   Int_t           ev_fnucl;
   Int_t           ev_hadr;
   Int_t           ev_earmaccept;
   Int_t           ev_harmaccept;
   Double_t        Harm_LACScint_det_esum;
   Int_t           Harm_LACScint_hit_nhits;
   vector<int>     *Harm_LACScint_hit_row;
   vector<int>     *Harm_LACScint_hit_col;
   vector<int>     *Harm_LACScint_hit_cell;
   vector<int>     *Harm_LACScint_hit_plane;
   vector<int>     *Harm_LACScint_hit_wire;
   vector<double>  *Harm_LACScint_hit_xcell;
   vector<double>  *Harm_LACScint_hit_ycell;
   vector<double>  *Harm_LACScint_hit_zcell;
   vector<double>  *Harm_LACScint_hit_xcellg;
   vector<double>  *Harm_LACScint_hit_ycellg;
   vector<double>  *Harm_LACScint_hit_zcellg;
   vector<double>  *Harm_LACScint_hit_xhit;
   vector<double>  *Harm_LACScint_hit_yhit;
   vector<double>  *Harm_LACScint_hit_zhit;
   vector<double>  *Harm_LACScint_hit_xhitg;
   vector<double>  *Harm_LACScint_hit_yhitg;
   vector<double>  *Harm_LACScint_hit_zhitg;
   vector<double>  *Harm_LACScint_hit_sumedep;
   vector<double>  *Harm_LACScint_hit_tavg;
   vector<double>  *Harm_LACScint_hit_trms;
   vector<double>  *Harm_LACScint_hit_tmin;
   vector<double>  *Harm_LACScint_hit_tmax;
   Int_t           Harm_RICH_hit_nhits;
   vector<int>     *Harm_RICH_hit_PMT;
   vector<int>     *Harm_RICH_hit_row;
   vector<int>     *Harm_RICH_hit_col;
   vector<double>  *Harm_RICH_hit_xpmt;
   vector<double>  *Harm_RICH_hit_ypmt;
   vector<double>  *Harm_RICH_hit_zpmt;
   vector<double>  *Harm_RICH_hit_xgpmt;
   vector<double>  *Harm_RICH_hit_ygpmt;
   vector<double>  *Harm_RICH_hit_zgpmt;
   vector<int>     *Harm_RICH_hit_NumPhotoelectrons;
   vector<double>  *Harm_RICH_hit_Time_avg;
   vector<double>  *Harm_RICH_hit_Time_rms;
   vector<double>  *Harm_RICH_hit_Time_min;
   vector<double>  *Harm_RICH_hit_Time_max;
   vector<int>     *Harm_RICH_hit_mTrackNo;
   vector<double>  *Harm_RICH_hit_xhit;
   vector<double>  *Harm_RICH_hit_yhit;
   vector<double>  *Harm_RICH_hit_zhit;
   vector<double>  *Harm_RICH_hit_pxhit;
   vector<double>  *Harm_RICH_hit_pyhit;
   vector<double>  *Harm_RICH_hit_pzhit;
   vector<double>  *Harm_RICH_hit_pvx;
   vector<double>  *Harm_RICH_hit_pvy;
   vector<double>  *Harm_RICH_hit_pvz;
   vector<double>  *Harm_RICH_hit_ppx;
   vector<double>  *Harm_RICH_hit_ppy;
   vector<double>  *Harm_RICH_hit_ppz;
   vector<int>     *Harm_RICH_hit_volume_flag;
   Int_t           Harm_SBSGEM_hit_nhits;
   vector<int>     *Harm_SBSGEM_hit_plane;
   vector<int>     *Harm_SBSGEM_hit_strip;
   vector<double>  *Harm_SBSGEM_hit_x;
   vector<double>  *Harm_SBSGEM_hit_y;
   vector<double>  *Harm_SBSGEM_hit_z;
   vector<double>  *Harm_SBSGEM_hit_polx;
   vector<double>  *Harm_SBSGEM_hit_poly;
   vector<double>  *Harm_SBSGEM_hit_polz;
   vector<double>  *Harm_SBSGEM_hit_t;
   vector<double>  *Harm_SBSGEM_hit_trms;
   vector<double>  *Harm_SBSGEM_hit_tmin;
   vector<double>  *Harm_SBSGEM_hit_tmax;
   vector<double>  *Harm_SBSGEM_hit_tx;
   vector<double>  *Harm_SBSGEM_hit_ty;
   vector<double>  *Harm_SBSGEM_hit_xin;
   vector<double>  *Harm_SBSGEM_hit_yin;
   vector<double>  *Harm_SBSGEM_hit_zin;
   vector<double>  *Harm_SBSGEM_hit_xout;
   vector<double>  *Harm_SBSGEM_hit_yout;
   vector<double>  *Harm_SBSGEM_hit_zout;
   vector<double>  *Harm_SBSGEM_hit_txp;
   vector<double>  *Harm_SBSGEM_hit_typ;
   vector<double>  *Harm_SBSGEM_hit_xg;
   vector<double>  *Harm_SBSGEM_hit_yg;
   vector<double>  *Harm_SBSGEM_hit_zg;
   vector<int>     *Harm_SBSGEM_hit_trid;
   vector<int>     *Harm_SBSGEM_hit_mid;
   vector<int>     *Harm_SBSGEM_hit_pid;
   vector<double>  *Harm_SBSGEM_hit_vx;
   vector<double>  *Harm_SBSGEM_hit_vy;
   vector<double>  *Harm_SBSGEM_hit_vz;
   vector<double>  *Harm_SBSGEM_hit_p;
   vector<double>  *Harm_SBSGEM_hit_edep;
   vector<double>  *Harm_SBSGEM_hit_beta;
   Int_t           Harm_SBSGEM_Track_ntracks;
   vector<int>     *Harm_SBSGEM_Track_TID;
   vector<int>     *Harm_SBSGEM_Track_PID;
   vector<int>     *Harm_SBSGEM_Track_MID;
   vector<int>     *Harm_SBSGEM_Track_NumHits;
   vector<int>     *Harm_SBSGEM_Track_NumPlanes;
   vector<int>     *Harm_SBSGEM_Track_NDF;
   vector<double>  *Harm_SBSGEM_Track_Chi2fit;
   vector<double>  *Harm_SBSGEM_Track_Chi2true;
   vector<double>  *Harm_SBSGEM_Track_X;
   vector<double>  *Harm_SBSGEM_Track_Y;
   vector<double>  *Harm_SBSGEM_Track_Xp;
   vector<double>  *Harm_SBSGEM_Track_Yp;
   vector<double>  *Harm_SBSGEM_Track_T;
   vector<double>  *Harm_SBSGEM_Track_P;
   vector<double>  *Harm_SBSGEM_Track_Sx;
   vector<double>  *Harm_SBSGEM_Track_Sy;
   vector<double>  *Harm_SBSGEM_Track_Sz;
   vector<double>  *Harm_SBSGEM_Track_Xfit;
   vector<double>  *Harm_SBSGEM_Track_Yfit;
   vector<double>  *Harm_SBSGEM_Track_Xpfit;
   vector<double>  *Harm_SBSGEM_Track_Ypfit;
   Double_t        SBS_mTPC_det_esum;
   Int_t           SBS_mTPC_hit_nhits;
   vector<int>     *SBS_mTPC_hit_cell;
   vector<double>  *SBS_mTPC_hit_xcell;
   vector<double>  *SBS_mTPC_hit_ycell;
   vector<double>  *SBS_mTPC_hit_zcell;
   vector<double>  *SBS_mTPC_hit_xcellg;
   vector<double>  *SBS_mTPC_hit_ycellg;
   vector<double>  *SBS_mTPC_hit_zcellg;
   vector<double>  *SBS_mTPC_hit_xhit;
   vector<double>  *SBS_mTPC_hit_yhit;
   vector<double>  *SBS_mTPC_hit_zhit;
   vector<double>  *SBS_mTPC_hit_xhitg;
   vector<double>  *SBS_mTPC_hit_yhitg;
   vector<double>  *SBS_mTPC_hit_zhitg;
   vector<double>  *SBS_mTPC_hit_sumedep;
   vector<double>  *SBS_mTPC_hit_tavg;
   vector<double>  *SBS_mTPC_hit_trms;
   vector<double>  *SBS_mTPC_hit_tmin;
   vector<double>  *SBS_mTPC_hit_tmax;
   Int_t           SBS_mTPC_npart_mTPC;
   vector<int>     *SBS_mTPC_ihit;
   vector<double>  *SBS_mTPC_x;
   vector<double>  *SBS_mTPC_y;
   vector<double>  *SBS_mTPC_z;
   vector<double>  *SBS_mTPC_t;
   vector<double>  *SBS_mTPC_E;
   vector<double>  *SBS_mTPC_dt;
   vector<double>  *SBS_mTPC_L;
   vector<double>  *SBS_mTPC_vx;
   vector<double>  *SBS_mTPC_vy;
   vector<double>  *SBS_mTPC_vz;
   vector<int>     *SBS_mTPC_trid;
   vector<int>     *SBS_mTPC_mid;
   vector<int>     *SBS_mTPC_pid;
   vector<double>  *SBS_mTPC_p;
   vector<double>  *SBS_mTPC_px;
   vector<double>  *SBS_mTPC_py;
   vector<double>  *SBS_mTPC_pz;
   vector<double>  *SBS_mTPC_edep;
   vector<double>  *SBS_mTPC_ztravel;
   vector<int>     *SBS_mTPC_nstrips;
   Int_t           SBS_mTPC_part_npart;
   vector<int>     *SBS_mTPC_part_PID;
   vector<int>     *SBS_mTPC_part_MID;
   vector<int>     *SBS_mTPC_part_TID;
   vector<int>     *SBS_mTPC_part_nbounce;
   vector<int>     *SBS_mTPC_part_hitindex;
   vector<double>  *SBS_mTPC_part_vx;
   vector<double>  *SBS_mTPC_part_vy;
   vector<double>  *SBS_mTPC_part_vz;
   vector<double>  *SBS_mTPC_part_px;
   vector<double>  *SBS_mTPC_part_py;
   vector<double>  *SBS_mTPC_part_pz;
   Double_t        SBS_mTPCGEM_det_esum;
   Int_t           SBS_mTPCGEM_hit_nhits;
   vector<int>     *SBS_mTPCGEM_hit_row;
   vector<int>     *SBS_mTPCGEM_hit_col;
   vector<int>     *SBS_mTPCGEM_hit_cell;
   vector<int>     *SBS_mTPCGEM_hit_plane;
   vector<int>     *SBS_mTPCGEM_hit_wire;
   vector<double>  *SBS_mTPCGEM_hit_xcell;
   vector<double>  *SBS_mTPCGEM_hit_ycell;
   vector<double>  *SBS_mTPCGEM_hit_zcell;
   vector<double>  *SBS_mTPCGEM_hit_xcellg;
   vector<double>  *SBS_mTPCGEM_hit_ycellg;
   vector<double>  *SBS_mTPCGEM_hit_zcellg;
   vector<double>  *SBS_mTPCGEM_hit_xhit;
   vector<double>  *SBS_mTPCGEM_hit_yhit;
   vector<double>  *SBS_mTPCGEM_hit_zhit;
   vector<double>  *SBS_mTPCGEM_hit_xhitg;
   vector<double>  *SBS_mTPCGEM_hit_yhitg;
   vector<double>  *SBS_mTPCGEM_hit_zhitg;
   vector<double>  *SBS_mTPCGEM_hit_sumedep;
   vector<double>  *SBS_mTPCGEM_hit_tavg;
   vector<double>  *SBS_mTPCGEM_hit_trms;
   vector<double>  *SBS_mTPCGEM_hit_tmin;
   vector<double>  *SBS_mTPCGEM_hit_tmax;
   Double_t        SBS_mTPCHV_det_esum;
   Int_t           SBS_mTPCHV_hit_nhits;
   vector<int>     *SBS_mTPCHV_hit_row;
   vector<int>     *SBS_mTPCHV_hit_col;
   vector<int>     *SBS_mTPCHV_hit_cell;
   vector<int>     *SBS_mTPCHV_hit_plane;
   vector<int>     *SBS_mTPCHV_hit_wire;
   vector<double>  *SBS_mTPCHV_hit_xcell;
   vector<double>  *SBS_mTPCHV_hit_ycell;
   vector<double>  *SBS_mTPCHV_hit_zcell;
   vector<double>  *SBS_mTPCHV_hit_xcellg;
   vector<double>  *SBS_mTPCHV_hit_ycellg;
   vector<double>  *SBS_mTPCHV_hit_zcellg;
   vector<double>  *SBS_mTPCHV_hit_xhit;
   vector<double>  *SBS_mTPCHV_hit_yhit;
   vector<double>  *SBS_mTPCHV_hit_zhit;
   vector<double>  *SBS_mTPCHV_hit_xhitg;
   vector<double>  *SBS_mTPCHV_hit_yhitg;
   vector<double>  *SBS_mTPCHV_hit_zhitg;
   vector<double>  *SBS_mTPCHV_hit_sumedep;
   vector<double>  *SBS_mTPCHV_hit_tavg;
   vector<double>  *SBS_mTPCHV_hit_trms;
   vector<double>  *SBS_mTPCHV_hit_tmin;
   vector<double>  *SBS_mTPCHV_hit_tmax;
   Double_t        SBS_mTPCReadout_det_esum;
   Int_t           SBS_mTPCReadout_hit_nhits;
   vector<int>     *SBS_mTPCReadout_hit_row;
   vector<int>     *SBS_mTPCReadout_hit_col;
   vector<int>     *SBS_mTPCReadout_hit_cell;
   vector<int>     *SBS_mTPCReadout_hit_plane;
   vector<int>     *SBS_mTPCReadout_hit_wire;
   vector<double>  *SBS_mTPCReadout_hit_xcell;
   vector<double>  *SBS_mTPCReadout_hit_ycell;
   vector<double>  *SBS_mTPCReadout_hit_zcell;
   vector<double>  *SBS_mTPCReadout_hit_xcellg;
   vector<double>  *SBS_mTPCReadout_hit_ycellg;
   vector<double>  *SBS_mTPCReadout_hit_zcellg;
   vector<double>  *SBS_mTPCReadout_hit_xhit;
   vector<double>  *SBS_mTPCReadout_hit_yhit;
   vector<double>  *SBS_mTPCReadout_hit_zhit;
   vector<double>  *SBS_mTPCReadout_hit_xhitg;
   vector<double>  *SBS_mTPCReadout_hit_yhitg;
   vector<double>  *SBS_mTPCReadout_hit_zhitg;
   vector<double>  *SBS_mTPCReadout_hit_sumedep;
   vector<double>  *SBS_mTPCReadout_hit_tavg;
   vector<double>  *SBS_mTPCReadout_hit_trms;
   vector<double>  *SBS_mTPCReadout_hit_tmin;
   vector<double>  *SBS_mTPCReadout_hit_tmax;

   // List of branches
   TBranch        *b_ev;   //!
   TBranch        *b_Harm_LACScint_det_esum;   //!
   TBranch        *b_Harm_LACScint_hit_nhits;   //!
   TBranch        *b_Harm_LACScint_hit_row;   //!
   TBranch        *b_Harm_LACScint_hit_col;   //!
   TBranch        *b_Harm_LACScint_hit_cell;   //!
   TBranch        *b_Harm_LACScint_hit_plane;   //!
   TBranch        *b_Harm_LACScint_hit_wire;   //!
   TBranch        *b_Harm_LACScint_hit_xcell;   //!
   TBranch        *b_Harm_LACScint_hit_ycell;   //!
   TBranch        *b_Harm_LACScint_hit_zcell;   //!
   TBranch        *b_Harm_LACScint_hit_xcellg;   //!
   TBranch        *b_Harm_LACScint_hit_ycellg;   //!
   TBranch        *b_Harm_LACScint_hit_zcellg;   //!
   TBranch        *b_Harm_LACScint_hit_xhit;   //!
   TBranch        *b_Harm_LACScint_hit_yhit;   //!
   TBranch        *b_Harm_LACScint_hit_zhit;   //!
   TBranch        *b_Harm_LACScint_hit_xhitg;   //!
   TBranch        *b_Harm_LACScint_hit_yhitg;   //!
   TBranch        *b_Harm_LACScint_hit_zhitg;   //!
   TBranch        *b_Harm_LACScint_hit_sumedep;   //!
   TBranch        *b_Harm_LACScint_hit_tavg;   //!
   TBranch        *b_Harm_LACScint_hit_trms;   //!
   TBranch        *b_Harm_LACScint_hit_tmin;   //!
   TBranch        *b_Harm_LACScint_hit_tmax;   //!
   TBranch        *b_Harm_RICH_hit_nhits;   //!
   TBranch        *b_Harm_RICH_hit_PMT;   //!
   TBranch        *b_Harm_RICH_hit_row;   //!
   TBranch        *b_Harm_RICH_hit_col;   //!
   TBranch        *b_Harm_RICH_hit_xpmt;   //!
   TBranch        *b_Harm_RICH_hit_ypmt;   //!
   TBranch        *b_Harm_RICH_hit_zpmt;   //!
   TBranch        *b_Harm_RICH_hit_xgpmt;   //!
   TBranch        *b_Harm_RICH_hit_ygpmt;   //!
   TBranch        *b_Harm_RICH_hit_zgpmt;   //!
   TBranch        *b_Harm_RICH_hit_NumPhotoelectrons;   //!
   TBranch        *b_Harm_RICH_hit_Time_avg;   //!
   TBranch        *b_Harm_RICH_hit_Time_rms;   //!
   TBranch        *b_Harm_RICH_hit_Time_min;   //!
   TBranch        *b_Harm_RICH_hit_Time_max;   //!
   TBranch        *b_Harm_RICH_hit_mTrackNo;   //!
   TBranch        *b_Harm_RICH_hit_xhit;   //!
   TBranch        *b_Harm_RICH_hit_yhit;   //!
   TBranch        *b_Harm_RICH_hit_zhit;   //!
   TBranch        *b_Harm_RICH_hit_pxhit;   //!
   TBranch        *b_Harm_RICH_hit_pyhit;   //!
   TBranch        *b_Harm_RICH_hit_pzhit;   //!
   TBranch        *b_Harm_RICH_hit_pvx;   //!
   TBranch        *b_Harm_RICH_hit_pvy;   //!
   TBranch        *b_Harm_RICH_hit_pvz;   //!
   TBranch        *b_Harm_RICH_hit_ppx;   //!
   TBranch        *b_Harm_RICH_hit_ppy;   //!
   TBranch        *b_Harm_RICH_hit_ppz;   //!
   TBranch        *b_Harm_RICH_hit_volume_flag;   //!
   TBranch        *b_Harm_SBSGEM_hit_nhits;   //!
   TBranch        *b_Harm_SBSGEM_hit_plane;   //!
   TBranch        *b_Harm_SBSGEM_hit_strip;   //!
   TBranch        *b_Harm_SBSGEM_hit_x;   //!
   TBranch        *b_Harm_SBSGEM_hit_y;   //!
   TBranch        *b_Harm_SBSGEM_hit_z;   //!
   TBranch        *b_Harm_SBSGEM_hit_polx;   //!
   TBranch        *b_Harm_SBSGEM_hit_poly;   //!
   TBranch        *b_Harm_SBSGEM_hit_polz;   //!
   TBranch        *b_Harm_SBSGEM_hit_t;   //!
   TBranch        *b_Harm_SBSGEM_hit_trms;   //!
   TBranch        *b_Harm_SBSGEM_hit_tmin;   //!
   TBranch        *b_Harm_SBSGEM_hit_tmax;   //!
   TBranch        *b_Harm_SBSGEM_hit_tx;   //!
   TBranch        *b_Harm_SBSGEM_hit_ty;   //!
   TBranch        *b_Harm_SBSGEM_hit_xin;   //!
   TBranch        *b_Harm_SBSGEM_hit_yin;   //!
   TBranch        *b_Harm_SBSGEM_hit_zin;   //!
   TBranch        *b_Harm_SBSGEM_hit_xout;   //!
   TBranch        *b_Harm_SBSGEM_hit_yout;   //!
   TBranch        *b_Harm_SBSGEM_hit_zout;   //!
   TBranch        *b_Harm_SBSGEM_hit_txp;   //!
   TBranch        *b_Harm_SBSGEM_hit_typ;   //!
   TBranch        *b_Harm_SBSGEM_hit_xg;   //!
   TBranch        *b_Harm_SBSGEM_hit_yg;   //!
   TBranch        *b_Harm_SBSGEM_hit_zg;   //!
   TBranch        *b_Harm_SBSGEM_hit_trid;   //!
   TBranch        *b_Harm_SBSGEM_hit_mid;   //!
   TBranch        *b_Harm_SBSGEM_hit_pid;   //!
   TBranch        *b_Harm_SBSGEM_hit_vx;   //!
   TBranch        *b_Harm_SBSGEM_hit_vy;   //!
   TBranch        *b_Harm_SBSGEM_hit_vz;   //!
   TBranch        *b_Harm_SBSGEM_hit_p;   //!
   TBranch        *b_Harm_SBSGEM_hit_edep;   //!
   TBranch        *b_Harm_SBSGEM_hit_beta;   //!
   TBranch        *b_Harm_SBSGEM_Track_ntracks;   //!
   TBranch        *b_Harm_SBSGEM_Track_TID;   //!
   TBranch        *b_Harm_SBSGEM_Track_PID;   //!
   TBranch        *b_Harm_SBSGEM_Track_MID;   //!
   TBranch        *b_Harm_SBSGEM_Track_NumHits;   //!
   TBranch        *b_Harm_SBSGEM_Track_NumPlanes;   //!
   TBranch        *b_Harm_SBSGEM_Track_NDF;   //!
   TBranch        *b_Harm_SBSGEM_Track_Chi2fit;   //!
   TBranch        *b_Harm_SBSGEM_Track_Chi2true;   //!
   TBranch        *b_Harm_SBSGEM_Track_X;   //!
   TBranch        *b_Harm_SBSGEM_Track_Y;   //!
   TBranch        *b_Harm_SBSGEM_Track_Xp;   //!
   TBranch        *b_Harm_SBSGEM_Track_Yp;   //!
   TBranch        *b_Harm_SBSGEM_Track_T;   //!
   TBranch        *b_Harm_SBSGEM_Track_P;   //!
   TBranch        *b_Harm_SBSGEM_Track_Sx;   //!
   TBranch        *b_Harm_SBSGEM_Track_Sy;   //!
   TBranch        *b_Harm_SBSGEM_Track_Sz;   //!
   TBranch        *b_Harm_SBSGEM_Track_Xfit;   //!
   TBranch        *b_Harm_SBSGEM_Track_Yfit;   //!
   TBranch        *b_Harm_SBSGEM_Track_Xpfit;   //!
   TBranch        *b_Harm_SBSGEM_Track_Ypfit;   //!
   TBranch        *b_SBS_mTPC_det_esum;   //!
   TBranch        *b_SBS_mTPC_hit_nhits;   //!
   TBranch        *b_SBS_mTPC_hit_cell;   //!
   TBranch        *b_SBS_mTPC_hit_xcell;   //!
   TBranch        *b_SBS_mTPC_hit_ycell;   //!
   TBranch        *b_SBS_mTPC_hit_zcell;   //!
   TBranch        *b_SBS_mTPC_hit_xcellg;   //!
   TBranch        *b_SBS_mTPC_hit_ycellg;   //!
   TBranch        *b_SBS_mTPC_hit_zcellg;   //!
   TBranch        *b_SBS_mTPC_hit_xhit;   //!
   TBranch        *b_SBS_mTPC_hit_yhit;   //!
   TBranch        *b_SBS_mTPC_hit_zhit;   //!
   TBranch        *b_SBS_mTPC_hit_xhitg;   //!
   TBranch        *b_SBS_mTPC_hit_yhitg;   //!
   TBranch        *b_SBS_mTPC_hit_zhitg;   //!
   TBranch        *b_SBS_mTPC_hit_sumedep;   //!
   TBranch        *b_SBS_mTPC_hit_tavg;   //!
   TBranch        *b_SBS_mTPC_hit_trms;   //!
   TBranch        *b_SBS_mTPC_hit_tmin;   //!
   TBranch        *b_SBS_mTPC_hit_tmax;   //!
   TBranch        *b_SBS_mTPC_npart_mTPC;   //!
   TBranch        *b_SBS_mTPC_ihit;   //!
   TBranch        *b_SBS_mTPC_x;   //!
   TBranch        *b_SBS_mTPC_y;   //!
   TBranch        *b_SBS_mTPC_z;   //!
   TBranch        *b_SBS_mTPC_t;   //!
   TBranch        *b_SBS_mTPC_E;   //!
   TBranch        *b_SBS_mTPC_dt;   //!
   TBranch        *b_SBS_mTPC_L;   //!
   TBranch        *b_SBS_mTPC_vx;   //!
   TBranch        *b_SBS_mTPC_vy;   //!
   TBranch        *b_SBS_mTPC_vz;   //!
   TBranch        *b_SBS_mTPC_trid;   //!
   TBranch        *b_SBS_mTPC_mid;   //!
   TBranch        *b_SBS_mTPC_pid;   //!
   TBranch        *b_SBS_mTPC_p;   //!
   TBranch        *b_SBS_mTPC_px;   //!
   TBranch        *b_SBS_mTPC_py;   //!
   TBranch        *b_SBS_mTPC_pz;   //!
   TBranch        *b_SBS_mTPC_edep;   //!
   TBranch        *b_SBS_mTPC_ztravel;   //!
   TBranch        *b_SBS_mTPC_nstrips;   //!
   TBranch        *b_SBS_mTPC_part_npart;   //!
   TBranch        *b_SBS_mTPC_part_PID;   //!
   TBranch        *b_SBS_mTPC_part_MID;   //!
   TBranch        *b_SBS_mTPC_part_TID;   //!
   TBranch        *b_SBS_mTPC_part_nbounce;   //!
   TBranch        *b_SBS_mTPC_part_hitindex;   //!
   TBranch        *b_SBS_mTPC_part_vx;   //!
   TBranch        *b_SBS_mTPC_part_vy;   //!
   TBranch        *b_SBS_mTPC_part_vz;   //!
   TBranch        *b_SBS_mTPC_part_px;   //!
   TBranch        *b_SBS_mTPC_part_py;   //!
   TBranch        *b_SBS_mTPC_part_pz;   //!
   TBranch        *b_SBS_mTPCGEM_det_esum;   //!
   TBranch        *b_SBS_mTPCGEM_hit_nhits;   //!
   TBranch        *b_SBS_mTPCGEM_hit_row;   //!
   TBranch        *b_SBS_mTPCGEM_hit_col;   //!
   TBranch        *b_SBS_mTPCGEM_hit_cell;   //!
   TBranch        *b_SBS_mTPCGEM_hit_plane;   //!
   TBranch        *b_SBS_mTPCGEM_hit_wire;   //!
   TBranch        *b_SBS_mTPCGEM_hit_xcell;   //!
   TBranch        *b_SBS_mTPCGEM_hit_ycell;   //!
   TBranch        *b_SBS_mTPCGEM_hit_zcell;   //!
   TBranch        *b_SBS_mTPCGEM_hit_xcellg;   //!
   TBranch        *b_SBS_mTPCGEM_hit_ycellg;   //!
   TBranch        *b_SBS_mTPCGEM_hit_zcellg;   //!
   TBranch        *b_SBS_mTPCGEM_hit_xhit;   //!
   TBranch        *b_SBS_mTPCGEM_hit_yhit;   //!
   TBranch        *b_SBS_mTPCGEM_hit_zhit;   //!
   TBranch        *b_SBS_mTPCGEM_hit_xhitg;   //!
   TBranch        *b_SBS_mTPCGEM_hit_yhitg;   //!
   TBranch        *b_SBS_mTPCGEM_hit_zhitg;   //!
   TBranch        *b_SBS_mTPCGEM_hit_sumedep;   //!
   TBranch        *b_SBS_mTPCGEM_hit_tavg;   //!
   TBranch        *b_SBS_mTPCGEM_hit_trms;   //!
   TBranch        *b_SBS_mTPCGEM_hit_tmin;   //!
   TBranch        *b_SBS_mTPCGEM_hit_tmax;   //!
   TBranch        *b_SBS_mTPCHV_det_esum;   //!
   TBranch        *b_SBS_mTPCHV_hit_nhits;   //!
   TBranch        *b_SBS_mTPCHV_hit_row;   //!
   TBranch        *b_SBS_mTPCHV_hit_col;   //!
   TBranch        *b_SBS_mTPCHV_hit_cell;   //!
   TBranch        *b_SBS_mTPCHV_hit_plane;   //!
   TBranch        *b_SBS_mTPCHV_hit_wire;   //!
   TBranch        *b_SBS_mTPCHV_hit_xcell;   //!
   TBranch        *b_SBS_mTPCHV_hit_ycell;   //!
   TBranch        *b_SBS_mTPCHV_hit_zcell;   //!
   TBranch        *b_SBS_mTPCHV_hit_xcellg;   //!
   TBranch        *b_SBS_mTPCHV_hit_ycellg;   //!
   TBranch        *b_SBS_mTPCHV_hit_zcellg;   //!
   TBranch        *b_SBS_mTPCHV_hit_xhit;   //!
   TBranch        *b_SBS_mTPCHV_hit_yhit;   //!
   TBranch        *b_SBS_mTPCHV_hit_zhit;   //!
   TBranch        *b_SBS_mTPCHV_hit_xhitg;   //!
   TBranch        *b_SBS_mTPCHV_hit_yhitg;   //!
   TBranch        *b_SBS_mTPCHV_hit_zhitg;   //!
   TBranch        *b_SBS_mTPCHV_hit_sumedep;   //!
   TBranch        *b_SBS_mTPCHV_hit_tavg;   //!
   TBranch        *b_SBS_mTPCHV_hit_trms;   //!
   TBranch        *b_SBS_mTPCHV_hit_tmin;   //!
   TBranch        *b_SBS_mTPCHV_hit_tmax;   //!
   TBranch        *b_SBS_mTPCReadout_det_esum;   //!
   TBranch        *b_SBS_mTPCReadout_hit_nhits;   //!
   TBranch        *b_SBS_mTPCReadout_hit_row;   //!
   TBranch        *b_SBS_mTPCReadout_hit_col;   //!
   TBranch        *b_SBS_mTPCReadout_hit_cell;   //!
   TBranch        *b_SBS_mTPCReadout_hit_plane;   //!
   TBranch        *b_SBS_mTPCReadout_hit_wire;   //!
   TBranch        *b_SBS_mTPCReadout_hit_xcell;   //!
   TBranch        *b_SBS_mTPCReadout_hit_ycell;   //!
   TBranch        *b_SBS_mTPCReadout_hit_zcell;   //!
   TBranch        *b_SBS_mTPCReadout_hit_xcellg;   //!
   TBranch        *b_SBS_mTPCReadout_hit_ycellg;   //!
   TBranch        *b_SBS_mTPCReadout_hit_zcellg;   //!
   TBranch        *b_SBS_mTPCReadout_hit_xhit;   //!
   TBranch        *b_SBS_mTPCReadout_hit_yhit;   //!
   TBranch        *b_SBS_mTPCReadout_hit_zhit;   //!
   TBranch        *b_SBS_mTPCReadout_hit_xhitg;   //!
   TBranch        *b_SBS_mTPCReadout_hit_yhitg;   //!
   TBranch        *b_SBS_mTPCReadout_hit_zhitg;   //!
   TBranch        *b_SBS_mTPCReadout_hit_sumedep;   //!
   TBranch        *b_SBS_mTPCReadout_hit_tavg;   //!
   TBranch        *b_SBS_mTPCReadout_hit_trms;   //!
   TBranch        *b_SBS_mTPCReadout_hit_tmin;   //!
   TBranch        *b_SBS_mTPCReadout_hit_tmax;   //!

   g4sbstree_mtpc(TTree *tree=0);
   virtual ~g4sbstree_mtpc();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef g4sbstree_mtpc_cxx
g4sbstree_mtpc::g4sbstree_mtpc(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("g4sbs_out_protons_60-400MeV_30-70degrees_1kevents.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("g4sbs_out_protons_60-400MeV_30-70degrees_1kevents.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

g4sbstree_mtpc::~g4sbstree_mtpc()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t g4sbstree_mtpc::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t g4sbstree_mtpc::LoadTree(Long64_t entry)
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

void g4sbstree_mtpc::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Harm_LACScint_hit_row = 0;
   Harm_LACScint_hit_col = 0;
   Harm_LACScint_hit_cell = 0;
   Harm_LACScint_hit_plane = 0;
   Harm_LACScint_hit_wire = 0;
   Harm_LACScint_hit_xcell = 0;
   Harm_LACScint_hit_ycell = 0;
   Harm_LACScint_hit_zcell = 0;
   Harm_LACScint_hit_xcellg = 0;
   Harm_LACScint_hit_ycellg = 0;
   Harm_LACScint_hit_zcellg = 0;
   Harm_LACScint_hit_xhit = 0;
   Harm_LACScint_hit_yhit = 0;
   Harm_LACScint_hit_zhit = 0;
   Harm_LACScint_hit_xhitg = 0;
   Harm_LACScint_hit_yhitg = 0;
   Harm_LACScint_hit_zhitg = 0;
   Harm_LACScint_hit_sumedep = 0;
   Harm_LACScint_hit_tavg = 0;
   Harm_LACScint_hit_trms = 0;
   Harm_LACScint_hit_tmin = 0;
   Harm_LACScint_hit_tmax = 0;
   Harm_RICH_hit_PMT = 0;
   Harm_RICH_hit_row = 0;
   Harm_RICH_hit_col = 0;
   Harm_RICH_hit_xpmt = 0;
   Harm_RICH_hit_ypmt = 0;
   Harm_RICH_hit_zpmt = 0;
   Harm_RICH_hit_xgpmt = 0;
   Harm_RICH_hit_ygpmt = 0;
   Harm_RICH_hit_zgpmt = 0;
   Harm_RICH_hit_NumPhotoelectrons = 0;
   Harm_RICH_hit_Time_avg = 0;
   Harm_RICH_hit_Time_rms = 0;
   Harm_RICH_hit_Time_min = 0;
   Harm_RICH_hit_Time_max = 0;
   Harm_RICH_hit_mTrackNo = 0;
   Harm_RICH_hit_xhit = 0;
   Harm_RICH_hit_yhit = 0;
   Harm_RICH_hit_zhit = 0;
   Harm_RICH_hit_pxhit = 0;
   Harm_RICH_hit_pyhit = 0;
   Harm_RICH_hit_pzhit = 0;
   Harm_RICH_hit_pvx = 0;
   Harm_RICH_hit_pvy = 0;
   Harm_RICH_hit_pvz = 0;
   Harm_RICH_hit_ppx = 0;
   Harm_RICH_hit_ppy = 0;
   Harm_RICH_hit_ppz = 0;
   Harm_RICH_hit_volume_flag = 0;
   Harm_SBSGEM_hit_plane = 0;
   Harm_SBSGEM_hit_strip = 0;
   Harm_SBSGEM_hit_x = 0;
   Harm_SBSGEM_hit_y = 0;
   Harm_SBSGEM_hit_z = 0;
   Harm_SBSGEM_hit_polx = 0;
   Harm_SBSGEM_hit_poly = 0;
   Harm_SBSGEM_hit_polz = 0;
   Harm_SBSGEM_hit_t = 0;
   Harm_SBSGEM_hit_trms = 0;
   Harm_SBSGEM_hit_tmin = 0;
   Harm_SBSGEM_hit_tmax = 0;
   Harm_SBSGEM_hit_tx = 0;
   Harm_SBSGEM_hit_ty = 0;
   Harm_SBSGEM_hit_xin = 0;
   Harm_SBSGEM_hit_yin = 0;
   Harm_SBSGEM_hit_zin = 0;
   Harm_SBSGEM_hit_xout = 0;
   Harm_SBSGEM_hit_yout = 0;
   Harm_SBSGEM_hit_zout = 0;
   Harm_SBSGEM_hit_txp = 0;
   Harm_SBSGEM_hit_typ = 0;
   Harm_SBSGEM_hit_xg = 0;
   Harm_SBSGEM_hit_yg = 0;
   Harm_SBSGEM_hit_zg = 0;
   Harm_SBSGEM_hit_trid = 0;
   Harm_SBSGEM_hit_mid = 0;
   Harm_SBSGEM_hit_pid = 0;
   Harm_SBSGEM_hit_vx = 0;
   Harm_SBSGEM_hit_vy = 0;
   Harm_SBSGEM_hit_vz = 0;
   Harm_SBSGEM_hit_p = 0;
   Harm_SBSGEM_hit_edep = 0;
   Harm_SBSGEM_hit_beta = 0;
   Harm_SBSGEM_Track_TID = 0;
   Harm_SBSGEM_Track_PID = 0;
   Harm_SBSGEM_Track_MID = 0;
   Harm_SBSGEM_Track_NumHits = 0;
   Harm_SBSGEM_Track_NumPlanes = 0;
   Harm_SBSGEM_Track_NDF = 0;
   Harm_SBSGEM_Track_Chi2fit = 0;
   Harm_SBSGEM_Track_Chi2true = 0;
   Harm_SBSGEM_Track_X = 0;
   Harm_SBSGEM_Track_Y = 0;
   Harm_SBSGEM_Track_Xp = 0;
   Harm_SBSGEM_Track_Yp = 0;
   Harm_SBSGEM_Track_T = 0;
   Harm_SBSGEM_Track_P = 0;
   Harm_SBSGEM_Track_Sx = 0;
   Harm_SBSGEM_Track_Sy = 0;
   Harm_SBSGEM_Track_Sz = 0;
   Harm_SBSGEM_Track_Xfit = 0;
   Harm_SBSGEM_Track_Yfit = 0;
   Harm_SBSGEM_Track_Xpfit = 0;
   Harm_SBSGEM_Track_Ypfit = 0;
   SBS_mTPC_hit_cell = 0;
   SBS_mTPC_hit_xcell = 0;
   SBS_mTPC_hit_ycell = 0;
   SBS_mTPC_hit_zcell = 0;
   SBS_mTPC_hit_xcellg = 0;
   SBS_mTPC_hit_ycellg = 0;
   SBS_mTPC_hit_zcellg = 0;
   SBS_mTPC_hit_xhit = 0;
   SBS_mTPC_hit_yhit = 0;
   SBS_mTPC_hit_zhit = 0;
   SBS_mTPC_hit_xhitg = 0;
   SBS_mTPC_hit_yhitg = 0;
   SBS_mTPC_hit_zhitg = 0;
   SBS_mTPC_hit_sumedep = 0;
   SBS_mTPC_hit_tavg = 0;
   SBS_mTPC_hit_trms = 0;
   SBS_mTPC_hit_tmin = 0;
   SBS_mTPC_hit_tmax = 0;
   SBS_mTPC_ihit = 0;
   SBS_mTPC_x = 0;
   SBS_mTPC_y = 0;
   SBS_mTPC_z = 0;
   SBS_mTPC_t = 0;
   SBS_mTPC_E = 0;
   SBS_mTPC_dt = 0;
   SBS_mTPC_L = 0;
   SBS_mTPC_vx = 0;
   SBS_mTPC_vy = 0;
   SBS_mTPC_vz = 0;
   SBS_mTPC_trid = 0;
   SBS_mTPC_mid = 0;
   SBS_mTPC_pid = 0;
   SBS_mTPC_p = 0;
   SBS_mTPC_px = 0;
   SBS_mTPC_py = 0;
   SBS_mTPC_pz = 0;
   SBS_mTPC_edep = 0;
   SBS_mTPC_ztravel = 0;
   SBS_mTPC_nstrips = 0;
   SBS_mTPC_part_PID = 0;
   SBS_mTPC_part_MID = 0;
   SBS_mTPC_part_TID = 0;
   SBS_mTPC_part_nbounce = 0;
   SBS_mTPC_part_hitindex = 0;
   SBS_mTPC_part_vx = 0;
   SBS_mTPC_part_vy = 0;
   SBS_mTPC_part_vz = 0;
   SBS_mTPC_part_px = 0;
   SBS_mTPC_part_py = 0;
   SBS_mTPC_part_pz = 0;
   SBS_mTPCGEM_hit_row = 0;
   SBS_mTPCGEM_hit_col = 0;
   SBS_mTPCGEM_hit_cell = 0;
   SBS_mTPCGEM_hit_plane = 0;
   SBS_mTPCGEM_hit_wire = 0;
   SBS_mTPCGEM_hit_xcell = 0;
   SBS_mTPCGEM_hit_ycell = 0;
   SBS_mTPCGEM_hit_zcell = 0;
   SBS_mTPCGEM_hit_xcellg = 0;
   SBS_mTPCGEM_hit_ycellg = 0;
   SBS_mTPCGEM_hit_zcellg = 0;
   SBS_mTPCGEM_hit_xhit = 0;
   SBS_mTPCGEM_hit_yhit = 0;
   SBS_mTPCGEM_hit_zhit = 0;
   SBS_mTPCGEM_hit_xhitg = 0;
   SBS_mTPCGEM_hit_yhitg = 0;
   SBS_mTPCGEM_hit_zhitg = 0;
   SBS_mTPCGEM_hit_sumedep = 0;
   SBS_mTPCGEM_hit_tavg = 0;
   SBS_mTPCGEM_hit_trms = 0;
   SBS_mTPCGEM_hit_tmin = 0;
   SBS_mTPCGEM_hit_tmax = 0;
   SBS_mTPCHV_hit_row = 0;
   SBS_mTPCHV_hit_col = 0;
   SBS_mTPCHV_hit_cell = 0;
   SBS_mTPCHV_hit_plane = 0;
   SBS_mTPCHV_hit_wire = 0;
   SBS_mTPCHV_hit_xcell = 0;
   SBS_mTPCHV_hit_ycell = 0;
   SBS_mTPCHV_hit_zcell = 0;
   SBS_mTPCHV_hit_xcellg = 0;
   SBS_mTPCHV_hit_ycellg = 0;
   SBS_mTPCHV_hit_zcellg = 0;
   SBS_mTPCHV_hit_xhit = 0;
   SBS_mTPCHV_hit_yhit = 0;
   SBS_mTPCHV_hit_zhit = 0;
   SBS_mTPCHV_hit_xhitg = 0;
   SBS_mTPCHV_hit_yhitg = 0;
   SBS_mTPCHV_hit_zhitg = 0;
   SBS_mTPCHV_hit_sumedep = 0;
   SBS_mTPCHV_hit_tavg = 0;
   SBS_mTPCHV_hit_trms = 0;
   SBS_mTPCHV_hit_tmin = 0;
   SBS_mTPCHV_hit_tmax = 0;
   SBS_mTPCReadout_hit_row = 0;
   SBS_mTPCReadout_hit_col = 0;
   SBS_mTPCReadout_hit_cell = 0;
   SBS_mTPCReadout_hit_plane = 0;
   SBS_mTPCReadout_hit_wire = 0;
   SBS_mTPCReadout_hit_xcell = 0;
   SBS_mTPCReadout_hit_ycell = 0;
   SBS_mTPCReadout_hit_zcell = 0;
   SBS_mTPCReadout_hit_xcellg = 0;
   SBS_mTPCReadout_hit_ycellg = 0;
   SBS_mTPCReadout_hit_zcellg = 0;
   SBS_mTPCReadout_hit_xhit = 0;
   SBS_mTPCReadout_hit_yhit = 0;
   SBS_mTPCReadout_hit_zhit = 0;
   SBS_mTPCReadout_hit_xhitg = 0;
   SBS_mTPCReadout_hit_yhitg = 0;
   SBS_mTPCReadout_hit_zhitg = 0;
   SBS_mTPCReadout_hit_sumedep = 0;
   SBS_mTPCReadout_hit_tavg = 0;
   SBS_mTPCReadout_hit_trms = 0;
   SBS_mTPCReadout_hit_tmin = 0;
   SBS_mTPCReadout_hit_tmax = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ev", &ev_count, &b_ev);
   fChain->SetBranchAddress("Harm.LACScint.det.esum", &Harm_LACScint_det_esum, &b_Harm_LACScint_det_esum);
   fChain->SetBranchAddress("Harm.LACScint.hit.nhits", &Harm_LACScint_hit_nhits, &b_Harm_LACScint_hit_nhits);
   fChain->SetBranchAddress("Harm.LACScint.hit.row", &Harm_LACScint_hit_row, &b_Harm_LACScint_hit_row);
   fChain->SetBranchAddress("Harm.LACScint.hit.col", &Harm_LACScint_hit_col, &b_Harm_LACScint_hit_col);
   fChain->SetBranchAddress("Harm.LACScint.hit.cell", &Harm_LACScint_hit_cell, &b_Harm_LACScint_hit_cell);
   fChain->SetBranchAddress("Harm.LACScint.hit.plane", &Harm_LACScint_hit_plane, &b_Harm_LACScint_hit_plane);
   fChain->SetBranchAddress("Harm.LACScint.hit.wire", &Harm_LACScint_hit_wire, &b_Harm_LACScint_hit_wire);
   fChain->SetBranchAddress("Harm.LACScint.hit.xcell", &Harm_LACScint_hit_xcell, &b_Harm_LACScint_hit_xcell);
   fChain->SetBranchAddress("Harm.LACScint.hit.ycell", &Harm_LACScint_hit_ycell, &b_Harm_LACScint_hit_ycell);
   fChain->SetBranchAddress("Harm.LACScint.hit.zcell", &Harm_LACScint_hit_zcell, &b_Harm_LACScint_hit_zcell);
   fChain->SetBranchAddress("Harm.LACScint.hit.xcellg", &Harm_LACScint_hit_xcellg, &b_Harm_LACScint_hit_xcellg);
   fChain->SetBranchAddress("Harm.LACScint.hit.ycellg", &Harm_LACScint_hit_ycellg, &b_Harm_LACScint_hit_ycellg);
   fChain->SetBranchAddress("Harm.LACScint.hit.zcellg", &Harm_LACScint_hit_zcellg, &b_Harm_LACScint_hit_zcellg);
   fChain->SetBranchAddress("Harm.LACScint.hit.xhit", &Harm_LACScint_hit_xhit, &b_Harm_LACScint_hit_xhit);
   fChain->SetBranchAddress("Harm.LACScint.hit.yhit", &Harm_LACScint_hit_yhit, &b_Harm_LACScint_hit_yhit);
   fChain->SetBranchAddress("Harm.LACScint.hit.zhit", &Harm_LACScint_hit_zhit, &b_Harm_LACScint_hit_zhit);
   fChain->SetBranchAddress("Harm.LACScint.hit.xhitg", &Harm_LACScint_hit_xhitg, &b_Harm_LACScint_hit_xhitg);
   fChain->SetBranchAddress("Harm.LACScint.hit.yhitg", &Harm_LACScint_hit_yhitg, &b_Harm_LACScint_hit_yhitg);
   fChain->SetBranchAddress("Harm.LACScint.hit.zhitg", &Harm_LACScint_hit_zhitg, &b_Harm_LACScint_hit_zhitg);
   fChain->SetBranchAddress("Harm.LACScint.hit.sumedep", &Harm_LACScint_hit_sumedep, &b_Harm_LACScint_hit_sumedep);
   fChain->SetBranchAddress("Harm.LACScint.hit.tavg", &Harm_LACScint_hit_tavg, &b_Harm_LACScint_hit_tavg);
   fChain->SetBranchAddress("Harm.LACScint.hit.trms", &Harm_LACScint_hit_trms, &b_Harm_LACScint_hit_trms);
   fChain->SetBranchAddress("Harm.LACScint.hit.tmin", &Harm_LACScint_hit_tmin, &b_Harm_LACScint_hit_tmin);
   fChain->SetBranchAddress("Harm.LACScint.hit.tmax", &Harm_LACScint_hit_tmax, &b_Harm_LACScint_hit_tmax);
   fChain->SetBranchAddress("Harm.RICH.hit.nhits", &Harm_RICH_hit_nhits, &b_Harm_RICH_hit_nhits);
   fChain->SetBranchAddress("Harm.RICH.hit.PMT", &Harm_RICH_hit_PMT, &b_Harm_RICH_hit_PMT);
   fChain->SetBranchAddress("Harm.RICH.hit.row", &Harm_RICH_hit_row, &b_Harm_RICH_hit_row);
   fChain->SetBranchAddress("Harm.RICH.hit.col", &Harm_RICH_hit_col, &b_Harm_RICH_hit_col);
   fChain->SetBranchAddress("Harm.RICH.hit.xpmt", &Harm_RICH_hit_xpmt, &b_Harm_RICH_hit_xpmt);
   fChain->SetBranchAddress("Harm.RICH.hit.ypmt", &Harm_RICH_hit_ypmt, &b_Harm_RICH_hit_ypmt);
   fChain->SetBranchAddress("Harm.RICH.hit.zpmt", &Harm_RICH_hit_zpmt, &b_Harm_RICH_hit_zpmt);
   fChain->SetBranchAddress("Harm.RICH.hit.xgpmt", &Harm_RICH_hit_xgpmt, &b_Harm_RICH_hit_xgpmt);
   fChain->SetBranchAddress("Harm.RICH.hit.ygpmt", &Harm_RICH_hit_ygpmt, &b_Harm_RICH_hit_ygpmt);
   fChain->SetBranchAddress("Harm.RICH.hit.zgpmt", &Harm_RICH_hit_zgpmt, &b_Harm_RICH_hit_zgpmt);
   fChain->SetBranchAddress("Harm.RICH.hit.NumPhotoelectrons", &Harm_RICH_hit_NumPhotoelectrons, &b_Harm_RICH_hit_NumPhotoelectrons);
   fChain->SetBranchAddress("Harm.RICH.hit.Time_avg", &Harm_RICH_hit_Time_avg, &b_Harm_RICH_hit_Time_avg);
   fChain->SetBranchAddress("Harm.RICH.hit.Time_rms", &Harm_RICH_hit_Time_rms, &b_Harm_RICH_hit_Time_rms);
   fChain->SetBranchAddress("Harm.RICH.hit.Time_min", &Harm_RICH_hit_Time_min, &b_Harm_RICH_hit_Time_min);
   fChain->SetBranchAddress("Harm.RICH.hit.Time_max", &Harm_RICH_hit_Time_max, &b_Harm_RICH_hit_Time_max);
   fChain->SetBranchAddress("Harm.RICH.hit.mTrackNo", &Harm_RICH_hit_mTrackNo, &b_Harm_RICH_hit_mTrackNo);
   fChain->SetBranchAddress("Harm.RICH.hit.xhit", &Harm_RICH_hit_xhit, &b_Harm_RICH_hit_xhit);
   fChain->SetBranchAddress("Harm.RICH.hit.yhit", &Harm_RICH_hit_yhit, &b_Harm_RICH_hit_yhit);
   fChain->SetBranchAddress("Harm.RICH.hit.zhit", &Harm_RICH_hit_zhit, &b_Harm_RICH_hit_zhit);
   fChain->SetBranchAddress("Harm.RICH.hit.pxhit", &Harm_RICH_hit_pxhit, &b_Harm_RICH_hit_pxhit);
   fChain->SetBranchAddress("Harm.RICH.hit.pyhit", &Harm_RICH_hit_pyhit, &b_Harm_RICH_hit_pyhit);
   fChain->SetBranchAddress("Harm.RICH.hit.pzhit", &Harm_RICH_hit_pzhit, &b_Harm_RICH_hit_pzhit);
   fChain->SetBranchAddress("Harm.RICH.hit.pvx", &Harm_RICH_hit_pvx, &b_Harm_RICH_hit_pvx);
   fChain->SetBranchAddress("Harm.RICH.hit.pvy", &Harm_RICH_hit_pvy, &b_Harm_RICH_hit_pvy);
   fChain->SetBranchAddress("Harm.RICH.hit.pvz", &Harm_RICH_hit_pvz, &b_Harm_RICH_hit_pvz);
   fChain->SetBranchAddress("Harm.RICH.hit.ppx", &Harm_RICH_hit_ppx, &b_Harm_RICH_hit_ppx);
   fChain->SetBranchAddress("Harm.RICH.hit.ppy", &Harm_RICH_hit_ppy, &b_Harm_RICH_hit_ppy);
   fChain->SetBranchAddress("Harm.RICH.hit.ppz", &Harm_RICH_hit_ppz, &b_Harm_RICH_hit_ppz);
   fChain->SetBranchAddress("Harm.RICH.hit.volume_flag", &Harm_RICH_hit_volume_flag, &b_Harm_RICH_hit_volume_flag);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.nhits", &Harm_SBSGEM_hit_nhits, &b_Harm_SBSGEM_hit_nhits);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.plane", &Harm_SBSGEM_hit_plane, &b_Harm_SBSGEM_hit_plane);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.strip", &Harm_SBSGEM_hit_strip, &b_Harm_SBSGEM_hit_strip);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.x", &Harm_SBSGEM_hit_x, &b_Harm_SBSGEM_hit_x);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.y", &Harm_SBSGEM_hit_y, &b_Harm_SBSGEM_hit_y);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.z", &Harm_SBSGEM_hit_z, &b_Harm_SBSGEM_hit_z);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.polx", &Harm_SBSGEM_hit_polx, &b_Harm_SBSGEM_hit_polx);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.poly", &Harm_SBSGEM_hit_poly, &b_Harm_SBSGEM_hit_poly);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.polz", &Harm_SBSGEM_hit_polz, &b_Harm_SBSGEM_hit_polz);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.t", &Harm_SBSGEM_hit_t, &b_Harm_SBSGEM_hit_t);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.trms", &Harm_SBSGEM_hit_trms, &b_Harm_SBSGEM_hit_trms);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.tmin", &Harm_SBSGEM_hit_tmin, &b_Harm_SBSGEM_hit_tmin);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.tmax", &Harm_SBSGEM_hit_tmax, &b_Harm_SBSGEM_hit_tmax);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.tx", &Harm_SBSGEM_hit_tx, &b_Harm_SBSGEM_hit_tx);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.ty", &Harm_SBSGEM_hit_ty, &b_Harm_SBSGEM_hit_ty);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.xin", &Harm_SBSGEM_hit_xin, &b_Harm_SBSGEM_hit_xin);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.yin", &Harm_SBSGEM_hit_yin, &b_Harm_SBSGEM_hit_yin);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.zin", &Harm_SBSGEM_hit_zin, &b_Harm_SBSGEM_hit_zin);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.xout", &Harm_SBSGEM_hit_xout, &b_Harm_SBSGEM_hit_xout);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.yout", &Harm_SBSGEM_hit_yout, &b_Harm_SBSGEM_hit_yout);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.zout", &Harm_SBSGEM_hit_zout, &b_Harm_SBSGEM_hit_zout);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.txp", &Harm_SBSGEM_hit_txp, &b_Harm_SBSGEM_hit_txp);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.typ", &Harm_SBSGEM_hit_typ, &b_Harm_SBSGEM_hit_typ);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.xg", &Harm_SBSGEM_hit_xg, &b_Harm_SBSGEM_hit_xg);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.yg", &Harm_SBSGEM_hit_yg, &b_Harm_SBSGEM_hit_yg);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.zg", &Harm_SBSGEM_hit_zg, &b_Harm_SBSGEM_hit_zg);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.trid", &Harm_SBSGEM_hit_trid, &b_Harm_SBSGEM_hit_trid);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.mid", &Harm_SBSGEM_hit_mid, &b_Harm_SBSGEM_hit_mid);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.pid", &Harm_SBSGEM_hit_pid, &b_Harm_SBSGEM_hit_pid);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.vx", &Harm_SBSGEM_hit_vx, &b_Harm_SBSGEM_hit_vx);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.vy", &Harm_SBSGEM_hit_vy, &b_Harm_SBSGEM_hit_vy);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.vz", &Harm_SBSGEM_hit_vz, &b_Harm_SBSGEM_hit_vz);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.p", &Harm_SBSGEM_hit_p, &b_Harm_SBSGEM_hit_p);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.edep", &Harm_SBSGEM_hit_edep, &b_Harm_SBSGEM_hit_edep);
   fChain->SetBranchAddress("Harm.SBSGEM.hit.beta", &Harm_SBSGEM_hit_beta, &b_Harm_SBSGEM_hit_beta);
   fChain->SetBranchAddress("Harm.SBSGEM.Track.ntracks", &Harm_SBSGEM_Track_ntracks, &b_Harm_SBSGEM_Track_ntracks);
   fChain->SetBranchAddress("Harm.SBSGEM.Track.TID", &Harm_SBSGEM_Track_TID, &b_Harm_SBSGEM_Track_TID);
   fChain->SetBranchAddress("Harm.SBSGEM.Track.PID", &Harm_SBSGEM_Track_PID, &b_Harm_SBSGEM_Track_PID);
   fChain->SetBranchAddress("Harm.SBSGEM.Track.MID", &Harm_SBSGEM_Track_MID, &b_Harm_SBSGEM_Track_MID);
   fChain->SetBranchAddress("Harm.SBSGEM.Track.NumHits", &Harm_SBSGEM_Track_NumHits, &b_Harm_SBSGEM_Track_NumHits);
   fChain->SetBranchAddress("Harm.SBSGEM.Track.NumPlanes", &Harm_SBSGEM_Track_NumPlanes, &b_Harm_SBSGEM_Track_NumPlanes);
   fChain->SetBranchAddress("Harm.SBSGEM.Track.NDF", &Harm_SBSGEM_Track_NDF, &b_Harm_SBSGEM_Track_NDF);
   fChain->SetBranchAddress("Harm.SBSGEM.Track.Chi2fit", &Harm_SBSGEM_Track_Chi2fit, &b_Harm_SBSGEM_Track_Chi2fit);
   fChain->SetBranchAddress("Harm.SBSGEM.Track.Chi2true", &Harm_SBSGEM_Track_Chi2true, &b_Harm_SBSGEM_Track_Chi2true);
   fChain->SetBranchAddress("Harm.SBSGEM.Track.X", &Harm_SBSGEM_Track_X, &b_Harm_SBSGEM_Track_X);
   fChain->SetBranchAddress("Harm.SBSGEM.Track.Y", &Harm_SBSGEM_Track_Y, &b_Harm_SBSGEM_Track_Y);
   fChain->SetBranchAddress("Harm.SBSGEM.Track.Xp", &Harm_SBSGEM_Track_Xp, &b_Harm_SBSGEM_Track_Xp);
   fChain->SetBranchAddress("Harm.SBSGEM.Track.Yp", &Harm_SBSGEM_Track_Yp, &b_Harm_SBSGEM_Track_Yp);
   fChain->SetBranchAddress("Harm.SBSGEM.Track.T", &Harm_SBSGEM_Track_T, &b_Harm_SBSGEM_Track_T);
   fChain->SetBranchAddress("Harm.SBSGEM.Track.P", &Harm_SBSGEM_Track_P, &b_Harm_SBSGEM_Track_P);
   fChain->SetBranchAddress("Harm.SBSGEM.Track.Sx", &Harm_SBSGEM_Track_Sx, &b_Harm_SBSGEM_Track_Sx);
   fChain->SetBranchAddress("Harm.SBSGEM.Track.Sy", &Harm_SBSGEM_Track_Sy, &b_Harm_SBSGEM_Track_Sy);
   fChain->SetBranchAddress("Harm.SBSGEM.Track.Sz", &Harm_SBSGEM_Track_Sz, &b_Harm_SBSGEM_Track_Sz);
   fChain->SetBranchAddress("Harm.SBSGEM.Track.Xfit", &Harm_SBSGEM_Track_Xfit, &b_Harm_SBSGEM_Track_Xfit);
   fChain->SetBranchAddress("Harm.SBSGEM.Track.Yfit", &Harm_SBSGEM_Track_Yfit, &b_Harm_SBSGEM_Track_Yfit);
   fChain->SetBranchAddress("Harm.SBSGEM.Track.Xpfit", &Harm_SBSGEM_Track_Xpfit, &b_Harm_SBSGEM_Track_Xpfit);
   fChain->SetBranchAddress("Harm.SBSGEM.Track.Ypfit", &Harm_SBSGEM_Track_Ypfit, &b_Harm_SBSGEM_Track_Ypfit);
   fChain->SetBranchAddress("SBS.mTPC.det.esum", &SBS_mTPC_det_esum, &b_SBS_mTPC_det_esum);
   fChain->SetBranchAddress("SBS.mTPC.hit.nhits", &SBS_mTPC_hit_nhits, &b_SBS_mTPC_hit_nhits);
   fChain->SetBranchAddress("SBS.mTPC.hit.cell", &SBS_mTPC_hit_cell, &b_SBS_mTPC_hit_cell);
   fChain->SetBranchAddress("SBS.mTPC.hit.xcell", &SBS_mTPC_hit_xcell, &b_SBS_mTPC_hit_xcell);
   fChain->SetBranchAddress("SBS.mTPC.hit.ycell", &SBS_mTPC_hit_ycell, &b_SBS_mTPC_hit_ycell);
   fChain->SetBranchAddress("SBS.mTPC.hit.zcell", &SBS_mTPC_hit_zcell, &b_SBS_mTPC_hit_zcell);
   fChain->SetBranchAddress("SBS.mTPC.hit.xcellg", &SBS_mTPC_hit_xcellg, &b_SBS_mTPC_hit_xcellg);
   fChain->SetBranchAddress("SBS.mTPC.hit.ycellg", &SBS_mTPC_hit_ycellg, &b_SBS_mTPC_hit_ycellg);
   fChain->SetBranchAddress("SBS.mTPC.hit.zcellg", &SBS_mTPC_hit_zcellg, &b_SBS_mTPC_hit_zcellg);
   fChain->SetBranchAddress("SBS.mTPC.hit.xhit", &SBS_mTPC_hit_xhit, &b_SBS_mTPC_hit_xhit);
   fChain->SetBranchAddress("SBS.mTPC.hit.yhit", &SBS_mTPC_hit_yhit, &b_SBS_mTPC_hit_yhit);
   fChain->SetBranchAddress("SBS.mTPC.hit.zhit", &SBS_mTPC_hit_zhit, &b_SBS_mTPC_hit_zhit);
   fChain->SetBranchAddress("SBS.mTPC.hit.xhitg", &SBS_mTPC_hit_xhitg, &b_SBS_mTPC_hit_xhitg);
   fChain->SetBranchAddress("SBS.mTPC.hit.yhitg", &SBS_mTPC_hit_yhitg, &b_SBS_mTPC_hit_yhitg);
   fChain->SetBranchAddress("SBS.mTPC.hit.zhitg", &SBS_mTPC_hit_zhitg, &b_SBS_mTPC_hit_zhitg);
   fChain->SetBranchAddress("SBS.mTPC.hit.sumedep", &SBS_mTPC_hit_sumedep, &b_SBS_mTPC_hit_sumedep);
   fChain->SetBranchAddress("SBS.mTPC.hit.tavg", &SBS_mTPC_hit_tavg, &b_SBS_mTPC_hit_tavg);
   fChain->SetBranchAddress("SBS.mTPC.hit.trms", &SBS_mTPC_hit_trms, &b_SBS_mTPC_hit_trms);
   fChain->SetBranchAddress("SBS.mTPC.hit.tmin", &SBS_mTPC_hit_tmin, &b_SBS_mTPC_hit_tmin);
   fChain->SetBranchAddress("SBS.mTPC.hit.tmax", &SBS_mTPC_hit_tmax, &b_SBS_mTPC_hit_tmax);
   fChain->SetBranchAddress("SBS.mTPC.npart_mTPC", &SBS_mTPC_npart_mTPC, &b_SBS_mTPC_npart_mTPC);
   fChain->SetBranchAddress("SBS.mTPC.ihit", &SBS_mTPC_ihit, &b_SBS_mTPC_ihit);
   fChain->SetBranchAddress("SBS.mTPC.x", &SBS_mTPC_x, &b_SBS_mTPC_x);
   fChain->SetBranchAddress("SBS.mTPC.y", &SBS_mTPC_y, &b_SBS_mTPC_y);
   fChain->SetBranchAddress("SBS.mTPC.z", &SBS_mTPC_z, &b_SBS_mTPC_z);
   fChain->SetBranchAddress("SBS.mTPC.t", &SBS_mTPC_t, &b_SBS_mTPC_t);
   fChain->SetBranchAddress("SBS.mTPC.E", &SBS_mTPC_E, &b_SBS_mTPC_E);
   fChain->SetBranchAddress("SBS.mTPC.dt", &SBS_mTPC_dt, &b_SBS_mTPC_dt);
   fChain->SetBranchAddress("SBS.mTPC.L", &SBS_mTPC_L, &b_SBS_mTPC_L);
   fChain->SetBranchAddress("SBS.mTPC.vx", &SBS_mTPC_vx, &b_SBS_mTPC_vx);
   fChain->SetBranchAddress("SBS.mTPC.vy", &SBS_mTPC_vy, &b_SBS_mTPC_vy);
   fChain->SetBranchAddress("SBS.mTPC.vz", &SBS_mTPC_vz, &b_SBS_mTPC_vz);
   fChain->SetBranchAddress("SBS.mTPC.trid", &SBS_mTPC_trid, &b_SBS_mTPC_trid);
   fChain->SetBranchAddress("SBS.mTPC.mid", &SBS_mTPC_mid, &b_SBS_mTPC_mid);
   fChain->SetBranchAddress("SBS.mTPC.pid", &SBS_mTPC_pid, &b_SBS_mTPC_pid);
   fChain->SetBranchAddress("SBS.mTPC.p", &SBS_mTPC_p, &b_SBS_mTPC_p);
   fChain->SetBranchAddress("SBS.mTPC.px", &SBS_mTPC_px, &b_SBS_mTPC_px);
   fChain->SetBranchAddress("SBS.mTPC.py", &SBS_mTPC_py, &b_SBS_mTPC_py);
   fChain->SetBranchAddress("SBS.mTPC.pz", &SBS_mTPC_pz, &b_SBS_mTPC_pz);
   fChain->SetBranchAddress("SBS.mTPC.edep", &SBS_mTPC_edep, &b_SBS_mTPC_edep);
   fChain->SetBranchAddress("SBS.mTPC.ztravel", &SBS_mTPC_ztravel, &b_SBS_mTPC_ztravel);
   fChain->SetBranchAddress("SBS.mTPC.nstrips", &SBS_mTPC_nstrips, &b_SBS_mTPC_nstrips);
   fChain->SetBranchAddress("SBS.mTPC.part.npart", &SBS_mTPC_part_npart, &b_SBS_mTPC_part_npart);
   fChain->SetBranchAddress("SBS.mTPC.part.PID", &SBS_mTPC_part_PID, &b_SBS_mTPC_part_PID);
   fChain->SetBranchAddress("SBS.mTPC.part.MID", &SBS_mTPC_part_MID, &b_SBS_mTPC_part_MID);
   fChain->SetBranchAddress("SBS.mTPC.part.TID", &SBS_mTPC_part_TID, &b_SBS_mTPC_part_TID);
   fChain->SetBranchAddress("SBS.mTPC.part.nbounce", &SBS_mTPC_part_nbounce, &b_SBS_mTPC_part_nbounce);
   fChain->SetBranchAddress("SBS.mTPC.part.hitindex", &SBS_mTPC_part_hitindex, &b_SBS_mTPC_part_hitindex);
   fChain->SetBranchAddress("SBS.mTPC.part.vx", &SBS_mTPC_part_vx, &b_SBS_mTPC_part_vx);
   fChain->SetBranchAddress("SBS.mTPC.part.vy", &SBS_mTPC_part_vy, &b_SBS_mTPC_part_vy);
   fChain->SetBranchAddress("SBS.mTPC.part.vz", &SBS_mTPC_part_vz, &b_SBS_mTPC_part_vz);
   fChain->SetBranchAddress("SBS.mTPC.part.px", &SBS_mTPC_part_px, &b_SBS_mTPC_part_px);
   fChain->SetBranchAddress("SBS.mTPC.part.py", &SBS_mTPC_part_py, &b_SBS_mTPC_part_py);
   fChain->SetBranchAddress("SBS.mTPC.part.pz", &SBS_mTPC_part_pz, &b_SBS_mTPC_part_pz);
   fChain->SetBranchAddress("SBS.mTPCGEM.det.esum", &SBS_mTPCGEM_det_esum, &b_SBS_mTPCGEM_det_esum);
   fChain->SetBranchAddress("SBS.mTPCGEM.hit.nhits", &SBS_mTPCGEM_hit_nhits, &b_SBS_mTPCGEM_hit_nhits);
   fChain->SetBranchAddress("SBS.mTPCGEM.hit.row", &SBS_mTPCGEM_hit_row, &b_SBS_mTPCGEM_hit_row);
   fChain->SetBranchAddress("SBS.mTPCGEM.hit.col", &SBS_mTPCGEM_hit_col, &b_SBS_mTPCGEM_hit_col);
   fChain->SetBranchAddress("SBS.mTPCGEM.hit.cell", &SBS_mTPCGEM_hit_cell, &b_SBS_mTPCGEM_hit_cell);
   fChain->SetBranchAddress("SBS.mTPCGEM.hit.plane", &SBS_mTPCGEM_hit_plane, &b_SBS_mTPCGEM_hit_plane);
   fChain->SetBranchAddress("SBS.mTPCGEM.hit.wire", &SBS_mTPCGEM_hit_wire, &b_SBS_mTPCGEM_hit_wire);
   fChain->SetBranchAddress("SBS.mTPCGEM.hit.xcell", &SBS_mTPCGEM_hit_xcell, &b_SBS_mTPCGEM_hit_xcell);
   fChain->SetBranchAddress("SBS.mTPCGEM.hit.ycell", &SBS_mTPCGEM_hit_ycell, &b_SBS_mTPCGEM_hit_ycell);
   fChain->SetBranchAddress("SBS.mTPCGEM.hit.zcell", &SBS_mTPCGEM_hit_zcell, &b_SBS_mTPCGEM_hit_zcell);
   fChain->SetBranchAddress("SBS.mTPCGEM.hit.xcellg", &SBS_mTPCGEM_hit_xcellg, &b_SBS_mTPCGEM_hit_xcellg);
   fChain->SetBranchAddress("SBS.mTPCGEM.hit.ycellg", &SBS_mTPCGEM_hit_ycellg, &b_SBS_mTPCGEM_hit_ycellg);
   fChain->SetBranchAddress("SBS.mTPCGEM.hit.zcellg", &SBS_mTPCGEM_hit_zcellg, &b_SBS_mTPCGEM_hit_zcellg);
   fChain->SetBranchAddress("SBS.mTPCGEM.hit.xhit", &SBS_mTPCGEM_hit_xhit, &b_SBS_mTPCGEM_hit_xhit);
   fChain->SetBranchAddress("SBS.mTPCGEM.hit.yhit", &SBS_mTPCGEM_hit_yhit, &b_SBS_mTPCGEM_hit_yhit);
   fChain->SetBranchAddress("SBS.mTPCGEM.hit.zhit", &SBS_mTPCGEM_hit_zhit, &b_SBS_mTPCGEM_hit_zhit);
   fChain->SetBranchAddress("SBS.mTPCGEM.hit.xhitg", &SBS_mTPCGEM_hit_xhitg, &b_SBS_mTPCGEM_hit_xhitg);
   fChain->SetBranchAddress("SBS.mTPCGEM.hit.yhitg", &SBS_mTPCGEM_hit_yhitg, &b_SBS_mTPCGEM_hit_yhitg);
   fChain->SetBranchAddress("SBS.mTPCGEM.hit.zhitg", &SBS_mTPCGEM_hit_zhitg, &b_SBS_mTPCGEM_hit_zhitg);
   fChain->SetBranchAddress("SBS.mTPCGEM.hit.sumedep", &SBS_mTPCGEM_hit_sumedep, &b_SBS_mTPCGEM_hit_sumedep);
   fChain->SetBranchAddress("SBS.mTPCGEM.hit.tavg", &SBS_mTPCGEM_hit_tavg, &b_SBS_mTPCGEM_hit_tavg);
   fChain->SetBranchAddress("SBS.mTPCGEM.hit.trms", &SBS_mTPCGEM_hit_trms, &b_SBS_mTPCGEM_hit_trms);
   fChain->SetBranchAddress("SBS.mTPCGEM.hit.tmin", &SBS_mTPCGEM_hit_tmin, &b_SBS_mTPCGEM_hit_tmin);
   fChain->SetBranchAddress("SBS.mTPCGEM.hit.tmax", &SBS_mTPCGEM_hit_tmax, &b_SBS_mTPCGEM_hit_tmax);
   fChain->SetBranchAddress("SBS.mTPCHV.det.esum", &SBS_mTPCHV_det_esum, &b_SBS_mTPCHV_det_esum);
   fChain->SetBranchAddress("SBS.mTPCHV.hit.nhits", &SBS_mTPCHV_hit_nhits, &b_SBS_mTPCHV_hit_nhits);
   fChain->SetBranchAddress("SBS.mTPCHV.hit.row", &SBS_mTPCHV_hit_row, &b_SBS_mTPCHV_hit_row);
   fChain->SetBranchAddress("SBS.mTPCHV.hit.col", &SBS_mTPCHV_hit_col, &b_SBS_mTPCHV_hit_col);
   fChain->SetBranchAddress("SBS.mTPCHV.hit.cell", &SBS_mTPCHV_hit_cell, &b_SBS_mTPCHV_hit_cell);
   fChain->SetBranchAddress("SBS.mTPCHV.hit.plane", &SBS_mTPCHV_hit_plane, &b_SBS_mTPCHV_hit_plane);
   fChain->SetBranchAddress("SBS.mTPCHV.hit.wire", &SBS_mTPCHV_hit_wire, &b_SBS_mTPCHV_hit_wire);
   fChain->SetBranchAddress("SBS.mTPCHV.hit.xcell", &SBS_mTPCHV_hit_xcell, &b_SBS_mTPCHV_hit_xcell);
   fChain->SetBranchAddress("SBS.mTPCHV.hit.ycell", &SBS_mTPCHV_hit_ycell, &b_SBS_mTPCHV_hit_ycell);
   fChain->SetBranchAddress("SBS.mTPCHV.hit.zcell", &SBS_mTPCHV_hit_zcell, &b_SBS_mTPCHV_hit_zcell);
   fChain->SetBranchAddress("SBS.mTPCHV.hit.xcellg", &SBS_mTPCHV_hit_xcellg, &b_SBS_mTPCHV_hit_xcellg);
   fChain->SetBranchAddress("SBS.mTPCHV.hit.ycellg", &SBS_mTPCHV_hit_ycellg, &b_SBS_mTPCHV_hit_ycellg);
   fChain->SetBranchAddress("SBS.mTPCHV.hit.zcellg", &SBS_mTPCHV_hit_zcellg, &b_SBS_mTPCHV_hit_zcellg);
   fChain->SetBranchAddress("SBS.mTPCHV.hit.xhit", &SBS_mTPCHV_hit_xhit, &b_SBS_mTPCHV_hit_xhit);
   fChain->SetBranchAddress("SBS.mTPCHV.hit.yhit", &SBS_mTPCHV_hit_yhit, &b_SBS_mTPCHV_hit_yhit);
   fChain->SetBranchAddress("SBS.mTPCHV.hit.zhit", &SBS_mTPCHV_hit_zhit, &b_SBS_mTPCHV_hit_zhit);
   fChain->SetBranchAddress("SBS.mTPCHV.hit.xhitg", &SBS_mTPCHV_hit_xhitg, &b_SBS_mTPCHV_hit_xhitg);
   fChain->SetBranchAddress("SBS.mTPCHV.hit.yhitg", &SBS_mTPCHV_hit_yhitg, &b_SBS_mTPCHV_hit_yhitg);
   fChain->SetBranchAddress("SBS.mTPCHV.hit.zhitg", &SBS_mTPCHV_hit_zhitg, &b_SBS_mTPCHV_hit_zhitg);
   fChain->SetBranchAddress("SBS.mTPCHV.hit.sumedep", &SBS_mTPCHV_hit_sumedep, &b_SBS_mTPCHV_hit_sumedep);
   fChain->SetBranchAddress("SBS.mTPCHV.hit.tavg", &SBS_mTPCHV_hit_tavg, &b_SBS_mTPCHV_hit_tavg);
   fChain->SetBranchAddress("SBS.mTPCHV.hit.trms", &SBS_mTPCHV_hit_trms, &b_SBS_mTPCHV_hit_trms);
   fChain->SetBranchAddress("SBS.mTPCHV.hit.tmin", &SBS_mTPCHV_hit_tmin, &b_SBS_mTPCHV_hit_tmin);
   fChain->SetBranchAddress("SBS.mTPCHV.hit.tmax", &SBS_mTPCHV_hit_tmax, &b_SBS_mTPCHV_hit_tmax);
   fChain->SetBranchAddress("SBS.mTPCReadout.det.esum", &SBS_mTPCReadout_det_esum, &b_SBS_mTPCReadout_det_esum);
   fChain->SetBranchAddress("SBS.mTPCReadout.hit.nhits", &SBS_mTPCReadout_hit_nhits, &b_SBS_mTPCReadout_hit_nhits);
   fChain->SetBranchAddress("SBS.mTPCReadout.hit.row", &SBS_mTPCReadout_hit_row, &b_SBS_mTPCReadout_hit_row);
   fChain->SetBranchAddress("SBS.mTPCReadout.hit.col", &SBS_mTPCReadout_hit_col, &b_SBS_mTPCReadout_hit_col);
   fChain->SetBranchAddress("SBS.mTPCReadout.hit.cell", &SBS_mTPCReadout_hit_cell, &b_SBS_mTPCReadout_hit_cell);
   fChain->SetBranchAddress("SBS.mTPCReadout.hit.plane", &SBS_mTPCReadout_hit_plane, &b_SBS_mTPCReadout_hit_plane);
   fChain->SetBranchAddress("SBS.mTPCReadout.hit.wire", &SBS_mTPCReadout_hit_wire, &b_SBS_mTPCReadout_hit_wire);
   fChain->SetBranchAddress("SBS.mTPCReadout.hit.xcell", &SBS_mTPCReadout_hit_xcell, &b_SBS_mTPCReadout_hit_xcell);
   fChain->SetBranchAddress("SBS.mTPCReadout.hit.ycell", &SBS_mTPCReadout_hit_ycell, &b_SBS_mTPCReadout_hit_ycell);
   fChain->SetBranchAddress("SBS.mTPCReadout.hit.zcell", &SBS_mTPCReadout_hit_zcell, &b_SBS_mTPCReadout_hit_zcell);
   fChain->SetBranchAddress("SBS.mTPCReadout.hit.xcellg", &SBS_mTPCReadout_hit_xcellg, &b_SBS_mTPCReadout_hit_xcellg);
   fChain->SetBranchAddress("SBS.mTPCReadout.hit.ycellg", &SBS_mTPCReadout_hit_ycellg, &b_SBS_mTPCReadout_hit_ycellg);
   fChain->SetBranchAddress("SBS.mTPCReadout.hit.zcellg", &SBS_mTPCReadout_hit_zcellg, &b_SBS_mTPCReadout_hit_zcellg);
   fChain->SetBranchAddress("SBS.mTPCReadout.hit.xhit", &SBS_mTPCReadout_hit_xhit, &b_SBS_mTPCReadout_hit_xhit);
   fChain->SetBranchAddress("SBS.mTPCReadout.hit.yhit", &SBS_mTPCReadout_hit_yhit, &b_SBS_mTPCReadout_hit_yhit);
   fChain->SetBranchAddress("SBS.mTPCReadout.hit.zhit", &SBS_mTPCReadout_hit_zhit, &b_SBS_mTPCReadout_hit_zhit);
   fChain->SetBranchAddress("SBS.mTPCReadout.hit.xhitg", &SBS_mTPCReadout_hit_xhitg, &b_SBS_mTPCReadout_hit_xhitg);
   fChain->SetBranchAddress("SBS.mTPCReadout.hit.yhitg", &SBS_mTPCReadout_hit_yhitg, &b_SBS_mTPCReadout_hit_yhitg);
   fChain->SetBranchAddress("SBS.mTPCReadout.hit.zhitg", &SBS_mTPCReadout_hit_zhitg, &b_SBS_mTPCReadout_hit_zhitg);
   fChain->SetBranchAddress("SBS.mTPCReadout.hit.sumedep", &SBS_mTPCReadout_hit_sumedep, &b_SBS_mTPCReadout_hit_sumedep);
   fChain->SetBranchAddress("SBS.mTPCReadout.hit.tavg", &SBS_mTPCReadout_hit_tavg, &b_SBS_mTPCReadout_hit_tavg);
   fChain->SetBranchAddress("SBS.mTPCReadout.hit.trms", &SBS_mTPCReadout_hit_trms, &b_SBS_mTPCReadout_hit_trms);
   fChain->SetBranchAddress("SBS.mTPCReadout.hit.tmin", &SBS_mTPCReadout_hit_tmin, &b_SBS_mTPCReadout_hit_tmin);
   fChain->SetBranchAddress("SBS.mTPCReadout.hit.tmax", &SBS_mTPCReadout_hit_tmax, &b_SBS_mTPCReadout_hit_tmax);
   Notify();
}

Bool_t g4sbstree_mtpc::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void g4sbstree_mtpc::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t g4sbstree_mtpc::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef g4sbstree_mtpc_cxx
