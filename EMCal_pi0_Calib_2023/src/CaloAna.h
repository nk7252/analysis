#ifndef CALOANA_H__
#define CALOANA_H__

#include <fun4all/SubsysReco.h>
#include <random>
#include <vector>
// Forward declarations
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TNtuple;
class TTree;
class TH3F;
class TH2F;
class TH1F;
class TH1;
class TH3;
class TF1;
class TProfile2D;
class TProfile;
class TLorentzVector;
class TRandom3;

class CaloAna : public SubsysReco
{
 private:
  std::mt19937 rng;  // Mersenne Twister random number generator
 public:
  //! constructor
  CaloAna(const std::string& name = "CaloAna", const std::string& fname = "MyNtuple.root");

  //! destructor
  virtual ~CaloAna();

  //! full initialization
  int Init(PHCompositeNode*) override;

  //! event processing method
  int process_event(PHCompositeNode*) override;

  //! end of run method
  int End(PHCompositeNode*) override;

  int process_g4hits(PHCompositeNode*);
  int process_g4cells(PHCompositeNode*);
  int process_towers(PHCompositeNode*);
  int process_clusters(PHCompositeNode*);

  void Detector(const std::string& name) { detector = name; }
  void set_timing_cut_width(const int& t) { _range = t; }
  void set_vertex_cut(const float& v) { _vz = v; }
  void apply_vertex_cut(bool Vtx_cut) { m_vtxCut = Vtx_cut; }

  float getWeight(int ieta, float pt);
  float getSPMCpTspectrum(float pt);

  TF1* fitHistogram(TH1* h);
  void fitEtaSlices(const std::string& infile, const std::string& fitOutFile, const std::string& cdbFile);

  double generateRandomNumber();
  TLorentzVector SmearPhoton4vector(TLorentzVector sourcephoton, double smearfactor);

 protected:
  std::string detector;
  std::string outfilename;
  int Getpeaktime(TH1* h);
  Fun4AllHistoManager* hm = nullptr;
  TFile* outfile = nullptr;
  TFile* frw= nullptr;
  TH1F* h_cutCounter = nullptr;

  TH2F* h_emcal_mbd_correlation = nullptr;
  TH2F* h_ohcal_mbd_correlation = nullptr;
  TH2F* h_ihcal_mbd_correlation = nullptr;
  TH2F* h_emcal_hcal_correlation = nullptr;
  TH2F* h_emcal_zdc_correlation = nullptr;
  TH2F* h_clusE_nTow{nullptr};

  TH1F* h_InvMass = nullptr;
  TH1F* h_InvMass_weighted = nullptr;
  TH1F* h_InvMass_w = nullptr;
  TH1F* h_InvMassMix = nullptr;
  TH2F* h_InvMass_2d = nullptr;
  TH2F* h_truthmatched_mass1_2d = nullptr;
  TH2F* h_truthmatched_mass2_2d = nullptr;
  TH2F* h_truthmatched_mass3_2d = nullptr;
  TH1F* h_InvMass_smear = nullptr;
  TH2F* h_InvMass_smear_2d = nullptr;
  TH1F* h_InvMass_smear_weighted = nullptr;
  TH2F* h_InvMass_smear_weighted_2d = nullptr;

  // should add truth histos too?
  // h_etaphidist_InvMass_over200M->Draw("LEGO")
  // h_etaphidist_InvMass_under200M->Draw("LEGO")
  // gPad->SetLogz()

  TH2F* h_cemc_etaphi = nullptr;
  TH2F* h_hcalin_etaphi = nullptr;
  TH2F* h_hcalout_etaphi = nullptr;
  TH2F* h_cemc_etaphi_wQA = nullptr;
  TH2F* h_hcalin_etaphi_wQA = nullptr;
  TH2F* h_hcalout_etaphi_wQA = nullptr;
  TH1* h_delR_recTrth = nullptr;
  TH1* h_delR_pionrecTrth = nullptr;
  TH1* h_totalzdc_e;

  TProfile2D* h_cemc_etaphi_time = nullptr;
  TProfile2D* h_hcalin_etaphi_time = nullptr;
  TProfile2D* h_hcalout_etaphi_time = nullptr;

  TProfile2D* h_cemc_etaphi_badChi2 = nullptr;
  TProfile2D* h_hcalin_etaphi_badChi2 = nullptr;
  TProfile2D* h_hcalout_etaphi_badChi2 = nullptr;

  TH1* hzdctime;
  TH1* hmbdtime;
  TH1* hemcaltime;
  TH1* hihcaltime;
  TH1* hohcaltime;

  TH1* hzdctime_cut;
  TH1* hmbdtime_cut;
  TH1* hemcaltime_cut;
  TH1* hihcaltime_cut;
  TH1* hohcaltime_cut;
  TH1* h_tower_e;

  TH1* hvtx_z_raw;
  TH1* hvtx_z_cut;

  TH1* hzdcSouthraw;
  TH1* hzdcNorthraw;
  TH1* hzdcSouthcalib;
  TH1* hzdcNorthcalib;

  TH1F* h_clusE;
  TH2F* h_etaphi_clus;

  TNtuple* g4hitntuple = nullptr;
  TNtuple* g4cellntuple = nullptr;
  TTree* towerntuple = nullptr;
  TNtuple* clusterntuple = nullptr;
  std::vector<float> m_energy;
  std::vector<int> m_etabin;
  std::vector<int> m_phibin;
  std::vector<int> m_time;

  std::vector<float> m_hcalin_energy;
  std::vector<int> m_hcalin_etabin;
  std::vector<int> m_hcalin_phibin;
  std::vector<int> m_hcalin_time;

  std::vector<float> m_hcalout_energy;
  std::vector<int> m_hcalout_etabin;
  std::vector<int> m_hcalout_phibin;
  std::vector<int> m_hcalout_time;

  std::vector<float> m_zdc_energy;
  std::vector<int> m_zdc_index;
  std::vector<int> m_zdc_side;

  std::vector<float> m_bbc_energy;
  std::vector<int> m_bbc_type;
  std::vector<int> m_bbc_side;

  int _eventcounter;
  int funkyCaloStuffcounter;
  int VertexMapFailcounter;

  // bool pileup = false;
  std::string clustcontainerstring;
  std::string calotowerinfostring;

  // cuts
  // float maxDr = 1.1;
  // float maxAlpha = 0.6;
  // float clus_chisq_cut = 4;
  // float nClus_ptCut = 0.5;
  // int max_nClusCount = 37;
  // float ptMaxCut = 7;  // 7 in data? ** keep this in mind. 3 may make more sense, but 7 is
  // float pt1ClusCut = 1.3;  // centrality dependence cuts 2.2 for both // 1.3
  // float pt2ClusCut = 0.7;  // // 0.7
  // float pi0ptcut = 1.22 * (pt1ClusCut + pt2ClusCut);

  int _range = 1;
  float _vz = 0.;
  bool m_vtxCut = false;
  bool dynMaskClus = false;
  bool getVtx = true;
  bool cutson = true;
  bool clust_waveform = true;  // if on will use  WAVEFORM_CEMC for towers, superseeds TOWERINFO_CALIB_CEMC
  bool recluster = false;      // if on will use CLUSTERINFO_CEMC2 for clusterContainer
  bool poscor = false;         // if on will use position correction for clusterContainer, superseeds CLUSTER_CEMC
  bool pp_rawcluster = true;   // if on will use GetEVec, superseeds GetECoreVec. ECore is for AuAu
  bool debug = false;
  bool Cluster_Debug = true; // cluster etaphi histograms
  bool Cluster_Debug2 = false; // true if you want to add the 2nd cluster for every pair, rather than just individual clusters
  //so if this is true, you will fill every time a cluster shows up past that cut, rather than just the first time
  bool etabyeta = false; // mass and pt eta histograms, 96 for each histogram(4)
  //cuts
  bool eTCutbool = false;       // replace pt1 and pt2 cuts with et cuts
  bool etaCutbool = false;      // do an eta cut or not
  bool clusterprobcut = false;  // if true use cluster prob cut, if false use chisq cut
  bool zvtxcut_bool = false;
  // SPMC
  bool SPMC_bool = true;
  float badcalibsmear;
  bool Pythia_weight = false;
  bool eta_weight = true;
  // NEVER USE LEADING ZEROS IN DECIMALS. IT WILL BE INTERPRETED AS OCTAL
  bool additionalsmearing = true;  // should be on for spmc in all cases. if you want 0 use 0 for smearint
  int badcalibsmearint = 0;      // thousandths. note that if pos corr is on I found 130 to be right. without I found 125.
  // gen MC: pythia, *should remove smeared/weighted histograms.*
  bool matchmctruth = true;  // these two might be redundant
  bool filltruthspectrum = true;
  bool missingprimarypions = true; // if true, will combine truth photons until the mass is within 0.001 of the pi0 mass. these pions will be added to the pion spectrum

  TH2F* h_reco_etaphi;
  TH1F* h_pt1;
  TH1F* h_pt2;
  TH1F* h_nclusters;
  TH1F* h_mass_eta_lt[96];
  TH1F* h_pt_eta[96];
  TH1F* h_mass_eta_lt_rw[96];
  TH1F* h_pt_eta_rw[96];
  TH1F* h_pt_rw[96];
  TH1* h_matched_res;
  TH1* h_res_e;
  TH3* h_res_e_phi;
  TH3* h_res_e_eta;
  TH3* h_res_e_eta_pdc;
  TH1* h_res;
  TH2F* h_res_ptTr;
  TH3* h_m_pt_eta;
  TH3* h_m_ptTr_eta;
  TH3* h_m_ptTr_eta_trKin;
  TH3* h_delPhi_e_eta;
  TH3* h_delEta_e_eta;
  TH3* h_delPhi_e_phi;
  TProfile* pr_eta_shower;
  TProfile* pr_phi_shower;
  TH2F* h_vert_xy;
  TH1* h_truthE;

  TH2F* h_reco_etaphi_cuts[14];
  TH2F* h_cluster_etaphi_cuts[14];
  TH1F* h_vtxmap_fail;
  TH1F* h_zvtx;
  TH2F* h_vert_xy;
  TH1F* h_nevents;
  TH1F* h_sp_pt_rw;
  TH1F* h_emcal_e_eta;
  TH1F* h_truth_eta;
  TH1F* h_truth_phi;
  TH1F* h_truth_e;
  TH1F* h_truth_pt;
  TH1F* h_FullTruth_eta;
  TH1F* h_FullTruth_e;
  TH1F* h_FullTruth_pt;
  TH1F* h_FullTruth_p;
  TH1F* h_truth_spectrum1;
  TH1F* h_truth_spectrum2;
  TH1F* h_truth_spectrum3;
  TH1F* h_truth_spectrum4;
  TH1F* h_truth_spectrum5;
  TH1F* h_truth_etaspectrum;
  TH1F* h_truthmatched_mass1;
  TH1F* h_truthmatched_mass2;
  TH1F* h_truthmatched_mass3;
  TH1F* h_truthmatched_mass;
  TH2F* h_truthmatched_mass_2d;
  TH2F* h_truthmatched_mass_eta_2d;
  TH2F* h_truthmatched_mass_etameson_weighted_2d;
  TH1F* h_truthmatched_mass_etameson_weighted;
  TH3F* h_truthmatched_mass_etameson_weighted_eta_3d;
  TH3F* h_truthmatched_mass_etameson_weighted_etabin_3d;
  TH1F* h_pion_pt;
  TH1F* h_pion_pt_weight;
  TH1F* h_truth_pid_p;
  TH1F* h_truth_pid_s;
  TH2F* h_inv_yield;
  TH2F* h_yield;
  TH1F* h_smear_pi0E;
  TH1F* h_nosmear_pi0E;
  TH1F* h_smear_pi0E_weighted;
  TH1F* h_nosmear_pi0E_weighted;
  TH1F* h_smear_nosmear_pi0E;
  TH1F* h_smear_nosmear_pi0E_weighted;
  TH3F* h_InvMass_smear_weighted_asymmetry_3d = nullptr;
  TH3F* h_InvMass_photonE_smear_weighted_3d = nullptr;
  // mass vs eta
  TH3F* h_InvMass_smear_eta_3d = nullptr;
  TH3F* h_InvMass_smear_weighted_eta_3d = nullptr;
  TH2F* h_InvMass_smear_eta_2d = nullptr;
  TH2F* h_InvMass_smear_weighted_eta_2d = nullptr;

  TH1F* h_truthmatched_photon1E;
  TH1F* h_truthmatched_photon2E;
  TH1F* h_truthmatched_AllphotonE;
  TH1F* h_truth_ALLphotonE;
  TH1F* h_truth_ALLphotonE_weighted;
  TH1F* h_truth_ALLphotonpt;
  TH1F* h_truth_ALLphotonp;
  TH1F* h_truthmatched_photon1E_weighted;
  TH1F* h_truthmatched_photon2E_weighted;
  TH1F* h_truthmatched_AllphotonE_weighted;

  TH1F* h_reco_photon1E;
  TH1F* h_reco_photon2E;
  TH2F* h_reco_photon1E_2d;
  TH2F* h_reco_photon2E_2d;
  TH2F* h_reco_ALLphotonE_2d;
  TH1F* h_reco_ALLphotonE;
  TH1F* h_reco_photon1E_weighted;
  TH1F* h_reco_photon2E_weighted;
  TH1F* h_reco_ALLphotonE_weighted;
  TH1F* h_truthmatched_Photon_delR;

  TH2F* h_pi0_ELoss_2d;
  TH2F* h_pi0_ERatio_2d;
  TH2F* h_eta_ELoss_2d;
  TH2F* h_eta_ERatio_2d;
  TH2F* h_clus_ELoss_2d;
  TH2F* h_clus_ERatio_2d;
  TH1F* h_clusmultimatch;
  TH1F* h_ndecayphotons;

  TH1F* h_temp_pion_pt;
  TH1F* h_temp_pion_eta;
  TH1F* h_temp_pion_phi;
  TH1F* h_temp_pion_mass;
  TH1F* h_temp_pion_multimatch;
  TH1F* h_primaryphotonpair_massdiff;
  TH1F* h_primaryphotonpair_massdiff2; 

  float target_pi0_mass = 0.145;
  double truth_pt;
  double WeightScale = 1;  // e+14
  double inv_yield;
  TRandom3* rnd;
};

#endif
