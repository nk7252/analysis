#include "CaloAna.h"

// vertex includes
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>

// Fun4All includes
#include "CLHEP/Vector/ThreeVector.h"  // for Hep3Vector
// #include "CLHEP/Vector/ThreeVector.icc"      // for Hep3Vector::Hep3Vector
#include "fun4all/SubsysReco.h"  // for SubsysReco
#include "phool/phool.h"         // for PHWHERE
// #include <Event/Event.h>
// #include <Event/packet.h>
#include <cdbobjects/CDBTTree.h>  // for CDBTTree
// #include <ffamodules/CDBInterface.h>
// #include <ffamodules/FlagHandler.h>
// #include <ffamodules/HeadReco.h>
// #include <ffamodules/SyncReco.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
// #include <g4main/PHG4Hit.h>
// #include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
// #include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
// #include <phool/recoConsts.h>
#include "g4main/PHG4VtxPoint.h"

/// Jet includes
#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>

// Calorimeter/Cluster includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
// #include <calobase/RawTower.h>
// #include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
// #include <calobase/TowerInfoContainerv1.h>
// #include <calobase/TowerInfoContainerv2.h>
// #include <calobase/TowerInfoContainerv3.h>
#include <calobase/TowerInfoDefs.h>
// #include <calobase/TowerInfov1.h>
// #include <calobase/TowerInfov2.h>

// MBD
// #include <mbd/BbcGeom.h>
// #include <mbd/MbdPmtContainerV1.h>
// #include <mbd/MbdPmtHit.h>

// ROOT includes
// #include <Math/Vector4D.h>  // for ROOT::Math::PtEtaPhiMVector
// #include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
// #include <TH1F.h>
#include <TH2.h>
#include <TH3.h>
// #include <TH3F.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TNtuple.h>
// #include <TProfile.h>
#include <TAxis.h>  // for TAxis
#include <TRandom3.h>
#include <TString.h>  // for Form
#include <TTree.h>

// general includes
// #include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <set>
// #include <sstream>
#include <Rtypes.h>   // for R__LOAD_LIBRARY
#include <stdlib.h>   // for exit, size_t
#include <iterator>   // for next
#include <stdexcept>  // for runtime_error
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

/// HEPMC truth includes
// #pragma GCC diagnostic push
// #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
// #include <HepMC/GenEvent.h>
// #include <HepMC/GenVertex.h>
// #pragma GCC diagnostic pop
// #include <phhepmc/PHHepMCGenEvent.h>
// #include <phhepmc/PHHepMCGenEventMap.h>

R__LOAD_LIBRARY(libLiteCaloEvalTowSlope.so)

using namespace std;

CaloAna::CaloAna(const std::string& name, const std::string& filename)
  : SubsysReco(name)
  , detector("CEMC")  // HCALIN
  , outfilename(filename)
{
  _eventcounter = 0;
  std::random_device rd;  // Obtain a random number from hardware
  rng.seed(rd());         // Seed the generator
}

CaloAna::~CaloAna()
{
  delete hm;
  delete g4hitntuple;
  delete g4cellntuple;
  delete towerntuple;
  delete clusterntuple;
}

int CaloAna::Init(PHCompositeNode*)
{
  if (debug) std::cout << " " << "Init Start " << std::endl;
  hm = new Fun4AllHistoManager(Name());
  // create and register your histos (all types) here
  //////////////////////////
  // pT rewieghting
  if (SPMC_bool)
  {
    TH1F* h_original = nullptr;

    if (eta_weight)
    {
      frw = TFile::Open("/sphenix/user/nkumar/analysis/EMCal_pi0_Calib_2023/macros/seta_spectrum.root", "READ");
      h_original = (TH1F*) frw->Get("seta_pt_spectrum");
      std::cout << "Using eta spectrum" << std::endl;
    }
    else
    {
      frw = TFile::Open("/sphenix/user/nkumar/analysis/EMCal_pi0_Calib_2023/macros/spi0_spectrum.root", "READ");
      h_original = (TH1F*) frw->Get("spi0_pt_spectrum");
      std::cout << "Using pi0 spectrum" << std::endl;
    }

    if (h_original == nullptr)
    {
      std::cerr << "Error: Original histogram not found!" << std::endl;
    }
    else
    {
      h_sp_pt_rw = (TH1F*) h_original->Clone("h_sp_pt_rw_clone");
      if (h_sp_pt_rw == nullptr)
      {
        std::cerr << "Error: Cloning of histogram failed!" << std::endl;
      }
      else
      {
        h_sp_pt_rw->SetDirectory(0);
      }
    }

    // Close the file and clean up
    frw->Close();
    delete frw;
  }

  outfile = new TFile(outfilename.c_str(), "RECREATE");
  // if (SPMC_bool) h_sp_pt_rw->SetDirectory(outfile);  // Attach pt weight to outfile
  //  cutQA
  h_reco_etaphi = new TH2F("h_reco_etaphi", "Reco etaphi clusters", 256, -1 * TMath::Pi(), TMath::Pi(), 96, -1.2, 1.2);
  std::string vtxtype;
  if (MBDvtx)
    vtxtype = "MBD";
  else
    vtxtype = "Global";
  h_vtxmap_fail = new TH1F(Form("h_vtxmap_fail_%s", vtxtype.c_str()), Form("%s Vertex Map Failure type", vtxtype.c_str()), 3, -0.5, 2.5);
  h_vtxmap_fail->GetXaxis()->SetBinLabel(0, "Missing Vertex Map");
  h_vtxmap_fail->GetXaxis()->SetBinLabel(1, "Empty Vertex Map");
  h_vtxmap_fail->GetXaxis()->SetBinLabel(2, "No Vertex");
  h_zvtx = new TH1F(Form("h_zvtx_%s", vtxtype.c_str()), Form("%s Z Vertex; z vtx (cm)", vtxtype.c_str()), 100, -100, 100);
  h_vert_xy = new TH2F("h_vert_xy", "Vertex XY", 500, -120, 120, 500, -120, 120);
  h_nevents = new TH1F("h_nevents", "Number of events", 1, 0.5, 1.5);
  h_nevents->GetXaxis()->SetBinLabel(1, "N Events");

  // if (Cluster_Debug)
  //{
  for (int i = 0; i < 14; i++)
  {
    h_reco_etaphi_cuts[i] = new TH2F(Form("h_reco_etaphi_cuts%d", i), Form("h_reco_etaphi_cuts%d", i), 96, -1.2, 1.2, 256, -1 * TMath::Pi(), TMath::Pi());
    h_cluster_etaphi_cuts[i] = new TH2F(Form("h_cluster_etaphi_cuts%d", i), Form("h_cluster_etaphi_cuts%d", i), 96, -1.2, 1.2, 256, -1 * TMath::Pi(), TMath::Pi());
  }
  //}

  h_cutCounter = new TH1F("h_cutCounter", "Cut Counter", 13, 0.5, 13.5);
  // list of cuts
  //  clus1 chi2, clus1 cuts, tower eta>95,hotclus1,clus1=clus2,clus2 chi2, clus2 cuts, assym, Dr, pi0pt cut, hotclus2

  // correlation plots
  if (etabyeta)
  {
    for (int i = 0; i < 96; i++)
    {
      h_mass_eta_lt[i] = new TH1F(Form("h_mass_eta_lt%d", i), "", 50, 0, 0.5);
      h_mass_eta_lt_rw[i] = new TH1F(Form("h_mass_eta_lt_rw%d", i), "", 50, 0, 0.5);
      h_pt_eta[i] = new TH1F(Form("h_pt_eta%d", i), "", 100, 0, 10);
      h_pt_eta_rw[i] = new TH1F(Form("h_pt_eta_rw%d", i), "", 100, 0, 10);
    }
  }

  h_cemc_etaphi = new TH2F("h_cemc_etaphi", "", 96, 0, 96, 256, 0, 256);

  // 1D distributions
  // h_InvMass = new TH1F("h_InvMass", "Invariant Mass", 600, 0, 1.2);
  h_InvMass_w = new TH1F("h_InvMass_w", "Invariant Mass", 600, 0, 1.2);
  h_InvMassMix = new TH1F("h_InvMassMix", "Invariant Mass", 600, 0, 1.2);

  // cluster QA
  h_etaphi_clus = new TH2F("h_etaphi_clus", "", 256, -1 * TMath::Pi(), TMath::Pi(), 96, -1.2, 1.2);
  h_clusE = new TH1F("h_clusE", "", 100, 0, 50);
  h_emcal_e_eta = new TH1F("h_emcal_e_eta", "", 96, 0, 96);
  // h_pt1 = new TH1F("h_pt1", "", 100, 0, 20);
  // h_pt2 = new TH1F("h_pt2", "", 100, 0, 20);
  h_pion_pt = new TH1F("h_pion_pt", "", 100, 0, 50);
  h_pion_pt_weight = new TH1F("h_pion_pt_weight", "", 100, 0, 50);
  h_nclusters = new TH1F("h_nclusters", "", 1000, 0, 1000);

  // Truth histos
  h_truth_eta = new TH1F("h_truth_eta", "", 100, -1.2, 1.2);
  h_truth_e = new TH1F("h_truth_e", "", 100, 0, 50);
  h_truth_pt = new TH1F("h_truth_pt", "", 100, 0, 50);
  h_truth_spectrum1 = new TH1F("h_truth_spectrum1", "", 5 * 10, 0, 50);
  h_truth_spectrum2 = new TH1F("h_truth_spectrum2", "", 5 * 10, 0, 50);
  h_truth_spectrum3 = new TH1F("h_truth_spectrum3", "", 5 * 10, 0, 50);
  h_truth_spectrum4 = new TH1F("h_truth_spectrum4", "primary and secondary pions", 5 * 10, 0, 50);
  h_truth_spectrum5 = new TH1F("h_truth_spectrum5", "primary + matched (2 truth photons) truth pions", 5 * 10, 0, 50);
  h_truth_etaspectrum = new TH1F("h_truth_etaspectrum", "", 5 * 10, 0, 50);
  h_truth_pid_p = new TH1F("h_truth_pid_p", "Primary particle PIDs", 1000, -500, 500);
  h_truth_pid_s = new TH1F("h_truth_pid_s", "Secondary particle PIDs", 1000, -500, 500);
  h_delR_recTrth = new TH1F("h_delR_recTrth", "", 1000, 0, 5);
  h_delR_pionrecTrth = new TH1F("h_delR_pionrecTrth", "", 5000, 0, 5);
  h_matched_res = new TH2F("h_matched_res", "", 100, 0, 1.5, 20, -1, 1);
  h_FullTruth_eta = new TH1F("h_FullTruth_eta", "Full Truth eta", 100, -1.2, 1.2);
  h_FullTruth_e = new TH1F("h_FullTruth_e", "Full Truth e", 10 * 10, 0, 50);
  h_FullTruth_pt = new TH1F("h_FullTruth_pt", "Full Truth pt", 10 * 10, 0, 50);
  h_FullTruth_p = new TH1F("h_FullTruth_p", "Full Truth p", 10 * 10, 0, 50);

  h_temp_pion_pt = new TH1F("h_temp_pion_pt", "missing primary truth pion candidates, pT", 5 * 10, 0, 50);
  h_temp_pion_eta = new TH1F("h_temp_pion_eta", "missing primary truth pion candidates, Eta", 96, -1.2, 1.2);
  h_temp_pion_phi = new TH1F("h_temp_pion_phi", "missing primary truth pion candidates, Phi", 256, -1 * TMath::Pi(), TMath::Pi());
  h_temp_pion_mass = new TH1F("h_temp_pion_mass", "missing primary truth pion candidates, Mass", 600, 0, 1.2);
  h_temp_pion_multimatch = new TH1F("h_temp_pion_multimatch", "missing primary truth pion candidates, multiple matches", 3, 0, 3);
  h_primaryphotonpair_massdiff = new TH1F("h_primaryphotonpair_massdiff", "mass difference between primary photon pairs", 100, -0.00003, 1);
  h_primaryphotonpair_massdiff2 = new TH1F("h_primaryphotonpair_massdiff2", "mass difference between primary photon pairs, zoomed in", 100, -0.00003, 0.00003);

  // histograms to extract MC photon resolution
  h_truthmatched_AllphotonE = new TH1F("h_truthmatched_AllphotonE", "All Photon Energy", 20 * 10, 0, 50);
  h_truth_ALLphotonE = new TH1F("h_truth_ALLphotonE", "All Photon Energy", 20 * 10, 0, 50);
  h_truthmatched_AllphotonE_weighted = new TH1F("h_truthmatched_AllphotonE_weighted", "All Photon Energy, weighted", 20 * 10, 0, 50);
  // h_truth_ALLphotonE_weighted = new TH1F("h_truth_ALLphotonE_weighted", "All Photon Energy, weighted", 20 * 10, 0, 50);
  // reco photons
  h_reco_ALLphotonE_2d = new TH2F("h_reco_ALLphotonE_2d", "pT vs All Reco Photon Energy", 20 * 10, 0, 50, 20 * 10, 0, 50);
  h_reco_ALLphotonE = new TH1F("h_reco_ALLphotonE", "All Reco Photon Energy", 20 * 10, 0, 50);
  // All truth photons(secondaries)
  h_truth_ALLphotonpt = new TH1F("h_truth_ALLphotonpt", "All Photon pt", 5 * 10, 0, 50);
  h_truth_ALLphotonp = new TH1F("h_truth_ALLphotonp", "All Photon p", 5 * 10, 0, 50);
  // truthmatching histograms from blair
  ///*
  h_res_e = new TH2F("h_res_e", "", 750, 0, 1.5, 50, 0, 50);
  // h_res_e_phi = new TH3F("h_res_e_phi","",100,0,1.5,10,0,20,256,0,256);
  // h_res_e_eta = new TH3F("h_res_e_eta","",300,0,1.5,40,0,20,96,0,96);
  // h_m_pt_eta = new TH3F("h_m_pt_eta","",70,0,0.7,10,0,10,96,0,96);
  h_m_ptTr_eta = new TH3F("h_m_ptTr_eta", "", 70, 0, 0.7, 10, 0, 10, 96, 0, 96);
  h_m_ptTr_eta_trKin = new TH3F("h_m_ptTr_eta_trKin", "", 70, 0, 0.7, 10, 0, 10, 96, 0, 96);
  h_res = new TH1F("h_res", "", 750, 0, 1.5);
  h_res_ptTr = new TH2F("h_res_ptTr", "", 20 * 10, 0, 50, 750, 0, 1.5);
  // h_delEta_e_eta = new TH3F("h_delEta_e_eta","",100,-0.1,0.1,10,0,20,96,0,96);
  // h_delPhi_e_eta = new TH3F("h_delPhi_e_eta","",100,-0.3,0.3,20,0,20,96,0,96);
  // h_delPhi_e_phi = new TH3F("h_delPhi_e_phi","",100,-0.1,0.1,20,0,20,256,0,256);
  // pr_eta_shower = new TProfile("pr_eta_shower","",96,-48.5,47.5, -1,1.5);
  // pr_phi_shower = new TProfile("pr_phi_shower","",256,-128.5,127.5, -1,1.5);
  // h_vert_xy = new TH2F("h_vert_xy","",500,-120,120,500,-120,120);
  h_truthE = new TH1F("h_truthE", "", 10000, 0, 30);
  h_pi0_ELoss_2d = new TH2F("h_pi0_ELoss_2d", "pT vs Eloss;Truth Meson pT (GeV);Truth_E - Reco_E(GeV)", 20 * 10, 0, 50, 2000, -10, 10);
  h_pi0_ERatio_2d = new TH2F("h_pi0_ERatio_2d", "pT vs reco/truth ERatio;Truth Meson pT (GeV);Meson ERatio(GeV)", 20 * 10, 0, 50, 1000, 0, 10);
  h_eta_ELoss_2d = new TH2F("h_eta_ELoss_2d", "pT vs Eloss;Truth Meson pT (GeV);Truth_E - Reco_E(GeV)", 20 * 10, 0, 50, 2000, -10, 10);
  h_eta_ERatio_2d = new TH2F("h_eta_ERatio_2d", "pT vs reco/truth ERatio;Truth Meson pT (GeV);Meson ERatio(GeV)", 20 * 10, 0, 50, 1000, 0, 10);
  h_clus_ELoss_2d = new TH2F("h_clus_ELoss_2d", "Cluster Eloss;Truth CLuster pT (GeV);Truth_E - Reco_E(GeV)", 20 * 10, 0, 50, 2000, -10, 10);
  h_clus_ELoss_weighted_2d = new TH2F("h_clus_ELoss_weighted_2d", "Cluster Eloss, weighted;Truth CLuster pT (GeV);Truth_E - Reco_E(GeV)", 20 * 10, 0, 50, 2000, -10, 10);
  h_clus_ERatio_2d = new TH2F("h_clus_ERatio_2d", "reco/truth Cluster ERatio;Truth CLuster pT (GeV);Cluster ERatio(GeV)", 20 * 10, 0, 50, 1000, 0, 100);
  h_clus_ERatio_weighted_2d = new TH2F("h_clus_ERatio_weighted_2d", "reco/truth Cluster ERatio, weighted;Truth CLuster pT (GeV);Cluster ERatio(GeV)", 20 * 10, 0, 50, 1000, 0, 100);

  h_clusmultimatch = new TH1F("h_clusmultimatch", "Cluster multi match", 11, -0.5, 10.5);
  h_ndecayphotons = new TH1F("h_ndecayphotons", "Number of decay photons", 5, -0.5, 4.5);
  //*/

  h_reco_photon1E_weighted = new TH1F("h_reco_photon1E_weighted", "Reco Photon 1 Energy, weighted", 20 * 10, 0, 50);
  h_reco_photon2E_weighted = new TH1F("h_reco_photon2E_weighted", "Reco Photon 2 Energy, weighted", 20 * 10, 0, 50);
  h_reco_ALLphotonE_weighted = new TH1F("h_reco_ALLphotonE_weighted", "All Reco Photon Energy, weighted", 20 * 10, 0, 50);

  // pT differential Inv Mass
  h_InvMass = new TH1F("h_InvMass", "Invariant Mass", 600, 0, 1.2);
  h_InvMass_2d = new TH2F("h_InvMass_2d", "pT vs Invariant Mass", 20 * 10, 0, 50, 600, 0, 1.2);
  h_InvMass_weighted = new TH1F("h_InvMass_weighted", "Invariant Mass, weighted WSHP", 600, 0, 1.2);

  h_inv_yield = new TH2F("h_inv_yield", "Invariant Yield distribution", 20 * 10, 0, 50, 100, 0, 20);
  h_yield = new TH2F("h_yield", "Yield distribution", 20 * 10, 0, 50, 100, 0, 70);
  h_truthmatched_mass = new TH1F("h_truthmatched_mass", "Invariant Mass, ph truth matched", 600, 0, 1.2);
  h_truthmatched_mass_2d = new TH2F("h_truthmatched_mass_2d", "pT vs Invariant Mass, ph truth matched", 20 * 10, 0, 50, 600, 0, 1.2);
  h_truthmatched_mass_eta_2d = new TH2F("h_truthmatched_mass_eta_2d", "eta vs Invariant Mass, ph truth matched", 24, -1.2, 1.2, 600, 0, 1.2);
  h_truthmatched_mass1 = new TH1F("h_truthmatched_mass1", "Invariant Mass, truth matched(delR<0.015)", 600, 0, 1.2);
  h_truthmatched_mass2 = new TH1F("h_truthmatched_mass2", "Invariant Mass, truth matched(delR<0.1)", 600, 0, 1.2);
  h_truthmatched_mass3 = new TH1F("h_truthmatched_mass3", "Invariant Mass, truth matched(delR<0.2)", 600, 0, 1.2);
  h_truthmatched_mass1_2d = new TH2F("h_truthmatched_mass1_2d", "pT vs Invariant Mass, truth matched(delR<0.015)", 20 * 10, 0, 50, 600, 0, 1.2);
  h_truthmatched_mass2_2d = new TH2F("h_truthmatched_mass2_2d", "pT vs Invariant Mass, truth matched(delR<0.1)", 20 * 10, 0, 50, 600, 0, 1.2);
  h_truthmatched_mass3_2d = new TH2F("h_truthmatched_mass3_2d", "pT vs Invariant Mass, truth matched(delR<0.2)", 20 * 10, 0, 50, 600, 0, 1.2);
  h_truthmatched_mass_weighted = new TH1F("h_truthmatched_mass_weighted", "Invariant Mass, truth matched, weighted", 600, 0, 1.2);
  h_truthmatched_mass_weighted_2d = new TH2F("h_truthmatched_mass_weighted_2d", "pT vs Invariant Mass, truth matched, weighted", 20 * 10, 0, 50, 600, 0, 1.2);
  h_truthmatched_mass_etameson_weighted_2d = new TH2F("h_truthmatched_mass_etameson_weighted_2d", "pT vs Invariant Mass, truth matched, weighted", 20 * 10, 0, 50, 600, 0, 1.2);
  h_truthmatched_mass_etameson_weighted = new TH1F("h_truthmatched_mass_etameson_weighted", "Invariant Mass, truth matched, weighted", 600, 0, 1.2);
  // h_truthmatched_mass_etameson_weighted_etabin_3d = new TH3F("h_truthmatched_mass_etameson_weighted_etabin_3d", "pT vs Invariant Mass vs eta (bin), truth matched, weighted", 20 * 10, 0, 50, 600, 0, 1.2, 96, 0, 96);
  h_truthmatched_mass_etameson_weighted_eta_3d = new TH3F("h_truthmatched_mass_etameson_weighted_eta_3d", "pT vs Invariant Mass vs eta, truth matched, weighted", 20 * 10, 0, 50, 600, 0, 1.2, 96, -1.2, 1.2);

  // 3d histogram to check for corelation between photon/cluster energies and invariant mass.
  // h_InvMass_photonE_smear_weighted_3d = new TH3F(Form("h_InvMass_smear%d_weighted_photonE_3d",badcalibsmearint ), Form("Photon Energies vs Invariant Mass, smear, weighted: %f percent",badcalibsmearint / 10.0f), 40, 0, 20, 40, 0, 20, 60, 0, 1.2);
  // 3d histogram to check for for corelation between pt invariant mass and asymmetry
  // h_InvMass_smear_weighted_asymmetry_3d = new TH3F(Form("h_InvMass_smear%d_weighted_asymmetry_3d", badcalibsmearint ),Form("pT vs Invariant Mass vs asymmetry + smear, weighted: %f percent", badcalibsmearint / 10.0f), 20 * 10, 0, 50, 60, 0, 1.2, 10, 0, 1);

  // 3d histogram to check for corelation between eta, pt and invariant mass
  h_InvMass_smear_eta_3d = new TH3F(Form("h_InvMass_smear%d_eta_3d", badcalibsmearint), Form("pT vs Invariant Mass vs eta + smear: %f percent", badcalibsmearint / 10.0f), 20 * 10, 0, 50, 60, 0, 1.2, 24, -1.2, 1.2);
  h_InvMass_smear_weighted_eta_3d = new TH3F(Form("h_InvMass_smear%d_weighted_eta_3d", badcalibsmearint), Form("pT vs Invariant Mass vs eta + smear, weighted: %f percent", badcalibsmearint / 10.0f), 20 * 10, 0, 50, 60, 0, 1.2, 24, -1.2, 1.2);
  // 2d histogram to check for corelation between eta, and invariant mass
  h_InvMass_smear_eta_2d = new TH2F(Form("h_InvMass_smear%d_eta_2d", badcalibsmearint), Form("eta vs Invariant Mass+ smear: %f percent", badcalibsmearint / 10.0f), 24, -1.2, 1.2, 120, 0, 1.2);
  h_InvMass_smear_weighted_eta_2d = new TH2F(Form("h_InvMass_smear%d_weighted_eta_2d", badcalibsmearint), Form("eta vs Invariant Mass+ smear, weighted: %f percent", badcalibsmearint / 10.0f), 24, -1.2, 1.2, 120, 0, 1.2);

  // for (int i = 0; i < 96; i++) h_pt_rw[i] = (TH1F*) frw->Get(Form("h_pt_eta%d", i));
  // these histograms below were not working when using SPMC bool
  rnd = new TRandom3();
  // smearing added SPMC
  badcalibsmear = static_cast<float>(badcalibsmearint) / 1000.0f;

  h_InvMass_smear = new TH1F(Form("h_InvMass_smear_%d", badcalibsmearint), Form("Invariant Mass + const smear: %f percent", badcalibsmearint / 10.0f), 600, 0, 1.2);

  h_InvMass_smear_2d = new TH2F(Form("h_InvMass_smear_2d_%d", badcalibsmearint), Form("pT vs Invariant Mass + const smear: %f percent", badcalibsmearint / 10.0f), 20 * 10, 0, 50, 600, 0, 1.2);

  // weighted variants
  h_InvMass_smear_weighted = new TH1F(Form("h_InvMass_smear_weighted_%d", badcalibsmearint), Form("Invariant Mass + const smear, weighted: %f percent", badcalibsmearint / 10.0f), 600, 0, 1.2);

  h_InvMass_smear_weighted_2d = new TH2F(Form("h_InvMass_smear_weighted_2d_%d", badcalibsmearint), Form("pT vs Invariant Mass + const smear, weighted: %f percent", badcalibsmearint / 10.0f), 20 * 10, 0, 50, 600, 0, 1.2);

  // towers node selection
  calotowerinfostring = (clust_waveform == true) ? "WAVEFORM_CEMC" : "TOWERINFO_CALIB_CEMC";
  // clustercontainer node selection//&& SPMC_bool == false
  clustcontainerstring = (poscor == true) ? "CLUSTER_POS_COR_CEMC" : ((recluster == true) ? "CLUSTERINFO_CEMC2" : ((clust_waveform == true) ? "CLUSTERINFO_CEMC" : "CLUSTERINFO_CEMC"));
  std::cout << "clustcontainerstring: " << clustcontainerstring << std::endl;
  std::cout << "calotowerinfostring: " << calotowerinfostring << std::endl;

  funkyCaloStuffcounter = 0;
  VertexMapFailcounter = 0;
  if (additionalsmearing == false) std::cout << "additional smearing is not being added" << std::endl;
  if (additionalsmearing == true) std::cout << "additional smearing is being added" << std::endl;
  if (debug) std::cout << " " << "CaloAna: End Init  " << std::endl;
  // return 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloAna::process_event(PHCompositeNode* topNode)
{
  _eventcounter++;
  if (debug) std::cout << " " << "Process Event: Start  " << std::endl;
  h_nevents->Fill(1);
  process_towers(topNode);
  if (debug) std::cout << " " << "Process Event: End  " << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloAna::process_towers(PHCompositeNode* topNode)
{
  if ((_eventcounter % 1000) == 0) std::cout << _eventcounter << std::endl;

  // Declare tracking sets for each cut
  std::unordered_set<RawCluster*> filledClustersAfterCut1;
  std::unordered_set<RawCluster*> filledClustersAfterCut2;
  std::unordered_set<RawCluster*> filledClustersAfterCut3;
  std::unordered_set<RawCluster*> filledClustersAfterCut4;
  std::unordered_set<RawCluster*> filledClustersAfterCut5;
  std::unordered_set<RawCluster*> filledClustersAfterCut6;
  std::unordered_set<RawCluster*> filledClustersAfterCut7;
  std::unordered_set<RawCluster*> filledClustersAfterCut8;
  std::unordered_set<RawCluster*> filledClustersAfterCut9;
  std::unordered_set<RawCluster*> filledClustersAfterCut10;
  std::unordered_set<RawCluster*> filledClustersAfterCut11;
  std::unordered_set<RawCluster*> filledClustersAfterCut12;

  // cuts
  if (debug) std::cout << " " << "Cuts " << std::endl;
  float maxDr = 100;          // cone cut, 100 is effectively off
  //float maxAlpha = 0.6;       // asymmetry cut
  //float clus_chisq_cut = 10;  // normally 4
  //float clusterprob = 0.1;    // replacing chisqr cut
  float nClus_ptCut = 0.0;    // 0.5 normally
  float pi0ptcutfactor = 0;   // seto to 0 to effectively disable it
  float ptMaxCut = 100;       // no cut in data, as far as I know. so I set it to a value it is unlikely to reach
  float pt1ClusCut = cluspTcut.first;//= 1.0;     // centrality dependence cuts 2.2 for both // 1.3
  float pt2ClusCut = cluspTcut.second;//0.6;     // 0.7
  float etcut = 1.0;          // cluster ET cut
  //float etacutval = 0.6;      // cluster pseudo-rapidity cut
  //float zvtx_cut_val = 30;    // z vertex cut value

  /*
  if (nClusCount > 30)
  {
    pt1ClusCut += 1.4 * (nClusCount - 29) / 200.0;
    pt2ClusCut += 1.4 * (nClusCount - 29) / 200.0;
  }
  //*/

  float pi0ptcut = pi0ptcutfactor * (pt1ClusCut + pt2ClusCut);

  // int max_nClusCount = 75;

  //-----------------------get vertex----------------------------------------//

  float vtx_z = 0;
  float vtx_x = 0;
  float vtx_y = 0;
  MbdVertexMap* mbdvertexmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
  GlobalVertexMap* vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (getVtx)
  {
    if ((!MBDvtx && !vertexmap) || (MBDvtx && !mbdvertexmap))
    {
      // if (debug) std::cout << PHWHERE << " Fatal Error - GlobalVertexMap node is missing"<< std::endl;
      if (!MBDvtx) std::cout << "CaloAna GlobalVertexMap node is missing" << std::endl;
      if (MBDvtx) std::cout << "CaloAna MbdVertexMap node is missing" << std::endl;
      h_vtxmap_fail->Fill(0);
      VertexMapFailcounter++;
    }
    if (MBDvtx)
    {
      if (mbdvertexmap && !mbdvertexmap->empty())
      {
        MbdVertex* mbdvtx = mbdvertexmap->begin()->second;
        if (mbdvtx)
        {
          vtx_z = mbdvtx->get_z();
          h_zvtx->Fill(vtx_z);
          if (debug) std::cout << "vtx_z: " << vtx_z << std::endl;
        }
        else
        {
          if (debug) std::cout << "CaloAna MBDVertex node returns no vtx" << std::endl;
          h_vtxmap_fail->Fill(2);
          VertexMapFailcounter++;
        }
      }
      else
      {
        if (debug) std::cout << "CaloAna GlobalVertexMap node is empty" << std::endl;  // if (debug)
        h_vtxmap_fail->Fill(1);
        VertexMapFailcounter++;
      }
    }
    else
    {
      if (vertexmap && !vertexmap->empty())
      {
        GlobalVertex* vtx = vertexmap->begin()->second;
        if (vtx)
        {
          vtx_z = vtx->get_z();
          h_zvtx->Fill(vtx_z);
          if (debug) std::cout << "vtx_z: " << vtx_z << std::endl;
        }
        else
        {
          if (debug) std::cout << "CaloAna GlobalVertex node returns no vtx" << std::endl;
          h_vtxmap_fail->Fill(2);
          VertexMapFailcounter++;
        }
      }
      else
      {
        if (debug) std::cout << "CaloAna GlobalVertexMap node is empty" << std::endl;  // if (debug)
        h_vtxmap_fail->Fill(1);
        VertexMapFailcounter++;
      }
    }
  }

  if (zvtxcut_bool && abs(vtx_z) > zvtx_cut_val)
  {
    h_cutCounter->Fill(12);
    return Fun4AllReturnCodes::EVENT_OK;
  }

  if (SPMC_bool)
  {
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    PHG4TruthInfoContainer::VtxRange vtxrange = truthinfo->GetVtxRange();
    int n_vertex = 0;
    float vertex_x[1000] = {0};
    float vertex_y[1000] = {0};
    float vertex_z[1000] = {0};
    float vertex_id[1000] = {0};
    for (PHG4TruthInfoContainer::ConstVtxIterator iter = vtxrange.first; iter != vtxrange.second; ++iter)
    {
      PHG4VtxPoint* vtx = iter->second;
      // if ( n_vertex > 0) continue;
      vertex_x[n_vertex] = vtx->get_x();
      vertex_y[n_vertex] = vtx->get_y();
      vertex_z[n_vertex] = vtx->get_z();
      vertex_id[n_vertex] = vtx->get_id();
      // if (n_vertex < 10) std::cout << "vx=" << vertex_x[n_vertex] << "  vy=" << vertex_y[n_vertex] << "   vz=" << vertex_z[n_vertex] << "  id=" << vertex_id[n_vertex] << std::endl;
      if (vertex_id[n_vertex] == 1)
      {
        // if (false) std::cout << "vx=" << vertex_x[n_vertex] << "  vy=" << vertex_y[n_vertex] << "   vz=" << vertex_z[n_vertex] << "  id=" << vertex_id[n_vertex] << std::endl;
        //  Pvtx_id = n_vertex;
        // std::cout << "truth vertex: " << vertex_x[n_vertex] << " " << vertex_y[n_vertex] << " " << vertex_z[n_vertex] << ", ID:" << vertex_id[n_vertex] << std::endl;
        h_vert_xy->Fill(vertex_x[n_vertex], vertex_y[n_vertex]);
        vtx_x = vertex_x[n_vertex];
        vtx_y = vertex_y[n_vertex];
        vtx_z = vertex_z[n_vertex];
      }
      n_vertex++;
      if (n_vertex >= 100000) break;
    }
    // std::cout << "truth vertex: " << vtx_x << " " << vtx_y << " " << vtx_z << ", ID:" << vertex_id[0] << std::endl;
  }

  //////////////////////////////////////////////
  //         towers
  vector<float> ht_eta;
  vector<float> ht_phi;

  // if (!m_vtxCut || abs(vtx_z) > _vz)  return Fun4AllReturnCodes::EVENT_OK;

  TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, Form("%s", calotowerinfostring.c_str()));

  RawClusterContainer* clusterContainer = findNode::getClass<RawClusterContainer>(topNode, Form("%s", clustcontainerstring.c_str()));
  // changed from CLUSTERINFO_CEMC2
  // Blair using "CLUSTER_POS_COR_CEMC" now. change from CLUSTER_CEMC
  // RawClusterContainer* clusterContainer = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_POS_COR_CEMC");
  if (!clusterContainer)
  {
    if (debug) std::cout << PHWHERE << "funkyCaloStuff::process_event - Fatal Error - RawClusterContainer node is missing. " << std::endl;
    funkyCaloStuffcounter++;
    return 0;
  }

  //////////////////////////////////////////
  // geometry for hot tower/cluster masking
  std::string towergeomnodename = "TOWERGEOM_CEMC";
  RawTowerGeomContainer* m_geometry = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodename);
  if (!m_geometry)
  {
    std::cout << Name() << "::"
              << "CreateNodeTree"
              << ": Could not find node " << towergeomnodename << std::endl;
    throw std::runtime_error("failed to find TOWERGEOM node in RawClusterDeadHotMask::CreateNodeTree");
  }


  
  /// Get reco jets
  /*
  JetContainer* reco_jets = findNode::getClass<JetContainer>(topNode, "AntiKt_Tower_r04");  //"AntiKt_Tower_r04_Sub1"recoJetName
  if (!reco_jets)
  {
    std::cout
        << "CaloAna::process_event - Error: can not find Reco JetContainer node:" << "AntiKt_Tower_r04" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
    //exit(-1);
  }
  //*/

  /// Get truth jets
  JetContainer *truth_jets = findNode::getClass<JetContainer>(topNode, "AntiKt_Truth_r04"); // truthJetName
  if (!truth_jets)
  {
    std::cout
        << "MyJetAnalysis::process_event - Error: can not find Truth JetContainer node: "
        << "AntiKt_Truth_r04" << std::endl;
    //exit(-1);
    return Fun4AllReturnCodes::EVENT_OK;
  }



  bool range_efficient = false;
  // loop over jets and check if pT range is efficient
  if (pythiajets)
  {
    for (auto jet : *truth_jets)
    {
      if (jet->get_pt() > efficiencyrange.first && jet->get_pt() < efficiencyrange.second)
      {
        range_efficient = true;
      }
    }
    if (!range_efficient)
    {
      if(debug) std::cout << "CaloAna::process_event - Error: pT range not efficient for Jet Sample" << std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }

  RawClusterContainer::ConstRange clusterEnd = clusterContainer->getClusters();
  RawClusterContainer::ConstIterator clusterIter;
  RawClusterContainer::ConstIterator clusterIter2;
  int nClusCount = 0;
  for (clusterIter = clusterEnd.first; clusterIter != clusterEnd.second; clusterIter++)
  {
    RawCluster* recoCluster = clusterIter->second;

    CLHEP::Hep3Vector vertex(vtx_x, vtx_y, vtx_z);
    CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*recoCluster, vertex);

    float clus_pt = E_vec_cluster.perp();
    float clus_chisq = recoCluster->get_chi2();
    float prob = recoCluster->get_prob();

    if (clus_pt < nClus_ptCut && cutson) continue;
    if (clusterprobcut && prob < clusterprob && cutson)
      continue;
    else if (clus_chisq > clus_chisq_cut && cutson && !clusterprobcut)
      continue;
    nClusCount++;
  }

  h_nclusters->Fill(nClusCount);

  // if (nClusCount > max_nClusCount && cutson) return Fun4AllReturnCodes::EVENT_OK;

  vector<float> save_pt;
  vector<float> save_eta;
  vector<float> save_phi;
  vector<float> save_e;

  vector<TLorentzVector> truth_photons;
  vector<TLorentzVector> truth_pi0_photons;

  if (filltruthspectrum)
  {
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    if (truthinfo)
    {
      if (matchmctruth)
      {  // primaries
        std::unordered_set<const PHG4Particle*> truth_Prim_photons;
        std::unordered_set<const PHG4Particle*> used_photons;
        std::vector<std::pair<const PHG4Particle*, const PHG4Particle*>> primary_reco_pions;

        PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
        for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
        {
          const PHG4Particle* truth = iter->second;
          if (!truthinfo->is_primary(truth)) continue;
          if (truth->get_pid() == 111 && !eta_weight)
          {
            float pion_pt = sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py());
            float pion_p = sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py() + truth->get_pz() * truth->get_pz());
            float pion_e = truth->get_e();
            float pion_phi = atan2(truth->get_py(), truth->get_px());
            float pion_eta = atanh(truth->get_pz() / sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py() + truth->get_pz() * truth->get_pz()));
            TLorentzVector truthpi0 = TLorentzVector();
            truthpi0.SetPtEtaPhiE(pion_pt, pion_eta, pion_phi, pion_e);
            h_truth_spectrum2->Fill(truthpi0.Pt());
            h_truth_spectrum4->Fill(truthpi0.Pt());
            h_truth_spectrum5->Fill(truthpi0.Pt());
            if (abs(pion_eta) <= 1.0) h_truth_spectrum3->Fill(pion_pt);
            h_FullTruth_e->Fill(pion_e);
            h_FullTruth_eta->Fill(pion_eta);
            h_FullTruth_pt->Fill(pion_pt);
            h_FullTruth_p->Fill(pion_p);
            if (debug) std::cout << "primary pid=" << truth->get_pid() << "   E=" << pion_e << "  pt=" << pion_pt << "  eta=" << pion_eta << "  phi=" << pion_phi << std::endl;
            // truth_photons.push_back(truthpi0);
          }
          if (truth->get_pid() == 221)
          {
            h_truth_etaspectrum->Fill(sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py()));
            if (debug) std::cout << "Primary eta found at vertex: " << truth->get_vtx_id() << std::endl;
            if (eta_weight)
            {
              float eta_pt = sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py());
              float eta_p = sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py() + truth->get_pz() * truth->get_pz());
              float eta_e = truth->get_e();
              float eta_phi = atan2(truth->get_py(), truth->get_px());
              float eta_prapid = atanh(truth->get_pz() / sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py() + truth->get_pz() * truth->get_pz()));
              TLorentzVector truthetameson = TLorentzVector();
              truthetameson.SetPtEtaPhiE(eta_pt, eta_prapid, eta_phi, eta_e);
              h_truth_spectrum2->Fill(eta_pt);
              if (abs(eta_pt) <= 1.0) h_truth_spectrum3->Fill(eta_pt);
              h_FullTruth_e->Fill(eta_e);
              h_FullTruth_eta->Fill(eta_prapid);
              h_FullTruth_pt->Fill(eta_pt);
              h_FullTruth_p->Fill(eta_p);
              // truth_photons.push_back(truthetameson);
            }
          }
          // photon loop
          if (truth->get_pid() == 22)
          {
            truth_Prim_photons.insert(truth);
            float photon_e = truth->get_e();
            h_truth_ALLphotonE->Fill(photon_e);
            // what should the photon weight be?
            // h_truth_ALLphotonE_weighted->Fill(photon_e, inv_yield);
          }
        }

        if (missingprimarypions)
        {
          // Iterate over unique pairs of photons
          for (auto it1 = truth_Prim_photons.begin(); it1 != truth_Prim_photons.end(); ++it1)
          {
            for (auto it2 = std::next(it1); it2 != truth_Prim_photons.end(); ++it2)
            {
              const PHG4Particle* p1 = *it1;
              const PHG4Particle* p2 = *it2;
              TLorentzVector temp_pion = TLorentzVector();
              temp_pion.SetPxPyPzE(p1->get_px() + p2->get_px(), p1->get_py() + p2->get_py(), p1->get_pz() + p2->get_pz(), p1->get_e() + p2->get_e());
              // bool multicount = false;
              //  Skip pairs with particles that have already been used
              if (used_photons.count(p1) > 0 || used_photons.count(p2) > 0)
              {
                // multicount = true;
                h_temp_pion_multimatch->Fill(1);
                continue;
              }
              float massdiff = fabs(temp_pion.M() - 0.135);
              h_primaryphotonpair_massdiff->Fill(massdiff);
              if (massdiff < 1) h_primaryphotonpair_massdiff2->Fill(massdiff);
              // Check if the pair's mass is near the target mass, accounting for floating point error
              if (massdiff < 0.00003)  //( temp_pion.M() == 0.135)//0.13497
              {
                primary_reco_pions.emplace_back(p1, p2);  // Store the pair
                used_photons.insert(p1);                  // Mark photons as used
                used_photons.insert(p2);

                h_temp_pion_pt->Fill(temp_pion.Pt());
                h_temp_pion_eta->Fill(temp_pion.Eta());
                h_temp_pion_phi->Fill(temp_pion.Phi());
                h_temp_pion_mass->Fill(temp_pion.M());
                h_truth_spectrum5->Fill(temp_pion.Pt());
                h_temp_pion_multimatch->Fill(2);
                // break;
              }
              else
                h_temp_pion_multimatch->Fill(0);
            }
          }
        }

        // secondaries
        PHG4TruthInfoContainer::Range second_range = truthinfo->GetSecondaryParticleRange();
        for (PHG4TruthInfoContainer::ConstIterator siter = second_range.first; siter != second_range.second; ++siter)
        {
          const PHG4Particle* truth = siter->second;
          int id = truth->get_pid();
          h_truth_pid_s->Fill(id);
          if (filltruthspectrum)
          {
            if (id == 22)
            {
              float photon_e = truth->get_e();
              float photon_pt = sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py());
              float photon_p = sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py() + truth->get_pz() * truth->get_pz());
              h_truth_ALLphotonE->Fill(photon_e);
              h_truth_ALLphotonpt->Fill(photon_pt);
              h_truth_ALLphotonp->Fill(photon_p);

              PHG4Particle* parent = truthinfo->GetParticle(truth->get_parent_id());
              if (eta_weight && parent->get_pid() == 221)
              {
                if (photon_pt < 0.1) continue;
                if (debug) std::cout << "Parent eta found at vertex: " << parent->get_vtx_id() << std::endl;
                float phot_phi = atan2(truth->get_py(), truth->get_px());
                float phot_eta = atanh(truth->get_pz() / sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py() + truth->get_pz() * truth->get_pz()));
                TLorentzVector myVector;
                myVector.SetXYZT(truth->get_px(), truth->get_py(), truth->get_pz(), truth->get_e());
                truth_pi0_photons.push_back(myVector);
                if (debug) std::cout << "2nd photons  pt=" << photon_pt << " e=" << photon_e << " phi=" << phot_phi << " eta=" << phot_eta << endl;
              }
              else if (!eta_weight && parent->get_pid() == 111)
              {
                if (photon_pt < 0.1) continue;
                float phot_phi = atan2(truth->get_py(), truth->get_px());
                float phot_eta = atanh(truth->get_pz() / sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py() + truth->get_pz() * truth->get_pz()));
                TLorentzVector myVector;
                myVector.SetXYZT(truth->get_px(), truth->get_py(), truth->get_pz(), truth->get_e());
                truth_pi0_photons.push_back(myVector);
                if (debug) std::cout << "2nd photons  pt=" << photon_pt << " e=" << photon_e << " phi=" << phot_phi << " eta=" << phot_eta << endl;
              }
            }

            if (id == 111)
            {
              h_truth_spectrum4->Fill(sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py()));
            }
          }
        }
      }
    }
  }

  // float smear = 0.00;
  bool match1 = false;
  bool match2 = false;
  int multimatchint = 0;
  if (debug) std::cout << " " << "Cluster Loop: 1 " << std::endl;
  for (clusterIter = clusterEnd.first; clusterIter != clusterEnd.second; clusterIter++)
  {
    RawCluster* recoCluster = clusterIter->second;

    CLHEP::Hep3Vector vertex(vtx_x, vtx_y, vtx_z);
    CLHEP::Hep3Vector E_vec_cluster;
    if (pp_rawcluster)
    {
      E_vec_cluster = RawClusterUtility::GetEVec(*recoCluster, vertex);
    }
    else if (!pp_rawcluster)  // i.e if AuAu
    {
      E_vec_cluster = RawClusterUtility::GetECoreVec(*recoCluster, vertex);
    }
    // GetEVec? look at rawclusterv1 maybe?
    std::vector<TLorentzVector> pi0smearvec(3);

    float clusE = E_vec_cluster.mag();
    float clus_eta = E_vec_cluster.pseudoRapidity();
    float clus_phi = E_vec_cluster.phi();
    float clus_pt = E_vec_cluster.perp();
    float prob = recoCluster->get_prob();
    // clus_pt *= rnd->Gaus(1, smear);
    float clus_chisq = recoCluster->get_chi2();
    h_cluster_etaphi_cuts[1]->Fill(clus_eta, clus_phi);

    // int lt_eta = recoCluster->get_lead_tower().first;
    // int lt_phi = recoCluster->get_lead_tower().second;

    //*
    if (clusterprobcut)
    {
      if (prob < clusterprob && cutson)
      {
        h_cutCounter->Fill(1);
        continue;
      }
    }
    else if (clus_chisq > clus_chisq_cut && cutson)
    {
      h_cutCounter->Fill(1);
      continue;
    }

    //*/
    TLorentzVector photon1;
    photon1.SetPtEtaPhiE(clus_pt, clus_eta, clus_phi, clusE);
    pi0smearvec[0] = SmearPhoton4vector(photon1, badcalibsmear);

    // auto& photons = (SPMC_bool) ? truth_photons : truth_pi0_photons;
    auto& photons = truth_pi0_photons;
    // std::cout << "truth_pi0_photons.size() = " << truth_pi0_photons.size() << std::endl;
    h_ndecayphotons->Fill(truth_pi0_photons.size());

    for (auto tr_phot : photons)
    {
      float delR = pi0smearvec[0].DeltaR(tr_phot);
      h_delR_recTrth->Fill(delR);
      float res = pi0smearvec[0].E() / tr_phot.E();
      float delPhi = pi0smearvec[0].Phi() - tr_phot.Phi();
      // if (delPhi > TMath::TwoPi()) delPhi -= TMath::TwoPi();
      // if (delPhi < -TMath::TwoPi()) delPhi += TMath::TwoPi();
      if (delR < 0.02)
      {
        if (debug) std::cout << "match clusE=" << pi0smearvec[0].E() << "  truthE=" << tr_phot.E() << " delPhi=" << delPhi << std::endl;
        h_matched_res->Fill(res, pi0smearvec[0].Eta());
        h_res_e->Fill(res, pi0smearvec[0].E());
        h_res_ptTr->Fill(tr_phot.Pt(), res);
        // h_res_e_eta->Fill(res, tr_phot.E(), lt_eta);
        // h_res_e_phi->Fill(res, tr_phot.E(), lt_phi);
        h_res->Fill(res);
        // h_delEta_e_eta->Fill(pi0smearvec[0].Eta() - tr_phot.Eta(), tr_phot.E(), lt_eta);
        // h_delPhi_e_eta->Fill(delPhi, tr_phot.E(), lt_eta);
        // h_delPhi_e_phi->Fill(delPhi, tr_phot.E(), lt_phi);
        h_truthE->Fill(tr_phot.E());
        h_clus_ELoss_2d->Fill(tr_phot.Pt(), tr_phot.E() - clusE);
        h_clus_ERatio_2d->Fill(tr_phot.Pt(), clusE / tr_phot.E());
        multimatchint++;
      }
    }

    TLorentzVector ph1_trEtaPhi;
    ph1_trEtaPhi.SetPtEtaPhiE(0, 0, 0, 0);
    for (auto tr_phot : photons)
    {
      float delR = pi0smearvec[0].DeltaR(tr_phot);
      // float res = pi0smearvec[0].E() / tr_phot.E();
      if (delR < 0.02)
      {
        ph1_trEtaPhi.SetPtEtaPhiE(tr_phot.E() / TMath::CosH(tr_phot.Eta()), tr_phot.Eta(), tr_phot.Phi(), tr_phot.E());
        if (debug) std::cout << "match  eta=" << ph1_trEtaPhi.Eta() << " E=" << ph1_trEtaPhi.E() << std::endl;
        match1 = true;
        break;
      }
    }

    if (Cluster_Debug)
    {
      if (filledClustersAfterCut1.insert(recoCluster).second)
      {
        h_cluster_etaphi_cuts[2]->Fill(clus_eta, clus_phi);
      }
    }

    /*
    if (eTCutbool)
    {
      if (pi0smearvec[0].Et() < etcut && cutson)
      {
        h_cutCounter->Fill(2);
        continue;
      }
    }
    else if ((pi0smearvec[0].Pt() < pt1ClusCut || pi0smearvec[0].Pt() > ptMaxCut) && cutson)
    {
      h_cutCounter->Fill(2);
      continue;
    }
    if (Cluster_Debug)
    {
      if (filledClustersAfterCut1.insert(recoCluster).second)
      {
        h_cluster_etaphi_cuts[2]->Fill(pi0smearvec[0].Eta(), pi0smearvec[0].Phi());
      }
    }
    */

    if (additionalsmearing)
    {
      if (Cluster_Debug) h_cluster_etaphi_cuts[3]->Fill(pi0smearvec[0].Eta(), pi0smearvec[0].Phi());

      /*
      if (clusterprobcut)
      {
        if (prob < clusterprob && cutson)
        {
          h_cutCounter->Fill(1);
          continue;
        }
      }
      else if (clus_chisq > clus_chisq_cut && cutson)
      {
        h_cutCounter->Fill(1);
        continue;
      }
      //*/

      if (eTCutbool)
      {
        if (pi0smearvec[0].Et() < etcut && cutson)
        {
          h_cutCounter->Fill(2);
          continue;
        }
      }
      else if ((pi0smearvec[0].Pt() < pt1ClusCut || pi0smearvec[0].Pt() > ptMaxCut) && cutson)
      {
        h_cutCounter->Fill(2);
        continue;
      }

      if (Cluster_Debug)
      {
        if (filledClustersAfterCut2.insert(recoCluster).second)
        {
          h_cluster_etaphi_cuts[4]->Fill(pi0smearvec[0].Eta(), pi0smearvec[0].Phi());
        }
      }
    }
    else if (!additionalsmearing)
    {
      if (Cluster_Debug) h_cluster_etaphi_cuts[3]->Fill(clus_eta, clus_phi);

      if (eTCutbool)
      {
        if (photon1.Et() < etcut && cutson)
        {
          h_cutCounter->Fill(2);
          continue;
        }
      }
      else if ((photon1.Pt() < pt1ClusCut || photon1.Pt() > ptMaxCut) && cutson)
      {
        h_cutCounter->Fill(2);
        continue;
      }

      if (Cluster_Debug)
      {
        if (filledClustersAfterCut2.insert(recoCluster).second)
        {
          h_cluster_etaphi_cuts[4]->Fill(clus_eta, clus_phi);
        }
      }
    }
    h_clusE->Fill(clusE);
    // if (clusE < 0.2) continue;
    //  std::cout << "clusE = " << clusE <<  " clus_eta = " << clus_eta <<  " clus_phi = " << clus_phi <<  " clus_pt = " << clus_pt <<  " clus_chisq = " << clus_chisq << std::endl;

    // loop over the towers in the cluster
    RawCluster::TowerConstRange towerCR = recoCluster->get_towers();
    RawCluster::TowerConstIterator toweriter;
    // lt_-> lead tower related
    float lt_e = -1000;
    unsigned int lt_eta = -1;
    bool hotClus = false;
    for (toweriter = towerCR.first; toweriter != towerCR.second; ++toweriter)
    {
      int towereta = m_geometry->get_tower_geometry(toweriter->first)->get_bineta();
      int towerphi = m_geometry->get_tower_geometry(toweriter->first)->get_binphi();
      unsigned int key = TowerInfoDefs::encode_emcal(towereta, towerphi);
      unsigned int channel = towers->decode_key(key);
      float energy = towers->get_tower_at_channel(channel)->get_energy();
      if (energy > lt_e)
      {
        lt_e = energy;
        lt_eta = towereta;
      }
    }

    if (lt_eta > 95)
    {
      h_cutCounter->Fill(3);
      continue;
    }
    if (Cluster_Debug && filledClustersAfterCut3.insert(recoCluster).second)
    {
      if (additionalsmearing)
        h_cluster_etaphi_cuts[5]->Fill(pi0smearvec[0].Eta(), pi0smearvec[0].Phi());
      else if (!additionalsmearing)
        h_cluster_etaphi_cuts[5]->Fill(clus_eta, clus_phi);
    }

    if (etabyeta) h_pt_eta[lt_eta]->Fill(clus_pt);

    if (dynMaskClus && hotClus == true && cutson)
    {
      h_cutCounter->Fill(4);
      continue;
    }
    if (Cluster_Debug && filledClustersAfterCut4.insert(recoCluster).second)
    {
      if (additionalsmearing)
        h_cluster_etaphi_cuts[6]->Fill(pi0smearvec[0].Eta(), pi0smearvec[0].Phi());
      else if (!additionalsmearing)
        h_cluster_etaphi_cuts[6]->Fill(clus_eta, clus_phi);
    }
    if (debug) std::cout << " " << "Cluster Loop: 2 " << std::endl;
    // for (clusterIter2 = clusterEnd.first; clusterIter2 != clusterEnd.second; clusterIter2++)
    for (clusterIter2 = std::next(clusterIter); clusterIter2 != clusterEnd.second; ++clusterIter2)
    {
      if (clusterIter == clusterIter2)
      {
        h_cutCounter->Fill(5);
        continue;
      }
      if (Cluster_Debug && filledClustersAfterCut5.insert(recoCluster).second)
      {
        if (additionalsmearing)
          h_cluster_etaphi_cuts[7]->Fill(pi0smearvec[0].Eta(), pi0smearvec[0].Phi());
        else if (!additionalsmearing)
          h_cluster_etaphi_cuts[7]->Fill(clus_eta, clus_phi);
      }

      RawCluster* recoCluster2 = clusterIter2->second;

      CLHEP::Hep3Vector E_vec_cluster2;
      if (pp_rawcluster)
      {
        E_vec_cluster2 = RawClusterUtility::GetEVec(*recoCluster2, vertex);
      }
      else if (!pp_rawcluster)  // i.e if AuAu
      {
        E_vec_cluster2 = RawClusterUtility::GetECoreVec(*recoCluster2, vertex);
      }

      float clus2E = E_vec_cluster2.mag();
      float clus2_eta = E_vec_cluster2.pseudoRapidity();
      float clus2_phi = E_vec_cluster2.phi();
      float clus2_pt = E_vec_cluster2.perp();
      float clus2_chisq = recoCluster2->get_chi2();
      float prob2 = recoCluster2->get_prob();
      if (clusterprobcut)
      {
        if (prob2 < clusterprob && cutson)
        {
          h_cutCounter->Fill(6);
          continue;
        }
      }
      else if (clus2_chisq > clus_chisq_cut && cutson)
      {
        h_cutCounter->Fill(6);
        continue;
      }
      if (Cluster_Debug && filledClustersAfterCut6.insert(recoCluster).second)
      {
        if (additionalsmearing)
          h_cluster_etaphi_cuts[8]->Fill(pi0smearvec[0].Eta(), pi0smearvec[0].Phi());
        else if (!additionalsmearing)
          h_cluster_etaphi_cuts[8]->Fill(clus_eta, clus_phi);
      }
      TLorentzVector photon2;
      photon2.SetPtEtaPhiE(clus2_pt, clus2_eta, clus2_phi, clus2E);
      pi0smearvec[1] = SmearPhoton4vector(photon2, badcalibsmear);

      TLorentzVector ph2_trEtaPhi;
      ph2_trEtaPhi.SetPtEtaPhiE(0, 0, 0, 0);
      for (auto tr_phot : photons)
      {
        float delR = pi0smearvec[1].DeltaR(tr_phot);
        // float res = pi0smearvec[1].E() / tr_phot.E();
        if (delR < 0.02)
        {
          ph2_trEtaPhi.SetPtEtaPhiE(tr_phot.E() / TMath::CosH(tr_phot.Eta()), tr_phot.Eta(), tr_phot.Phi(), tr_phot.E());
          if (debug) std::cout << "match  eta=" << ph2_trEtaPhi.Eta() << " E=" << ph2_trEtaPhi.E() << std::endl;
          if (match1) match2 = true;
        }
      }

      TLorentzVector pi0;
      TLorentzVector pi0_trKin = ph1_trEtaPhi + ph2_trEtaPhi;
      pi0smearvec[2] = pi0smearvec[0] + pi0smearvec[1];

      // set pi0 to the smeared version if additional smearing is added
      if (additionalsmearing)
      {
        pi0 = pi0smearvec[2];
      }
      else if (!additionalsmearing)
      {
        pi0 = photon1 + photon2;
      }

      h_reco_etaphi_cuts[1]->Fill(pi0smearvec[2].Eta(), pi0smearvec[2].Phi());

      if (additionalsmearing)
      {
        if (eTCutbool)
        {
          if (pi0smearvec[1].Et() < etcut && cutson)
          {
            h_cutCounter->Fill(7);
            continue;
          }
        }
        else if ((pi0smearvec[1].Pt() < pt2ClusCut || pi0smearvec[1].Pt() > ptMaxCut) && cutson)
        {
          h_cutCounter->Fill(7);
          continue;
        }

        if (Cluster_Debug && filledClustersAfterCut7.insert(recoCluster).second)
        {
          h_reco_etaphi_cuts[2]->Fill(pi0smearvec[2].Eta(), pi0smearvec[2].Phi());
          h_cluster_etaphi_cuts[9]->Fill(pi0smearvec[0].Eta(), pi0smearvec[0].Phi());
          if (filledClustersAfterCut7.insert(recoCluster2).second || Cluster_Debug2) h_cluster_etaphi_cuts[9]->Fill(pi0smearvec[1].Eta(), pi0smearvec[1].Phi());
        }

        if (fabs(pi0smearvec[0].E() - pi0smearvec[1].E()) / (pi0smearvec[0].E() + pi0smearvec[1].E()) > maxAlpha && cutson)
        {
          h_cutCounter->Fill(8);
          continue;
        }

        if (Cluster_Debug && filledClustersAfterCut8.insert(recoCluster).second)
        {
          h_reco_etaphi_cuts[3]->Fill(pi0smearvec[2].Eta(), pi0smearvec[2].Phi());
          h_cluster_etaphi_cuts[10]->Fill(pi0smearvec[0].Eta(), pi0smearvec[0].Phi());
          if (filledClustersAfterCut8.insert(recoCluster2).second || Cluster_Debug2) h_cluster_etaphi_cuts[10]->Fill(pi0smearvec[1].Eta(), pi0smearvec[1].Phi());
        }

        if (pi0smearvec[0].DeltaR(pi0smearvec[1]) > maxDr && cutson)
        {
          h_cutCounter->Fill(9);
          continue;
        }

        if (Cluster_Debug && filledClustersAfterCut9.insert(recoCluster).second)
        {
          h_reco_etaphi_cuts[4]->Fill(pi0smearvec[2].Eta(), pi0smearvec[2].Phi());
          h_cluster_etaphi_cuts[11]->Fill(pi0smearvec[0].Eta(), pi0smearvec[0].Phi());
          if (filledClustersAfterCut9.insert(recoCluster2).second || Cluster_Debug2) h_cluster_etaphi_cuts[11]->Fill(pi0smearvec[1].Eta(), pi0smearvec[1].Phi());
        }

        if (pi0smearvec[2].Pt() < pi0ptcut)
        {
          h_cutCounter->Fill(10);
          continue;
        }

        if (Cluster_Debug && filledClustersAfterCut10.insert(recoCluster).second)
        {
          h_reco_etaphi_cuts[5]->Fill(pi0smearvec[2].Eta(), pi0smearvec[2].Phi());
          h_cluster_etaphi_cuts[12]->Fill(pi0smearvec[0].Eta(), pi0smearvec[0].Phi());
          if (filledClustersAfterCut10.insert(recoCluster2).second || Cluster_Debug2) h_cluster_etaphi_cuts[12]->Fill(pi0smearvec[1].Eta(), pi0smearvec[1].Phi());
        }

        if (etaCutbool && abs(pi0smearvec[2].Eta()) > etacutval)
        {
          h_cutCounter->Fill(13);
          continue;
        }

        if (Cluster_Debug && filledClustersAfterCut11.insert(recoCluster).second)
        {
          h_reco_etaphi_cuts[6]->Fill(pi0smearvec[2].Eta(), pi0smearvec[2].Phi());
          h_cluster_etaphi_cuts[13]->Fill(pi0smearvec[0].Eta(), pi0smearvec[0].Phi());
          if (filledClustersAfterCut11.insert(recoCluster2).second || Cluster_Debug2) h_cluster_etaphi_cuts[13]->Fill(pi0smearvec[1].Eta(), pi0smearvec[1].Phi());
        }
      }
      else if (!additionalsmearing)
      {
        if (eTCutbool)
        {
          if (photon2.Et() < etcut && cutson)
          {
            h_cutCounter->Fill(7);
            continue;
          }
        }
        else if ((photon2.Pt() < pt2ClusCut || photon2.Pt() > ptMaxCut) && cutson)
        {
          h_cutCounter->Fill(7);
          continue;
        }

        if (Cluster_Debug && filledClustersAfterCut7.insert(recoCluster).second)
        {
          h_reco_etaphi_cuts[2]->Fill(pi0.Eta(), pi0.Phi());
          h_cluster_etaphi_cuts[9]->Fill(clus_eta, clus_phi);
          if (filledClustersAfterCut7.insert(recoCluster2).second || Cluster_Debug2) h_cluster_etaphi_cuts[9]->Fill(clus2_eta, clus2_phi);
        }

        if (fabs(photon1.E() - photon2.E()) / (photon1.E() + photon2.E()) > maxAlpha && cutson)
        {
          h_cutCounter->Fill(8);
          continue;
        }

        if (Cluster_Debug && filledClustersAfterCut8.insert(recoCluster).second)
        {
          h_reco_etaphi_cuts[3]->Fill(pi0.Eta(), pi0.Phi());
          h_cluster_etaphi_cuts[10]->Fill(clus_eta, clus_phi);
          if (filledClustersAfterCut8.insert(recoCluster2).second || Cluster_Debug2) h_cluster_etaphi_cuts[10]->Fill(clus2_eta, clus2_phi);
        }

        if (photon1.DeltaR(photon2) > maxDr && cutson)
        {
          h_cutCounter->Fill(9);
          continue;
        }

        if (Cluster_Debug && filledClustersAfterCut9.insert(recoCluster).second)
        {
          h_reco_etaphi_cuts[4]->Fill(pi0.Eta(), pi0.Phi());
          h_cluster_etaphi_cuts[11]->Fill(clus_eta, clus_phi);
          if (filledClustersAfterCut9.insert(recoCluster2).second || Cluster_Debug2) h_cluster_etaphi_cuts[11]->Fill(clus2_eta, clus2_phi);
        }

        if (pi0.Pt() < pi0ptcut)
        {
          h_cutCounter->Fill(10);
          continue;
        }

        if (Cluster_Debug && filledClustersAfterCut10.insert(recoCluster).second)
        {
          h_reco_etaphi_cuts[5]->Fill(pi0.Eta(), pi0.Phi());
          h_cluster_etaphi_cuts[12]->Fill(clus_eta, clus_phi);
          if (filledClustersAfterCut10.insert(recoCluster2).second || Cluster_Debug2) h_cluster_etaphi_cuts[12]->Fill(clus2_eta, clus2_phi);
        }
        if (etaCutbool && abs(pi0.Eta()) > etacutval)
        {
          h_cutCounter->Fill(13);
          continue;
        }

        if (Cluster_Debug && filledClustersAfterCut11.insert(recoCluster).second)
        {
          h_reco_etaphi_cuts[6]->Fill(pi0.Eta(), pi0.Phi());
          h_cluster_etaphi_cuts[13]->Fill(clus_eta, clus_phi);
          if (filledClustersAfterCut11.insert(recoCluster2).second || Cluster_Debug2) h_cluster_etaphi_cuts[13]->Fill(clus2_eta, clus2_phi);
        }
      }

      // loop over the towers in the cluster
      RawCluster::TowerConstRange towerCR2 = recoCluster2->get_towers();
      RawCluster::TowerConstIterator toweriter2;
      bool hotClus2 = false;
      for (toweriter2 = towerCR2.first; toweriter2 != towerCR2.second; ++toweriter2)
      {
        int towereta = m_geometry->get_tower_geometry(toweriter2->first)->get_bineta();
        int towerphi = m_geometry->get_tower_geometry(toweriter2->first)->get_binphi();

        for (size_t i = 0; i < ht_eta.size(); i++)
        {
          if (towerphi == ht_phi[i] && towereta == ht_phi[i]) hotClus2 = true;
        }
      }
      h_etaphi_clus->Fill(clus_phi, clus_eta);

      if (dynMaskClus && hotClus2 == true && cutson)
      {
        h_cutCounter->Fill(11);
        continue;
      }

      if (Cluster_Debug && filledClustersAfterCut12.insert(recoCluster).second)
      {
        if (additionalsmearing)
          h_reco_etaphi_cuts[7]->Fill(pi0smearvec[2].Eta(), pi0smearvec[2].Phi());
        else if (!additionalsmearing)
          h_reco_etaphi_cuts[7]->Fill(pi0.Eta(), pi0.Phi());
      }

      h_reco_etaphi->Fill(pi0.Phi(), pi0.Eta());  // pi0 is the same as the smeared version if adding smearing

      /////////////////////////////////////////////////
      //// Truth info
      // float weight = 1;
      PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
      // vector<TLorentzVector> truth_pions;
      if (debug) std::cout << " " << "truth: Primary Loop " << std::endl;
      if (truthinfo)
      {
        PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
        for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
        {
          // Get truth particle
          // PHG4Particle* particle = truthinfo->GetParticle(1);  // primary for the SPMC
          const PHG4Particle* truth = iter->second;
          if (!truthinfo->is_primary(truth)) continue;
          TLorentzVector myVector;
          // LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> myVector;
          if (!eta_weight)
          {
            myVector.SetXYZM(truth->get_px(), truth->get_py(), truth->get_pz(), 0.13497);
            // myVector.SetPxPyPzE(truth->get_px(), truth->get_py(), truth->get_pz(), truth->get_e());
          }
          else if (eta_weight)
          {
            myVector.SetXYZM(truth->get_px(), truth->get_py(), truth->get_pz(), 0.54786);
            // myVector.SetPxPyPzE(truth->get_px(), truth->get_py(), truth->get_pz(), truth->get_e());
          }
          float energy = myVector.E();
          // weight = myVector.Pt() * TMath::Exp(-3 * myVector.Pt());

          int id = truth->get_pid();
          h_truth_pid_p->Fill(id);

          double Pt = myVector.Pt();
          double weight_function;
          if (Pythia_weight)
          {
            double t = 4.5229;
            double w = 0.0912;
            double A = 528.197;
            double B = 141.505;
            double n = 9.86803;
            double m_param = 7.48221;
            double p0 = 0.592415;
            //"((1 / (1 + exp((x - [0]) / [1]))) * [2] / pow(1 + x / [3], [4]) + (1 - (1 / (1 + exp((x - [0]) / [1])))) * [5] / (pow(x, [6])))"
            weight_function = ((1 / (1 + exp((Pt - t) / w))) * A / pow(1 + Pt / p0, m_param) + (1 - (1 / (1 + exp((Pt - t) / w)))) * B / (pow(Pt, n)));
          }
          else
          {
            //--------------------Alternative paramaterization, woods saxon + hagedorn + power law
            double t = 4.5;
            double w = 0.114;
            double A = 229.6;
            double B = 14.43;
            double n = 8.1028;
            double m_param = 10.654;
            double p0 = 1.466;
            weight_function = ((1 / (1 + exp((Pt - t) / w))) * A / pow(1 + Pt / p0, m_param) + (1 - (1 / (1 + exp((Pt - t) / w)))) * B / (pow(Pt, n)));

            if (eta_weight)
            {
              inv_yield *= 0.5 * pow((1.2 + sqrt(pow(0.54786, 2) + pow(Pt, 2))) / (1.2 + sqrt(pow(0.1349768, 2) + pow(Pt, 2))), -10);  // mT scaling
              weight_function *= 0.5 * pow((1.2 + sqrt(pow(0.54786, 2) + pow(Pt, 2))) / (1.2 + sqrt(pow(0.1349768, 2) + pow(Pt, 2))), -10);
            }
          }

          inv_yield = WeightScale * Pt * weight_function;  //
          // h_pion_pt_weight->Fill(pi0.Pt(), inv_yield);

          /*
          if (SPMC_bool && inv_yield != 0)
          {
            float spectrum_value = getSPMCpTspectrum(static_cast<float>(Pt));
            if (spectrum_value != 0)
            {
              inv_yield = inv_yield / spectrum_value;
            }
            else
            {
              std::cerr << "Warning: Spectrum value is zero at Pt = " << Pt << std::endl;
              inv_yield = 0;
            }
          }
          */

          h_inv_yield->Fill(Pt, inv_yield);
          h_yield->Fill(Pt, weight_function);
          h_InvMass_weighted->Fill(pi0.M(), inv_yield);
          h_InvMass_smear_weighted->Fill(pi0smearvec[2].M(), inv_yield);
          h_InvMass_smear_weighted_2d->Fill(pi0smearvec[2].Pt(), pi0smearvec[2].M(), inv_yield);
          // h_InvMass_photonE_smear_weighted_3d->Fill(pi0smearvec[0].Pt(), pi0smearvec[1].Pt(), pi0smearvec[2].M(), inv_yield);
          // h_InvMass_smear_weighted_asymmetry_3d->Fill(pi0smearvec[2].Pt(), pi0smearvec[2].M(), fabs(pi0smearvec[0].E() - pi0smearvec[1].E()) / (pi0smearvec[0].E() + pi0smearvec[1].E()), inv_yield);
          h_truth_e->Fill(energy, inv_yield);
          h_truth_eta->Fill(myVector.Eta(), inv_yield);
          h_truth_pt->Fill(myVector.Pt(), inv_yield);

          h_InvMass_smear_weighted_eta_3d->Fill(pi0smearvec[2].Pt(), pi0smearvec[2].M(), pi0smearvec[2].Eta(), inv_yield);
          h_InvMass_smear_weighted_eta_2d->Fill(pi0smearvec[2].Eta(), pi0smearvec[2].M(), inv_yield);
          h_reco_etaphi_cuts[8]->Fill(pi0.Eta(), pi0.Phi(), inv_yield);

          if (filltruthspectrum)
          {
            float delR = pi0smearvec[2].DeltaR(myVector);
            if ((id == 111 || (eta_weight && id == 221)) && delR < 0.015)
            {
              h_truth_spectrum1->Fill(myVector.Pt());
              h_truthmatched_mass1_2d->Fill(pi0.Pt(), pi0.M(), inv_yield);
              h_truthmatched_mass1->Fill(pi0.M(), inv_yield);
            }
            else if ((id == 111 || (eta_weight && id == 221)) && delR < 0.1)
            {
              h_truthmatched_mass2_2d->Fill(pi0.Pt(), pi0.M(), inv_yield);
              h_truthmatched_mass2->Fill(pi0.M(), inv_yield);
            }
            else if ((id == 111 || (eta_weight && id == 221)) && delR < 0.2)
            {
              h_truthmatched_mass3_2d->Fill(pi0.Pt(), pi0.M(), inv_yield);
              h_truthmatched_mass3->Fill(pi0.M(), inv_yield);
            }
            // photon truthmatching
          }

          if (debug) std::cout << "truth pt=" << Pt << "   weight function=" << weight_function << "  inv_yield=" << inv_yield << std::endl;
          if (debug) std::cout << "M=" << myVector.M() << "   E=" << energy << "  pt=" << myVector.Pt() << "  eta=" << myVector.Eta() << std::endl;
        }

        // truth secondary loops
        if (debug) std::cout << " " << "truth: Secondary Loop " << std::endl;
        PHG4TruthInfoContainer::Range second_range = truthinfo->GetSecondaryParticleRange();
        float m_g4 = 0;
        for (PHG4TruthInfoContainer::ConstIterator siter = second_range.first; siter != second_range.second; ++siter)
        {
          if (m_g4 >= 19999) break;
          // Get photons from pi0 decays
          const PHG4Particle* truth = siter->second;
          int id = truth->get_pid();
          h_truth_pid_s->Fill(id);
          //*/
        }
      }

      //*
      //*/
      if (debug) std::cout << " " << "truth: Loops Done " << std::endl;
      if (additionalsmearing)
      {
        h_InvMass_smear->Fill(pi0smearvec[2].M());
        h_InvMass_smear_2d->Fill(pi0smearvec[2].Pt(), pi0smearvec[2].M());
        h_InvMass_smear_eta_3d->Fill(pi0smearvec[2].Pt(), pi0smearvec[2].M(), pi0smearvec[2].Eta());
        h_InvMass_smear_eta_2d->Fill(pi0smearvec[2].Eta(), pi0smearvec[2].M());
      }
      h_InvMass_2d->Fill(pi0.Pt(), pi0.M());
      h_pion_pt->Fill(pi0.Pt());
      h_InvMass->Fill(pi0.M());
      h_reco_ALLphotonE->Fill(photon1.E());
      h_reco_ALLphotonE->Fill(photon2.E());
      h_reco_ALLphotonE_2d->Fill(photon1.Pt(), photon1.E());
      h_reco_ALLphotonE_2d->Fill(photon2.Pt(), photon2.E());
      for (auto tr_phot : photons)
      {
        float delR = pi0smearvec[0].DeltaR(tr_phot);
        if (delR < 0.02)
        {
          h_clus_ELoss_weighted_2d->Fill(tr_phot.Pt(), tr_phot.E() - clusE);
          h_clus_ERatio_weighted_2d->Fill(tr_phot.Pt(), clusE / tr_phot.E());
        }
      }
      if (match2 && pi0_trKin.M() > 0.001)
      {
        h_delR_pionrecTrth->Fill(pi0smearvec[2].DeltaR(pi0_trKin));
        h_truthmatched_mass->Fill(pi0smearvec[2].M());
        h_truthmatched_mass_2d->Fill(pi0smearvec[2].Pt(), pi0smearvec[2].M());
        h_truthmatched_mass_weighted->Fill(pi0smearvec[2].M(), inv_yield);
        h_truthmatched_mass_weighted_2d->Fill(pi0smearvec[2].Pt(), pi0smearvec[2].M(), inv_yield);
        h_truthmatched_mass_eta_2d->Fill(pi0smearvec[2].Eta(), pi0smearvec[2].M());
        h_reco_etaphi_cuts[9]->Fill(pi0smearvec[2].Eta(), pi0smearvec[2].Phi());
        h_reco_etaphi_cuts[10]->Fill(pi0smearvec[2].Eta(), pi0smearvec[2].Phi(), inv_yield);
        h_pi0_ELoss_2d->Fill(pi0_trKin.Pt(), pi0_trKin.E() - pi0smearvec[2].E());
        if (pi0_trKin.E() != 0) h_pi0_ERatio_2d->Fill(pi0_trKin.Pt(), pi0smearvec[2].E() / pi0_trKin.E());
        // std::cout << pi0_trKin.M() << std::endl;//
        if (eta_weight)  //&&pi0_trKin.M() >= 0.4 && pi0_trKin.M() <= 0.8
        {
          h_truthmatched_mass_etameson_weighted->Fill(pi0smearvec[2].M(), inv_yield);
          h_truthmatched_mass_etameson_weighted_2d->Fill(pi0smearvec[2].Pt(), pi0smearvec[2].M(), inv_yield);
          h_reco_etaphi_cuts[11]->Fill(pi0smearvec[2].Eta(), pi0smearvec[2].Phi(), inv_yield);
          h_truthmatched_mass_etameson_weighted_eta_3d->Fill(pi0smearvec[2].Pt(), pi0smearvec[2].M(), pi0smearvec[2].Eta(), inv_yield);
          h_eta_ELoss_2d->Fill(pi0_trKin.Pt(), pi0_trKin.E() - pi0smearvec[2].E());
          if (pi0_trKin.E() != 0) h_eta_ERatio_2d->Fill(pi0_trKin.Pt(), pi0smearvec[2].E() / pi0_trKin.E());
        }
      }
    }  // clusterIter2

    h_clusmultimatch->Fill(multimatchint);
    multimatchint = 0;
  }  // clusteriter1 loop

  return Fun4AllReturnCodes::EVENT_OK;
  //}
}

int CaloAna::End(PHCompositeNode* /*topNode*/)
{
  outfile->cd();

  outfile->Write();
  outfile->Close();
  delete outfile;
  hm->dumpHistos(outfilename, "UPDATE");
  std::cout << "funkycounter: " << funkyCaloStuffcounter << ", VertexMap Failure:" << VertexMapFailcounter << std::endl;

  return 0;
}

float CaloAna::getWeight(int ieta, float pt)
{
  if (ieta < 0 || ieta > 95) return 0;
  float val = h_pt_rw[ieta]->GetBinContent(h_pt_rw[ieta]->FindBin(pt));
  if (val == 0) return 0;
  return 1 / val;
}

float CaloAna::getSPMCpTspectrum(float pt)
{
  if (h_sp_pt_rw == nullptr)
  {
    std::cerr << "Error: h_sp_pt_rw is nullptr in getSPMCpTspectrum!" << std::endl;
    return 0.0f;
  }

  int bin = h_sp_pt_rw->FindBin(pt);
  float val = h_sp_pt_rw->GetBinContent(bin);
  if (debug) std::cout << "Histogram value at pt=" << pt << " is " << val << std::endl;
  return val;
}

TF1* CaloAna::fitHistogram(TH1* h)
{
  TF1* fitFunc = new TF1("fitFunc", "[0]/[2]/2.5*exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x + [5]*x^2 + [6]*x^3", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());

  fitFunc->SetParameter(0, 5);
  fitFunc->SetParameter(1, target_pi0_mass);
  fitFunc->SetParameter(2, 0.01);
  fitFunc->SetParameter(3, 0.0);
  fitFunc->SetParameter(4, 0.0);
  fitFunc->SetParameter(5, 0.0);
  fitFunc->SetParameter(6, 100);

  fitFunc->SetParLimits(0, 0, 10);
  fitFunc->SetParLimits(1, 0.113, 0.25);
  fitFunc->SetParLimits(2, 0.01, 0.04);
  fitFunc->SetParLimits(3, -2, 1);
  fitFunc->SetParLimits(4, 0, 40);
  fitFunc->SetParLimits(5, -150, 50);
  fitFunc->SetParLimits(6, 0, 200);

  fitFunc->SetRange(0.05, 0.7);

  // Perform the fit
  h->Fit("fitFunc", "QN");

  return fitFunc;
}

void CaloAna::fitEtaSlices(const std::string& infile, const std::string& fitOutFile, const std::string& cdbFile)
{
  TFile* fin = new TFile(infile.c_str());
  std::cout << "getting hists" << std::endl;
  TH1F* h_peak_eta = new TH1F("h_peak_eta", "", 96, 0, 96);
  TH1F* h_sigma_eta = new TH1F("h_sigma_eta", "", 96, 0, 96);
  TH1F* h_p3_eta = new TH1F("h_p3_eta", "", 96, 0, 96);
  TH1F* h_p4_eta = new TH1F("h_p4_eta", "", 96, 0, 96);
  TH1F* h_p5_eta = new TH1F("h_p5_eta", "", 96, 0, 96);
  TH1F* h_p6_eta = new TH1F("h_p6_eta", "", 96, 0, 96);
  TH1F* h_p0_eta = new TH1F("h_p0_eta", "", 96, 0, 96);
  if (!fin)
  {
    std::cout << "pi0EtaByEta::fitEtaSlices null fin" << std::endl;
    exit(1);
  }
  TH1F* h_M_eta[96];

  if (etabyeta)
  {
    for (int i = 0; i < 96; i++)
    {
      h_M_eta[i] = (TH1F*) fin->Get(Form("h_mass_eta_lt_rw%d", i));
      h_M_eta[i]->Scale(1. / h_M_eta[i]->Integral(), "width");
    }
  }

  TF1* fitFunOut[96];
  for (int i = 0; i < 96; i++)
  {
    if (!h_M_eta[i])
    {
      std::cout << "pi0EtaByEta::fitEtaSlices null hist" << std::endl;
    }

    fitFunOut[i] = fitHistogram(h_M_eta[i]);
    fitFunOut[i]->SetName(Form("f_pi0_eta%d", i));
    float mass_val_out = fitFunOut[i]->GetParameter(1);
    float mass_err_out = fitFunOut[i]->GetParError(1);
    h_peak_eta->SetBinContent(i + 1, mass_val_out);
    if (isnan(h_M_eta[i]->GetEntries()))
    {
      h_peak_eta->SetBinError(i + 1, 0);
      continue;
    }
    h_peak_eta->SetBinError(i + 1, mass_err_out);
    h_sigma_eta->SetBinContent(i + 1, fitFunOut[i]->GetParameter(2));
    h_sigma_eta->SetBinError(i + 1, fitFunOut[i]->GetParError(2));
    h_p3_eta->SetBinContent(i + 1, fitFunOut[i]->GetParameter(3));
    h_p3_eta->SetBinError(i + 1, fitFunOut[i]->GetParError(3));
    h_p4_eta->SetBinContent(i + 1, fitFunOut[i]->GetParameter(4));
    h_p4_eta->SetBinError(i + 1, fitFunOut[i]->GetParError(4));
    h_p5_eta->SetBinContent(i + 1, fitFunOut[i]->GetParameter(5));
    h_p5_eta->SetBinError(i + 1, fitFunOut[i]->GetParError(5));
    h_p6_eta->SetBinContent(i + 1, fitFunOut[i]->GetParameter(6));
    h_p6_eta->SetBinError(i + 1, fitFunOut[i]->GetParError(6));
    h_p0_eta->SetBinContent(i + 1, fitFunOut[i]->GetParameter(0));
    h_p0_eta->SetBinError(i + 1, fitFunOut[i]->GetParError(0));
  }

  CDBTTree* cdbttree1 = new CDBTTree(cdbFile.c_str());
  CDBTTree* cdbttree2 = new CDBTTree(cdbFile.c_str());

  std::string m_fieldname = "Femc_datadriven_qm1_correction";

  for (int i = 0; i < 96; i++)
  {
    for (int j = 0; j < 256; j++)
    {
      float correction = target_pi0_mass / h_peak_eta->GetBinContent(i + 1);
      unsigned int key = TowerInfoDefs::encode_emcal(i, j);
      float val1 = cdbttree1->GetFloatValue(key, m_fieldname);
      cdbttree2->SetFloatValue(key, m_fieldname, val1 * correction);
    }
  }

  cdbttree2->Commit();
  cdbttree2->WriteCDBTTree();
  delete cdbttree2;
  delete cdbttree1;

  TFile* fit_out = new TFile(fitOutFile.c_str(), "recreate");
  fit_out->cd();
  for (auto& i : h_M_eta)
  {
    i->Write();
    delete i;
  }
  for (auto& i : fitFunOut)
  {
    i->Write();
    delete i;
  }

  h_p3_eta->Write();
  h_p4_eta->Write();
  h_p5_eta->Write();
  h_p6_eta->Write();
  h_p0_eta->Write();
  h_sigma_eta->Write();
  h_peak_eta->Write();
  fin->Close();

  std::cout << "finish fitting suc" << std::endl;

  return;
}

// function to generate random numbers
double CaloAna::generateRandomNumber()
{
  std::normal_distribution<double> dist(0.0, 1.0);  // mean 0, sigma 1
  float rand = dist(rng);
  // cout << "test rnd gen code: " << rand << endl;
  return rand;
}

// function to smear photon 4 vectors

TLorentzVector CaloAna::SmearPhoton4vector(const TLorentzVector& sourcephoton, double smearfactor)
{
  double smear = generateRandomNumber() * smearfactor + 1;
  TLorentzVector smearedphoton = sourcephoton * smear;
  return smearedphoton;
}