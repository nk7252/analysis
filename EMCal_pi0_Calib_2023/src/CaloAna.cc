#include "CaloAna.h"

// Global vertex includes
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

// Fun4All includes
#include <Event/Event.h>
#include <Event/packet.h>
#include <cdbobjects/CDBTTree.h>  // for CDBTTree
#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/recoConsts.h>

// G4Cells includes
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>

// Calorimeter/Cluster includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfoContainerv3.h>
#include <calobase/TowerInfoDefs.h>
#include <calobase/TowerInfov1.h>
#include <calobase/TowerInfov2.h>

// MBD
#include <mbd/BbcGeom.h>
#include <mbd/MbdPmtContainerV1.h>
#include <mbd/MbdPmtHit.h>

// ROOT includes
#include <Math/Vector4D.h>  // for ROOT::Math::PtEtaPhiMVector
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH3.h>
#include <TH3F.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TTree.h>

// general includes
#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

/// HEPMC truth includes
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#pragma GCC diagnostic pop
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

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
    }
    else
    {
      frw = TFile::Open("/sphenix/user/nkumar/analysis/EMCal_pi0_Calib_2023/macros/spi0_spectrum.root", "READ");
      h_original = (TH1F*) frw->Get("spi0_pt_spectrum");
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
  h_vtxmap_fail = new TH1F("h_vtxmap_fail", "Vtxmap Fail", 2, 0, 2);
  h_zvtx = new TH1F("h_zvtx", "Zvtx", 1000, -500, 500);

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
  h_clusE = new TH1F("h_clusE", "", 100, 0, 20);
  h_emcal_e_eta = new TH1F("h_emcal_e_eta", "", 96, 0, 96);
  h_pt1 = new TH1F("h_pt1", "", 100, 0, 20);
  h_pt2 = new TH1F("h_pt2", "", 100, 0, 20);
  h_pion_pt = new TH1F("h_pion_pt", "", 100, 0, 20);
  h_pion_pt_weight = new TH1F("h_pion_pt_weight", "", 100, 0, 20);
  h_nclusters = new TH1F("h_nclusters", "", 1000, 0, 1000);

  // Truth histos
  h_truth_eta = new TH1F("h_truth_eta", "", 100, -1.2, 1.2);
  h_truth_e = new TH1F("h_truth_e", "", 100, 0, 20);
  h_truth_pt = new TH1F("h_truth_pt", "", 100, 0, 20);
  h_truth_spectrum1 = new TH1F("h_truth_spectrum1", "", 8 * 10, 0, 20);
  h_truth_spectrum2 = new TH1F("h_truth_spectrum2", "", 8 * 10, 0, 20);
  h_truth_etaspectrum = new TH1F("h_truth_etaspectrum", "", 8 * 10, 0, 20);
  h_truth_pid_p = new TH1F("h_truth_pid_p", "Primary particle PIDs", 1000, -500, 500);
  h_truth_pid_s = new TH1F("h_truth_pid_s", "Secondary particle PIDs", 1000, -500, 500);
  h_delR_recTrth = new TH1F("h_delR_recTrth", "", 1000, 0, 5);
  h_delR_pionrecTrth = new TH1F("h_delR_recTrth", "", 5000, 0, 5);
  h_matched_res = new TH2F("h_matched_res", "", 100, 0, 1.5, 20, -1, 1);
  // histograms to extract MC photon resolution
  h_truthmatched_photon1E = new TH1F("h_truthmatched_photon1E", "matched Photon 1 Energy", 8 * 10, 0, 20);
  h_truthmatched_photon2E = new TH1F("h_truthmatched_photon2E", "matchedPhoton 2 Energy", 8 * 10, 0, 20);
  h_truthmatched_AllphotonE = new TH1F("h_truthmatched_AllphotonE", "All Photon Energy", 8 * 10, 0, 20);
  h_truth_ALLphotonE = new TH1F("h_truth_ALLphotonE", "All Photon Energy", 8 * 10, 0, 20);
  h_truthmatched_photon1E_weighted = new TH1F("h_truthmatched_photon1E_weighted", "matched Photon 1 Energy, weighted", 8 * 10, 0, 20);
  h_truthmatched_photon2E_weighted = new TH1F("h_truthmatched_photon2E_weighted", "matchedPhoton 2 Energy, weighted", 8 * 10, 0, 20);
  h_truthmatched_AllphotonE_weighted = new TH1F("h_truthmatched_AllphotonE_weighted", "All Photon Energy, weighted", 8 * 10, 0, 20);
  h_truthmatched_Photon_delR = new TH1F("h_truthmatched_Photon_delR", "Photon delR", 10000, 0, 5);
  // h_truth_ALLphotonE_weighted = new TH1F("h_truth_ALLphotonE_weighted", "All Photon Energy, weighted", 8 * 10, 0, 20);
  // reco photons
  h_reco_photon1E = new TH1F("h_reco_photon1E", "Reco Photon 1 Energy", 8 * 10, 0, 20);
  h_reco_photon2E = new TH1F("h_reco_photon2E", "Reco Photon 2 Energy", 8 * 10, 0, 20);
  h_reco_photon1E_2d = new TH2F("h_reco_photon1E_2d", "pT vs Reco Photon 1 Energy", 8 * 10, 0, 20, 8 * 10, 0, 20);
  h_reco_photon2E_2d = new TH2F("h_reco_photon2E_2d", "pT vs Reco Photon 2 Energy", 8 * 10, 0, 20, 8 * 10, 0, 20);
  h_reco_ALLphotonE_2d = new TH2F("h_reco_ALLphotonE_2d", "pT vs All Reco Photon Energy", 8 * 10, 0, 20, 8 * 10, 0, 20);
  h_reco_ALLphotonE = new TH1F("h_reco_ALLphotonE", "All Reco Photon Energy", 8 * 10, 0, 20);
  // All truth photons(secondaries)
  h_truth_ALLphotonpt = new TH1F("h_truth_ALLphotonpt", "All Photon pt", 100, 0, 20);
  h_truth_ALLphotonp = new TH1F("h_truth_ALLphotonp", "All Photon p", 100, 0, 20);

  h_reco_photon1E_weighted = new TH1F("h_reco_photon1E_weighted", "Reco Photon 1 Energy, weighted", 8 * 10, 0, 20);
  h_reco_photon2E_weighted = new TH1F("h_reco_photon2E_weighted", "Reco Photon 2 Energy, weighted", 8 * 10, 0, 20);
  h_reco_ALLphotonE_weighted = new TH1F("h_reco_ALLphotonE_weighted", "All Reco Photon Energy, weighted", 8 * 10, 0, 20);

  // pT differential Inv Mass
  h_InvMass = new TH1F("h_InvMass", "Invariant Mass", 600, 0, 1.2);
  h_InvMass_2d = new TH2F("h_InvMass_2d", "pT vs Invariant Mass", 8 * 10, 0, 20, 600, 0, 1.2);
  h_InvMass_weighted = new TH1F("h_InvMass_weighted", "Invariant Mass, weighted WSHP", 600, 0, 1.2);

  h_inv_yield = new TH2F("h_inv_yield", "Invariant Yield distribution", 8 * 10, 0, 20, 100, 0, 20);
  h_yield = new TH2F("h_yield", "Yield distribution", 8 * 10, 0, 20, 100, 0, 70);
  h_truthmatched_mass1 = new TH1F("h_truthmatched_mass1", "Invariant Mass, truth matched(delR<0.015)", 600, 0, 1.2);
  h_truthmatched_mass2 = new TH1F("h_truthmatched_mass2", "Invariant Mass, truth matched(delR<0.1)", 600, 0, 1.2);
  h_truthmatched_mass3 = new TH1F("h_truthmatched_mass3", "Invariant Mass, truth matched(delR<0.2)", 600, 0, 1.2);
  h_truthmatched_mass1_2d = new TH2F("h_truthmatched_mass1_2d", "pT vs Invariant Mass, truth matched(delR<0.015)", 8 * 10, 0, 20, 600, 0, 1.2);
  h_truthmatched_mass2_2d = new TH2F("h_truthmatched_mass2_2d", "pT vs Invariant Mass, truth matched(delR<0.1)", 8 * 10, 0, 20, 600, 0, 1.2);
  h_truthmatched_mass3_2d = new TH2F("h_truthmatched_mass3_2d", "pT vs Invariant Mass, truth matched(delR<0.2)", 8 * 10, 0, 20, 600, 0, 1.2);
  h_FullTruth_eta = new TH1F("h_FullTruth_eta", "Full Truth eta", 100, -1.2, 1.2);
  h_FullTruth_e = new TH1F("h_FullTruth_e", "Full Truth e", 100, 0, 20);
  h_FullTruth_pt = new TH1F("h_FullTruth_pt", "Full Truth pt", 100, 0, 20);
  h_FullTruth_p = new TH1F("h_FullTruth_p", "Full Truth p", 100, 0, 20);

  // 3d histogram to check for corelation between photon/cluster energies and invariant mass.
  // h_InvMass_photonE_smear_weighted_3d = new TH3F(Form("h_InvMass_smear%d_weighted_photonE_3d",badcalibsmearint ), Form("Photon Energies vs Invariant Mass, smear, weighted: %f percent",badcalibsmearint / 10.0f), 40, 0, 20, 40, 0, 20, 60, 0, 1.2);
  // 3d histogram to check for for corelation between pt invariant mass and asymmetry
  // h_InvMass_smear_weighted_asymmetry_3d = new TH3F(Form("h_InvMass_smear%d_weighted_asymmetry_3d", badcalibsmearint ),Form("pT vs Invariant Mass vs asymmetry + smear, weighted: %f percent", badcalibsmearint / 10.0f), 8 * 10, 0, 20, 60, 0, 1.2, 10, 0, 1);

  // 3d histogram to check for corelation between eta, pt and invariant mass
  h_InvMass_smear_eta_3d = new TH3F(Form("h_InvMass_smear%d_eta_3d", badcalibsmearint), Form("pT vs Invariant Mass vs eta + smear: %f percent", badcalibsmearint / 10.0f), 8 * 10, 0, 20, 60, 0, 1.2, 24, -1.2, 1.2);
  h_InvMass_smear_weighted_eta_3d = new TH3F(Form("h_InvMass_smear%d_weighted_eta_3d", badcalibsmearint), Form("pT vs Invariant Mass vs eta + smear, weighted: %f percent", badcalibsmearint / 10.0f), 8 * 10, 0, 20, 60, 0, 1.2, 24, -1.2, 1.2);
  // 2d histogram to check for corelation between eta, and invariant mass
  h_InvMass_smear_eta_2d = new TH2F(Form("h_InvMass_smear%d_eta_2d", badcalibsmearint), Form("eta vs Invariant Mass+ smear: %f percent", badcalibsmearint / 10.0f), 24, -1.2, 1.2, 120, 0, 1.2);
  h_InvMass_smear_weighted_eta_2d = new TH2F(Form("h_InvMass_smear%d_weighted_eta_2d", badcalibsmearint), Form("eta vs Invariant Mass+ smear, weighted: %f percent", badcalibsmearint / 10.0f), 24, -1.2, 1.2, 120, 0, 1.2);

  // for (int i = 0; i < 96; i++) h_pt_rw[i] = (TH1F*) frw->Get(Form("h_pt_eta%d", i));
  // these histograms below were not working when using SPMC bool
  rnd = new TRandom3();
  // smearing added SPMC
  badcalibsmear = static_cast<float>(badcalibsmearint) / 1000.0f;

  h_InvMass_smear = new TH1F(Form("h_InvMass_smear_%d", badcalibsmearint), Form("Invariant Mass + const smear: %f percent", badcalibsmearint / 10.0f), 600, 0, 1.2);

  h_InvMass_smear_2d = new TH2F(Form("h_InvMass_smear_2d_%d", badcalibsmearint), Form("pT vs Invariant Mass + const smear: %f percent", badcalibsmearint / 10.0f), 8 * 10, 0, 20, 600, 0, 1.2);

  // weighted variants
  h_InvMass_smear_weighted = new TH1F(Form("h_InvMass_smear_weighted_%d", badcalibsmearint), Form("Invariant Mass + const smear, weighted: %f percent", badcalibsmearint / 10.0f), 600, 0, 1.2);

  h_InvMass_smear_weighted_2d = new TH2F(Form("h_InvMass_smear_weighted_2d_%d", badcalibsmearint), Form("pT vs Invariant Mass + const smear, weighted: %f percent", badcalibsmearint / 10.0f), 8 * 10, 0, 20, 600, 0, 1.2);

  // towers node selection
  calotowerinfostring = (clust_waveform == true) ? "WAVEFORM_CEMC" : "TOWERINFO_CALIB_CEMC";
  // clustercontainer node selection
  clustcontainerstring = (poscor == true) ? "CLUSTER_POS_COR_CEMC" : ((clust_waveform == true) ? "CLUSTERINFO_CEMC" : "CLUSTER_CEMC");

  funkyCaloStuffcounter = 0;
  VertexMapFailcounter = 0;
  if (additionalsmearing == false) std::cout << "additional smearing is not being added" << std::endl;
  if (additionalsmearing == true) std::cout << "additional smearing is being added" << std::endl;
  // return 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloAna::process_event(PHCompositeNode* topNode)
{
  _eventcounter++;
  if (debug) std::cout << " " << "Process Event: Start  " << std::endl;
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
  float maxDr = 1.1;
  float maxAlpha = 0.5;      // asymmetry cut
  float clus_chisq_cut = 4;  // normally 4
  float clusterprob = 0.1;   // replacing chisqr cut
  float nClus_ptCut = 0.0;   // 0.5 normally
  float pi0ptcutfactor = 0;  // seto to 0 to effectively disable it
  float ptMaxCut = 100;      // no cut in data, as far as I know. so I set it to a value it is unlikely to reach
  float pt1ClusCut = 1.3;    // centrality dependence cuts 2.2 for both // 1.3
  float pt2ClusCut = 0.7;    // 0.7
  float etcut = 1.0;         // cluster ET cut
  float etacutval = 0.6;     // cluster pseudo-rapidity cut
  float zvtx_cut_val = 30;   // z vertex cut value

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
  if (getVtx)
  {
    GlobalVertexMap* vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
    if (!vertexmap)
    {
      // if (debug) std::cout << PHWHERE << " Fatal Error - GlobalVertexMap node is missing"<< std::endl;
      if (debug) std::cout << "CaloAna GlobalVertexMap node is missing" << std::endl;
      h_vtxmap_fail->Fill(1);
      VertexMapFailcounter++;
      // return 0;
      //  return Fun4AllReturnCodes::ABORTRUN;
    }
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
        h_vtxmap_fail->Fill(1);
        VertexMapFailcounter++;
      }
    }
    else
    {
      if (debug) std::cout << "CaloAna GlobalVertexMap node is empty" << std::endl;
      h_vtxmap_fail->Fill(1);
      VertexMapFailcounter++;
    }
  }
  if (zvtxcut_bool && abs(vtx_z) > zvtx_cut_val)
  {
    h_cutCounter->Fill(12);
    return Fun4AllReturnCodes::EVENT_OK;
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

  RawClusterContainer::ConstRange clusterEnd = clusterContainer->getClusters();
  RawClusterContainer::ConstIterator clusterIter;
  RawClusterContainer::ConstIterator clusterIter2;
  int nClusCount = 0;
  for (clusterIter = clusterEnd.first; clusterIter != clusterEnd.second; clusterIter++)
  {
    RawCluster* recoCluster = clusterIter->second;

    CLHEP::Hep3Vector vertex(0, 0, vtx_z);
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

  // float smear = 0.00;
  if (debug) std::cout << " " << "Cluster Loop: 1 " << std::endl;
  for (clusterIter = clusterEnd.first; clusterIter != clusterEnd.second; clusterIter++)
  {
    RawCluster* recoCluster = clusterIter->second;

    CLHEP::Hep3Vector vertex(0, 0, vtx_z);
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

    TLorentzVector photon1;
    photon1.SetPtEtaPhiE(clus_pt, clus_eta, clus_phi, clusE);
    pi0smearvec[0] = SmearPhoton4vector(photon1, badcalibsmear);
    
    if (eTCutbool)
    {
      if (pi0smearvec[0].Et() < etcut && cutson)
      {
        h_cutCounter->Fill(2);
        continue;
      }
    }

    if (additionalsmearing)
    {
      if (Cluster_Debug) h_cluster_etaphi_cuts[3]->Fill(pi0smearvec[0].Eta(), pi0smearvec[0].Phi());

      //switched in prob cut for etcut
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

      if (Cluster_Debug)
      {
        if (filledClustersAfterCut1.insert(recoCluster).second)
        {
          h_cluster_etaphi_cuts[2]->Fill(clus_eta, clus_phi);
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
      TLorentzVector pi0;
      // TLorentzVector pi0 = photon1 + photon2;
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
      vector<TLorentzVector> truth_photons;
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
          h_reco_photon1E_weighted->Fill(photon1.E(), inv_yield);
          h_reco_photon2E_weighted->Fill(photon2.E(), inv_yield);
          h_reco_ALLphotonE_weighted->Fill(photon1.E(), inv_yield);
          h_reco_ALLphotonE_weighted->Fill(photon2.E(), inv_yield);

          h_InvMass_smear_weighted_eta_3d->Fill(pi0smearvec[2].Pt(), pi0smearvec[2].M(), pi0smearvec[2].Eta(), inv_yield);
          h_InvMass_smear_weighted_eta_2d->Fill(pi0smearvec[2].Eta(), pi0smearvec[2].M(), inv_yield);
          h_reco_etaphi_cuts[8]->Fill(pi0.Eta(), pi0.Phi(), inv_yield);

          if (filltruthspectrum && (matchmctruth))
          {
            float delR = pi0.DeltaR(myVector);
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
          if (filltruthspectrum && (matchmctruth))
          {
            if (id == 22)
            {
              TLorentzVector myPhotonVector;
              float photonE = truth->get_e();
              myPhotonVector.SetPxPyPzE(truth->get_px(), truth->get_py(), truth->get_pz(), photonE);
              if (additionalsmearing)
              {
                float delRP1 = pi0smearvec[0].DeltaR(myPhotonVector);
                float delRP2 = pi0smearvec[1].DeltaR(myPhotonVector);
                h_truthmatched_Photon_delR->Fill(delRP1);
                h_truthmatched_Photon_delR->Fill(delRP2);
                if (pi0smearvec[0].DeltaR(myPhotonVector) < 0.015)
                {
                  h_truthmatched_photon1E->Fill(photonE);
                  h_truthmatched_photon1E_weighted->Fill(photonE, inv_yield);
                  h_truthmatched_AllphotonE->Fill(photonE);
                  h_truthmatched_AllphotonE_weighted->Fill(photonE, inv_yield);
                }
                if (pi0smearvec[1].DeltaR(myPhotonVector) < 0.015)
                {
                  h_truthmatched_photon2E->Fill(photonE);
                  h_truthmatched_photon2E_weighted->Fill(photonE, inv_yield);
                  h_truthmatched_AllphotonE->Fill(photonE);
                  h_truthmatched_AllphotonE_weighted->Fill(photonE, inv_yield);
                }
              }
              else if (!additionalsmearing)
              {
                float delRP1 = photon1.DeltaR(myPhotonVector);
                float delRP2 = photon2.DeltaR(myPhotonVector);
                h_truthmatched_Photon_delR->Fill(delRP1);
                h_truthmatched_Photon_delR->Fill(delRP2);
                if (photon1.DeltaR(myPhotonVector) < 0.015)
                {
                  h_truthmatched_photon1E->Fill(photonE);
                  h_truthmatched_photon1E_weighted->Fill(photonE, inv_yield);
                  h_truthmatched_AllphotonE->Fill(photonE);
                  h_truthmatched_AllphotonE_weighted->Fill(photonE, inv_yield);
                }
                if (photon2.DeltaR(myPhotonVector) < 0.015)
                {
                  h_truthmatched_photon2E->Fill(photonE);
                  h_truthmatched_photon2E_weighted->Fill(photonE, inv_yield);
                  h_truthmatched_AllphotonE->Fill(photonE);
                  h_truthmatched_AllphotonE_weighted->Fill(photonE, inv_yield);
                }
              }
            }
          }

          //*/
        }
      }

      //*
      for (auto tr_phot : truth_photons)
      {
        float delR = photon1.DeltaR(tr_phot);
        if (debug) std::cout << delR << " " << std::endl;
        h_delR_recTrth->Fill(delR);
        if (delR < 0.0101)
        {  // choose this value based on looking at delR distribution
          float res = photon1.E() / tr_phot.E();
          h_matched_res->Fill(res, photon1.Eta());
        }
      }
      //*/
      if (debug) std::cout << " " << "truth: Loops Done " << std::endl;
      if (additionalsmearing)
      {
        h_InvMass_smear->Fill(pi0smearvec[2].M());
        h_InvMass_smear_2d->Fill(pi0smearvec[2].Pt(), pi0smearvec[2].M());
        h_InvMass_smear_eta_3d->Fill(pi0smearvec[2].Pt(), pi0smearvec[2].M(), pi0smearvec[2].Eta());
        h_InvMass_smear_eta_2d->Fill(pi0smearvec[2].Eta(), pi0smearvec[2].M());
      }
      h_pt1->Fill(photon1.Pt());
      h_pt2->Fill(photon2.Pt());
      h_InvMass_2d->Fill(pi0.Pt(), pi0.M());
      h_pion_pt->Fill(pi0.Pt());
      h_InvMass->Fill(pi0.M());
      h_reco_photon1E->Fill(photon1.E());
      h_reco_photon1E_2d->Fill(photon1.Pt(), photon1.E());
      h_reco_photon2E->Fill(photon2.E());
      h_reco_photon2E_2d->Fill(photon2.Pt(), photon2.E());
      h_reco_ALLphotonE->Fill(photon1.E());
      h_reco_ALLphotonE->Fill(photon2.E());
      h_reco_ALLphotonE_2d->Fill(photon1.Pt(), photon1.E());
      h_reco_ALLphotonE_2d->Fill(photon2.Pt(), photon2.E());

    }  // clusterIter2
  }  // clusteriter1 loop

  if (filltruthspectrum)
  {
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    if (truthinfo)
    {
      if (matchmctruth)
      {  // primaries
        PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
        for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
        {
          const PHG4Particle* truth = iter->second;
          if (!truthinfo->is_primary(truth)) continue;
          if (truth->get_pid() == 111 || (eta_weight && truth->get_pid() == 221))
          {
            float pion_pt = sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py());
            float pion_p = sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py() + truth->get_pz() * truth->get_pz());
            float pion_e = truth->get_e();
            float pion_phi = atan2(truth->get_py(), truth->get_px());
            float pion_eta = atanh(truth->get_pz() / sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py() + truth->get_pz() * truth->get_pz()));
            TLorentzVector truthpi0 = TLorentzVector();
            truthpi0.SetPtEtaPhiE(pion_pt, pion_eta, pion_phi, pion_e);
            h_truth_spectrum2->Fill(truthpi0.Pt());
            h_FullTruth_e->Fill(pion_e);
            h_FullTruth_eta->Fill(pion_eta);
            h_FullTruth_pt->Fill(pion_pt);
            h_FullTruth_p->Fill(pion_p);
          }
          if (truth->get_pid() == 221)
          {
            h_truth_etaspectrum->Fill(sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py()));
          }
          // photon loop
          if (truth->get_pid() == 22)
          {
            float photon_e = truth->get_e();
            h_truth_ALLphotonE->Fill(photon_e);
            // what should the photon weight be?
            // h_truth_ALLphotonE_weighted->Fill(photon_e, inv_yield);
          }
        }

        // secondaries
        PHG4TruthInfoContainer::Range second_range = truthinfo->GetSecondaryParticleRange();
        for (PHG4TruthInfoContainer::ConstIterator siter = second_range.first; siter != second_range.second; ++siter)
        {
          const PHG4Particle* truth = siter->second;
          int id = truth->get_pid();
          h_truth_pid_s->Fill(id);
          if (filltruthspectrum && (matchmctruth))
          {
            if (id == 22)
            {
              float photon_e = truth->get_e();
              float photon_pt = sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py());
              float photon_p = sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py() + truth->get_pz() * truth->get_pz()); /*

               float photon_phi = atan2(truth->get_py(), truth->get_px());
               float photon_eta = atanh(truth->get_pz() / sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py() + truth->get_pz() * truth->get_pz()));
               TLorentzVector truthphoton = TLorentzVector();
               truthphoton.SetPtEtaPhiE(photon_pt, photon_eta, photon_phi, photon_e);
               */
              h_truth_ALLphotonE->Fill(photon_e);
              h_truth_ALLphotonpt->Fill(photon_pt);
              h_truth_ALLphotonp->Fill(photon_p);
              // what should the photon weight be?
              // h_truth_ALLphotonE_weighted->Fill(photon_e, inv_yield);
            }
          }
        }
      }
    }
  }

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

TLorentzVector CaloAna::SmearPhoton4vector(TLorentzVector sourcephoton, double smearfactor)
{
  double smear = generateRandomNumber() * smearfactor + 1;
  TLorentzVector smearedphoton = sourcephoton * smear;
  return smearedphoton;
}
