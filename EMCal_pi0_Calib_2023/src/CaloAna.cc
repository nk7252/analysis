#include "CaloAna.h"

// G4Hits includes
#include <TLorentzVector.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>

// G4Cells includes
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>

// Tower includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

// Cluster includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

// MBD
#include <mbd/BbcGeom.h>
#include <mbd/MbdPmtContainerV1.h>
#include <mbd/MbdPmtHit.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TNtuple.h>
#include <TProfile2D.h>
#include <TTree.h>

#include <Event/Event.h>
#include <Event/packet.h>
#include <cassert>
#include <sstream>
#include <string>

#include <TMath.h>
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TRandom3.h"

#include <cdbobjects/CDBTTree.h>  // for CDBTTree
#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>

#include <g4main/PHG4TruthInfoContainer.h>

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

  outfile = new TFile(outfilename.c_str(), "RECREATE");
  // cutQA
  h_cutCounter = new TH1F("h_cutCounter", "Cut Counter", 13, 0.5, 13.5);
  // list of cuts
  //  clus1 chi2, clus1 cuts, tower eta>95,hotclus1,clus1=clus2,clus2 chi2, clus2 cuts, assym, Dr, pi0pt cut, hotclus2

  // correlation plots
  for (int i = 0; i < 96; i++)
  {
    h_mass_eta_lt[i] = new TH1F(Form("h_mass_eta_lt%d", i), "", 50, 0, 0.5);
    h_mass_eta_lt_rw[i] = new TH1F(Form("h_mass_eta_lt_rw%d", i), "", 50, 0, 0.5);
    h_pt_eta[i] = new TH1F(Form("h_pt_eta%d", i), "", 100, 0, 10);
    h_pt_eta_rw[i] = new TH1F(Form("h_pt_eta_rw%d", i), "", 100, 0, 10);
  }

  h_cemc_etaphi = new TH2F("h_cemc_etaphi", "", 96, 0, 96, 256, 0, 256);

  // 1D distributions
  // h_InvMass = new TH1F("h_InvMass", "Invariant Mass", 600, 0, 1.2);
  h_InvMass_w = new TH1F("h_InvMass_w", "Invariant Mass", 600, 0, 1.2);
  h_InvMassMix = new TH1F("h_InvMassMix", "Invariant Mass", 120, 0, 1.2);

  // cluster QA
  h_etaphi_clus = new TH2F("h_etaphi_clus", "", 140, -1.2, 1.2, 64, -1 * TMath::Pi(), TMath::Pi());
  h_clusE = new TH1F("h_clusE", "", 100, 0, 10);
  h_emcal_e_eta = new TH1F("h_emcal_e_eta", "", 96, 0, 96);

  h_pt1 = new TH1F("h_pt1", "", 100, 0, 5);
  h_pt2 = new TH1F("h_pt2", "", 100, 0, 5);
  h_pion_pt = new TH1F("h_pion_pt", "", 100, 0, 10);
  h_pion_pt_weight = new TH1F("h_pion_pt_weight", "", 100, 0, 10);

  h_nclusters = new TH1F("h_nclusters", "", 1000, 0, 1000);
  // Truth histos
  h_truth_eta = new TH1F("h_truth_eta", "", 100, -1.2, 1.2);
  h_truth_e = new TH1F("h_truth_e", "", 100, 0, 10);
  h_truth_pt = new TH1F("h_truth_pt", "", 100, 0, 10);
  h_truth_spectrum1 = new TH1F("h_truth_spectrum1", "", 10000, 0, 10);
  h_truth_spectrum2 = new TH1F("h_truth_spectrum", "", 10000, 0, 10);
  h_truth_pid_p = new TH1F("h_truth_pid_p", "Primary particle PIDs", 400, -200, 200);
  h_truth_pid_s = new TH1F("h_truth_pid_s", "Secondary particle PIDs", 400, -200, 200);
  h_delR_recTrth = new TH1F("h_delR_recTrth", "", 1000, 0, 5);
  h_delR_pionrecTrth = new TH1F("h_delR_recTrth", "", 5000, 0, 5);
  h_matched_res = new TH2F("h_matched_res", "", 100, 0, 1.5, 20, -1, 1);

  // pT differential Inv Mass
  h_InvMass = new TH1F("h_InvMass", "Invariant Mass", 600, 0, 1.2);
  h_InvMass_weighted = new TH1F("h_InvMass_weighted", "Invariant Mass, weighted WSHP", 600, 0, 1.2);
  h_inv_yield = new TH1F("h_inv_yield", "Invariant Yield distribution", 100, 0, 1e13);
  h_truthmatched_mass = new TH1F("h_truthmatched_mass", "Invariant Mass, truth matched(delR<0.015)", 600, 0, 1.2);
  h_truthmatched_mass2 = new TH1F("h_truthmatched_mass2", "Invariant Mass, truth matched(delR<0.1)", 600, 0, 1.2);
  h_truthmatched_mass3 = new TH1F("h_truthmatched_mass3", "Invariant Mass, truth matched(delR<0.2)", 600, 0, 1.2);
  h_InvMass_2d = new TH2F("h_InvMass_2d", "pT vs Invariant Mass", 4 * 10, 0, 10, 600, 0, 1.2);
  h_truthmatched_mass2_2d = new TH2F("h_truthmatched_mass2_2d", "pT vs Invariant Mass, truth matched(delR<0.1)", 4 * 10, 0, 10, 600, 0, 1.2);
  h_truthmatched_mass3_2d = new TH2F("h_truthmatched_mass3_2d", "pT vs Invariant Mass, truth matched(delR<0.25)", 4 * 10, 0, 10, 600, 0, 1.2);

  // high mass tail diagnostic
  std::vector<std::string> HistList = {"photon1", "photon2", "all photons", "pions"};
  for (int i = 0; i < 4; i++)
  {
    // eta and phi distributions
    h_phidist_InvMass_under200M[i] = new TH1F(Form("h_phidist_%s_InvMass_under200M", HistList[i].c_str()), Form("Phi dist for Inv mass under 200 MeV:%s", HistList[i].c_str()), 140, -7 / 8 * TMath::Pi(), 7 / 8 * TMath::Pi());
    h_phidist_InvMass_over200M[i] = new TH1F(Form("h_phidist_%s_InvMass_over200M", HistList[i].c_str()), Form("Phi dist for Inv mass over 200 MeV:%s", HistList[i].c_str()), 140, -7 / 8 * TMath::Pi(), 7 / 8 * TMath::Pi());

    h_etadist_InvMass_under200M[i] = new TH1F(Form("h_etadist_%s_InvMass_under200M", HistList[i].c_str()), Form("Eta dist for Inv mass under 200 MeV:%s", HistList[i].c_str()), 140, -1.2, 1.2);
    h_etadist_InvMass_over200M[i] = new TH1F(Form("h_etadist_%s_InvMass_over200M", HistList[i].c_str()), Form("Eta dist for Inv mass over 200 MeV:%s", HistList[i].c_str()), 140, -1.2, 1.2);

    // eta-phi distributions
    h_etaphidist_InvMass_under200M[i] = new TH2F(Form("h_etaphidist_%s_InvMass_under200M", HistList[i].c_str()), Form("Eta-Phi dist for Inv mass under 200 MeV:%s", HistList[i].c_str()), 24, -1.2, 1.2, 64, -1 * TMath::Pi(), TMath::Pi());
    h_etaphidist_InvMass_over200M[i] = new TH2F(Form("h_etaphidist_%s_InvMass_over200M", HistList[i].c_str()), Form("Eta-Phi dist for Inv mass over 200 MeV:%s", HistList[i].c_str()), 24, -1.2, 1.2, 64, -1 * TMath::Pi(), TMath::Pi());  // eta used to be 140
  }

  // dist of deta or dphi
  h_Dphidist_InvMass_under200M = new TH1F("h_Dphidist_InvMass_under200M", "Delta Phi dist for Inv mass under 200 MeV", 90, -1 * TMath::Pi() / 4, 1 * TMath::Pi() / 4);
  h_Dphidist_InvMass_over200M = new TH1F("h_Dphidist_InvMass_over200M", "Delta Phi dist for Inv mass over 200 MeV", 90, -1 * TMath::Pi() / 4, 1 * TMath::Pi() / 4);
  h_Detadist_InvMass_under200M = new TH1F("h_Detadist_InvMass_under200M", "Delta Eta dist for Inv mass under 200 MeV", 140, -1.2, 1.2);
  h_Detadist_InvMass_over200M = new TH1F("h_Detadist_InvMass_over200M", "Delta Eta dist for Inv mass over 200 MeV", 140, -1.2, 1.2);

  pidcuts = {0.001, 0.005, 0.01, 0.05, 0.1, 1};  // GeV? pretty sure that is the case

  for (int i = 0; i < 6; i++)
  {
    h_truth_pid_cuts[i] = new TH1F(Form("h_truth_pid_cut_%f", pidcuts[i]), Form("truth pid cut at %f MeV", pidcuts[i]), 150, -30, 120);
  }

  //////////////////////////
  // pT rewieghting
  // frw = new TFile("/sphenix/u/bseidlitz/work/analysis/EMCal_pi0_Calib_2023/macros/rw_pt.root");
  // for (int i = 0; i < 96; i++) h_pt_rw[i] = (TH1F*) frw->Get(Form("h_pt_eta%d", i));

  rnd = new TRandom3();
  // smearing added SPMC
  badcalibsmear = static_cast<float>(badcalibsmearint) / 1000.0f;

  h_InvMass_smear = new TH1F(Form("h_InvMass_smear_%d", badcalibsmearint), Form("Invariant Mass + const smear: %f percent", badcalibsmearint / 10.0f), 600, 0, 1.2);

  h_InvMass_smear_weighted = new TH1F(Form("h_InvMass_smear_weighted_%d", badcalibsmearint), Form("Invariant Mass + const smear, weighted: %f percent", badcalibsmearint / 10.0f), 600, 0, 1.2);

  // 2d variations
  h_InvMass_smear_weighted_2d = new TH2F(Form("h_InvMass_smear_weighted_2d_%d", badcalibsmearint), Form("pT vs Invariant Mass + const smear, weighted: %f percent", badcalibsmearint / 10.0f), 4 * 10, 0, 10, 600, 0, 1.2);

  h_InvMass_smear_2d = new TH2F(Form("h_InvMass_smear_2d_%d", badcalibsmearint), Form("pT vs Invariant Mass + const smear: %f percent", badcalibsmearint / 10.0f), 4 * 10, 0, 10, 600, 0, 1.2);

  std::vector<std::string> RestrictEtaCuts = {"Low_Eta", "Mid_Eta", "High_Eta"};
  etaRanges = {
      {0.0, 0.2},  // pair1
      {0.2, 0.4},  // pair2
      {0.4, 1.1}   // pair3
  };
  for (int i = 0; i < 3; i++)
  {
    h_InvMass_smear_risingpt[i] = new TH1F(
        Form("h_InvMass_smear_risingpt_%d_%s", badcalibsmearint, RestrictEtaCuts[i].c_str()),
        Form("Invariant Mass, rising_pt+%s+smearing: %f percent", RestrictEtaCuts[i].c_str(), badcalibsmearint / 10.0f), 120, 0, 0.6);

    h_InvMass_smear_fallingpt[i] = new TH1F(
        Form("h_InvMass_smear_fallingpt_%d_%s", badcalibsmearint, RestrictEtaCuts[i].c_str()),
        Form("Invariant Mass, falling_pt+%s+smearing: %f percent", RestrictEtaCuts[i].c_str(), badcalibsmearint / 10.0f), 120, 0, 0.6);

    h_InvMass_smear_flatpt[i] = new TH1F(
        Form("h_InvMass_smear_flatpt_%d_%s", badcalibsmearint, RestrictEtaCuts[i].c_str()),
        Form("Invariant Mass, flat_pt+%s+smearing: %f percent", RestrictEtaCuts[i].c_str(), badcalibsmearint / 10.0f), 120, 0, 0.6);
  }

  if (poscor == true)
  {
    clustposcorstring = "CLUSTER_POS_COR_CEMC";
  }
  else
  {
    clustposcorstring = "CLUSTER_CEMC";
  }

  funkyCaloStuffcounter = 0;
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

  // cuts
  if (debug) std::cout << " " << "Cuts " << std::endl;
  float maxDr = 1.1;
  float maxAlpha = 0.6;
  float clus_chisq_cut = 4;
  float nClus_ptCut = 0.5;  // 0.5
  float pi0ptcutfactor = 0;
  float ptMaxCut = 20;     // 7 in data? ** keep this in mind. 3 may make more sense, but 7 is
  float pt1ClusCut = 1.3;  // centrality dependence cuts 2.2 for both // 1.3
  float pt2ClusCut = 0.7;  // // 0.7

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
      // std::cout << PHWHERE << " Fatal Error - GlobalVertexMap node is missing"<< std::endl;
      std::cout << "CaloAna GlobalVertexMap node is missing" << std::endl;
      // return Fun4AllReturnCodes::ABORTRUN;
    }
    if (vertexmap && !vertexmap->empty())
    {
      GlobalVertex* vtx = vertexmap->begin()->second;
      if (vtx)
      {
        vtx_z = vtx->get_z();
      }
    }
  }

  //////////////////////////////////////////////
  //         towers
  vector<float> ht_eta;
  vector<float> ht_phi;

  // if (!m_vtxCut || abs(vtx_z) > _vz)  return Fun4AllReturnCodes::EVENT_OK;

  TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");

  RawClusterContainer* clusterContainer = findNode::getClass<RawClusterContainer>(topNode, Form("%s", clustposcorstring.c_str()));
  // changed from CLUSTERINFO_CEMC2
  // Blair using "CLUSTER_POS_COR_CEMC" now. change from CLUSTER_CEMC
  // RawClusterContainer* clusterContainer = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_POS_COR_CEMC");
  if (!clusterContainer)
  {
    std::cout << PHWHERE << "funkyCaloStuff::process_event - Fatal Error - CLUSTER_CEMC node is missing. " << std::endl;
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

    if (clus_pt < nClus_ptCut && cutson) continue;
    if (clus_chisq > clus_chisq_cut && cutson) continue;

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
    CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*recoCluster, vertex);
    std::vector<TLorentzVector> pi0smearvec(3);

    float clusE = E_vec_cluster.mag();
    float clus_eta = E_vec_cluster.pseudoRapidity();
    float clus_phi = E_vec_cluster.phi();
    float clus_pt = E_vec_cluster.perp();
    // clus_pt *= rnd->Gaus(1, smear);
    float clus_chisq = recoCluster->get_chi2();
    if (clus_chisq > clus_chisq_cut && cutson)
    {
      h_cutCounter->Fill(1);
      continue;
    }
    TLorentzVector photon1;
    photon1.SetPtEtaPhiE(clus_pt, clus_eta, clus_phi, clusE);

    pi0smearvec[0] = SmearPhoton4vector(photon1, badcalibsmear);

    if (additionalsmearing)
    {
      if ((pi0smearvec[0].Pt() < pt1ClusCut || pi0smearvec[0].Pt() > ptMaxCut) && cutson)
      {
        h_cutCounter->Fill(2);
        continue;
      }
    }
    else if (!additionalsmearing)
    {
      if ((photon1.Pt() < pt1ClusCut || photon1.Pt() > ptMaxCut) && cutson)
      {
        h_cutCounter->Fill(2);
        continue;
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
    h_pt_eta[lt_eta]->Fill(clus_pt);

    if (dynMaskClus && hotClus == true && cutson)
    {
      h_cutCounter->Fill(4);
      continue;
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
      RawCluster* recoCluster2 = clusterIter2->second;

      CLHEP::Hep3Vector E_vec_cluster2 = RawClusterUtility::GetECoreVec(*recoCluster2, vertex);

      float clus2E = E_vec_cluster2.mag();
      float clus2_eta = E_vec_cluster2.pseudoRapidity();
      float clus2_phi = E_vec_cluster2.phi();
      float clus2_pt = E_vec_cluster2.perp();
      float clus2_chisq = recoCluster2->get_chi2();

      if (clus2_chisq > clus_chisq_cut && cutson)
      {
        h_cutCounter->Fill(6);
        continue;
      }
      TLorentzVector photon2;
      photon2.SetPtEtaPhiE(clus2_pt, clus2_eta, clus2_phi, clus2E);
      pi0smearvec[1] = SmearPhoton4vector(photon2, badcalibsmear);
      TLorentzVector pi0 = photon1 + photon2;
      pi0smearvec[2] = pi0smearvec[0] + pi0smearvec[1];

      if (additionalsmearing)
      {
        if ((pi0smearvec[1].Pt() < pt2ClusCut || pi0smearvec[1].Pt() > ptMaxCut) && cutson)
        {
          h_cutCounter->Fill(7);
          continue;
        }
        if (fabs(pi0smearvec[0].E() - pi0smearvec[1].E()) / (pi0smearvec[0].E() + pi0smearvec[1].E()) > maxAlpha && cutson)
        {
          h_cutCounter->Fill(8);
          continue;
        }
        if (pi0smearvec[0].DeltaR(pi0smearvec[1]) > maxDr && cutson)
        {
          h_cutCounter->Fill(9);
          continue;
        }
        if (pi0smearvec[2].Pt() < pi0ptcut)
        {
          h_cutCounter->Fill(10);
          continue;
        }
      }
      else if (!additionalsmearing)
      {
        if ((photon2.Pt() < pt2ClusCut || photon2.Pt() > ptMaxCut) && cutson)
        {
          h_cutCounter->Fill(7);
          continue;
        }
        if (fabs(photon1.E() - photon2.E()) / (photon1.E() + photon2.E()) > maxAlpha && cutson)
        {
          h_cutCounter->Fill(8);
          continue;
        }
        if (photon1.DeltaR(photon2) > maxDr && cutson)
        {
          h_cutCounter->Fill(9);
          continue;
        }
        if (pi0.Pt() < pi0ptcut)
        {
          h_cutCounter->Fill(10);
          continue;
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
      h_etaphi_clus->Fill(clus_eta, clus_phi);

      if (dynMaskClus && hotClus2 == true && cutson)
      {
        h_cutCounter->Fill(11);
        continue;
      }

      /////////////////////////////////////////////////
      //// Truth info
      float weight = 1;
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
          myVector.SetXYZM(truth->get_px(), truth->get_py(), truth->get_pz(), 0.13497);
          float energy = myVector.E();
          weight = myVector.Pt() * TMath::Exp(-3 * myVector.Pt());
          h_truth_e->Fill(energy, weight);
          h_truth_eta->Fill(myVector.Eta(), weight);
          // h_truth_pt->Fill(myVector.Pt());
          h_truth_pt->Fill(myVector.Pt(), weight);
          int id = truth->get_pid();
          h_truth_pid_p->Fill(id);
          //--------------------Alternative paramaterization, woods saxon + hagedorn + power law

          double t = 4.5;
          double w = 0.114;
          double A = 229.6;
          double B = 14.43;
          double n = 8.1028;
          double m_param = 10.654;
          double p0 = 1.466;
          double Pt = myVector.Pt();
          double weight_function = ((1 / (1 + exp((Pt - t) / w))) * A / pow(1 + Pt / p0, m_param) + (1 - (1 / (1 + exp((Pt - t) / w)))) * B / (pow(Pt, n)));
          inv_yield = WeightScale * Pt * weight_function;  //
          // h_pion_pt_weight->Fill(pi0.Pt(), inv_yield);
          h_inv_yield->Fill(inv_yield);
          h_InvMass_weighted->Fill(pi0.M(), inv_yield);
          h_InvMass_smear_weighted->Fill(pi0smearvec[2].M(), inv_yield);
          h_InvMass_smear_weighted_2d->Fill(pi0smearvec[2].Pt(), pi0smearvec[2].M(), inv_yield);
          if (debug) std::cout << "truth pt=" << Pt << "   weight function=" << weight_function << "  inv_yield=" << inv_yield << std::endl;
          if (additionalsmearing)
          {
            for (size_t i = 0; i < 3; i++)  // break up inv mass spectrum if debugging.
            {
              if (std::abs(myVector.Eta()) > etaRanges[i].first && std::abs(myVector.Eta()) < etaRanges[i].second)
              {
                if (Pt < 4.75)
                {
                  h_InvMass_smear_risingpt[i]->Fill(pi0smearvec[2].M());
                }
                else if (4.75 < Pt && Pt < 5.25)
                {
                  h_InvMass_smear_flatpt[i]->Fill(pi0smearvec[2].M());
                }
                else
                {
                  h_InvMass_smear_fallingpt[i]->Fill(pi0smearvec[2].M());
                }
              }
            }
          }

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

          ///*if (filltruthspectrum && matchmctruth))
          {
            TLorentzVector truthpi0 = TLorentzVector();
            float pion_pt = sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py());
            // float pion_pz = truth->get_pz();
            float pion_e = truth->get_e();
            float pion_phi = atan2(truth->get_py(), truth->get_px());
            float pion_eta = atanh(truth->get_pz() / sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py() + truth->get_pz() * truth->get_pz()));
            truthpi0.SetPtEtaPhiE(pion_pt, pion_eta, pion_phi, pion_e);
            float delR = pi0.DeltaR(truthpi0);
            h_delR_pionrecTrth->Fill(delR);
            if (delR < 0.015)
            {
              h_truth_spectrum1->Fill(truthpi0.Pt());
            }
            //
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
      }
      h_pt1->Fill(photon1.Pt());
      h_pt2->Fill(photon2.Pt());
      h_InvMass_2d->Fill(pi0.Pt(), pi0.M());
      h_pion_pt->Fill(pi0.Pt());
      h_InvMass->Fill(pi0.M());
    }  // clusterIter2
  }  // clusteriter1 loop

  if ((filltruthspectrum && !matchmctruth) || (filltruthspectrum && Fillanyways))
  {
    PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    if (truthinfo)
    {
      // secondaries
      PHG4TruthInfoContainer::Range second_range = truthinfo->GetSecondaryParticleRange();
      for (PHG4TruthInfoContainer::ConstIterator siter = second_range.first; siter != second_range.second; ++siter)
      {
        const PHG4Particle* truth = siter->second;
        //int id = truth->get_pid();
        if (truth->get_pid() == 111)
        {
          float pion_pt = sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py());
          float pion_e = truth->get_e();
          float pion_phi = atan2(truth->get_py(), truth->get_px());
          float pion_eta = atanh(truth->get_pz() / sqrt(truth->get_px() * truth->get_px() + truth->get_py() * truth->get_py() + truth->get_pz() * truth->get_pz()));
          TLorentzVector truthpi0 = TLorentzVector();
          truthpi0.SetPtEtaPhiE(pion_pt, pion_eta, pion_phi, pion_e);
          h_truth_spectrum2->Fill(truthpi0.Pt());
          // h_truth_pion_pt->Fill(truthpi0.Pt());
          // h_truth_pion_eta->Fill(truthpi0.Eta());
          // h_truth_pion_e->Fill(truthpi0.E());
          // h_truth_pion_phi->Fill(truthpi0.Phi());
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
  std::cout << "funkycounter: " << funkyCaloStuffcounter << std::endl;
  return 0;
}

float CaloAna::getWeight(int ieta, float pt)
{
  if (ieta < 0 || ieta > 95) return 0;
  float val = h_pt_rw[ieta]->GetBinContent(h_pt_rw[ieta]->FindBin(pt));
  if (val == 0) return 0;
  return 1 / val;
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
  for (int i = 0; i < 96; i++)
  {
    h_M_eta[i] = (TH1F*) fin->Get(Form("h_mass_eta_lt_rw%d", i));
    h_M_eta[i]->Scale(1. / h_M_eta[i]->Integral(), "width");
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
