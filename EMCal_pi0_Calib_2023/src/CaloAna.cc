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

#include <iostream>
#include <utility>
#include <vector>
#include <string>
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include <TMath.h>

#include <cdbobjects/CDBTTree.h>  // for CDBTTree
#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>

#include <g4main/PHG4TruthInfoContainer.h>

R__LOAD_LIBRARY(libLiteCaloEvalTowSlope.so)

using namespace std;

CaloAna::CaloAna(const std::string& name, const std::string& filename)
  : SubsysReco(name)
  , detector("HCALIN")
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
  hm = new Fun4AllHistoManager(Name());
  // create and register your histos (all types) here

  outfile = new TFile(outfilename.c_str(), "RECREATE");

  // correlation plots
  for (int i = 0; i < 96; i++)
  {
    h_mass_eta_lt[i] = new TH1F(Form("h_mass_eta_lt%d", i), "", 50, 0, 0.5);
    h_pt_eta[i] = new TH1F(Form("h_pt_eta%d", i), "", 100, 0, 10);
  }

  h_cemc_etaphi = new TH2F("h_cemc_etaphi", "", 96, 0, 96, 256, 0, 256);

  // 1D distributions
  h_InvMass = new TH1F("h_InvMass", "Invariant Mass", 120, 0, 1.2);
  h_InvMassMix = new TH1F("h_InvMassMix", "Invariant Mass", 120, 0, 1.2);

  // cluster QA
  h_etaphi_clus = new TH2F("h_etaphi_clus", "", 140, -1.2, 1.2, 64, -1 * TMath::Pi(), TMath::Pi());
  h_clusE = new TH1F("h_clusE", "", 100, 0, 10);

  h_emcal_e_eta = new TH1F("h_emcal_e_eta", "", 96, 0, 96);

  h_pt1 = new TH1F("h_pt1", "", 100, 0, 5);
  h_pt2 = new TH1F("h_pt2", "", 100, 0, 5);

  h_nclusters = new TH1F("h_nclusters", "", 1000, 0, 1000);
  // Truth histos
  h_truth_eta = new TH1F("h_truth_eta", "", 100, -1.2, 1.2);
  h_truth_e = new TH1F("h_truth_e", "", 100, 0, 10);
  h_truth_pt = new TH1F("h_truth_pt", "", 100, 0, 10);
  h_truth_pid = new TH1F("h_truth_pid", "", 150, -30, 120);

  // pT differential Inv Mass
  h_InvMass = new TH1F("h_InvMass", "Invariant Mass", 120, 0, 0.6);
  h_pTdiff_InvMass = new TH2F("h_pTdiff_InvMass", "Invariant Mass", 2 * 64, 0, 64, 100, 0, 1.2);

  // vector for bad calib smearing.
  badcalibsmear={1,3,5,10};
  // high mass tail diagnostic
  std::vector<std::string> HistList={"photon1","photon2","all photons","pions"};
  for(int i=0; i<4;i++){
    h_InvMass_badcalib_smear[i]= new TH1F(Form("h_InvMass_badcalib_smear_%d",badcalibsmear[i]), Form("Invariant Mass with 'bad calibration' smearing applied:%d",badcalibsmear[i]), 120, 0, 0.6);

    // eta and phi distributions
    h_phidist_InvMass_under200M[i] = new TH1F(Form("h_phidist_%s_InvMass_under200M",HistList[i].c_str()), Form("Phi dist for Inv mass under 200 MeV:%s",HistList[i].c_str()) , 140, -7/8 * TMath::Pi(), 7/8 * TMath::Pi());
    h_phidist_InvMass_over200M[i] = new TH1F(Form("h_phidist_%s_InvMass_over200M",HistList[i].c_str()), Form("Phi dist for Inv mass over 200 MeV:%s",HistList[i].c_str()), 140, -7/8 * TMath::Pi(), 7/8 * TMath::Pi());

    h_etadist_InvMass_under200M[i] = new TH1F(Form("h_etadist_%s_InvMass_under200M",HistList[i].c_str()),Form("Eta dist for Inv mass under 200 MeV:%s",HistList[i].c_str()) , 140, -1.2, 1.2);
    h_etadist_InvMass_over200M[i] = new TH1F(Form("h_etadist_%s_InvMass_over200M",HistList[i].c_str()), Form("Eta dist for Inv mass over 200 MeV:%s",HistList[i].c_str()), 140, -1.2, 1.2);

    // eta-phi distributions
    h_etaphidist_InvMass_under200M[i] = new TH2F(Form("h_etaphidist_%s_InvMass_under200M",HistList[i].c_str()),Form("Eta-Phi dist for Inv mass under 200 MeV:%s",HistList[i].c_str()) , 24, -1.2, 1.2, 64, -1 * TMath::Pi(), TMath::Pi());
    h_etaphidist_InvMass_over200M[i] = new TH2F(Form("h_etaphidist_%s_InvMass_over200M",HistList[i].c_str()),Form("Eta-Phi dist for Inv mass over 200 MeV:%s",HistList[i].c_str()) , 24, -1.2, 1.2, 64, -1 * TMath::Pi(), TMath::Pi());  // eta used to be 140
  }
  
  // dist of deta or dphi
  h_Dphidist_InvMass_under200M = new TH1F("h_Dphidist_InvMass_under200M","Delta Phi dist for Inv mass under 200 MeV", 90, -1 * TMath::Pi() / 4, 1 * TMath::Pi() / 4);
  h_Dphidist_InvMass_over200M = new TH1F("h_Dphidist_InvMass_over200M","Delta Phi dist for Inv mass over 200 MeV", 90, -1 * TMath::Pi() / 4, 1 * TMath::Pi() / 4);
  h_Detadist_InvMass_under200M = new TH1F("h_Detadist_InvMass_under200M","Delta Eta dist for Inv mass under 200 MeV", 140, -1.2, 1.2);
  h_Detadist_InvMass_over200M = new TH1F("h_Detadist_InvMass_over200M","Delta Eta dist for Inv mass over 200 MeV", 140, -1.2, 1.2);

  std::vector<float> pidcuts ={0.5,1,5,10,20,50};
  for(int i=0; i<6; i++){
    h_truth_pid_cuts[i]= new TH1F(Form("h_truth_pid_cut_%f",pidcuts[i]), Form("truth pid cut at %f MeV",pidcuts[i]), 150, -30, 120); 
  }

  funkyCaloStuffcounter = 0;
  return 0;
}

int CaloAna::process_event(PHCompositeNode* topNode)
{
  _eventcounter++;

  process_towers(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloAna::process_towers(PHCompositeNode* topNode)
{
  if ((_eventcounter % 1000) == 0) std::cout << _eventcounter << std::endl;

  // float emcaldownscale = 1000000 / 800;

  float emcal_hit_threshold = 0.5;  // GeV

  // cuts
  bool cutson = true;
  // if(cutson){std::cout << "Cuts are on" << std::endl;}else{std::cout << "Cuts are off" << std::endl;}
  float maxDr = 1.1;
  float maxAlpha = 0.6;
  float clus_chisq_cut = 4;
  float nClus_ptCut = 0.5;
  int max_nClusCount = 75;

  //----------------------------------get vertex------------------------------------------------------//

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

  vector<float> ht_eta;
  vector<float> ht_phi;

  // if (!m_vtxCut || abs(vtx_z) > _vz)  return Fun4AllReturnCodes::EVENT_OK;

  TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
  if (towers)
  {
    int size = towers->size();  // online towers should be the same!
    for (int channel = 0; channel < size; channel++)
    {
      TowerInfo* tower = towers->get_tower_at_channel(channel);
      float offlineenergy = tower->get_energy();
      unsigned int towerkey = towers->encode_key(channel);
      int ieta = towers->getTowerEtaBin(towerkey);
      int iphi = towers->getTowerPhiBin(towerkey);
      bool isGood = !(tower->get_isBadChi2());
      if (!isGood && offlineenergy > 0.2)
      {
        ht_eta.push_back(ieta);
        ht_phi.push_back(iphi);
      }
      if (isGood) h_emcal_e_eta->Fill(ieta, offlineenergy);
      if (offlineenergy > emcal_hit_threshold)
      {
        h_cemc_etaphi->Fill(ieta, iphi);
      }
    }
  }

  RawClusterContainer* clusterContainer = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_POS_COR_CEMC");  // changed from CLUSTERINFO_CEMC2
  // Blair using "CLUSTER_POS_COR_CEMC" now. change from CLUSTER_CEMC
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

  if (nClusCount > max_nClusCount && cutson) return Fun4AllReturnCodes::EVENT_OK;

  float ptMaxCut = 3;  // 7 in data? ** keep this in mind.

  // float ptClusMax = 7;
  float pt1ClusCut = 1.3;  // 1.3
  float pt2ClusCut = 0.7;  // 0.7

  // if (nClusCount > 30)// no cluster dependent cut.
  //{
  //   pt1ClusCut += 1.4 * (nClusCount - 29) / 200.0;
  //   pt2ClusCut += 1.4 * (nClusCount - 29) / 200.0;
  // }

  float pi0ptcut = 1.22 * (pt1ClusCut + pt2ClusCut);

  vector<float> save_pt;
  vector<float> save_eta;
  vector<float> save_phi;
  vector<float> save_e;

  for (clusterIter = clusterEnd.first; clusterIter != clusterEnd.second; clusterIter++)
  {
    RawCluster* recoCluster = clusterIter->second;

    CLHEP::Hep3Vector vertex(0, 0, vtx_z);
    CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*recoCluster, vertex);
    std::vector<TLorentzVector> pi0gammavec(3);
    std::vector<TLorentzVector> pi0smearvec(4);// only filled with pions. each is a different level of smearing. smearing level is defined in init(?)


    float clusE = E_vec_cluster.mag();
    float clus_eta = E_vec_cluster.pseudoRapidity();
    float clus_phi = E_vec_cluster.phi();
    float clus_pt = E_vec_cluster.perp();
    float clus_chisq = recoCluster->get_chi2();
    h_clusE->Fill(clusE);
    // std::cout << "clusE = " << clusE <<  " clus_eta = " << clus_eta <<  " clus_phi = " << clus_phi <<  " clus_pt = " << clus_pt <<  " clus_chisq = " << clus_chisq << std::endl;

    if (clus_chisq > clus_chisq_cut && cutson) continue;

    // loop over the towers in the cluster
    RawCluster::TowerConstRange towerCR = recoCluster->get_towers();
    RawCluster::TowerConstIterator toweriter;
    bool hotClus = false;
    float lt_e = -1000;
    unsigned int lt_eta = -1;
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

      for (size_t i = 0; i < ht_eta.size(); i++)
        if (towerphi == ht_phi[i] && towereta == ht_eta[i])
          hotClus = true;
    }

    if (lt_eta > 95) continue;

    h_pt_eta[lt_eta]->Fill(clus_pt);

    if (dynMaskClus && hotClus == true && cutson) continue;

    h_etaphi_clus->Fill(clus_eta, clus_phi);

    TLorentzVector photon1;
    photon1.SetPtEtaPhiE(clus_pt, clus_eta, clus_phi, clusE);

    if ((clus_pt < pt1ClusCut || clus_pt > ptMaxCut) && cutson) continue;
    // || clus_pt > ptClusMax this was in the cuts to data.
    for (clusterIter2 = clusterEnd.first; clusterIter2 != clusterEnd.second; clusterIter2++)
    {
      if (clusterIter == clusterIter2)
      {
        continue;
      }
      RawCluster* recoCluster2 = clusterIter2->second;

      CLHEP::Hep3Vector E_vec_cluster2 = RawClusterUtility::GetECoreVec(*recoCluster2, vertex);

      float clus2E = E_vec_cluster2.mag();
      float clus2_eta = E_vec_cluster2.pseudoRapidity();
      float clus2_phi = E_vec_cluster2.phi();
      float clus2_pt = E_vec_cluster2.perp();
      float clus2_chisq = recoCluster2->get_chi2();

      if ((clus2_pt < pt2ClusCut || clus2_pt > ptMaxCut) && cutson) continue;
      // || clus2_pt > ptClusMax is found in the data.
      if (clus2_chisq > clus_chisq_cut && cutson) continue;
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

      if (dynMaskClus && hotClus2 == true && cutson) continue;

      TLorentzVector photon2;
      photon2.SetPtEtaPhiE(clus2_pt, clus2_eta, clus2_phi, clus2E);

      if (fabs(clusE - clus2E) / (clusE + clus2E) > maxAlpha && cutson) continue;

      if (photon1.DeltaR(photon2) > maxDr && cutson) continue;

      TLorentzVector pi0 = photon1 + photon2;
      pi0gammavec[0]=photon1;//photon1
      pi0gammavec[1]=photon2;//photon2
      pi0gammavec[2]=pi0;//pion

      if (pi0.Pt() < pi0ptcut) continue;

      // maybe need two more histograms for safety.
      // is it possible for them to be bent more than 180 from each other?
      //  does it matter here?
      float dphi = clus_phi - clus2_phi;  // for now ensure it is -pi to pi
      if (dphi > 3.14159) dphi -= 2 * 3.14159;
      if (dphi < -3.14159) dphi += 2 * 3.14159;

      float deta = clus_eta - clus2_eta;

      if (pi0.M() > 0.2)
      {
        // need to fix the fill loops 
        h_Detadist_InvMass_over200M->Fill(deta);
        h_Dphidist_InvMass_over200M->Fill(dphi);


        for (size_t i = 0; i < 4; ++i) {//h_phidist_InvMass_under200M.size()
          // For the first three histograms, fill with the corresponding TLorentzVector
          if (i < 3) {
            //histograms[i]->Fill(tlorentzVectors[i].Phi()); // Example property
            h_phidist_InvMass_over200M[i]->Fill(pi0gammavec[i].Phi()); 
            h_etadist_InvMass_over200M[i]->Fill(pi0gammavec[i].Eta()); 
            // eta-phi distributions
            h_etaphidist_InvMass_over200M[i]->Fill(pi0gammavec[i].Eta(),pi0gammavec[i].Phi()); 
          } 
          // Special case: 4th histogram gets filled with both the 1st and 2nd TLorentzVectors
          else if (i == 3) {
            //ph1
            h_phidist_InvMass_over200M[i]->Fill(pi0gammavec[0].Phi()); 
            h_etadist_InvMass_over200M[i]->Fill(pi0gammavec[0].Eta()); 
            // eta-phi distributions
            h_etaphidist_InvMass_over200M[i]->Fill(pi0gammavec[0].Eta(),pi0gammavec[0].Phi()); 
            //ph2
            h_phidist_InvMass_over200M[i]->Fill(pi0gammavec[1].Phi()); 
            h_etadist_InvMass_over200M[i]->Fill(pi0gammavec[1].Eta()); 
            // eta-phi distributions
            h_etaphidist_InvMass_over200M[i]->Fill(pi0gammavec[1].Eta(),pi0gammavec[1].Phi()); 
          }
        }
      }
      else
      {
        h_Detadist_InvMass_under200M->Fill(deta);
        h_Dphidist_InvMass_under200M->Fill(dphi);

       for (size_t i = 0; i < 4; ++i) {
          // For the first three histograms, fill with the corresponding TLorentzVector
          if (i < 3) {
            //histograms[i]->Fill(tlorentzVectors[i].Phi()); // Example property
            h_phidist_InvMass_under200M[i]->Fill(pi0gammavec[i].Phi()); 
            h_etadist_InvMass_under200M[i]->Fill(pi0gammavec[i].Eta()); 
            // eta-phi distributions
            h_etaphidist_InvMass_under200M[i]->Fill(pi0gammavec[i].Eta(),pi0gammavec[i].Phi()); 
          } 
          // Special case: 4th histogram gets filled with both the 1st and 2nd TLorentzVectors
          else if (i == 3) {
            //ph1
            h_phidist_InvMass_under200M[i]->Fill(pi0gammavec[0].Phi()); 
            h_etadist_InvMass_under200M[i]->Fill(pi0gammavec[0].Eta()); 
            // eta-phi distributions
            h_etaphidist_InvMass_under200M[i]->Fill(pi0gammavec[0].Eta(),pi0gammavec[0].Phi()); 
            //ph2
            h_phidist_InvMass_under200M[i]->Fill(pi0gammavec[1].Phi()); 
            h_etadist_InvMass_under200M[i]->Fill(pi0gammavec[1].Eta()); 
            // eta-phi distributions
            h_etaphidist_InvMass_under200M[i]->Fill(pi0gammavec[1].Eta(),pi0gammavec[1].Phi()); 
          }
        }
      }

      h_pt1->Fill(photon1.Pt());
      h_pt2->Fill(photon2.Pt());
      h_pTdiff_InvMass->Fill(pi0.Pt(), pi0.M());
      h_InvMass->Fill(pi0.M());
      for(int i=0; i<4; i++){
        pi0smearvec[i]=photon1*(generateRandomNumber()*badcalibsmear[i]+1)+photon2*(generateRandomNumber()*badcalibsmear[i]+1);
        h_InvMass_badcalib_smear[i]->Fill(pi0smearvec[i].M());
      }
      h_mass_eta_lt[lt_eta]->Fill(pi0.M());
    }
  }  // clus1 loop

  /////////////////////////////////////////////////
  //// Truth info
  float wieght = 1;
  PHG4TruthInfoContainer* truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (truthinfo)
  {
    PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
    for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
    {
      // Get truth particle
      const PHG4Particle* truth = iter->second;
      if (!truthinfo->is_primary(truth)) continue;// continue if it is not the primary? turn off for now. and see what secondaries there are.
      TLorentzVector myVector;
      myVector.SetXYZM(truth->get_px(), truth->get_py(), truth->get_pz(), 0.13497);

      float energy = myVector.E();
      h_truth_eta->Fill(myVector.Eta());
      h_truth_e->Fill(energy, wieght);
      h_truth_pt->Fill(myVector.Pt());

      int id =  truth->get_pid();
      h_truth_pid->Fill(id);
      //std::cout << "id=" << id << "   E=" << energy << "  pt=" << myVector.Pt() << "  eta=" << myVector.Eta() << std::endl;
    }
  // try to see secondaries
  PHG4TruthInfoContainer::Range ranges = truthinfo->GetSecondaryParticleRange();
    for (PHG4TruthInfoContainer::ConstIterator iters = ranges.first; iters != ranges.second; ++iters)
    {
      // Get truth particle
      const PHG4Particle* truths = iters->second;
      TLorentzVector myVector;
      myVector.SetXYZM(truth->get_px(), truth->get_py(), truth->get_pz(), 0.13497);
       
      float energy = myVector.E();
      //h_truth_eta->Fill(myVector.Eta());
      //h_truth_e->Fill(energy, wieght);
      //h_truth_pt->Fill(myVector.Pt());

      int id =  truths->get_pid();
      h_truth_pid->Fill(id);

      for(int i=0; i<6; i++){
        if(energy>pidcuts[i]){
          h_truth_pid_cuts[i]->Fill(id);
        }
        else{
          break;
        }
      }
      //std::cout << "id=" << id << "   E=" << energy << "  pt=" << myVector.Pt() << "  eta=" << myVector.Eta() << std::endl;
    }


  }

  ht_phi.clear();
  ht_eta.clear();

  return Fun4AllReturnCodes::EVENT_OK;
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

std::pair<double, double> CaloAna::fitHistogram(TH1F* h)
{
  TF1* fitFunc = new TF1("fitFunc", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x + [5]*x^2 + [6]*x^3", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());

  fitFunc->SetParameter(0, h->GetMaximum());
  fitFunc->SetParameter(1, target_pi0_mass);
  fitFunc->SetParameter(2, 0.01);
  fitFunc->SetParameter(3, 0.0);
  fitFunc->SetParameter(4, 0.0);
  fitFunc->SetParameter(5, 0.0);
  fitFunc->SetParameter(6, 0.0);

  fitFunc->SetParLimits(1, 0.1, 0.2);

  // Perform the fit
  h->Fit("fitFunc", "QN");

  // Get the mean and its error
  double mean = fitFunc->GetParameter(1);
  double errorOnMean = fitFunc->GetParError(1);

  // Create a pair to store the results
  std::pair<double, double> result(mean, errorOnMean);

  return result;
}

void CaloAna::fitEtaSlices(std::string infile, std::string fitOutFile, std::string cdbFile)
{
  TFile* fin = new TFile(infile.c_str());
  cout << "getting hists" << endl;
  TH1F* h_peak_eta = new TH1F("h_peak_eta", "", 96, 0, 96);
  if (!fin) cout << "CaloAna::fitEtaSlices null fin" << endl;
  TH1F* h_M_eta[96];
  for (int i = 0; i < 96; i++) h_M_eta[i] = (TH1F*) fin->Get(Form("h_mass_eta_lt%d", i));

  for (int i = 0; i < 96; i++)
  {
    if (!h_M_eta[i]) cout << "CaloAna::fitEtaSlices null hist" << endl;
    std::pair<double, double> result = fitHistogram(h_M_eta[i]);
    h_peak_eta->SetBinContent(i + 1, result.first);
    h_peak_eta->SetBinError(i + 1, result.second);
  }
  cdbFile = "";

  CDBTTree* cdbttree1 = new CDBTTree(cdbFile.c_str());
  CDBTTree* cdbttree2 = new CDBTTree(cdbFile.c_str());

  string m_fieldname = "Femc_datadriven_qm1_correction";

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
  h_peak_eta->Write();
  for (int i = 0; i < 96; i++)
  {
    h_M_eta[i]->Write();
    delete h_M_eta[i];
  }

  fin->Close();

  cout << "finish fitting suc" << endl;

  return;
}

// Method to generate random numbers
double generateRandomNumber() {
    std::normal_distribution<double> dist(0.0, 1.0);// mean 0, sigma 1
    return dist(rng);
}
