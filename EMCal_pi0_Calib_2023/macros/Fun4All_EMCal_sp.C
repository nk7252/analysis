#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
//#ifndef FUN4ALL_YEAR1_C
//#define FUN4ALL_YEAR1_C
#include <caloreco/CaloTowerCalib.h>
#include <caloreco/RawClusterBuilderTemplate.h>
#include <caloreco/RawClusterPositionCorrection.h>
//#include <caloreco/RawClusterDeadHotMask.h>
//#include <caloreco/TowerInfoDeadHotMask.h>

#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>

#include <fun4allraw/Fun4AllPrdfInputManager.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/InputFileHandler.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/SubsysReco.h>

#include <caloreco/CaloGeomMapping.h>
#include <caloreco/CaloTowerBuilder.h>
#include <caloreco/CaloTowerCalib.h>
#include <caloreco/CaloWaveformProcessing.h>
#include <caloreco/CaloTowerStatus.h>
#include <caloreco/RawClusterBuilderTemplate.h>
#include <caloreco/RawClusterPositionCorrection.h>

#include <calowaveformsim/CaloWaveformSim.h>

//#include <globalvertex/GlobalVertexReco.h>
//#include <mbd/MbdReco.h>

#include <phool/recoConsts.h>


#include <ffamodules/CDBInterface.h>
//#include <cdbobjects/CDBTTree.h>  // for CDBTTree
#include <GlobalVariables.C>

//#include <litecaloeval/LiteCaloEval.h>

R__LOAD_LIBRARY(libcdbobjects)

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libffamodules.so)
//R__LOAD_LIBRARY(libmbd.so)
//R__LOAD_LIBRARY(libglobalvertex.so)
//R__LOAD_LIBRARY(libLiteCaloEvalTowSlope.so)

#include <caloana/CaloAna.h>
R__LOAD_LIBRARY(libcaloana.so)
#endif
//void createLocalEMCalCalibFile(const string fname, int runNumber);

// void Fun4All_EMCal(int nevents = 0, const std::string &fname = "inputdata.txt",int iter = 0, const std::string &calib_fname="local_calib_copy.root")
void Fun4All_EMCal_sp(int nevents = 10000, const std::string &fname = "inputdata_sp.txt", const std::string &fname_truth = "g4hits.list", const std::string &fnameglobal = "dst_global.list")
{
  //bool enableMasking = 0;

  //bool doFit = 0;
  //bool doHistMake = 0;

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  // se->Verbosity(verbosity);
  recoConsts *rc = recoConsts::instance();

  ifstream file(fname);
  string first_file;
  getline(file, first_file);

  //===============
  // conditions DB flags
  //===============
  pair<int, int> runseg = Fun4AllUtils::GetRunSegment(first_file);
  int runnumber = runseg.first;
  cout << "run number = " << runnumber << endl;

  // global tag
  // rc->set_StringFlag("CDB_GLOBALTAG","PROD_ana.401");//not sure if this is useful for my purposes?
  rc->set_StringFlag("CDB_GLOBALTAG", "MDC2");
  // // 64 bit timestamp
  rc->set_uint64Flag("TIMESTAMP", runnumber);
  cout << "tag, timestamp done" << endl;
  //===============
  // add input managers
  //===============
  Fun4AllInputManager *in = new Fun4AllDstInputManager("DST_TOWERS");
  Fun4AllInputManager *intruth = new Fun4AllDstInputManager("DST_TRUTH");
  Fun4AllInputManager *inglobal = new Fun4AllDstInputManager("DST_GLOBAL");
  cout << "add listfiles to input manager" << endl;

  in->AddListFile(fname,1);
  intruth->AddListFile(fname_truth,1);
  inglobal->AddListFile(fnameglobal,1);
  
  cout << "files added" << endl;

  se->registerInputManager(in);
  se->registerInputManager(intruth);
  se->registerInputManager(inglobal);
  // add this pedestal file for waveform simulation
  Fun4AllInputManager *hitsin = new Fun4AllNoSyncDstInputManager("DST2");
  hitsin->AddFile("pedestal-00046796.root");
  hitsin->Repeat();
  se->registerInputManager(hitsin);

  //* calo waveform simulation
  CEMC_Towers();
  CEMC_Clusters();
  CaloWaveformSim* caloWaveformSim = new CaloWaveformSim();
  caloWaveformSim->set_detector_type(CaloTowerDefs::CEMC);
  caloWaveformSim->set_detector("CEMC");
  caloWaveformSim->set_nsamples(12);
  caloWaveformSim->set_pedestalsamples(12);
  caloWaveformSim->set_timewidth(0.2);
  caloWaveformSim->set_peakpos(6);
  caloWaveformSim->set_calibName("cemc_pi0_twrSlope_v1_default");

  //caloWaveformSim->set_noise_type(CaloWaveformSim::NOISE_NONE);
  //Emma recomended commenting this out. It seems to work regardless of that.
  //caloWaveformSim->get_light_collection_model().load_data_file(string(getenv("CALIBRATIONROOT")) + string("/CEMC/LightCollection/Prototype3Module.xml"),"data_grid_light_guide_efficiency", "data_grid_fiber_trans");
  caloWaveformSim->Verbosity(2);
  se->registerSubsystem(caloWaveformSim);


  CaloTowerBuilder* ca2 = new CaloTowerBuilder();
  ca2->set_detector_type(CaloTowerDefs::CEMC);
  ca2->set_nsamples(12);
  ca2->set_dataflag(false);
  ca2->set_processing_type(CaloWaveformProcessing::TEMPLATE);
  ca2->set_builder_type(CaloTowerDefs::kWaveformTowerv2);
  //match our current ZS threshold ~14ADC for emcal
  //ca2->set_tbt_softwarezerosuppression("/sphenix/user/nkumar/NK_Work_2024/waveformZS/test/src/_tbt_CEMC_zs_x2.root");
  ca2->set_tbt_softwarezerosuppression("/sphenix/user/nkumar/NK_Work_2024/waveformZS/test/src/_tbt_CEMC_zs.root");
  //ca2->set_softwarezerosuppression(true, 14);
  se->registerSubsystem(ca2);

  CaloTowerStatus *statusEMC = new CaloTowerStatus("CEMCSTATUS");
  statusEMC->set_detector_type(CaloTowerDefs::CEMC);
  statusEMC->set_time_cut(1);
  se->registerSubsystem(statusEMC);

  std::cout << "Calibrating EMCal" << std::endl;
  CaloTowerCalib *calibEMC = new CaloTowerCalib("CEMCCALIB");
  calibEMC->set_detector_type(CaloTowerDefs::CEMC);
  calibEMC->set_outputNodePrefix("TOWERINFO_CALIB_");
  se->registerSubsystem(calibEMC);

  std::cout << "Building clusters" << std::endl;
  RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
  ClusterBuilder->Detector("CEMC");
  ClusterBuilder->set_threshold_energy(0.030);  // for when using basic calibration
  std::string emc_prof = getenv("CALIBRATIONROOT");
  emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
  ClusterBuilder->LoadProfile(emc_prof);
  ClusterBuilder->set_UseTowerInfo(1);  // to use towerinfo objects rather than old RawTower
  se->registerSubsystem(ClusterBuilder);
  //*/

  cout << "input manager registered" << endl;

  // Fun4AllInputManager *in2 = new Fun4AllDstInputManager("DST_TOWERS2");
  // in2->AddListFile(fnamehits);
  // se->registerInputManager(in2);

  std::string filename = first_file.substr(first_file.find_last_of("/\\") + 1);
  std::string OutFile = Form("OUTHIST_iter_%s", filename.c_str());
  // std::string OutFile = Form("OUTHIST_iter%d_%s",iter , filename.c_str());

  /*
  if (iter == 0)
  {
    cout << "creating emcal calib" << endl;
    createLocalEMCalCalibFile(calib_fname.c_str(), runnumber);
    cout << "creating " << calib_fname.c_str() << " and exiting" << endl;
    return;
  }

  ////////////////////
  // Calibrate towers

  std::cout << "Calibrating EMCal" << std::endl;
  CaloTowerCalib *calibEMC = new CaloTowerCalib("CEMCCALIB");
  calibEMC->set_detector_type(CaloTowerDefs::CEMC);
  calibEMC->set_directURL(calib_fname.c_str());
  se->registerSubsystem(calibEMC);


  //////////////////
  // Clusters
  std::cout << "Building clusters" << std::endl;
  RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
  ClusterBuilder->Detector("CEMC");
  ClusterBuilder->set_threshold_energy(0.20);  // for when using basic calibration
  std::string emc_prof = getenv("CALIBRATIONROOT");
  emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
  ClusterBuilder->LoadProfile(emc_prof);
  ClusterBuilder->set_UseTowerInfo(1);  // to use towerinfo objects rather than old RawTower
  ClusterBuilder->setOutputClusterNodeName("CLUSTERINFO_CEMC2");
  se->registerSubsystem(ClusterBuilder);
  */
  /*
    std::cout << "Applying Position Dependent Correction" << std::endl;
    RawClusterPositionCorrection *clusterCorrection = new RawClusterPositionCorrection("CEMC");
    clusterCorrection->set_UseTowerInfo(1);  // to use towerinfo objects rather than old RawTower
    se->registerSubsystem(clusterCorrection);
  */

  ///////////////////
  // analysis modules
  /*
  if (iter==1){
  LiteCaloEval *eval7e = new LiteCaloEval("CEMCEVALUATOR2", "CEMC",OutFile);
  eval7e->CaloType(LiteCaloEval::CEMC);
  eval7e->setInputTowerNodeName("TOWERINFO_CALIB_CEMC");
  se->registerSubsystem(eval7e);
  }
  */

  // if (iter>1){
  CaloAna *ca = new CaloAna("calomodulename", OutFile);
  
  ca->set_timing_cut_width(16);
  ca->apply_vertex_cut(false);
  ca->set_vertex_cut(20.);
  se->registerSubsystem(ca);
  std::cout << "Subsystems registered" << std::endl;
  //}

  se->run(nevents);
  se->End();
  se->PrintTimer();
  delete se;
  std::cout << "All done!" << std::endl;

  TFile *f_done_signal = new TFile("DONE.root", "recreate");
  std::cout << "All done!" << std::endl;
  gSystem->Exit(0);
}

/*
void createLocalEMCalCalibFile(const string fname, int runNumber)
{
  string default_time_independent_calib = "cemc_pi0_twrSlope_v1_default";
  string m_calibName = "cemc_pi0_twrSlope_v1";
  std::cout << "calib:GET URL" << std::endl;
  string calibdir = CDBInterface::instance()->getUrl(m_calibName);
  string filePath;
  std::cout << "calib:check calibdir" << std::endl;
  if (!calibdir.empty())
  {
    std::cout << "calib:calibdir not empty" << std::endl;
    filePath = calibdir;
    // cdbttree = new CDBTTree(calibdir);
  }
  else
  {
    std::cout << "calib:Calibdir maybe empty" << std::endl;
    calibdir = CDBInterface::instance()->getUrl(default_time_independent_calib);

    if (calibdir.empty())
    {
      std::cout << "No EMCal Calibration NOT even a default" << std::endl;
      exit(1);
    }
    filePath = calibdir;
    // cdbttree = new CDBTTree(calibdir);
    std::cout << "No specific file for " << m_calibName << " found, using default calib " << default_time_independent_calib << std::endl;
  }

  TFile *f_cdb = new TFile(filePath.c_str());
  f_cdb->Cp(fname.c_str());

  std::cout << "created local Calib file for run " << runNumber << " named " << fname << std::endl;

  delete f_cdb;
}
//#endif
*/

