
#ifndef MACRO_FUN4ALL_READDST_C
#define MACRO_FUN4ALL_READDST_C

#include <tpcrawtottree/TPCRawDataTree.h>

#include <fun4allraw/Fun4AllPrdfInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libTPCRawDataTree.so)

int Fun4All_TPC_UnpackPRDF(const int nEvents = 100,
                           const string &inputFile = "/sphenix/lustre01/sphnxpro/rawdata/commissioning/tpc/beam/TPC_ebdc*_beam-00011012-0000.prdf"  //
)
{
  //---------------
  // Fun4All server
  //---------------
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(1);

  string outDir = "/sphenix/u/jamesj3j3/tpc/sPHENIXProjects/beam-run-11012/";

  string fileName = inputFile;
  size_t pos = fileName.find("TPC_ebdc");
  fileName.erase(fileName.begin(),fileName.begin()+pos);
  
  TPCRawDataTree *r2tree = new TPCRawDataTree(outDir + fileName + "_TPCRawDataTree_skip100.root");/////////////////////////////

  // add all possible TPC packet that we need to analyze
  for (int packet = 4000; packet<=4230; packet+=10)
  {
    r2tree->AddPacket(packet);
    r2tree->AddPacket(packet+1);
  }

  se->registerSubsystem(r2tree);

  Fun4AllPrdfInputManager *in1 = new Fun4AllPrdfInputManager("PRDF1");
  in1->AddFile(inputFile);
  se->registerInputManager(in1);

  se->skip(100);/////////////////////////////
  se->run(nEvents);

  se->End();

  delete se;
  std::cout << "All done processing" << std::endl;
  gSystem->Exit(0);
  return 0;
}
#endif
