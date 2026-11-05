#include <GlobalVariables.C>

#include <DisplayOn.C>
#include <G4Setup_sPHENIX.C>
#include <G4_Global.C>
#include <G4_ActsGeom.C>

#include <Trkr_RecoInit.C>
#include <Trkr_Clustering.C>
#include <Trkr_LaserClustering.C>

#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllUtils.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <g4tpc/PHG4TpcCentralMembrane.h>

#include <tpccentralmembraneplotter/TpcCentralMembranePlotter.h>

void Fun4All_CentralMembranePlotter(
    const int nEvents = 1000,
    const std::string filelist = "run79507.list")
{
  auto se = Fun4AllServer::instance();
  se->Verbosity(0);

  TRACKING::tpc_zero_supp = true;
  G4TPC::ENABLE_CENTRAL_MEMBRANE_CLUSTERING = true;
  Enable::MVTX_APPLYMISALIGNMENT = true;
  ACTSGEOM::mvtx_applymisalignment = Enable::MVTX_APPLYMISALIGNMENT;

  auto rc = recoConsts::instance();
  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");

  std::ifstream inputList(filelist);
  std::string firstInputFile;
  std::getline(inputList, firstInputFile);
  if (firstInputFile.empty())
  {
    std::cout << "Fun4All_LaserClusterGarfieldPropagation - empty input file list: " << filelist << std::endl;
    gSystem->Exit(1);
  }

  const auto runSegment = Fun4AllUtils::GetRunSegment(firstInputFile);
  const int runnumber = runSegment.first;
  rc->set_IntFlag("RUNNUMBER", runnumber);
  rc->set_uint64Flag("TIMESTAMP", runnumber);

  auto hitsin = new Fun4AllDstInputManager("InputManager");
  hitsin->AddListFile(filelist);
  se->registerInputManager(hitsin);

  const std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");
  auto ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  G4TPC::ENABLE_MODULE_EDGE_CORRECTIONS = true;
  G4TPC::ENABLE_STATIC_CORRECTIONS = false;
  G4TPC::USE_PHI_AS_RAD_STATIC_CORRECTIONS = false;
  G4TPC::ENABLE_AVERAGE_CORRECTIONS = false;
  G4TPC::laser_adc_threshold = 100;

  TrackingInit();
  ACTSGEOM::ActsGeomInit();

  auto centralmembraneplotter = new TpcCentralMembranePlotter;
  centralmembraneplotter->Verbosity(0);
  centralmembraneplotter->setOutputfile("centralmembraneclusters.root");
  se->registerSubsystem(centralmembraneplotter);

  se->run(nEvents);
  se->End();

  gSystem->Exit(0);
}