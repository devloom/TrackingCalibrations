#ifndef MACRO_TRKRLASERCLUS_C
#define MACRO_TRKRLASERCLUS_C

#include <GlobalVariables.C>

#include <G4_ActsGeom.C>
#include <G4_TrkrVariables.C>

#include <tpc/LaserClusterizer.h>

#include <tpccalib/TpcDirectLaserReconstruction.h>
#include <tpccalib/TpcLaminationFitting.h>

#include <fun4all/Fun4AllServer.h>

R__LOAD_LIBRARY(libtpc.so)
R__LOAD_LIBRARY(libtpccalib.so)

void TPC_LaserClustering()
{
  // central membrane reconstruction
  if (G4TPC::ENABLE_CENTRAL_MEMBRANE_CLUSTERING)
  {
    int verbosity = std::max(Enable::VERBOSITY, Enable::TPC_VERBOSITY);
    ACTSGEOM::ActsGeomInit();
    Fun4AllServer *se = Fun4AllServer::instance();

    LaserClusterizer *laserClusterizer = new LaserClusterizer;
    laserClusterizer->Verbosity(verbosity);
    laserClusterizer->set_max_time_samples(TRACKING::reco_tpc_maxtime_sample);
    laserClusterizer->set_adc_threshold(G4TPC::laser_adc_threshold);
    laserClusterizer->set_do_sequential(G4TPC::LaserClusteringSequential);
    laserClusterizer->set_do_fitting(G4TPC::laserClusterFitting);
    se->registerSubsystem(laserClusterizer);
  }

  // direct laser reconstruction
  if (G4TPC::ENABLE_DIRECT_LASER_HITS)
  {
    Fun4AllServer *se = Fun4AllServer::instance();
    auto *directLaserReconstruction = new TpcDirectLaserReconstruction;
    directLaserReconstruction->set_outputfile(G4TPC::DIRECT_LASER_ROOTOUTPUT_FILENAME);
    directLaserReconstruction->set_savehistograms(G4TPC::DIRECT_LASER_SAVEHISTOGRAMS);
    directLaserReconstruction->set_histogram_outputfile(G4TPC::DIRECT_LASER_HISTOGRAMOUTPUT_FILENAME);
    se->registerSubsystem(directLaserReconstruction);
  }
}

void TPC_LaminationClustering()
{
  // central membrane reconstruction
  if (G4TPC::ENABLE_CENTRAL_MEMBRANE_CLUSTERING)
  {
    int verbosity = std::max(Enable::VERBOSITY, Enable::TPC_VERBOSITY);
    ACTSGEOM::ActsGeomInit();
    Fun4AllServer *se = Fun4AllServer::instance();

    LaserClusterizer *laminationClusterizer = new LaserClusterizer;
    laminationClusterizer->Verbosity(verbosity);
    laminationClusterizer->set_max_time_samples(TRACKING::reco_tpc_maxtime_sample);
    laminationClusterizer->set_adc_threshold(G4TPC::laser_adc_threshold);
    laminationClusterizer->set_do_sequential(G4TPC::LaserClusteringSequential);
    laminationClusterizer->set_do_fitting(G4TPC::laserClusterFitting);
    laminationClusterizer->set_lamination(true);
    se->registerSubsystem(laminationClusterizer);
  }
}

void TPC_LaminationFitting()
{
  if (G4TPC::ENABLE_CENTRAL_MEMBRANE_CLUSTERING)
  {
    int verbosity = std::max(Enable::VERBOSITY, Enable::TPC_VERBOSITY);
    ACTSGEOM::ActsGeomInit();
    Fun4AllServer *se = Fun4AllServer::instance();

    TpcLaminationFitting *laminations = new TpcLaminationFitting;
    laminations->Verbosity(verbosity);
    // laminations->set_nLayerCut(1);
    laminations->set_nLayerCut(2);
    laminations->set_fieldOff(G4TPC::BFieldOff);
    // laminations->set_stripePatternFile("/sphenix/u/bkimelman/CMStripePattern_Feb6_2026.root");
  
    // use ideal stripe pattern
    laminations->set_stripePatternFile("/sphenix/u/bkimelman/CMStripePattern_ideal.root");

    // use stripe pattern from field off beam off data
    // laminations->set_stripePatternFile("/sphenix/user/dloomis/DistortionsAnalysis/stripePatterns/CMStripePattern_Mar25_2026_alignedcoordinates.root");
    // laminations->set_stripePatternFile("/sphenix/user/dloomis/DistortionsAnalysis/stripePatterns/CMStripePattern_Apr15_2026_idealcoordinates_78242fieldoffbeamoff_updatedstripefits.root");

    // use stripe pattern from field on beam off data with static corrections    
    // laminations->set_stripePatternFile("/sphenix/user/dloomis/DistortionsAnalysis/stripePatterns/CMStripePattern_Apr10_2026_idealcoordinates_83319fieldonbeamoff_wstaticcorrections.root");
    // laminations->set_stripePatternFile("/sphenix/user/dloomis/DistortionsAnalysis/stripePatterns/CMStripePattern_Apr11_2026_alignedcoordinates_83313fieldonbeamoff_wstaticcorrections.root");
    
    // laminations->set_stripePatternFile("/sphenix/user/dloomis/DistortionsAnalysis/stripePatterns/CMStripePattern_Apr13_2026_idealcoordinates_83319fieldonbeamoff_wstaticcorrections_updatedstaticmap.root");
    // laminations->set_stripePatternFile("/sphenix/user/dloomis/DistortionsAnalysis/stripePatterns/CMStripePattern_Apr14_2026_idealcoordinates_83319fieldonbeamoff_wstaticcorrections_updatedstaticmap_updatedstripefits.root");

    // use stripe pattern from field on beam off data without static corrections    
    // laminations->set_stripePatternFile("/sphenix/user/dloomis/DistortionsAnalysis/stripePatterns/CMStripePattern_Mar31_2026_alignedcoordinates_75074fieldonbeamoff.root");



    // use stripe pattern from field on beam on data with static corrections
    // use stripe pattern from field on beam off data with static corrections
    laminations->set_saveAllLaminationHistograms(true);
    // laminations->set_ppMode(TRACKING::pp_mode);
    laminations->set_ppMode(TRACKING::streaming_mode);
    laminations->setOutputfile(G4TPC::LaminationOutputName);
    if (!G4TPC::LaminationQAName.empty())
    {
      laminations->set_QAFileName(G4TPC::LaminationQAName);
    }
    se->registerSubsystem(laminations);
  }
}

#endif
