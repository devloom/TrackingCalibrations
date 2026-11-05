#include "TpcLaserClusterGarfieldPropagator.h"

#include <centralmembraneanalysisv2/StripeComparison.h>
#include <centralmembraneanalysisv2/StripeDetector.h>
#include <centralmembraneanalysisv2/helpers.h>
#include <centralmembraneanalysisv2/StripeMatchingTypes.h>
#include <centralmembraneanalysisv2/parameters.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phgarfield/PHGarfield.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <trackbase/LaserCluster.h>
#include <trackbase/LaserClusterContainer.h>
#include <trackbase/TpcDefs.h>

#include <TAxis.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH2D.h>
#include <TObject.h>
#include <TPolyLine3D.h>
#include <TGraph.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

static TH2* LoadDetachedTH2(TFile* inputFile, const std::string& histogramName, const std::string& cloneName)
{
  if (!inputFile || inputFile->IsZombie())
  {
    return nullptr;
  }

  auto* histogram = dynamic_cast<TH2*>(inputFile->Get(histogramName.c_str()));
  if (!histogram)
  {
    std::cout << "LoadDetachedTH2 - missing histogram: " << histogramName << std::endl;
    return nullptr;
  }

  auto* clone = dynamic_cast<TH2*>(histogram->Clone(cloneName.c_str()));
  if (clone)
  {
    clone->SetDirectory(nullptr);
    clone->SetStats(false);
  }
  return clone;
}

static TH2D* BuildXYStripeHistogram(const TH2* phiRHistogram, const std::string& name, const std::string& title)
{
  if (!phiRHistogram)
  {
    return nullptr;
  }

  const double rMax = std::max(std::abs(phiRHistogram->GetYaxis()->GetXmin()), std::abs(phiRHistogram->GetYaxis()->GetXmax()));
  auto* xyHistogram = new TH2D(name.c_str(), title.c_str(),
                               2000, -rMax, rMax,
                               2000, -rMax, rMax);
  xyHistogram->SetDirectory(nullptr);
  xyHistogram->SetStats(false);

  for (int phiBin = 1; phiBin <= phiRHistogram->GetNbinsX(); ++phiBin)
  {
    const double phi = phiRHistogram->GetXaxis()->GetBinCenter(phiBin);
    for (int rBin = 1; rBin <= phiRHistogram->GetNbinsY(); ++rBin)
    {
      const double content = phiRHistogram->GetBinContent(phiBin, rBin);
      if (content == 0.0)
      {
        continue;
      }

      const double r = phiRHistogram->GetYaxis()->GetBinCenter(rBin);
      xyHistogram->Fill(r * std::cos(phi), r * std::sin(phi), content);
    }
  }

  return xyHistogram;
}

static TH2* BuildCleanedStripeHistogram(TH2* histogram, const std::string& cleanedName, std::vector<double>& radialGapCenters)
{
  radialGapCenters.clear();
  if (!histogram)
  {
    return nullptr;
  }

  TH2* cleanedHistogram = nullptr;
  if (!StripeDetectionParameters::applyPreStripeLaminationMask)
  {
    cleanedHistogram = dynamic_cast<TH2*>(histogram->Clone(cleanedName.c_str()));
  }
  else
  {
    auto cleaning = StripeMatchingHelpers::clean_laminations(histogram);
    cleanedHistogram = cleaning.cleanedHistogram;
    radialGapCenters = std::move(cleaning.radialGapCenters);
    if (cleanedHistogram)
    {
      cleanedHistogram->SetName(cleanedName.c_str());
    }
  }

  if (cleanedHistogram)
  {
    cleanedHistogram->SetDirectory(nullptr);
    cleanedHistogram->SetStats(false);
  }
  return cleanedHistogram;
}

static bool DetectStripesWithMatchingCode(TH2* histogram, const std::string& cleanedName, std::vector<StripeCluster>& stripes, std::vector<double>& radialGapCenters, TFile& outputFile)
{
  stripes.clear();
  TH2* cleanedHistogram = BuildCleanedStripeHistogram(histogram, cleanedName, radialGapCenters);
  if (!cleanedHistogram)
  {
    std::cout << "DetectStripesWithMatchingCode - failed to build cleaned histogram " << cleanedName << std::endl;
    return false;
  }

  outputFile.cd();
  cleanedHistogram->Write("", TObject::kOverwrite);

  StripeDetector detector;
  const bool ok = detector.detect(cleanedHistogram, stripes);
  delete cleanedHistogram;
  return ok;
}

static TGraph* BuildStripeGraph(const std::vector<StripeCluster>& stripes,
                                  const std::string& name, const std::string& title,
                                  bool xy)
{
  auto* graph = new TGraph(static_cast<int>(stripes.size()));
  graph->SetName(name.c_str());
  graph->SetTitle(title.c_str());
  for (std::size_t i = 0; i < stripes.size(); ++i)
  {
    if (xy)
    {
      graph->SetPoint(static_cast<int>(i), stripes[i].r * std::cos(stripes[i].phi), stripes[i].r * std::sin(stripes[i].phi));
    }
    else
    {
      graph->SetPoint(static_cast<int>(i), stripes[i].phi, stripes[i].r);
    }
  }
  return graph;
}



TpcLaserClusterGarfieldPropagator::TpcLaserClusterGarfieldPropagator(const std::string& name)
  : SubsysReco(name)
{
}

TpcLaserClusterGarfieldPropagator::~TpcLaserClusterGarfieldPropagator()
{
  for (int side = 0; side < 2; ++side)
  {
    delete m_hRaw[side];
    delete m_hGarfield[side];
    delete m_hRawXY[side];
    delete m_hGarfieldXY[side];
    delete m_hDeltaR[side];
    delete m_hDeltaPhi[side];
    delete m_hDeltaRInterpolated[side];
    delete m_hDeltaPhiInterpolated[side];
    delete m_hDeltaRInterpolatedXY[side];
    delete m_hDeltaPhiInterpolatedXY[side];
  }
  delete m_hDriftStats;
}

void TpcLaserClusterGarfieldPropagator::setOutputfile(const std::string& outputfile)
{
  m_outputfile = outputfile;
}

void TpcLaserClusterGarfieldPropagator::setPHGarfield(PHGarfield* phgarfield)
{
  m_phgarfield = phgarfield;
}

void TpcLaserClusterGarfieldPropagator::setZLaunch(double zAbsCm)
{
  m_zAbsLaunchCm = std::abs(zAbsCm);
}

void TpcLaserClusterGarfieldPropagator::setHistogramBinning(int nPhi, int nR, double rMin, double rMax)
{
  if (nPhi <= 0 || nR <= 0 || rMin >= rMax)
  {
    std::cout << "TpcLaserClusterGarfieldPropagator::setHistogramBinning - invalid binning requested" << std::endl;
    return;
  }

  m_phiBins = nPhi;
  m_rBins = nR;
  m_rMin = rMin;
  m_rMax = rMax;
}

void TpcLaserClusterGarfieldPropagator::setDeltaMapBinning(int nPhi, int nR)
{
  if (nPhi <= 0 || nR <= 0)
  {
    std::cout << "TpcLaserClusterGarfieldPropagator::setDeltaMapBinning - invalid binning requested" << std::endl;
    return;
  }

  m_deltaMapPhiBins = nPhi;
  m_deltaMapRBins = nR;
}

void TpcLaserClusterGarfieldPropagator::setFieldOffBeamOffStripeFile(const std::string& filename)
{
  m_fieldOffBeamOffStripeFile = filename;
}

void TpcLaserClusterGarfieldPropagator::setStripeComparisonBranchWeights(double residualWeight, double unmatchedWeight, double cleaningWeight)
{
  m_branchResidualWeight = residualWeight;
  m_branchUnmatchedWeight = unmatchedWeight;
  m_branchCleaningWeight = cleaningWeight;
}

int TpcLaserClusterGarfieldPropagator::InitRun(PHCompositeNode* /*topNode*/)
{
  if (!m_phgarfield)
  {
    std::cout << "TpcLaserClusterGarfieldPropagator::InitRun - PHGarfield pointer is null" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  CreateHistograms();
  for (int side = 0; side < 2; ++side)
  {
    if (!m_hRaw[side] || !m_hGarfield[side] || !m_hRawXY[side] || !m_hGarfieldXY[side] || !m_hDeltaR[side] || !m_hDeltaPhi[side] || !m_hDeltaRInterpolated[side] || !m_hDeltaPhiInterpolated[side] || !m_hDeltaRInterpolatedXY[side] || !m_hDeltaPhiInterpolatedXY[side])
    {
      std::cout << "TpcLaserClusterGarfieldPropagator::InitRun - failed to create histograms" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcLaserClusterGarfieldPropagator::process_event(PHCompositeNode* topNode)
{
  const int nodeStatus = GetNodes(topNode);
  if (nodeStatus != Fun4AllReturnCodes::EVENT_OK)
  {
    return nodeStatus;
  }

  ++m_processedEvents;
  const auto clusterRange = m_laserClusterContainer->getClusters();
  for (auto clusterIter = clusterRange.first; clusterIter != clusterRange.second; ++clusterIter)
  {
    const auto& clusterKey = clusterIter->first;
    const auto* cluster = clusterIter->second;
    ++m_seenClusters;

    if (!cluster)
    {
      continue;
    }

    if (cluster->getNLayers() < 2 || cluster->getSDWeightedLayer() >= 0.5)
    {
      continue;
    }

    const int side = TpcDefs::getSide(clusterKey) ? 1 : 0;
    const double x = cluster->getX();
    const double y = cluster->getY();
    const double phi = NormalizePhi(std::atan2(y, x));
    const double r = std::hypot(x, y);

    m_hRaw[side]->Fill(phi, r);
    m_hRawXY[side]->Fill(x, y);
    ++m_acceptedClusters[side];
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcLaserClusterGarfieldPropagator::End(PHCompositeNode* /*topNode*/)
{
  if (!m_phgarfield)
  {
    std::cout << "TpcLaserClusterGarfieldPropagator::End - PHGarfield pointer is null" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PropagateOccupiedBins();
  FillFullDeltaMaps();
  WriteOutput();

  std::cout << "TpcLaserClusterGarfieldPropagator processed " << m_processedEvents
            << " events, accepted clusters negz=" << m_acceptedClusters[0]
            << " posz=" << m_acceptedClusters[1] << std::endl;
  std::cout << "Garfield drifted occupied bins attempted=" << m_attemptedBins
            << " success=" << m_successfulBins
            << " null_paths=" << m_nullPaths
            << " short_paths=" << m_shortPaths
            << " no_cm_point=" << m_noCentralMembranePoint
            << " out_of_range=" << m_outOfRangePropagated << std::endl;
  std::cout << "Garfield full delta maps attempted=" << m_fullMapAttemptedBins
            << " success=" << m_fullMapSuccessfulBins << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcLaserClusterGarfieldPropagator::GetNodes(PHCompositeNode* topNode)
{
  m_laserClusterContainer = findNode::getClass<LaserClusterContainer>(topNode, "LASER_CLUSTER");
  if (!m_laserClusterContainer)
  {
    std::cout << "TpcLaserClusterGarfieldPropagator::GetNodes - LASER_CLUSTER node missing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void TpcLaserClusterGarfieldPropagator::CreateHistograms()
{
  const char* sideNames[2] = {"negz", "posz"};
  for (int side = 0; side < 2; ++side)
  {
    delete m_hRaw[side];
    delete m_hGarfield[side];
    delete m_hRawXY[side];
    delete m_hGarfieldXY[side];
    delete m_hDeltaR[side];
    delete m_hDeltaPhi[side];
    delete m_hDeltaRInterpolated[side];
    delete m_hDeltaPhiInterpolated[side];
    delete m_hDeltaRInterpolatedXY[side];
    delete m_hDeltaPhiInterpolatedXY[side];

    const std::string suffix = sideNames[side];
    m_hRaw[side] = new TH2D(("hLaserCluster_raw_" + suffix).c_str(),
                            ("Raw laser clusters - " + suffix + ";#phi [rad];R [cm]").c_str(),
                            m_phiBins, 0.0, 2.0 * M_PI, m_rBins, m_rMin, m_rMax);
    m_hGarfield[side] = new TH2D(("hLaserCluster_garfield_" + suffix).c_str(),
                                 ("Laser clusters reverse-drifted with PHGarfield to central membrane - " + suffix + ";#phi [rad];R [cm]").c_str(),
                                 m_phiBins, 0.0, 2.0 * M_PI, m_rBins, m_rMin, m_rMax);
    m_hRawXY[side] = new TH2D(("hLaserCluster_raw_xy_" + suffix).c_str(),
                              ("Raw laser clusters - " + suffix + ";x [cm];y [cm]").c_str(),
                              2000, -m_rMax, m_rMax, 2000, -m_rMax, m_rMax);
    m_hGarfieldXY[side] = new TH2D(("hLaserCluster_garfield_xy_" + suffix).c_str(),
                                   ("Laser clusters reverse-drifted with PHGarfield to central membrane - " + suffix + ";x [cm];y [cm]").c_str(),
                                   2000, -m_rMax, m_rMax, 2000, -m_rMax, m_rMax);
    m_hDeltaR[side] = new TH2D(("hLaserCluster_garfield_deltaR_" + suffix).c_str(),
                               ("PHGarfield readout minus central-membrane R - " + suffix + ";#phi_{readout} [rad];R_{readout} [cm];#DeltaR [cm]").c_str(),
                               m_phiBins, 0.0, 2.0 * M_PI, m_rBins, m_rMin, m_rMax);
    m_hDeltaPhi[side] = new TH2D(("hLaserCluster_garfield_deltaPhi_" + suffix).c_str(),
                                 ("PHGarfield readout minus central-membrane #phi - " + suffix + ";#phi_{readout} [rad];R_{readout} [cm];#Delta#phi [rad]").c_str(),
                                 m_phiBins, 0.0, 2.0 * M_PI, m_rBins, m_rMin, m_rMax);
    m_hDeltaRInterpolated[side] = new TH2D(("hLaserCluster_garfield_deltaR_interpolated_" + suffix).c_str(),
                                           ("PHGarfield full-grid interpolated readout minus central-membrane R - " + suffix + ";#phi_{readout} [rad];R_{readout} [cm];#DeltaR [cm]").c_str(),
                                           m_deltaMapPhiBins, 0.0, 2.0 * M_PI, m_deltaMapRBins, m_rMin, m_rMax);
    m_hDeltaPhiInterpolated[side] = new TH2D(("hLaserCluster_garfield_deltaPhi_interpolated_" + suffix).c_str(),
                                             ("PHGarfield full-grid interpolated readout minus central-membrane #phi - " + suffix + ";#phi_{readout} [rad];R_{readout} [cm];#Delta#phi [rad]").c_str(),
                                             m_deltaMapPhiBins, 0.0, 2.0 * M_PI, m_deltaMapRBins, m_rMin, m_rMax);
    m_hDeltaRInterpolatedXY[side] = new TH2D(("hLaserCluster_garfield_deltaR_interpolated_xy_" + suffix).c_str(),
                                             ("PHGarfield full-grid interpolated readout minus central-membrane R - " + suffix + ";x_{readout} [cm];y_{readout} [cm];#DeltaR [cm]").c_str(),
                                             m_deltaMapRBins, -m_rMax, m_rMax, m_deltaMapRBins, -m_rMax, m_rMax);
    m_hDeltaPhiInterpolatedXY[side] = new TH2D(("hLaserCluster_garfield_deltaPhi_interpolated_xy_" + suffix).c_str(),
                                               ("PHGarfield full-grid interpolated readout minus central-membrane #phi - " + suffix + ";x_{readout} [cm];y_{readout} [cm];#Delta#phi [rad]").c_str(),
                                               m_deltaMapRBins, -m_rMax, m_rMax, m_deltaMapRBins, -m_rMax, m_rMax);

    m_hRaw[side]->SetDirectory(nullptr);
    m_hGarfield[side]->SetDirectory(nullptr);
    m_hRawXY[side]->SetDirectory(nullptr);
    m_hGarfieldXY[side]->SetDirectory(nullptr);
    m_hDeltaR[side]->SetDirectory(nullptr);
    m_hDeltaPhi[side]->SetDirectory(nullptr);
    m_hDeltaRInterpolated[side]->SetDirectory(nullptr);
    m_hDeltaPhiInterpolated[side]->SetDirectory(nullptr);
    m_hDeltaRInterpolatedXY[side]->SetDirectory(nullptr);
    m_hDeltaPhiInterpolatedXY[side]->SetDirectory(nullptr);
    m_hRaw[side]->SetStats(false);
    m_hGarfield[side]->SetStats(false);
    m_hRawXY[side]->SetStats(false);
    m_hGarfieldXY[side]->SetStats(false);
    m_hDeltaRInterpolated[side]->SetStats(false);
    m_hDeltaPhiInterpolated[side]->SetStats(false);
    m_hDeltaRInterpolatedXY[side]->SetStats(false);
    m_hDeltaPhiInterpolatedXY[side]->SetStats(false);
  }

  delete m_hDriftStats;
  m_hDriftStats = new TH1D("hGarfieldDriftStats", "PHGarfield drift counters", 12, 0.5, 12.5);
  m_hDriftStats->SetDirectory(nullptr);
  m_hDriftStats->GetXaxis()->SetBinLabel(1, "processed_events");
  m_hDriftStats->GetXaxis()->SetBinLabel(2, "seen_clusters");
  m_hDriftStats->GetXaxis()->SetBinLabel(3, "accepted_negz");
  m_hDriftStats->GetXaxis()->SetBinLabel(4, "accepted_posz");
  m_hDriftStats->GetXaxis()->SetBinLabel(5, "attempted_bins");
  m_hDriftStats->GetXaxis()->SetBinLabel(6, "successful_bins");
  m_hDriftStats->GetXaxis()->SetBinLabel(7, "null_paths");
  m_hDriftStats->GetXaxis()->SetBinLabel(8, "short_paths");
  m_hDriftStats->GetXaxis()->SetBinLabel(9, "no_cm_point");
  m_hDriftStats->GetXaxis()->SetBinLabel(10, "out_of_range");
  m_hDriftStats->GetXaxis()->SetBinLabel(11, "full_map_attempted");
  m_hDriftStats->GetXaxis()->SetBinLabel(12, "full_map_success");
}

void TpcLaserClusterGarfieldPropagator::PropagateOccupiedBins()
{
  for (int side = 0; side < 2; ++side)
  {
    if (!m_hRaw[side])
    {
      continue;
    }

    for (int phiBin = 1; phiBin <= m_hRaw[side]->GetNbinsX(); ++phiBin)
    {
      for (int rBin = 1; rBin <= m_hRaw[side]->GetNbinsY(); ++rBin)
      {
        const double weight = m_hRaw[side]->GetBinContent(phiBin, rBin);
        if (weight <= 0.0)
        {
          continue;
        }

        PropagateBin(side, phiBin, rBin, weight);
      }
    }
  }
}

void TpcLaserClusterGarfieldPropagator::FillFullDeltaMaps()
{
  for (int side = 0; side < 2; ++side)
  {
    if (!m_hDeltaRInterpolated[side] || !m_hDeltaPhiInterpolated[side] || !m_hDeltaRInterpolatedXY[side] || !m_hDeltaPhiInterpolatedXY[side])
    {
      continue;
    }

    for (int phiBin = 1; phiBin <= m_hDeltaRInterpolated[side]->GetNbinsX(); ++phiBin)
    {
      for (int rBin = 1; rBin <= m_hDeltaRInterpolated[side]->GetNbinsY(); ++rBin)
      {
        FillFullDeltaMapBin(side, phiBin, rBin);
      }
    }

    for (int xBin = 1; xBin <= m_hDeltaRInterpolatedXY[side]->GetNbinsX(); ++xBin)
    {
      for (int yBin = 1; yBin <= m_hDeltaRInterpolatedXY[side]->GetNbinsY(); ++yBin)
      {
        FillFullDeltaMapXYBin(side, xBin, yBin);
      }
    }
  }
}

bool TpcLaserClusterGarfieldPropagator::PropagateBin(int side, int phiBin, int rBin, double weight)
{
  if (side < 0 || side >= 2 || weight <= 0.0)
  {
    return false;
  }

  ++m_attemptedBins;

  const double phiReadout = m_hRaw[side]->GetXaxis()->GetBinCenter(phiBin);
  const double rReadout = m_hRaw[side]->GetYaxis()->GetBinCenter(rBin);
  const double xReadout = rReadout * std::cos(phiReadout);
  const double yReadout = rReadout * std::sin(phiReadout);
  const double zLaunch = side == 1 ? m_zAbsLaunchCm : -m_zAbsLaunchCm;

  std::unique_ptr<TPolyLine3D> path(m_phgarfield->ReverseDrift(xReadout, yReadout, zLaunch));
  if (!path)
  {
    ++m_nullPaths;
    return false;
  }
  if (path->GetN() < 2)
  {
    ++m_shortPaths;
    return false;
  }

  double rCm = 0.0;
  double phiCm = 0.0;
  if (!InterpolateOrClosestAtZ(path.get(), 0.0, rCm, phiCm))
  {
    ++m_noCentralMembranePoint;
    return false;
  }

  const int phiCmBin = m_hGarfield[side]->GetXaxis()->FindBin(phiCm);
  const int rCmBin = m_hGarfield[side]->GetYaxis()->FindBin(rCm);
  if (phiCmBin < 1 || phiCmBin > m_phiBins || rCmBin < 1 || rCmBin > m_rBins)
  {
    ++m_outOfRangePropagated;
    return false;
  }

  double deltaPhi = phiReadout - phiCm;
  if (deltaPhi > M_PI)
  {
    deltaPhi -= 2.0 * M_PI;
  }
  else if (deltaPhi < -M_PI)
  {
    deltaPhi += 2.0 * M_PI;
  }

  m_hGarfield[side]->Fill(phiCm, rCm, weight);
  m_hGarfieldXY[side]->Fill(rCm * std::cos(phiCm), rCm * std::sin(phiCm), weight);
  m_hDeltaR[side]->SetBinContent(phiBin, rBin, rReadout - rCm);
  m_hDeltaPhi[side]->SetBinContent(phiBin, rBin, deltaPhi);
  ++m_successfulBins;

  return true;
}

bool TpcLaserClusterGarfieldPropagator::FillFullDeltaMapBin(int side, int phiBin, int rBin)
{
  if (side < 0 || side >= 2 || !m_hDeltaRInterpolated[side] || !m_hDeltaPhiInterpolated[side])
  {
    return false;
  }

  ++m_fullMapAttemptedBins;

  const double phiReadout = m_hDeltaRInterpolated[side]->GetXaxis()->GetBinCenter(phiBin);
  const double rReadout = m_hDeltaRInterpolated[side]->GetYaxis()->GetBinCenter(rBin);
  const double xReadout = rReadout * std::cos(phiReadout);
  const double yReadout = rReadout * std::sin(phiReadout);
  const double zLaunch = side == 1 ? m_zAbsLaunchCm : -m_zAbsLaunchCm;

  std::unique_ptr<TPolyLine3D> path(m_phgarfield->ReverseDrift(xReadout, yReadout, zLaunch));
  if (!path || path->GetN() < 2)
  {
    return false;
  }

  double rCm = 0.0;
  double phiCm = 0.0;
  if (!InterpolateOrClosestAtZ(path.get(), 0.0, rCm, phiCm))
  {
    return false;
  }

  double deltaPhi = phiReadout - phiCm;
  if (deltaPhi > M_PI)
  {
    deltaPhi -= 2.0 * M_PI;
  }
  else if (deltaPhi < -M_PI)
  {
    deltaPhi += 2.0 * M_PI;
  }

  const double deltaR = rReadout - rCm;
  m_hDeltaRInterpolated[side]->SetBinContent(phiBin, rBin, deltaR);
  m_hDeltaPhiInterpolated[side]->SetBinContent(phiBin, rBin, deltaPhi);
  ++m_fullMapSuccessfulBins;

  return true;
}

bool TpcLaserClusterGarfieldPropagator::FillFullDeltaMapXYBin(int side, int xBin, int yBin)
{
  if (side < 0 || side >= 2 || !m_hDeltaRInterpolatedXY[side] || !m_hDeltaPhiInterpolatedXY[side])
  {
    return false;
  }

  ++m_fullMapAttemptedBins;

  const double xReadout = m_hDeltaRInterpolatedXY[side]->GetXaxis()->GetBinCenter(xBin);
  const double yReadout = m_hDeltaRInterpolatedXY[side]->GetYaxis()->GetBinCenter(yBin);
  const double rReadout = std::hypot(xReadout, yReadout);
  if (rReadout < m_rMin || rReadout > m_rMax)
  {
    return false;
  }

  const double phiReadout = NormalizePhi(std::atan2(yReadout, xReadout));
  const double zLaunch = side == 1 ? m_zAbsLaunchCm : -m_zAbsLaunchCm;

  std::unique_ptr<TPolyLine3D> path(m_phgarfield->ReverseDrift(xReadout, yReadout, zLaunch));
  if (!path || path->GetN() < 2)
  {
    return false;
  }

  double rCm = 0.0;
  double phiCm = 0.0;
  if (!InterpolateOrClosestAtZ(path.get(), 0.0, rCm, phiCm))
  {
    return false;
  }

  double deltaPhi = phiReadout - phiCm;
  if (deltaPhi > M_PI)
  {
    deltaPhi -= 2.0 * M_PI;
  }
  else if (deltaPhi < -M_PI)
  {
    deltaPhi += 2.0 * M_PI;
  }

  m_hDeltaRInterpolatedXY[side]->SetBinContent(xBin, yBin, rReadout - rCm);
  m_hDeltaPhiInterpolatedXY[side]->SetBinContent(xBin, yBin, deltaPhi);
  ++m_fullMapSuccessfulBins;

  return true;
}

bool TpcLaserClusterGarfieldPropagator::InterpolateOrClosestAtZ(TPolyLine3D* path, double zTarget, double& rOut, double& phiOut) const
{
  if (!path || path->GetN() < 1)
  {
    return false;
  }

  const int nPoints = path->GetN();
  const float* points = path->GetP();
  double bestDistance = std::numeric_limits<double>::max();
  int bestIndex = -1;

  for (int index = 0; index < nPoints - 1; ++index)
  {
    const double x1 = points[3 * index + 0];
    const double y1 = points[3 * index + 1];
    const double z1 = points[3 * index + 2];
    const double x2 = points[3 * (index + 1) + 0];
    const double y2 = points[3 * (index + 1) + 1];
    const double z2 = points[3 * (index + 1) + 2];

    const double zMin = std::min(z1, z2);
    const double zMax = std::max(z1, z2);
    if (zTarget >= zMin && zTarget <= zMax)
    {
      const double fraction = z2 == z1 ? 0.0 : (zTarget - z1) / (z2 - z1);
      const double x = x1 + fraction * (x2 - x1);
      const double y = y1 + fraction * (y2 - y1);
      rOut = std::hypot(x, y);
      phiOut = NormalizePhi(std::atan2(y, x));
      return true;
    }

    const double distance = std::abs(z1 - zTarget);
    if (distance < bestDistance)
    {
      bestDistance = distance;
      bestIndex = index;
    }
  }

  const double lastDistance = std::abs(points[3 * (nPoints - 1) + 2] - zTarget);
  if (lastDistance < bestDistance)
  {
    bestIndex = nPoints - 1;
  }

  if (bestIndex < 0)
  {
    return false;
  }

  const double x = points[3 * bestIndex + 0];
  const double y = points[3 * bestIndex + 1];
  rOut = std::hypot(x, y);
  phiOut = NormalizePhi(std::atan2(y, x));
  return true;
}

double TpcLaserClusterGarfieldPropagator::NormalizePhi(double phi) const
{
  while (phi < 0.0)
  {
    phi += 2.0 * M_PI;
  }
  while (phi >= 2.0 * M_PI)
  {
    phi -= 2.0 * M_PI;
  }
  return phi;
}

void TpcLaserClusterGarfieldPropagator::WriteOutput()
{
  if (m_hDriftStats)
  {
    m_hDriftStats->SetBinContent(1, m_processedEvents);
    m_hDriftStats->SetBinContent(2, m_seenClusters);
    m_hDriftStats->SetBinContent(3, m_acceptedClusters[0]);
    m_hDriftStats->SetBinContent(4, m_acceptedClusters[1]);
    m_hDriftStats->SetBinContent(5, m_attemptedBins);
    m_hDriftStats->SetBinContent(6, m_successfulBins);
    m_hDriftStats->SetBinContent(7, m_nullPaths);
    m_hDriftStats->SetBinContent(8, m_shortPaths);
    m_hDriftStats->SetBinContent(9, m_noCentralMembranePoint);
    m_hDriftStats->SetBinContent(10, m_outOfRangePropagated);
    m_hDriftStats->SetBinContent(11, m_fullMapAttemptedBins);
    m_hDriftStats->SetBinContent(12, m_fullMapSuccessfulBins);
  }

  TFile outputFile(m_outputfile.c_str(), "RECREATE");
  if (outputFile.IsZombie())
  {
    std::cout << "TpcLaserClusterGarfieldPropagator::WriteOutput - could not open output file: " << m_outputfile << std::endl;
    return;
  }

  outputFile.cd();
  for (int side = 0; side < 2; ++side)
  {
    if (m_hRaw[side])
    {
      m_hRaw[side]->Write();
    }
    if (m_hGarfield[side])
    {
      m_hGarfield[side]->Write();
    }
    if (m_hRawXY[side])
    {
      m_hRawXY[side]->Write();
    }
    if (m_hGarfieldXY[side])
    {
      m_hGarfieldXY[side]->Write();
    }
    if (m_hDeltaR[side])
    {
      m_hDeltaR[side]->Write();
    }
    if (m_hDeltaPhi[side])
    {
      m_hDeltaPhi[side]->Write();
    }
    if (m_hDeltaRInterpolated[side])
    {
      m_hDeltaRInterpolated[side]->Write();
    }
    if (m_hDeltaPhiInterpolated[side])
    {
      m_hDeltaPhiInterpolated[side]->Write();
    }
    if (m_hDeltaRInterpolatedXY[side])
    {
      m_hDeltaRInterpolatedXY[side]->Write();
    }
    if (m_hDeltaPhiInterpolatedXY[side])
    {
      m_hDeltaPhiInterpolatedXY[side]->Write();
    }
  }
  if (m_hDriftStats)
  {
    m_hDriftStats->Write();
  }
  WriteFieldOffBeamOffStripeHistograms(outputFile, m_fieldOffBeamOffStripeFile,
                                       m_hGarfield[0], m_hGarfield[1],
                                       m_hDeltaRInterpolated[0], m_hDeltaPhiInterpolated[0],
                                       m_hDeltaRInterpolatedXY[0], m_hDeltaPhiInterpolatedXY[0],
                                       m_hDeltaRInterpolated[1], m_hDeltaPhiInterpolated[1],
                                       m_hDeltaRInterpolatedXY[1], m_hDeltaPhiInterpolatedXY[1]);
  outputFile.Close();
}
