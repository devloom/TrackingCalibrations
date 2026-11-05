#ifndef TPC_LASER_CLUSTER_GARFIELD_PROPAGATOR_H
#define TPC_LASER_CLUSTER_GARFIELD_PROPAGATOR_H

#include <fun4all/SubsysReco.h>

#include <string>

class LaserClusterContainer;
class PHCompositeNode;
class PHGarfield;
class TH1D;
class TH2D;
class TPolyLine3D;

class TpcLaserClusterGarfieldPropagator : public SubsysReco
{
 public:
  explicit TpcLaserClusterGarfieldPropagator(const std::string& name = "TpcLaserClusterGarfieldPropagator");
  ~TpcLaserClusterGarfieldPropagator() override;

  void setOutputfile(const std::string& outputfile);
  void setPHGarfield(PHGarfield* phgarfield);
  void setZLaunch(double zAbsCm = 102.0);
  void setHistogramBinning(int nPhi, int nR, double rMin, double rMax);
  void setDeltaMapBinning(int nPhi, int nR);
  void setFieldOffBeamOffStripeFile(const std::string& filename);
  void setStripeComparisonBranchWeights(double residualWeight, double unmatchedWeight, double cleaningWeight);

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

 private:
  int GetNodes(PHCompositeNode* topNode);
  void CreateHistograms();
  void PropagateOccupiedBins();
  void FillFullDeltaMaps();
  bool PropagateBin(int side, int phiBin, int rBin, double weight);
  bool FillFullDeltaMapBin(int side, int phiBin, int rBin);
  bool FillFullDeltaMapXYBin(int side, int xBin, int yBin);
  bool InterpolateOrClosestAtZ(TPolyLine3D* path, double zTarget, double& rOut, double& phiOut) const;
  double NormalizePhi(double phi) const;
  void WriteOutput();

  std::string m_outputfile{"LaserClusterGarfieldPropagation.root"};
  std::string m_fieldOffBeamOffStripeFile;
  PHGarfield* m_phgarfield{nullptr};
  LaserClusterContainer* m_laserClusterContainer{nullptr};

  int m_phiBins{8000};
  int m_rBins{500};
  int m_deltaMapPhiBins{48};
  int m_deltaMapRBins{48};
  double m_rMin{30.0};
  double m_rMax{80.0};
  double m_zAbsLaunchCm{102.0};
  double m_branchResidualWeight{1.0};
  double m_branchUnmatchedWeight{1.0};
  double m_branchCleaningWeight{1.0};

  TH2D* m_hRaw[2]{nullptr, nullptr};
  TH2D* m_hGarfield[2]{nullptr, nullptr};
  TH2D* m_hRawXY[2]{nullptr, nullptr};
  TH2D* m_hGarfieldXY[2]{nullptr, nullptr};
  TH2D* m_hDeltaR[2]{nullptr, nullptr};
  TH2D* m_hDeltaPhi[2]{nullptr, nullptr};
  TH2D* m_hDeltaRInterpolated[2]{nullptr, nullptr};
  TH2D* m_hDeltaPhiInterpolated[2]{nullptr, nullptr};
  TH2D* m_hDeltaRInterpolatedXY[2]{nullptr, nullptr};
  TH2D* m_hDeltaPhiInterpolatedXY[2]{nullptr, nullptr};
  TH1D* m_hDriftStats{nullptr};

  long long m_processedEvents{0};
  long long m_seenClusters{0};
  long long m_acceptedClusters[2]{0, 0};
  long long m_attemptedBins{0};
  long long m_successfulBins{0};
  long long m_nullPaths{0};
  long long m_shortPaths{0};
  long long m_noCentralMembranePoint{0};
  long long m_outOfRangePropagated{0};
  long long m_fullMapAttemptedBins{0};
  long long m_fullMapSuccessfulBins{0};
};

#endif
