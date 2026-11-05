#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)

std::pair<double,double> SetCommonYRange(const std::vector<TH1*>& histograms);

void plot1D_Zbin(TH3* h3, TH1* h1, float z, bool convert_P_2_RP)
{
  int xbin = h3->GetXaxis()->FindBin((4.55503+4.85652)/2.);
  int zbin = h3->GetZaxis()->FindBin(z);
  for (int i = 1; i <= h3->GetNbinsY(); i++)
  {
      double value = h3->GetBinContent(xbin, i, zbin);
      double error = h3->GetBinError(xbin, i, zbin);
      if (convert_P_2_RP)
      {
        double rcenter = h3->GetYaxis()->GetBinCenter(i);
	//cout<<"rcenter = "<<rcenter<<" , i = "<<i<<endl;
        value *= rcenter;
        error *= rcenter;
      }
      h1->SetBinContent(i, value);
      h1->SetBinError(i, error);
  }
}

void plot1D(TH2* h2, TH1* h1, bool convert_P_2_RP)
{
  int xbin = h2->GetXaxis()->FindBin(4.55073);
  for (int i = 1; i <= h2->GetNbinsY(); i++)
  {
      int xminbin = h2->GetXaxis()->FindBin(4.55073);
      int xmaxbin = h2->GetXaxis()->FindBin(4.84711);
      float bincontent=0, binerror=0;
      int nbin=0;
      for (int j = xminbin; j <= xmaxbin; j++)
      {
        bincontent+=h2->GetBinContent(j, i);
        binerror+=h2->GetBinError(j, i);
	      nbin++;
      }
      bincontent/=nbin;
      binerror/=nbin;
      if (convert_P_2_RP)
        {
          bincontent *= h2->GetYaxis()->GetBinCenter(i);
          binerror *= h2->GetYaxis()->GetBinCenter(i);
        }
      h1->SetBinContent(i, bincontent);
      h1->SetBinError(i, binerror);
  }
}

void draw1D_r_from3D()
{
gStyle->SetOptStat(0);
TGaxis::SetMaxDigits(3);

const int nrun = 1;
//int runs[nrun] = {79507,79509,79510,79511,79512,79513,79514,79515,79516};
// int runs[nrun] = {79507,79516};
int runs[nrun] = {79516};
// int runs[nrun] = {79507};

//get map from cdb lamination study
std::string cdbfilename[nrun];
TFile* cdbfile[nrun];
TH2 *hcdb_N_pr_pos[nrun], *hcdb_R_pr_pos[nrun], *hcdb_P_pr_pos[nrun], *hcdb_Z_pr_pos[nrun];
TH2 *hcdb_N_pr_neg[nrun], *hcdb_R_pr_neg[nrun], *hcdb_P_pr_neg[nrun], *hcdb_Z_pr_neg[nrun];
TH1 *hcdb_N_pos[nrun], *hcdb_R_pos[nrun], *hcdb_P_pos[nrun], *hcdb_Z_pos[nrun];
TH1 *hcdb_N_neg[nrun], *hcdb_R_neg[nrun], *hcdb_P_neg[nrun], *hcdb_Z_neg[nrun];

for (int i = 0; i < nrun; i++)
{
  // auto rc = recoConsts::instance();
  // rc->set_uint64Flag("TIMESTAMP", runs[i]);
  // rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");
  // cdbfilename[i] = CDBInterface::instance()->getUrl("TPC_LAMINATION_FIT_CORRECTION");
  // cout<<"CDB file path: "<<cdbfilename[i].c_str()<<endl;
  // cdbfile[i] = new TFile(cdbfilename[i].c_str(),"");
  

  // ************ If changing from CDB to local files, make sure to set third argument of plot1D to false to avoid multiplying by R again ******************* //
  // cdbfile[i] = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationfits/79516_fieldonbeamon/globalalignedcoordinates_mariiaalignment_81000reference/Laminations_run3pp_new_nocdbtag_v001-00079516.root", "READ");
  // cdbfile[i] = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationfits/79507_fieldonbeamon/globalalignedcoordinates_mariiaalignment_81000reference/Laminations_run3pp_new_nocdbtag_v001-00079507.root", "READ");
  // cdbfile[i] = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationfits/79516_fieldonbeamon/idealtpccoordinates/Laminations_run3pp_new_nocdbtag_v001-00079516.root", "READ");
  

  cdbfile[i] = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationfits/79516_fieldonbeamon/idealtpccoordinates_83319fieldonbeamoffwstaticcorrectionsstripereference/Laminations_run3pp_new_nocdbtag_v001-00079516.root", "READ");

  hcdb_R_pr_pos[i] = (TH2*) cdbfile[i]->Get("hIntDistortionR_posz");
  hcdb_P_pr_pos[i] = (TH2*) cdbfile[i]->Get("hIntDistortionP_posz");
  hcdb_Z_pr_pos[i] = (TH2*) cdbfile[i]->Get("hIntDistortionZ_posz");
  hcdb_N_pr_pos[i] = (TH2*) cdbfile[i]->Get("hEntries_negz");
  hcdb_R_pr_neg[i] = (TH2*) cdbfile[i]->Get("hIntDistortionR_negz");
  hcdb_P_pr_neg[i] = (TH2*) cdbfile[i]->Get("hIntDistortionP_negz");
  hcdb_Z_pr_neg[i] = (TH2*) cdbfile[i]->Get("hIntDistortionZ_negz");
  hcdb_N_pr_neg[i] = (TH2*) cdbfile[i]->Get("hEntries_negz");

  hcdb_N_pos[i] = new TH1F(Form("hentries_posz_%d_cdb",runs[i]),Form("N;R (cm);N"),hcdb_N_pr_pos[i]->GetYaxis()->GetNbins(),hcdb_N_pr_pos[i]->GetYaxis()->GetXmin(),hcdb_N_pr_pos[i]->GetYaxis()->GetXmax());
  hcdb_R_pos[i] = new TH1F(Form("hIntDistortionR_posz_%d_cdb",runs[i]),Form("dR;R (cm);dR (cm)"),hcdb_R_pr_pos[i]->GetYaxis()->GetNbins(),hcdb_R_pr_pos[i]->GetYaxis()->GetXmin(),hcdb_R_pr_pos[i]->GetYaxis()->GetXmax());
  hcdb_P_pos[i] = new TH1F(Form("hIntDistortionP_posz_%d_cdb",runs[i]),Form("Rdphi;R (cm);Rdphi (cm)"),hcdb_P_pr_pos[i]->GetYaxis()->GetNbins(),hcdb_P_pr_pos[i]->GetYaxis()->GetXmin(),hcdb_P_pr_pos[i]->GetYaxis()->GetXmax());
  hcdb_Z_pos[i] = new TH1F(Form("hIntDistortionZ_posz_%d_cdb",runs[i]),Form("dz;R (cm);dz (cm)"),hcdb_Z_pr_pos[i]->GetYaxis()->GetNbins(),hcdb_Z_pr_pos[i]->GetYaxis()->GetXmin(),hcdb_Z_pr_pos[i]->GetYaxis()->GetXmax());
  hcdb_N_neg[i] = new TH1F(Form("hentries_negz_%d_cdb",runs[i]),Form("N;R (cm);N"),hcdb_N_pr_neg[i]->GetYaxis()->GetNbins(),hcdb_N_pr_neg[i]->GetYaxis()->GetXmin(),hcdb_N_pr_neg[i]->GetYaxis()->GetXmax());
  hcdb_R_neg[i] = new TH1F(Form("hIntDistortionR_negz_%d_cdb",runs[i]),Form("dR;R (cm);dR (cm)"),hcdb_R_pr_neg[i]->GetYaxis()->GetNbins(),hcdb_R_pr_neg[i]->GetYaxis()->GetXmin(),hcdb_R_pr_neg[i]->GetYaxis()->GetXmax());
  hcdb_P_neg[i] = new TH1F(Form("hIntDistortionP_negz_%d_cdb",runs[i]),Form("Rdphi;R (cm);Rdphi (cm)"),hcdb_P_pr_neg[i]->GetYaxis()->GetNbins(),hcdb_P_pr_neg[i]->GetYaxis()->GetXmin(),hcdb_P_pr_neg[i]->GetYaxis()->GetXmax());
  hcdb_Z_neg[i] = new TH1F(Form("hIntDistortionZ_negz_%d_cdb",runs[i]),Form("dz;R (cm);dz (cm)"),hcdb_Z_pr_neg[i]->GetYaxis()->GetNbins(),hcdb_Z_pr_neg[i]->GetYaxis()->GetXmin(),hcdb_Z_pr_neg[i]->GetYaxis()->GetXmax());
  plot1D(hcdb_N_pr_pos[i],hcdb_N_pos[i], false);
  plot1D(hcdb_R_pr_pos[i],hcdb_R_pos[i], false);
  plot1D(hcdb_P_pr_pos[i],hcdb_P_pos[i], true);
  plot1D(hcdb_Z_pr_pos[i],hcdb_Z_pos[i], false);
  plot1D(hcdb_N_pr_neg[i],hcdb_N_neg[i], false);
  plot1D(hcdb_R_pr_neg[i],hcdb_R_neg[i], false);
  plot1D(hcdb_P_pr_neg[i],hcdb_P_neg[i], true);
  plot1D(hcdb_Z_pr_neg[i],hcdb_Z_neg[i], false);
  hcdb_N_pos[i]->SetLineColor(1); hcdb_N_pos[i]->SetLineWidth(1); hcdb_N_pos[i]->SetFillColor(0); hcdb_N_pos[i]->SetMarkerColor(1);
  hcdb_R_pos[i]->SetLineColor(1); hcdb_R_pos[i]->SetLineWidth(1); hcdb_R_pos[i]->SetFillColor(0); hcdb_R_pos[i]->SetMarkerColor(1);
  hcdb_P_pos[i]->SetLineColor(1); hcdb_P_pos[i]->SetLineWidth(1); hcdb_P_pos[i]->SetFillColor(0); hcdb_P_pos[i]->SetMarkerColor(1);
  hcdb_Z_pos[i]->SetLineColor(1); hcdb_Z_pos[i]->SetLineWidth(1); hcdb_Z_pos[i]->SetFillColor(0); hcdb_Z_pos[i]->SetMarkerColor(1);
  hcdb_N_neg[i]->SetLineColor(1); hcdb_N_neg[i]->SetLineWidth(1); hcdb_N_neg[i]->SetFillColor(0); hcdb_R_neg[i]->SetMarkerColor(1);
  hcdb_R_neg[i]->SetLineColor(1); hcdb_R_neg[i]->SetLineWidth(1); hcdb_R_neg[i]->SetFillColor(0); hcdb_R_neg[i]->SetMarkerColor(1);
  hcdb_P_neg[i]->SetLineColor(1); hcdb_P_neg[i]->SetLineWidth(1); hcdb_P_neg[i]->SetFillColor(0); hcdb_P_neg[i]->SetMarkerColor(1);
  hcdb_Z_neg[i]->SetLineColor(1); hcdb_Z_neg[i]->SetLineWidth(1); hcdb_Z_neg[i]->SetFillColor(0); hcdb_Z_neg[i]->SetMarkerColor(1);
}

std::vector<float> selectZs={5, 10, 15, 30, 60, 80};
//std::vector<float> selectZs={5};
for (const auto& selectZ : selectZs)
{

TFile* file_3D_map_method1[nrun];
TH3 *h_N_prz_pos_method1[nrun], *h_R_prz_pos_method1[nrun], *h_P_prz_pos_method1[nrun], *h_Z_prz_pos_method1[nrun];
TH3 *h_N_prz_neg_method1[nrun], *h_R_prz_neg_method1[nrun], *h_P_prz_neg_method1[nrun], *h_Z_prz_neg_method1[nrun];
TH1 *h_N_pos_method1[nrun], *h_R_pos_method1[nrun], *h_P_pos_method1[nrun], *h_Z_pos_method1[nrun];
TH1 *h_N_neg_method1[nrun], *h_R_neg_method1[nrun], *h_P_neg_method1[nrun], *h_Z_neg_method1[nrun];

TFile* file_3D_map_method2[nrun];
TH3 *h_N_prz_pos_method2[nrun], *h_R_prz_pos_method2[nrun], *h_P_prz_pos_method2[nrun], *h_Z_prz_pos_method2[nrun];
TH3 *h_N_prz_neg_method2[nrun], *h_R_prz_neg_method2[nrun], *h_P_prz_neg_method2[nrun], *h_Z_prz_neg_method2[nrun];
TH1 *h_N_pos_method2[nrun], *h_R_pos_method2[nrun], *h_P_pos_method2[nrun], *h_Z_pos_method2[nrun];
TH1 *h_N_neg_method2[nrun], *h_R_neg_method2[nrun], *h_P_neg_method2[nrun], *h_Z_neg_method2[nrun];

TFile* file_3D_map_method3[nrun];
TH3 *h_N_prz_pos_method3[nrun], *h_R_prz_pos_method3[nrun], *h_P_prz_pos_method3[nrun], *h_Z_prz_pos_method3[nrun];
TH3 *h_N_prz_neg_method3[nrun], *h_R_prz_neg_method3[nrun], *h_P_prz_neg_method3[nrun], *h_Z_prz_neg_method3[nrun];
TH1 *h_N_pos_method3[nrun], *h_R_pos_method3[nrun], *h_P_pos_method3[nrun], *h_Z_pos_method3[nrun];
TH1 *h_N_neg_method3[nrun], *h_R_neg_method3[nrun], *h_P_neg_method3[nrun], *h_Z_neg_method3[nrun];

for (int i = 0; i < nrun; i++)
{
  file_3D_map_method1[i] = new TFile(Form("/sphenix/u/xyu3/workarea/TPCdistortion/Si_TPOT_fit/run3pp/jobB/Rootfiles/Distortions_full_mm_%d.root",runs[i]),"");
  h_R_prz_pos_method1[i] = (TH3*) file_3D_map_method1[i]->Get("hIntDistortionR_posz");
  h_P_prz_pos_method1[i] = (TH3*) file_3D_map_method1[i]->Get("hIntDistortionP_posz");
  h_Z_prz_pos_method1[i] = (TH3*) file_3D_map_method1[i]->Get("hIntDistortionZ_posz");
  h_N_prz_pos_method1[i] = (TH3*) file_3D_map_method1[i]->Get("hentries_posz");
  h_R_prz_neg_method1[i] = (TH3*) file_3D_map_method1[i]->Get("hIntDistortionR_negz");
  h_P_prz_neg_method1[i] = (TH3*) file_3D_map_method1[i]->Get("hIntDistortionP_negz");
  h_Z_prz_neg_method1[i] = (TH3*) file_3D_map_method1[i]->Get("hIntDistortionZ_negz");
  h_N_prz_neg_method1[i] = (TH3*) file_3D_map_method1[i]->Get("hentries_negz");

  h_N_pos_method1[i] = new TH1F(Form("hentries_posz_%d",runs[i]),Form("N @ Z=%d cm;R (cm);N",(int)selectZ),h_N_prz_pos_method1[i]->GetYaxis()->GetNbins(),h_N_prz_pos_method1[i]->GetYaxis()->GetXmin(),h_N_prz_pos_method1[i]->GetYaxis()->GetXmax());
  h_R_pos_method1[i] = new TH1F(Form("hIntDistortionR_posz_%d",runs[i]),Form("dR @ Z=%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_pos_method1[i]->GetYaxis()->GetNbins(),h_R_prz_pos_method1[i]->GetYaxis()->GetXmin(),h_R_prz_pos_method1[i]->GetYaxis()->GetXmax());
  h_P_pos_method1[i] = new TH1F(Form("hIntDistortionP_posz_%d",runs[i]),Form("Rdphi @ Z=%d cm;R (cm);Rdphi (cm)",(int)selectZ),h_P_prz_pos_method1[i]->GetYaxis()->GetNbins(),h_P_prz_pos_method1[i]->GetYaxis()->GetXmin(),h_P_prz_pos_method1[i]->GetYaxis()->GetXmax());
  h_Z_pos_method1[i] = new TH1F(Form("hIntDistortionZ_posz_%d",runs[i]),Form("dz @ Z=%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_pos_method1[i]->GetYaxis()->GetNbins(),h_Z_prz_pos_method1[i]->GetYaxis()->GetXmin(),h_Z_prz_pos_method1[i]->GetYaxis()->GetXmax());
  h_N_neg_method1[i] = new TH1F(Form("hentries_negz_%d",runs[i]),Form("N @ Z=-%d cm;R (cm);N",(int)selectZ),h_N_prz_neg_method1[i]->GetYaxis()->GetNbins(),h_N_prz_neg_method1[i]->GetYaxis()->GetXmin(),h_N_prz_neg_method1[i]->GetYaxis()->GetXmax());
  h_R_neg_method1[i] = new TH1F(Form("hIntDistortionR_negz_%d",runs[i]),Form("dR @ Z=-%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_neg_method1[i]->GetYaxis()->GetNbins(),h_R_prz_neg_method1[i]->GetYaxis()->GetXmin(),h_R_prz_neg_method1[i]->GetYaxis()->GetXmax());
  h_P_neg_method1[i] = new TH1F(Form("hIntDistortionP_negz_%d",runs[i]),Form("Rdphi @ Z=-%d cm;R (cm);Rdphi (cm)",(int)selectZ),h_P_prz_neg_method1[i]->GetYaxis()->GetNbins(),h_P_prz_neg_method1[i]->GetYaxis()->GetXmin(),h_P_prz_neg_method1[i]->GetYaxis()->GetXmax());
  h_Z_neg_method1[i] = new TH1F(Form("hIntDistortionZ_negz_%d",runs[i]),Form("dz @ Z=-%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_neg_method1[i]->GetYaxis()->GetNbins(),h_Z_prz_neg_method1[i]->GetYaxis()->GetXmin(),h_Z_prz_neg_method1[i]->GetYaxis()->GetXmax());
  plot1D_Zbin(h_N_prz_pos_method1[i],h_N_pos_method1[i],selectZ, false);
  plot1D_Zbin(h_R_prz_pos_method1[i],h_R_pos_method1[i],selectZ, false);
  plot1D_Zbin(h_P_prz_pos_method1[i],h_P_pos_method1[i],selectZ, true);
  plot1D_Zbin(h_Z_prz_pos_method1[i],h_Z_pos_method1[i],selectZ, false);
  plot1D_Zbin(h_N_prz_neg_method1[i],h_N_neg_method1[i],-selectZ, false);
  plot1D_Zbin(h_R_prz_neg_method1[i],h_R_neg_method1[i],-selectZ, false);
  plot1D_Zbin(h_P_prz_neg_method1[i],h_P_neg_method1[i],-selectZ, true);
  plot1D_Zbin(h_Z_prz_neg_method1[i],h_Z_neg_method1[i],-selectZ, false);
  h_N_pos_method1[i]->SetLineColor(2); h_N_pos_method1[i]->SetLineWidth(1); h_N_pos_method1[i]->SetFillColor(0); h_N_pos_method1[i]->SetMarkerColor(2);
  h_R_pos_method1[i]->SetLineColor(2); h_R_pos_method1[i]->SetLineWidth(1); h_R_pos_method1[i]->SetFillColor(0); h_R_pos_method1[i]->SetMarkerColor(2);
  h_P_pos_method1[i]->SetLineColor(2); h_P_pos_method1[i]->SetLineWidth(1); h_P_pos_method1[i]->SetFillColor(0); h_P_pos_method1[i]->SetMarkerColor(2);
  h_Z_pos_method1[i]->SetLineColor(2); h_Z_pos_method1[i]->SetLineWidth(1); h_Z_pos_method1[i]->SetFillColor(0); h_Z_pos_method1[i]->SetMarkerColor(2);
  h_N_neg_method1[i]->SetLineColor(2); h_N_neg_method1[i]->SetLineWidth(1); h_N_neg_method1[i]->SetFillColor(0); h_R_neg_method1[i]->SetMarkerColor(2);
  h_R_neg_method1[i]->SetLineColor(2); h_R_neg_method1[i]->SetLineWidth(1); h_R_neg_method1[i]->SetFillColor(0); h_R_neg_method1[i]->SetMarkerColor(2);
  h_P_neg_method1[i]->SetLineColor(2); h_P_neg_method1[i]->SetLineWidth(1); h_P_neg_method1[i]->SetFillColor(0); h_P_neg_method1[i]->SetMarkerColor(2);
  h_Z_neg_method1[i]->SetLineColor(2); h_Z_neg_method1[i]->SetLineWidth(1); h_Z_neg_method1[i]->SetFillColor(0); h_Z_neg_method1[i]->SetMarkerColor(2);

  file_3D_map_method2[i] = new TFile(Form("/sphenix/u/xyu3/workarea/TPCdistortion/Si_TPOT_fit/run3pp_newclusErrorPara/jobB/Rootfiles/Distortions_full_mm_%d.root",runs[i]),"");
  h_R_prz_pos_method2[i] = (TH3*) file_3D_map_method2[i]->Get("hIntDistortionR_posz");
  h_P_prz_pos_method2[i] = (TH3*) file_3D_map_method2[i]->Get("hIntDistortionP_posz");
  h_Z_prz_pos_method2[i] = (TH3*) file_3D_map_method2[i]->Get("hIntDistortionZ_posz");
  h_N_prz_pos_method2[i] = (TH3*) file_3D_map_method2[i]->Get("hentries_posz");
  h_R_prz_neg_method2[i] = (TH3*) file_3D_map_method2[i]->Get("hIntDistortionR_negz");
  h_P_prz_neg_method2[i] = (TH3*) file_3D_map_method2[i]->Get("hIntDistortionP_negz");
  h_Z_prz_neg_method2[i] = (TH3*) file_3D_map_method2[i]->Get("hIntDistortionZ_negz");
  h_N_prz_neg_method2[i] = (TH3*) file_3D_map_method2[i]->Get("hentries_negz");

  h_N_pos_method2[i] = new TH1F(Form("hentries_posz_%d",runs[i]),Form("N @ Z=%d cm;R (cm);N",(int)selectZ),h_N_prz_pos_method2[i]->GetYaxis()->GetNbins(),h_N_prz_pos_method2[i]->GetYaxis()->GetXmin(),h_N_prz_pos_method2[i]->GetYaxis()->GetXmax());
  h_R_pos_method2[i] = new TH1F(Form("hIntDistortionR_posz_%d",runs[i]),Form("dR @ Z=%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_pos_method2[i]->GetYaxis()->GetNbins(),h_R_prz_pos_method2[i]->GetYaxis()->GetXmin(),h_R_prz_pos_method2[i]->GetYaxis()->GetXmax());
  h_P_pos_method2[i] = new TH1F(Form("hIntDistortionP_posz_%d",runs[i]),Form("Rdphi @ Z=%d cm;R (cm);Rdphi (cm)",(int)selectZ),h_P_prz_pos_method2[i]->GetYaxis()->GetNbins(),h_P_prz_pos_method2[i]->GetYaxis()->GetXmin(),h_P_prz_pos_method2[i]->GetYaxis()->GetXmax());
  h_Z_pos_method2[i] = new TH1F(Form("hIntDistortionZ_posz_%d",runs[i]),Form("dz @ Z=%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_pos_method2[i]->GetYaxis()->GetNbins(),h_Z_prz_pos_method2[i]->GetYaxis()->GetXmin(),h_Z_prz_pos_method2[i]->GetYaxis()->GetXmax());
  h_N_neg_method2[i] = new TH1F(Form("hentries_negz_%d",runs[i]),Form("N @ Z=-%d cm;R (cm);N",(int)selectZ),h_N_prz_neg_method2[i]->GetYaxis()->GetNbins(),h_N_prz_neg_method2[i]->GetYaxis()->GetXmin(),h_N_prz_neg_method2[i]->GetYaxis()->GetXmax());
  h_R_neg_method2[i] = new TH1F(Form("hIntDistortionR_negz_%d",runs[i]),Form("dR @ Z=-%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_neg_method2[i]->GetYaxis()->GetNbins(),h_R_prz_neg_method2[i]->GetYaxis()->GetXmin(),h_R_prz_neg_method2[i]->GetYaxis()->GetXmax());
  h_P_neg_method2[i] = new TH1F(Form("hIntDistortionP_negz_%d",runs[i]),Form("Rdphi @ Z=-%d cm;R (cm);Rdphi (cm)",(int)selectZ),h_P_prz_neg_method2[i]->GetYaxis()->GetNbins(),h_P_prz_neg_method2[i]->GetYaxis()->GetXmin(),h_P_prz_neg_method2[i]->GetYaxis()->GetXmax());
  h_Z_neg_method2[i] = new TH1F(Form("hIntDistortionZ_negz_%d",runs[i]),Form("dz @ Z=-%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_neg_method2[i]->GetYaxis()->GetNbins(),h_Z_prz_neg_method2[i]->GetYaxis()->GetXmin(),h_Z_prz_neg_method2[i]->GetYaxis()->GetXmax());
  plot1D_Zbin(h_N_prz_pos_method2[i],h_N_pos_method2[i],selectZ, false);
  plot1D_Zbin(h_R_prz_pos_method2[i],h_R_pos_method2[i],selectZ, false);
  plot1D_Zbin(h_P_prz_pos_method2[i],h_P_pos_method2[i],selectZ, true);
  plot1D_Zbin(h_Z_prz_pos_method2[i],h_Z_pos_method2[i],selectZ, false);
  plot1D_Zbin(h_N_prz_neg_method2[i],h_N_neg_method2[i],-selectZ, false);
  plot1D_Zbin(h_R_prz_neg_method2[i],h_R_neg_method2[i],-selectZ, false);
  plot1D_Zbin(h_P_prz_neg_method2[i],h_P_neg_method2[i],-selectZ, true);
  plot1D_Zbin(h_Z_prz_neg_method2[i],h_Z_neg_method2[i],-selectZ, false);
  h_N_pos_method2[i]->SetLineColor((4)); h_N_pos_method2[i]->SetLineWidth(1); h_N_pos_method2[i]->SetFillColor(0); h_N_pos_method2[i]->SetMarkerColor((4));
  h_R_pos_method2[i]->SetLineColor((4)); h_R_pos_method2[i]->SetLineWidth(1); h_R_pos_method2[i]->SetFillColor(0); h_R_pos_method2[i]->SetMarkerColor((4));
  h_P_pos_method2[i]->SetLineColor((4)); h_P_pos_method2[i]->SetLineWidth(1); h_P_pos_method2[i]->SetFillColor(0); h_P_pos_method2[i]->SetMarkerColor((4));
  h_Z_pos_method2[i]->SetLineColor((4)); h_Z_pos_method2[i]->SetLineWidth(1); h_Z_pos_method2[i]->SetFillColor(0); h_Z_pos_method2[i]->SetMarkerColor((4));
  h_N_neg_method2[i]->SetLineColor((4)); h_N_neg_method2[i]->SetLineWidth(1); h_N_neg_method2[i]->SetFillColor(0); h_R_neg_method2[i]->SetMarkerColor((4));
  h_R_neg_method2[i]->SetLineColor((4)); h_R_neg_method2[i]->SetLineWidth(1); h_R_neg_method2[i]->SetFillColor(0); h_R_neg_method2[i]->SetMarkerColor((4));
  h_P_neg_method2[i]->SetLineColor((4)); h_P_neg_method2[i]->SetLineWidth(1); h_P_neg_method2[i]->SetFillColor(0); h_P_neg_method2[i]->SetMarkerColor((4));
  h_Z_neg_method2[i]->SetLineColor((4)); h_Z_neg_method2[i]->SetLineWidth(1); h_Z_neg_method2[i]->SetFillColor(0); h_Z_neg_method2[i]->SetMarkerColor((4));

  file_3D_map_method3[i] = new TFile(Form("/sphenix/u/xyu3/workarea/TPCdistortion/Si_TPOT_fit/run3pp_newclusErrorParaTPOT/jobB/Rootfiles/Distortions_full_mm_%d.root",runs[i]),"");
  h_R_prz_pos_method3[i] = (TH3*) file_3D_map_method3[i]->Get("hIntDistortionR_posz");
  h_P_prz_pos_method3[i] = (TH3*) file_3D_map_method3[i]->Get("hIntDistortionP_posz");
  h_Z_prz_pos_method3[i] = (TH3*) file_3D_map_method3[i]->Get("hIntDistortionZ_posz");
  h_N_prz_pos_method3[i] = (TH3*) file_3D_map_method3[i]->Get("hentries_posz");
  h_R_prz_neg_method3[i] = (TH3*) file_3D_map_method3[i]->Get("hIntDistortionR_negz");
  h_P_prz_neg_method3[i] = (TH3*) file_3D_map_method3[i]->Get("hIntDistortionP_negz");
  h_Z_prz_neg_method3[i] = (TH3*) file_3D_map_method3[i]->Get("hIntDistortionZ_negz");
  h_N_prz_neg_method3[i] = (TH3*) file_3D_map_method3[i]->Get("hentries_negz");

  h_N_pos_method3[i] = new TH1F(Form("hentries_posz_%d",runs[i]),Form("N @ Z=%d cm;R (cm);N",(int)selectZ),h_N_prz_pos_method3[i]->GetYaxis()->GetNbins(),h_N_prz_pos_method3[i]->GetYaxis()->GetXmin(),h_N_prz_pos_method3[i]->GetYaxis()->GetXmax());
  h_R_pos_method3[i] = new TH1F(Form("hIntDistortionR_posz_%d",runs[i]),Form("dR @ Z=%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_pos_method3[i]->GetYaxis()->GetNbins(),h_R_prz_pos_method3[i]->GetYaxis()->GetXmin(),h_R_prz_pos_method3[i]->GetYaxis()->GetXmax());
  h_P_pos_method3[i] = new TH1F(Form("hIntDistortionP_posz_%d",runs[i]),Form("Rdphi @ Z=%d cm;R (cm);Rdphi (cm)",(int)selectZ),h_P_prz_pos_method3[i]->GetYaxis()->GetNbins(),h_P_prz_pos_method3[i]->GetYaxis()->GetXmin(),h_P_prz_pos_method3[i]->GetYaxis()->GetXmax());
  h_Z_pos_method3[i] = new TH1F(Form("hIntDistortionZ_posz_%d",runs[i]),Form("dz @ Z=%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_pos_method3[i]->GetYaxis()->GetNbins(),h_Z_prz_pos_method3[i]->GetYaxis()->GetXmin(),h_Z_prz_pos_method3[i]->GetYaxis()->GetXmax());
  h_N_neg_method3[i] = new TH1F(Form("hentries_negz_%d",runs[i]),Form("N @ Z=-%d cm;R (cm);N",(int)selectZ),h_N_prz_neg_method3[i]->GetYaxis()->GetNbins(),h_N_prz_neg_method3[i]->GetYaxis()->GetXmin(),h_N_prz_neg_method3[i]->GetYaxis()->GetXmax());
  h_R_neg_method3[i] = new TH1F(Form("hIntDistortionR_negz_%d",runs[i]),Form("dR @ Z=-%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_neg_method3[i]->GetYaxis()->GetNbins(),h_R_prz_neg_method3[i]->GetYaxis()->GetXmin(),h_R_prz_neg_method3[i]->GetYaxis()->GetXmax());
  h_P_neg_method3[i] = new TH1F(Form("hIntDistortionP_negz_%d",runs[i]),Form("Rdphi @ Z=-%d cm;R (cm);Rdphi (cm)",(int)selectZ),h_P_prz_neg_method3[i]->GetYaxis()->GetNbins(),h_P_prz_neg_method3[i]->GetYaxis()->GetXmin(),h_P_prz_neg_method3[i]->GetYaxis()->GetXmax());
  h_Z_neg_method3[i] = new TH1F(Form("hIntDistortionZ_negz_%d",runs[i]),Form("dz @ Z=-%d cm;R (cm);dz (cm)",(int)selectZ),h_Z_prz_neg_method3[i]->GetYaxis()->GetNbins(),h_Z_prz_neg_method3[i]->GetYaxis()->GetXmin(),h_Z_prz_neg_method3[i]->GetYaxis()->GetXmax());
  plot1D_Zbin(h_N_prz_pos_method3[i],h_N_pos_method3[i],selectZ, false);
  plot1D_Zbin(h_R_prz_pos_method3[i],h_R_pos_method3[i],selectZ, false);
  plot1D_Zbin(h_P_prz_pos_method3[i],h_P_pos_method3[i],selectZ, true);
  plot1D_Zbin(h_Z_prz_pos_method3[i],h_Z_pos_method3[i],selectZ, false);
  plot1D_Zbin(h_N_prz_neg_method3[i],h_N_neg_method3[i],-selectZ, false);
  plot1D_Zbin(h_R_prz_neg_method3[i],h_R_neg_method3[i],-selectZ, false);
  plot1D_Zbin(h_P_prz_neg_method3[i],h_P_neg_method3[i],-selectZ, true);
  plot1D_Zbin(h_Z_prz_neg_method3[i],h_Z_neg_method3[i],-selectZ, false);
  h_N_pos_method3[i]->SetLineColor((3)); h_N_pos_method3[i]->SetLineWidth(1); h_N_pos_method3[i]->SetFillColor(0); h_N_pos_method3[i]->SetMarkerColor((3));
  h_R_pos_method3[i]->SetLineColor((3)); h_R_pos_method3[i]->SetLineWidth(1); h_R_pos_method3[i]->SetFillColor(0); h_R_pos_method3[i]->SetMarkerColor((3));
  h_P_pos_method3[i]->SetLineColor((3)); h_P_pos_method3[i]->SetLineWidth(1); h_P_pos_method3[i]->SetFillColor(0); h_P_pos_method3[i]->SetMarkerColor((3));
  h_Z_pos_method3[i]->SetLineColor((3)); h_Z_pos_method3[i]->SetLineWidth(1); h_Z_pos_method3[i]->SetFillColor(0); h_Z_pos_method3[i]->SetMarkerColor((3));
  h_N_neg_method3[i]->SetLineColor((3)); h_N_neg_method3[i]->SetLineWidth(1); h_N_neg_method3[i]->SetFillColor(0); h_R_neg_method3[i]->SetMarkerColor((3));
  h_R_neg_method3[i]->SetLineColor((3)); h_R_neg_method3[i]->SetLineWidth(1); h_R_neg_method3[i]->SetFillColor(0); h_R_neg_method3[i]->SetMarkerColor((3));
  h_P_neg_method3[i]->SetLineColor((3)); h_P_neg_method3[i]->SetLineWidth(1); h_P_neg_method3[i]->SetFillColor(0); h_P_neg_method3[i]->SetMarkerColor((3));
  h_Z_neg_method3[i]->SetLineColor((3)); h_Z_neg_method3[i]->SetLineWidth(1); h_Z_neg_method3[i]->SetFillColor(0); h_Z_neg_method3[i]->SetMarkerColor((3));

  std::pair<double,double> yrange_R = SetCommonYRange({h_R_neg_method1[i], h_R_pos_method1[i], h_R_neg_method2[i], h_R_pos_method2[i], h_R_neg_method3[i], h_R_pos_method3[i], hcdb_R_neg[i], hcdb_R_pos[i]});
  std::pair<double,double> yrange_P = SetCommonYRange({h_P_neg_method1[i], h_P_pos_method1[i], h_P_neg_method2[i], h_P_pos_method2[i], h_P_neg_method3[i], h_P_pos_method3[i], hcdb_P_neg[i], hcdb_P_pos[i]});
  std::pair<double,double> yrange_Z = SetCommonYRange({h_Z_neg_method1[i], h_Z_pos_method1[i], h_Z_neg_method2[i], h_Z_pos_method2[i], h_Z_neg_method3[i], h_Z_pos_method3[i], hcdb_Z_neg[i], hcdb_Z_pos[i]});
}

TFile* ofile = new TFile(Form("79516/Rootfiles/hist_dr_1Dr_Z%d.root",(int)selectZ),"recreate");
ofile->cd();
for (int i=0; i<nrun; i++)
{
  h_R_pos_method1[i]->Write();
  h_R_pos_method2[i]->Write();
  h_R_pos_method3[i]->Write();
  hcdb_R_pos[i]->Write();
  h_R_neg_method1[i]->Write();
  h_R_neg_method2[i]->Write();
  h_R_neg_method3[i]->Write();
  hcdb_R_neg[i]->Write();
}
ofile->Write();

TFile* ofile2 = new TFile(Form("79516/Rootfiles/hist_rdphi_1Dr_Z%d.root",(int)selectZ),"recreate");
ofile2->cd();
for (int i=0; i<nrun; i++)
{
  h_P_pos_method1[i]->Write();
  h_P_pos_method2[i]->Write();
  h_P_pos_method3[i]->Write();
  hcdb_P_pos[i]->Write();
  h_P_neg_method1[i]->Write();
  h_P_neg_method2[i]->Write();
  h_P_neg_method3[i]->Write();
  hcdb_P_neg[i]->Write();
}
ofile2->Write();

TFile* ofile3 = new TFile(Form("79516/Rootfiles/hist_dz_1Dr_Z%d.root",(int)selectZ),"recreate");
ofile3->cd();
for (int i=0; i<nrun; i++)
{
  h_Z_pos_method1[i]->Write();
  h_Z_pos_method2[i]->Write();
  h_Z_pos_method3[i]->Write();
  hcdb_Z_pos[i]->Write();
  h_Z_neg_method1[i]->Write();
  h_Z_neg_method2[i]->Write();
  h_Z_neg_method3[i]->Write();
  hcdb_Z_neg[i]->Write();
}
ofile3->Write();



TCanvas *c = new TCanvas("c","",800,600);
c->Divide(1,2);
c->cd(1);
gPad->SetLogy(0);
// h_P_pos_method1[i]->Draw("hist,e,same");
// h_P_pos_method2[i]->Draw("hist,e,same");
h_P_pos_method3[0]->Draw("hist,e,same");
h_P_pos_method3[0]->SetTitle("North");
hcdb_P_pos[0]->Draw("hist,same");
TLegend *legend = new TLegend(0.4, 0.6, 0.8, 0.8);
legend->SetBorderSize(0);
legend->AddEntry(h_P_pos_method3[0], Form("Track based @ Z =#pm 5 cm"), "l");
legend->AddEntry(hcdb_P_pos[0], Form("Lamination"), "l");
legend->Draw();
c->cd(2);
gPad->SetLogy(0);
h_P_neg_method3[0]->SetTitle("South");
h_P_neg_method3[0]->Draw("hist,e,same");
hcdb_P_neg[0]->Draw("hist,same");

c->SaveAs("79516/figure/rdphi_distortions.png");






for (int i=0; i<nrun; i++)
{
TCanvas* can = new TCanvas("can","",3200,1200);
can->Divide(4,2);
can->cd(1);
gPad->SetLogy(0);
// h_P_pos_method1[i]->Draw("hist,e,same");
// h_P_pos_method2[i]->Draw("hist,e,same");
h_P_pos_method3[i]->Draw("hist,e,same");
hcdb_P_pos[i]->Draw("hist,same");
can->cd(2);
gPad->SetLogy(0);
// h_R_pos_method1[i]->Draw("hist,e,same");
// h_R_pos_method2[i]->Draw("hist,e,same");
h_R_pos_method3[i]->Draw("hist,e,same");
hcdb_R_pos[i]->Draw("hist,same");
can->cd(3);
gPad->SetLogy(0);
// h_Z_pos_method1[i]->Draw("hist,e,same");
// h_Z_pos_method2[i]->Draw("hist,e,same");
h_Z_pos_method3[i]->Draw("hist,e,same");
hcdb_Z_pos[i]->Draw("hist,same");
can->cd(4);
gPad->SetLogy(1);
// h_N_pos_method1[i]->Draw("hist,e,same");
// h_N_pos_method2[i]->Draw("hist,e,same");
h_N_pos_method3[i]->Draw("hist,e,same");
//hcdb_N_pos[i]->Draw("hist,same");
can->cd(5);
gPad->SetLogy(0);
// h_P_neg_method1[i]->Draw("hist,e,same");
// h_P_neg_method2[i]->Draw("hist,e,same");
h_P_neg_method3[i]->Draw("hist,e,same");
hcdb_P_neg[i]->Draw("hist,same");
can->cd(6);
gPad->SetLogy(0);
// h_R_neg_method1[i]->Draw("hist,e,same");
// h_R_neg_method2[i]->Draw("hist,e,same");
h_R_neg_method3[i]->Draw("hist,e,same");
hcdb_R_neg[i]->Draw("hist,same");
can->cd(7);
gPad->SetLogy(0);
// h_Z_neg_method1[i]->Draw("hist,e,same");
// h_Z_neg_method2[i]->Draw("hist,e,same");
h_Z_neg_method3[i]->Draw("hist,e,same");
hcdb_Z_neg[i]->Draw("hist,same");
can->cd(8);
gPad->SetLogy(1);
// h_N_neg_method1[i]->Draw("hist,e,same");
// h_N_neg_method2[i]->Draw("hist,e,same");
h_N_neg_method3[i]->Draw("hist,e,same");
//hcdb_N_neg[i]->Draw("hist,same");

gPad->RedrawAxis();

can->Update();
can->SaveAs(Form("79516/figure/resid_vsR_from3D_atZ%d_%d.pdf",(int)selectZ,runs[i]));

TLegend *legend = new TLegend(0.1, 0.1, 0.9, 0.9);
legend->SetHeader(Form("Run %d",runs[i]));
legend->AddEntry(h_P_pos_method1[i], Form("repo cluster error para"), "l");
legend->AddEntry(h_P_pos_method2[i], Form("new cluster error para (w/o TPOT)"), "l");
legend->AddEntry(h_P_pos_method3[i], Form("new cluster error para (w/ TPOT)"), "l");
legend->AddEntry(hcdb_P_pos[i], Form("Lamination fitting in CDB"), "l");
legend->Draw();
TCanvas* can_leg = new TCanvas("can_leg","",1600,1200);
legend->Draw();
can_leg->SaveAs(Form("79516/figure/resid_vsR_from3D_leg_%d.pdf",runs[i]));

delete can;
delete can_leg;
}

}

}

std::pair<double,double> SetCommonYRange(const std::vector<TH1*>& histograms)
{
    Double_t yMin = TMath::Infinity();
    Double_t yMax = -TMath::Infinity();

    for (TH1* h : histograms) {
        if (!h) continue;
        yMin = TMath::Min(yMin, h->GetMinimum());
        yMax = TMath::Max(yMax, h->GetMaximum());
    }

    if (yMin == TMath::Infinity()) yMin = 0;
    if (yMax == -TMath::Infinity()) yMax = 1;

    for (TH1* h : histograms) {
        if (h) {
            h->SetMinimum(yMin);
            h->SetMaximum(yMax);
        }
    }
    return {yMin, yMax};
}
