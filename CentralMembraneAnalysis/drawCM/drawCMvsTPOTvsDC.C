#include <ffamodules/CDBInterface.h>
#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)

std::pair<double,double> SetCommonYRange(const std::vector<TH1*>& histograms);

void plot1D_Zbin(TH3* h3, TH1* h1, float z, bool convert_P_2_RP, bool convert_RP_2_P)
{
  // int xbin = h3->GetXaxis()->FindBin((4.55503+4.85652)/2.);
  int zbin = h3->GetZaxis()->FindBin(z);
  for (int i = 1; i <= h3->GetNbinsY(); i++)
  {
      int xminbin = h3->GetXaxis()->FindBin(0.0);
      int xmaxbin = h3->GetXaxis()->FindBin(6.28);
      // int xminbin = h3->GetXaxis()->FindBin(2.0);
      // int xmaxbin = h3->GetXaxis()->FindBin(3.0);
      float bincontent=0, binerror=0;
      int nbin=0;
      for (int j = xminbin; j <= xmaxbin; j++)
      {
        bincontent+=h3->GetBinContent(j, i, zbin);
        binerror+=h3->GetBinError(j, i, zbin);
	      nbin++;
      }
      bincontent/=nbin;
      binerror/=nbin;
      // double value = h3->GetBinContent(xbin, i, zbin);
      // double error = h3->GetBinError(xbin, i, zbin);
      if (convert_P_2_RP)
      {
        double rcenter = h3->GetYaxis()->GetBinCenter(i);
	//cout<<"rcenter = "<<rcenter<<" , i = "<<i<<endl;
        bincontent *= rcenter;
        binerror *= rcenter;
      }
      if (convert_RP_2_P)
      {
        double rcenter = h3->GetYaxis()->GetBinCenter(i);
	//cout<<"rcenter = "<<rcenter<<" , i = "<<i<<endl;
        bincontent /= rcenter;
        binerror /= rcenter;
      }
      h1->SetBinContent(i, bincontent);
      h1->SetBinError(i, binerror);
  }
}

void plot1D(TH2* h2, TH1* h1, bool convert_P_2_RP, bool convert_RP_2_P)
{
  // int xbin = h2->GetXaxis()->FindBin(4.55073);
  for (int i = 1; i <= h2->GetNbinsY(); i++)
  {
      // int xminbin = h2->GetXaxis()->FindBin(4.55073);
      // int xmaxbin = h2->GetXaxis()->FindBin(4.84711);
      int xminbin = h2->GetXaxis()->FindBin(0.0);
      int xmaxbin = h2->GetXaxis()->FindBin(6.28);
      // int xminbin = h2->GetXaxis()->FindBin(2.0);
      // int xmaxbin = h2->GetXaxis()->FindBin(3.0);
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

      if (convert_RP_2_P)
        {
          bincontent /= h2->GetYaxis()->GetBinCenter(i);
          binerror /= h2->GetYaxis()->GetBinCenter(i);
        }
      h1->SetBinContent(i, bincontent);
      h1->SetBinError(i, binerror);
  }
}

void drawCMvsTPOTvsDC()
{
  gStyle->SetOptStat(0);
  gROOT->SetBatch(kTRUE);

  gROOT->LoadMacro("/sphenix/user/dloomis/Tools/sPhenixStyle.C");
  gROOT->ProcessLine("SetsPhenixStyle();");

  TGaxis::SetMaxDigits(3);

  const int nrun = 1;
  //int runs[nrun] = {79507,79509,79510,79511,79512,79513,79514,79515,79516};
  // int runs[nrun] = {79507,79516};
  int runs[nrun] = {79516};
  // int runs[nrun] = {79507};

  // ********* Central Membrane Distortions From CDB ************ //
  std::string cdbfilename[nrun];
  TFile* cdbfile[nrun];
  TH2 *hcdb_R_pr_pos[nrun], *hcdb_P_pr_pos[nrun];
  TH2 *hcdb_R_pr_neg[nrun], *hcdb_P_pr_neg[nrun];
  TH1 *hcdb_R_pos[nrun], *hcdb_P_pos[nrun];
  TH1 *hcdb_R_neg[nrun], *hcdb_P_neg[nrun];

  for (int i = 0; i < nrun; i++)
  {
    auto rc = recoConsts::instance();
    rc->set_uint64Flag("TIMESTAMP", runs[i]);
    rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");
    cdbfilename[i] = CDBInterface::instance()->getUrl("TPC_LAMINATION_FIT_CORRECTION");
    cout<<"CDB file path: "<<cdbfilename[i].c_str()<<endl;
    cdbfile[i] = new TFile(cdbfilename[i].c_str(),"");
    

    // ************ If changing from CDB to local files, make sure to set third argument of plot1D to false to avoid multiplying by R again ******************* //
    // cdbfile[i] = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationfits/79516_fieldonbeamon/globalalignedcoordinates_mariiaalignment_81000reference/Laminations_run3pp_new_nocdbtag_v001-00079516.root", "READ");
    // cdbfile[i] = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationfits/79507_fieldonbeamon/globalalignedcoordinates_mariiaalignment_81000reference/Laminations_run3pp_new_nocdbtag_v001-00079507.root", "READ");
    // cdbfile[i] = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationfits/79516_fieldonbeamon/idealtpccoordinates/Laminations_run3pp_new_nocdbtag_v001-00079516.root", "READ");
    

    // cdbfile[i] = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationfits/79516_fieldonbeamon/idealtpccoordinates_78242fieldoffbeamoffidealstripereference_updatedreferencestripefits/Laminations_run3pp_new_nocdbtag_v001-00079516.root", "READ");
    // cdbfile[i] = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationfits/79516_fieldonbeamon/idealtpccoordinates_78242fieldoffbeamoffidealstripereference_updatedreferencestripefits_updatedstripefits/Laminations_run3pp_new_nocdbtag_v001-00079516.root", "READ");
    // cdbfile[i] = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationfits/79516_fieldonbeamon/idealtpccoordinates_83319fieldonbeamoffwupdatedstaticcorrectionsidealstripereference_updatedreferencestripefits/Laminations_run3pp_new_nocdbtag_v001-00079516.root", "READ");
    
    // cdbfile[i] = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationfits/79516_fieldonbeamon/idealtpccoordinates_83319fieldonbeamoffwupdatedstaticcorrectionsidealstripereference_updatedreferencestripefits_updatedstripefits/Laminations_run3pp_new_nocdbtag_v001-00079516.root", "READ");
    // cdbfile[i] = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationfits/79507_fieldonbeamon/idealtpccoordinates_83319fieldonbeamoffwupdatedstaticcorrectionsidealstripereference_updatedreferencestripefits_updatedstripefits/Laminations_run3pp_new_nocdbtag_v001-00079507.root", "READ");

    // cdbfile[i] = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationfits/79516_fieldonbeamon/globalalignedcoordinates_mariiaalignment_83313fieldonbeamoffwstaticcorrectionsalignedstripereference/Laminations_run3pp_new_nocdbtag_v001-00079516.root", "READ");
    cdbfile[i] = new TFile("../CentralMembraneStripeMatching/macros/output/CMDistortionCorrections-00079513_reference83319.root", "READ");


    // cdbfile[i] = new TFile("/sphenix/user/dloomis/TrackingAnalysis/CentralMembraneAnalysisv2/macros/CMDistortionCorrections_run3pp_new_nocdbtag_v001-00079507_reference83319fieldonbeamoff_0to20segs.root", "READ");

    hcdb_R_pr_pos[i] = (TH2*) cdbfile[i]->Get("hIntDistortionR_posz");
    hcdb_P_pr_pos[i] = (TH2*) cdbfile[i]->Get("hIntDistortionP_posz");
    hcdb_R_pr_neg[i] = (TH2*) cdbfile[i]->Get("hIntDistortionR_negz");
    hcdb_P_pr_neg[i] = (TH2*) cdbfile[i]->Get("hIntDistortionP_negz");

    hcdb_R_pos[i] = new TH1F(Form("hIntDistortionR_posz_%d_cdb",runs[i]),Form("dR;R (cm);dR (cm)"),hcdb_R_pr_pos[i]->GetYaxis()->GetNbins(),hcdb_R_pr_pos[i]->GetYaxis()->GetXmin(),hcdb_R_pr_pos[i]->GetYaxis()->GetXmax());
    hcdb_P_pos[i] = new TH1F(Form("hIntDistortionP_posz_%d_cdb",runs[i]),Form("Rdphi;R (cm);Rdphi (cm)"),hcdb_P_pr_pos[i]->GetYaxis()->GetNbins(),hcdb_P_pr_pos[i]->GetYaxis()->GetXmin(),hcdb_P_pr_pos[i]->GetYaxis()->GetXmax());
    hcdb_R_neg[i] = new TH1F(Form("hIntDistortionR_negz_%d_cdb",runs[i]),Form("dR;R (cm);dR (cm)"),hcdb_R_pr_neg[i]->GetYaxis()->GetNbins(),hcdb_R_pr_neg[i]->GetYaxis()->GetXmin(),hcdb_R_pr_neg[i]->GetYaxis()->GetXmax());
    hcdb_P_neg[i] = new TH1F(Form("hIntDistortionP_negz_%d_cdb",runs[i]),Form("Rdphi;R (cm);Rdphi (cm)"),hcdb_P_pr_neg[i]->GetYaxis()->GetNbins(),hcdb_P_pr_neg[i]->GetYaxis()->GetXmin(),hcdb_P_pr_neg[i]->GetYaxis()->GetXmax());
    plot1D(hcdb_R_pr_pos[i],hcdb_R_pos[i], false, false);
    plot1D(hcdb_P_pr_pos[i],hcdb_P_pos[i], false, false);
    plot1D(hcdb_R_pr_neg[i],hcdb_R_neg[i], false, false);
    plot1D(hcdb_P_pr_neg[i],hcdb_P_neg[i], false, false);
    hcdb_R_pos[i]->SetLineColor(1); hcdb_R_pos[i]->SetLineWidth(1); hcdb_R_pos[i]->SetFillColor(0); hcdb_R_pos[i]->SetMarkerColor(1);
    hcdb_P_pos[i]->SetLineColor(1); hcdb_P_pos[i]->SetLineWidth(1); hcdb_P_pos[i]->SetFillColor(0); hcdb_P_pos[i]->SetMarkerColor(1);
    hcdb_R_neg[i]->SetLineColor(1); hcdb_R_neg[i]->SetLineWidth(1); hcdb_R_neg[i]->SetFillColor(0); hcdb_R_neg[i]->SetMarkerColor(1);
    hcdb_P_neg[i]->SetLineColor(1); hcdb_P_neg[i]->SetLineWidth(1); hcdb_P_neg[i]->SetFillColor(0); hcdb_P_neg[i]->SetMarkerColor(1);

    hcdb_P_pos[i]->SetTitle("dphi");
    hcdb_P_neg[i]->SetTitle("dphi");
    hcdb_P_pos[i]->GetYaxis()->SetTitle("dphi (rad)");
    hcdb_P_neg[i]->GetYaxis()->SetTitle("dphi (rad)");
    hcdb_P_pos[i]->GetYaxis()->SetRangeUser(-0.02,0.04);
    hcdb_P_neg[i]->GetYaxis()->SetRangeUser(-0.02,0.04);
    hcdb_R_pos[i]->GetYaxis()->SetRangeUser(-1.0,1.0);
    hcdb_R_neg[i]->GetYaxis()->SetRangeUser(-1.0,1.0);
  }

  // ************************ //


  // ************** Updated central membrane ************** //
  std::string cmfilename[nrun];
  TFile* cmfile[nrun];
  TH2 *hcm_R_pr_pos[nrun], *hcm_P_pr_pos[nrun];
  TH2 *hcm_R_pr_neg[nrun], *hcm_P_pr_neg[nrun];
  TH1 *hcm_R_pos[nrun], *hcm_P_pos[nrun];
  TH1 *hcm_R_neg[nrun], *hcm_P_neg[nrun];

  for (int i = 0; i < nrun; i++)
  {
  
    // ************ If changing from CDB to local files, make sure to set third argument of plot1D to false to avoid multiplying by R again ******************* //

    // cmfile[i] = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationfits/79516_fieldonbeamon/idealtpccoordinates_78242fieldoffbeamoffidealstripereference_updatedreferencestripefits/Laminations_run3pp_new_nocdbtag_v001-00079516.root", "READ");

    // cmfile[i] = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationfits/79516_fieldonbeamon/idealtpccoordinates_83319fieldonbeamoffwupdatedstaticcorrectionsidealstripereference_updatedreferencestripefits_updatedstripefits/Laminations_run3pp_new_nocdbtag_v001-00079516.root", "READ");
    // cmfile[i] = new TFile("/sphenix/user/dloomis/DistortionsAnalysis/laminationfits/79507_fieldonbeamon/idealtpccoordinates_83319fieldonbeamoffwupdatedstaticcorrectionsidealstripereference_updatedreferencestripefits_updatedstripefits/Laminations_run3pp_new_nocdbtag_v001-00079507.root", "READ");

    // just stripe matching for dphi and dR (no lamination fits)
    // cmfile[i] = new TFile("/sphenix/user/dloomis/TrackingAnalysis/CentralMembraneAnalysisv2/macros/CMDistortionCorrections_run3pp_new_nocdbtag_v001-00079507_reference75073fieldonbeamoff.root", "READ");

    // cmfile[i] = new TFile("/sphenix/user/dloomis/TrackingAnalysis/CentralMembraneAnalysisv2/macros/CMDistortionCorrections_run3pp_new_nocdbtag_v001-00083319_staticmodelclosure.root", "READ");
    // cmfile[i] = new TFile("/sphenix/user/dloomis/TrackingAnalysis/CentralMembraneAnalysisv2/macros/CMDistortionCorrections_run3pp_new_nocdbtag_v001-00079507_reference83319fieldonbeamoff_20to40segs.root", "READ");
    
    // cmfile[i] = new TFile("/sphenix/user/dloomis/TrackingAnalysis/CentralMembraneAnalysisv2/macros/output/CMDistortionCorrections-00075073.root", "READ");

    // cmfile[i] = new TFile("/sphenix/user/dloomis/TrackingAnalysis/CentralMembraneAnalysisv2/macros/output/CMDistortionCorrections_run3pp_new_nocdbtag_v001-00083319.root", "READ");
    // cmfile[i] = new TFile("../CentralMembraneStripeMatching/macros/output/CMDistortionCorrections-00079513_reference83319.root", "READ");
    cmfile[i] = new TFile("../CentralMembraneStripeMatching/macros/output/CMDistortionCorrections-00079513_reference75073.root", "READ");
    


    hcm_R_pr_pos[i] = (TH2*) cmfile[i]->Get("hIntDistortionR_posz");
    hcm_P_pr_pos[i] = (TH2*) cmfile[i]->Get("hIntDistortionP_posz");
    hcm_R_pr_neg[i] = (TH2*) cmfile[i]->Get("hIntDistortionR_negz");
    hcm_P_pr_neg[i] = (TH2*) cmfile[i]->Get("hIntDistortionP_negz");

    hcm_R_pos[i] = new TH1F(Form("hIntDistortionR_posz_%d_cm",runs[i]),Form("dR;R (cm);dR (cm)"),hcm_R_pr_pos[i]->GetYaxis()->GetNbins(),hcm_R_pr_pos[i]->GetYaxis()->GetXmin(),hcm_R_pr_pos[i]->GetYaxis()->GetXmax());
    hcm_P_pos[i] = new TH1F(Form("hIntDistortionP_posz_%d_cm",runs[i]),Form("Rdphi;R (cm);Rdphi (cm)"),hcm_P_pr_pos[i]->GetYaxis()->GetNbins(),hcm_P_pr_pos[i]->GetYaxis()->GetXmin(),hcm_P_pr_pos[i]->GetYaxis()->GetXmax());
    hcm_R_neg[i] = new TH1F(Form("hIntDistortionR_negz_%d_cm",runs[i]),Form("dR;R (cm);dR (cm)"),hcm_R_pr_neg[i]->GetYaxis()->GetNbins(),hcm_R_pr_neg[i]->GetYaxis()->GetXmin(),hcm_R_pr_neg[i]->GetYaxis()->GetXmax());
    hcm_P_neg[i] = new TH1F(Form("hIntDistortionP_negz_%d_cm",runs[i]),Form("Rdphi;R (cm);Rdphi (cm)"),hcm_P_pr_neg[i]->GetYaxis()->GetNbins(),hcm_P_pr_neg[i]->GetYaxis()->GetXmin(),hcm_P_pr_neg[i]->GetYaxis()->GetXmax());
    plot1D(hcm_R_pr_pos[i],hcm_R_pos[i], false, false);
    plot1D(hcm_P_pr_pos[i],hcm_P_pos[i], false, false);
    plot1D(hcm_R_pr_neg[i],hcm_R_neg[i], false, false);
    plot1D(hcm_P_pr_neg[i],hcm_P_neg[i], false, false);
    hcm_R_pos[i]->SetLineColor(4); hcm_R_pos[i]->SetLineWidth(1); hcm_R_pos[i]->SetFillColor(0); hcm_R_pos[i]->SetMarkerColor(4);
    hcm_P_pos[i]->SetLineColor(4); hcm_P_pos[i]->SetLineWidth(1); hcm_P_pos[i]->SetFillColor(0); hcm_P_pos[i]->SetMarkerColor(4);
    hcm_R_neg[i]->SetLineColor(4); hcm_R_neg[i]->SetLineWidth(1); hcm_R_neg[i]->SetFillColor(0); hcm_R_neg[i]->SetMarkerColor(4);
    hcm_P_neg[i]->SetLineColor(4); hcm_P_neg[i]->SetLineWidth(1); hcm_P_neg[i]->SetFillColor(0); hcm_P_neg[i]->SetMarkerColor(4);

    hcm_P_pos[i]->SetTitle("dphi");
    hcm_P_neg[i]->SetTitle("dphi");
    hcm_P_pos[i]->GetYaxis()->SetTitle("dphi (rad)");
    hcm_P_neg[i]->GetYaxis()->SetTitle("dphi (rad)");
    hcm_P_pos[i]->GetYaxis()->SetRangeUser(-0.04,0.04);
    hcm_P_neg[i]->GetYaxis()->SetRangeUser(-0.04,0.04);
    // hcm_R_pos[i]->GetYaxis()->SetRangeUser(-2.5,2.5);
    // hcm_R_neg[i]->GetYaxis()->SetRangeUser(-2.5,2.5);
    hcm_R_pos[i]->GetYaxis()->SetRangeUser(-1,1);
    hcm_R_neg[i]->GetYaxis()->SetRangeUser(-1,1);
  }

  
  // ************ Track based (TPOT) distortions ************ //
  // float selectZ = 0.0;

  // TFile* file_3D_map_method3[nrun];
  // TH3 *h_R_prz_pos_method3[nrun], *h_P_prz_pos_method3[nrun];
  // TH3 *h_R_prz_neg_method3[nrun], *h_P_prz_neg_method3[nrun];
  // TH1 *h_R_pos_method3[nrun], *h_P_pos_method3[nrun];
  // TH1 *h_R_neg_method3[nrun], *h_P_neg_method3[nrun];

  // for (int i = 0; i < nrun; i++)
  // {
  //   file_3D_map_method3[i] = new TFile(Form("/sphenix/u/xyu3/workarea/TPCdistortion/Si_TPOT_fit/run3pp_newclusErrorParaTPOT/jobB/Rootfiles/Distortions_full_mm_%d.root",runs[i]),"");
  //   h_R_prz_pos_method3[i] = (TH3*) file_3D_map_method3[i]->Get("hIntDistortionR_posz");
  //   h_P_prz_pos_method3[i] = (TH3*) file_3D_map_method3[i]->Get("hIntDistortionP_posz");
  //   h_R_prz_neg_method3[i] = (TH3*) file_3D_map_method3[i]->Get("hIntDistortionR_negz");
  //   h_P_prz_neg_method3[i] = (TH3*) file_3D_map_method3[i]->Get("hIntDistortionP_negz");

  //   h_R_pos_method3[i] = new TH1F(Form("hIntDistortionR_posz_%d",runs[i]),Form("dR @ Z=%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_pos_method3[i]->GetYaxis()->GetNbins(),h_R_prz_pos_method3[i]->GetYaxis()->GetXmin(),h_R_prz_pos_method3[i]->GetYaxis()->GetXmax());
  //   h_P_pos_method3[i] = new TH1F(Form("hIntDistortionP_posz_%d",runs[i]),Form("Rdphi @ Z=%d cm;R (cm);Rdphi (cm)",(int)selectZ),h_P_prz_pos_method3[i]->GetYaxis()->GetNbins(),h_P_prz_pos_method3[i]->GetYaxis()->GetXmin(),h_P_prz_pos_method3[i]->GetYaxis()->GetXmax());
    
  //   h_R_neg_method3[i] = new TH1F(Form("hIntDistortionR_negz_%d",runs[i]),Form("dR @ Z=-%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_neg_method3[i]->GetYaxis()->GetNbins(),h_R_prz_neg_method3[i]->GetYaxis()->GetXmin(),h_R_prz_neg_method3[i]->GetYaxis()->GetXmax());
  //   h_P_neg_method3[i] = new TH1F(Form("hIntDistortionP_negz_%d",runs[i]),Form("Rdphi @ Z=-%d cm;R (cm);Rdphi (cm)",(int)selectZ),h_P_prz_neg_method3[i]->GetYaxis()->GetNbins(),h_P_prz_neg_method3[i]->GetYaxis()->GetXmin(),h_P_prz_neg_method3[i]->GetYaxis()->GetXmax());
    
  //   plot1D_Zbin(h_R_prz_pos_method3[i],h_R_pos_method3[i],selectZ, false);
  //   plot1D_Zbin(h_P_prz_pos_method3[i],h_P_pos_method3[i],selectZ, false);
  //   plot1D_Zbin(h_R_prz_neg_method3[i],h_R_neg_method3[i],-selectZ, false);
  //   plot1D_Zbin(h_P_prz_neg_method3[i],h_P_neg_method3[i],-selectZ, false);
  //   h_R_pos_method3[i]->SetLineColor((3)); h_R_pos_method3[i]->SetLineWidth(1); h_R_pos_method3[i]->SetFillColor(0); h_R_pos_method3[i]->SetMarkerColor((3));
  //   h_P_pos_method3[i]->SetLineColor((3)); h_P_pos_method3[i]->SetLineWidth(1); h_P_pos_method3[i]->SetFillColor(0); h_P_pos_method3[i]->SetMarkerColor((3));
  //   h_R_neg_method3[i]->SetLineColor((3)); h_R_neg_method3[i]->SetLineWidth(1); h_R_neg_method3[i]->SetFillColor(0); h_R_neg_method3[i]->SetMarkerColor((3));
  //   h_P_neg_method3[i]->SetLineColor((3)); h_P_neg_method3[i]->SetLineWidth(1); h_P_neg_method3[i]->SetFillColor(0); h_P_neg_method3[i]->SetMarkerColor((3));

  //   h_P_pos_method3[i]->GetYaxis()->SetRangeUser(-0.01,0.04);
  //   h_P_neg_method3[i]->GetYaxis()->SetRangeUser(-0.01,0.04);

  // }


  // ************************ //


  // ************ Digital Current distortions ************ //
  float selectZ = 0.0;
  TFile* file_digitalcurrents_spacecharge[nrun];

  TH3F *h_R_prz_pos_DC[nrun], *h_P_prz_pos_DC[nrun];
  TH3F *h_R_prz_neg_DC[nrun], *h_P_prz_neg_DC[nrun];

  TH1 *h_R_pos_DC[nrun], *h_P_pos_DC[nrun];
  TH1 *h_R_neg_DC[nrun], *h_P_neg_DC[nrun];

  for (int i=0; i<nrun; i++)
  {

    // Dhanush's digital current maps with space charge distortions and static distortions subtracted out
    file_digitalcurrents_spacecharge[i] = new TFile(Form("/sphenix/user/dloomis/Distortions/DistortionMaps_DigitalCurrent/distortion_spacecharge_run%d_mean.distortion_map.hist.root",runs[i]),"");
    // file_digitalcurrents_spacecharge[i] = new TFile(Form("/sphenix/user/dloomis/Distortions/DistortionMaps_DigitalCurrent/distortion_static_spacecharge_run%d_mean.distortion_map.hist.root",runs[i]),"");
    
    // These are static distortions from Ross's model
    // file_digitalcurrents_spacecharge[i] = new TFile("/sphenix/user/dloomis/Distortions/DistortionMaps_Static/static_only_inverted_10-new_hugoupdate.root");
    // file_digitalcurrents_spacecharge[i] = new TFile("/sphenix/user/dloomis/Distortions/DistortionMaps_Static/GarfieldStaticCorrectionMaps.root");
    

    // file_digitalcurrents_spacecharge[i] = new TFile("/sphenix/user/dloomis/TrackingAnalysis/CentralMembraneAnalysisv2/macros/CMDistortionCorrections_run3pp_new_nocdbtag_v001-00079507_reference83319fieldonbeamoff_20to40segs.root");
    
    h_R_prz_pos_DC[i] = (TH3F*) file_digitalcurrents_spacecharge[i]->Get("hIntDistortionR_posz");
    h_P_prz_pos_DC[i] = (TH3F*) file_digitalcurrents_spacecharge[i]->Get("hIntDistortionP_posz");
    h_R_prz_neg_DC[i] = (TH3F*) file_digitalcurrents_spacecharge[i]->Get("hIntDistortionR_negz");
    h_P_prz_neg_DC[i] = (TH3F*) file_digitalcurrents_spacecharge[i]->Get("hIntDistortionP_negz");


    h_R_pos_DC[i] = new TH1F(Form("hIntDistortionR_posz_%d",runs[i]),Form("dR @ Z=%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_pos_DC[i]->GetYaxis()->GetNbins(),h_R_prz_pos_DC[i]->GetYaxis()->GetXmin(),h_R_prz_pos_DC[i]->GetYaxis()->GetXmax());
    h_P_pos_DC[i] = new TH1F(Form("hIntDistortionP_posz_%d",runs[i]),Form("Rdphi @ Z=%d cm;R (cm);Rdphi (cm)",(int)selectZ),h_P_prz_pos_DC[i]->GetYaxis()->GetNbins(),h_P_prz_pos_DC[i]->GetYaxis()->GetXmin(),h_P_prz_pos_DC[i]->GetYaxis()->GetXmax());
    
    h_R_neg_DC[i] = new TH1F(Form("hIntDistortionR_negz_%d",runs[i]),Form("dR @ Z=-%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_neg_DC[i]->GetYaxis()->GetNbins(),h_R_prz_neg_DC[i]->GetYaxis()->GetXmin(),h_R_prz_neg_DC[i]->GetYaxis()->GetXmax());
    h_P_neg_DC[i] = new TH1F(Form("hIntDistortionP_negz_%d",runs[i]),Form("Rdphi @ Z=-%d cm;R (cm);Rdphi (cm)",(int)selectZ),h_P_prz_neg_DC[i]->GetYaxis()->GetNbins(),h_P_prz_neg_DC[i]->GetYaxis()->GetXmin(),h_P_prz_neg_DC[i]->GetYaxis()->GetXmax());

    plot1D_Zbin(h_R_prz_pos_DC[i],h_R_pos_DC[i],selectZ, false, false);
    // plot1D_Zbin(h_P_prz_pos_DC[i],h_P_pos_DC[i],selectZ, false, true);
    plot1D_Zbin(h_P_prz_pos_DC[i],h_P_pos_DC[i],selectZ, false, false);
    plot1D_Zbin(h_R_prz_neg_DC[i],h_R_neg_DC[i],-selectZ, false, false);
    // plot1D_Zbin(h_P_prz_neg_DC[i],h_P_neg_DC[i],-selectZ, false, true);
    plot1D_Zbin(h_P_prz_neg_DC[i],h_P_neg_DC[i],-selectZ, false, false);
    


    h_P_pos_DC[i]->SetLineColor((2)); h_P_pos_DC[i]->SetLineWidth(1); h_P_pos_DC[i]->SetFillColor(0); h_P_pos_DC[i]->SetMarkerColor((2));
    h_R_pos_DC[i]->SetLineColor((2)); h_R_pos_DC[i]->SetLineWidth(1); h_R_pos_DC[i]->SetFillColor(0); h_R_pos_DC[i]->SetMarkerColor((2));
    h_P_neg_DC[i]->SetLineColor((2)); h_P_neg_DC[i]->SetLineWidth(1); h_P_neg_DC[i]->SetFillColor(0); h_P_neg_DC[i]->SetMarkerColor((2));
    h_R_neg_DC[i]->SetLineColor((2)); h_R_neg_DC[i]->SetLineWidth(1); h_R_neg_DC[i]->SetFillColor(0); h_R_neg_DC[i]->SetMarkerColor((2));

    // h_P_pos_DC[i]->GetYaxis()->SetRangeUser(-0.04,0.04);
    // h_P_neg_DC[i]->GetYaxis()->SetRangeUser(-0.04,0.04);


    // std::pair<double,double> yrange_R = SetCommonYRange({h_R_neg_method3[i], h_R_pos_method3[i], hcdb_R_neg[i], hcdb_R_pos[i]});
    // std::pair<double,double> yrange_P = SetCommonYRange({h_P_neg_method3[i], h_P_pos_method3[i], hcdb_P_neg[i], hcdb_P_pos[i]});

    // std::pair<double,double> yrange_R = SetCommonYRange({h_R_neg_DC[i], h_R_pos_DC[i], hcdb_R_neg[i], hcdb_R_pos[i], h_R_neg_method3[i], h_R_pos_method3[i]});
    // std::pair<double,double> yrange_P = SetCommonYRange({h_P_neg_DC[i], h_P_pos_DC[i], hcdb_P_neg[i], hcdb_P_pos[i], h_R_neg_method3[i], h_R_pos_method3[i]});


  }



  // ************ Static model distortions ************ //
  TFile* file_model_static[nrun];

  TH3F *h_R_prz_pos_staticmodel[nrun], *h_P_prz_pos_staticmodel[nrun];
  TH3F *h_R_prz_neg_staticmodel[nrun], *h_P_prz_neg_staticmodel[nrun];

  TH1 *h_R_pos_staticmodel[nrun], *h_P_pos_staticmodel[nrun];
  TH1 *h_R_neg_staticmodel[nrun], *h_P_neg_staticmodel[nrun];

  for (int i=0; i<nrun; i++)
  {

    // These are static distortions from Ross's model
    // file_digitalcurrents_spacecharge[i] = new TFile("/sphenix/user/dloomis/Distortions/DistortionMaps_Static/static_only_inverted_10-new_hugoupdate.root");
    // file_model_static[i] = new TFile("/sphenix/user/dloomis/Distortions/DistortionMaps_Static/GarfieldStaticCorrectionMaps.root");
    file_model_static[i] = new TFile("/sphenix/user/dloomis/Distortions/DistortionMaps_Static/GarfieldStaticCorrectionMaps_Aug_tpcReal_magReal.root");
    
    h_R_prz_pos_staticmodel[i] = (TH3F*) file_model_static[i]->Get("hIntDistortionR_posz");
    h_P_prz_pos_staticmodel[i] = (TH3F*) file_model_static[i]->Get("hIntDistortionP_posz");
    h_R_prz_neg_staticmodel[i] = (TH3F*) file_model_static[i]->Get("hIntDistortionR_negz");
    h_P_prz_neg_staticmodel[i] = (TH3F*) file_model_static[i]->Get("hIntDistortionP_negz");


    h_R_pos_staticmodel[i] = new TH1F(Form("hIntDistortionR_posz_%d",runs[i]),Form("dR @ Z=%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_pos_staticmodel[i]->GetYaxis()->GetNbins(),h_R_prz_pos_staticmodel[i]->GetYaxis()->GetXmin(),h_R_prz_pos_staticmodel[i]->GetYaxis()->GetXmax());
    h_P_pos_staticmodel[i] = new TH1F(Form("hIntDistortionP_posz_%d",runs[i]),Form("Rdphi @ Z=%d cm;R (cm);Rdphi (cm)",(int)selectZ),h_P_prz_pos_staticmodel[i]->GetYaxis()->GetNbins(),h_P_prz_pos_staticmodel[i]->GetYaxis()->GetXmin(),h_P_prz_pos_staticmodel[i]->GetYaxis()->GetXmax());
    
    h_R_neg_staticmodel[i] = new TH1F(Form("hIntDistortionR_negz_%d",runs[i]),Form("dR @ Z=-%d cm;R (cm);dR (cm)",(int)selectZ),h_R_prz_neg_staticmodel[i]->GetYaxis()->GetNbins(),h_R_prz_neg_staticmodel[i]->GetYaxis()->GetXmin(),h_R_prz_neg_staticmodel[i]->GetYaxis()->GetXmax());
    h_P_neg_staticmodel[i] = new TH1F(Form("hIntDistortionP_negz_%d",runs[i]),Form("Rdphi @ Z=-%d cm;R (cm);Rdphi (cm)",(int)selectZ),h_P_prz_neg_staticmodel[i]->GetYaxis()->GetNbins(),h_P_prz_neg_staticmodel[i]->GetYaxis()->GetXmin(),h_P_prz_neg_staticmodel[i]->GetYaxis()->GetXmax());

    plot1D_Zbin(h_R_prz_pos_staticmodel[i],h_R_pos_staticmodel[i],selectZ, false, false);
    // plot1D_Zbin(h_P_prz_pos_staticmodel[i],h_P_pos_staticmodel[i],selectZ, false, true);
    plot1D_Zbin(h_P_prz_pos_staticmodel[i],h_P_pos_staticmodel[i],selectZ, false, false);
    plot1D_Zbin(h_R_prz_neg_staticmodel[i],h_R_neg_staticmodel[i],-selectZ, false, false);
    // plot1D_Zbin(h_P_prz_neg_staticmodel[i],h_P_neg_staticmodel[i],-selectZ, false, true);
    plot1D_Zbin(h_P_prz_neg_staticmodel[i],h_P_neg_staticmodel[i],-selectZ, false, false);
    


    h_P_pos_staticmodel[i]->SetLineColor((3)); h_P_pos_staticmodel[i]->SetLineWidth(1); h_P_pos_staticmodel[i]->SetFillColor(0); h_P_pos_staticmodel[i]->SetMarkerColor((3));
    h_R_pos_staticmodel[i]->SetLineColor((3)); h_R_pos_staticmodel[i]->SetLineWidth(1); h_R_pos_staticmodel[i]->SetFillColor(0); h_R_pos_staticmodel[i]->SetMarkerColor((3));
    h_P_neg_staticmodel[i]->SetLineColor((3)); h_P_neg_staticmodel[i]->SetLineWidth(1); h_P_neg_staticmodel[i]->SetFillColor(0); h_P_neg_staticmodel[i]->SetMarkerColor((3));
    h_R_neg_staticmodel[i]->SetLineColor((3)); h_R_neg_staticmodel[i]->SetLineWidth(1); h_R_neg_staticmodel[i]->SetFillColor(0); h_R_neg_staticmodel[i]->SetMarkerColor((3));

    // h_P_pos_staticmodel[i]->GetYaxis()->SetRangeUser(-0.04,0.04);
    // h_P_neg_staticmodel[i]->GetYaxis()->SetRangeUser(-0.04,0.04);


    // std::pair<double,double> yrange_R = SetCommonYRange({h_R_neg_method3[i], h_R_pos_method3[i], hcdb_R_neg[i], hcdb_R_pos[i]});
    // std::pair<double,double> yrange_P = SetCommonYRange({h_P_neg_method3[i], h_P_pos_method3[i], hcdb_P_neg[i], hcdb_P_pos[i]});

    // std::pair<double,double> yrange_R = SetCommonYRange({h_R_neg_DC[i], h_R_pos_DC[i], hcdb_R_neg[i], hcdb_R_pos[i], h_R_neg_method3[i], h_R_pos_method3[i]});
    // std::pair<double,double> yrange_P = SetCommonYRange({h_P_neg_DC[i], h_P_pos_DC[i], hcdb_P_neg[i], hcdb_P_pos[i], h_R_neg_method3[i], h_R_pos_method3[i]});


  }


  for (int i=0; i<nrun; i++)
  {
    TCanvas* can = new TCanvas("can","",1200,1200);
    can->Divide(2,2);
    can->cd(1);
    gPad->SetLogy(0);
    // hcdb_P_pos[i]->Draw("hist,same");
    hcm_P_pos[i]->Draw("hist,same");
    // h_P_pos_method3[i]->Draw("hist,e,same");
    h_P_pos_DC[i]->Draw("hist,same");
    // h_P_pos_staticmodel[i]->Draw("hist,same");

    TLegend *legend = new TLegend(0.45, 0.65, 0.8, 0.9);
    // legend->SetHeader(Form("Run %d",runs[i]));
    legend->SetHeader("Run 79513");
    // legend->SetHeader("Static Distortions");
    legend->SetTextSize(0.035);
    // legend->AddEntry(h_P_pos_method3[i], Form("TPOT z=5 cm"), "l");
    // legend->AddEntry(hcdb_P_pos[i], Form("CDB"), "l");
    legend->AddEntry(hcm_P_pos[i], Form("CM Stripe Matching"), "l");
    legend->AddEntry(h_P_pos_DC[i], Form("Digital Current - 79516"), "l");
    // legend->AddEntry(h_P_pos_DC[i], Form("Static Model"), "l");
    // legend->AddEntry(h_P_pos_staticmodel[i], Form("Static Model (w/ tilt)"), "l");

    legend->Draw();   

    can->cd(2);
    gPad->SetLogy(0);
    // hcdb_R_pos[i]->Draw("hist,same");
    hcm_R_pos[i]->Draw("hist,same");
    // h_R_pos_method3[i]->Draw("hist,e,same");
    h_R_pos_DC[i]->Draw("hist,same");
    // h_R_pos_staticmodel[i]->Draw("hist,same");
    can->cd(3);
    gPad->SetLogy(0);
    // hcdb_P_neg[i]->Draw("hist,same");
    hcm_P_neg[i]->Draw("hist,same");
    // h_P_neg_method3[i]->Draw("hist,e,same");
    h_P_neg_DC[i]->Draw("hist,same");
    // h_P_neg_staticmodel[i]->Draw("hist,same");
    can->cd(4);
    gPad->SetLogy(0);
    // hcdb_R_neg[i]->Draw("hist,same");
    hcm_R_neg[i]->Draw("hist,same");
    // h_R_neg_method3[i]->Draw("hist,e,same");
    h_R_neg_DC[i]->Draw("hist,same");
    // h_R_neg_staticmodel[i]->Draw("hist,same");


    gPad->RedrawAxis();

    can->Update();
    can->SaveAs(Form("79507/figure/resid_vsR_from3D_atZ%d_%d_ref75073.pdf",(int)selectZ,runs[i]));

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
