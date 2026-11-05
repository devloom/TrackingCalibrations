#include <iostream>
#include <string>
#include <utility>
#include <vector>

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

  struct CMPlotConfig
  {
    std::string label;
    std::string filename;
    int color;
  };

  struct CMPlot
  {
    TFile* file = nullptr;
    TH2* R_pr_pos = nullptr;
    TH2* P_pr_pos = nullptr;
    TH2* R_pr_neg = nullptr;
    TH2* P_pr_neg = nullptr;
    TH1* R_pos = nullptr;
    TH1* P_pos = nullptr;
    TH1* R_neg = nullptr;
    TH1* P_neg = nullptr;
  };

  std::vector<CMPlotConfig> cmPlots = {
    {
      "79507",
      "../CentralMembraneStripeMatching/macros/output/CMDistortionCorrections-00079507_reference83319.root",
      kBlack
    },
    {
      "79508",
      "../CentralMembraneStripeMatching/macros/output/CMDistortionCorrections-00079508_reference83319.root",
      kRed + 1
    },
    {
      "79509",
      "../CentralMembraneStripeMatching/macros/output/CMDistortionCorrections-00079509_reference83319.root",
      kOrange + 1
    },
    {
      "79510",
      "../CentralMembraneStripeMatching/macros/output/CMDistortionCorrections-00079510_reference83319.root",
      kYellow + 1
    },
    {
      "79511",
      "../CentralMembraneStripeMatching/macros/output/CMDistortionCorrections-00079511_reference83319.root",
      kMagenta + 1
    },
    {
      "79512",
      "../CentralMembraneStripeMatching/macros/output/CMDistortionCorrections-00079512_reference83319.root",
      kGreen + 1
    },
    {
      "79513",
      "../CentralMembraneStripeMatching/macros/output/CMDistortionCorrections-00079513_reference83319.root",
      kCyan + 1
    },
    {
      "79514",
      "../CentralMembraneStripeMatching/macros/output/CMDistortionCorrections-00079514_reference83319.root",
      kTeal + 1
    },
    {
      "79515",
      "../CentralMembraneStripeMatching/macros/output/CMDistortionCorrections-00079515_reference83319.root",
      kAzure + 1
    },
    {
      "79516",
      "../CentralMembraneStripeMatching/macros/output/CMDistortionCorrections-00079516_reference83319.root",
      kBlue + 1
    },
  };

  std::vector<CMPlot> cmHists(cmPlots.size());

  for (size_t icm = 0; icm < cmPlots.size(); icm++)
  {
    CMPlot& plot = cmHists[icm];
    const CMPlotConfig& config = cmPlots[icm];

    plot.file = new TFile(config.filename.c_str(), "READ");
    if (!plot.file || plot.file->IsZombie())
    {
      std::cout << "Could not open CM file: " << config.filename << std::endl;
      continue;
    }

    plot.R_pr_pos = (TH2*) plot.file->Get("hIntDistortionR_posz");
    plot.P_pr_pos = (TH2*) plot.file->Get("hIntDistortionP_posz");
    plot.R_pr_neg = (TH2*) plot.file->Get("hIntDistortionR_negz");
    plot.P_pr_neg = (TH2*) plot.file->Get("hIntDistortionP_negz");

    if (!plot.R_pr_pos || !plot.P_pr_pos || !plot.R_pr_neg || !plot.P_pr_neg)
    {
      std::cout << "Missing one or more CM histograms in: " << config.filename << std::endl;
      continue;
    }

    plot.R_pos = new TH1F(Form("hIntDistortionR_posz_cm_%zu", icm),
                          "dR;R (cm);dR (cm)",
                          plot.R_pr_pos->GetYaxis()->GetNbins(),
                          plot.R_pr_pos->GetYaxis()->GetXmin(),
                          plot.R_pr_pos->GetYaxis()->GetXmax());
    plot.P_pos = new TH1F(Form("hIntDistortionP_posz_cm_%zu", icm),
                          "dphi;R (cm);dphi (rad)",
                          plot.P_pr_pos->GetYaxis()->GetNbins(),
                          plot.P_pr_pos->GetYaxis()->GetXmin(),
                          plot.P_pr_pos->GetYaxis()->GetXmax());
    plot.R_neg = new TH1F(Form("hIntDistortionR_negz_cm_%zu", icm),
                          "dR;R (cm);dR (cm)",
                          plot.R_pr_neg->GetYaxis()->GetNbins(),
                          plot.R_pr_neg->GetYaxis()->GetXmin(),
                          plot.R_pr_neg->GetYaxis()->GetXmax());
    plot.P_neg = new TH1F(Form("hIntDistortionP_negz_cm_%zu", icm),
                          "dphi;R (cm);dphi (rad)",
                          plot.P_pr_neg->GetYaxis()->GetNbins(),
                          plot.P_pr_neg->GetYaxis()->GetXmin(),
                          plot.P_pr_neg->GetYaxis()->GetXmax());

    plot1D(plot.R_pr_pos, plot.R_pos, false, false);
    plot1D(plot.P_pr_pos, plot.P_pos, false, false);
    plot1D(plot.R_pr_neg, plot.R_neg, false, false);
    plot1D(plot.P_pr_neg, plot.P_neg, false, false);

    std::vector<TH1*> hists = {plot.R_pos, plot.P_pos, plot.R_neg, plot.P_neg};
    for (TH1* hist : hists)
    {
      hist->SetLineColor(config.color);
      hist->SetLineWidth(2);
      hist->SetFillColor(0);
      hist->SetMarkerColor(config.color);
    }
  }

  std::vector<TH1*> pPosHists;
  std::vector<TH1*> rPosHists;
  std::vector<TH1*> pNegHists;
  std::vector<TH1*> rNegHists;

  for (CMPlot& plot : cmHists)
  {
    if (!plot.P_pos) continue;
    pPosHists.push_back(plot.P_pos);
    rPosHists.push_back(plot.R_pos);
    pNegHists.push_back(plot.P_neg);
    rNegHists.push_back(plot.R_neg);
  }

  SetCommonYRange(pPosHists);
  SetCommonYRange(rPosHists);
  SetCommonYRange(pNegHists);
  SetCommonYRange(rNegHists);

  TCanvas* can = new TCanvas("can_compareCMplots", "", 1200, 1200);
  can->Divide(2, 2);

  TLegend *legend = new TLegend(0.6, 0.5, 0.9, 0.9);
//   legend->SetEntrySeparation(0.05);
  legend->SetTextSize(0.05);

  for (size_t icm = 0; icm < cmHists.size(); icm++)
  {
    CMPlot& plot = cmHists[icm];
    if (plot.P_pos) legend->AddEntry(plot.P_pos, cmPlots[icm].label.c_str(), "l");
  }

  can->cd(1);
  gPad->SetLogy(0);
  for (size_t ihist = 0; ihist < pPosHists.size(); ihist++)
  {
    pPosHists[ihist]->Draw(ihist == 0 ? "hist" : "hist,same");
  }
  legend->Draw();

  can->cd(2);
  gPad->SetLogy(0);
  for (size_t ihist = 0; ihist < rPosHists.size(); ihist++)
  {
    rPosHists[ihist]->Draw(ihist == 0 ? "hist" : "hist,same");
  }

  can->cd(3);
  gPad->SetLogy(0);
  for (size_t ihist = 0; ihist < pNegHists.size(); ihist++)
  {
    pNegHists[ihist]->Draw(ihist == 0 ? "hist" : "hist,same");
  }

  can->cd(4);
  gPad->SetLogy(0);
  for (size_t ihist = 0; ihist < rNegHists.size(); ihist++)
  {
    rNegHists[ihist]->Draw(ihist == 0 ? "hist" : "hist,same");
  }

  for (int ipad = 1; ipad <= 4; ipad++)
  {
    can->cd(ipad);
    gPad->RedrawAxis();
  }

  can->Update();
  can->SaveAs("79507/figure/resid_vsR_compareCMplots_ref83319.pdf");

}

void drawMultipleCM()
{
  drawCMvsTPOTvsDC();
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
            h->SetMinimum(yMin-(abs(yMin)*0.5));
            h->SetMaximum(yMax+(abs(yMax)*0.5));
        }
    }
    return {yMin, yMax};
}
