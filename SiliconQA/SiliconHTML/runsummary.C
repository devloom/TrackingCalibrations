////// Taken from Hao-Ren: /sphenix/user/hjheng/sPHENIXRepo/MVTX/runcheck/runsummary.C
/////  For testing the implementation of raw hit MVTX QA with previously implemented cluster QA


#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <sstream>
#include <string>
#include <boost/algorithm/string.hpp>
#include <memory> // For smart pointers

constexpr int kMarkerStyle = 25;
constexpr float kMarkerSize = 0.6f;

std::map<std::string, std::pair<int, int>> map_beamcrossing_run()
{
    std::map<std::string, std::pair<int, int>> m_beamcrossing_run;
    m_beamcrossing_run["0 mrad"] = std::make_pair(46619, 49034);
    m_beamcrossing_run["+1.5 mrad"] = std::make_pair(49034, 49054);
    m_beamcrossing_run["0 mrad"] = std::make_pair(49054, 51242);
    m_beamcrossing_run["+1.5 mrad"] = std::make_pair(51242, std::numeric_limits<int>::quiet_NaN());

    return m_beamcrossing_run;
}

std::tuple<float, float, float, float> fitClusPhi(TH1F *hm)
{
    // get the average number of clusters per bin
    float avgc = hm->Integral() / hm->GetNbinsX();
    TF1 *f1 = new TF1("f1", "pol0(0)+[1]*cos(x+[2])", -3.14159, 3.14159);
    // set the initial parameters
    f1->SetParameter(0, avgc);
    f1->SetParameter(1, avgc * 0.07);
    f1->SetParameter(2, 9.5);
    f1->SetParLimits(0, avgc*0.5, avgc*1.5);
    f1->SetParLimits(1, avgc*0.001, avgc*0.5);
    f1->SetParLimits(2, 9.3, 9.7);
    f1->SetLineColor(kRed);
    hm->Fit("f1", "R P M");
    return std::make_tuple(f1->GetParameter(0), f1->GetParameter(1), f1->GetParameter(2), (f1->GetChisquare() / f1->GetParameter(0)) / f1->GetNDF());
}

float rawHitAcceptance(TH2F *hlay)
{
    int totalchips = hlay->GetNbinsX();
    int totalstaves = hlay->GetNbinsY();
    int emptybins = 0;

    for (int x = 0; x < totalchips; x++)
    {
      for (int y = 0; y < totalstaves; y++)
      {
        if (hlay->GetBinContent(x+1,y+1) == 0)
        {
          emptybins++;
        }
      }
    }

    return (float(emptybins)/(totalchips*totalstaves));
}

void makesummarytable(std::map<std::string, std::vector<int>> map_runstatuslist)
{
    ofstream out("runstatus_table_v2.tex");
    out << "\\documentclass{article}" << endl;
    out << "\\usepackage{hyperref}" << endl;
    out << "\\usepackage{graphicx}" << endl;
    out << "\\usepackage[dvipsnames]{xcolor}" << endl;
    out << "\\usepackage[yyyymmdd,hhmmss]{datetime}" << endl;
    out << "\\usepackage{geometry}" << endl;
    out << "\\geometry{a4paper, landscape, margin=0.5in, top=0.5in}" << endl;
    out << "\\pagenumbering{gobble}" << endl;
    out << "\\makeatletter" << endl;
    out << "\\title{MVTX Run Summary}" << endl;
    out << "\\author{Hao-Ren Jheng}" << endl;
    out << "\\date{\\today} " << endl;
    out << "\\makeatother" << endl;

    out << "\\begin{document}" << endl;
    out << "\\maketitle" << endl;
    // Create the table first
    out << "\\begin{table}[h!]" << endl;
    out << "\\centering" << endl;
    out << "\\begin{tabular}{|c|c|c|}" << endl;
    out << "\\hline" << endl;
    out << "Run Status & Number of Runs & Fraction \\\\" << endl;
    out << "\\hline" << endl;
    out << "Total number of runs & "
        << map_runstatuslist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}"].size() + map_runstatuslist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}"].size() +
               map_runstatuslist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}"].size() + map_runstatuslist["\\textcolor{Purple}{No Cluster, Hit QA}"].size() + 
               map_runstatuslist["Raw hit acceptance: \\textcolor{Red}{BAD}"].size()
        << " & - \\\\" << endl;
    out << "\\hline" << endl;
    for (const auto &elem : map_runstatuslist)
    {
        if (elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}" && 
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}" && 
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}" &&
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}; Run length: \\textcolor{Orange}{Short}" &&
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}; Run length: \\textcolor{ForestGreen}{Normal}" &&
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}; Run length: \\textcolor{Orange}{Short}" &&
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}; Run length: \\textcolor{ForestGreen}{Normal}" &&
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}; Run length: \\textcolor{Orange}{Short}" &&
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}; Run length: \\textcolor{ForestGreen}{Normal}")
        {
            float frac =
                (float)elem.second.size() / (map_runstatuslist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}"].size() + 
                                             map_runstatuslist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}"].size() +
                                             map_runstatuslist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}"].size() + 
                                             map_runstatuslist["\\textcolor{Purple}{No Cluster, Hit QA}"].size() +
                                             map_runstatuslist["Raw hit acceptance: \\textcolor{Red}{BAD}"].size());
            // Format the fraction to percentage and 3 significant figures
            stringstream ss;
            ss << setprecision(3) << frac * 100.;
            string fracstr = ss.str();
            out << elem.first << " & " << elem.second.size() << " & " << fracstr << "\\% \\\\" << endl;
        }
    }
    out << "\\hline" << endl;
    out << "\\end{tabular}" << endl;
    out << "\\end{table}" << endl;
    out << endl;

    // bullet list, each item includes status, number of runs, and list of runs with proper line break
    out << "\\begin{itemize}" << endl;
    for (const auto &elem : map_runstatuslist)
    {
        if (elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}" && 
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}" && 
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}" &&
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}; Run length: \\textcolor{Orange}{Short}" &&
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}; Run length: \\textcolor{ForestGreen}{Normal}" &&
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}; Run length: \\textcolor{Orange}{Short}" &&
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}; Run length: \\textcolor{ForestGreen}{Normal}" &&
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}; Run length: \\textcolor{Orange}{Short}" &&
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}; Run length: \\textcolor{ForestGreen}{Normal}")
        {
            out << "\\item \\textbf{" << elem.first << "} \\newline" << endl;
            out << "Run list: ";
            for (const auto &run : elem.second)
            {
                if (run != elem.second.back())
                    out << run << ", ";
                else
                    out << run << endl;
            }
            out << endl;
        }
    }
    out << "\\end{itemize}" << endl;

    out << "\\begin{figure}[p]" << endl;
    out << "\\centering" << endl;
    out << "\\includegraphics[width=0.99\\textwidth]{Chi2NdfGraph.pdf}" << endl;
    out << "\\caption{$\\log_{10}$ of scaled $\\chi^{2}/$ndf}" << endl;
    out << "\\end{figure}" << endl;

    out << "\\begin{figure}[p]" << endl;
    out << "\\centering" << endl;
    out << "\\includegraphics[width=0.99\\textwidth]{AvgnclusGraph.pdf}" << endl;
    out << "\\caption{Average number of clusters per bin (constant offset A)}" << endl;
    out << "\\end{figure}" << endl;

    out << "\\begin{figure}[p]" << endl;
    out << "\\centering" << endl;
    out << "\\includegraphics[width=0.99\\textwidth]{BoverAGraph.pdf}" << endl;
    out << "\\caption{Modulation amplitude ratio (B/A)}" << endl;
    out << "\\end{figure}" << endl;

    out << "\\begin{figure}[p]" << endl;
    out << "\\centering" << endl;
    out << "\\includegraphics[width=0.99\\textwidth]{ChipOccupancyOverflowGraph.pdf}" << endl;
    out << "\\caption{$\\log_{10}$ of probability of chip occupancy $\\geq$ 0.3\\%)}" << endl;
    out << "\\end{figure}" << endl;

    //out << "\\newpage" << endl;
    //out << "Plots from individual runs can be found in \\href{https://sphenix-intra.sdcc.bnl.gov/WWW/user/hjheng/Run2024/MvtxQA-Summary/}{https://sphenix-intra.sdcc.bnl.gov/WWW/user/hjheng/Run2024/MvtxQA-Summary/}." << endl;

    out << "\\end{document}" << endl;
    out.close();
}

std::string modifystr4html(std::string s)
{
    // replace the latex $\chi^{2}$ with html code
    boost::replace_all(s, "$\\chi^{2}$", "&chi;<sup>2</sup>");
    // replace the latex \textcolor{ForestGreen}{GOOD} with html code
    boost::replace_all(s, "\\textcolor{ForestGreen}{GOOD}", "<font color=\"#11cf17\">GOOD</font>");
    // replace the latex \textcolor{Orange}{OK} with html code
    boost::replace_all(s, "\\textcolor{Orange}{OK}", "<font color=\"#ff7f00\">OK</font>");
    // replace the latex \textcolor{Red}{BAD} with html code
    boost::replace_all(s, "\\textcolor{Red}{BAD}", "<font color=\"#ff0000\">BAD</font>");
    // replace the latex \textcolor{Purple}{No Cluster QA} with html code
    boost::replace_all(s, "\\textcolor{Purple}{No Cluster QA}", "<font color=\"#800080\">No Cluster QA</font>");
    // replace the latex \textcolor{Purple}{No Cluster, Hit QA} with html code
    boost::replace_all(s, "\\textcolor{Purple}{No Cluster, Hit QA}", "<font color=\"#800080\">No Cluster, Hit QA</font>");
    // replace the latex \textcolor{Purple}{No Hit QA} with html code
    boost::replace_all(s, "\\textcolor{Purple}{No Hit QA}", "<font color=\"#800080\">No Hit QA</font>");
    // replace the latex \textcolor{Orange}{Short} with html code
    boost::replace_all(s, "\\textcolor{Orange}{Short}", "<font color=\"#ff7f00\">Short</font>");
    // replace the latex \textcolor{ForestGreen}{Normal} with html code
    boost::replace_all(s, "\\textcolor{ForestGreen}{Normal}", "<font color=\"#11cf17\">Normal</font>");
    // replace B/A to modulation amplitude-to-average ratio
    boost::replace_all(s, "B/A ratio", "modulation amplitude-to-average (B/A) ratio");

    return s;
}

void makehtml(std::map<std::string, std::vector<int>> map_runstatuslist)
{
    ofstream out("runstatus_v2.html");
    out << "<!DOCTYPE html>\n";
    out << "<html lang=\"en\">\n";
    out << "<head>\n";
    out << "<meta charset=\"UTF-8\">\n";
    out << "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n";
    out << "<title>MVTX Run Summary</title>\n";
    out << "<style>\n";
    out << "  body {\n";
    out << "    font-family: 'Arial', sans-serif;\n";
    out << "    margin: 0;\n";
    out << "    padding: 20px;\n";
    out << "    background-color: #f4f7f9;\n";
    out << "    color: #333;\n";
    out << "  }\n";
    out << "  h2 {\n";
    out << "    text-align: center;\n";
    out << "    font-size: 2em;\n";
    out << "    color: #4CAF50;\n";
    out << "    margin-bottom: 20px;\n";
    out << "  }\n";
    out << "  button {\n";
    out << "    background-color: #4CAF50;\n";
    out << "    color: white;\n";
    out << "    padding: 10px 20px;\n";
    out << "    border: none;\n";
    out << "    border-radius: 5px;\n";
    out << "    cursor: pointer;\n";
    out << "    font-size: 1em;\n";
    out << "    margin: 5px 0;\n";
    out << "    transition: background-color 0.3s ease;\n";
    out << "  }\n";
    out << "  button:hover {\n";
    out << "    background-color: #45a049;\n";
    out << "  }\n";
    out << "  button.close-btn {\n";
    out << "    background-color: #f44336;\n";
    out << "  }\n";
    out << "  button.close-btn:hover {\n";
    out << "    background-color: #d32f2f;\n";
    out << "  }\n";
    out << "  table {\n";
    out << "    width: 100%;\n";
    out << "    border-collapse: collapse;\n";
    out << "    margin-bottom: 20px;\n";
    out << "    background-color: white;\n";
    out << "    box-shadow: 0 0 10px rgba(0, 0, 0, 0.1);\n";
    out << "  }\n";
    out << "  th, td {\n";
    out << "    padding: 12px 15px;\n";
    out << "    border: 1px solid #dddddd;\n";
    out << "    text-align: left;\n";
    out << "    font-size: 1em;\n";
    out << "  }\n";
    out << "  th {\n";
    out << "    background-color: #4CAF50;\n";
    out << "    color: white;\n";
    out << "    text-transform: uppercase;\n";
    out << "    letter-spacing: 0.1em;\n";
    out << "  }\n";
    out << "  tr:nth-child(even) {\n";
    out << "    background-color: #f9f9f9;\n";
    out << "  }\n";
    out << "  td {\n";
    out << "    font-size: 0.95em;\n";
    out << "  }\n";
    out << "  td:last-child {\n";
    out << "    text-align: right;\n";
    out << "  }\n";
    out << "  ul {\n";
    out << "    padding-left: 20px;\n";
    out << "  }\n";
    out << "  ul li {\n";
    out << "    font-size: 1em;\n";
    out << "    margin-bottom: 10px;\n";
    out << "  }\n";
    out << "  ul li strong {\n";
    out << "    display: block;\n";
    out << "    margin-bottom: 5px;\n";
    out << "  }\n";
    out << "  @media (max-width: 768px) {\n";
    out << "    h2 {\n";
    out << "      font-size: 1.5em;\n";
    out << "    }\n";
    out << "    table, th, td {\n";
    out << "      font-size: 0.85em;\n";
    out << "    }\n";
    out << "    button {\n";
    out << "      font-size: 0.9em;\n";
    out << "    }\n";
    out << "  }\n";
    out << "</style>\n";
    out << "</head>\n";
    out << "<body>\n";
    out << "<h2>MVTX Run Summary</h2>" << endl;
    //out << "<button style=\"background-color: #4CAF50; color: white; padding: 10px 20px; border: none; cursor: pointer;\" onclick=\"window.print()\">Download as PDF</button>" << endl; 
    //out << "<button style=\"background-color: #4CAF50; color: white; padding: 10px 20px; border: none; cursor: pointer;\" onclick=\"window.open('https://sphenix-intra.sdcc.bnl.gov/WWW/user/hjheng/Run2024/MvtxQA-Summary/SummaryPlots.pdf')\">Download Summary Plots</button>" << endl;
    //out << "<button style=\"background-color: #4CAF50; color: white; padding: 10px 20px; border: none; cursor: pointer;\" onclick=\"window.open('https://sphenix-intra.sdcc.bnl.gov/WWW/user/hjheng/Run2024/MvtxQA-Summary/AddedRuns_ClusPhi.pdf')\">Cluster Phi (Since Last Update)</button>" << endl;
    //out << "<button style=\"background-color: #4CAF50; color: white; padding: 10px 20px; border: none; cursor: pointer;\" onclick=\"window.open('https://sphenix-intra.sdcc.bnl.gov/WWW/user/hjheng/Run2024/MvtxQA-Summary/AddedRuns_Occupancy.pdf')\">Chip Occupancy (Since Last Update)</button>" << endl;
    out << "<br>" << endl;
    out << "<br>" << endl;
    out << "<button style=\"background-color: #f44336; color: white; padding: 10px 20px; border: none; cursor: pointer;\" onclick=\"window.close()\">Close</button>" << endl;
    out << "<br>" << endl;
    out << "<br>" << endl;
    out << "<table>" << endl;
    out << "  <tr>" << endl;
    out << "    <th>Run Status</th>" << endl;
    out << "    <th>Number of Runs</th>" << endl;
    out << "    <th>Fraction</th>" << endl;
    out << "  </tr>" << endl;
    out << "  <tr>" << endl;
    out << "    <td>Total number of runs</td>" << endl;
    out << "    <td>" << map_runstatuslist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}"].size() + map_runstatuslist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}"].size() +
               map_runstatuslist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}"].size() + map_runstatuslist["\\textcolor{Purple}{No Cluster, Hit QA}"].size() +
                map_runstatuslist["Raw hit acceptance: \\textcolor{Red}{BAD}"].size()
        << "</td>" << endl;
    out << "    <td>-</td>" << endl;
    out << "  </tr>" << endl;
    for (const auto &elem : map_runstatuslist)
    {
        if (elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}" && 
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}" && 
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}" &&
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}; Run length: \\textcolor{Orange}{Short}" &&
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}; Run length: \\textcolor{ForestGreen}{Normal}" &&
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}; Run length: \\textcolor{Orange}{Short}" &&
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}; Run length: \\textcolor{ForestGreen}{Normal}" &&
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}; Run length: \\textcolor{Orange}{Short}" &&
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}; Run length: \\textcolor{ForestGreen}{Normal}")
        {
            float frac =
                (float)elem.second.size() / (map_runstatuslist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}"].size() + 
                                             map_runstatuslist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}"].size() +
                                             map_runstatuslist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}"].size() + 
                                             map_runstatuslist["\\textcolor{Purple}{No Cluster, Hit QA}"].size() +
                                             map_runstatuslist["Raw hit acceptance: \\textcolor{Red}{BAD}"].size());
            // modify the string for html
            std::string ms = modifystr4html(elem.first);
            // Format the fraction to percentage and 3 significant figures
            stringstream ss;
            ss << setprecision(3) << frac * 100.;
            string fracstr = ss.str();
            out << "  <tr>" << endl;
            out << "    <td>" << ms << "</td>" << endl;
            out << "    <td>" << elem.second.size() << "</td>" << endl;
            out << "    <td>" << fracstr << "\%</td>"
                << "</td>" << endl;
            out << "  </tr>" << endl;
        }
    }
    out << "</table>" << endl;
    out << "<ul>" << endl;
    for (const auto &elem : map_runstatuslist)
    {
        if (elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}" && 
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}" && 
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}" &&
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}; Run length: \\textcolor{Orange}{Short}" &&
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}; Run length: \\textcolor{ForestGreen}{Normal}" &&
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}; Run length: \\textcolor{Orange}{Short}" &&
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}; Run length: \\textcolor{ForestGreen}{Normal}" &&
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}; Run length: \\textcolor{Orange}{Short}" &&
            elem.first != "Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}; Run length: \\textcolor{ForestGreen}{Normal}")
        {
            // modify the string for html
            std::string ms = modifystr4html(elem.first);
            out << "<li><strong>" << ms << "</strong><br>" << endl;
            out << "Run list: ";
            for (const auto &run : elem.second)
            {
                if (run != elem.second.back())
                    out << run << ", ";
                else
                    out << run << "<br>";
            }
            out << "</li>" << endl;
        }
    }
    out << "</ul>" << endl;
    // out << "<button style=\"background-color: #4CAF50; color: white; padding: 10px 20px; border: none; cursor: pointer;\" onclick=\"window.print()\">Download as PDF</button>" << endl; 
    // out << "<button style=\"background-color: #4CAF50; color: white; padding: 10px 20px; border: none; cursor: pointer;\" onclick=\"window.open('https://sphenix-intra.sdcc.bnl.gov/WWW/user/hjheng/Run2024/MvtxQA-Summary/SummaryPlots.pdf')\">Download Summary Plots</button>" << endl;
    // out << "<button style=\"background-color: #f44336; color: white; padding: 10px 20px; border: none; cursor: pointer;\" onclick=\"window.close()\">Close</button>" << endl;
    out << "</body>" << endl;
    out << "</html>" << endl;

    out.close();
}

std::vector<int> getrunlist_processed()
{
    std::vector<int> v_runlist;
    std::string inputfile = "run_processed.txt";
    std::ifstream in(inputfile);
    if (!in)
    {
        std::cerr << "Cannot open the input file " << inputfile << std::endl;
        exit(1);
    }
    std::string line;
    while (std::getline(in, line))
    {
        std::vector<std::string> v;
        boost::split(v, line, boost::is_any_of("\t"));
        v_runlist.push_back(::atoi(v.at(0).c_str()));
    }
    in.close();
    return v_runlist;
}

void addrun2processed(std::vector<int> v_newruns)
{
    std::string inputfile = "run_processed.txt";
    std::ofstream out(inputfile, std::ios::app);
    if (!out)
    {
        std::cerr << "Cannot open the output file " << inputfile << std::endl;
        exit(1);
    }
    for (const auto &run : v_newruns)
    {
        out << run << std::endl;
    }
    out.close();
}

void runsummary()
{
    system("mkdir -p ./ClusPhiChecker");
    system("mkdir -p ./ChipOccupancyChecker");
    std::string inputbasdir = "/sphenix/data/data02/sphnxpro/QAhtml/aggregated";
    std::map<std::string, std::tuple<float, float, float, float>> map_run_fitresult;
    std::map<std::string, float> map_run_chipoccupancy_overflow;

    // *********** NEW MAP FOR RAW HIT QA ***********//
    std::map<std::string, std::tuple<float, float, float>> map_run_deadstaveschips;
    // *********************** //

    std::vector<int> v_processedrun = getrunlist_processed();
    std::vector<int> v_addedrun;

    // Adjust the thresholds here
    float scaledchi2ndf_good = 0.0015;
    float scaledchi2ndf_bad = 0.01;
    float scaledchi2ndf_good_b49961 = 0.0015;
    float scaledchi2ndf_bad_b49961 = 0.01;
    float scaledchi2ndf_good_a49961 = 0.015;
    float scaledchi2ndf_bad_a49961 = 0.025;
    float avgclus_threshold = -1E9; // will be set later
    float bovera_high = 0.28;
    float bovera_low = 0.24;
    float bovera_high_b49961 = 0.28;
    float bovera_low_b49961 = 0.24;
    float bovera_high_a49961 = 0.1;
    float bovera_low_a49961 = 0.05;

    // *********** DEAD AREA THRESHOLDS FOR RAW HIT QA ***********//
    float layer0hitthresh = 0.05;
    float layer1hitthresh = 0.05;
    float layer2hitthresh = 0.05;
    // *********************** //

    for (const auto &entry : std::filesystem::directory_iterator(inputbasdir))
    {
        std::string inputfile = entry.path().string();
        
        if (inputfile.find("HIST_DST_TRKR_CLUSTER") == std::string::npos)
        {
            continue;
        }
              

        std::string runnumber = inputfile.substr(inputfile.find("-000") + 4, 5);

        std::cout << "Processing " << inputfile << " " << runnumber << std::endl;
        TFile *f = new TFile(inputfile.c_str(), "READ");
        f->cd();
        TH1F *h_MvtxClusterQA_clusterPhi_incl = dynamic_cast<TH1F *>(f->Get("h_MvtxClusterQA_clusterPhi_incl"));
        TH1F *h_MvtxClusterQA_chipOccupancy = dynamic_cast<TH1F *>(f->Get("h_MvtxClusterQA_chipOccupancy"));

        if (!h_MvtxClusterQA_clusterPhi_incl && !h_MvtxClusterQA_chipOccupancy)
        {
            std::cerr << "One of the histograms h_MvtxClusterQA_clusterPhi_incl or h_MvtxClusterQA_chipOccupancy is missing" << std::endl;
            continue;
        }

        // *********** 2D HISTOGRAMS FOR RAW HIT QA ***********//
        auto pos = inputfile.find("HIST_DST_TRKR_CLUSTER");
        std::string inputfile_hit = inputfile;
        inputfile_hit.replace(pos, std::string("HIST_DST_TRKR_CLUSTER").size(), "HIST_DST_TRKR_HIT");

        if (!std::filesystem::exists(inputfile_hit))
        {
            continue;
        }

        TFile *f_hit = new TFile(inputfile_hit.c_str(), "READ");
        f_hit->cd();
        TH2F *h_MvtxRawHitQA_nhits_stave_chip_layer0 = dynamic_cast<TH2F *>(f_hit->Get("h_MvtxRawHitQA_nhits_stave_chip_layer0"));
        TH2F *h_MvtxRawHitQA_nhits_stave_chip_layer1 = dynamic_cast<TH2F *>(f_hit->Get("h_MvtxRawHitQA_nhits_stave_chip_layer1"));
        TH2F *h_MvtxRawHitQA_nhits_stave_chip_layer2 = dynamic_cast<TH2F *>(f_hit->Get("h_MvtxRawHitQA_nhits_stave_chip_layer2"));
        
        if (!h_MvtxRawHitQA_nhits_stave_chip_layer0 || !h_MvtxRawHitQA_nhits_stave_chip_layer1 || !h_MvtxRawHitQA_nhits_stave_chip_layer2)
        {
            std::cerr << "One of the histograms h_MvtxRawHitQA_nhits_stave_chip_layer\\%d is missing" << std::endl;
            continue;
        }
        // *********************** //

               

        if (std::find(v_processedrun.begin(), v_processedrun.end(), ::atoi(runnumber.c_str())) == v_processedrun.end())
            v_addedrun.push_back(::atoi(runnumber.c_str()));

        std::tuple<float, float, float, float> fitresult;
        if (h_MvtxClusterQA_clusterPhi_incl)
            fitresult = fitClusPhi(h_MvtxClusterQA_clusterPhi_incl);
        if (h_MvtxClusterQA_clusterPhi_incl)
            map_run_fitresult[runnumber] = fitresult;


        // *********** DO RAW HIT QA (GET FRACTION OF DEAD STAVES/CHIPS) ***********//
        std::tuple<float, float, float> deadstaveschips;
        if (h_MvtxRawHitQA_nhits_stave_chip_layer0)
            std::get<0>(deadstaveschips) = rawHitAcceptance(h_MvtxRawHitQA_nhits_stave_chip_layer0);
        if (h_MvtxRawHitQA_nhits_stave_chip_layer1)
            std::get<1>(deadstaveschips) = rawHitAcceptance(h_MvtxRawHitQA_nhits_stave_chip_layer1);
        if (h_MvtxRawHitQA_nhits_stave_chip_layer2)
            std::get<2>(deadstaveschips) = rawHitAcceptance(h_MvtxRawHitQA_nhits_stave_chip_layer2);

        if (h_MvtxRawHitQA_nhits_stave_chip_layer0 && h_MvtxRawHitQA_nhits_stave_chip_layer1 && h_MvtxRawHitQA_nhits_stave_chip_layer2)
            map_run_deadstaveschips[runnumber] = deadstaveschips;
         // *********************** //
        
        


        TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
        gPad->SetRightMargin(0.08);
        gPad->SetTopMargin(0.08);
        if (h_MvtxClusterQA_clusterPhi_incl)
            h_MvtxClusterQA_clusterPhi_incl->Draw();
        TLegend *leg = new TLegend(0.2, 0.2, 0.88, 0.45);
        leg->SetTextSize(0.03);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetHeader(Form("Run %s Fit Result: A+B#timescos(#phi+C)", runnumber.c_str()));
        leg->SetNColumns(2);
        leg->AddEntry((TObject *)0, Form("#chi^{2}/ndf=%g", std::get<3>(fitresult) * std::get<0>(fitresult)), "");
        leg->AddEntry((TObject *)0, Form("A=%g", std::get<0>(fitresult)), "");
        leg->AddEntry((TObject *)0, Form("Scaled #chi^{2}/ndf=%g", std::get<3>(fitresult)), "");
        leg->AddEntry((TObject *)0, Form("B=%g", std::get<1>(fitresult)), "");
        leg->AddEntry((TObject *)0, Form("log_{10}(scaled #chi^{2}/ndf)=%g", log10(std::get<3>(fitresult))), "");
        leg->AddEntry((TObject *)0, Form("B/A=%g", std::get<1>(fitresult) / std::get<0>(fitresult)), "");
        leg->AddEntry((TObject *)0, "", "");
        leg->AddEntry((TObject *)0, Form("C=%g", std::get<2>(fitresult)), "");
        leg->Draw();
        c1->SaveAs(Form("./ClusPhiChecker/ClusPhiChecker_%s.pdf", runnumber.c_str()));
        c1->SaveAs(Form("./ClusPhiChecker/ClusPhiChecker_%s.png", runnumber.c_str()));
        //c1->SaveAs(Form("/sphenix/WWW/user/hjheng/Run2024/MvtxQA-Summary/ClusPhi/ClusPhiChecker_%s.pdf", runnumber.c_str()));
        //c1->SaveAs(Form("/sphenix/WWW/user/hjheng/Run2024/MvtxQA-Summary/ClusPhi/ClusPhiChecker_%s.png", runnumber.c_str()));
        c1->cd();
        if (h_MvtxClusterQA_chipOccupancy->Integral(-1, -1) > 0.)
            c1->SetLogy(1);
        h_MvtxClusterQA_chipOccupancy->Scale(1. / h_MvtxClusterQA_chipOccupancy->Integral(-1, -1));
        h_MvtxClusterQA_chipOccupancy->GetXaxis()->SetRange(1, h_MvtxClusterQA_chipOccupancy->GetNbinsX() + 1);
        h_MvtxClusterQA_chipOccupancy->GetXaxis()->SetTitle("Chip occupancy [%]");
        h_MvtxClusterQA_chipOccupancy->GetYaxis()->SetTitle("Normalized Entries");
        // fix the y axis range
        h_MvtxClusterQA_chipOccupancy->GetYaxis()->SetRangeUser(1E-8, 3);
        h_MvtxClusterQA_chipOccupancy->Draw();
        TLegend *leg_runnum = new TLegend(0.7, 0.8, 0.88, 0.85);
        leg_runnum->SetTextSize(0.04);
        leg_runnum->SetBorderSize(0);
        leg_runnum->SetFillStyle(0);
        leg_runnum->SetHeader(Form("Run %s", runnumber.c_str()));
        leg_runnum->Draw();
        c1->SaveAs(Form("./ChipOccupancyChecker/ChipOccupancy_%s.pdf", runnumber.c_str()));
        c1->SaveAs(Form("./ChipOccupancyChecker/ChipOccupancy_%s.png", runnumber.c_str()));
        //c1->SaveAs(Form("/sphenix/WWW/user/hjheng/Run2024/MvtxQA-Summary/ChipOccupancy/ChipOccupancy_%s.pdf", runnumber.c_str()));
        //c1->SaveAs(Form("/sphenix/WWW/user/hjheng/Run2024/MvtxQA-Summary/ChipOccupancy/ChipOccupancy_%s.png", runnumber.c_str()));

        // sum the entries for bins > 0.3%
        map_run_chipoccupancy_overflow[runnumber] = h_MvtxClusterQA_chipOccupancy->Integral(h_MvtxClusterQA_chipOccupancy->FindBin(0.3), h_MvtxClusterQA_chipOccupancy->GetNbinsX() + 1);
        cout << "Sum of entries larger than 0.3%: " << map_run_chipoccupancy_overflow[runnumber] << endl;

        delete f;
        delete f_hit;
        delete c1;
    }

    // check if the two maps have the same size
    // if (map_run_fitresult.size() != map_run_chipoccupancy_overflow.size())
    // {
    //     std::cerr << "The two maps have different sizes, exiting..." << std::endl;
    //     exit(1);
    // }

    

    TH1F *hM_chi2ndf = new TH1F("hM_chi2ndf", "Chi2/ndf of the cluster phi fit", 150, 0, 0.015);
    TH1F *hM_BoverA = new TH1F("hM_BoverA", "B/A ratio of the cluster phi fit", 120, 0.15, 0.45);
    TH1F *hM_logoccupancy = new TH1F("hM_logoccupancy", "Log10 of the sum of chip occupancy", 120, -6, -3);
    TGraph *gM_chi2ndf = new TGraph();
    TGraph *gM_avgnclus = new TGraph();
    TGraph *gM_BoverA = new TGraph();
    TGraph *gM_chipoccupancyoverflow = new TGraph();
    TH2F *hM_logchi2ndf_BoverA = new TH2F("hM_logchi2ndf_BoverA", "Log10 of chi2/ndf vs B/A", 100, -3.5, 1.5, 100, 0.15, 0.35);
    std::vector<float> v_avgclus;

    for (const auto &elem : map_run_fitresult)
    {
        // int index = std::distance(map_run_fitresult.begin(), map_run_fitresult.find(elem.first)); // NOT USED; KEEP FOR RECORD
        hM_chi2ndf->Fill(std::get<3>(elem.second));
        hM_BoverA->Fill(std::get<1>(elem.second) / std::get<0>(elem.second));
        gM_chi2ndf->AddPoint(::atof(elem.first.c_str()), log10(std::get<3>(elem.second)));
        gM_avgnclus->AddPoint(::atof(elem.first.c_str()), std::get<0>(elem.second));
        gM_BoverA->AddPoint(::atof(elem.first.c_str()), std::get<1>(elem.second) / std::get<0>(elem.second));
        v_avgclus.push_back(std::get<0>(elem.second));
        hM_logchi2ndf_BoverA->Fill(log10(std::get<3>(elem.second)), std::get<1>(elem.second) / std::get<0>(elem.second));
    }
    
     

    for (const auto &elem : map_run_chipoccupancy_overflow)
    {
        // int index = std::distance(map_run_chipoccupancy_overflow.begin(), map_run_chipoccupancy_overflow.find(elem.first)); // NOT USED; KEEP FOR RECORD
        if (elem.second == 0 || std::isnan(elem.second))
            continue;

        hM_logoccupancy->Fill(log10(elem.second));
        gM_chipoccupancyoverflow->AddPoint(::atof(elem.first.c_str()), log10(elem.second));
    }

    

    // get the threshold for which 10% quantile of the average number of clusters
    std::sort(v_avgclus.begin(), v_avgclus.end());
    // remove 0 entries
    v_avgclus.erase(std::remove(v_avgclus.begin(), v_avgclus.end(), 0), v_avgclus.end());
    
    if (!v_avgclus.empty())
    {
        avgclus_threshold = v_avgclus.at(v_avgclus.size() * 0.1);
    }
    else
    {
        avgclus_threshold = 0;
    }

    cout << "10% of the runs have average number of clusters less than " << avgclus_threshold << endl;

    

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    c2->SetTopMargin(0.2);
    c2->SetRightMargin(0.05);
    c2->cd();
    hM_chi2ndf->GetXaxis()->SetTitle("Scaled #chi^{2}/ndf");
    hM_chi2ndf->GetYaxis()->SetTitle("Counts");
    float maxwooverflow = hM_chi2ndf->GetMaximum();
    // Draw overflow bin
    hM_chi2ndf->GetXaxis()->SetRange(1, hM_chi2ndf->GetNbinsX() + 1);
    hM_chi2ndf->Draw("hist");
    gPad->Update();
    TH1F *dummy_good = new TH1F("dummy_good", "dummy_good", 1, 0, scaledchi2ndf_good);
    dummy_good->SetFillColorAlpha(kGreen + 3, 0.2);
    dummy_good->SetLineWidth(0);
    // Maximum range of the canvas
    dummy_good->SetBinContent(1, gPad->GetFrame()->GetY2());
    dummy_good->Draw("same");
    TH1F *dummy_ok = new TH1F("dummy_bad", "dummy_bad", 1, scaledchi2ndf_good, scaledchi2ndf_bad);
    dummy_ok->SetFillColorAlpha(kOrange - 3, 0.2);
    dummy_ok->SetLineWidth(0);
    dummy_ok->SetBinContent(1, gPad->GetFrame()->GetY2());
    dummy_ok->Draw("same");
    TH1F *dummy_bad = new TH1F("dummy_bad", "dummy_bad", 1, scaledchi2ndf_bad, hM_chi2ndf->GetXaxis()->GetXmax() + hM_chi2ndf->GetBinWidth(1));
    dummy_bad->SetFillColorAlpha(kRed, 0.2);
    dummy_bad->SetLineWidth(0);
    dummy_bad->SetBinContent(1, gPad->GetFrame()->GetY2());
    dummy_bad->Draw("same");
    TLine *l1 = new TLine(scaledchi2ndf_good, 0, scaledchi2ndf_good, gPad->GetFrame()->GetY2());
    l1->SetLineColor(kOrange);
    l1->SetLineWidth(2);
    l1->Draw();
    TLine *l2 = new TLine(scaledchi2ndf_bad, 0, scaledchi2ndf_bad, gPad->GetFrame()->GetY2());
    l2->SetLineColor(kRed);
    l2->SetLineWidth(2);
    l2->Draw();
    hM_chi2ndf->Draw("histsame");
    TLegend *l = new TLegend(0.2, 0.85, 0.85, 0.95);
    l->SetTextSize(0.04);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->AddEntry(l1, Form("Threshold for GOOD (scaled #chi^{2}/ndf #leq %g)", scaledchi2ndf_good), "l");
    l->AddEntry(l2, Form("Threshold for BAD (scaled #chi^{2}/ndf > %g)", scaledchi2ndf_bad), "l");
    l->Draw();
    c2->SaveAs("./Chi2NdfHist.pdf");
    c2->SaveAs("./Chi2NdfHist.png");

    c2->cd();
    c2->SetTopMargin(0.05);
    hM_BoverA->GetXaxis()->SetTitle("B/A");
    hM_BoverA->GetYaxis()->SetTitle("Counts");
    hM_BoverA->Draw("hist");
    // Get 5% and 95% quantiles and draw TLines
    Double_t BoverA5, BoverA95, q10 = 0.1, q90 = 0.9;
    hM_BoverA->GetQuantiles(1, &BoverA5, &q10);
    hM_BoverA->GetQuantiles(1, &BoverA95, &q90);
    TLine *l11 = new TLine(BoverA5, hM_BoverA->GetMinimum(), BoverA5, hM_BoverA->GetMaximum());
    l11->SetLineColorAlpha(kBlue, 0.6);
    l11->SetLineStyle(2);
    l11->SetLineWidth(2);
    l11->Draw();
    TLine *l12 = new TLine(BoverA95, hM_BoverA->GetMinimum(), BoverA95, hM_BoverA->GetMaximum());
    l12->SetLineColorAlpha(kBlue, 0.6);
    l12->SetLineStyle(2);
    l12->SetLineWidth(2);
    l12->Draw();
    c2->SaveAs("./BoverAHist.pdf");
    c2->SaveAs("./BoverAHist.png");

    c2->cd();
    c2->SetTopMargin(0.05);
    // c2->SetLogy(1);
    hM_logoccupancy->GetXaxis()->SetTitle("Log_{10}(Prob. of chip occu.#geq 0.3%)");
    hM_logoccupancy->GetYaxis()->SetTitle("Counts");
    hM_logoccupancy->Draw("hist");
    c2->SaveAs("./ChipOccupancyHist.pdf");
    c2->SaveAs("./ChipOccupancyHist.png");

    c2->SetTopMargin(0.05);
    TH1F *hM_logoccupancy_clone = (TH1F *)hM_logoccupancy->Clone();
    hM_logoccupancy_clone->Scale(1. / hM_logoccupancy_clone->Integral());
    TH1F *hM_logoccupancy_cumulative = (TH1F *)hM_logoccupancy_clone->GetCumulative();
    hM_logoccupancy_cumulative->GetXaxis()->SetTitle("Log_{10}(Prob. of chip occu.#geq 0.3%)");
    hM_logoccupancy_cumulative->GetYaxis()->SetTitle("Cumulative fraction");
    hM_logoccupancy_cumulative->Draw("hist");

    // get the X value when the cumulative fraction is 0.95
    float x95 = 0;
    for (int i = 1; i < hM_logoccupancy_cumulative->GetNbinsX(); i++)
    {
        if (hM_logoccupancy_cumulative->GetBinContent(i) > 0.95)
        {
            x95 = hM_logoccupancy_cumulative->GetBinCenter(i);
            break;
        }
    }
    cout << "95% of the runs have chip occupancy overflow less than " << pow(10, x95) << endl;
    TLine *l6 = new TLine(x95, hM_logoccupancy_cumulative->GetMinimum(), x95, hM_logoccupancy_cumulative->GetMaximum());
    l6->SetLineColorAlpha(kOrange + 10, 0.6);
    l6->SetLineStyle(2);
    l6->SetLineWidth(2);
    l6->Draw();
    c2->SaveAs("./ChipOccupancyCumulativeHist.pdf");
    c2->SaveAs("./ChipOccupancyCumulativeHist.png");

    c2->cd();
    c2->SetTopMargin(0.05);
    c2->SetRightMargin(0.13);
    hM_logchi2ndf_BoverA->GetXaxis()->SetTitle("Log_{10}(scaled #chi^{2}/ndf)");
    hM_logchi2ndf_BoverA->GetYaxis()->SetTitle("B/A");
    hM_logchi2ndf_BoverA->Draw("colz");
    gPad->Update();
    // Add TLines for the cut values
    TLine *l8 = new TLine(log10(scaledchi2ndf_good), gPad->GetFrame()->GetY1(), log10(scaledchi2ndf_good), gPad->GetFrame()->GetY2());
    l8->SetLineColorAlpha(kOrange, 0.6);
    l8->SetLineStyle(2);
    l8->SetLineWidth(2);
    l8->Draw();
    TLine *l9 = new TLine(log10(scaledchi2ndf_bad), gPad->GetFrame()->GetY1(), log10(scaledchi2ndf_bad), gPad->GetFrame()->GetY2());
    l9->SetLineColorAlpha(kRed, 0.6);
    l9->SetLineStyle(2);
    l9->SetLineWidth(2);
    l9->Draw();
    TLine *l10 = new TLine(gPad->GetFrame()->GetX1(), 0.24, gPad->GetFrame()->GetX2(), 0.24);
    l10->SetLineColorAlpha(kRed, 0.6);
    l10->SetLineStyle(2);
    l10->SetLineWidth(2);
    l10->Draw();
    c2->SaveAs("./LogChi2NdfBoverA.pdf");
    c2->SaveAs("./LogChi2NdfBoverA.png");
    delete c2;

    std::map<std::string, std::vector<int>> map_runstatus_runlist;

    for (const auto &elem : map_run_fitresult)
    {

        auto elem_rawhit = map_run_deadstaveschips.find(elem.first);

        if (std::isnan(std::get<3>(elem.second)) || elem_rawhit == map_run_deadstaveschips.end())
        {
            map_runstatus_runlist["\\textcolor{Purple}{No Cluster, Hit QA}"].push_back(::atoi(elem.first.c_str()));
        }

        else if (std::get<0>(elem_rawhit->second) < layer0hitthresh && std::get<1>(elem_rawhit->second) < layer1hitthresh && std::get<2>(elem_rawhit->second) < layer2hitthresh)
        {
                // else if (std::get<3>(elem.second) <= scaledchi2ndf_good)
            // run-dependent thresholds
            if ((::atoi(elem.first.c_str()) <= 49961 && std::get<3>(elem.second) <= scaledchi2ndf_good_b49961) ||
                    (::atoi(elem.first.c_str()) > 49961 && std::get<3>(elem.second) <= scaledchi2ndf_good_a49961))
            {
                map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}"].push_back(::atoi(elem.first.c_str()));
                if (std::get<0>(elem.second) < avgclus_threshold)
                {
                    map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}; Run length: \\textcolor{Orange}{Short}"].push_back(::atoi(elem.first.c_str()));
                    // if ((std::get<1>(elem.second) / std::get<0>(elem.second)) < bovera_low || (std::get<1>(elem.second) / std::get<0>(elem.second)) > bovera_high)
                    if ((::atoi(elem.first.c_str()) <= 49961 && (std::get<1>(elem.second) / std::get<0>(elem.second)) < bovera_low_b49961) ||
                        (::atoi(elem.first.c_str()) > 49961 && (std::get<1>(elem.second) / std::get<0>(elem.second)) < bovera_low_a49961))
                    {
                        map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}; Run length: \\textcolor{Orange}{Short}; B/A ratio: \\textcolor{Red}{BAD}"].push_back(
                            ::atoi(elem.first.c_str()));
                    }
                    else
                        map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}; Run length: \\textcolor{Orange}{Short}; B/A ratio: \\textcolor{ForestGreen}{GOOD}"].push_back(
                            ::atoi(elem.first.c_str()));
                }
                else
                {
                    map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}; Run length: \\textcolor{ForestGreen}{Normal}"].push_back(::atoi(elem.first.c_str()));
                    // if ((std::get<1>(elem.second) / std::get<0>(elem.second)) < bovera_low || (std::get<1>(elem.second) / std::get<0>(elem.second)) > bovera_high)
                    if ((::atoi(elem.first.c_str()) <= 49961 && (std::get<1>(elem.second) / std::get<0>(elem.second)) < bovera_low_b49961) ||
                        (::atoi(elem.first.c_str()) > 49961 && (std::get<1>(elem.second) / std::get<0>(elem.second)) < bovera_low_a49961))
                    {
                        map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}; Run length: \\textcolor{ForestGreen}{Normal}; B/A ratio: \\textcolor{Red}{BAD}"].push_back(
                            ::atoi(elem.first.c_str()));
                    }
                    else
                        map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}; Run length: \\textcolor{ForestGreen}{Normal}; B/A ratio: \\textcolor{ForestGreen}{GOOD}"].push_back(
                            ::atoi(elem.first.c_str()));
                }
            }
            // else if (std::get<3>(elem.second) > scaledchi2ndf_good && std::get<3>(elem.second) <= scaledchi2ndf_bad)
            else if ((::atoi(elem.first.c_str()) <= 49961 && std::get<3>(elem.second) > scaledchi2ndf_good_b49961 && std::get<3>(elem.second) <= scaledchi2ndf_bad_b49961) ||
                    (::atoi(elem.first.c_str()) > 49961 && std::get<3>(elem.second) > scaledchi2ndf_good_a49961 && std::get<3>(elem.second) <= scaledchi2ndf_bad_a49961))
            {
                map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}"].push_back(::atoi(elem.first.c_str()));
                if (std::get<0>(elem.second) < avgclus_threshold)
                {
                    map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}; Run length: \\textcolor{Orange}{Short}"].push_back(::atoi(elem.first.c_str()));
                    // if ((std::get<1>(elem.second) / std::get<0>(elem.second)) < bovera_low || (std::get<1>(elem.second) / std::get<0>(elem.second)) > bovera_high)
                    if ((::atoi(elem.first.c_str()) <= 49961 && (std::get<1>(elem.second) / std::get<0>(elem.second)) < bovera_low_b49961) ||
                        (::atoi(elem.first.c_str()) > 49961 && (std::get<1>(elem.second) / std::get<0>(elem.second)) < bovera_low_a49961))
                    {
                        map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}; Run length: \\textcolor{Orange}{Short}; B/A ratio: \\textcolor{Red}{BAD}"].push_back(::atoi(elem.first.c_str()));
                    }
                    else
                        map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}; Run length: \\textcolor{Orange}{Short}; B/A ratio: \\textcolor{ForestGreen}{GOOD}"].push_back(
                            ::atoi(elem.first.c_str()));
                }
                else
                {
                    map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}; Run length: \\textcolor{ForestGreen}{Normal}"].push_back(::atoi(elem.first.c_str()));
                    // if ((std::get<1>(elem.second) / std::get<0>(elem.second)) < bovera_low || (std::get<1>(elem.second) / std::get<0>(elem.second)) > bovera_high)
                    if ((::atoi(elem.first.c_str()) <= 49961 && (std::get<1>(elem.second) / std::get<0>(elem.second)) < bovera_low_b49961) ||
                        (::atoi(elem.first.c_str()) > 49961 && (std::get<1>(elem.second) / std::get<0>(elem.second)) < bovera_low_a49961))
                    {
                        map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}; Run length: \\textcolor{ForestGreen}{Normal}; B/A ratio: \\textcolor{Red}{BAD}"].push_back(
                            ::atoi(elem.first.c_str()));
                    }
                    else
                        map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}; Run length: \\textcolor{ForestGreen}{Normal}; B/A ratio: \\textcolor{ForestGreen}{GOOD}"].push_back(
                            ::atoi(elem.first.c_str()));
                }
            }
            else
            {
                map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}"].push_back(::atoi(elem.first.c_str()));
                if (std::get<0>(elem.second) < avgclus_threshold)
                {
                    map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}; Run length: \\textcolor{Orange}{Short}"].push_back(::atoi(elem.first.c_str()));
                    // if ((std::get<1>(elem.second) / std::get<0>(elem.second)) < bovera_low || (std::get<1>(elem.second) / std::get<0>(elem.second)) > bovera_high)
                    if ((::atoi(elem.first.c_str()) <= 49961 && (std::get<1>(elem.second) / std::get<0>(elem.second)) < bovera_low_b49961) ||
                        (::atoi(elem.first.c_str()) > 49961 && (std::get<1>(elem.second) / std::get<0>(elem.second)) < bovera_low_a49961))
                    {
                        map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}; Run length: \\textcolor{Orange}{Short}; B/A ratio: \\textcolor{Red}{BAD}"].push_back(::atoi(elem.first.c_str()));
                    }
                    else
                        map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}; Run length: \\textcolor{Orange}{Short}; B/A ratio: \\textcolor{ForestGreen}{GOOD}"].push_back(
                            ::atoi(elem.first.c_str()));
                }
                else
                {
                    map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}; Run length: \\textcolor{ForestGreen}{Normal}"].push_back(::atoi(elem.first.c_str()));
                    // if ((std::get<1>(elem.second) / std::get<0>(elem.second)) < bovera_low || (std::get<1>(elem.second) / std::get<0>(elem.second)) > bovera_high)
                    if ((::atoi(elem.first.c_str()) <= 49961 && (std::get<1>(elem.second) / std::get<0>(elem.second)) < bovera_low_b49961) ||
                        (::atoi(elem.first.c_str()) > 49961 && (std::get<1>(elem.second) / std::get<0>(elem.second)) < bovera_low_a49961))
                    {
                        map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}; Run length: \\textcolor{ForestGreen}{Normal}; B/A ratio: \\textcolor{Red}{BAD}"].push_back(
                            ::atoi(elem.first.c_str()));
                    }
                    else
                        map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}; Run length: \\textcolor{ForestGreen}{Normal}; B/A ratio: \\textcolor{ForestGreen}{GOOD}"].push_back(
                            ::atoi(elem.first.c_str()));
                }
            }
        }
        else
        {
            map_runstatus_runlist["Raw hit acceptance: \\textcolor{Red}{BAD}"].push_back(::atoi(elem.first.c_str()));
        }
        
    }

    
    makesummarytable(map_runstatus_runlist);
    makehtml(map_runstatus_runlist);

    TCanvas *c3 = new TCanvas("c3", "c3", 1500, 600);
    c3->SetRightMargin(0.02);
    c3->SetLeftMargin(0.1);
    c3->cd();
    gM_chi2ndf->SetMarkerStyle(kMarkerStyle);
    gM_chi2ndf->SetLineWidth(2);
    gM_chi2ndf->SetMarkerSize(kMarkerSize);
    gM_chi2ndf->SetMarkerColor(kBlue);
    gM_chi2ndf->GetXaxis()->SetTitle("Run number");
    gM_chi2ndf->GetYaxis()->SetTitle("Log_{10}(scaled #chi^{2}/ndf)");
    // gM_chi2ndf->GetHistogram()->SetMaximum(0.03);
    // gM_chi2ndf->GetHistogram()->SetMinimum(0);
    gM_chi2ndf->GetYaxis()->SetTitleOffset(0.8);
    gM_chi2ndf->Draw("AP");
    // Add TLines for the cut values
    // TLine *l3 = new TLine(gM_chi2ndf->GetXaxis()->GetXmin(), log10(scaledchi2ndf_good), gM_chi2ndf->GetXaxis()->GetXmax(), log10(scaledchi2ndf_good));
    TLine *l3 = new TLine(gM_chi2ndf->GetXaxis()->GetXmin(), log10(scaledchi2ndf_good_b49961), 49961, log10(scaledchi2ndf_good_b49961));
    l3->SetLineColorAlpha(kOrange - 3, 0.6);
    l3->SetLineStyle(2);
    l3->SetLineWidth(2);
    l3->Draw();
    // TLine *l4 = new TLine(gM_chi2ndf->GetXaxis()->GetXmin(), log10(scaledchi2ndf_bad), gM_chi2ndf->GetXaxis()->GetXmax(), log10(scaledchi2ndf_bad));
    TLine *l4 = new TLine(gM_chi2ndf->GetXaxis()->GetXmin(), log10(scaledchi2ndf_bad_b49961), 49961, log10(scaledchi2ndf_bad_b49961));
    l4->SetLineColorAlpha(kRed, 0.6);
    l4->SetLineStyle(2);
    l4->SetLineWidth(2);
    l4->Draw();
    TLine *l17 = new TLine(49961, log10(scaledchi2ndf_good_a49961), gM_chi2ndf->GetXaxis()->GetXmax(), log10(scaledchi2ndf_good_a49961));
    l17->SetLineColorAlpha(kOrange - 3, 0.6);
    l17->SetLineStyle(2);
    l17->SetLineWidth(2);
    l17->Draw();
    TLine *l18 = new TLine(49961, log10(scaledchi2ndf_bad_a49961), gM_chi2ndf->GetXaxis()->GetXmax(), log10(scaledchi2ndf_bad_a49961));
    l18->SetLineColorAlpha(kRed, 0.6);
    l18->SetLineStyle(2);
    l18->SetLineWidth(2);
    l18->Draw();
    // Add some stats
    TLatex *runstat = new TLatex();
    float runstat_x0 = 0.15;
    runstat->SetTextSize(0.05);
    runstat->DrawLatexNDC(runstat_x0, 0.87, Form("Number of runs processed: %lu", map_run_fitresult.size()));
    runstat->SetTextSize(0.04);
    runstat->SetTextColor(kGreen + 3);
    runstat->DrawLatexNDC(runstat_x0, 0.82,
                          Form("GOOD (scaled #chi^{2}/ndf #leq %g): %lu runs (%2.1f%%)", scaledchi2ndf_good, map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}"].size(),
                               map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{ForestGreen}{GOOD}"].size() / (float)map_run_fitresult.size() * 100.));
    runstat->SetTextColor(kOrange - 3);
    runstat->DrawLatexNDC(runstat_x0, 0.77,
                          Form("OK (scaled #chi^{2}/ndf #leq %g): %lu runs (%2.1f%%)", scaledchi2ndf_bad, map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}"].size(),
                               map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Orange}{OK}"].size() / (float)map_run_fitresult.size() * 100.));
    runstat->SetTextColor(kRed);
    runstat->DrawLatexNDC(runstat_x0, 0.72,
                          Form("BAD (scaled #chi^{2}/ndf > %g): %lu runs (%2.1f%%)", scaledchi2ndf_bad, map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}"].size(),
                               map_runstatus_runlist["Raw hit acceptance: \\textcolor{ForestGreen}{GOOD}; $\\chi^{2}$/ndf: \\textcolor{Red}{BAD}"].size() / (float)map_run_fitresult.size() * 100.));
    runstat->SetTextColor(kMagenta + 1);
    runstat->DrawLatexNDC(runstat_x0, 0.67,
                          Form("No MVTX cluster QA: %lu runs (%2.1f%%)", map_runstatus_runlist["\\textcolor{Purple}{No Cluster QA}"].size(),
                               map_runstatus_runlist["\\textcolor{Purple}{No Cluster QA}"].size() / (float)map_run_fitresult.size() * 100.));
    c3->SaveAs("./Chi2NdfGraph.pdf");
    c3->SaveAs("./Chi2NdfGraph.png");

    c3->cd();
    c3->SetLogy(1);
    gM_avgnclus->SetMarkerStyle(kMarkerStyle);
    gM_avgnclus->SetLineWidth(2);
    gM_avgnclus->SetMarkerSize(kMarkerSize);
    gM_avgnclus->SetMarkerColor(kBlue);
    gM_avgnclus->GetXaxis()->SetTitle("Run number");
    gM_avgnclus->GetYaxis()->SetTitle("Average number of clusters (per bin)");
    gM_avgnclus->GetYaxis()->SetTitleOffset(0.8);
    gM_avgnclus->GetHistogram()->SetMaximum(gM_avgnclus->GetHistogram()->GetMaximum() * 10);
    gM_avgnclus->GetHistogram()->SetMinimum(0.1);
    gM_avgnclus->Draw("AP");
    TLine *l7 = new TLine(gM_avgnclus->GetXaxis()->GetXmin(), avgclus_threshold, gM_avgnclus->GetXaxis()->GetXmax(), avgclus_threshold);
    l7->SetLineColorAlpha(kOrange + 10, 0.6);
    l7->SetLineStyle(2);
    l7->SetLineWidth(2);
    l7->Draw();
    // Add text for the Tline
    TText *t = new TText();
    t->SetTextSize(0.03);
    t->SetTextAlign(33);
    t->SetTextColor(kOrange + 10);
    t->DrawText(gM_avgnclus->GetXaxis()->GetXmax() - 100, avgclus_threshold, "10th percentile");
    c3->SaveAs("./AvgnclusGraph.pdf");
    c3->SaveAs("./AvgnclusGraph.png");

    c3->cd();
    c3->SetLogy(0);
    gM_BoverA->SetMarkerStyle(kMarkerStyle);
    gM_BoverA->SetLineWidth(2);
    gM_BoverA->SetMarkerSize(kMarkerSize);
    gM_BoverA->SetMarkerColor(kBlue);
    gM_BoverA->GetXaxis()->SetTitle("Run number");
    gM_BoverA->GetYaxis()->SetTitle("Modulation amplitude-to-average (B/A) ratio");
    gM_BoverA->GetYaxis()->SetTitleOffset(0.8);
    gM_BoverA->GetHistogram()->SetMaximum(0.5);
    gM_BoverA->GetHistogram()->SetMinimum(0.);
    gM_BoverA->Draw("AP");
    // TLine *l13 = new TLine(gM_BoverA->GetXaxis()->GetXmin(), bovera_low, gM_BoverA->GetXaxis()->GetXmax(), bovera_low);
    TLine *l13 = new TLine(gM_BoverA->GetXaxis()->GetXmin(), bovera_low_b49961, 49961, bovera_low);
    l13->SetLineColorAlpha(kOrange + 10, 0.6);
    l13->SetLineStyle(2);
    l13->SetLineWidth(2);
    l13->Draw();
    // TLine *l14 = new TLine(gM_BoverA->GetXaxis()->GetXmin(), bovera_high, gM_BoverA->GetXaxis()->GetXmax(), bovera_high);
    TLine *l14 = new TLine(gM_BoverA->GetXaxis()->GetXmin(), bovera_high_b49961, 49961, bovera_high);
    l14->SetLineColorAlpha(kOrange + 10, 0.6);
    l14->SetLineStyle(2);
    l14->SetLineWidth(2);
    l14->Draw();
    TLine *l15 = new TLine(49961, bovera_low_a49961, gM_BoverA->GetXaxis()->GetXmax(), bovera_low_a49961);
    l15->SetLineColorAlpha(kOrange + 10, 0.6);
    l15->SetLineStyle(2);
    l15->SetLineWidth(2);
    l15->Draw();
    TLine *l16 = new TLine(49961, bovera_high_a49961, gM_BoverA->GetXaxis()->GetXmax(), bovera_high_a49961);
    l16->SetLineColorAlpha(kOrange + 10, 0.6);
    l16->SetLineStyle(2);
    l16->SetLineWidth(2);
    l16->Draw();
    c3->SaveAs("./BoverAGraph.pdf");
    c3->SaveAs("./BoverAGraph.png");

    c3->cd();
    c3->SetLogy(0);
    gM_chipoccupancyoverflow->SetMarkerStyle(kMarkerStyle);
    gM_chipoccupancyoverflow->SetLineWidth(2);
    gM_chipoccupancyoverflow->SetMarkerSize(kMarkerSize);
    gM_chipoccupancyoverflow->SetMarkerColor(kBlue);
    gM_chipoccupancyoverflow->GetXaxis()->SetTitle("Run number");
    gM_chipoccupancyoverflow->GetYaxis()->SetTitle("Log_{10}(Prob. of chip occu.#geq 0.3%)");
    gM_chipoccupancyoverflow->GetYaxis()->SetTitleOffset(0.8);
    gM_chipoccupancyoverflow->Draw("AP");
    TLine *l5 = new TLine(gM_chipoccupancyoverflow->GetXaxis()->GetXmin(), x95, gM_chipoccupancyoverflow->GetXaxis()->GetXmax(), x95);
    l5->SetLineColorAlpha(kOrange + 10, 0.6);
    l5->SetLineStyle(2);
    l5->SetLineWidth(2);
    l5->Draw();
    TText *t2 = new TText();
    t2->SetTextSize(0.03);
    t2->SetTextAlign(11);
    t2->SetTextColor(kOrange + 10);
    t2->DrawText(gM_chipoccupancyoverflow->GetXaxis()->GetXmin() + 125, x95 + 0.02, "95th percentile");
    c3->SaveAs("./ChipOccupancyOverflowGraph.pdf");
    c3->SaveAs("./ChipOccupancyOverflowGraph.png");
    delete c3;

    
    std::sort(v_addedrun.begin(), v_addedrun.end());
    /*
    if (v_addedrun.size() > 0)
    {
        std::string pdfmergecmd_clusphi = "gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=/sphenix/WWW/user/hjheng/Run2024/MvtxQA-Summary/AddedRuns_ClusPhi.pdf -dBATCH ";
        std::string pdfmergecmd_occupancy = "gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=/sphenix/WWW/user/hjheng/Run2024/MvtxQA-Summary/AddedRuns_Occupancy.pdf -dBATCH ";
        for (const auto &elem : v_addedrun)
        {
            pdfmergecmd_clusphi += Form("/sphenix/WWW/user/hjheng/Run2024/MvtxQA-Summary/ClusPhi/ClusPhiChecker_%d.pdf ", elem);
            pdfmergecmd_occupancy += Form("/sphenix/WWW/user/hjheng/Run2024/MvtxQA-Summary/ChipOccupancy/ChipOccupancy_%d.pdf ", elem);
        }
        system(pdfmergecmd_clusphi.c_str());
        system(pdfmergecmd_occupancy.c_str());
    }
    */
    

    addrun2processed(v_addedrun);
    
}