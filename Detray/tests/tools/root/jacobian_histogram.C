/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// ROOT include(s).
#include <TCanvas.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TStyle.h>

#include <ROOT/RCsvDS.hxx>
#include <ROOT/RDataFrame.hxx>

// System include(s).
#include <array>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

const double x_pos1 = -13.4;
const double x_pos2 = 0.06;
const double y_pos1 = 2e5;
const double y_pos2 = 6.75e4;

// const double title_x = x_pos;
// const double title_y = 0.8197;

const int label_font = 132;
const double label_font_size = 0.055;
const double header_text_size = 0.055;
const double geom_text_size = 0.0434028;
const double titleX_font_size = 0.05;
const double titleY_font_size = 0.055;
const double x_title_offset = 1.75;
const double y_title_offset = 1.34;
const double x_label_offset = 0.015;
const double y_label_offset = 0.015;
const int fill_style = 1001;
const double color_alpha = 0.5;

const std::array<float, 2> cdim{700, 600};
const double maximum = 1e6;

const double pad_x0 = 0.00;
const double pad_x1 = 1.;
const double pad_y0 = 0.00;
const double pad_y1 = 1.;

const double bin_width = 0.2;
}  // namespace

void draw_text(double x1, double y1, double y2, double s1, double s2,
               std::string t1, std::string t2) {
    TLatex* ttext1 = new TLatex(x1, y1, t1.c_str());
    TLatex* ttext2 = new TLatex(x1, y2, t2.c_str());
    ttext1->SetTextFont(22);
    ttext1->SetTextSize(s1);
    ttext2->SetTextFont(132);
    ttext2->SetTextSize(s2);

    ttext1->Draw();
    ttext2->Draw();
}

void histo_setup(TH1D* histo) {
    histo->GetXaxis()->SetLabelFont(label_font);
    histo->GetYaxis()->SetLabelFont(label_font);
    histo->GetXaxis()->SetLabelSize(label_font_size);
    histo->GetYaxis()->SetLabelSize(label_font_size);
    histo->GetXaxis()->SetTitleSize(titleX_font_size);
    histo->GetYaxis()->SetTitleSize(titleY_font_size);
    histo->GetYaxis()->SetTitleOffset(y_title_offset);
    histo->GetXaxis()->SetTitleOffset(x_title_offset + 0.1);
    histo->GetXaxis()->SetLabelOffset(x_label_offset);
    histo->GetYaxis()->SetLabelOffset(y_label_offset);
    histo->GetYaxis()->SetNdivisions(504);
    histo->GetYaxis()->SetMaxDigits(1);
    histo->GetXaxis()->CenterTitle(true);
    histo->GetYaxis()->CenterTitle(true);
    histo->GetXaxis()->SetTitleFont(132);
    histo->GetYaxis()->SetTitleFont(132);
}

void set_yaxis_title(TH1D* h, const double text_size) {
    double bin_width = h->GetBinWidth(0u);
    std::string str = std::to_string(bin_width);
    str.erase(str.find_last_not_of('0') + 1, std::string::npos);
    str.erase(str.find_last_not_of('.') + 1, std::string::npos);
    std::string y_axis_title = "Counts / (" + str + ")";
    h->GetYaxis()->SetTitle(y_axis_title.c_str());
    h->GetYaxis()->SetTitleSize(text_size);
}

void draw_histogram(const std::string root_name, const int num) {

    const std::string rect_title = "Bound-to-bound transport";
    const std::string geom_title = "RKN with the ODD magnetic field and CsI";

    TFile* f = TFile::Open(root_name.c_str(), "read");
    TTree* t = (TTree*)f->Get("inhom_rect_material");

    auto dthetadl0_canvas =
        new TCanvas("dthetadl0_canvas", "dthetadl0_canvas", cdim[0], cdim[1]);
    dthetadl0_canvas->SetLogy();

    TPad* pad1 = new TPad("pad1", "pad1", pad_x0, pad_y0, pad_x1, pad_y1);
    pad1->Draw();
    pad1->cd();
    pad1->SetLeftMargin(110. / pad1->GetWw());
    pad1->SetBottomMargin(125. / pad1->GetWh());
    pad1->SetLogy();

    t->Draw("log10(abs(dthetadl0_E)) >> htemp(100,-14,-4)");
    TH1D* dthetadl0_hist = (TH1D*)gPad->GetPrimitive("htemp");
    histo_setup(dthetadl0_hist);
    set_yaxis_title(dthetadl0_hist, titleY_font_size);
    dthetadl0_hist->GetXaxis()->SetNdivisions(505);
    dthetadl0_hist->GetXaxis()->SetTitle(
        "log_{10}( #left|#frac{#partial#theta_{F}}{#partiall_{0I}}#right| "
        "[mm^{-1}] )");
    dthetadl0_hist->SetMaximum(maximum);
    dthetadl0_hist->SetLineColor(kOrange + 3);
    dthetadl0_hist->SetFillStyle(fill_style);
    dthetadl0_hist->SetFillColorAlpha(kOrange + 2, color_alpha);

    draw_text(x_pos1, y_pos1, y_pos2, header_text_size, geom_text_size,
              rect_title.c_str(), geom_title.c_str());

    dthetadl0_canvas->Draw();
    dthetadl0_canvas->SaveAs("bound_to_bound_dthetadl0_E_histo.pdf");

    auto dqopdqop_canvas =
        new TCanvas("dqopdqop_canvas", "dqopdqop_canvas", cdim[0], cdim[1]);
    dqopdqop_canvas->SetLogy();

    TPad* pad2 = new TPad("pad2", "pad2", pad_x0, pad_y0, pad_x1, pad_y1);
    pad2->Draw();
    pad2->cd();
    pad2->SetLeftMargin(110. / pad2->GetWw());
    pad2->SetBottomMargin(125. / pad2->GetWh());
    pad2->SetLogy();

    t->Draw("log10(dqopdqop_E) >> htemp2(100,0,1)");
    TH1D* dqopdqop_hist = (TH1D*)gPad->GetPrimitive("htemp2");
    histo_setup(dqopdqop_hist);
    set_yaxis_title(dqopdqop_hist, titleY_font_size);

    histo_setup(dqopdqop_hist);
    dqopdqop_hist->GetXaxis()->SetNdivisions(505);
    dqopdqop_hist->GetXaxis()->SetTitle(
        "log_{10}( "
        "#left|#frac{#partial#lambda_{F}}{#partial#lambda_{I}}#right| )");
    dqopdqop_hist->SetMaximum(maximum);
    dqopdqop_hist->SetLineColor(kGreen + 3);
    dqopdqop_hist->SetFillStyle(fill_style);
    dqopdqop_hist->SetFillColorAlpha(kGreen + 2, color_alpha);

    draw_text(x_pos2, y_pos1, y_pos2, header_text_size, geom_text_size,
              rect_title.c_str(), geom_title.c_str());

    dqopdqop_canvas->Draw();
    dqopdqop_canvas->SaveAs("bound_to_bound_dqopdqop_E_histo.pdf");
}

void jacobian_histogram(int num) {
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    const std::string name = "inhom_rect_material_" + std::to_string(num);
    const std::string csv_name = name + ".csv";
    const std::string root_name = "residual_histogram.root";

    std::clog << "Processing file: " << csv_name << std::endl;

    auto rdf = ROOT::RDF::FromCSV(csv_name);
    // Create root file
    rdf.Snapshot("inhom_rect_material", root_name.c_str());

    draw_histogram(root_name, num);
}
