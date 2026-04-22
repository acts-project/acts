/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// ROOT include(s).
#include <TCanvas.h>
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
const double step_title_x_offset = 1.4;
const double step_title_x = 0.175;
const double step_title_y = 0.835;
const double step_ygap = -0.0505;
const double x_label_offset = 0.015;
const double x_margin = 1.;
const double y_label_offset = 0.015;
const double rk_title_x_offset = 0.9;
const double rk_title_offset_fraction = 0.0218;
const double rk_title_y = 14.69;
const double rk_ygap = -0.83;
const double rk_header_text_size = 0.046;
const double rk_geom_text_size = 0.0362903;
const double step_header_text_size = 0.055;
const double step_geom_text_size = 0.0434028;
const double label_font_size_step = 0.055;
const double title_font_size_step = 0.055;
const double label_font_size_rk_tol = 0.046;
const double title_font_size_rk_tol = 0.046;
const int title_font = 132;
const int label_font = 132;
const int legend_font = 132;
const double pad_x0 = 0.005;
const double pad_x1 = 1.;
const double pad_y0 = 0.005;
const double pad_y1 = 1.;
const double ymin = -2.;
const double ymax = 4.;
std::array<double, 4u> ldim{0.20985, 0.533, 0.87997, 0.939};
}  // namespace

std::vector<std::string> create_labels() {

    std::vector<std::string> varI = {
        "{#partiall_{0I}}", "{#partiall_{1I}}", "{#partial#phi_{I}}",
        "{#partial#theta_{I}}", "{#partial#lambda_{I}}"};

    std::vector<std::string> varF = {
        "{#partiall_{0F}}", "{#partiall_{1F}}", "{#partial#phi_{F}}",
        "{#partial#theta_{F}}", "{#partial#lambda_{F}}"};

    std::vector<std::string> labels;

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            std::string head = "#frac";
            labels.push_back(head + varF[i] + varI[j]);
        }
    }

    return labels;
}

void draw_text(double x1, double y1, double y_delta, double s1, double s2,
               std::string t1, std::string t2) {
    TLatex* ttext1 = new TLatex(x1, y1, t1.c_str());
    TLatex* ttext2 = new TLatex(x1, y1 + y_delta, t2.c_str());
    ttext1->SetTextFont(22);
    ttext1->SetTextSize(s1);
    ttext2->SetTextFont(132);
    ttext2->SetTextSize(s2);

    ttext1->Draw();
    ttext2->Draw();
}

std::vector<std::string> create_columns(const std::string& tag) {
    std::vector<std::string> varI = {"dl0", "dl1", "dphi", "dtheta", "dqop"};

    std::vector<std::string> varF = {"dl0", "dl1", "dphi", "dtheta", "dqop"};

    std::vector<std::string> columns;

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            columns.push_back(varF[i] + varI[j] + "_" + tag);
        }
    }

    return columns;
}

std::array<double, 25u> get_means(ROOT::RDataFrame& rdf) {

    std::array<double, 25u> ret;

    auto col_names_R = create_columns("R");

    for (std::size_t i = 0u; i < 25u; i++) {
        const double residual_mean = *rdf.Mean<double>(col_names_R[i]);
        ret[i] = residual_mean;
    }

    return ret;
}

std::map<std::string, std::vector<double>> get_means(
    const std::vector<std::string> labels, const std::string tag, const int min,
    const int max, std::vector<double>& mean_step_sizes) {

    std::map<std::string, std::vector<double>> ret;

    int num = min;
    while (num <= max + 1e-3) {

        const std::string name = tag + "_" + std::to_string(num);
        const std::string csv_name = name + ".csv";

        const std::string sign = num >= 0 ? "p" : "m";
        const std::string root_name = name + ".root";

        std::clog << "Processing file: " << csv_name << std::endl;

        auto rdf = ROOT::RDF::FromCSV(csv_name);

        const std::array<double, 25u> means = get_means(rdf);

        for (unsigned int i = 0; i < 25u; i++) {
            ret[labels[i]].push_back(TMath::Log10(means[i]));
        }

        mean_step_sizes.push_back(
            TMath::Log10(*rdf.Mean<double>("average_step_size")));

        // Create root file
        rdf.Snapshot(tag.c_str(), root_name.c_str());

        num = num + 2;
    }

    return ret;
}

std::vector<double> get_x_vector(const int min, const int max) {
    std::vector<double> ret;

    int num = min;
    while (num <= max + 1e-3) {
        ret.push_back(num);
        num = num + 2;
    }

    return ret;
}

void draw_graphs(const std::string header_title, const std::string geom_title,
                 const std::vector<std::string> labels,
                 const std::vector<double> x_vec,
                 std::map<std::string, std::vector<double>> means) {

    TPad* gr_pad = new TPad("gr_pad", "gr_pad", 0, 0, 1, 1);
    gr_pad->Draw();
    gr_pad->cd();
    gr_pad->SetLeftMargin(122. / gr_pad->GetWw());
    gr_pad->SetTopMargin(50. / gr_pad->GetWh());
    gr_pad->SetBottomMargin(120. / gr_pad->GetWh());

    TGraph* gr[25];
    TMultiGraph* mg = new TMultiGraph();

    const std::array<int, 5u> marker_styles = {7, 2, 5, 27, 25};
    const std::array<double, 5u> marker_sizes = {2.135, 2.135, 2.135, 2.135,
                                                 1.49};
    const std::array<int, 5u> line_styles = {1, 3, 2, 7, 4};
    const std::array<int, 5u> hues = {kOrange + 2, kPink + 5, kBlue + 2,
                                      kCyan + 2, kGreen + 2};

    auto legend = new TLegend(ldim[0u], ldim[1u], ldim[2u], ldim[3u]);

    // legend->SetHeader(header_title.c_str());
    legend->SetHeader("");
    legend->SetNColumns(5);
    legend->SetColumnSeparation(-0.3);
    legend->SetFillStyle(0);
    legend->SetMargin(0.55);
    legend->SetTextFont(legend_font);

    for (int i = 0; i < 25; i++) {
        gr[i] = new TGraph(x_vec.size(), &x_vec[0], &means[labels[i]][0]);

        const int n = i / 5;
        const int m = i % 5;

        gr[i]->SetMarkerStyle(marker_styles[n]);
        gr[i]->SetMarkerSize(marker_sizes[n]);
        gr[i]->SetLineStyle(line_styles[m]);
        gr[i]->SetMarkerColor(hues[m]);
        gr[i]->SetLineColor(hues[m]);

        mg->Add(gr[i]);
        legend->AddEntry(gr[i], labels[i].c_str(), "lp");
    }

    mg->GetXaxis()->SetTitle("log_{10}(#font[12]{#tau} [mm])");
    mg->GetXaxis()->SetLabelOffset(-0.005);
    mg->GetXaxis()->SetTitleOffset(rk_title_x_offset);
    mg->GetXaxis()->SetTitleSize(title_font_size_rk_tol);
    mg->GetXaxis()->SetLabelSize(label_font_size_rk_tol);
    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);

    double yaxis_min = -10;
    double yaxis_max = 15;
    double yaxis_margin = 1.;

    std::clog << "Vec size: " << x_vec.size() << std::endl;
    mg->GetYaxis()->SetTitle("log_{10}(Mean of #font[12]{#Omega_{R}})");
    mg->GetYaxis()->SetTitleOffset(1.48);
    mg->GetYaxis()->SetTitleSize(title_font_size_rk_tol);
    mg->GetYaxis()->SetLabelSize(label_font_size_rk_tol);
    mg->GetYaxis()->SetRangeUser(yaxis_min - yaxis_margin,
                                 yaxis_max + yaxis_margin);
    mg->GetYaxis()->SetLabelOffset(0.01);
    mg->GetYaxis()->SetTitleFont(title_font);
    mg->GetXaxis()->SetTitleFont(title_font);

    double x_min = x_vec.front();
    double x_max = x_vec.back();

    if (x_vec.size() > 10) {
        if (x_vec.size() == 13) {
            x_margin = 2;
        }

        x_min = x_min - x_margin;
        x_max = x_max + x_margin;
        mg->GetXaxis()->SetLimits(x_min, x_max);
        mg->GetXaxis()->SetLabelSize(0);
        mg->GetXaxis()->SetTickLength(0);

        mg->Draw("APL");

        if (x_vec.size() == 11) {
            auto ga = new TGaxis(x_vec.front(), yaxis_min - yaxis_margin,
                                 x_vec.back(), yaxis_min - yaxis_margin,
                                 x_vec.front(), x_vec.back(), 405, "N");
            ga->SetLabelFont(label_font);
            ga->SetLabelSize(label_font_size_rk_tol);
            ga->SetLabelOffset(-0.0065);
            ga->Draw();
        } else if (x_vec.size() == 13) {
            auto ga = new TGaxis(x_vec.front(), yaxis_min - yaxis_margin,
                                 x_vec.back(), yaxis_min - yaxis_margin,
                                 x_vec.front(), x_vec.back(), 304, "N");
            ga->SetLabelFont(label_font);
            ga->SetLabelSize(label_font_size_rk_tol);
            ga->SetLabelOffset(-0.0065);
            ga->Draw();
        }

        mg->GetYaxis()->SetLabelSize(0);
        mg->GetYaxis()->SetTickLength(0);
        auto ga_y = new TGaxis(x_min, yaxis_min, x_min, yaxis_max, yaxis_min,
                               yaxis_max, 505, "N");
        ga_y->SetLabelFont(label_font);
        ga_y->SetLabelOffset(0.02);
        ga_y->SetLabelSize(label_font_size_rk_tol);
        ga_y->Draw();

    } else {
        mg->Draw("APL");
    }

    TLegendEntry* header =
        (TLegendEntry*)legend->GetListOfPrimitives()->First();
    header->SetTextFont(22);
    header->SetTextSize(.033);

    double rk_title_deltaX =
        (x_vec.back() - x_vec.front()) * rk_title_offset_fraction;

    if (x_vec.size() == 13) {
        rk_title_deltaX = -0.2;
    }

    legend->Draw();
    draw_text(rk_title_deltaX + x_vec.front(), rk_title_y, rk_ygap,
              rk_header_text_size, rk_geom_text_size, header_title, geom_title);
}

void draw_mean_step_size(const std::string header_title,
                         const std::string geom_title,
                         const std::vector<double>& x_vec,
                         const std::vector<double>& means) {
    TPad* step_pad =
        new TPad("step_pad", "step_pad", pad_x0, pad_y0, pad_x1, pad_y1);
    step_pad->Draw();
    step_pad->cd();

    step_pad->SetLeftMargin(93. / step_pad->GetWw());
    step_pad->SetBottomMargin(90. / step_pad->GetWh());
    step_pad->SetTopMargin(50. / step_pad->GetWh());

    TGraph* gr = new TGraph(x_vec.size(), &x_vec[0], &means[0]);

    gr->SetMarkerStyle(4);
    gr->SetMarkerSize(1.2);
    gr->GetXaxis()->SetTitle("log_{10}(#font[12]{#tau} [mm])");
    gr->GetYaxis()->SetTitle("log_{10}(Mean of avg. step size [mm])");
    gr->GetXaxis()->SetLimits(x_vec.front() - 0.5, x_vec.back() + 0.5);

    if (x_vec.size() == 13) {
        ymin = -4;
    }

    gr->GetYaxis()->SetRangeUser(ymin, ymax);
    gr->GetYaxis()->SetNdivisions(505);
    gr->GetXaxis()->SetLabelSize(label_font_size_step);
    gr->GetYaxis()->SetLabelSize(label_font_size_step);
    gr->GetXaxis()->SetTitleSize(title_font_size_step);
    gr->GetYaxis()->SetTitleSize(title_font_size_step);
    gr->GetXaxis()->SetTitleOffset(step_title_x_offset);
    gr->GetYaxis()->SetTitleOffset(1.15);
    gr->GetYaxis()->SetLabelOffset(y_label_offset);
    gr->GetYaxis()->SetTitleFont(title_font);
    gr->GetXaxis()->SetTitleFont(title_font);
    gr->GetXaxis()->SetLabelOffset(x_label_offset);
    gr->GetXaxis()->CenterTitle(true);
    gr->GetYaxis()->CenterTitle(true);
    gr->GetXaxis()->SetLabelFont(label_font);
    gr->GetYaxis()->SetLabelFont(label_font);

    if (x_vec.size() > 10) {
        if (x_vec.size() == 13) {
            x_margin = 2;
        }

        gr->GetXaxis()->SetLimits(x_vec.front() - x_margin,
                                  x_vec.back() + x_margin);
        gr->GetXaxis()->SetLabelSize(0);
        gr->GetXaxis()->SetTickLength(0);

        gr->Draw("APL");

        if (x_vec.size() == 11) {
            auto ga = new TGaxis(x_vec.front(), ymin, x_vec.back(), ymin,
                                 x_vec.front(), x_vec.back(), 405, "N");
            ga->SetLabelFont(label_font);
            ga->SetLabelSize(label_font_size_step);
            ga->SetLabelOffset(x_label_offset);
            ga->Draw();
        } else if (x_vec.size() == 13) {
            auto ga = new TGaxis(x_vec.front(), ymin, x_vec.back(), ymin,
                                 x_vec.front(), x_vec.back(), 304, "N");
            ga->SetLabelFont(label_font);
            ga->SetLabelSize(label_font_size_step);
            ga->SetLabelOffset(x_label_offset);
            ga->Draw();
        }

    } else {
        gr->Draw();
    }

    TPad* text_pad =
        new TPad("text_pad", "text_pad", pad_x0, pad_y0, pad_x1, pad_y1);
    text_pad->SetFillStyle(4000);
    text_pad->Draw();
    text_pad->cd();

    draw_text(step_title_x, step_title_y, step_ygap, step_header_text_size,
              step_geom_text_size, header_title, geom_title);
}

// ROOT Script for jacboain file reading
void rk_tolerance_comparison(int min, int max) {
    gStyle->SetOptTitle(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendTextSize(0.0333);

    const std::array<float, 2> cdim1{800, 1300};
    const std::array<float, 2> cdim2{700, 600};

    auto labels = create_labels();

    const auto x_vec = get_x_vector(min, max);

    const std::string rect_header = "Bound-to-bound transport";
    const std::string wire_header = "Perigee-to-perigee transport";
    const std::string geom_header = "RKN with the ODD magnetic field and CsI";

    /************************
     *  Rectangular
     * **********************/

    const std::string rect_pdf = "bound_to_bound_rk_tolerance.pdf";

    auto rect_canvas =
        new TCanvas("rect_canvas", "rect_canvas", cdim1[0], cdim1[1]);

    std::vector<double> rect_mean_step_sizes;
    const auto rect_y_means = get_means(labels, "inhom_rect_material", min, max,
                                        rect_mean_step_sizes);
    draw_graphs(rect_header, geom_header, labels, x_vec, rect_y_means);

    rect_canvas->SaveAs(rect_pdf.c_str());

    auto rect_canvas2 =
        new TCanvas("rect_canvas2", "rect_canvas2", cdim2[0], cdim2[1]);
    const std::string rect_mean_step_pdf = "bound_to_bound_mean_step_size.pdf";

    draw_mean_step_size(rect_header, geom_header, x_vec, rect_mean_step_sizes);
    rect_canvas2->SaveAs(rect_mean_step_pdf.c_str());

    /************************
     *  Wire
     * **********************/

    const std::string wire_pdf = "perigee_to_perigee_rk_tolerance.pdf";

    auto wire_canvas =
        new TCanvas("wire_canvas", "wire_canvas", cdim1[0], cdim1[1]);
    std::vector<double> wire_mean_step_sizes;
    const auto wire_y_means = get_means(labels, "inhom_wire_material", min, max,
                                        wire_mean_step_sizes);
    draw_graphs(wire_header, geom_header, labels, x_vec, wire_y_means);

    wire_canvas->SaveAs(wire_pdf.c_str());

    auto wire_canvas2 =
        new TCanvas("wire_canvas2", "wire_canvas2", cdim2[0], cdim2[1]);
    const std::string wire_mean_step_pdf =
        "perigee_to_perigee_mean_step_size.pdf";

    draw_mean_step_size(wire_header, geom_header, x_vec, wire_mean_step_sizes);
    wire_canvas2->SaveAs(wire_mean_step_pdf.c_str());
}
