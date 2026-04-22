/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// ROOT include(s).
#include <Math/ProbFuncMathCore.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>

#include <ROOT/RCsvDS.hxx>
#include <ROOT/RDataFrame.hxx>

// System include(s).
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

namespace {
const double x_pos = 0.205;
const double title_x = x_pos;
const double title_y = 0.8197;
const double y_gap = -0.0505;
const double header_text_size = 0.055;
const double geom_text_size = 0.0434028;

const double pull_fit_title_x = x_pos;
const double pull_fit_title_y = 0.700;
const double pval_fit_title_x = x_pos;
const double pval_fit_title_y = 0.700;
const double gaus_fit_par_x = x_pos;
const double number_offset = 0.125;
const double gaus_fit_par_y = pull_fit_title_y - 0.065;
const double const_fit_par_x = x_pos;
const double const_fit_par_y = pval_fit_title_y - 0.0459;
const double tolerance_x = 0.7;
const double tolerance_y = 0.67;
const double pull_text_size = 0.0434028;
const double pval_text_size = 0.0434028;
const double pad_x0 = 0.00;
const double pad_x1 = 1.;
const double pad_y0 = 0.00;
const double pad_y1 = 1.;
const int label_font = 132;
const double label_font_size = 0.055;
const double titleX_font_size = 0.055;
const double titleY_font_size = 0.055;
const double x_title_offset = 1.25;
const double y_title_offset = 1.34;
const double y_title_offset_pval = 0.9;
const double x_label_offset = 0.015;
const double y_label_offset = 0.015;
const double pull_min = -6.5;
const double pull_max = 6.5;

}  // namespace

auto get_tree(std::string name) {

    const std::string csv_name = name + ".csv";
    const std::string root_name = name + ".root";

    auto rdf = ROOT::RDF::FromCSV(csv_name);

    // Create root file
    rdf.Snapshot(name, root_name);

    auto f = TFile::Open(root_name.c_str(), "update");
    auto t = (TTree*)f->Get(name.c_str());

    return t;
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

std::pair<std::array<double, 3u>, std::array<double, 3u>> fit_pull(
    TH1D* h_pull, std::array<double, 14u>& arr_pull) {

    // Function used for the fit.
    TF1 gaus{"gaus", "gaus", pull_min, pull_max};
    double fit_par[3];
    double fit_par_error[3];

    // Set the mean seed to 0
    gaus.SetParameters(1, 0);
    gaus.SetParLimits(1, -1., 1.);
    // Set the standard deviation seed to 1
    gaus.SetParameters(2, 1.0);
    gaus.SetParLimits(2, 0.5, 2.);

    auto res = h_pull->Fit("gaus", "Q0S");
    gaus.GetParameters(&fit_par[0]);

    std::array<double, 3u> par{fit_par[0], fit_par[1], fit_par[2]};
    std::array<double, 3u> error;
    error[0] = gaus.GetParError(0);
    error[1] = gaus.GetParError(1);
    error[2] = gaus.GetParError(2);

    arr_pull[0u] = fit_par[0u];
    arr_pull[1u] = fit_par[1u];
    arr_pull[2u] = fit_par[2u];
    arr_pull[3u] = error[0u];
    arr_pull[4u] = error[1u];
    arr_pull[5u] = error[2u];
    arr_pull[6u] = res->Ndf();
    arr_pull[7u] = res->Chi2();
    arr_pull[8u] = arr_pull[6u] / arr_pull[7u];
    arr_pull[9u] = ROOT::Math::chisquared_cdf_c(arr_pull[7u], arr_pull[6u]);

    return {par, error};
}

std::pair<double, double> fit_pval(TH1D* h_pval) {

    // Function used for the fit.
    TF1 unif{"uniform", "[0]", 0.f, 1.f};
    double fit_par[1];

    auto res = h_pval->Fit("uniform", "Q0S");
    unif.GetParameters(&fit_par[0]);
    double error = unif.GetParError(0);

    return {fit_par[0], error};
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

void set_xaxis_title(TH1D* h, const double text_size) {

    std::string x_axis_title;

    const TString h_name = h->GetName();

    if (h_name.Contains("l0")) {
        x_axis_title = "PL(l_{0F})";
    } else if (h_name.Contains("l1")) {
        x_axis_title = "PL(l_{1F})";
    } else if (h_name.Contains("phi")) {
        x_axis_title = "PL(#phi_{F})";
    } else if (h_name.Contains("theta")) {
        x_axis_title = "PL(#theta_{F})";
    } else if (h_name.Contains("qop")) {
        x_axis_title = "PL(#lambda_{F})";
    } else if (h_name.Contains("pval")) {
        x_axis_title = "p-value";
    }

    h->GetXaxis()->SetTitle(x_axis_title.c_str());
    h->GetXaxis()->SetTitleSize(text_size);
}

void draw_fit_title(const std::string title, const double x, const double y,
                    const double text_size) {

    TLatex* ttext = new TLatex(x, y, title.c_str());
    ttext->SetTextFont(22);
    ttext->SetTextSize(text_size);
    ttext->Draw();
    gPad->cd();
}

void draw_gaus_fit_par(const std::array<double, 3u>& fit_par,
                       const std::array<double, 3u>& fit_par_error,
                       const double x, const double y, const double text_size) {
    TLatex* ttext = new TLatex(x, y, "#splitline{Mean}{Std Dev}");
    ttext->SetTextFont(132);
    ttext->SetTextSize(text_size);
    ttext->Draw();

    std::stringstream mean_stream;
    mean_stream << std::fixed << std::setprecision(3) << fit_par[1] << " #pm "
                << fit_par_error[1];
    std::stringstream sigma_stream;
    sigma_stream << std::fixed << std::setprecision(3) << fit_par[2] << " #pm "
                 << fit_par_error[2];

    TLatex* ttext2 = new TLatex(x + number_offset, y,
                                "#splitline{" + TString(mean_stream.str()) +
                                    "}{" + TString(sigma_stream.str()) + "}");
    ttext2->SetTextFont(132);
    ttext2->SetTextSize(text_size);
    ttext2->Draw();
}

void draw_const_fit_par(const double fit_par, const double fit_par_error,
                        const double x, const double y,
                        const double text_size) {
    std::stringstream val_stream;
    val_stream << std::fixed << std::setprecision(3) << fit_par << " #pm "
               << fit_par_error;

    TLatex* ttext = new TLatex(x, y, "Value  " + TString(val_stream.str()));
    ttext->SetTextFont(132);
    ttext->SetTextSize(text_size);
    ttext->Draw();
}

void draw_tolerance(const double log10_rk_tolerance, const double x,
                    const double y) {
    std::stringstream val_stream;
    val_stream << "#tau = 10^{" << int(log10_rk_tolerance) << "} mm";

    TLatex* ttext = new TLatex(x, y, TString(val_stream.str()));
    ttext->SetTextFont(132);
    ttext->SetTextSize(0.052);
    ttext->Draw();
}

void draw_pull(TH1D* h_pull, const std::string& header_text,
               const std::string& geom_text, const double log10_rk_tol,
               std::array<double, 14u>& arr_pull) {

    TPad* pull_pad =
        new TPad("pull_pad", "pull_pad", pad_x0, pad_y0, pad_x1, pad_y1);
    pull_pad->SetLogy();
    pull_pad->Draw();
    pull_pad->cd();

    pull_pad->SetLeftMargin(110. / pull_pad->GetWw());
    pull_pad->SetBottomMargin(95. / pull_pad->GetWh());

    auto fit_res = fit_pull(h_pull, arr_pull);
    auto fit_pars = fit_res.first;
    auto fit_errors = fit_res.second;

    set_xaxis_title(h_pull, pull_text_size);
    set_yaxis_title(h_pull, pull_text_size);
    const double y_axis_max = h_pull->GetEntries() * 50.;
    h_pull->GetYaxis()->SetRangeUser(0.5f, y_axis_max);
    h_pull->GetXaxis()->SetLabelFont(label_font);
    h_pull->GetYaxis()->SetLabelFont(label_font);
    h_pull->GetXaxis()->SetLabelSize(label_font_size);
    h_pull->GetYaxis()->SetLabelSize(label_font_size);
    h_pull->GetXaxis()->SetTitleSize(titleX_font_size);
    h_pull->GetYaxis()->SetTitleSize(titleY_font_size);
    h_pull->GetYaxis()->SetTitleOffset(y_title_offset);
    h_pull->GetXaxis()->SetTitleOffset(x_title_offset + 0.1);
    h_pull->GetXaxis()->SetLabelOffset(x_label_offset);
    h_pull->GetYaxis()->SetLabelOffset(y_label_offset);
    h_pull->GetYaxis()->SetNdivisions(504);
    h_pull->GetYaxis()->SetMaxDigits(1);
    h_pull->GetXaxis()->CenterTitle(true);
    h_pull->GetYaxis()->CenterTitle(true);
    h_pull->GetXaxis()->SetTitleFont(132);
    h_pull->GetYaxis()->SetTitleFont(132);
    h_pull->Draw();
    TF1* gaus = new TF1(h_pull->GetName(), "gaus", pull_min, pull_max);
    gaus->SetParameters(fit_pars[0], fit_pars[1], fit_pars[2]);
    gaus->Draw("same");

    TPad* text_pad = new TPad("gaus_text_pad", "gaus_text_pad", pad_x0, pad_y0,
                              pad_x1, pad_y1);
    text_pad->SetFillStyle(4000);
    text_pad->Draw();
    text_pad->cd();

    draw_text(title_x, title_y, y_gap, header_text_size, geom_text_size,
              header_text.c_str(), geom_text.c_str());
    draw_fit_title("Gaussian fit", pull_fit_title_x, pull_fit_title_y,
                   pull_text_size);
    draw_gaus_fit_par(fit_pars, fit_errors, gaus_fit_par_x, gaus_fit_par_y,
                      pull_text_size);
    // draw_tolerance(log10_rk_tol, tolerance_x, tolerance_y);
}

void draw_pval(TH1D* h_pval, const std::string& header_text,
               const std::string& geom_text, const double log10_rk_tol) {

    TPad* pval_pad =
        new TPad("pval_pad", "pval_pad", pad_x0, pad_y0, pad_x1, pad_y1);
    pval_pad->Draw();
    pval_pad->cd();

    pval_pad->SetLeftMargin(110. / pval_pad->GetWw());
    pval_pad->SetBottomMargin(95. / pval_pad->GetWh());

    auto fit_res = fit_pval(h_pval);
    auto fit_par = fit_res.first;
    auto fit_error = fit_res.second;
    set_xaxis_title(h_pval, pval_text_size);
    set_yaxis_title(h_pval, pval_text_size);
    const double y_axis_max = 2. * h_pval->GetEntries() / h_pval->GetNbinsX();
    h_pval->GetYaxis()->SetRangeUser(0.f, y_axis_max);
    h_pval->GetXaxis()->SetLabelFont(label_font);
    h_pval->GetYaxis()->SetLabelFont(label_font);
    h_pval->GetXaxis()->SetLabelSize(label_font_size);
    h_pval->GetYaxis()->SetLabelSize(label_font_size);
    h_pval->GetXaxis()->SetTitleSize(titleX_font_size);
    h_pval->GetYaxis()->SetTitleSize(titleY_font_size);
    h_pval->GetYaxis()->SetTitleOffset(y_title_offset_pval);
    h_pval->GetXaxis()->SetTitleOffset(x_title_offset);
    h_pval->GetYaxis()->SetLabelOffset(y_label_offset);
    h_pval->GetXaxis()->SetLabelOffset(x_label_offset);
    h_pval->GetXaxis()->SetNdivisions(505);
    h_pval->GetYaxis()->SetMaxDigits(2);
    h_pval->GetYaxis()->SetNdivisions(505);
    h_pval->GetYaxis()->SetDecimals();
    h_pval->GetXaxis()->CenterTitle(true);
    h_pval->GetYaxis()->CenterTitle(true);
    h_pval->GetXaxis()->SetTitleFont(132);
    h_pval->GetYaxis()->SetTitleFont(132);

    h_pval->Draw();
    TLine* line = new TLine(0.f, fit_par, 1.f, fit_par);
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    line->Draw();

    TPad* text_pad = new TPad("const_text_pad", "const_text_pad", 0, 0, 1, 1);
    text_pad->SetFillStyle(4000);
    text_pad->Draw();
    text_pad->cd();

    draw_text(title_x, title_y, y_gap, header_text_size, geom_text_size,
              header_text.c_str(), geom_text.c_str());
    draw_fit_title("Constant fit", pval_fit_title_x, pval_fit_title_y,
                   pval_text_size);
    draw_const_fit_par(fit_par, fit_error, const_fit_par_x, const_fit_par_y,
                       pval_text_size);
    // draw_tolerance(log10_rk_tol, tolerance_x, tolerance_y);
}

std::string to_pdf(const std::string& name) {
    return name + ".pdf";
}

void read_tree(TTree* t, const std::string& tag,
               const std::string& header_title, const std::string& geom_title) {
    const std::array<float, 2> cdim1{700, 600};
    const std::array<float, 2> cdim2{700, 600};

    int n_bins = (pull_max - pull_min) / 0.1;

    double pull_l0;
    double pull_l1;
    double pull_phi;
    double pull_theta;
    double pull_qop;
    double chi2;
    double log10_rk_tolerance;
    double p_value;

    t->SetBranchAddress("pull_l0", &pull_l0);
    t->SetBranchAddress("pull_l1", &pull_l1);
    t->SetBranchAddress("pull_phi", &pull_phi);
    t->SetBranchAddress("pull_theta", &pull_theta);
    t->SetBranchAddress("pull_qop", &pull_qop);
    t->SetBranchAddress("chi2", &chi2);
    t->SetBranchAddress("log10_rk_tolerance", &log10_rk_tolerance);
    auto b_pval = t->Branch("p_value", &p_value, "p_value/D");

    std::string l0_name = tag + "_pull_l0";
    std::string l1_name = tag + "_pull_l1";
    std::string phi_name = tag + "_pull_phi";
    std::string theta_name = tag + "_pull_theta";
    std::string qop_name = tag + "_pull_qop";
    std::string chi2_name = tag + "_chi2";
    std::string pval_name = tag + "_pval";

    TH1D* h_l0 =
        new TH1D(l0_name.c_str(), l0_name.c_str(), n_bins, pull_min, pull_max);
    TH1D* h_l1 =
        new TH1D(l1_name.c_str(), l1_name.c_str(), n_bins, pull_min, pull_max);
    TH1D* h_phi = new TH1D(phi_name.c_str(), phi_name.c_str(), n_bins, pull_min,
                           pull_max);
    TH1D* h_theta = new TH1D(theta_name.c_str(), theta_name.c_str(), n_bins,
                             pull_min, pull_max);
    TH1D* h_qop = new TH1D(qop_name.c_str(), qop_name.c_str(), n_bins, pull_min,
                           pull_max);
    TH1D* h_chi2 =
        new TH1D(chi2_name.c_str(), chi2_name.c_str(), 50, 0.f, 50.f);
    TH1D* h_pval =
        new TH1D(pval_name.c_str(), pval_name.c_str(), 50, 0.f, 1.0f);

    // N, mean, stddev, N_error, mean_error, stddev_error, ndf, chi2, ndf/chi2,
    // pval, N_4sigma, N_4sigma_fraction, N_expected_4sigma,
    // N_4sigma/N_expected_4sigma
    std::array<double, 14u> finfo_l0;
    std::array<double, 14u> finfo_l1;
    std::array<double, 14u> finfo_phi;
    std::array<double, 14u> finfo_theta;
    std::array<double, 14u> finfo_qop;
    finfo_l0[10u] = 0;
    finfo_l1[10u] = 0;
    finfo_phi[10u] = 0;
    finfo_theta[10u] = 0;
    finfo_qop[10u] = 0;

    unsigned int n_outliers{0u};

    // Fill the histograms
    for (int i = 0; i < t->GetEntries(); i++) {
        t->GetEntry(i);
        h_l0->Fill(pull_l0);
        h_l1->Fill(pull_l1);
        h_phi->Fill(pull_phi);
        h_theta->Fill(pull_theta);
        h_qop->Fill(pull_qop);
        h_chi2->Fill(chi2);
        p_value = ROOT::Math::chisquared_cdf_c(chi2, 5.f);
        h_pval->Fill(p_value);
        b_pval->Fill();

        if ((pull_l0 < pull_min) || (pull_l0 > pull_max) ||
            (pull_l1 < pull_min) || (pull_l1 > pull_max) ||
            (pull_phi < pull_min) || (pull_phi > pull_max) ||
            (pull_theta < pull_min) || (pull_theta > pull_max) ||
            (pull_qop < pull_min) || (pull_qop > pull_max)) {
            n_outliers++;
        }

        if (abs(pull_l0) > 4) {
            finfo_l0[10u]++;
        }
        if (abs(pull_l1) > 4) {
            finfo_l1[10u]++;
        }
        if (abs(pull_phi) > 4) {
            finfo_phi[10u]++;
        }
        if (abs(pull_theta) > 4) {
            finfo_theta[10u]++;
        }
        if (abs(pull_qop) > 4) {
            finfo_qop[10u]++;
        }
    }

    t->Write("", TObject::kOverwrite);

    std::clog << tag << " n outliers: " << n_outliers << std::endl;

    finfo_l0[11u] = double(finfo_l0[10u]) / double(t->GetEntries());
    finfo_l1[11u] = double(finfo_l1[10u]) / double(t->GetEntries());
    finfo_phi[11u] = double(finfo_phi[10u]) / double(t->GetEntries());
    finfo_theta[11u] = double(finfo_theta[10u]) / double(t->GetEntries());
    finfo_qop[11u] = double(finfo_qop[10u]) / double(t->GetEntries());

    finfo_l0[12u] = t->GetEntries() * 0.000063342484;
    finfo_l1[12u] = t->GetEntries() * 0.000063342484;
    finfo_phi[12u] = t->GetEntries() * 0.000063342484;
    finfo_theta[12u] = t->GetEntries() * 0.000063342484;
    finfo_qop[12u] = t->GetEntries() * 0.000063342484;

    finfo_l0[13u] = finfo_l0[10u] / finfo_l0[12u];
    finfo_l1[13u] = finfo_l1[10u] / finfo_l1[12u];
    finfo_phi[13u] = finfo_phi[10u] / finfo_phi[12u];
    finfo_theta[13u] = finfo_theta[10u] / finfo_theta[12u];
    finfo_qop[13u] = finfo_qop[10u] / finfo_qop[12u];

    auto c_l0 =
        new TCanvas(h_l0->GetName(), h_l0->GetName(), cdim1[0], cdim1[1]);
    c_l0->SetLogy();

    draw_pull(h_l0, header_title, geom_title, log10_rk_tolerance, finfo_l0);

    c_l0->SaveAs(to_pdf(l0_name).c_str());

    auto c_l1 =
        new TCanvas(h_l1->GetName(), h_l1->GetName(), cdim1[0], cdim1[1]);
    c_l1->SetLogy();
    draw_pull(h_l1, header_title, geom_title, log10_rk_tolerance, finfo_l1);
    c_l1->SaveAs(to_pdf(l1_name).c_str());

    auto c_phi =
        new TCanvas(h_phi->GetName(), h_phi->GetName(), cdim1[0], cdim1[1]);
    c_phi->SetLogy();
    draw_pull(h_phi, header_title, geom_title, log10_rk_tolerance, finfo_phi);
    c_phi->SaveAs(to_pdf(phi_name).c_str());

    auto c_theta =
        new TCanvas(h_theta->GetName(), h_theta->GetName(), cdim1[0], cdim1[1]);
    c_theta->SetLogy();

    draw_pull(h_theta, header_title, geom_title, log10_rk_tolerance,
              finfo_theta);
    c_theta->SaveAs(to_pdf(theta_name).c_str());

    auto c_qop =
        new TCanvas(h_qop->GetName(), h_qop->GetName(), cdim1[0], cdim1[1]);
    c_qop->SetLogy();
    draw_pull(h_qop, header_title, geom_title, log10_rk_tolerance, finfo_qop);
    c_qop->SaveAs(to_pdf(qop_name).c_str());

    auto c_chi2 =
        new TCanvas(h_chi2->GetName(), h_chi2->GetName(), cdim1[0], cdim1[1]);
    h_chi2->Draw();
    c_chi2->SaveAs(to_pdf(chi2_name).c_str());

    auto c_pval =
        new TCanvas(h_pval->GetName(), h_pval->GetName(), cdim2[0], cdim2[1]);
    draw_pval(h_pval, header_title, geom_title, log10_rk_tolerance);
    c_pval->SaveAs(to_pdf(pval_name).c_str());

    std::ofstream f;

    std::string fname = tag + ".txt";
    f.open(fname.c_str());

    f << "N, mean, stddev, N_error, mean_error, stddev_error, ndf, chi2, "
         "ndf/chi2, pval, N_4sigma, N_4sigma_fraction, N_4sigma_expected, "
         "N_4sigma/N_4sigma_expected";
    f << std::endl;

    for (const auto& n : finfo_l0) {
        f << n << ",";
    }
    f << std::endl;

    for (const auto& n : finfo_l1) {
        f << n << ",";
    }
    f << std::endl;

    for (const auto& n : finfo_phi) {
        f << n << ",";
    }
    f << std::endl;

    for (const auto& n : finfo_theta) {
        f << n << ",";
    }
    f << std::endl;

    for (const auto& n : finfo_qop) {
        f << n << ",";
    }
    f << std::endl;

    f.close();
}

// ROOT Script for covariance file reading
void covariance_validation() {

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    const std::string rect_title = "Bound-to-bound transport";
    const std::string wire_title = "Perigee-to-perigee transport";
    const std::string geom_title = "RKN with the ODD magnetic field and CsI";

    /************************
     *  Rectangular
     * **********************/

    std::string rect_name = "rect_cov_transport";
    auto rect_tree = get_tree(rect_name);
    read_tree(rect_tree, "bound_to_bound", rect_title, geom_title);

    /************************
     *  Wire
     * **********************/

    std::string wire_name = "wire_cov_transport";
    auto wire_tree = get_tree(wire_name);
    read_tree(wire_tree, "perigee_to_perigee", wire_title, geom_title);
}
