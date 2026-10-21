#include <fstream>
#include <iostream>

#include "TEfficiency.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TProfile.h"

void process_teff(TEfficiency *eff, std::string tree, std::string out_dir) {
    const TH1 *hist = eff->GetTotalHistogram();

    int num_bins = hist->GetNbinsX();

    std::ofstream myfile;
    myfile.open(out_dir + "/" + tree + ".csv");

    myfile << "bin_left,bin_right,ntotal,efficiency,err_low,err_high"
           << std::endl;

    for (int i = 1; i <= num_bins; ++i) {
        myfile << hist->GetBinLowEdge(i) << ","
               << (hist->GetBinLowEdge(i) + hist->GetBinWidth(i)) << ","
               << hist->GetBinContent(i) << "," << eff->GetEfficiency(i) << ","
               << eff->GetEfficiencyErrorLow(i) << ","
               << eff->GetEfficiencyErrorUp(i) << std::endl;
    }

    myfile.close();
}

void process_profile(TProfile *hist, std::string tree, std::string out_dir) {
    int num_bins = hist->GetNbinsX();

    std::ofstream myfile;
    myfile.open(out_dir + "/" + tree + ".csv");

    myfile << "bin_left,bin_right,value,err_low,err_high" << std::endl;

    for (int i = 1; i <= num_bins; ++i) {
        myfile << hist->GetBinLowEdge(i) << ","
               << (hist->GetBinLowEdge(i) + hist->GetBinWidth(i)) << ","
               << hist->GetBinContent(i) << "," << hist->GetBinErrorLow(i)
               << "," << hist->GetBinErrorUp(i) << std::endl;
    }

    myfile.close();
}

void process_th1(TH1 *hist, std::string tree, std::string out_dir) {
    int num_bins = hist->GetNbinsX();

    std::ofstream myfile;
    myfile.open(out_dir + "/" + tree + ".csv");

    myfile << "bin_left,bin_right,ntotal,err_low,err_high" << std::endl;

    for (int i = 1; i <= num_bins; ++i) {
        myfile << hist->GetBinLowEdge(i) << ","
               << (hist->GetBinLowEdge(i) + hist->GetBinWidth(i)) << ","
               << hist->GetBinContent(i) << "," << hist->GetBinErrorLow(i)
               << "," << hist->GetBinErrorUp(i) << std::endl;
    }

    myfile.close();
}

void read_and_dump(TFile *file, std::string tree, std::string out_dir) {
    TEfficiency *eff = dynamic_cast<TEfficiency *>(file->Get(tree.c_str()));

    if (eff != nullptr) {
        std::cout << tree << " is a TEfficiency" << std::endl;
        process_teff(eff, tree, out_dir);
        return;
    }

    TProfile *prof = dynamic_cast<TProfile *>(file->Get(tree.c_str()));

    if (prof != nullptr) {
        std::cout << tree << " is a TProfile" << std::endl;
        process_profile(prof, tree, out_dir);
        return;
    }

    TH1 *th1 = dynamic_cast<TH1 *>(file->Get(tree.c_str()));

    if (th1 != nullptr) {
        std::cout << tree << " is a TH1" << std::endl;
        process_th1(th1, tree, out_dir);
        return;
    }
}

int main(int argc, char **argv) {
    if (argc < 3) {
        throw std::runtime_error("Too few args.");
    }

    std::unique_ptr<TFile> file(TFile::Open(argv[1]));

    std::vector<std::string> trees{"seeding_trackeff_vs_pT",
                                   "seeding_trackeff_vs_eta",
                                   "seeding_trackeff_vs_phi",
                                   "finding_trackeff_vs_pT",
                                   "finding_trackeff_vs_eta",
                                   "finding_trackeff_vs_phi",
                                   "pval",
                                   "ndf",
                                   "completeness",
                                   "purity",
                                   "finding_nDuplicated_vs_eta",
                                   "finding_nFakeTracks_vs_eta",
                                   "res_d0",
                                   "res_z0",
                                   "res_phi",
                                   "res_qop",
                                   "res_qopT",
                                   "res_qopz",
                                   "res_theta",
                                   "pull_d0",
                                   "pull_z0",
                                   "pull_phi",
                                   "pull_qop",
                                   "pull_qopT",
                                   "pull_qopz",
                                   "pull_theta"};

    for (const auto &tree : trees) {
        read_and_dump(file.get(), tree, argv[2]);
    }

    return 0;
}
