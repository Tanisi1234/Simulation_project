#include <TApplication.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <iostream>
#include <vector>
#include <numeric>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>

#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentAnsys123.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/ViewField.hh"

using namespace Garfield;

int main(int argc, char* argv[]) {
    TApplication app("app", &argc, argv);

    // Gas Setup
    MediumMagboltz gas("ar", 90., "co2", 10.);
    gas.SetTemperature(293.15);
    gas.SetPressure(760.);
    gas.Initialise(true);
    gas.EnablePenningTransfer(0.51, 0., "ar");
    gas.LoadIonMobility("IonMobility_Ar+_Ar.txt");

    // Field Map
    ComponentAnsys123 fm;
    fm.Initialise("ELIST.lis", "NLIST.lis", "MPLIST.lis", "PRNSOL.lis", "mm");
    fm.EnableMirrorPeriodicityX();
    fm.EnableMirrorPeriodicityY();
    fm.SetGas(&gas);

    constexpr double pitch = 0.014;

    // Field + Mesh Plot
    ViewField fieldView(&fm);
    ViewFEMesh meshView(&fm);

    TCanvas* cf = new TCanvas("cf", "Electric Potential", 750, 650);
    fieldView.SetCanvas(cf);
    fieldView.SetPlane(0, -1, 0, 0, 0, 0);
    fieldView.SetArea(-0.5 * pitch, -0.02, 0.5 * pitch, 0.02);
    fieldView.SetVoltageRange(-160., 160.);
    fieldView.PlotContour();

    meshView.SetCanvas(cf);
    meshView.SetPlane(0, -1, 0, 0, 0, 0);
    meshView.SetArea(-0.5 * pitch, -0.02, 0.5 * pitch, 0.02);
    meshView.SetFillMesh(true);
    meshView.SetColor(2, kGray);
    meshView.Plot(true);

    // Sensor and Simulation Setup
    Sensor sensor(&fm);
    sensor.SetArea(-5 * pitch, -5 * pitch, -0.01, 5 * pitch, 5 * pitch, 0.025);  \\ the x0 y0 randomly generated 
    AvalancheMicroscopic aval(&sensor);
    AvalancheMC drift(&sensor);
    drift.SetDistanceSteps(2.e-4);

    ViewDrift driftView;
    aval.EnablePlotting(&driftView, 10);

    // Statistics
    std::vector<double> absGains, effGains;
    std::size_t nTotalIons = 0;
    std::size_t nBFIons = 0;

    constexpr std::size_t nEvents = 100;

    // Main Simulation Loop
    for (std::size_t i = 0; i < nEvents; ++i) {
        const double x0 = -0.5 * pitch + RndmUniform() * pitch;
        const double y0 = -0.5 * pitch + RndmUniform() * pitch;

        aval.AvalancheElectron(x0, y0, 0.02, 0., 0.1, 0., 0., 0.);  \\start point of electons= 0.2cm 

        // Absolute Gain
        int ne = 0, ni = 0;
        aval.GetAvalancheSize(ne, ni);
        absGains.push_back(static_cast<double>(ne));

        // Effective Gain + IBF
        int nEff = 0;
        for (const auto& electron : aval.GetElectrons()) {
            if (electron.path.back().z < -0.009) {   \\ readout threshold 
                
                ++nEff;
            }

            const auto& birth = electron.path[0];
            drift.DriftIon(birth.x, birth.y, birth.z, birth.t);  \\ birth point of ions 
                        ++nTotalIons;

            if (!drift.GetIons().empty()) {
                if (drift.GetIons().front().path.back().z > 0.012) {  \\ backflow ion threshold 
                    ++nBFIons;
                }
            }
        }
        effGains.push_back(static_cast<double>(nEff));
    }

    // Calculate Statistics
    auto meanStd = [](const std::vector<double>& v) {
        if (v.empty()) return std::make_pair(0.0, 0.0);
        double sum = std::accumulate(v.begin(), v.end(), 0.0);
        double mean = sum / v.size();
        double sqsum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
        double stdev = std::sqrt(std::abs(sqsum / v.size() - mean * mean));
        return std::make_pair(mean, stdev);
    };

    auto [absMean, absStd] = meanStd(absGains);
    auto [effMean, effStd] = meanStd(effGains);
    
    // Calculate Statistical Error (Sigma / sqrt(N))
    double absError = absStd / std::sqrt(nEvents);
    double effError = effStd / std::sqrt(nEvents);
    double ibfFrac = nTotalIons > 0 ? static_cast<double>(nBFIons) / nTotalIons : 0.0;

    // Print minimal results
    std::cout << "GEM Simulation Results\n";
    std::cout << "Events          : " << nEvents << "\n";
    std::cout << "Absolute Gain   : " << absMean << " ± " << absError << "\n";
    std::cout << "Effective Gain  : " << effMean << " ± " << effError << "\n";
    std::cout << "IBF Fraction    : " << ibfFrac << "\n";




    // gain plots
    TCanvas* cGain = new TCanvas("cGain", "Gain Distribution", 1200, 500);
    cGain->Divide(2, 1);
    gStyle->SetOptStat(1110); // Entries, Mean, RMS Underflow Over flows


    cGain->cd(1);
    TH1D* hAbs = new TH1D("hAbs", "Absolute Gain;Total Electrons;Events", 30, 0, 30);
    for (double g : absGains) hAbs->Fill(g);
    hAbs->SetFillColor(kBlue-9);
    hAbs->SetLineColor(kBlue+2);
    hAbs->Draw();

    cGain->cd(2);
    TH1D* hEff = new TH1D("hEff", "Effective Gain;Collected Electrons;Events", 30, 0, 30);
    for (double g : effGains) hEff->Fill(g);
    hEff->SetFillColor(kRed-9);
    hEff->SetLineColor(kRed+2);
    hEff->Draw();

    // ==========================================================
    // DATA FROM SCANS (0% to 90% Ar)
    // ==========================================================
    const int n = 10; 
    double x_ar[]   = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90};
    // Absolute Gain values from your results
    double y_gain[] = {1.52, 1.41, 1.63, 1.79, 2.51, 3.29, 3.69, 5.07, 8.17, 22.42};
    // Statistical errors (sigma/sqrt(N))
    double e_gain[] = {0.08, 0.09, 0.12, 0.14, 0.19, 0.24, 0.31, 0.44, 0.80, 1.75};
    // IBF values from your results
    double y_ibf[]  = {0.79, 0.74, 0.69, 0.64, 0.56, 0.46, 0.46, 0.40, 0.33, 0.26};

    // ==========================================================
    // TREND CURVES
    // ==========================================================
    TCanvas *cTrends = new TCanvas("cTrends", "GEM Performance Trends", 1200, 500);
    cTrends->Divide(2, 1);

    // Left Plot: Gain vs Argon % (Log Scale for Exponential Rise)
    cTrends->cd(1);
    gPad->SetGrid();
    gPad->SetLogy(); // Log scale is critical for Gain
    TGraphErrors *grGain = new TGraphErrors(n, x_ar, y_gain, nullptr, e_gain);
    grGain->SetTitle("Absolute Gain Trend;Argon Fraction (%);Absolute Gain");
    grGain->SetMarkerStyle(21);
    grGain->SetMarkerColor(kBlue+1);
    grGain->SetLineColor(kBlue+1);
    grGain->SetLineWidth(2);
    grGain->Draw("APC"); // A=Axes, P=Points, C=Smooth Spline Curve

    // Right Plot: IBF vs Argon % (Linear scale for decline)
    cTrends->cd(2);
    gPad->SetGrid();
    TGraph *grIBF = new TGraph(n, x_ar, y_ibf);
    grIBF->SetTitle("IBF Trend;Argon Fraction (%);Ion Backflow Fraction");
    grIBF->SetMarkerStyle(20);
    grIBF->SetMarkerColor(kRed+1);
    grIBF->SetLineColor(kRed+1);
    grIBF->SetLineWidth(2);
    grIBF->Draw("APC");

    std::cout << "Trend curves generated successfully." << std::endl;
    app.Run();
    return 0;

   
}