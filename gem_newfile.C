#include <TApplication.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>

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

    MediumMagboltz gas("ar", 90., "co2", 10.);
    gas.SetTemperature(293.15);
    gas.SetPressure(760.);
    gas.Initialise(true);
    gas.EnablePenningTransfer(0.51, 0., "ar");
    gas.LoadIonMobility("IonMobility_Ar+_Ar.txt");


    ComponentAnsys123 fm;
    fm.Initialise("ELIST.lis", "NLIST.lis", "MPLIST.lis", "PRNSOL.lis", "mm");
    fm.EnableMirrorPeriodicityX();
    fm.EnableMirrorPeriodicityY();
    fm.SetGas(&gas);

    constexpr double pitch = 0.014;

    ViewField fieldView(&fm);
    ViewFEMesh meshView(&fm);

    TCanvas* cf = new TCanvas("cf", "Electric Potential & Geometry", 750, 650);
    fieldView.SetCanvas(cf);
    fieldView.SetPlane(0, -1, 0, 0, 0, 0); // XZ Projection Slice
    fieldView.SetArea(-0.5 * pitch, -0.02, 0.5 * pitch, 0.02);
    fieldView.SetVoltageRange(-160., 160.);
    fieldView.PlotContour();

    meshView.SetCanvas(cf);
    meshView.SetPlane(0, -1, 0, 0, 0, 0);
    meshView.SetArea(-0.5 * pitch, -0.02, 0.5 * pitch, 0.02);
    meshView.SetFillMesh(true);
    meshView.SetColor(2, kGray);
    meshView.Plot(true); // Overlay mesh structure over fields

    
    // SIMULATION CONFIGURATION
    Sensor sensor(&fm);
    sensor.SetArea(-5 * pitch, -5 * pitch, -0.01, 5 * pitch, 5 * pitch, 0.025);
    
    AvalancheMicroscopic aval(&sensor);
    AvalancheMC drift(&sensor);
    drift.SetDistanceSteps(2.e-4);

    ViewDrift driftView;
    aval.EnablePlotting(&driftView, 10); // Display avalanche tracks


    // STATISTICS & RAW METRIC COUNTERS
    std::size_t totalStartingElectrons   = 0;
    std::size_t totalMultipliedElectrons = 0;
    std::size_t totalCollectedElectrons  = 0;
    std::size_t nTotalIons               = 0;
    std::size_t nBFIons                  = 0;

    constexpr std::size_t nEvents = 100;

    
    // MAIN SIMULATION LOOP
    for (std::size_t i = 0; i < nEvents; ++i) {
        const double x0 = -0.5 * pitch + RndmUniform() * pitch;
        const double y0 = -0.5 * pitch + RndmUniform() * pitch;

        // Seeding initial primary electron at z = 0.02 cm
        aval.AvalancheElectron(x0, y0, 0.02, 0., 0.1, 0., 0., 0.);
        totalStartingElectrons += 1; 

        // 1. Get raw amount after avalanche multiplication
        int ne = 0, ni = 0;
        aval.GetAvalancheSize(ne, ni);
        totalMultipliedElectrons += ne;

        // 2. Loop tracks for readout collection and ion drift
        for (const auto& electron : aval.GetElectrons()) {
            // Check if electron reached the bottom readout threshold
            if (electron.path.back().z < -0.009) {
                ++totalCollectedElectrons;
            }

            // Track accompanying positive ion
            const auto& birth = electron.path[0];
            drift.DriftIon(birth.x, birth.y, birth.z, birth.t);
            ++nTotalIons;

            // Check if ion backflow drifted back past upper boundary
            if (!drift.GetIons().empty()) {
                if (drift.GetIons().front().path.back().z > 0.012) {
                    ++nBFIons;
                }
            }
        }
    }


    double ibfFrac = nTotalIons > 0 ? static_cast<double>(nBFIons) / nTotalIons : 0.0;

  
    std::cout << "Total Simulated Events          : " << nEvents << "\n";
    std::cout << "Total Starting Primary e-       : " << totalStartingElectrons << "\n";
    std::cout << "Total Multiplied e- (Absolute)  : " << totalMultipliedElectrons << "\n";
    std::cout << "Total Collected e- (Readout)    : " << totalCollectedElectrons << "\n";
    std::cout << "Total Ions Generated            : " << nTotalIons << "\n";
    std::cout << "Total Ions Backflowed           : " << nBFIons << "\n";
    std::cout << "Calculated IBF Fraction         : " << ibfFrac << "\n";

    TCanvas* cD = new TCanvas("cD", "Avalanche Tracks Mesh Layout", 600, 600);
    ViewFEMesh meshDriftView(&fm);
    
    meshDriftView.SetCanvas(cD);
    meshDriftView.SetPlane(0, -1, 0, 0, 0, 0); 
    meshDriftView.SetArea(-2 * pitch, -0.02, 2 * pitch, 0.02);
    meshDriftView.SetFillMesh(true);
    meshDriftView.SetColor(2, kYellow + 3); 
    meshDriftView.EnableAxes();
    
    meshDriftView.SetViewDrift(&driftView);
    meshDriftView.Plot();

    app.Run();
    return 0;
}
