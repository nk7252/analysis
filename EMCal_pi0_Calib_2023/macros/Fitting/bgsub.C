#include <Math/MinimizerOptions.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TStyle.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "sPhenixStyle.C"
#include "sPhenixStyle.h"

// Combined function for Gaussian + Polynomial
double combinedFunctionDoubleGauss(double* x, double* par)
{
    double gauss1 = par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));
    double gauss2 = par[8] * exp(-0.5 * pow((x[0] - par[9]) / par[10], 2));
    double poly = par[3] + par[4] * x[0] + par[5] * x[0] * x[0] + par[6] * x[0] * x[0] * x[0] + par[7] * x[0] * x[0] * x[0] * x[0];
    return gauss1 + gauss2 + poly;
}

// Combined function for Gaussian + Double Polynomial
double combinedFunctionDoubleGaussDoublePoly(double* x, double* par)
{
    double gauss1 = par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));
    double gauss2 = par[7] * exp(-0.5 * pow((x[0] - par[8]) / par[9], 2));
    double poly1 = par[3] + par[4] * x[0] + par[5] * x[0] * x[0] + par[6] * x[0] * x[0] * x[0];
    double poly2 = par[10] + par[11] * x[0] + par[12] * x[0] * x[0];
    return gauss1 + gauss2 + poly1 + poly2;
}

// Combined function for Gaussian + Logarithm
double combinedFunctionDoubleGaussLog(double* x, double* par)
{
    double gauss1 = par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));
    double gauss2 = par[3] * exp(-0.5 * pow((x[0] - par[4]) / par[5], 2));
    double logBg = par[6] * log(x[0]) + par[7];
    return gauss1 + gauss2 + logBg;
}

// Background fit function with double polynomial
double doublePolyBG(double* x, double* par)
{
    double poly1 = 0;
    if (x[0] >= par[4] && x[0] <= par[5])
    {
        poly1 = par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0];
    }
    double poly2 = 0;
    if (x[0] >= par[9] && x[0] <= par[10])
    {
        poly2 = par[6] + par[7] * x[0] + par[8] * x[0] * x[0];
    }
    return poly1 + poly2;
}

double poly4bg(double* x, double* par)
{
    double gauss_min1 = par[5];
    double gauss_max1 = par[6];
    double gauss_min2 = par[7];
    double gauss_max2 = par[8];
    if ((x[0] >= gauss_min1 && x[0] <= gauss_max1) || (x[0] >= gauss_min2 && x[0] <= gauss_max2))
    {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0] + par[4] * x[0] * x[0] * x[0] * x[0];
}

class FitManager
{
public:
    FitManager(double scale_factor, std::vector<float>& limits, int startBin, int endBin, int projectionBins, int rebinFactor, bool dynamic_left, int background_scheme)
        : scale_factor(scale_factor)
        , limits(limits)
        , startBin(startBin)
        , endBin(endBin)
        , projectionBins(projectionBins)
        , rebinFactor(rebinFactor)
        , dynamic_left(dynamic_left)
        , background_scheme(background_scheme)
    {
        initialize();
    }

    void fitHistogram()
    {
        SetsPhenixStyle();
        std::vector<double> pionPt, pionPeak, pionRes, pionPtErr, pionPeakErr, pionResErr;
        std::vector<double> etaPeak, etaRes, etaPeakErr, etaResErr;
        std::vector<double> PeakRatio, PeakRatioErr;

        for (int i = startBin; i <= endBin; i += projectionBins)
        {
            int lastBin = std::min(i + projectionBins - 1, nBinsX);
            TH1D* hist = hist2D->ProjectionY(Form("proj_%d", i), i, lastBin);
            if (hist->GetEntries() < 1000)
            {
                delete hist;
                continue;
            }

            if (rebinFactor > 1)
            {
                hist->Rebin(rebinFactor);
            }

            hist->Scale(1 / (hist->GetBinLowEdge(2) - hist->GetBinLowEdge(1)));

            if (dynamic_left)
            {
                float leftmost_limit = 0;
                for (int bin = 1; bin <= hist->GetNbinsX(); ++bin)
                {
                    if (hist->GetBinContent(bin) > 0)
                    {
                        leftmost_limit = hist->GetBinLowEdge(bin);
                        limits[0] = leftmost_limit;
                        break;
                    }
                }
            }

            double pt_min = hist2D->GetXaxis()->GetBinLowEdge(i);
            double pt_max = hist2D->GetXaxis()->GetBinUpEdge(lastBin);
            TString ptRange = Form("pt_%.2f-%.2f_GeV", pt_min, pt_max);

            TF1* Backgroundonly;
            if (background_scheme == 0)
            {
                Backgroundonly = new TF1("Backgroundonly", poly4bg, limits[0], limits[1], 9);
            }
            else if (background_scheme == 1)
            {
                Backgroundonly = new TF1("Backgroundonly", doublePolyBG, limits[0], limits[1], 11);
            }
            else if (background_scheme == 2)
            {
                Backgroundonly = new TF1("Backgroundonly", "log(x)", limits[0], limits[1]);
            }

            hist->Fit(Backgroundonly, "RME");

            TF1* combinedFit = createCombinedFunction(hist, Backgroundonly);

            hist->Fit(combinedFit, "RME");

            saveFitResults(hist, combinedFit, pionPt, pionPtErr, pionPeak, pionPeakErr, pionRes, pionResErr, etaPeak, etaPeakErr, etaRes, etaResErr, PeakRatio, PeakRatioErr, ptRange);
            saveGraphs(pionPt, pionPtErr, pionPeak, pionPeakErr, pionRes, pionResErr, etaPeak, etaPeakErr, etaRes, etaResErr, PeakRatio, PeakRatioErr);

            delete hist;
            delete Backgroundonly;
            delete combinedFit;
        }
    }

private:
    double scale_factor;
    std::vector<float>& limits;
    int startBin;
    int endBin;
    int projectionBins;
    int rebinFactor;
    bool dynamic_left;
    int background_scheme;
    TFile* file;
    TH2F* hist2D;
    int nBinsX;

    void initialize()
    {
        ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
        file = new TFile("/sphenix/u/nkumar/analysis/EMCal_pi0_Calib_2023/macros/condor/output/merged_file.root");
        hist2D = (TH2F*) file->Get("h_InvMass_2d");
        nBinsX = hist2D->GetNbinsX();
        if (endBin == -1) endBin = nBinsX;
    }

    TF1* createCombinedFunction(TH1D* hist, TF1* Backgroundonly)
    {
        TF1* combinedFit;
        if (background_scheme == 0)
        {
            combinedFit = new TF1("combinedFit", combinedFunctionDoubleGauss, limits[0], limits[1], 11);
            setInitialParameters(combinedFit, hist, Backgroundonly, 3, 8);
        }
        else if (background_scheme == 1)
        {
            combinedFit = new TF1("combinedFit", combinedFunctionDoubleGaussDoublePoly, limits[0], limits[1], 13);
            setInitialParameters(combinedFit, hist, Backgroundonly, 3, 7);
            setInitialParameters(combinedFit, hist, Backgroundonly, 10, 3);
        }
        else if (background_scheme == 2)
        {
            combinedFit = new TF1("combinedFit", combinedFunctionDoubleGaussLog, limits[0], limits[1], 8);
            setInitialParameters(combinedFit, hist, Backgroundonly, 3, 2);
        }
        return combinedFit;
    }

    void setInitialParameters(TF1* combinedFit, TH1D* hist, TF1* Backgroundonly, int startIdx, int numParams)
    {
        for (int i = 0; i < 3; ++i) combinedFit->SetParameter(i, hist->GetFunction("gausFit")->GetParameter(i));
        for (int i = startIdx; i < startIdx + numParams; ++i) combinedFit->SetParameter(i, Backgroundonly->GetParameter(i - startIdx));
    }

    void saveFitResults(TH1D* hist, TF1* combinedFit, std::vector<double>& pionPt, std::vector<double>& pionPtErr, std::vector<double>& pionPeak, std::vector<double>& pionPeakErr, std::vector<double>& pionRes, std::vector<double>& pionResErr, std::vector<double>& etaPeak, std::vector<double>& etaPeakErr, std::vector<double>& etaRes, std::vector<double>& etaResErr, std::vector<double>& PeakRatio, std::vector<double>& PeakRatioErr, TString ptRange)
    {
        double pion_pt = (hist->GetXaxis()->GetXmin() + hist->GetXaxis()->GetXmax()) / 2.0;
        double pion_peak = combinedFit->GetParameter(1);
        double pion_peak_err = combinedFit->GetParError(1);
        double pion_res = combinedFit->GetParameter(2) / combinedFit->GetParameter(1);
        double pion_res_err = pion_res * sqrt(pow(combinedFit->GetParError(2) / combinedFit->GetParameter(2), 2) + pow(pion_peak_err / pion_peak, 2));

        pionPt.push_back(pion_pt);
        pionPtErr.push_back(0); // Adjust as necessary
        pionPeak.push_back(pion_peak);
        pionPeakErr.push_back(pion_peak_err);
        pionRes.push_back(pion_res);
        pionResErr.push_back(pion_res_err);

        double eta_peak = combinedFit->GetParameter(9);
        double eta_peak_err = combinedFit->GetParError(9);
        double eta_res = combinedFit->GetParameter(10) / combinedFit->GetParameter(9);
        double eta_res_err = eta_res * sqrt(pow(combinedFit->GetParError(10) / combinedFit->GetParameter(10), 2) + pow(eta_peak_err / eta_peak, 2));

        etaPeak.push_back(eta_peak);
        etaPeakErr.push_back(eta_peak_err);
        etaRes.push_back(eta_res);
        etaResErr.push_back(eta_res_err);
        PeakRatio.push_back(pion_peak / eta_peak);
        PeakRatioErr.push_back(sqrt(pow(eta_peak_err / eta_peak, 2) + pow(pion_peak_err / pion_peak, 2)));
    }

    void saveGraphs(std::vector<double>& pionPt, std::vector<double>& pionPtErr, std::vector<double>& pionPeak, std::vector<double>& pionPeakErr, std::vector<double>& pionRes, std::vector<double>& pionResErr, std::vector<double>& etaPeak, std::vector<double>& etaPeakErr, std::vector<double>& etaRes, std::vector<double>& etaResErr, std::vector<double>& PeakRatio, std::vector<double>& PeakRatioErr)
    {
        int nPoints = pionPt.size();
        TGraphErrors* gPionPeak = new TGraphErrors(nPoints, &pionPt[0], &pionPeak[0], &pionPtErr[0], &pionPeakErr[0]);
        TGraphErrors* gPionRes = new TGraphErrors(nPoints, &pionPt[0], &pionRes[0], &pionPtErr[0], &pionResErr[0]);
        TGraphErrors* gEtaPeak = new TGraphErrors(nPoints, &pionPt[0], &etaPeak[0], &pionPtErr[0], &etaPeakErr[0]);
        TGraphErrors* gEtaRes = new TGraphErrors(nPoints, &pionPt[0], &etaRes[0], &pionPtErr[0], &etaResErr[0]);
        TGraphErrors* gPeakRatio = new TGraphErrors(nPoints, &pionPt[0], &PeakRatio[0], &pionPtErr[0], &PeakRatioErr[0]);

        TCanvas* cPionPeak = new TCanvas("cPionPeak", "Pion Peak Position", 800, 600);
        gPionPeak->Draw("ALP");
        cPionPeak->Print("2D_Histogram_Fits.pdf");

        TCanvas* cPionRes = new TCanvas("cPionRes", "Pion Relative Resolution", 800, 600);
        gPionRes->Draw("ALP");
        cPionRes->Print("2D_Histogram_Fits.pdf");

        TCanvas* cEtaPeak = new TCanvas("cEtaPeak", "Eta Peak Position", 800, 600);
        gEtaPeak->Draw("ALP");
        cEtaPeak->Print("2D_Histogram_Fits.pdf");

        TCanvas* cEtaRes = new TCanvas("cEtaRes", "Eta Relative Resolution", 800, 600);
        gEtaRes->Draw("ALP");
        cEtaRes->Print("2D_Histogram_Fits.pdf");

        TCanvas* cPeakRatio = new TCanvas("cPeakRatio", "Pion/Eta Mass Ratio", 800, 600);
        gPeakRatio->Draw("ALP");
        cPeakRatio->Print("2D_Histogram_Fits.pdf");

        delete cPionPeak;
        delete cPionRes;
        delete gPionPeak;
        delete gPionRes;
        delete cEtaPeak;
        delete cEtaRes;
        delete gEtaPeak;
        delete gEtaRes;
        delete cPeakRatio;
        delete gPeakRatio;
    }
};

int main()
{
    double scale_factor = 1.0;
    std::vector<float> limits = {0.05, 1.0, 0.11, 0.19, 0.11, 0.19, 0.5, 0.7};  // Example limits
    int startBin = 1;
    int endBin = -1;
    int projectionBins = 1;
    int rebinFactor = 1;
    bool dynamic_left = true;
    int background_scheme = 2;  // Set to 2 for logarithm background

    FitManager fitManager(scale_factor, limits, startBin, endBin, projectionBins, rebinFactor, dynamic_left, background_scheme);
    fitManager.fitHistogram();

    return 0;
}