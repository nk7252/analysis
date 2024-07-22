#pragma once
// root includes
#include <Math/MinimizerOptions.h>
#include <Math/Minimizer.h>
#include <Math/Functor.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TString.h>
#include <TStyle.h>
// general includes
#include <fstream>
#include <iostream>
#include <vector>
// sphenix includes
#include "sPhenixStyle.C"
#include "sPhenixStyle.h"

//global bin var
std::vector<double> nuBins = {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.36, 0.40, 0.44, 0.48, 0.50, 0.52, 0.54, 0.56, 0.58, 0.60, 0.62, 0.64, 0.66, 0.68, 0.70, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 1.0, 1.04, 1.08, 1.12, 1.16, 1.2};
// Combined function for Gaussian + Polynomial
double combinedFunction(double *x, double *par)
{
  // Gaussian part
  double gauss = par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));

  // Polynomial part
  double poly = par[3] + par[4] * x[0] + par[5] * x[0] * x[0] + par[6] * x[0] * x[0] * x[0] + par[7] * x[0] * x[0] * x[0] * x[0];

  return gauss + poly;
}

double combinedFunctionDoubleGauss(double *x, double *par)
{
  // First Gaussian part (e.g., pion peak)
  double gauss1 = par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));

  // Second Gaussian part (e.g., eta peak)
  double gauss2 = par[8] * exp(-0.5 * pow((x[0] - par[9]) / par[10], 2));

  // Polynomial part
  double poly = par[3] + par[4] * x[0] + par[5] * x[0] * x[0] + par[6] * x[0] * x[0] * x[0] + par[7] * x[0] * x[0] * x[0] * x[0];

  return gauss1 + gauss2 + poly;
}

double combinedFunctionDoubleGaussDoublePoly(double *x, double *par)
{
  // First Gaussian part (e.g., pion peak)
  double gauss1 = par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));

  // Second Gaussian part (e.g., eta peak)
  double gauss2 = par[3] * exp(-0.5 * pow((x[0] - par[4]) / par[5], 2));

  // Polynomial part
  double poly1 = 0;
  double poly2 = 0;
  double boundary = 0.35;
  if (x[0] <= boundary)
  {
    poly1 = par[6] + par[7] * x[0] + par[8] * x[0] * x[0] + par[9] * x[0] * x[0] * x[0];
  }
  else
  {
    // Calculate poly1 at the boundary
    double poly1_boundary = par[6] + par[7] * boundary + par[8] * boundary * boundary + par[9] * boundary * boundary * boundary;
    // Calculate the derivative of poly1 at the boundary
    double poly1_derivative = par[7] + 2 * par[8] * boundary + 3 * par[9] * boundary * boundary;

    // Set the coefficients of poly2 to match the value and derivative at the boundary
    poly2 = poly1_boundary + poly1_derivative * (x[0] - boundary) + par[10] * (x[0] - boundary) * (x[0] - boundary);
  }

  return gauss1 + gauss2 + poly1 + poly2;
}

double combinedFunctionDoubleGaussLog(double *x, double *par)
{
  // First Gaussian part (e.g., pion peak)
  double gauss1 = 0;
  if (x[0] >= 0.09 && x[0] <= 0.21)
  {  // Check if x is in the range of the first Gaussian
    double gauss1 = par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));
  }

  // Second Gaussian part (e.g., eta peak)
  double gauss2 = 0;
  if (x[0] >= 0.5 && x[0] <= 0.7)
  {  // Check if x is in the range of the first Gaussian
    double gauss2 = par[3] * exp(-0.5 * pow((x[0] - par[4]) / par[5], 2));
  }
  // Logarithm background part
  double logBg = par[6] * log(x[0]) + par[7];

  return gauss1 + gauss2 + logBg;
}

double combinedFunctionDoubleGaussPoly2(double *x, double *par)
{
  // First Gaussian part (e.g., pion peak)
  double gauss1 = par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));

  // Second Gaussian part (e.g., eta peak)
  double gauss2 = par[3] * exp(-0.5 * pow((x[0] - par[4]) / par[5], 2));

  // Polynomial part (2nd degree)
  double poly = par[6] + par[7] * x[0] + par[8] * x[0] * x[0];

  return gauss1 + gauss2 + poly;
}

double combinedFunctionDoubleGaussPoly3(double *x, double *par)
{
  // First Gaussian part (e.g., pion peak)
  double gauss1 = par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));

  // Second Gaussian part (e.g., eta peak)
  double gauss2 = par[3] * exp(-0.5 * pow((x[0] - par[4]) / par[5], 2));

  // Polynomial part (3rd degree)
  double poly = par[6] + par[7] * x[0] + par[8] * x[0] * x[0] + par[9] * x[0] * x[0] * x[0];

  return gauss1 + gauss2 + poly;
}

double combinedFunctionDoubleGaussPoly5(double *x, double *par)
{
  // First Gaussian part (e.g., pion peak)
  double gauss1 = par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));

  // Second Gaussian part (e.g., eta peak)
  double gauss2 = par[3] * exp(-0.5 * pow((x[0] - par[4]) / par[5], 2));

  // Polynomial part (5th degree)
  double poly = par[6] + par[7] * x[0] + par[8] * x[0] * x[0] + par[9] * x[0] * x[0] * x[0] + par[10] * x[0] * x[0] * x[0] * x[0] + par[11] * x[0] * x[0] * x[0] * x[0] * x[0];

  return gauss1 + gauss2 + poly;
}

double combinedFunctionDoubleGaussPoly6(double *x, double *par)
{
  // First Gaussian part (e.g., pion peak)
  double gauss1 = par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));

  // Second Gaussian part (e.g., eta peak)
  double gauss2 = par[3] * exp(-0.5 * pow((x[0] - par[4]) / par[5], 2));

  // Polynomial part (5th degree)
  double poly = par[6] + par[7] * x[0] + par[8] * x[0] * x[0] + par[9] * x[0] * x[0] * x[0] + par[10] * x[0] * x[0] * x[0] * x[0] + par[11] * x[0] * x[0] * x[0] * x[0] * x[0] + par[12] * x[0] * x[0] * x[0] * x[0] * x[0] * x[0];

  return gauss1 + gauss2 + poly;
}

double doubleGauss(double *x, double *par)
{
  // First Gaussian part (e.g., pion peak)
  double gauss1 = 0;
  double gauss2 = 0;
  if (x[0] >= par[6] && x[0] <= par[7])
  {  // Check if x is in the range of the first Gaussian
    gauss1 = par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));
  }
  else if(x[0] >= par[8] && x[0] <= par[9])
  {  // Check if x is in the range of the second Gaussian
    gauss2 = par[3] * exp(-0.5 * pow((x[0] - par[4]) / par[5], 2));
  }
  else 
  {
    TF1::RejectPoint();
    return 0;
  } 
  return gauss1 + gauss2;
}

double doublePolyBG(double *x, double *par)
{
  double poly1 = 0;
  double poly2 = 0;
  double boundary = 0.35;

  if (x[0] <= boundary)
  {
    poly1 = par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0];
  }
  else
  {
    // Calculate poly1 at the boundary
    double poly1_boundary = par[0] + par[1] * boundary + par[2] * boundary * boundary + par[3] * boundary * boundary * boundary;
    // Calculate the derivative of poly1 at the boundary
    double poly1_derivative = par[1] + 2 * par[2] * boundary + 3 * par[3] * boundary * boundary;

    // Set the coefficients of poly2 to match the value and derivative at the boundary
    poly2 = poly1_boundary + poly1_derivative * (x[0] - boundary) + par[4] * (x[0] - boundary) * (x[0] - boundary);
  }

  // Check if x is in the range of any Gaussian fit
  if ((x[0] >= 0.1 && x[0] <= 0.2) || (x[0] >= 0.5 && x[0] <= 0.7))
  {
    TF1::RejectPoint();
    return 0;
  }

  return poly1 + poly2;
}

double ONLYdoublePolyBG(double *x, double *par)
{
  // double  poly1 = par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0];
  // double  poly2 = par[4] + par[5] * x[0] + par[6] * x[0] * x[0];
  double poly1 = 0;
  double boundary = 0.35;
  if (x[0] >= 0.05 && x[0] <= 0.35)
  {  // Check if x is in the range of the first Gaussian
    poly1 = par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0];
  }
  else
  {
    // Calculate poly1 at the boundary
    double poly1_boundary = par[0] + par[1] * boundary + par[2] * boundary * boundary + par[3] * boundary * boundary * boundary;
    // Calculate the derivative of poly1 at the boundary
    double poly1_derivative = par[1] + 2 * par[2] * boundary + 3 * par[3] * boundary * boundary;

    poly1 = poly1_boundary + poly1_derivative * (x[0] - boundary) + par[4] * (x[0] - boundary) * (x[0] - boundary);  // par[4] + par[5] * x[0] + par[6] * x[0] * x[0];
  }
  return poly1;  // + poly2;
}

double LogBG(double *x, double *par)
{
  // First Gaussian part (e.g., pion peak)
  // double logBg = 0;

  double logBg = par[0] * log(x[0]) + par[1];

  // Check if x is in the range of any Gaussian fit
  if ((x[0] >= 0.1 && x[0] <= 0.2) || (x[0] >= 0.52 && x[0] <= 0.68))
  {
    TF1::RejectPoint();
    return 0;
  }

  return logBg;  // + poly2;
}

double poly2BG(double *x, double *par)
{
  // 2nd degree polynomial background
  // Check if x is in the range of any Gaussian fit
  if ((x[0] >= 0.1 && x[0] <= 0.2) || (x[0] >= 0.52 && x[0] <= 0.68))
  {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0];
}

double poly3BG(double *x, double *par)
{
  // 3rd degree polynomial background
  // Check if x is in the range of any Gaussian fit
  if ((x[0] >= 0.1 && x[0] <= 0.2) || (x[0] >= 0.52 && x[0] <= 0.68))
  {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0];
}

double poly5BG(double *x, double *par)
{
  // 5th degree polynomial background
  // Check if x is in the range of any Gaussian fit
  if ((x[0] >= 0.1 && x[0] <= 0.2) || (x[0] >= 0.52 && x[0] <= 0.68))
  {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0] + par[4] * x[0] * x[0] * x[0] * x[0] + par[5] * x[0] * x[0] * x[0] * x[0] * x[0];
}

double poly6BG(double *x, double *par)
{
  // 6th degree polynomial background
  // Check if x is in the range of any Gaussian fit
  if ((x[0] >= 0.1 && x[0] <= 0.2) || (x[0] >= 0.52 && x[0] <= 0.68))
  {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0] + par[4] * x[0] * x[0] * x[0] * x[0] + par[5] * x[0] * x[0] * x[0] * x[0] * x[0] + par[6] * x[0] * x[0] * x[0] * x[0] * x[0] * x[0];
}

// leftRightPolynomial function to optionally exclude two Gaussian regions
double leftRightPolynomial(double *x, double *par)
{
  // Define the range for the first Gaussian fit (pion peak)
  double gauss_min1 = par[5];
  double gauss_max1 = par[6];

  // Define the range for the second Gaussian fit (eta peak)
  double gauss_min2 = par[7];
  double gauss_max2 = par[8];

  // Check if x is in the range of any Gaussian fit
  if ((x[0] >= gauss_min1 && x[0] <= gauss_max1) || (x[0] >= gauss_min2 && x[0] <= gauss_max2))
  {
    TF1::RejectPoint();
    return 0;
  }

  // Polynomial (4th degree) calculation
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0] + par[4] * x[0] * x[0] * x[0] * x[0];
}

TH1D *rebinHistogram(TH1D *h, const std::vector<double> &binEdges)
{
  const int Nbins = binEdges.size() - 1;
  double bins[Nbins + 1];
  for (int i = 0; i < Nbins + 1; ++i)
  {
    bins[i] = binEdges[i];
  }
  return (TH1D *) h->Rebin(Nbins, "hrb", bins);
}

// scale the histogram's error bars
void scale_histogram_errors(TH1D *hist_error_scale, double scale_factor)
{
  for (int i = 1; i <= hist_error_scale->GetNbinsX(); ++i)
  {
    // Get the current error
    double current_error = hist_error_scale->GetBinError(i);

    // Scale the error by the scale factor
    hist_error_scale->SetBinError(i, current_error * scale_factor);
    // std::cout << "orig bin cont: " << hist_error_scale->GetBinContent(i) << " . bin error: " << current_error << " . New bin error: " << hist_error_scale->GetBinError(i) <<std::endl;
  }
}

void appendtextfile(TF1 *fitFunc, const std::string &fitName, double scale_factor)
{
  // Open a text file in append mode
  std::ofstream outfile;
  outfile.open("/sphenix/u/nkumar/analysis/EMCal_pi0_Calib_2023/macros/Fitting/fit_parameters.txt", std::ios_base::app);

  // Check if the file is open (and thus valid)
  if (outfile.is_open())
  {
    // Check if the file is empty
    outfile.seekp(0, std::ios::end);
    size_t size = outfile.tellp();

    if (size == 0)
    {
      // File is new, write header
      outfile << "Fit Name          Scale Factor          Mean          Sigma          Sigma/Mean          chi2          NDF           chi2/NDF" << std::endl;
    }

    // Write the parameters to the file
    outfile << fitName << "     " << scale_factor << "     " << fitFunc->GetParameter(1) << "     " << fitFunc->GetParameter(2) << "     " << fitFunc->GetParameter(2) / fitFunc->GetParameter(1) << "     " << fitFunc->GetChisquare() << "     " << fitFunc->GetNDF() << "     " << fitFunc->GetChisquare() / fitFunc->GetNDF() << std::endl;

    outfile.close();
  }
  else
  {
    std::cerr << "Error: Unable to open file." << std::endl;
  }
}

void fit_2d_histogram(Double_t scale_factor, std::vector<float> &limits, bool fitEtaPeak = false, int startBin = 1, int endBin = -1, int projectionBins = 1, int rebinFactor = 1, bool dynamic_left = false, int background_scheme = 0, bool var_bins = false)//std::vector<double> &rebinEdges
{
  // more thorough minimizer for fit

  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Fumili2");  //,"Simplex", "Migrad", "Combined"
  // ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
  // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("GSLMultiMin");//, "ConjugateFR"
  //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("GSLMultiFit", "LevenbergMarquardt");
  // ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(3);
  //   Set the global fit strategy
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
  ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(10000);
  //ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.0001);
  //ROOT::Math::MinimizerOptions::SetDefaultPrecision(1e-8);
  SetsPhenixStyle();

  // Ensure the limits vector has the correct size
  if (limits.size() < 6 || (fitEtaPeak && limits.size() < 10))
  {
    std::cerr << "Insufficient limits provided. Expected 6 (or 10 if fitting eta peak) values." << std::endl;
    return;
  }

  // Open the ROOT file and get the 2D histogram
  // TFile *file = new TFile("/sphenix/u/nkumar/analysis/EMCal_pi0_Calib_2023/macros/condor/output/merged_file.root");
  TFile *file = new TFile("/sphenix/u/nkumar/analysis/EMCal_pi0_Calib_2023/macros/condor/output/OUTHIST_iter_DST_CALO_CLUSTER_pythia8_pp_mb_3MHz_0000000011_merged.root");
  TH2F *hist2D = (TH2F *) file->Get("h_InvMass_2d");

  int nBinsX = hist2D->GetNbinsX();
  if (endBin == -1) endBin = nBinsX;  // Default to the last bin if not specified

  // Create a PDF to save the canvases
  TCanvas *dummyCanvas = new TCanvas();
  dummyCanvas->Print("2D_Histogram_Fits.pdf[");

  // Vectors to store fit results
  std::vector<double> pionPt, pionPeak, pionRes, pionPtErr, pionPeakErr, pionResErr;
  std::vector<double> etaPeak, etaRes, etaPeakErr, etaResErr;
  std::vector<double> PeakRatio, PeakRatioErr;

  for (int i = startBin; i <= endBin; i += projectionBins)
  {
    // Project the histogram along the Y-axis
    int lastBin = std::min(i + projectionBins - 1, nBinsX);
    TH1D *hist = hist2D->ProjectionY(Form("proj_%d", i), i, lastBin);

    // Check if the projection has enough entries to perform a fit
    if (hist->GetEntries() < 1000)
    {  // Adjust the threshold as needed
      delete hist;
      continue;
    }

    // Rebin the projected histogram if needed
    TH1D *histF = (TH1D *) hist;
    if (var_bins && !nuBins.empty())
    {
      std::cout << "Rebinning histogram with non-uniform edges" << std::endl;
      histF = rebinHistogram(histF, nuBins);//nuBins
    }
    else if (rebinFactor > 1)
    {
      histF->Rebin(rebinFactor);
    }

    // Rebin the projected histogram if needed
    // if (rebinFactor > 1)
    //{
    // hist->Rebin(rebinFactor);
    //}

    // Normalize the histogram
    histF->Scale(1. / 2, "width");
    // histF->Scale(1 / (histF->GetBinLowEdge(2) - histF->GetBinLowEdge(1)));  // divide bin content by width so it becomes dN/dM

    // Determine the leftmost point with a value in the projection histograms
    if (dynamic_left)
    {
      float leftmost_limit = 0;
      for (int bin = 1; bin <= histF->GetNbinsX(); ++bin)
      {
        if (histF->GetBinContent(bin) > 0)
        {
          leftmost_limit = histF->GetBinLowEdge(bin);
          limits[0] = leftmost_limit;
          break;
        }
      }
    }

    // Get the pt range for the current slice
    double pt_min = hist2D->GetXaxis()->GetBinLowEdge(i);
    double pt_max = hist2D->GetXaxis()->GetBinUpEdge(lastBin);
    TString ptRange = Form("pt_%.2f-%.2f_GeV", pt_min, pt_max);

    // Set histogram range and scale errors
    // histF->GetXaxis()->SetRangeUser(0, 1.0);
    scale_histogram_errors(histF, scale_factor);

    // Fit left and right regions with a polynomial, excluding Gaussian regions
    TF1 *leftRightFit;
    if (fitEtaPeak)
    {
      if (background_scheme == 0)  // poly4
      {
        leftRightFit = new TF1("leftRightFit", leftRightPolynomial, limits[0], limits[1], 9);
        leftRightFit->SetParameter(5, limits[4]);
        leftRightFit->SetParameter(6, limits[5]);
        leftRightFit->SetParameter(7, limits[8]);
        leftRightFit->SetParameter(8, limits[9]);
      }
      else if (background_scheme == 1)  // poly3+poly2
      {
        leftRightFit = new TF1("leftRightFit", doublePolyBG, limits[0], limits[1], 5);
      }
      else if (background_scheme == 2)
      {
        leftRightFit = new TF1("leftRightFit", LogBG, limits[0], limits[1], 2);
      }
      else if (background_scheme == 3)
      {
        leftRightFit = new TF1("leftRightFit", poly2BG, limits[0], limits[1], 3);
      }
      else if (background_scheme == 4)
      {
        leftRightFit = new TF1("leftRightFit", poly3BG, limits[0], limits[1], 4);
      }
      else if (background_scheme == 5)
      {
        leftRightFit = new TF1("leftRightFit", poly5BG, limits[0], limits[1], 6);
      }
      else if (background_scheme == 6)
      {
        leftRightFit = new TF1("leftRightFit", poly6BG, limits[0], limits[1], 7);
      }
    }
    else
    {
      leftRightFit = new TF1("leftRightFit", leftRightPolynomial, limits[0], limits[1], 7);
      leftRightFit->SetParameter(5, limits[4]);
      leftRightFit->SetParameter(6, limits[5]);
    }
    histF->Fit(leftRightFit, "RE");

    // Fit first Gaussian in the specified range
    TF1 *gausFit = new TF1("gausFit", "gaus", limits[2], limits[3]);
    gausFit->SetParLimits(1, 0.11, 0.19);
    gausFit->SetParLimits(2, 0.05, 0.25);
    histF->Fit(gausFit, "RE");

    // Combined Gaussian + Polynomial fit
    TF1 *combinedFit;
    if (fitEtaPeak)
    {
      if (background_scheme == 0)  // poly4
      {
        combinedFit = new TF1("combinedFit", combinedFunctionDoubleGauss, limits[0], limits[1], 11);  // 2 Gaussians + 1 polynomial = 3 + 3 + 5
      }
      else if (background_scheme == 1)  // poly3+poly2
      {
        combinedFit = new TF1("combinedFit", combinedFunctionDoubleGaussDoublePoly, limits[0], limits[1], 11);  // 2 Gaussians + 1 poly3 +1poly2 = 3 + 3 + 4 + 3=13; mathing at boundaries removes 2 free parameters
      }
      else if (background_scheme == 2)
      {
        combinedFit = new TF1("combinedFit", combinedFunctionDoubleGaussLog, limits[0], limits[1], 8);
      }
      else if (background_scheme == 3)
      {
        combinedFit = new TF1("combinedFit", combinedFunctionDoubleGaussPoly2, limits[0], limits[1], 9);
      }
      else if (background_scheme == 4)
      {
        combinedFit = new TF1("combinedFit", combinedFunctionDoubleGaussPoly3, limits[0], limits[1], 10);
      }
      else if (background_scheme == 5)
      {
        combinedFit = new TF1("combinedFit", combinedFunctionDoubleGaussPoly5, limits[0], limits[1], 12);
      }
      else if (background_scheme == 6)
      {
        combinedFit = new TF1("combinedFit", combinedFunctionDoubleGaussPoly6, limits[0], limits[1], 13);
      }
    }
    else
    {
      combinedFit = new TF1("combinedFit", combinedFunction, limits[0], limits[1], 8);
    }

    // Set initial parameters from previous fits
    for (int j = 0; j < 3; ++j)
    {
      combinedFit->SetParameter(j, gausFit->GetParameter(j));
    }
    combinedFit->SetParLimits(0, 0, gausFit->GetParameter(0) * 1.25);  // gausFit->GetParameter(0) *0.5
    combinedFit->SetParLimits(1, 0.11, 0.19);
    combinedFit->SetParLimits(2, 0.01, 0.25);
    // for (int j = 3; j < 8; ++j) combinedFit->SetParameter(j, leftRightFit->GetParameter(j - 3));

    // Fit second Gaussian in the specified range
    TF1 *gausFit2 = new TF1("gausFit2", "gaus", limits[6], limits[7]);
    gausFit2->SetParLimits(1, 0.55, 0.63);
    gausFit2->SetParLimits(2, 0.03, 0.25);
    histF->Fit(gausFit2, "RE");
    if (fitEtaPeak)
    {
      if (background_scheme == 0)  // poly4
      {
        for (int j = 0; j < 3; ++j) combinedFit->SetParameter(j + 8, gausFit2->GetParameter(j));
        for (int j = 0; j < 5; ++j) combinedFit->SetParameter(j + 3, leftRightFit->GetParameter(j));
        // combinedFit->SetParLimits(8, 10, gausFit->GetParameter(0) / 6);
        combinedFit->SetParLimits(9, 0.55, 0.63);
        combinedFit->SetParLimits(10, 0.05, 0.25);
      }
      else if (background_scheme == 1)  // poly3+poly2
      {
        /*
        for (int j = 0; j < 4; ++j)  // poly3
        {
          combinedFit->SetParameter(j + 6, leftRightFit->GetParameter(j));
        }
        for (int j = 0; j < 3; ++j)  // poly2
        {
          combinedFit->SetParameter(j + 10, leftRightFit->GetParameter(j + 4));
        }*/
        for (int j = 0; j < 4; ++j) combinedFit->SetParameter(j + 6, leftRightFit->GetParameter(j));
        combinedFit->SetParameter(10, leftRightFit->GetParameter(4));

        // combinedFit->SetParLimits(6, leftRightFit->GetParameter(0) * 0.95, leftRightFit->GetParameter(0) * 1.05);
        // combinedFit->SetParLimits(7, leftRightFit->GetParameter(1) * 0.99, leftRightFit->GetParameter(1) * 1.0);
        // combinedFit->SetParLimits(8, leftRightFit->GetParameter(2) * 0.95, leftRightFit->GetParameter(2) * 1.05);
        // combinedFit->SetParLimits(9, leftRightFit->GetParameter(3) * 0.99, leftRightFit->GetParameter(3) * 1.0);
        // combinedFit->SetParLimits(10, leftRightFit->GetParameter(4) * 0.95, leftRightFit->GetParameter(4) * 1.05);
        // combinedFit->SetParLimits(11, leftRightFit->GetParameter(5) * 0.95, leftRightFit->GetParameter(5) * 1.05);
        // combinedFit->SetParLimits(12, leftRightFit->GetParameter(6) * 0.95, leftRightFit->GetParameter(6) * 1.05);
        for (int j = 0; j < 3; ++j)
        {
          combinedFit->SetParameter(j + 3, gausFit2->GetParameter(j));
        }
        combinedFit->SetParLimits(3, 0, gausFit2->GetParameter(0) * 1.05);  // gausFit2->GetParameter(0) *0.5
        combinedFit->SetParLimits(4, 0.5, 0.63);
        combinedFit->SetParLimits(5, 0.03, 0.25);
      }
      else if (background_scheme == 2)
      {
        for (int j = 0; j < 2; ++j) combinedFit->SetParameter(j + 6, leftRightFit->GetParameter(j));
        for (int j = 0; j < 3; ++j) combinedFit->SetParameter(j + 3, gausFit2->GetParameter(j));
        combinedFit->SetParLimits(4, 0.55, 0.63);
        combinedFit->SetParLimits(5, 0.03, 0.25);
      }
      else if (background_scheme == 3)
      {
        for (int j = 0; j < 3; ++j) combinedFit->SetParameter(j + 6, leftRightFit->GetParameter(j));
        for (int j = 0; j < 3; ++j) combinedFit->SetParameter(j + 3, gausFit2->GetParameter(j));
        combinedFit->SetParLimits(4, 0.55, 0.63);
        combinedFit->SetParLimits(5, 0.03, 0.25);
      }
      else if (background_scheme == 4)
      {
        for (int j = 0; j < 4; ++j) combinedFit->SetParameter(j + 6, leftRightFit->GetParameter(j));
        for (int j = 0; j < 3; ++j) combinedFit->SetParameter(j + 3, gausFit2->GetParameter(j));
        combinedFit->SetParLimits(4, 0.55, 0.63);
        combinedFit->SetParLimits(5, 0.03, 0.25);
      }
      else if (background_scheme == 5)
      {
        for (int j = 0; j < 6; ++j) combinedFit->SetParameter(j + 6, leftRightFit->GetParameter(j));
        for (int j = 0; j < 3; ++j) combinedFit->SetParameter(j + 3, gausFit2->GetParameter(j));
        combinedFit->SetParLimits(4, 0.55, 0.63);
        combinedFit->SetParLimits(5, 0.03, 0.25);
        combinedFit->SetParLimits(3, 0, gausFit2->GetParameter(0) * 1.05);  // gausFit2->GetParameter(0) *0.5
      }
      else if (background_scheme == 6)
      {
        for (int j = 0; j < 7; ++j) combinedFit->SetParameter(j + 6, leftRightFit->GetParameter(j));
        for (int j = 0; j < 3; ++j) combinedFit->SetParameter(j + 3, gausFit2->GetParameter(j));
        combinedFit->SetParLimits(4, 0.55, 0.63);
        combinedFit->SetParLimits(5, 0.03, 0.25);
        combinedFit->SetParLimits(3, 0, gausFit2->GetParameter(0) * 1.05);  // gausFit2->GetParameter(0) *0.5
      }
    }
    else
    {
      for (int j = 3; j < 8; ++j) combinedFit->SetParameter(j, leftRightFit->GetParameter(j - 3));
    }

    // Fit the combined function
    combinedFit->SetNpx(1000);
    histF->Fit(combinedFit, "REI");
    
    // After fitting
    std::cout << "Background only Fit Parameters:" << std::endl;
    for (int i = 0; i < leftRightFit->GetNpar(); ++i)
    {
      std::cout << "Param " << i << ": " << leftRightFit->GetParameter(i) << std::endl;
    }
    std::cout << "Combined Fit Parameters:" << std::endl;
    for (int i = 0; i < combinedFit->GetNpar(); ++i)
    {
      std::cout << "Param " << i << ": " << combinedFit->GetParameter(i) << std::endl;
    }
    // Create a new function for just the polynomial part
    TF1 *polyPart;
    if (background_scheme == 0)  // poly4
    {
      polyPart = new TF1("polyPart", "pol4", limits[0], limits[1]);
      for (int j = 0; j < 5; ++j) polyPart->SetParameter(j, combinedFit->GetParameter(j + 3));
    }
    if (background_scheme == 1)  // poly3+poly2
    {
      ///////////////combined fit is not working. just use leftrightfit params for now
      /*
      polyPart = new TF1("polyPart", ONLYdoublePolyBG, limits[0], limits[1], 7);
      for (int j = 0; j < 4; ++j)
      {
        polyPart->SetParameter(j, combinedFit->GetParameter(j + 6));  // 6,7,8,9
      }
      for (int k = 0; k < 3; k++)
      {
        polyPart->SetParameter(k + 4, combinedFit->GetParameter(k + 10));  // 10,11,12
      }
      std::cout << "Polynomial Parameters:" << std::endl;
      for (int i = 0; i < polyPart->GetNpar(); ++i)
      {
        std::cout << "Param " << i << ": " << polyPart->GetParameter(i) << std::endl;
      }
      */
      polyPart = new TF1("polyPart", ONLYdoublePolyBG, limits[0], limits[1], 5);
      for (int j = 0; j < 4; ++j) polyPart->SetParameter(j, leftRightFit->GetParameter(j));
      polyPart->SetParameter(4, leftRightFit->GetParameter(4));

      std::cout << "Polynomial Parameters:" << std::endl;
      for (int i = 0; i < polyPart->GetNpar(); ++i)
      {
        std::cout << "Param " << i << ": " << polyPart->GetParameter(i) << std::endl;
      }
    }
    else if (background_scheme == 2)
    {
      polyPart = new TF1("polyPart", "log(x)", limits[0], limits[1]);
      for (int j = 0; j < 2; ++j) polyPart->SetParameter(j, combinedFit->GetParameter(j + 6));
    }
    else if (background_scheme == 3)
    {
      polyPart = new TF1("polyPart", "pol2", limits[0], limits[1]);
      for (int j = 0; j < 3; ++j) polyPart->SetParameter(j, combinedFit->GetParameter(j + 6));
    }
    else if (background_scheme == 4)
    {
      polyPart = new TF1("polyPart", "pol3", limits[0], limits[1]);
      for (int j = 0; j < 4; ++j) polyPart->SetParameter(j, combinedFit->GetParameter(j + 6));
    }
    else if (background_scheme == 5)
    {
      polyPart = new TF1("polyPart", "pol5", limits[0], limits[1]);
      for (int j = 0; j < 6; ++j) polyPart->SetParameter(j, combinedFit->GetParameter(j + 6));
    }
    else if (background_scheme == 6)
    {
      polyPart = new TF1("polyPart", "pol6", limits[0], limits[1]);
      for (int j = 0; j < 7; ++j) polyPart->SetParameter(j, combinedFit->GetParameter(j + 6));
    }

    // Create a new histogram to store the subtracted data
    TH1F *histSubtracted = (TH1F *) histF->Clone(Form("histSubtracted_%d", i));
    for (int j = 1; j <= histF->GetNbinsX(); ++j)
    {
      double x = histF->GetBinCenter(j);
      if (histF->GetBinContent(j) == 0)
      {
        histSubtracted->SetBinContent(j, 0);
      }
      else
      {
        double y = histF->GetBinContent(j) - polyPart->Eval(x);
        histSubtracted->SetBinContent(j, y);
      }
    }

    // Fit the subtracted histogram with the double Gaussian function
    TF1 *doubleGaussFit = new TF1("doubleGaussFit", doubleGauss, limits[0], limits[1], 10);
    doubleGaussFit->SetNpx(1000);
    doubleGaussFit->SetParameter(0, combinedFit->GetParameter(0));
    doubleGaussFit->SetParameter(1, combinedFit->GetParameter(1));
    doubleGaussFit->SetParameter(2, combinedFit->GetParameter(2));
    doubleGaussFit->SetParameter(6, limits[2]);
    doubleGaussFit->SetParameter(7, limits[3]);
    doubleGaussFit->SetParLimits(1, 0.11, 0.19);
    doubleGaussFit->SetParLimits(2, 0.01, 0.25);
    if (fitEtaPeak)
    {
      if (background_scheme == 0)  // poly4
      {
        doubleGaussFit->SetParameter(3, combinedFit->GetParameter(8));
        doubleGaussFit->SetParameter(4, combinedFit->GetParameter(9));
        doubleGaussFit->SetParameter(5, combinedFit->GetParameter(10));
      }
      else if (background_scheme == 1)  // poly3+poly2
      {
        doubleGaussFit->SetParameter(3, combinedFit->GetParameter(7));
        doubleGaussFit->SetParameter(4, combinedFit->GetParameter(8));
        doubleGaussFit->SetParameter(5, combinedFit->GetParameter(9));
      }
      else if (background_scheme == 2 || background_scheme == 3 || background_scheme == 4 || background_scheme == 5 || background_scheme == 6)
      {
        doubleGaussFit->SetParameter(3, combinedFit->GetParameter(3));
        doubleGaussFit->SetParameter(4, combinedFit->GetParameter(4));
        doubleGaussFit->SetParameter(5, combinedFit->GetParameter(5));
      }
      doubleGaussFit->SetParameter(8, limits[6]);
      doubleGaussFit->SetParameter(9, limits[7]);
      doubleGaussFit->SetParLimits(4, 0.55, 0.63);
      doubleGaussFit->SetParLimits(5, 0.01, 0.25);
    }
    histSubtracted->Fit(doubleGaussFit, "RE");

    // Store pion peak position and resolution; instead of combined fit use doubleGaussFit
    double pion_pt = (pt_min + pt_max) / 2.0;
    double pion_peak = combinedFit->GetParameter(1);
    double pion_peak_err = combinedFit->GetParError(1);
    double pion_res = combinedFit->GetParameter(2) / combinedFit->GetParameter(1);
    double pion_res_err = pion_res * sqrt(pow(combinedFit->GetParError(2) / combinedFit->GetParameter(2), 2) + pow(pion_peak_err / pion_peak, 2));

    pionPt.push_back(pion_pt);
    pionPeak.push_back(pion_peak);
    pionPeakErr.push_back(pion_peak_err);  //
    pionRes.push_back(pion_res);
    pionResErr.push_back(pion_res_err);

    // Store eta peak position and resolution if fitting eta
    if (fitEtaPeak)
    {
      double eta_peak, eta_peak_err, eta_res, eta_res_err, peak_ratio_err;
      if (background_scheme == 0)  // poly4
      {
        eta_peak = combinedFit->GetParameter(9);
        eta_peak_err = combinedFit->GetParError(9);
        eta_res = combinedFit->GetParameter(10) / combinedFit->GetParameter(9);
        eta_res_err = eta_res * sqrt(pow(combinedFit->GetParError(10) / combinedFit->GetParameter(10), 2) + pow(eta_peak_err / eta_peak, 2));
        peak_ratio_err = sqrt(pow(eta_peak_err / eta_peak, 2) + pow(pion_peak_err / pion_peak, 2));
      }
      else if (background_scheme == 1)  // poly3+poly2
      {
        eta_peak = combinedFit->GetParameter(4);
        eta_peak_err = combinedFit->GetParError(4);
        eta_res = combinedFit->GetParameter(5) / combinedFit->GetParameter(4);
        eta_res_err = eta_res * sqrt(pow(combinedFit->GetParError(5) / combinedFit->GetParameter(5), 2) + pow(eta_peak_err / eta_peak, 2));
        peak_ratio_err = sqrt(pow(eta_peak_err / eta_peak, 2) + pow(pion_peak_err / pion_peak, 2));
      }
      else if (background_scheme == 2 || background_scheme == 3 || background_scheme == 4 || background_scheme == 5 || background_scheme == 6)
      {
        eta_peak = combinedFit->GetParameter(4);
        eta_peak_err = combinedFit->GetParError(4);
        eta_res = combinedFit->GetParameter(5) / combinedFit->GetParameter(4);
        eta_res_err = eta_res * sqrt(pow(combinedFit->GetParError(5) / combinedFit->GetParameter(5), 2) + pow(eta_peak_err / eta_peak, 2));
        peak_ratio_err = sqrt(pow(eta_peak_err / eta_peak, 2) + pow(pion_peak_err / pion_peak, 2));
      }

      etaPeak.push_back(eta_peak);
      etaPeakErr.push_back(eta_peak_err);
      etaRes.push_back(eta_res);
      etaResErr.push_back(eta_res_err);
      PeakRatio.push_back(pion_peak / eta_peak);
      PeakRatioErr.push_back(peak_ratio_err);
      // std::cout << "Combined Fit Error Parameters:" << std::endl;
      // for (int i = 0; i < combinedFit->GetNpar(); ++i) {
      // std::cout << "Param " << i << ": " << pion_peak_err << " , " << pion_res_err <<" , " << eta_peak_err <<" , " << eta_res_err <<" , " << peak_ratio_err <<std::endl;
      //}
    }
    else
    {
      etaPeak.push_back(0);
      etaPeakErr.push_back(0);
      etaRes.push_back(0);
      etaResErr.push_back(0);
    }
    // Draw the fits and subtracted histograms
    TCanvas *c1 = new TCanvas(Form("c1_%s", ptRange.Data()), "Fits", 800, 600);
    histF->SetTitle(Form("Combined Fit; #it{m}_{#gamma#gamma} (GeV); dN/d#it{m}_{#gamma#gamma}; pT: %s", ptRange.Data()));
    histF->Draw("E");
    histF->SetMinimum(0.0);
    polyPart->SetLineColor(kRed);
    polyPart->Draw("SAME");
    combinedFit->SetLineColor(kBlack);
    combinedFit->Draw("SAME");
    //leftRightFit->SetLineColor(kGreen);
    //leftRightFit->Draw("SAME");
    gausFit->SetLineColor(kMagenta);
    gausFit->Draw("SAME");
    gausFit2->SetLineColor(kMagenta);
    gausFit2->Draw("SAME");

    TLegend *leg1 = new TLegend(0.5, 0.5, 0.95, 0.95);
    leg1->SetFillStyle(0);
    leg1->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
    leg1->AddEntry("", "pythia: p+p #sqrt{s_{NN}} = 200 GeV", "");
    leg1->AddEntry(polyPart, "Background Fit");
    leg1->AddEntry(combinedFit, "Combined Fit");
    //leg1->AddEntry(leftRightFit, "originalBG");
    leg1->AddEntry(gausFit, "originalGauss");
    leg1->Draw();
    leg1->SetTextAlign(32);
    c1->Update();
    c1->Print("2D_Histogram_Fits.pdf");

    TCanvas *c2 = new TCanvas(Form("c2_%s", ptRange.Data()), "Subtracted Peak", 800, 600);
    histSubtracted->SetTitle(Form("Background Subtracted Peak; #it{m}_{#gamma#gamma} (GeV); dN/d#it{m}_{#gamma#gamma}; pT: %s", ptRange.Data()));
    // hist->GetYaxis()->SetTitleOffset(1.5); // Adjust this value as needed
    // hist->GetYaxis()->SetLabelOffset(0.02); // Adjust if the labels overlap
    histSubtracted->Draw();
    histSubtracted->SetMinimum(0.0);
    histSubtracted->GetYaxis()->SetTitleOffset(1.5);
    TLegend *leg = new TLegend(0.5, 0.8, 0.93, 0.93);
    leg->SetFillStyle(0);
    leg->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
    leg->AddEntry("", "pythia: p+p #sqrt{s_{NN}} = 200 GeV", "");
    leg->SetTextAlign(32);
    histSubtracted->SetStats(0);

    TPaveText *pt2 = new TPaveText(0.6, 0.5, 0.93, 0.78, "NDC");
    pt2->SetFillColor(0);
    pt2->SetFillStyle(0);
    pt2->AddText(Form("#chi^{2}/NDF = %.2f", doubleGaussFit->GetChisquare() / doubleGaussFit->GetNDF()));
    pt2->AddText(Form("Mean = %.4f", doubleGaussFit->GetParameter(1)));
    pt2->AddText(Form("Sigma = %.4f", doubleGaussFit->GetParameter(2)));
    pt2->AddText(Form("Relative Width: %.2f%%", doubleGaussFit->GetParameter(2) * 100.0f / doubleGaussFit->GetParameter(1)));
    if (fitEtaPeak)
    {
      pt2->AddText(Form("Eta Mean = %.4f", doubleGaussFit->GetParameter(4)));
      pt2->AddText(Form("Eta Sigma = %.4f", doubleGaussFit->GetParameter(5)));
      pt2->AddText(Form("Eta Relative Width: %.2f%%", doubleGaussFit->GetParameter(5) * 100.0f / doubleGaussFit->GetParameter(4)));
    }
    pt2->Draw("SAME");
    c2->Print("2D_Histogram_Fits.pdf");

    // Create a canvas to display the fit parameters
    TCanvas *c3 = new TCanvas(Form("fitInfo_%s", ptRange.Data()), "Fit Parameters", 800, 600);
    TPaveText *fitInfo = new TPaveText(0.1, 0.1, 0.9, 0.9, "blNDC");
    fitInfo->SetTextAlign(12);
    fitInfo->SetFillColor(0);
    fitInfo->AddText(Form("Fit for pT range: %s", ptRange.Data()));
    fitInfo->AddText("Fit Parameters:");
    fitInfo->AddText(Form("Combined Fit Range = %f to %f", limits[0], limits[1]));
    fitInfo->AddText(Form("Pion Mean = %f +/- %f", combinedFit->GetParameter(1), combinedFit->GetParError(1)));
    fitInfo->AddText(Form("Pion Sigma = %f +/- %f", combinedFit->GetParameter(2), combinedFit->GetParError(2)));
    if (background_scheme == 0)  // poly4
    {
      fitInfo->AddText(Form("Eta Mean = %f +/- %f", combinedFit->GetParameter(9), combinedFit->GetParError(9)));
      fitInfo->AddText(Form("Eta Sigma= %f +/- %f", combinedFit->GetParameter(10), combinedFit->GetParError(10)));
    }
    else if (background_scheme == 1)  // poly3+poly2
    {
      fitInfo->AddText(Form("Eta Mean = %f +/- %f", combinedFit->GetParameter(4), combinedFit->GetParError(4)));
      fitInfo->AddText(Form("Eta Sigma = %f +/- %f", combinedFit->GetParameter(5), combinedFit->GetParError(5)));
    }
    else if (background_scheme == 2 || background_scheme == 3 || background_scheme == 4 || background_scheme == 5 || background_scheme == 6)
    {
      fitInfo->AddText(Form("Eta Mean = %f +/- %f", combinedFit->GetParameter(4), combinedFit->GetParError(4)));
      fitInfo->AddText(Form("Eta Sigma = %f +/- %f", combinedFit->GetParameter(5), combinedFit->GetParError(5)));
    }
    fitInfo->AddText(Form("Combined Chi2/NDF = %f / %d = %f", combinedFit->GetChisquare(), combinedFit->GetNDF(), combinedFit->GetChisquare() / combinedFit->GetNDF()));
    fitInfo->AddText(Form("Background Subtracted Peak Fit = %f to %f", limits[2], limits[3]));
    fitInfo->AddText(Form("Pion Mean = %f +/- %f", doubleGaussFit->GetParameter(1), doubleGaussFit->GetParError(1)));
    fitInfo->AddText(Form("Pion Sigma = %f +/- %f", doubleGaussFit->GetParameter(2), doubleGaussFit->GetParError(2)));
    fitInfo->AddText(Form("Relative Width: %f", doubleGaussFit->GetParameter(2) * 100.0f / doubleGaussFit->GetParameter(1)));
    fitInfo->AddText(Form("Chi2/NDF = %f / %d = %f", doubleGaussFit->GetChisquare(), doubleGaussFit->GetNDF(), doubleGaussFit->GetChisquare() / doubleGaussFit->GetNDF()));
    if (fitEtaPeak)
    {
      fitInfo->AddText(Form("Eta Peak Mean = %f +/- %f", doubleGaussFit->GetParameter(4), doubleGaussFit->GetParError(4)));
      fitInfo->AddText(Form("Eta Peak Sigma = %f +/- %f", doubleGaussFit->GetParameter(5), doubleGaussFit->GetParError(5)));
      fitInfo->AddText(Form("Eta Relative Width: %f", doubleGaussFit->GetParameter(5) * 100.0f / doubleGaussFit->GetParameter(4)));
    }
    fitInfo->Draw();
    c3->Print("2D_Histogram_Fits.pdf");

    appendtextfile(combinedFit, Form("Combined Fit_%s", ptRange.Data()), scale_factor);
    appendtextfile(doubleGaussFit, Form("subpgaus fit_%s", ptRange.Data()), scale_factor);

    std::cout << "Before deleting canvases" << std::endl;
    /*
    delete hist;
    delete histF;
    delete leftRightFit;
    delete gausFit;
    delete combinedFit;
    delete polyPart;
    delete histSubtracted;
    delete doubleGaussFit;
    delete leg1;
    delete leg;
    delete pt2;
    delete fitInfo;
    delete c1;
    delete c2;
    delete c3;
    */
    std::cout << "After deleting canvases" << std::endl;
  }
  std::cout << "Before TGraphErrors initialization" << std::endl;
  // Create TGraphErrors for the collected fit results
  int nPoints = pionPt.size();
  TGraphErrors *gPionPeak = new TGraphErrors(nPoints, &pionPt[0], &pionPeak[0], &pionPtErr[0], &pionPeakErr[0]);
  TGraphErrors *gPionRes = new TGraphErrors(nPoints, &pionPt[0], &pionRes[0], &pionPtErr[0], &pionResErr[0]);
  std::cout << "After TGraphErrors initialization" << std::endl;

  gPionPeak->SetTitle("Pion Peak Position; #it{pT}_{#gamma#gamma} (GeV); Pion Peak Position (GeV)");
  gPionRes->SetTitle("Pion Relative Resolution; #it{pT}_{#gamma#gamma} (GeV); Pion Relative Resolution");
  gPionPeak->GetYaxis()->SetTitleOffset(1.5);
  gPionRes->GetYaxis()->SetTitleOffset(1.5);
  // Draw TGraphErrors and add them to the PDF
  TCanvas *cPionPeak = new TCanvas("cPionPeak", "Pion Peak Position", 800, 600);
  gPionPeak->Draw("ALP");
  cPionPeak->Print("2D_Histogram_Fits.pdf");

  TCanvas *cPionRes = new TCanvas("cPionRes", "Pion Relative Width", 800, 600);
  gPionRes->Draw("ALP");
  cPionRes->Print("2D_Histogram_Fits.pdf");

  if (fitEtaPeak)
  {
    TGraphErrors *gEtaPeak = new TGraphErrors(nPoints, &pionPt[0], &etaPeak[0], &pionPtErr[0], &etaPeakErr[0]);
    TGraphErrors *gEtaRes = new TGraphErrors(nPoints, &pionPt[0], &etaRes[0], &pionPtErr[0], &etaResErr[0]);
    TGraphErrors *gPeakRatio = new TGraphErrors(nPoints, &pionPt[0], &PeakRatio[0], &pionPtErr[0], &PeakRatioErr[0]);
    gEtaPeak->SetTitle("Eta Peak Position; #it{pT}_{#gamma#gamma} (GeV); Eta Peak Position (GeV)");
    gEtaRes->SetTitle("Eta Relative Width; #it{pT}_{#gamma#gamma} (GeV); Eta Relative Width");
    gPeakRatio->SetTitle("Pion/Eta Mass Ratio; #it{pT}_{#gamma#gamma} (GeV); Pion/Eta Mass");
    gEtaPeak->GetYaxis()->SetTitleOffset(1.5);
    gEtaRes->GetYaxis()->SetTitleOffset(1.5);
    gPeakRatio->GetYaxis()->SetTitleOffset(1.5);

    TCanvas *cEtaPeak = new TCanvas("cEtaPeak", "Eta Peak Position", 800, 600);
    gEtaPeak->Draw("ALPE");
    cEtaPeak->Print("2D_Histogram_Fits.pdf");

    TCanvas *cEtaRes = new TCanvas("cEtaRes", "Eta Relative Width", 800, 600);
    gEtaRes->Draw("ALPE");
    cEtaRes->Print("2D_Histogram_Fits.pdf");

    TCanvas *cPeakRatio = new TCanvas("cPeakRatio", "Pion/Eta Mass Ratio", 800, 600);
    gPeakRatio->Draw("ALPE");
    cPeakRatio->Print("2D_Histogram_Fits.pdf");

    // Clean up for eta graphs
    delete cEtaPeak;
    delete cEtaRes;
    delete gEtaPeak;
    delete gEtaRes;
    delete gPeakRatio;
    delete cPeakRatio;
  }

  // Close the PDF file
  dummyCanvas->Print("2D_Histogram_Fits.pdf]");

  // Clean up
  delete file;
  delete dummyCanvas;
  delete cPionPeak;
  delete cPionRes;
  delete gPionPeak;
  delete gPionRes;
}

void bgsub(double scale_factor = 1, float polyL = 0.05, float polygauss1L = 0.08, float gauss1L = 0.11, float gauss1R = 0.19, float polygauss1R = 0.3, float polygauss2L = 0.5, float gauss2L = 0.55, float gauss2R = 0.65, float polygauss2R = 0.7, float polyR = 1.0, int startBin = 1, int endBin = -1, int projectionBins = 1, int rebinFactor = 1, int background_scheme = 0, bool fitEtaPeak = false, bool dynamic_left = true, bool variable_bins = false)
{
  // code is run with a command like this:
  // root bgsub.C'(1.4,0.05,0.05,0.11,0.19,0.35,0.35,0.52,0.68,1.0,1.0,9,29,1,1,1,true,true,true)'
  // background scheme:0=poly4, 1= poly3 + poly2,2=log,3=poly2,4=poly3,5=poly5

  // Fit limits for the polynomial and Gaussian fits
  std::vector<float> limits = {
      polyL, polyR,              // 0,1 Polynomial fit range: left and right limits
      gauss1L, gauss1R,          // 2-3 First Gaussian fit range: left and right limits
      polygauss1L, polygauss1R,  // 4,5 Exclusion zone for left and right polynomials: first gaussian
      gauss2L, gauss2R,          // 6,7 Second Gaussian fit range (if fitting eta peak): left and right limits
      polygauss2L, polygauss2R   // 8,9 Exclusion zone for left and right polynomials: second gaussian
  };

  //if (variable_bins)
  //{
    // std::vector<double> nonUniformBins = {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.36, 0.40, 0.44, 0.48, 0.52, 0.56, 0.60, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 1.0};//, 1.04, 1.08, 1.12, 1.16, 1.2
    //std::vector<double> nonUniformBins=nuBins;
  fit_2d_histogram(scale_factor, limits, fitEtaPeak, startBin, endBin, projectionBins, rebinFactor, dynamic_left, background_scheme, variable_bins);
  //}
  //else
  //{
    //std::vector<double> nonUniformBins;
    //fit_2d_histogram(scale_factor, limits, fitEtaPeak, startBin, endBin, projectionBins, rebinFactor, dynamic_left, background_scheme, nonUniformBins);
  //}

  // fit_2d_histogram(scale_factor, limits, fitEtaPeak, startBin, endBin, projectionBins, rebinFactor, dynamic_left, background_scheme, nonUniformBins);
  //  return 0;
  gApplication->Terminate(0);
}
