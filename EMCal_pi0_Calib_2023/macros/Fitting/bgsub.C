// root includes
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
// general includes
#include <fstream>
#include <iostream>
#include <vector>
// sphenix includes
#include "sPhenixStyle.C"
#include "sPhenixStyle.h"

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
  //double gauss1 = par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));
  double gauss1 = 0;
  if (x[0] >= 0.1 && x[0] <= 0.2)
  {  // Check if x is in the range of the first Gaussian
  double gauss1 = par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));
  }

  // Second Gaussian part (e.g., eta peak)
  //double gauss2 = par[7] * exp(-0.5 * pow((x[0] - par[8]) / par[9], 2));
  double gauss2 = 0;
  if (x[0] >= 0.5 && x[0] <= 0.65)
  {  // Check if x is in the range of the first Gaussian
  double gauss2 = par[7] * exp(-0.5 * pow((x[0] - par[8]) / par[9], 2));
  }

  // Polynomial part
  // First Gaussian part (e.g., pion peak)
  double poly1 = 0;
  if (x[0] >= 0.05 && x[0] <= 0.3)
  {  // Check if x is in the range of the first Gaussian
  double poly1 = par[3] + par[4] * x[0] + par[5] * x[0] * x[0] + par[6] * x[0] * x[0] * x[0];
  }

  // Second Gaussian part (e.g., eta peak)
  double poly2 = 0;
  if (x[0] > 0.3 && x[0] <= 1.0)
  {
  double poly2 = par[10] + par[11] * x[0] + par[12] * x[0] * x[0];
  }

  return gauss1 + gauss2 + poly1 + poly2;
}

double doubleGauss(double *x, double *par)
{
  // First Gaussian part (e.g., pion peak)
  double gauss1 = 0;
  if (x[0] >= par[6] && x[0] <= par[7])
  {  // Check if x is in the range of the first Gaussian
    gauss1 = par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));
  }

  // Second Gaussian part (e.g., eta peak)
  double gauss2 = 0;
  if (x[0] >= par[8] && x[0] <= par[9])
  {  // Check if x is in the range of the second Gaussian
    gauss2 = par[3] * exp(-0.5 * pow((x[0] - par[4]) / par[5], 2));
  }

  return gauss1 + gauss2;
}

double doublePolyBG(double *x, double *par)
{
  // First Gaussian part (e.g., pion peak)
  double poly1 = 0;
  if (x[0] >= par[4] && x[0] <= par[5])
  {  // Check if x is in the range of the first Gaussian
    poly1 = par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0];
  }
  else
  {
    poly1 = par[6] + par[7] * x[0] + par[8] * x[0] * x[0];
  }
  // Check if x is in the range of any Gaussian fit
  if ((x[0] >= 0.11 && x[0] <= 0.19) || (x[0] >= 0.52 && x[0] <= 0.68))
  {
    TF1::RejectPoint();
    return 0;
  }

  // Second Gaussian part (e.g., eta peak)
  // double poly2 = 0;
  // if (x[0] >= par[9] && x[0] <= par[10])
  //{  // Check if x is in the range of the second Gaussian
  //  poly2 = par[6] + par[7] * x[0] + par[8] * x[0] * x[0];
  //}

  return poly1;  // + poly2;
}
double ONLYdoublePolyBG(double *x, double *par)
{
  //double  poly1 = par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0];
  //double  poly2 = par[4] + par[5] * x[0] + par[6] * x[0] * x[0];
  double poly1 = 0;
  if (x[0] >= 0.05 && x[0] <= 0.3)
  {  // Check if x is in the range of the first Gaussian
    poly1 = par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0];
  }
  else
  {
    poly1 = par[4] + par[5] * x[0] + par[6] * x[0] * x[0];
  }
  return poly1;// + poly2;
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

void fit_histogram(double scale_factor = 1, float leftmost_gauslimit = 0.05, float rightmost_gauslimit = 0.9, bool fitEtaPeak = false)
{
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
  SetsPhenixStyle();

  // Open the ROOT file and get the histogram
  TFile *file = new TFile("/sphenix/u/nkumar/analysis/EMCal_pi0_Calib_2023/macros/condor/output/merged_file.root");
  TH1F *hist = (TH1F *) file->Get("h_InvMass");

  // Rebin the histogram to have 'numBins' bins
  int numBins = 120;
  int currentNumBins = hist->GetNbinsX();
  double currentXMax = hist->GetXaxis()->GetXmax();
  int rebinFactor = currentNumBins / numBins;
  if (rebinFactor > 1)
  {
    std::cout << "current nbins: " << currentNumBins << " requested nbins: " << numBins << " rebin by: " << rebinFactor << std::endl;
    hist->Rebin(rebinFactor);
    std::cout << "new nbin check: " << hist->GetNbinsX() << std::endl;
  }

  // Overall limits
  float rightmost_limit = 0.9;  // fit range limit
  float leftmost_limit = 0.05;  // fit range limit. normally 0.05
  // Limits on gauss and poly
  float leftpolylim = 0.11;
  float rightpolylim = 0.19;

  hist->GetXaxis()->SetRangeUser(0, 1.0);
  // scale_histogram_errors(hist, scale_factor);

  // Fit left and right regions with a polynomial, excluding Gaussian regions
  TF1 *leftRightFit;
  if (fitEtaPeak)
  {
    leftRightFit = new TF1("leftRightFit", leftRightPolynomial, leftmost_limit, rightmost_limit, 9);
    leftRightFit->SetParameter(5, leftpolylim);
    leftRightFit->SetParameter(6, rightpolylim);
    leftRightFit->SetParameter(7, 0.5);  // Rough estimate for eta peak start
    leftRightFit->SetParameter(8, 0.7);  // Rough estimate for eta peak end
  }
  else
  {
    leftRightFit = new TF1("leftRightFit", leftRightPolynomial, leftmost_limit, rightmost_limit, 7);
    leftRightFit->SetParameter(5, leftpolylim);
    leftRightFit->SetParameter(6, rightpolylim);
  }
  hist->Fit(leftRightFit, "R");

  // Fit Gaussian in the specified range
  TF1 *gausFit = new TF1("gausFit", "gaus", leftpolylim, rightpolylim);  // leftpolylim, rightpolylim
  hist->Fit(gausFit, "R");

  // Combined Gaussian + Polynomial fit
  TF1 *combinedFit;
  if (fitEtaPeak)
  {
    combinedFit = new TF1("combinedFit", combinedFunctionDoubleGauss, leftmost_limit, rightmost_limit, 11);  // 2 Gaussians + 1 polynomial = 3 + 3 + 5
  }
  else
  {
    combinedFit = new TF1("combinedFit", combinedFunction, leftmost_limit, rightmost_limit, 8);
  }

  // Set initial parameters from previous fits
  for (int i = 0; i < 3; ++i) combinedFit->SetParameter(i, gausFit->GetParameter(i));
  for (int i = 3; i < 8; ++i) combinedFit->SetParameter(i, leftRightFit->GetParameter(i - 3));

  if (fitEtaPeak)
  {
    // Set initial guesses for the second Gaussian (eta peak)
    combinedFit->SetParameter(8, gausFit->GetParameter(0) / 3);  // Assume a smaller amplitude
    combinedFit->SetParameter(9, 0.6);                           // Rough guess for eta peak mean
    combinedFit->SetParameter(10, 0.05);                         // Rough guess for eta peak sigma
  }

  // Fit the combined function
  hist->Fit(combinedFit, "R");

  double chi2 = combinedFit->GetChisquare();
  double ndf = combinedFit->GetNDF();
  double chi2ndf = chi2 / ndf;

  std::cout << "Chi-squared: " << chi2 << std::endl;
  std::cout << "Number of Degrees of Freedom: " << ndf << std::endl;
  std::cout << "Chi-squared/NDF: " << chi2ndf << std::endl;

  // Create a new function for just the polynomial part
  TF1 *polyPart = new TF1("polyPart", "pol4", leftmost_limit, rightmost_limit);
  for (int i = 0; i < 5; ++i) polyPart->SetParameter(i, combinedFit->GetParameter(i + 3));

  // Create a new histogram to store the subtracted data
  TH1F *histSubtracted = (TH1F *) hist->Clone("histSubtracted");
  for (int i = 1; i <= hist->GetNbinsX(); ++i)
  {
    double x = hist->GetBinCenter(i);
    double y = hist->GetBinContent(i) - polyPart->Eval(x);
    histSubtracted->SetBinContent(i, y);
  }

  TF1 *gausFit2 = new TF1("gausFit2", "gaus", leftmost_gauslimit, rightmost_gauslimit);  // leftmost_limit, 0.25
  for (int i = 0; i < 3; ++i) gausFit2->SetParameter(i, combinedFit->GetParameter(i));
  histSubtracted->Fit(gausFit2, "R");

  double chi2_s = gausFit2->GetChisquare();
  double ndf_s = gausFit2->GetNDF();
  double chi2ndf_s = chi2_s / ndf_s;

  std::cout << "Chi-squared: " << chi2_s << std::endl;
  std::cout << "Number of Degrees of Freedom: " << ndf_s << std::endl;
  std::cout << "Chi-squared/NDF: " << chi2ndf_s << std::endl;
  std::cout << "Relative Width: " << 100 * gausFit2->GetParameter(2) / gausFit2->GetParameter(1) << " %" << std::endl;

  // Open PDF file
  TCanvas *dummyCanvas = new TCanvas();
  dummyCanvas->Print("fit_results.pdf[");

  // Draw everything
  TCanvas *c1 = new TCanvas("c1", "Fits", 800, 600);
  hist->Draw("E");
  hist->SetTitle("pi0 and eta Inv. Mass Peaks; Inv. Mass (GeV); Counts");

  polyPart->SetLineColor(kRed);
  polyPart->Draw("SAME");

  combinedFit->SetLineColor(kBlack);
  combinedFit->Draw("SAME");

  TLegend *leg1 = new TLegend(0.5, 0.5, 0.95, 0.95);
  leg1->SetFillStyle(0);
  leg1->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
  leg1->AddEntry("", "pythia: p+p #sqrt{s_{NN}} = 200 GeV", "");
  leg1->AddEntry(polyPart, "Background Fit");
  leg1->AddEntry(combinedFit, "Combined Fit");
  leg1->Draw();
  leg1->SetTextAlign(32);
  c1->Update();
  c1->Print("fit_results.pdf");

  TCanvas *c2 = new TCanvas("c2", "Subtracted Peak", 800, 600);
  histSubtracted->Draw();
  histSubtracted->SetMinimum(0.0);
  histSubtracted->SetTitle("Background Subtracted Peak; Inv. Mass (GeV); Counts (Background subtracted)");
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
  pt2->AddText(Form("#chi^{2}/NDF = %.2f", gausFit2->GetChisquare() / gausFit2->GetNDF()));
  pt2->AddText(Form("Mean = %.4f", gausFit2->GetParameter(1)));
  pt2->AddText(Form("Sigma = %.4f", gausFit2->GetParameter(2)));
  pt2->AddText(Form("Relative Width: %.2f%%", gausFit2->GetParameter(2) * 100.0f / gausFit2->GetParameter(1)));
  pt2->Draw("SAME");
  gPad->Modified();
  leg->Draw("Same");
  gPad->Update();
  c2->Print("fit_results.pdf");

  appendtextfile(combinedFit, "Combined Fit", scale_factor);
  appendtextfile(gausFit2, "subpgaus fit", scale_factor);

  // Fit Parameters Canvas
  TCanvas *c3 = new TCanvas("canvas2", "Fit Parameters", 800, 600);
  c3->cd();

  TPaveText *pt = new TPaveText(0.1, 0.1, 0.9, 0.9, "blNDC");
  pt->SetTextAlign(12);
  pt->SetFillColor(0);

  pt->AddText("Data Fit");
  pt->AddText("Fit Parameters:");
  pt->AddText(Form("Combined Fit Range = %f to %f", leftmost_limit, rightmost_limit));
  pt->AddText(Form("Pion Peak Mean = %f +/- %f", combinedFit->GetParameter(1), combinedFit->GetParError(1)));
  pt->AddText(Form("Pion Peak Sigma = %f +/- %f", combinedFit->GetParameter(2), combinedFit->GetParError(2)));
  if (fitEtaPeak)
  {
    pt->AddText(Form("Eta Peak Mean = %f +/- %f", combinedFit->GetParameter(9), combinedFit->GetParError(9)));
    pt->AddText(Form("Eta Peak Sigma = %f +/- %f", combinedFit->GetParameter(10), combinedFit->GetParError(10)));
    pt->AddText(Form("Eta Peak Relative Width = %f", combinedFit->GetParameter(10) / combinedFit->GetParameter(9)));
  }
  pt->AddText(Form("Background Subtracted Peak Fit = %f to %f", leftmost_gauslimit, rightmost_gauslimit));
  pt->AddText(Form("Pion Mean = %f +/- %f", gausFit2->GetParameter(1), gausFit2->GetParError(1)));
  pt->AddText(Form("Pion Sigma = %f +/- %f", gausFit2->GetParameter(2), gausFit2->GetParError(2)));
  pt->AddText(Form("Pion Relative Width: %f", gausFit2->GetParameter(2) * 100.0f / gausFit2->GetParameter(1)));
  pt->AddText(Form("Pion Chi2/NDF = %f / %d = %f", gausFit2->GetChisquare(), gausFit2->GetNDF(), gausFit2->GetChisquare() / gausFit2->GetNDF()));
  if (fitEtaPeak)
  {
    // pt->AddText(Form("Pion Mean = %f +/- %f", gausFit2->GetParameter(1), gausFit2->GetParError(1)));
    // pt->AddText(Form("Pion Sigma = %f +/- %f", gausFit2->GetParameter(2), gausFit2->GetParError(2)));
    // pt->AddText(Form("Pion Relative Width: %f", gausFit2->GetParameter(2) * 100.0f / gausFit2->GetParameter(1)));
    // pt->AddText(Form("Pion Chi2/NDF = %f / %d = %f", gausFit2->GetChisquare(), gausFit2->GetNDF(), gausFit2->GetChisquare() / gausFit2->GetNDF()));
  }

  pt->Draw();
  c3->Print("fit_results.pdf");

  // Close PDF file
  dummyCanvas->Print("fit_results.pdf]");

  delete file;
  delete c1;
  delete c2;
  delete c3;
  delete gausFit;
  delete gausFit2;
  delete polyPart;
  delete leftRightFit;
  delete combinedFit;
  delete leg;
  delete leg1;
  delete dummyCanvas;
}

void fit_2d_histogram(Double_t scale_factor, std::vector<float> &limits, bool fitEtaPeak = false, int startBin = 1, int endBin = -1, int projectionBins = 1, int rebinFactor = 1, bool dynamic_left = false, int background_scheme = 0)
{
  // more thorough minimizer for fit

  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2","Fumili");//,"Simplex"
  // Set the global fit strategy
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);
  SetsPhenixStyle();

  // Ensure the limits vector has the correct size
  if (limits.size() < 6 || (fitEtaPeak && limits.size() < 10))
  {
    std::cerr << "Insufficient limits provided. Expected 6 (or 10 if fitting eta peak) values." << std::endl;
    return;
  }

  // Open the ROOT file and get the 2D histogram
  TFile *file = new TFile("/sphenix/u/nkumar/analysis/EMCal_pi0_Calib_2023/macros/condor/output/merged_file.root");
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
    if (rebinFactor > 1)
    {
      hist->Rebin(rebinFactor);
    }

    hist->Scale(1 / (hist->GetBinLowEdge(2) - hist->GetBinLowEdge(1)));  // divide bin content by width so it becomes dN/dM

    // Determine the leftmost point with a value in the projection histograms
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

    // Get the pt range for the current slice
    double pt_min = hist2D->GetXaxis()->GetBinLowEdge(i);
    double pt_max = hist2D->GetXaxis()->GetBinUpEdge(lastBin);
    TString ptRange = Form("pt_%.2f-%.2f_GeV", pt_min, pt_max);

    // Set histogram range and scale errors
    // hist->GetXaxis()->SetRangeUser(0, 1.0);
    scale_histogram_errors(hist, scale_factor);

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
      if (background_scheme == 1)  // poly3+poly2
      {
        leftRightFit = new TF1("leftRightFit", doublePolyBG, limits[0], limits[1], 11);
        leftRightFit->SetParameter(4, limits[0]);  // left poly1 lim
        leftRightFit->SetParameter(5, 0.3);       // right poly2 lim
        // leftRightFit->SetParameter(9, 0.35);//left poly2 lim
        // leftRightFit->SetParameter(10, limits[1]);//right poly2 lim
      }
    }
    else
    {
      leftRightFit = new TF1("leftRightFit", leftRightPolynomial, limits[0], limits[1], 7);
      leftRightFit->SetParameter(5, limits[4]);
      leftRightFit->SetParameter(6, limits[5]);
    }
    hist->Fit(leftRightFit, "R");

    // Fit first Gaussian in the specified range
    TF1 *gausFit = new TF1("gausFit", "gaus", limits[2], limits[3]);
    gausFit->SetParLimits(1, 0.11, 0.19);
    gausFit->SetParLimits(2, 0.05, 0.25);
    hist->Fit(gausFit, "R");

    // Combined Gaussian + Polynomial fit
    TF1 *combinedFit;
    if (fitEtaPeak)
    {
      if (background_scheme == 0)  // poly4
      {
        combinedFit = new TF1("combinedFit", combinedFunctionDoubleGauss, limits[0], limits[1], 11);  // 2 Gaussians + 1 polynomial = 3 + 3 + 5
      }
      if (background_scheme == 1)  // poly3+poly2
      {
        combinedFit = new TF1("combinedFit", combinedFunctionDoubleGaussDoublePoly, limits[0], limits[1], 13);  // 2 Gaussians + 1 poly3 +1poly2 = 3 + 3 + 4 + 3=13
        // if using fit limits for parts you need 4 more paramms:13,14,15,16
      }
    }
    else
    {
      combinedFit = new TF1("combinedFit", combinedFunction, limits[0], limits[1], 8);
    }

    // Set initial parameters from previous fits
    for (int j = 0; j < 3; ++j){
      combinedFit->SetParameter(j, gausFit->GetParameter(j));
    } 
    combinedFit->SetParLimits(0, 0, gausFit->GetParameter(0) *1.25);//gausFit->GetParameter(0) *0.5
    combinedFit->SetParLimits(1, 0.11, 0.19);
    combinedFit->SetParLimits(2, 0.01, 0.11);
    // for (int j = 3; j < 8; ++j) combinedFit->SetParameter(j, leftRightFit->GetParameter(j - 3));

    // Fit second Gaussian in the specified range
    TF1 *gausFit2 = new TF1("gausFit2", "gaus", limits[6], limits[7]);
    gausFit2->SetParLimits(1, 0.55, 0.63);
    gausFit2->SetParLimits(2, 0.05, 0.25);
    hist->Fit(gausFit2, "R");
    if (fitEtaPeak)
    {
      if (background_scheme == 0)  // poly4
      {
        for (int j = 0; j < 3; ++j) combinedFit->SetParameter(j+8, gausFit2->GetParameter(j));
        for (int j = 0; j < 5; ++j) combinedFit->SetParameter(j+3, leftRightFit->GetParameter(j));
        // combinedFit->SetParLimits(8, 10, gausFit->GetParameter(0) / 6);
        combinedFit->SetParLimits(9, 0.55, 0.63);
        combinedFit->SetParLimits(10, 0.05, 0.25);
      }
      if (background_scheme == 1)  // poly3+poly2
      {
        for (int j = 0; j < 4; ++j){
          combinedFit->SetParameter(j + 3, leftRightFit->GetParameter(j));
        } 
        for (int j = 0; j < 3; ++j){
          combinedFit->SetParameter(j + 10, leftRightFit->GetParameter(j + 6));
        } 
        //combinedFit->SetParLimits(3, leftRightFit->GetParameter(0)*0.95, leftRightFit->GetParameter(0)*1.05);
        //combinedFit->SetParLimits(4, leftRightFit->GetParameter(1)*0.95, leftRightFit->GetParameter(1)*1.05);
        //combinedFit->SetParLimits(5, leftRightFit->GetParameter(2)*0.95, leftRightFit->GetParameter(2)*1.05);
        //combinedFit->SetParLimits(6, leftRightFit->GetParameter(3)*0.95, leftRightFit->GetParameter(3)*1.05);
        //combinedFit->SetParLimits(10, leftRightFit->GetParameter(6)*0.95, leftRightFit->GetParameter(6)*1.05);
        //combinedFit->SetParLimits(11, leftRightFit->GetParameter(7)*0.95, leftRightFit->GetParameter(7)*1.05);
        //combinedFit->SetParLimits(12, leftRightFit->GetParameter(8)*0.95, leftRightFit->GetParameter(8)*1.05);
        //combinedFit->SetParameter(13,limits[0]);
        //combinedFit->SetParameter(14,0.3);
        //combinedFit->SetParameter(15,0.3);
        //combinedFit->SetParameter(16,limits[1]);
        for (int j = 0; j < 3; ++j){
          combinedFit->SetParameter(j+7, gausFit2->GetParameter(j));
        } 
        //combinedFit->SetParLimits(7, 0, gausFit2->GetParameter(0) *1.05);//gausFit2->GetParameter(0) *0.5
        combinedFit->SetParLimits(8, 0.55, 0.63);
        combinedFit->SetParLimits(9, 0.01, 0.25);
      }
    }
    else
    {
      for (int j = 3; j < 8; ++j) combinedFit->SetParameter(j, leftRightFit->GetParameter(j - 3));
    }

    // Fit the combined function
    hist->Fit(combinedFit, "R");
  // After fitting
  std::cout << "Combined Fit Parameters:" << std::endl;
  for (int i = 0; i < combinedFit->GetNpar(); ++i) {
      std::cout << "Param " << i << ": " << combinedFit->GetParameter(i) << std::endl;
  }
    // Store pion peak position and resolution
    double pion_pt = (pt_min + pt_max) / 2.0;
    double pion_peak = combinedFit->GetParameter(1);
    double pion_peak_err = combinedFit->GetParError(1);
    double pion_res = combinedFit->GetParameter(2) / combinedFit->GetParameter(1);
    double pion_res_err = pion_res * sqrt(pow(combinedFit->GetParError(2) / combinedFit->GetParameter(2), 2) + pow(pion_peak_err / pion_peak, 2));

    pionPt.push_back(pion_pt);
    pionPeak.push_back(pion_peak);
    pionPeakErr.push_back(pion_peak_err);//
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
      if (background_scheme == 1)  // poly3+poly2
      {
        eta_peak = combinedFit->GetParameter(8);
        eta_peak_err = combinedFit->GetParError(8);
        eta_res = combinedFit->GetParameter(9) / combinedFit->GetParameter(8);
        eta_res_err = eta_res * sqrt(pow(combinedFit->GetParError(9) / combinedFit->GetParameter(9), 2) + pow(eta_peak_err / eta_peak, 2));
        peak_ratio_err = sqrt(pow(eta_peak_err / eta_peak, 2) + pow(pion_peak_err / pion_peak, 2));
      }

      etaPeak.push_back(eta_peak);
      etaPeakErr.push_back(eta_peak_err);
      etaRes.push_back(eta_res);
      etaResErr.push_back(eta_res_err);
      PeakRatio.push_back(pion_peak / eta_peak);
      PeakRatioErr.push_back(peak_ratio_err);
      std::cout << "Combined Fit Error Parameters:" << std::endl;
      for (int i = 0; i < combinedFit->GetNpar(); ++i) {
        std::cout << "Param " << i << ": " << pion_peak_err << " , " << pion_res_err <<" , " << eta_peak_err <<" , " << eta_res_err <<" , " << peak_ratio_err <<std::endl;
      }
    }
    else
    {
      etaPeak.push_back(0);
      etaPeakErr.push_back(0);
      etaRes.push_back(0);
      etaResErr.push_back(0);
    }

    // Create a new function for just the polynomial part
    TF1 *polyPart;
    //TF1 *gauspoly3 = new TF1("gpol3", "pol3", limits[0], 0.3);
    //TF1 *gauspoly2 = new TF1("gpol2", "pol2", 0.3, limits[1]);
    if (background_scheme == 0)  // poly4
    {
      polyPart = new TF1("polyPart", "pol4", limits[0], limits[1]);
      for (int j = 0; j < 5; ++j) polyPart->SetParameter(j, combinedFit->GetParameter(j + 3));
    }
    if (background_scheme == 1)  // poly3+poly2
    {
      //for (int j = 0; j < 4; ++j) gauspoly3->SetParameter(j, combinedFit->GetParameter(j + 3));       // 3,4,5,6
      //for (int j = 0; j < 3; ++j) gauspoly2->SetParameter(j, combinedFit->GetParameter(j + 10));  // 10,11,12
      //polyPart = new TF1("polyPart", "gpol3+gpol2", limits[0], limits[1]);
      polyPart = new TF1("polyPart", ONLYdoublePolyBG, limits[0], limits[1], 7);
      for (int j = 0; j < 4; ++j){
        polyPart->SetParameter(j, combinedFit->GetParameter(j + 3));       // 3,4,5,6
      } 
      for (int k = 0; k < 3; k++){
        polyPart->SetParameter(k + 4, combinedFit->GetParameter(k + 10));  // 10,11,12
      }
      std::cout << "Polynomial Parameters:" << std::endl;
      for (int i = 0; i < polyPart->GetNpar(); ++i) {
        std::cout << "Param " << i << ": " << polyPart->GetParameter(i) << std::endl;
      } 
    }

    // Create a new histogram to store the subtracted data
    TH1F *histSubtracted = (TH1F *) hist->Clone(Form("histSubtracted_%d", i));
    for (int j = 1; j <= hist->GetNbinsX(); ++j)
    {
      double x = hist->GetBinCenter(j);
      if (hist->GetBinContent(j) == 0)
      {
        histSubtracted->SetBinContent(j, 0);
      }
      else
      {
        double y = hist->GetBinContent(j) - polyPart->Eval(x);
        histSubtracted->SetBinContent(j, y);
      }
    }

    // Fit the subtracted histogram with the double Gaussian function
    TF1 *doubleGaussFit = new TF1("doubleGaussFit", doubleGauss, limits[0], limits[1], 10);
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
      if (background_scheme == 1)  // poly3+poly2
      {
        doubleGaussFit->SetParameter(3, combinedFit->GetParameter(7));
        doubleGaussFit->SetParameter(4, combinedFit->GetParameter(8));
        doubleGaussFit->SetParameter(5, combinedFit->GetParameter(9));
      }
      doubleGaussFit->SetParameter(8, limits[6]);
      doubleGaussFit->SetParameter(9, limits[7]);
      doubleGaussFit->SetParLimits(4, 0.55, 0.63);
      doubleGaussFit->SetParLimits(5, 0.01, 0.25);
    }
    histSubtracted->Fit(doubleGaussFit, "R");

    // Draw the fits and subtracted histograms
    TCanvas *c1 = new TCanvas(Form("c1_%s", ptRange.Data()), "Fits", 800, 600);
    hist->SetTitle(Form("Combined Fit; Inv. Mass (GeV); Counts; pT: %s", ptRange.Data()));
    hist->Draw("E");
    hist->SetMinimum(0.0);
    polyPart->SetLineColor(kRed);
    polyPart->Draw("SAME");
    combinedFit->SetLineColor(kBlack);
    combinedFit->Draw("SAME");
    leftRightFit->SetLineColor(kGreen);
    leftRightFit->Draw("SAME");
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
    leg1->AddEntry(leftRightFit, "originalBG");
    leg1->AddEntry(gausFit, "originalGauss");
    leg1->Draw();
    leg1->SetTextAlign(32);
    c1->Update();
    c1->Print("2D_Histogram_Fits.pdf");

    TCanvas *c2 = new TCanvas(Form("c2_%s", ptRange.Data()), "Subtracted Peak", 800, 600);
    histSubtracted->SetTitle(Form("Background Subtracted Peak; Inv. Mass (GeV); Counts (Background subtracted); pT: %s", ptRange.Data()));
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
    fitInfo->AddText(Form("Data Fit for pT range: %s", ptRange.Data()));
    fitInfo->AddText("Fit Parameters:");
    fitInfo->AddText(Form("Combined Fit Range = %f to %f", limits[0], limits[1]));
    fitInfo->AddText(Form("Pion Mean = %f +/- %f", combinedFit->GetParameter(1), combinedFit->GetParError(1)));
    fitInfo->AddText(Form("Pion Sigma = %f +/- %f", combinedFit->GetParameter(2), combinedFit->GetParError(2)));
    if (background_scheme == 0)  // poly4
    {
      fitInfo->AddText(Form("Eta Mean = %f +/- %f", combinedFit->GetParameter(9), combinedFit->GetParError(9)));
      fitInfo->AddText(Form("Eta Sigma= %f +/- %f", combinedFit->GetParameter(10), combinedFit->GetParError(10)));
    }
    if (background_scheme == 1)  // poly3+poly2
    {
      fitInfo->AddText(Form("Eta Mean = %f +/- %f", combinedFit->GetParameter(8), combinedFit->GetParError(8)));
      fitInfo->AddText(Form("Eta Sigma = %f +/- %f", combinedFit->GetParameter(9), combinedFit->GetParError(9)));
    }
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

    delete hist;
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
  }

  // Create TGraphErrors for the collected fit results
  int nPoints = pionPt.size();
  TGraphErrors *gPionPeak = new TGraphErrors(nPoints, &pionPt[0], &pionPeak[0], &pionPtErr[0], &pionPeakErr[0]);
  TGraphErrors *gPionRes = new TGraphErrors(nPoints, &pionPt[0], &pionRes[0], &pionPtErr[0], &pionResErr[0]);

  gPionPeak->SetTitle("Pion Peak Position; pT (GeV/c); Pion Peak Position (GeV/c^2)");
  gPionRes->SetTitle("Pion Relative Resolution; pT (GeV/c); Pion Relative Resolution");

  // Draw TGraphErrors and add them to the PDF
  TCanvas *cPionPeak = new TCanvas("cPionPeak", "Pion Peak Position", 800, 600);
  gPionPeak->Draw("ALP");
  cPionPeak->Print("2D_Histogram_Fits.pdf");

  TCanvas *cPionRes = new TCanvas("cPionRes", "Pion Relative Resolution", 800, 600);
  gPionRes->Draw("ALP");
  cPionRes->Print("2D_Histogram_Fits.pdf");

  if (fitEtaPeak)
  {
    TGraphErrors *gEtaPeak = new TGraphErrors(nPoints, &pionPt[0], &etaPeak[0], &pionPtErr[0], &etaPeakErr[0]);
    TGraphErrors *gEtaRes = new TGraphErrors(nPoints, &pionPt[0], &etaRes[0], &pionPtErr[0], &etaResErr[0]);
    TGraphErrors *gPeakRatio = new TGraphErrors(nPoints, &pionPt[0], &PeakRatio[0], &pionPtErr[0], &PeakRatioErr[0]);
    gEtaPeak->SetTitle("Eta Peak Position; pT (GeV/c); Eta Peak Position (GeV/c^2)");
    gEtaRes->SetTitle("Eta Relative Resolution; pT (GeV/c); Eta Relative Resolution");
    gPeakRatio->SetTitle("Pion/Eta Mass Ratio; pT (GeV/c); Pion/Eta Mass");
    TCanvas *cEtaPeak = new TCanvas("cEtaPeak", "Eta Peak Position", 800, 600);
    gEtaPeak->Draw("ALPE");
    cEtaPeak->Print("2D_Histogram_Fits.pdf");

    TCanvas *cEtaRes = new TCanvas("cEtaRes", "Eta Relative Resolution", 800, 600);
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

void bgsub(double scale_factor = 1, float polyL = 0.05, float polygauss1L = 0.08, float gauss1L = 0.11, float gauss1R = 0.19, float polygauss1R = 0.3, float polygauss2L = 0.5, float gauss2L = 0.55, float gauss2R = 0.65, float polygauss2R = 0.7, float polyR = 1.0, int startBin = 1, int endBin = -1, int projectionBins = 1, int rebinFactor = 1, int background_scheme = 0, bool fitEtaPeak = false, bool dynamic_left = true)
{
  // code is run with a command like this:
  // root bgsub.C'(1,0.07,0.1,0.11,0.19,0.20,0.48,0.48,0.72,0.72,0.99,16,-1,4,2,true)'
  //  Scale factor for histogram errors
  //  double scale_factor = 1.0;

  // background scheme:0=poly4, 1= poly3 + poly2

  // Fit limits for the polynomial and Gaussian fits
  std::vector<float> limits = {
      polyL, polyR,              //0,1 Polynomial fit range: left and right limits
      gauss1L, gauss1R,          //2-3 First Gaussian fit range: left and right limits
      polygauss1L, polygauss1R,  //4,5 Exclusion zone for left and right polynomials: first gaussian
      gauss2L, gauss2R,          //6,7 Second Gaussian fit range (if fitting eta peak): left and right limits
      polygauss2L, polygauss2R   //8,9 Exclusion zone for left and right polynomials: second gaussian
  };

  // Flag to indicate whether to fit the eta peak
  // bool fitEtaPeak = true;

  // Specify the start and end bins for the projections
  // int startBin = 1;
  // int endBin = 20;  //-1;  // Use -1 to indicate the last bin

  // Specify the number of bins for each projection and the rebin factor
  // int projectionBins = 1;  // Default to 1 bin
  // int rebinFactor = 1;     // Default to no rebinning

  // fit_histogram(scale_factor, leftmost_gauslimit, rightmost_gauslimit, true);
  fit_2d_histogram(scale_factor, limits, fitEtaPeak, startBin, endBin, projectionBins, rebinFactor, dynamic_left, background_scheme);
  // return 0;
  gApplication->Terminate(0);
}
