#include <Math/MinimizerOptions.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <fstream>
#include "sPhenixStyle.C"
#include "sPhenixStyle.h"

// Combined function for Gaussian + Polynomial
Double_t combinedFunction(Double_t *x, Double_t *par)
{
  // Gaussian part
  Double_t gauss = par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));

  // Polynomial part
  Double_t poly = par[3] + par[4] * x[0] + par[5] * x[0] * x[0] + par[6] * x[0] * x[0] * x[0] + par[7] * x[0] * x[0] * x[0] * x[0];

  return gauss + poly;
}

Double_t combinedFunctionDoubleGauss(Double_t *x, Double_t *par)
{
  // First Gaussian part (e.g., pion peak)
  Double_t gauss1 = par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));

  // Second Gaussian part (e.g., eta peak)
  Double_t gauss2 = par[8] * exp(-0.5 * pow((x[0] - par[9]) / par[10], 2));

  // Polynomial part
  Double_t poly = par[3] + par[4] * x[0] + par[5] * x[0] * x[0] + par[6] * x[0] * x[0] * x[0] + par[7] * x[0] * x[0] * x[0] * x[0];

  return gauss1 + gauss2 + poly;
}

// Custom function for simultaneous left and right polynomial fit, excluding the Gaussian region
Double_t leftRightPolynomial(Double_t *x, Double_t *par)
{
  // Define the range for the Gaussian fit
  Double_t gauss_min = 0.11;
  Double_t gauss_max = 0.19;

  // Check if x is in the range of the Gaussian fit
  if (x[0] >= gauss_min && x[0] <= gauss_max)
  {
    TF1::RejectPoint();
    return 0;
  }

  // Polynomial (4th degree) calculation
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0] + par[4] * x[0] * x[0] * x[0] * x[0];
}

// scale the histogram's error bars
void scale_histogram_errors(TH1F *hist_error_scale, Double_t scale_factor)
{
  for (int i = 1; i <= hist_error_scale->GetNbinsX(); ++i)
  {
    // Get the current error
    Double_t current_error = hist_error_scale->GetBinError(i);

    // Scale the error by the scale factor
    hist_error_scale->SetBinError(i, current_error * scale_factor);
    // std::cout << "orig bin cont: " << hist_error_scale->GetBinContent(i) << " . bin error: " << current_error << " . New bin error: " << hist_error_scale->GetBinError(i) <<std::endl;
  }
}
// scale the histogram's error bars
void replace_histogram_errors(TH1F *hist_error_replace, Double_t error_replace)
{
  for (int i = 1; i <= hist_error_replace->GetNbinsX(); ++i)
  {
    // Scale the error by the scale factor
    // std::cout << "orig bin cont: " << histCopy_ErrorScale->GetBinContent(i) << " . new bin cont: " << newBinContent_scaled << std::endl;
    hist_error_replace->SetBinError(i, error_replace);
  }
}

void appendtextfile(TF1 *fitFunc, const std::string &fitName, Double_t scale_factor)
{
  // Open a text file in append mode
  std::ofstream outfile;
  outfile.open("fit_parameters.txt", std::ios_base::app);

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

void fit_histogram(Double_t scale_factor = 1, float leftmost_gauslimit= 0.05, float rightmost_gauslimit= 0.3, bool fitEtaPeak = false)
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
  float rightmost_limit = 0.3;  // fit range limit
  float leftmost_limit = 0.05;  // fit range limit. normally 0.05
  // Limits on gauss and poly
  float leftpolylim = 0.11;
  float rightpolylim = 0.19;

  hist->GetXaxis()->SetRangeUser(0, 1.0);
  scale_histogram_errors(hist, scale_factor);

  // Fit left and right regions with a polynomial, excluding Gaussian region
  TF1 *leftRightFit = new TF1("leftRightFit", leftRightPolynomial, leftmost_limit, rightmost_limit, 5);
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
    combinedFit->SetParameter(9, 0.55);                          // Rough guess for eta peak mean
    combinedFit->SetParameter(10, 0.05);                         // Rough guess for eta peak sigma
  }

  // Fit the combined function
  hist->Fit(combinedFit, "RL");

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
  pt->AddText(Form("Peak Mean = %f +/- %f", combinedFit->GetParameter(1), combinedFit->GetParError(1)));
  pt->AddText(Form("Peak Sigma = %f +/- %f", combinedFit->GetParameter(2), combinedFit->GetParError(2)));
  if (fitEtaPeak)
  {
    pt->AddText(Form("Eta Peak Mean = %f +/- %f", combinedFit->GetParameter(9), combinedFit->GetParError(9)));
    pt->AddText(Form("Eta Peak Sigma = %f +/- %f", combinedFit->GetParameter(10), combinedFit->GetParError(10)));
  }
  pt->AddText(Form("Background Subtracted Peak Fit = %f to %f", leftmost_gauslimit, rightmost_gauslimit));
  pt->AddText(Form("Mean = %f +/- %f", gausFit2->GetParameter(1), gausFit2->GetParError(1)));
  pt->AddText(Form("Sigma = %f +/- %f", gausFit2->GetParameter(2), gausFit2->GetParError(2)));
  pt->AddText(Form("Relative Width: %f", gausFit2->GetParameter(2) * 100.0f / gausFit2->GetParameter(1)));
  pt->AddText(Form("Chi2/NDF = %f / %d = %f", gausFit2->GetChisquare(), gausFit2->GetNDF(), gausFit2->GetChisquare() / gausFit2->GetNDF()));

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

void fit_2d_histogram(Double_t scale_factor, float leftmost_gauslimit, float rightmost_gauslimit, bool fitEtaPeak = false, int startBin = 1, int endBin = -1)
{
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
  SetsPhenixStyle();

  // Open the ROOT file and get the 2D histogram
  TFile *file = new TFile("/sphenix/u/nkumar/analysis/EMCal_pi0_Calib_2023/macros/condor/output/merged_file.root");
  TH2F *hist2D = (TH2F *) file->Get("h_InvMass_2d");

  int nBinsX = hist2D->GetNbinsX();
  if (endBin == -1) endBin = nBinsX;  // Default to the last bin if not specified

  // Overall limits
  float rightmost_limit = 0.3;  // fit range limit
  float leftmost_limit = 0.05;  // fit range limit. normally 0.05

  // Limits on gauss and poly
  float leftpolylim = 0.11;
  float rightpolylim = 0.19;

  // Create a PDF to save the canvases
  TCanvas *dummyCanvas = new TCanvas();
  dummyCanvas->Print("2D_Histogram_Fits.pdf[");

  // Vectors to store fit results
  std::vector<double> pionPt, pionPeak, pionRes, etaPeak, etaRes;

  for (int i = startBin; i <= endBin; ++i)
  {
    // Project the histogram along the Y-axis
    TH1D *hist = hist2D->ProjectionY(Form("proj_%d", i), i, i);

    // Check if the projection has enough entries to perform a fit
    if (hist->GetEntries() < 100)
    {  // Adjust the threshold as needed
      delete hist;
      continue;
    }

    // Get the pt range for the current slice
    double pt_min = hist2D->GetXaxis()->GetBinLowEdge(i);
    double pt_max = hist2D->GetXaxis()->GetBinUpEdge(i);
    TString ptRange = Form("pt_%.2f-%.2f_GeV", pt_min, pt_max);

    // Rebin and scale the histogram
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

    // Set histogram range and scale errors
    hist->GetXaxis()->SetRangeUser(0, 1.0);
    // scale_histogram_errors(hist, scale_factor);

    // Fit left and right regions with a polynomial, excluding Gaussian region
    TF1 *leftRightFit = new TF1("leftRightFit", leftRightPolynomial, leftmost_limit, rightmost_limit, 5);
    hist->Fit(leftRightFit, "R");

    // Fit Gaussian in the specified range
    TF1 *gausFit = new TF1("gausFit", "gaus", leftpolylim, rightpolylim);
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
    for (int j = 0; j < 3; ++j) combinedFit->SetParameter(j, gausFit->GetParameter(j));
    for (int j = 3; j < 8; ++j) combinedFit->SetParameter(j, leftRightFit->GetParameter(j - 3));

    if (fitEtaPeak)
    {
      // Set initial guesses for the second Gaussian (eta peak)
      combinedFit->SetParameter(8, gausFit->GetParameter(0) / 3);  // Assume a smaller amplitude
      combinedFit->SetParameter(9, 0.55);                          // Rough guess for eta peak mean
      combinedFit->SetParameter(10, 0.05);                         // Rough guess for eta peak sigma
    }

    // Fit the combined function
    hist->Fit(combinedFit, "RL");

    // Store pion peak position and resolution
    pionPt.push_back((pt_min + pt_max) / 2.0);
    pionPeak.push_back(combinedFit->GetParameter(1));
    pionRes.push_back(combinedFit->GetParameter(2) / combinedFit->GetParameter(1));

    // Store eta peak position and resolution if fitting eta
    if (fitEtaPeak)
    {
      etaPeak.push_back(combinedFit->GetParameter(9));
      etaRes.push_back(combinedFit->GetParameter(10) / combinedFit->GetParameter(9));
    }
    else
    {
      etaPeak.push_back(0);
      etaRes.push_back(0);
    }

    // Create a new function for just the polynomial part
    TF1 *polyPart = new TF1("polyPart", "pol4", leftmost_limit, rightmost_limit);
    for (int j = 0; j < 5; ++j) polyPart->SetParameter(j, combinedFit->GetParameter(j + 3));

    // Create a new histogram to store the subtracted data
    TH1F *histSubtracted = (TH1F *) hist->Clone(Form("histSubtracted_%d", i));
    for (int j = 1; j <= hist->GetNbinsX(); ++j)
    {
      double x = hist->GetBinCenter(j);
      double y = hist->GetBinContent(j) - polyPart->Eval(x);
      histSubtracted->SetBinContent(j, y);
    }

    TF1 *gausFit2 = new TF1("gausFit2", "gaus", leftmost_gauslimit, rightmost_gauslimit);
    for (int j = 0; j < 3; ++j) gausFit2->SetParameter(j, combinedFit->GetParameter(j));
    histSubtracted->Fit(gausFit2, "R");

    // Draw the fits and subtracted histograms
    TCanvas *c1 = new TCanvas(Form("c1_%s", ptRange.Data()), "Fits", 800, 600);
    hist->Draw("E");
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
    c1->Print("2D_Histogram_Fits.pdf");

    TCanvas *c2 = new TCanvas(Form("c2_%s", ptRange.Data()), "Subtracted Peak", 800, 600);
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
    c2->Print("2D_Histogram_Fits.pdf");

    appendtextfile(combinedFit, Form("Combined Fit_%s", ptRange.Data()), scale_factor);
    appendtextfile(gausFit2, Form("subpgaus fit_%s", ptRange.Data()), scale_factor);

    delete hist;
    delete leftRightFit;
    delete gausFit;
    delete combinedFit;
    delete polyPart;
    delete histSubtracted;
    delete gausFit2;
    delete leg1;
    delete leg;
    delete pt2;
    delete c1;
    delete c2;
  }

  // Create TGraphs for the collected fit results
  int nPoints = pionPt.size();
  TGraph *gPionPeak = new TGraph(nPoints, &pionPt[0], &pionPeak[0]);
  TGraph *gPionRes = new TGraph(nPoints, &pionPt[0], &pionRes[0]);

  gPionPeak->SetTitle("Pion Peak Position; pT (GeV/c); Peak Position (GeV/c^2)");
  gPionRes->SetTitle("Pion Relative Resolution; pT (GeV/c); Relative Resolution");

  // Draw TGraphs and add them to the PDF
  TCanvas *cPionPeak = new TCanvas("cPionPeak", "Pion Peak Position", 800, 600);
  gPionPeak->Draw("ALP");
  cPionPeak->Print("2D_Histogram_Fits.pdf");

  TCanvas *cPionRes = new TCanvas("cPionRes", "Pion Relative Resolution", 800, 600);
  gPionRes->Draw("ALP");
  cPionRes->Print("2D_Histogram_Fits.pdf");

  if (fitEtaPeak)
  {
    TGraph *gEtaPeak = new TGraph(nPoints, &pionPt[0], &etaPeak[0]);
    TGraph *gEtaRes = new TGraph(nPoints, &pionPt[0], &etaRes[0]);
    gEtaPeak->SetTitle("Eta Peak Position; pT (GeV/c); Peak Position (GeV/c^2)");
    gEtaRes->SetTitle("Eta Relative Resolution; pT (GeV/c); Relative Resolution");
    TCanvas *cEtaPeak = new TCanvas("cEtaPeak", "Eta Peak Position", 800, 600);
    gEtaPeak->Draw("ALP");
    cEtaPeak->Print("2D_Histogram_Fits.pdf");

    TCanvas *cEtaRes = new TCanvas("cEtaRes", "Eta Relative Resolution", 800, 600);
    gEtaRes->Draw("ALP");
    cEtaRes->Print("2D_Histogram_Fits.pdf");

    // Clean up for eta graphs
    delete cEtaPeak;
    delete cEtaRes;
    delete gEtaPeak;
    delete gEtaRes;
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

void bgsub(Double_t scale_factor = 1, float leftmost_gauslimit = 0.05, float rightmost_gauslimit = 0.3)
{
  fit_histogram(scale_factor, leftmost_gauslimit, rightmost_gauslimit, false);
  fit_2d_histogram(scale_factor, leftmost_gauslimit, rightmost_gauslimit, false, 1, -1);
  // return 0;
}