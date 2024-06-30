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
// try ignoring the errors
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

void fit_histogram(Double_t scale_factor, float leftmost_gauslimit, float rightmost_gauslimit)
{
  // more thorough minimizer for fit
  // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  // Set the global fit strategy
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
  SetsPhenixStyle();

  // Open the ROOT file and get the histogram
  // old file->cluster dependent cuts. sig fraction implement that cut. cut on centrality?
  // TFile *file = new TFile("diClusMass_23726_23746_nomPi0CalibCuts.root");
  // new file. cluster dependent cut removed.
  TFile *file = new TFile("/sphenix/u/nkumar/analysis/EMCal_pi0_Calib_2023/macros/condor/output/merged_file.root");
  TH1F *hist = (TH1F *) file->Get("h_InvMass");

  // Rebin the histogram to have 'numBins' bins
  // First, calculate the rebin factor assuming the histogram's range is 0 to maxXRange
  int numBins = 120;
  int currentNumBins = hist->GetNbinsX();
  double currentXMax = hist->GetXaxis()->GetXmax();
  int rebinFactor = currentNumBins / numBins;
  if (rebinFactor > 1)
  {  // Only rebin if the factor is greater than 1
    std::cout << "current nbins: " << currentNumBins << " requested nbins: " << numBins << " rebin by: " << rebinFactor << std::endl;
    hist->Rebin(rebinFactor);
    std::cout << "new nbin check: " << hist->GetNbinsX() << std::endl;
  }

  // overall limits
  float rightmost_limit = 0.3;  // fit range limit
  float leftmost_limit = 0.05;  // fit range limit. normally 0.05
  // float rightmost_limit= 0.257;// fit range limit
  // float leftmost_limit= 0.07; // fit range limit
  //  limits on gauss and poly
  float leftpolylim = 0.11;
  float rightpolylim = 0.19;
  hist->GetXaxis()->SetRangeUser(0, 0.4);
  // Double_t scale_factor = 2.5; // Replace with the factor by which you want to scale the errors
  // Double_t error_replace= 0.1;
  scale_histogram_errors(hist, scale_factor);
  // replace_histogram_errors(hist, error_replace);

  // Fit left and right regions with a polynomial, excluding Gaussian region
  TF1 *leftRightFit = new TF1("leftRightFit", leftRightPolynomial, leftmost_limit, rightmost_limit, 5);
  hist->Fit(leftRightFit, "R");

  // Fit Gaussian in the specified range
  TF1 *gausFit = new TF1("gausFit", "gaus", leftpolylim, rightpolylim);  // leftpolylim, rightpolylim
  hist->Fit(gausFit, "R");

  // Combined Gaussian + Polynomial fit
  TF1 *combinedFit = new TF1("combinedFit", combinedFunction, leftmost_limit, rightmost_limit, 8);
  // Set initial parameters from previous fits
  for (int i = 0; i < 3; ++i) combinedFit->SetParameter(i, gausFit->GetParameter(i));
  for (int i = 3; i < 8; ++i) combinedFit->SetParameter(i, leftRightFit->GetParameter(i - 3));
  // try to improve the fit.
  hist->Fit(combinedFit, "RL");  // M
  double chi2 = combinedFit->GetChisquare();
  double ndf = combinedFit->GetNDF();
  double chi2ndf = chi2 / ndf;

  std::cout << "Chi-squared: " << chi2 << std::endl;
  std::cout << "Number of Degrees of Freedom: " << ndf << std::endl;
  std::cout << "Chi-squared/NDF: " << chi2ndf << std::endl;

  //-------------------------------------------show the poly4 part seperately

  // Create a new function for just the polynomial part
  TF1 *polyPart = new TF1("polyPart", "pol4", leftmost_limit, rightmost_limit);

  // Set the parameters of polyPart to those from the combined fit
  // Assuming the first 5 parameters of combinedFit are for the polynomial
  for (int i = 0; i < 5; ++i) polyPart->SetParameter(i, combinedFit->GetParameter(i + 3));

  // Create a new histogram to store the subtracted data
  TH1F *histSubtracted = (TH1F *) hist->Clone("histSubtracted");

  // Subtract the polynomial part
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

  // store 2 separate functions for visualization
  TF1 *fleft = new TF1("fleft", leftRightPolynomial, leftmost_limit, leftpolylim, 5);
  fleft->SetParameters(leftRightFit->GetParameters());
  // hist->GetListOfFunctions()->Add(fleft);
  // gROOT->GetListOfFunctions()->Remove(fleft);
  TF1 *fright = new TF1("fright", leftRightPolynomial, rightpolylim, rightmost_limit, 5);
  fright->SetParameters(leftRightFit->GetParameters());
  // hist->GetListOfFunctions()->Add(fright);
  // gROOT->GetListOfFunctions()->Remove(fright);

  // Draw everything
  //-------------------------------------------------------------------------------------------canvas 1
  TCanvas *c1 = new TCanvas("c1", "Fits", 800, 600);

  hist->Draw("E");
  hist->SetTitle("pi0 Inv. Mass Peak; Inv. Mass (GeV); Counts");
  // gausFit->SetLineColor(kRed);
  // gausFit->Draw("SAME");// draw the gaussian fit
  polyPart->SetLineColor(kRed);
  polyPart->Draw("SAME");

  // leftRightFit->SetLineColor(kBlue);
  fleft->SetLineColor(kBlue);
  fright->SetLineColor(kBlue);

  // fleft->Draw("SAME");
  // fright->Draw("SAME");// turn off to see the inflection better.
  // leftRightFit->Draw("SAME"); // Draw the left and right polynomial fits

  combinedFit->SetLineColor(kBlack);
  combinedFit->Draw("SAME");  // draw the combined fit

  // Add a legend
  TLegend *leg1 = new TLegend(0.5, 0.5, 0.95, 0.95);  // bot left x, bot left y, top right x, top right y
  //0.5, 0.2, 0.85, 0.5  bottom right
  leg1->SetFillStyle(0);
  leg1->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
  leg1->AddEntry("", "pythia: p+p #sqrt{s_{NN}} = 200 GeV", "");
  leg1->AddEntry(polyPart, "Background Fit");
  leg1->AddEntry(combinedFit, "Combined Fit");
  leg1->Draw();

  leg1->SetTextAlign(32);  // Center-left alignment
  // Update the canvas
  c1->Update();

  // Save the canvas
  c1->SaveAs("combined_fits.pdf");
  //-------------------------------------------------------------------------------------------canvas 2
  TCanvas *c2 = new TCanvas("c2", "Subtracted Peak", 800, 600);
  histSubtracted->Draw();
  histSubtracted->SetMinimum(0.0);
  histSubtracted->SetTitle("Background Subtracted Peak; Inv. Mass (GeV); Counts (Background subtracted)");
  // c2->SetLeftMargin(0.15); // Increase the left margin. Default is around 0.1.
  histSubtracted->GetYaxis()->SetTitleOffset(1.5);
  float xbleft = 0.5;
  float ybleft = 0.8;
  float xtright = 0.85;
  float ytright = 0.93;
  TLegend *leg = new TLegend(xbleft, ybleft, xtright, ytright);
  leg->SetFillStyle(0);
  leg->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
  leg->AddEntry("", "pythia:p+p #sqrt{s_{NN}} = 200 GeV", "");
  histSubtracted->SetStats(0);

  TPaveText *pt2 = new TPaveText(xbleft + .1, 0.5, xtright, 0.78, "NDC");  // Adjust coordinates as needed
  pt2->SetFillColor(0);                                                    // Set the fill color to 0 for transparency
  pt2->SetFillStyle(0);                                                    // Set fill style to 0 (solid) with color 0 for transparency
  pt2->AddText(Form("#chi^{2}/NDF = %.2f", gausFit2->GetChisquare() / gausFit2->GetNDF()));
  pt2->AddText(Form("Mean = %.4f", gausFit2->GetParameter(1)));
  pt2->AddText(Form("Sigma = %.4f", gausFit2->GetParameter(2)));
  pt2->AddText(Form("Relative Width: %.2f%%", gausFit2->GetParameter(2) * 100.0f / gausFit2->GetParameter(1)));
  pt2->Draw("SAME");
  gPad->Modified();  // Apply the changes to the pad

  leg->Draw("Same");
  gPad->Update();
  c2->SaveAs("Subtracted_Peak.pdf");

  // append // Append fit parameters to text file
  appendtextfile(combinedFit, "Combined Fit", scale_factor);
  appendtextfile(gausFit2, "subpgaus fit", scale_factor);

  // Second canvas: Custom list of fit results
  TCanvas *c3 = new TCanvas("canvas2", "Fit Parameters", 800, 600);
  c3->cd();

  TPaveText *pt = new TPaveText(0.1, 0.1, 0.9, 0.9, "blNDC");  // blNDC: borderless, normalized coordinates
  pt->SetTextAlign(12);                                        // Align text to the left
  pt->SetFillColor(0);                                         // Transparent background

  // Adding custom text entries
  pt->AddText("Data Fit");
  pt->AddText("Fit Parameters:");
  pt->AddText(Form("Combined Fit Range = %f to %f", leftmost_limit, rightmost_limit));
  pt->AddText(Form("Peak Mean = %f +/- %f", combinedFit->GetParameter(1), combinedFit->GetParError(1)));
  pt->AddText(Form("Peak Sigma = %f +/- %f", combinedFit->GetParameter(2), combinedFit->GetParError(2)));
  pt->AddText(Form("Background Subtracted Peak Fit = %f to %f", leftmost_gauslimit, rightmost_gauslimit));
  pt->AddText(Form("Mean = %f +/- %f", gausFit2->GetParameter(1), gausFit2->GetParError(1)));
  pt->AddText(Form("Sigma = %f +/- %f", gausFit2->GetParameter(2), gausFit2->GetParError(2)));
  pt->AddText(Form("Relative Width: %f", gausFit2->GetParameter(2) * 100.0f / gausFit2->GetParameter(1)));
  pt->AddText(Form("Chi2/NDF = %f / %d= %f", gausFit2->GetChisquare(), gausFit2->GetNDF(), gausFit2->GetChisquare() / gausFit2->GetNDF()));

  pt->Draw();
  c3->SaveAs("FitInfo.pdf");

  delete file;
  delete c1;
  delete c2;
  delete c3;
  delete gausFit;
  delete gausFit2;
  delete polyPart;
  delete leftRightFit;
  delete fleft;
  delete fright;
  delete combinedFit;
  delete leg;
  delete leg1;
}

void bgsub(Double_t scale_factor, float leftmost_gauslimit = 0.05, float rightmost_gauslimit = 0.3)
{
  fit_histogram(scale_factor, leftmost_gauslimit, rightmost_gauslimit);
  // return 0;
}