#include "sPhenixStyle.C"

void myText(Double_t x, Double_t y, Color_t color, const char* text, Doubl e_t tsize = 0.04)
{
  TLatex l;
  l.SetTextAlign(12);
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x, y, text);
}

void plot()
{
  SetsPhenixStyle();

  std::string file_path = "runList.txt";

  std::ifstream file(file_path);

  std::vector<int> run_numbers;
  std::string line;
  while (std::getline(file, line))
  {
    std::istringstream iss(line);
    int run_number;
    if (iss >> run_number)
    {
      run_numbers.push_back(run_number);
    }
    else
    {
      std::cerr << "Warning: Invalid integer on line: " << line << std ::endl;
    }
  }

  file.close();

  if (run_numbers.empty())
  {
    std::cerr << "Error: No valid run numbers found in the file." << std ::endl;
    return 1;
  }

  int first_run = run_numbers.front();
  int last_run = run_numbers.back();

  TFile* finpi0 = new TFile("combine_out/out_fin152.root");

  TH1F* h_InvMass = (TH1F*) finpi0->Get("h_InvMass");
  TH1F* h_InvMassMix = (TH1F*) finpi0->Get("h_InvMassMix");
  // h_InvMass->Rebin(2);
  // h_InvMassMix->Rebin(2);

  //////////////////////////////////
  // mix event
  /*
  int bin1 = h_InvMass->FindBin(0.7);
  int bin2 = h_InvMass->FindBin(1.19);
  h_InvMassMix->Scale(h_InvMass->Integral(bin1,bin2)/h_InvMassMix->Integra
  l(bin1,bin2));

  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  h_InvMass->Draw("ex0");
  h_InvMassMix->Draw("same ex0");
  h_InvMassMix->SetMarkerColor(kRed);
  h_InvMassMix->SetMarkerStyle(22);

  myText(0.2,0.98,1,Form("%d-%d",first_run,last_run));

  c1->SaveAs("figures/sig.pdf");

  TH1F* h_InvMass_sub = (TH1F*) h_InvMass->Clone("h_InvMass_sub");
  h_InvMass_sub->Add(h_InvMassMix,-1);

  TF1* fitFunc1 = new TF1("fitFunc1", "[0]+x*[1]",0,h_InvMassMix->GetXaxis
  ()->GetXmax());
  TF1* fitFunc2 = new TF1("fitFunc2", "[0]+x*[1]",0,h_InvMassMix->GetXaxis
  ()->GetXmax());
  h_InvMass->Fit(fitFunc1,"QRN","",0.4,0.8);
  h_InvMassMix->Fit(fitFunc2,"QRN","",0.4,0.8);

  TH1F* h_mix_adj = (TH1F*) h_InvMassMix->Clone("h_mix_adj");
  for(int i=1; i<h_mix_adj->GetNbinsX()+1; i++) {
      float x = h_mix_adj->GetBinCenter(i);
      float val = h_mix_adj->GetBinContent(i)*fitFunc1->Eval( x) /fitFunc2
  ->Eval(x);
      h_mix_adj->SetBinContent(i,val);
  }
  h_mix_adj->Scale(h_InvMass->Integral(bin1,bin2)/h_mix_adj->Integral(bin1
  ,bin2));
  TH1F* h_InvMass_sub2 = (TH1F*) h_InvMass->Clone("h_InvMass_sub2");
  h_InvMass_sub2->Divide(h_mix_adj);

  TCanvas* c2 = new TCanvas("c2","c2",600,600);
  h_InvMass->Draw();
  h_mix_adj->Draw("same ex0");
  h_mix_adj->SetMarkerColor(kBlue);
  h_mix_adj->SetMarkerStyle(22);

  c2->SaveAs("figures/sig_rwmix.pdf");


  TCanvas* c6 = new TCanvas("c6","c6",600,600);
  h_InvMass_sub2->Draw();

  c6->SaveAs("figures/sig_sub_rwmix.pdf");

*/

  /////////////////////////////////////
  // fitting mass hist
  // vector<double>vtemp = {0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,
  // 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.24, 0.2 6, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5, 0.52, 0.54, 0.5 6, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0. 86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06, 1.08, 1.1, 1.12, 1.14, 1. 16, 1.18, 1.2};
  vector<double> vtemp = {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0. 10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.36, 0.40, 0.44, 0.48, 0.52, 0.56, 0.60, 0.64, 0.68, 0.72, 0.76, 0. 8, 0.84, 0.88, 0.92, 0.96, 1.0, 1.04, 1.08, 1.12, 1.16, 1.2};
  const int Nbins = vtemp.size() - 1;

  double_t bins[100];
  for (int ib = 0; ib < Nbins + 1; ib++) bins[ib] = vtemp[ib];

  TH1F* h = (TH1F*) h_InvMass->Clone("h");
  // h->Rebin(2);

  // h->Rebin(5);
  TH1F* hrb = (TH1F*) h->Rebin(Nbins, "hrb", bins);
  hrb->Scale(1. / 2, "width");
  h = hrb;

  for (int ib = 1; ib < h->GetNbinsX() + 1; ib++) h->SetBinError(ib, h->GetBinError(ib) *1.414);

  TH1F* h_bkg = (TH1F*) hrb->Clone("h_bkg");
  TH1F* h_sig = (TH1F*) hrb->Clone("h_sig");

  //  removes bin with etas in bkg
  int b1 = h_bkg->GetXaxis()->FindBin(0.5);
  int b2 = h_bkg->GetXaxis()->FindBin(0.83);
  for (int ib = b1; ib < b2 + 1; ib++)
  {
    h_bkg->SetBinContent(ib, 0);
    h_bkg->SetBinError(ib, 100000000000000);
  }

  // fit background
  TF1* f_bkg2 = new TF1("fbkgg", "pol3", 0.05, 0.35);
  TF1* f_pi = new TF1("fpi", "gaus", 0.05, 0.35);
  f_pi->SetParLimits(1, 0.13, 0.19);
  f_pi->SetParLimits(2, 0.01, 0.25);
  // f_bkg2->SetParLimits(0,-100000,0);

  TF1* f_pi_fit = new TF1("f_pi_fit", "fbkgg+fpi", 0.05, 0.35);
  f_pi_fit->SetParLimits(6, 0.01, 0.25);
  f_pi_fit->SetParLimits(5, 0.13, 0.19);

  TF1* f_bkg = new TF1("fbkg", "pol2", 0.35, 1.2);
  TF1* f_eta = new TF1("feta", "gaus", 0.35, 1.2);
  f_eta->SetParLimits(1, 0.5, 0.75);

  h->Fit(f_pi_fit, "", "r", 0.05, 0.35);
  h_bkg->Fit(f_bkg, "", "r", 0.35, 1.15);

  TCanvas* c11 = new TCanvas("c11", "c11", 600, 600);
  h_bkg->Draw();

  h_sig->Reset();

  b1 = h_bkg->GetXaxis()->FindBin(0.37);
  b2 = h_bkg->GetXaxis()->FindBin(1.20);

  // subtract bkg
  for (int ib = b1; ib < b2 + 1; ib++)
  {
    float val = h->GetBinContent(ib) - f_bkg->Eval(h->GetBinCenter(ib));
    float err = h->GetBinError(ib);
    h_sig->SetBinContent(ib, val);
    h_sig->SetBinError(ib, err);
  }

  h_sig->Fit(f_eta, "", "r", 0.35, 1.1);

  TGraph* line0;
  line0 = new TGraph(2);
  line0->SetPoint(0, -100, 0);
  line0->SetPoint(1, 500, 0);
  line0->SetLineStyle(2);
  line0->SetLineColor(1);
  line0->SetLineWidth(1);

  f_bkg2->SetParameter(0, f_pi_fit->GetParameter(0));
  f_bkg2->SetParameter(1, f_pi_fit->GetParameter(1));
  f_bkg2->SetParameter(2, f_pi_fit->GetParameter(2));
  f_bkg2->SetParameter(3, f_pi_fit->GetParameter(3));

  TF1* f_fit_eta = new TF1("f_fit_eta", "fbkg+feta", 0.3, 1.25);
  TF1* fp_bkg = (TF1*) f_bkg->Clone("fp_bkg");

  TPad* pad[2];
  pad[0] = new TPad("pad0", "", 0, 0.3, 1.0, 1.0);
  pad[1] = new TPad("pad1", "", 0, 0.0, 1.0, 0.3);

  TCanvas* c3 = new TCanvas("c3", "c3", 600, 600);
  pad[0]->Draw();
  pad[1]->Draw();
  pad[0]->cd();

  h->Draw("ex0");
  h->GetXaxis()->SetRangeUser(0, 1.0);
  h->GetYaxis()->SetRangeUser(1, 9e5);
  h->SetXTitle("#it{M}_{#gamma#gamma}");
  h->SetYTitle("pairs / GeV");
  float bla = h->GetYaxis()->GetTitleSize();

  h->GetYaxis()->SetTitleSize(bla * 1.25);
  h->GetYaxis()->SetLabelSize(bla * 1.25);

  f_fit_eta->Draw("same");
  f_fit_eta->SetLineColor(kRed);
  f_fit_eta->GetXaxis()->SetRangeUser(0.3, 1.3);
  h->SetXTitle("#it{M}_{#gamma#gamma} [GeV]");

  f_pi_fit->Draw("same");
  f_pi_fit->SetLineColor(kRed);
  f_bkg2->Draw("same");
  f_bkg2->SetLineColor(kBlue);

  fp_bkg->Draw("same");
  fp_bkg->SetLineColor(kBlue);

  // gPad->SetLogy();

  float piM = f_pi_fit->GetParameter(5) * 1e3;
  float piW = f_pi_fit->GetParameter(6) * 1e3;
  float piErr = f_pi_fit->GetParError(5) * 1e3;
  float etaErr = f_eta->GetParError(1) * 1e3;
  float etaM = f_eta->GetParameter(1) * 1e3;
  float etaW = f_eta->GetParameter(2) * 1e3;

  TH1F* h_event = (TH1F*) finpi0->Get("h_event");

  myText(0.4, 0.88, 1, "#it{#bf{sPHENIX}} Preliminary", 0.06);
  myText(0.78, 0.88, 1, "June 7, 2024", 0.04);
  myText(0.4, 0.83 - 0.00, 1, "#it{p}+#it{p} 200 GeV", 0.06);
  myText(0.4, 0.77 - 0.00, 1, Form("%0.0fM Events", h_event->GetEntries() / 1e6), 0.06);
myText(0.4, 0.69-0.00, 1, "#it{p}_{T1,2} > 1.5 GeV, #it{A}_{#gamma#gamma}< 0.6 ",0.06);
myText(0.4, 0.62 - 0.00, 1, Form("#it{#mu}_{#pi}/#it{#mu}_{#eta} = %0.3f (
#it{m } _{#pi } / #it{m } _{#eta } = 0.246) ",piM/etaM),0.06);
myText(0.4, 0.55 - 0.00, 1, Form("#it{#sigma}_{#pi}/#it{#mu}_{#pi} = %0.3f",piW/piM),0.06);
myText(0.4, 0.47 - 0.00, 1, Form("#it{#sigma}_{#eta}/#it{#mu}_{#eta} = %0.3f ",etaW/etaM),0.06);

// myText(0.45, 0.69, 1, "#it{A}_{#gamma#gamma}<0.6 ",0.06);
// myText(0.5, 0.62, 1, "#it{p}_{T1}>1.3 GeV",0.06);
// myText(0.5, 0.55, 1, "#it{p}_{T2}>0.7 GeV",0.06);
// myText(0.45, 0.55, 1, Form("#it{#mu}_{#eta} = %0.1f #pm %0.1f MeV",etaM,etaErr),0.06);
// myText(0.45, 0.62, 1, Form("#it{#mu}_{#pi} = %0.1f #pm %0.1f MeV",piM,piErr),0.06);
// myText(0.45, 0.55, 1, Form("#it{#mu}_{#eta} = %0.1f #pm %0.1f MeV",etaM, etaErr),0.06);

// h->GetYaxis()->SetRangeUser(2e4,1e6);

gPad->SetBottomMargin(0.0);
gPad->SetLeftMargin(0.18);
gPad->SetTopMargin(0.07);

c3->SaveAs("figures/pi0_eta_mass_log.pdf");

pad[1]->cd();

gPad->SetLeftMargin(0.18);

gPad->SetTopMargin(0);
gPad->SetBottomMargin(0.3);
h_sig->Draw("ex0");
h_sig->SetXTitle("#it{m}_{#gamma#gamma} [GeV]  ");

h_sig->SetYTitle("data-bkg");
h_sig->GetXaxis()->SetRangeUser(0, 1.0);

h_sig->GetYaxis()->SetNdivisions(505, kTRUE);
h_sig->GetYaxis()->SetLabelSize(0.15);
h_sig->GetYaxis()->SetTitleSize(0.16);
h_sig->GetYaxis()->CenterTitle(kTRUE);
h_sig->GetYaxis()->SetTitleOffset(0.6);

h_sig->GetXaxis()->SetNdivisions(510, kTRUE);
h_sig->GetXaxis()->SetLabelSize(0.15);
h_sig->GetXaxis()->SetTickLength(0.10);
h_sig->GetXaxis()->SetTitleSize(0.16);
h_sig->GetXaxis()->SetTitleOffset(0.9);

line0->Draw("l");
f_eta->Draw("same");
f_eta->SetLineColor(kRed);

c3->SaveAs("figures/pi0_eta_mass.pdf");
}
