#include "TFeldmanCousins.h"

void spectrum()
{

  TFeldmanCousins fldc_err(0.68);
  
  // Opening Energy histogram of data root file
  TFile *fEdata = new TFile("/media/glados/Data_Storage_2/mdtax4/pass5/data/processing/hEdata.root","r");
  TH1D *hEdata = (TH1D*)fEdata->Get("hEdata");

  
  Int_t nxbins = hEdata->GetNbinsX();
  Double_t log10en_lo = hEdata->GetXaxis()->GetXmin();
  Double_t log10en_up = hEdata->GetXaxis()->GetXmax();

  
  // Opening Exposure histogram with fitted function from root file
  TFile *fExp = new TFile("/media/glados/Data_Storage_2/mdtax4/mc/processing/hExposure.root","r");
  TH1D *hExposure = (TH1D*)fExp->Get("hExposure");
  TF1 *func_hExposure = hExposure->GetFunction("fit_func");


  TGraphAsymmErrors *gE3J = new TGraphAsymmErrors(0);
  TGraphAsymmErrors *gJ = new TGraphAsymmErrors(0);

  Int_t npts = 0;
  for (int i=1 ; i <= nxbins ; i++)
    {
      // Defining variables for each bin
      Double_t log_bin_center = hEdata->GetBinCenter(i);
      Double_t log_bin_width = hEdata->GetBinWidth(i);
      Double_t log_bin_lower_edge = (log_bin_center - log_bin_width / 2.0);
      Double_t log_bin_upper_edge = (log_bin_center + log_bin_width / 2.0);
      Double_t bin_center = TMath::Power(10,log_bin_center);
      Double_t bin_lower_edge = TMath::Power(10,log_bin_lower_edge);
      Double_t bin_upper_edge = TMath::Power(10,log_bin_upper_edge);

      if(log_bin_center < 18.0)
	continue;
      
      if(log_bin_center > 20.5)
	continue;

      
      // Each piece of the spectrum calc for each bin
      Double_t nevent = (hEdata->GetBinContent(i));
      Double_t Ec3 = (TMath::Power(bin_center,3))/1e24;
      Double_t dE = (bin_upper_edge - bin_lower_edge);
      Double_t exposure = (hExposure->GetBinContent(i+10));
      
      Double_t n_e1, n_e2;

      if(nevent < 20.0 && nevent > 0)
	{
	  n_e1 = nevent-fldc_err.CalculateLowerLimit(nevent,0.0);
	  n_e2 = fldc_err.CalculateUpperLimit(nevent,0.0) - nevent;
	}
      else if (nevent >= 20)
	{
	  n_e1 = TMath::Sqrt(nevent);
          n_e2 = n_e1;
	}
      else
	{
	  n_e1 = 0;
	  n_e2 = 0;
	}

      gE3J->SetPoint(npts, log_bin_center, Ec3 * nevent / dE / exposure);
      gE3J->SetPointError(npts,0,0, Ec3 * n_e1 / dE / exposure, Ec3 * n_e2 / dE / exposure);

      gJ->SetPoint(npts, log_bin_center, nevent / dE / exposure);
      gJ->SetPointError(npts,0,0, n_e1 / dE / exposure, n_e2 / dE / exposure);     
      npts ++;
      
    }
  
  TCanvas *flux = new TCanvas("flux","flux",900,600);
  flux->SetLogy();
  
  gJ->SetMarkerStyle(21);
  gJ->SetTitle(";log_{10}(E/eV);J(E) (eV^{-1}m^{-2}sr^{-1}s^{-1})");
  gJ->GetYaxis()->SetRangeUser(1e-34,2e-30);
  gJ->Draw("a,e1p");

  TCanvas *spec = new TCanvas("spec","spec",900,600);
  spec->SetLogy();
  
  gE3J->SetMarkerStyle(21);
  gE3J->SetTitle(";log_{10}(E/eV);E^{3}*J(E) / 1x10^{24} (eV^{2}m^{-2}sr^{-1}s^{-1})");
  gE3J->GetYaxis()->SetRangeUser(0.1,10.0);
  gE3J->Draw("a,e1p");

  // const Int_t n = 14;
  // Float_t x[n]={18.55,18.65,18.75,18.85,18.95,19.05,19.15,19.25,19.35,19.45,19.55,19.65,19.75};
  // Float_t y[n]={1.09170730899202,0.311281306692546,1.19152030151854,1.21601951068668,0.934271821789413,0.32077797754299,0.448772377043857,0,3.67965860369104,2.68527233664669,0,5.89757210980019,4.42590037043264};
  // Float_t ey[n]={0.302785129259196,0.17971834621264,0.397173433839514,0.45961217352563,0.467135910894706,0.32077797754299,0.448772377043857,0,1.83982930184552,1.89877427857552,0,4.17021323137637,4.42590037043264};
  // TGraphAsymmErrors *gE3J_greg = new TGraphAsymmErrors(n,x,y,0,0,ey,ey);
  // gE3J_greg->SetMarkerColor(kBlue);
  // gE3J_greg->SetMarkerStyle(3);
  // gE3J_greg->SetLineColor(kBlue);
  // gE3J_greg->Draw("P");
  
  TFile *ta_comp = new TFile("ta_combined_icrc2019.root","r");
  TTree *compspec = (TTree*)ta_comp->Get("compspec");
  
  const Int_t nentries = compspec->GetEntries();
  Double_t E3J[51],BC[51];
  Double_t BinCenter,J,E3,E;
  Int_t num = 0;
  
  compspec->SetBranchAddress("BinCenter",&BinCenter);
  compspec->SetBranchAddress("J",&J);
  
  for (Int_t i=0; i < nentries ; i++){
      compspec->GetEntry(i);
      BC[i] = BinCenter;
      E = TMath::Power(10,BinCenter);
      E3 = TMath::Power(E,3)/1e24;
      E3J[i] = E3*J;
      num++;
    }
  TGraph *icrc_compspec = new TGraph(num,BC,E3J);
  icrc_compspec->SetMarkerColor(kRed);
  icrc_compspec->SetMarkerStyle(20);
  icrc_compspec->Draw("P");
}
